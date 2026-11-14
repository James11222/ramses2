#!/usr/bin/env python3
"""Estimate the GPU memory (VRAM) of a single-GPU RAMSES II run.

Usage:
    vram.py --ngridmax 2^27 [--npartmax 2^27] [--ncachemax 1e4] [--nvar 5]
            [--grav] [--mhd] [--amr] [--clump [--no-halo] [--tree]]
            [--npre 4|8] [--nsubgrid 1] [--no-lean] [--multimass] [--dep-algo 2]

Counts accept integers, floats (1.4e8) or powers (512^3, 2**27; quote the
latter in zsh, e.g. '2**27').
Values mirror the &AMR_PARAMS/&RUN_PARAMS namelist entries and the
NPRE, NSUBGRID, GRAV and MHD compile flags. Pass --nvar only for HYDRO=1
builds, with nvar = 5 + NENER + NPSCAL + NMETAL + NION. Pass --grav for
GRAV=1 builds with the Poisson solver on. --mhd counts ideal MHD only
(the split non-ideal MHD solver allocates more). Stars, sinks, turbulence,
cooling tables and LONGINT particle ids are not counted.

Pure N-body runs (--grav without --nvar, --mhd and --clump) use the lean mesh
of gpu_lean=.true.: rho, phi_old, the refinement flags and the oct swap tables
stop at ngridmax, and so does the multigrid mesh. Pass --no-lean for
gpu_lean=.false. With particles of a single mass and no mass_cut_refine, nref
is not allocated either; pass --multimass otherwise.

Examples:
    vram.py --ngridmax 1.4e8 --ncachemax 6e7 --npartmax 1024^3 --grav --amr
    vram.py --ngridmax 2^27 --npartmax 2^27 --grav --clump --npre 8
    vram.py --ngridmax 2^21 --nvar 5 --mhd --amr
"""
import argparse

def vram(ngridmax, npartmax, ncachemax, nvar=None, grav=False, mhd=False, amr=False, clump=False,
         halo=True, tree=False, npre=4, nsubgrid=1, lean=True, multimass=False, dep_algo=2):
    """
    This function estimates the GPU memory (VRAM) usage for RAMSES II simulations
    based on the number of grids, particles, cache, precision, number of hydro variables,
    usage of AMR, usage of the clump finder, and usage of merger trees.

    Args:
        ngridmax  : Maximum number of grids
        npartmax  : Maximum number of particles
        ncachemax : Maximum number of cache octs
        nvar      : Number of variables (optional, default is None).
                    If provided, it will be used to estimate additional memory usage
                    for storing hydro variables.
        grav      : Boolean indicating whether gravity is used (GRAV=1, Poisson solver on)
        mhd       : Boolean indicating whether ideal MHD is used (MHD=1)
        amr       : Boolean indicating whether AMR is used (True) or not (False)
        clump     : Boolean indicating whether the clump finder is used (clump_finder=.true.)
        halo      : Boolean indicating whether the clump finder runs in halo mode
                    (saddle_threshold > 0, rho_type_clump=1), only with clump=True
        tree      : Boolean indicating whether merger tree particles are used
                    (merger_tree=.true., only with clump=True, ntreemax = npartmax/100)
        npre      : Precision (4 for single, 8 for double)
        nsubgrid  : NSUBGRID compile flag (octs per neighbour block edge)
        lean      : gpu_lean namelist parameter. It only acts in pure N-body runs
                    (grav, no hydro, no MHD, no clump finder): the ghost octs of the
                    cache then have no rho, nref, phi_old and flags
        multimass : Particles of several masses, or mass_cut_refine set: nref is
                    kept in lean runs (it is always kept otherwise)
        dep_algo  : Largest *_dep_algo of the particle species: the small-stencil
                    deposit (3) keeps 8 more bytes per particle (idp_swap)

    Returns:
        Estimated VRAM usage in GiB (2^30 bytes) assuming 1 GPU [GiB]


    Written by James Sunseri (Princeton University, October 2026)
    """

    G = ngridmax
    C = ncachemax
    N = G + C
    Np = max(npartmax, 1e4)
    P = npre
    M = G/7 + max(C/7, 1e4)
    Nt = 0
    if clump and tree:
        Nt = npartmax/100
    Ns = max(Np, Nt)       # particle swap/sort buffers are sized for the largest species
    Nscan = max(N, Np, Nt) # prefix_sum: largest oct or particle scan (no more 8 cells/oct)


    # Lean mesh (gpu_lean): pure N-body runs keep rho, nref, phi_old, flag1, flag2,
    # the oct swap tables and the multigrid mesh out of the ghost octs of the cache
    lean = lean and grav and nvar is None and not mhd and not clump
    Nsrc = G if lean else N
    Mmg = G/7 if lean else M

    V = (56*N          # grid 32, hash 24
        + 8*Nsrc           # swap_local, swap_global
        + 4*(nsubgrid + 2)**3*G/nsubgrid**3  # nbor
        + 4*Nscan          # prefix_sum
        + (7*P + 8)*Np     # xp, vp, mp, levelp, sortp
        + (P + 4)*Ns       # xp_swap, isp_swap (the particle sort keeps its keys in xp_swap and prefix_sum)
        )
    if dep_algo == 3:
        V += 8*Ns          # idp_swap: Hilbert keys of the small-stencil deposit

    if grav:
        V += 32*P*N        # phi, f
        V += 16*P*Nsrc     # rho, phi_old
        if not lean or multimass:
            V += 8*P*Nsrc  # nref
        V += (56 + 32*P)*Mmg  # multigrid mesh: grid, hash, phi, f
        V += 4*(G + M) + 108*G/7  # multigrid mesh: father, nbor
    if mhd:
        V += 96*P*N  # bold + bnew (8 cells x 6 faces each, two copies)

    if nvar is not None:
        V += 16 * nvar * P * N # uold + unew (8 cells × nvar each, two copies)
                               # nvar = 5 + NENER + NPSCAL + NMETAL + NION
    if amr:
        V += 64 * Nsrc + 4 * N + 4 * Np  # flag1/flag2, father, idp (nlevelmax > levelmin)
    elif clump:
        V += 64 * N  # flag1/flag2 are also allocated in unigrid runs with the clump finder

    if clump:  # gpu_allocate_clump is persistent from startup
        nclumpmax = G/10
        V += 8*Np  # DM pid/hid
        V += (32 + P)*Nscan  # test-cell arrays, sized ntestmax like the scan buffers
        V += (40 + 21*P)*nclumpmax  # per-clump peak arrays
        if halo:
            V += (4 + 201*P)*nclumpmax  # ind_halo, halo_mass, mass_bin + phi_clump (nbin=100)
        if tree:
            V += (9*P + 28)*Nt  # tree xp, vp, mp, tp, tm, levelp, sortp, pid, hid, idt, idp, idm
            V += 8*nclumpmax  # ntree, min_tree_id

    return (V / (1024)**3) + 1.5  # Convert bytes to GiB + CUDA Context


def count(s):
    """Parse a count given as 1000, 1.4e8, 512^3 or 2**27."""
    s = s.replace("**","^")
    try:
        if "^" in s:
            base, exp = s.split("^",1)
            return int(float(base)**int(exp))
        return int(float(s))
    except ValueError:
        raise argparse.ArgumentTypeError("invalid count: "+s)

if __name__=="__main__":
    parser = argparse.ArgumentParser(
        description="Estimate the GPU memory (VRAM) of a single-GPU RAMSES II run",
        epilog="Counts accept integers, floats (1.4e8) or powers (512^3, 2**27).")
    parser.add_argument("--ngridmax", type=count, required=True,
                        help="maximum number of grids (octs)")
    parser.add_argument("--npartmax", type=count, default=0,
                        help="maximum number of particles (default: 0)")
    parser.add_argument("--ncachemax", type=count, default=10000,
                        help="maximum number of cache octs (default: 10000)")
    parser.add_argument("--nvar", type=int, default=None,
                        help="number of hydro variables, HYDRO=1 only (default: none)")
    parser.add_argument("--grav", action="store_true",
                        help="gravity on (GRAV=1 build with the Poisson solver)")
    parser.add_argument("--mhd", action="store_true",
                        help="ideal MHD on (MHD=1 build)")
    parser.add_argument("--amr", action="store_true",
                        help="AMR run (nlevelmax > levelmin)")
    parser.add_argument("--clump", action="store_true",
                        help="clump finder on (clump_finder=.true.)")
    parser.add_argument("--no-halo", dest="halo", action="store_false",
                        help="clump finder not in halo mode, only with --clump")
    parser.add_argument("--tree", action="store_true",
                        help="merger tree particles on (merger_tree=.true.), only with --clump")
    parser.add_argument("--npre", type=int, default=4, choices=[4,8],
                        help="float precision in bytes, NPRE (default: 4)")
    parser.add_argument("--nsubgrid", type=int, default=1,
                        help="NSUBGRID compile flag (default: 1)")
    parser.add_argument("--no-lean", dest="lean", action="store_false",
                        help="gpu_lean=.false.: pure N-body runs also keep rho, nref, phi_old and flags in ghost octs")
    parser.add_argument("--multimass", action="store_true",
                        help="particles of several masses or mass_cut_refine set: lean runs keep nref")
    parser.add_argument("--dep-algo", type=int, default=2, choices=[1,2,3],
                        help="largest *_dep_algo of the particle species (default: 2)")
    args = parser.parse_args()
    if (args.tree or not args.halo) and not args.clump:
        parser.error("--tree and --no-halo require --clump")

    v = vram(args.ngridmax, args.npartmax, args.ncachemax, nvar=args.nvar,
             grav=args.grav, mhd=args.mhd, amr=args.amr, clump=args.clump,
             halo=args.halo, tree=args.tree, npre=args.npre, nsubgrid=args.nsubgrid,
             lean=args.lean, multimass=args.multimass, dep_algo=args.dep_algo)
    print("ngridmax  = {:d}".format(args.ngridmax))
    print("npartmax  = {:d}".format(args.npartmax))
    print("ncachemax = {:d}".format(args.ncachemax))
    print("nvar      = {}".format(args.nvar if args.nvar is not None else "none"))
    print("npre      = {:d} ({} precision)".format(args.npre, "single" if args.npre==4 else "double"))
    print("nsubgrid  = {:d}".format(args.nsubgrid))
    print("grav      = {}".format(args.grav))
    print("mhd       = {}".format(args.mhd))
    print("amr       = {}".format(args.amr))
    print("clump     = {}".format(args.clump))
    lean = args.lean and args.grav and args.nvar is None and not args.mhd and not args.clump
    print("lean      = {}".format(lean))
    if lean:
        print("multimass = {}".format(args.multimass))
    print("dep_algo  = {:d}".format(args.dep_algo))
    if args.clump:
        print("halo      = {}".format(args.halo))
        print("tree      = {}".format(args.tree))
    print("Estimated VRAM: {:.2f} GiB (incl. 1.5 GiB CUDA context)".format(v))
