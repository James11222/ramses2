#!/usr/bin/env python3
"""Estimate the GPU memory (VRAM) of a single-GPU RAMSES II run.

Usage:
    vram.py --ngridmax 2^27 --npartmax 2^27 [--ncachemax 1e4] [--nvar 5]
            [--amr] [--clump [--no-halo] [--tree]] [--npre 4|8]

Counts accept integers, floats (1.4e8) or powers (512^3, 2**27).
Values mirror the &AMR_PARAMS/&RUN_PARAMS namelist entries and the
NPRE compile flag. Pass --nvar only for HYDRO=1 builds, with
nvar = 5 + NENER + NPSCAL + NMETAL + NION.

Examples:
    vram.py --ngridmax 1.4e8 --ncachemax 6e7 --npartmax 1024^3 --amr
    vram.py --ngridmax 2^27 --npartmax 2^27 --clump --npre 8
"""
import argparse

def vram(ngridmax, npartmax, ncachemax, nvar=None, amr=False, clump=False, halo=True, tree=False, npre=4):
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
        amr       : Boolean indicating whether AMR is used (True) or not (False)
        clump     : Boolean indicating whether the clump finder is used (clump_finder=.true.)
        halo      : Boolean indicating whether the clump finder runs in halo mode
                    (saddle_threshold > 0, rho_type_clump=1), only with clump=True
        tree      : Boolean indicating whether merger tree particles are used
                    (merger_tree=.true., only with clump=True, ntreemax = npartmax/100)
        npre      : Precision (4 for single, 8 for double)

    Returns:
        Estimated VRAM usage in GB assuming 1 GPU [GBs]


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

    V = ((96 + 56*P)*N + 108*G  # grid 64, hash 24, swap 8, gravity rho/phi/nref/phi_old/f
        + (92 + 32*P)*M + 4*G + 108*G/7  # multigrid mesh: grid, hash, phi/f, father, nbor
        + 4*Nscan          # prefix_sum
        + (7*P + 8)*Np     # xp, vp, mp, levelp, sortp
        + (P + 32)*Ns      # xp_swap, isp_swap, idp_swap, hkeyp + CUB radix sort workspace
        )

    if nvar is not None:
        V += 16 * nvar * P * N # uold + unew (8 cells × nvar each, two copies)
                               # nvar = 5 + NENER + NPSCAL + NMETAL + NION
    if amr:
        V += 68 * N + 4 * Np  # flag1/flag2, father, idp (nlevelmax > levelmin)
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

    return (V / (1024)**3) + 1.5  # Convert bytes to GB + CUDA Context


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
    args = parser.parse_args()
    if (args.tree or not args.halo) and not args.clump:
        parser.error("--tree and --no-halo require --clump")

    v = vram(args.ngridmax, args.npartmax, args.ncachemax, nvar=args.nvar,
             amr=args.amr, clump=args.clump, halo=args.halo, tree=args.tree,
             npre=args.npre)
    print("ngridmax  = {:d}".format(args.ngridmax))
    print("npartmax  = {:d}".format(args.npartmax))
    print("ncachemax = {:d}".format(args.ncachemax))
    print("nvar      = {}".format(args.nvar if args.nvar is not None else "none"))
    print("npre      = {:d} ({} precision)".format(args.npre, "single" if args.npre==4 else "double"))
    print("amr       = {}".format(args.amr))
    print("clump     = {}".format(args.clump))
    if args.clump:
        print("halo      = {}".format(args.halo))
        print("tree      = {}".format(args.tree))
    print("Estimated VRAM: {:.2f} GB (incl. 1.5 GB CUDA context)".format(v))
