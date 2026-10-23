#!/usr/bin/env python3
import yt
from yt.utilities.exceptions import YTFieldNotFound
import numpy as np
import ramses as ram
import argparse
import os
from scipy.io import FortranFile
from scipy.spatial import KDTree
from scipy.integrate import cumulative_trapezoid

#================================================
# This code convert an old ramses output folder
# into a new mini-ramses output folder.
# It creates the following files if relevant:
# params.bin, part.00001, star.00001,
# amr.00001, hydro.00001.
# Contributions from Corentin Cadiou.
# Romain Teyssier, Princeton, Oct 22 2025.
#================================================

#================================
# read old ramses parameter files
#================================
def rd_params(nout,**kwargs):
    path = kwargs.get("path","./")
    car1 = str(nout).zfill(5)
    filename = path+"/output_"+car1+"/amr_"+car1+".out00001"

    with FortranFile(filename, 'r') as f:
        ncpu, = f.read_ints('i')
        ndim, = f.read_ints('i')
        nx,ny,nz = f.read_ints('i')
        nlevelmax, = f.read_ints('i')
        ngridmax, = f.read_ints('i')
        nboundary, = f.read_ints('i')
        ngrid_current, = f.read_ints('i')
        boxlen, = f.read_reals('f8')

        noutput,iout,ifout = f.read_ints('i')
        tout = f.read_reals('f8')
        aout = f.read_reals('f8')
        t, = f.read_reals('f8')
        dtold = f.read_reals('f8')
        dtnew = f.read_reals('f8')
        nstep,nstep_coarse = f.read_ints('i')
        einit,mass_tot_0,rho_tot = f.read_reals('f8')
        omega_m,omega_l,omega_k,omega_b,h0,aexp_ini,boxlen_ini = f.read_reals('f8')
        aexp,hexp,aexp_old,epot_tot_int,epot_tot_old = f.read_reals('f8')
        mass_sph, = f.read_reals('f8')

    params = {"noutput":noutput,"iout":iout,"ifout":ifout,"tout":tout,"aout":aout,"t":t,
              "dtold":dtold,"dtnew":dtnew,"nstep":nstep,"nstep_coarse":nstep_coarse,
              "einit":einit,"mass_tot_0":mass_tot_0,"rho_tot":rho_tot,
              "aexp_ini":aexp_ini,"boxlen_ini":boxlen_ini,"aexp":aexp,"hexp":hexp,
              "aexp_old":aexp_old,"epot_tot_int":epot_tot_int,"epot_tot_old":epot_tot_old,
              "mass_sph":mass_sph}

    return params

#================================
# write new ramses header file
#================================
def wr_header(fileloc,npart,nfile,fields):
    with open(fileloc, 'w') as f:
        # Write header information
        f.write("Total number of particles\n")
        f.write(f"{npart}\n")
        f.write("Total number of files\n")
        f.write(f"{nfile}\n")
        # Particle fields: name, numpy data type, number of components
        f.write("Particle fields\n")
        for name,dtype,ncomp in fields:
            f.write(f"{name} {dtype} {ncomp}\n")

#================================
# conformal to proper time
#================================
def tau2t(tau,params1):
    # same integrals as friedman in init_time.f90,
    # both times are in units of 1/H0 and zero today
    om,ol,ok = params1["omega_m"],params1["omega_l"],params1["omega_k"]
    lna = np.linspace(np.log(1e-6),np.log(2.0),200001)
    a = np.exp(lna)
    dtau = a/np.sqrt(a**3*(om+ol*a**3+ok*a))
    dt = a/np.sqrt((om+ol*a**3+ok*a)/a)
    tau_tab = cumulative_trapezoid(dtau,lna,initial=0)
    t_tab = cumulative_trapezoid(dt,lna,initial=0)
    tau_tab = tau_tab - np.interp(0.0,lna,tau_tab)
    t_tab = t_tab - np.interp(0.0,lna,t_tab)
    return np.interp(tau,tau_tab,t_tab)

#================================
# write new ramses particle file
#================================
def wr_part(data,ptype,mask,ndim,prefix,star,conformal,params1,id_type):
    if mask is None:
        mask = slice(None)
    # yt code_length is unit_l*boxlen, new ramses positions are in [0,boxlen]
    boxlen = params1["boxlen"]
    x = (data[ptype,"particle_position_x"].to("code_length").v[mask]*boxlen).astype(np.float32)
    y = (data[ptype,"particle_position_y"].to("code_length").v[mask]*boxlen).astype(np.float32)
    z = (data[ptype,"particle_position_z"].to("code_length").v[mask]*boxlen).astype(np.float32)
    vx = data[ptype,"particle_velocity_x"].to("code_velocity").v[mask].astype(np.float32)
    vy = data[ptype,"particle_velocity_y"].to("code_velocity").v[mask].astype(np.float32)
    vz = data[ptype,"particle_velocity_z"].to("code_velocity").v[mask].astype(np.float32)
    m = data[ptype,"particle_mass"].to("code_mass").v[mask].astype(np.float32)
    level = data[ptype,"particle_level"].v[mask].astype(np.int32)
    idp = data[ptype,"particle_identity"].v[mask]
    if star:
        # fields are absent if the run has no stars
        try:
            metal = data[ptype,"particle_metallicity"].v[mask].astype(np.float32)
            # new ramses birth date is proper time
            if conformal:
                birth = tau2t(data[ptype,"conformal_birth_time"].v[mask],params1).astype(np.float32)
            else:
                birth = data[ptype,"particle_birth_time"].v[mask].astype(np.float32)
        except YTFieldNotFound:
            metal = np.zeros(len(x),dtype=np.float32)
            birth = np.zeros(len(x),dtype=np.float32)

    npart = len(x)
    print("Found ",npart," particles")

    fields = [("pos","f4",ndim),("vel","f4",ndim),("mass","f4",1)]
    if star:
        fields += [("metallicity","f4",1),("birth_date","f4",1)]
    fields += [("level","i4",1),("birth_id",np.dtype(id_type).str[1:],1)]
    wr_header(prefix+"_header.txt",npart,1,fields)

    with open(prefix+".00001", "wb") as f:
        np.array([ndim],dtype=np.int32).tofile(f)
        np.array([npart],dtype=np.int32).tofile(f)
        x.tofile(f)
        y.tofile(f)
        z.tofile(f)
        vx.tofile(f)
        vy.tofile(f)
        vz.tofile(f)
        m.tofile(f)
        if star:
            metal.tofile(f)
            birth.tofile(f)
        level.tofile(f)
        to_id(idp,id_type).tofile(f)

#================================
# check and cast particle ids
#================================
def to_id(idp,id_type):
    if len(idp) > 0 and np.abs(idp).max() > np.iinfo(id_type).max:
        raise SystemExit("Particle ids overflow int32: use --longint and compile with -DLONGINT")
    return idp.astype(id_type)

#================================
# write new ramses sink file
#================================
def wr_sink(csvfile,ndim,prefix,conformal,params1,id_type):
    # old sink csv: first comment line holds the column names
    with open(csvfile) as f:
        lines = f.readlines()
    names = lines[0].strip().lstrip("#").strip().split(",")
    col = {name:i for i,name in enumerate(names)}
    rows = [line for line in lines if not line.strip().startswith("#")]
    s = np.loadtxt(rows,delimiter=",",ndmin=2).reshape(-1,len(names))

    nsink = s.shape[0]
    print("Found ",nsink," sinks")
    # sink mass, the SMBH mass (mbh) is not used
    m = s[:,col["msink"]]
    tform = s[:,col["tform"]]
    if conformal:
        tform = tau2t(tform,params1)
    # new ramses stores the specific angular momentum
    jsink = s[:,[col["lx"],col["ly"],col["lz"]]]/np.maximum(m,1e-300)[:,None]

    fields = [("pos","f4",ndim),("vel","f4",ndim),("mass","f4",1),
              ("accel","f4",ndim),("angmom","f4",ndim),
              ("birth_date","f4",1),("merging_date","f4",1),("level","i4",1),
              ("birth_id",np.dtype(id_type).str[1:],1),("merging_id",np.dtype(id_type).str[1:],1)]
    wr_header(prefix+"_header.txt",nsink,1,fields)

    with open(prefix+".00001", "wb") as f:
        np.array([ndim],dtype=np.int32).tofile(f)
        np.array([nsink],dtype=np.int32).tofile(f)
        for c in ["x","y","z"][:ndim]:
            s[:,col[c]].astype(np.float32).tofile(f)
        for c in ["vx","vy","vz"][:ndim]:
            s[:,col[c]].astype(np.float32).tofile(f)
        m.astype(np.float32).tofile(f)
        # old force, overwritten at the first move step (as in input_part_ascii.f90)
        np.zeros(nsink*ndim,dtype=np.float32).tofile(f)
        for idim in range(0,ndim):
            jsink[:,idim].astype(np.float32).tofile(f)
        tform.astype(np.float32).tofile(f)
        # not merged yet, as in sink_formation.f90
        np.full(nsink,-1000,dtype=np.float32).tofile(f)
        s[:,col["level"]].astype(np.int32).tofile(f)
        to_id(s[:,col["id"]].astype(np.int64),id_type).tofile(f)
        np.zeros(nsink,dtype=id_type).tofile(f)

#================================
# write new ramses parameter file
#================================
def wr_params(params1,params2,gamma,cosmo,nout,**kwargs):
    path = kwargs.get("path","./")
    car1 = str(nout).zfill(5)
    filename = path+"/output_"+car1+"/params.bin"

    # proper time from conformal time, as in init_time.f90
    if cosmo:
        texp = float(tau2t(params2["t"],params1))
    else:
        texp = params2["t"]

    with open(filename, 'wb') as f:  # 'wb' = write binary mode
        int_array = np.array([1,1,params1["ndim"],params1["levelmin"],
                              params1["levelmax"]],dtype=np.int32)
        int_array.tofile(f)
        boxlen = np.array([params1["boxlen"]],dtype=np.float64)
        boxlen.tofile(f)
        int_array = np.array([params2["noutput"],params2["iout"],params2["ifout"],0],dtype=np.int32)
        int_array.tofile(f)
        dble_array = np.array(params2["tout"],dtype=np.float64)
        dble_array.tofile(f)
        dble_array = np.array(params2["aout"],dtype=np.float64)
        dble_array.tofile(f)
        dble_array = np.array(params2["t"],dtype=np.float64)
        dble_array.tofile(f)
        dble_array = np.array(params2["dtold"],dtype=np.float64)
        dble_array.tofile(f)
        dble_array = np.array(params2["dtnew"],dtype=np.float64)
        dble_array.tofile(f)
        int_array = np.array([params2["nstep"],params2["nstep_coarse"]],dtype=np.int32)
        int_array.tofile(f)
        dble_array = np.array([params2["einit"],params2["mass_tot_0"],params2["rho_tot"]],dtype=np.float64)
        dble_array.tofile(f)
        dble_array = np.array([params1["omega_m"],params1["omega_l"],params1["omega_k"],
                               params1["omega_b"],params1["H0"],params2["aexp_ini"],params2["boxlen_ini"]],dtype=np.float64)
        dble_array.tofile(f)
        dble_array = np.array([params2["aexp"],texp,params2["hexp"],
                               params2["aexp_old"],params2["epot_tot_int"],params2["epot_tot_old"],params2["mass_sph"],gamma],dtype=np.float64)
        dble_array.tofile(f)
        int_array = np.array([1,2,3,4,5,6],dtype=np.int64)
        int_array.tofile(f)

#================================
# write new ramses hydro header
#================================
def wr_hydro_header(fileloc,hydro_fields):
    names = ["density","velocity_x","velocity_y","velocity_z","thermal_pressure"]
    inener = 0
    iscal = 0
    for field in hydro_fields[5:]:
        if "non_thermal" in field:
            inener = inener + 1
            names.append(f"non_thermal_pressure_{inener}")
        else:
            iscal = iscal + 1
            if field == "Metallicity":
                names.append("metal_mass_fraction")
            else:
                names.append(f"scalar_{iscal}")
    with open(fileloc, 'w') as f:
        f.write(f"nvar        ={len(names):11d}\n")
        for ivar,name in enumerate(names):
            f.write(f"variable #{ivar+1:2d}: {name}\n")

#================================
# write new ramses info file
#================================
def wr_info(params1,params2,gamma,cosmo,nout,**kwargs):
    path = kwargs.get("path","./")
    car1 = str(nout).zfill(5)
    filename = path+"/output_"+car1+"/info.txt"

    # proper time from conformal time, as in init_time.f90
    if cosmo:
        texp = float(tau2t(params2["t"],params1))
    else:
        texp = params2["t"]

    with open(filename, 'w') as f:
        f.write(f"nfile       ={1:11d}\n")
        f.write(f"ncpu        ={1:11d}\n")
        f.write(f"ndim        ={params1['ndim']:11d}\n")
        f.write(f"levelmin    ={params1['levelmin']:11d}\n")
        f.write(f"levelmax    ={params1['levelmax']:11d}\n")
        f.write(f"ngridmax    ={params1['ngridmax']:11d}\n")
        f.write(f"nstep_coarse={params2['nstep_coarse']:11d}\n")
        f.write("\n")
        f.write(f"boxlen      ={params1['boxlen']:23.15E}\n")
        f.write(f"time        ={params2['t']:23.15E}\n")
        f.write(f"texp        ={texp:23.15E}\n")
        f.write(f"aexp        ={params2['aexp']:23.15E}\n")
        f.write(f"H0          ={params1['H0']:23.15E}\n")
        f.write(f"omega_m     ={params1['omega_m']:23.15E}\n")
        f.write(f"omega_l     ={params1['omega_l']:23.15E}\n")
        f.write(f"omega_k     ={params1['omega_k']:23.15E}\n")
        f.write(f"omega_b     ={params1['omega_b']:23.15E}\n")
        f.write(f"gamma       ={gamma:23.15E}\n")
        f.write(f"unit_l      ={params1['unit_l']:23.15E}\n")
        f.write(f"unit_d      ={params1['unit_d']:23.15E}\n")
        f.write(f"unit_t      ={params1['unit_t']:23.15E}\n")
        f.write("\n")

#================================
# main code
#================================

#================================
# argument parsing
#================================
parser = argparse.ArgumentParser()
parser.add_argument("nout", help="enter output number")
parser.add_argument("--inp", help="specify the input path")
parser.add_argument("--out", help="specifiy the output path")
parser.add_argument("--dmo", help="convert only DM particles",action="store_true")
parser.add_argument("--longint", help="write 8-byte particle ids (new ramses compiled with -DLONGINT)",action="store_true")
args = parser.parse_args()
id_type = np.int64 if args.longint else np.int32
nout = args.nout
print("Reading output number",nout)
path_in = args.inp
path_out = args.out
dm_only = args.dmo

if path_in==None:
    path_in="."

if dm_only==None:
    dmo_only=False

if path_out==None:
    path_out="new"

car1 = str(nout).zfill(5)
filename=path_in+"/output_"+car1+"/info_"+car1+".txt"

yt.set_log_level(50) # yt verbosity low

#================================
# load up general data set
#================================
params2 = rd_params(nout,path=path_in)
ds = yt.load(filename)
ndim = ds.parameters["ndim"]
levelmin = ds.parameters["levelmin"]
levelmax = ds.parameters["levelmax"]
print("levelmin=",levelmin)
print("levelmax=",levelmax)

# hydro fields in file order (ds.field_list is unordered)
fields = []
for fh in ds.index.domains[0].field_handlers:
    if fh.ftype == "ramses":
        fields = [field[1] for field in fh.field_list]
print("fields=",fields)
has_hydro = "Density" in fields

# new ramses order: density, velocity, pressure,
# non-thermal pressures, passive scalars.
# Magnetic fields are ignored.
hydro_fields = []
if has_hydro:
    if ndim != 3:
        raise SystemExit("Only 3D hydro outputs are supported")
    prim = ["Density","x-velocity","y-velocity","z-velocity","Pressure"]
    bfield = [f for f in fields if f.startswith("B_")]
    nener = [f for f in fields if "non_thermal" in f]
    scalars = [f for f in fields if f not in prim+bfield+nener]
    hydro_fields = prim + nener + scalars
    if len(bfield) > 0:
        print("Ignoring magnetic fields",bfield)
    print("Writing hydro fields",hydro_fields)
    print("New ramses must be compiled with NVAR=",len(hydro_fields)," NENER=",len(nener))
nvar = len(hydro_fields)

# gamma from the old hydro file header, set by yt
gamma = ds.gamma if has_hydro else 5/3
print("gamma=",gamma)

if not os.path.isdir(path_out):
    os.makedirs(path_out)
os.makedirs(path_out+"/output_"+car1,exist_ok=True)

wr_params(ds.parameters,params2,gamma,ds.cosmological_simulation,nout,path=path_out)
print("Parameter file written.")

wr_info(ds.parameters,params2,gamma,ds.cosmological_simulation,nout,path=path_out)
print("Info file written.")

#================================
# read and write particle files
#================================

data = ds.all_data()

ptypes = [field[0] for field in ds.field_list if field[1] == "particle_position_x"]
if ("io","particle_family") in ds.field_list:
    # particle families: use yt particle filters,
    # cloud and tracer particles are not converted
    dm_type, dm_mask = "DM", None
    star_type, star_mask = "star", None
elif "io" in ptypes:
    # no particle families: stars have a non-zero birth time
    if ("io","conformal_birth_time") in ds.field_list:
        birth = data["io","conformal_birth_time"].v
    elif ("io","particle_birth_time") in ds.field_list:
        birth = data["io","particle_birth_time"].v
    else:
        birth = np.zeros(len(data["io","particle_mass"]))
    dm_type, dm_mask = "io", birth == 0
    star_type, star_mask = "io", birth != 0
else:
    dm_type = None
    print("No particles found")

if dm_type is not None:
    print("Loading and writing DM particle data")
    wr_part(data,dm_type,dm_mask,ndim,path_out+"/output_"+car1+"/part",False,
            ds.use_conformal_time,ds.parameters,id_type)

    if not dm_only:
        print("Loading and writing star particle data")
        wr_part(data,star_type,star_mask,ndim,path_out+"/output_"+car1+"/star",True,
                ds.use_conformal_time,ds.parameters,id_type)

# old sinks are stored in a csv file, formation time is conformal in cosmo runs
file_sink_old = path_in+"/output_"+car1+"/sink_"+car1+".csv"
if os.path.exists(file_sink_old) and not dm_only:
    print("Loading and writing sink particle data")
    wr_sink(file_sink_old,ndim,path_out+"/output_"+car1+"/sink",
            ds.cosmological_simulation,ds.parameters,id_type)

#===================================
# read and write amr and hydro file
#===================================

file_amr = path_out+"/output_"+car1+"/amr.00001"
file_hydro = path_out+"/output_"+car1+"/hydro.00001"

ngrid = np.zeros(levelmax,dtype=np.int32)
ngridtot = np.sum(ngrid)

with open(file_amr, "wb") as f_amr:

    # write headers with zero grid counts for all levels,
    # actual counts are filled in level by level below
    np.array([ndim,levelmin,levelmax],dtype=np.int32).tofile(f_amr)
    np.zeros(levelmax-levelmin+1,dtype=np.int32).tofile(f_amr)
    if has_hydro:
        wr_hydro_header(path_out+"/output_"+car1+"/hydro_header.txt",hydro_fields)
        f_hydro = open(file_hydro, "wb")
        np.array([ndim,nvar,levelmin,levelmax],dtype=np.int32).tofile(f_hydro)
        np.zeros(levelmax-levelmin+1,dtype=np.int32).tofile(f_hydro)

    dxmin = 1/2**(np.linspace(0,levelmax,levelmax+1))
    for ilevel in range(levelmin,levelmax+1):
        print("Reading data for level=",ilevel)

        # read current level data
        ds = yt.load(filename,max_level=ilevel,max_level_convention="ramses")
        data = ds.all_data()
        if has_hydro:
            data.get_data([("ramses",f) for f in hydro_fields])
        # yt code_length is unit_l*boxlen: cells span [0,1] for any boxlen
        dx = data["index","dx"].to("code_length").v
        x = data["index","x"].to("code_length").v
        y = data["index","y"].to("code_length").v
        z = data["index","z"].to("code_length").v
        q = [data["ramses",f].v for f in hydro_fields]
        ind = np.where(dx < 1.1*dxmin[ilevel])
        ngrid[ilevel-1] = int(len(ind[0])/8)
        ngridtot = ngridtot + ngrid[ilevel-1]
        print("Level ",ilevel," has ",ngrid[ilevel-1]," grids.")
        if (ngrid[ilevel-1] == 0):
            print("Skipping loop, no more level...")
            break

        x = x[ind]
        y = y[ind]
        z = z[ind]
        q = [qq[ind] for qq in q]

        # read finer level data
        if (ilevel < levelmax):
            print("Reading data for finer level and computing KDTree")
            ds_fine = yt.load(filename,max_level=ilevel+1,max_level_convention="ramses")
            data_fine = ds_fine.all_data()
            dx_fine = data_fine["index","dx"].to("code_length").v
            x_fine = data_fine["index","x"].to("code_length").v
            y_fine = data_fine["index","y"].to("code_length").v
            z_fine = data_fine["index","z"].to("code_length").v
            ind_fine = np.where(dx_fine < 1.1*dxmin[ilevel+1])
            if len(ind_fine[0])>0:
                x_fine = x_fine[ind_fine]
                y_fine = y_fine[ind_fine]
                z_fine = z_fine[ind_fine]

                # compute refinement map
                xyz_fine = np.stack([x_fine, y_fine, z_fine], axis=-1)
                xyz_coarse = np.stack([x, y, z], axis=-1)
                tree_fine = KDTree(xyz_fine)
                distance, iii = tree_fine.query(xyz_coarse, p=np.inf, distance_upper_bound=dxmin[ilevel]/2, workers=-1)
                refined = np.isfinite(distance)
            else:
                refined = np.full(ngrid[ilevel-1]*2**ndim, False, dtype=bool)
        else:
            refined = np.full(ngrid[ilevel-1]*2**ndim, False, dtype=bool)

        print("Found ",len(np.where(refined == True)[0])," refined cells")

        # write number of grids in amr and hydro file
        f_amr.seek(12+4*(ilevel-levelmin))
        np.array([ngrid[ilevel-1]],dtype=np.int32).tofile(f_amr)
        if has_hydro:
            f_hydro.seek(16+4*(ilevel-levelmin))
            np.array([ngrid[ilevel-1]],dtype=np.int32).tofile(f_hydro)

        # compute cartesian keys sorted by hilbert index
        x = (x - dxmin[ilevel]/2)/dxmin[ilevel]
        y = (y - dxmin[ilevel]/2)/dxmin[ilevel]
        z = (z - dxmin[ilevel]/2)/dxmin[ilevel]
        x = np.round(x).astype(int)
        y = np.round(y).astype(int)
        z = np.round(z).astype(int)

        xg = x[::8]/2
        yg = y[::8]/2
        zg = z[::8]/2
        xg = xg.astype(int)
        yg = yg.astype(int)
        zg = zg.astype(int)

        print("Computing and sorting hilbert index")
        hkey = ram.hilbert3d(xg,yg,zg,ilevel-1)
        hkey = hkey.astype(int)
        ind = np.argsort(hkey)
        xg = xg[ind]
        yg = yg[ind]
        zg = zg[ind]
        print("Hilbert index done")

        # mapping C -> Fortran
        true_ind = [1,5,3,7,2,6,4,8]

        print("Writing amr file")

        # write amr data
        stride = ndim+1
        size = ngrid[ilevel-1]*stride
        out_array = np.empty(size,dtype=np.int32)
        out_array[0::stride] = xg
        out_array[1::stride] = yg
        out_array[2::stride] = zg
        refined_int = np.zeros(ngrid[ilevel-1],dtype=np.int32)
        for iind in range(0,2**ndim):
            jind = true_ind[iind]-1
            mask = refined[jind::2**ndim].astype(np.int32)
            refined_int = refined_int | (mask[ind] << iind)
        out_array[ndim::stride] = refined_int
        f_amr.seek(0,2)
        out_array.tofile(f_amr)

        if not has_hydro:
            continue

        print("Writing hydro file")

        # write hydro data
        size = ngrid[ilevel-1]*nvar*(2**ndim)
        out_array=np.empty(size,dtype=np.float32)
        for iind in range(0,8):
            jind = true_ind[iind]-1
            for ivar in range(0,nvar):
                out_1 = q[ivar][jind::8].astype(np.float32)
                out_array[ivar*8+iind::8*nvar] = out_1[ind]
        f_hydro.seek(0,2)
        out_array.tofile(f_hydro)

if has_hydro:
    f_hydro.close()
