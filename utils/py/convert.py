#!/usr/bin/env python3
import yt
import numpy as np
import ramses as ram
import argparse
import os
from scipy.io import FortranFile
from scipy.spatial import KDTree

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
def wr_header(fileloc,npart,nfile,star):
    with open(fileloc, 'w') as f:
        # Write header information
        f.write("Total number of particles\n")
        f.write(f"{npart}\n")
        f.write("Total number of files\n")
        f.write(f"{nfile}\n")
        # Particle fields: name, numpy data type, number of components
        f.write("Particle fields\n")
        f.write("pos f4 3\n")
        f.write("vel f4 3\n")
        f.write("mass f4 1\n")
        if star:
            f.write("metallicity f4 1\n")
            f.write("birth_date f4 1\n")
        f.write("level i4 1\n")
        f.write("birth_id i8 1\n")

#================================
# write new ramses parameter file
#================================
def wr_params(params1,params2,nout,**kwargs):
    path = kwargs.get("path","./")
    car1 = str(nout).zfill(5)
    filename = path+"/output_"+car1+"/params.bin"

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
        dble_array = np.array([params2["aexp"],params2["hexp"],params2["aexp_old"],
                               params2["aexp_old"],params2["epot_tot_int"],params2["epot_tot_old"],params2["mass_sph"],1.6667],dtype=np.float64)
        dble_array.tofile(f)
        int_array = np.array([1,2,3,4,5,6],dtype=np.int64)
        int_array.tofile(f)

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
args = parser.parse_args()
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
fields = [field[1] for field in ds.field_list if field[0] == 'ramses']
print("levelmin=",levelmin)
print("levelmax=",levelmax)
print("fields=",fields)
nvar = 5
if("Metallicity" in fields):
    nvar = nvar + 1
if("hydro_scalar_01" in fields):
    nvar = nvar + 1
if("hydro_scalar_02" in fields):
    nvar = nvar + 1
if("hydro_scalar_03" in fields):
    nvar = nvar + 1

if not os.path.isdir(path_out):
    os.makedirs(path_out)
os.makedirs(path_out+"/output_"+car1,exist_ok=True)

wr_params(ds.parameters,params2,nout,path=path_out)
print("Parameter file written.")

#================================
# read and write particle files
#================================

data = ds.all_data()

print("Loading and writing DM particle data")

file_part = path_out+"/output_"+car1+"/part.00001"

x = data["DM","particle_position_x"].to("code_length").v.astype(np.float32)
y = data["DM","particle_position_y"].to("code_length").v.astype(np.float32)
z = data["DM","particle_position_z"].to("code_length").v.astype(np.float32)
vx = data["DM","particle_velocity_x"].to("code_velocity").v.astype(np.float32)
vy = data["DM","particle_velocity_y"].to("code_velocity").v.astype(np.float32)
vz = data["DM","particle_velocity_z"].to("code_velocity").v.astype(np.float32)
m = data["DM","particle_mass"].to("code_mass").v.astype(np.float32)
level = data["DM","particle_level"].v.astype(np.int32)
idp = data["DM","particle_identity"].v.astype(np.int64)

npart = len(x)

header_part = path_out+"/output_"+car1+"/part_header.txt"

wr_header(header_part,npart,1,False)

with open(file_part, "wb") as f_part:
    np.array([ndim],dtype=np.int32).tofile(f_part)
    np.array([npart],dtype=np.int32).tofile(f_part)
    x.tofile(f_part)
    y.tofile(f_part)
    z.tofile(f_part)
    vx.tofile(f_part)
    vy.tofile(f_part)
    vz.tofile(f_part)
    m.tofile(f_part)
    level.tofile(f_part)
    idp.tofile(f_part)

if (dm_only):
    exit()

print("Loading and writing star particle data")

file_star = path_out+"/output_"+car1+"/star.00001"

x = data["star","particle_position_x"].to("code_length").v.astype(np.float32)
y = data["star","particle_position_y"].to("code_length").v.astype(np.float32)
z = data["star","particle_position_z"].to("code_length").v.astype(np.float32)
vx = data["star","particle_velocity_x"].to("code_velocity").v.astype(np.float32)
vy = data["star","particle_velocity_y"].to("code_velocity").v.astype(np.float32)
vz = data["star","particle_velocity_z"].to("code_velocity").v.astype(np.float32)
m = data["star","particle_mass"].to("code_mass").v.astype(np.float32)
level = data["star","particle_level"].v.astype(np.int32)
idp = data["star","particle_identity"].v.astype(np.int64)
metal = data["star","particle_metallicity"].v.astype(np.float32)
birth = data["star","particle_birth_time"].v.astype(np.float32)

npart = len(x)

header_star = path_out+"/output_"+car1+"/star_header.txt"

wr_header(header_star,npart,1,True)

with open(file_star, "wb") as f_star:
    np.array([ndim],dtype=np.int32).tofile(f_star)
    np.array([npart],dtype=np.int32).tofile(f_star)
    x.tofile(f_star)
    y.tofile(f_star)
    z.tofile(f_star)
    vx.tofile(f_star)
    vy.tofile(f_star)
    vz.tofile(f_star)
    m.tofile(f_star)
    metal.tofile(f_star)
    birth.tofile(f_star)
    level.tofile(f_star)
    idp.tofile(f_star)

#===================================
# read and write amr and hydro file
#===================================

file_amr = path_out+"/output_"+car1+"/amr.00001"
file_hydro = path_out+"/output_"+car1+"/hydro.00001"

ngrid = np.zeros(levelmax,dtype=np.int32)
ngridtot = np.sum(ngrid)

with open(file_amr, "wb") as f_amr, open(file_hydro, "wb") as f_hydro:
    np.array([ndim],dtype=np.int32).tofile(f_amr)
    np.array([levelmin],dtype=np.int32).tofile(f_amr)
    np.array([levelmax],dtype=np.int32).tofile(f_amr)
    for ilevel in range(levelmin,levelmax+1):
        np.array([ngrid[ilevel-1]],dtype=np.int32).tofile(f_amr)
    np.array([ndim],dtype=np.int32).tofile(f_hydro)
    np.array([nvar],dtype=np.int32).tofile(f_hydro)
    np.array([levelmin],dtype=np.int32).tofile(f_hydro)
    np.array([levelmax],dtype=np.int32).tofile(f_hydro)
    for ilevel in range(levelmin,levelmax+1):
        np.array([ngrid[ilevel-1]],dtype=np.int32).tofile(f_amr)

with open(file_amr, "wb") as f_amr, open(file_hydro, "wb") as f_hydro:

    dxmin = 1/2**(np.linspace(0,levelmax,levelmax+1))
    for ilevel in range(levelmin,levelmax+1):
        print("Reading data for level=",ilevel)

        # read current level data
        ds = yt.load(filename,max_level=ilevel,max_level_convention="ramses")
        data = ds.all_data()
        data.get_data([_ for _ in ds.field_list if _[0] == "ramses"])
        dx = data["ramses","dx"].to("code_length").v
        x = data["ramses","x"].to("code_length").v
        y = data["ramses","y"].to("code_length").v
        z = data["ramses","z"].to("code_length").v
        d = data["ramses","Density"].v
        p = data["ramses","Pressure"].v
        vx = data["ramses","x-velocity"].v
        vy = data["ramses","y-velocity"].v
        vz = data["ramses","z-velocity"].v
        if("Metallicity" in fields):
            metal = data["ramses","Metallicity"].v
        if("hydro_scalar_01" in fields):
            pscal1 = data["ramses","hydro_scalar_01"].v
        if("hydro_scalar_02" in fields):
            pscal2 = data["ramses","hydro_scalar_02"].v
        if("hydro_scalar_03" in fields):
            pscal3 = data["ramses","hydro_scalar_03"].v
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
        d = d[ind]
        p = p[ind]
        vx = vx[ind]
        vy = vy[ind]
        vz = vz[ind]
        if("Metallicity" in fields):
            metal = metal[ind]
        if("hydro_scalar_01" in fields):
            pscal1 = pscal1[ind]
        if("hydro_scalar_02" in fields):
            pscal2 = pscal2[ind]
        if("hydro_scalar_03" in fields):
            pscal3 = pscal3[ind]

        # read finer level data
        if (ilevel < levelmax):
            print("Reading data for finer level and computing KDTree")
            ds_fine = yt.load(filename,max_level=ilevel+1,max_level_convention="ramses")
            data_fine = ds_fine.all_data()
            # speed up file read
            data_fine.get_data([_ for _ in ds.field_list if _[0] == "ramses"])
            dx_fine = data_fine["ramses","dx"].to("code_length").v
            x_fine = data_fine["ramses","x"].to("code_length").v
            y_fine = data_fine["ramses","y"].to("code_length").v
            z_fine = data_fine["ramses","z"].to("code_length").v
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
        f_amr.seek(12+4*(ilevel-1))
        np.array([ngrid[ilevel-1]],dtype=np.int32).tofile(f_amr)
        f_hydro.seek(16+4*(ilevel-1))
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

        print("Writing hydro file")

        # write hydro data
        size = ngrid[ilevel-1]*nvar*(2**ndim)
        out_array=np.empty(size,dtype=np.float32)
        for iind in range(0,8):
            jind = true_ind[iind]-1
            out_1 = d[jind::8].astype(np.float32)
            out_array[0*8+iind::8*nvar] = out_1[ind]
            out_1 = vx[jind::8].astype(np.float32)
            out_array[1*8+iind::8*nvar] = out_1[ind]
            out_1 = vy[jind::8].astype(np.float32)
            out_array[2*8+iind::8*nvar] = out_1[ind]
            out_1 = vz[jind::8].astype(np.float32)
            out_array[3*8+iind::8*nvar] = out_1[ind]
            out_1 = p[jind::8].astype(np.float32)
            out_array[4*8+iind::8*nvar] = out_1[ind]
            if("Metallicity" in fields):
                out_1 = metal[jind::8].astype(np.float32)
                out_array[5*8+iind::8*nvar] = out_1[ind]
            if("hydro_scalar_01" in fields):
                out_1 = pscal1[jind::8].astype(np.float32)
                out_array[6*8+iind::8*nvar] = out_1[ind]
            if("hydro_scalar_02" in fields):
                out_1 = pscal2[jind::8].astype(np.float32)
                out_array[7*8+iind::8*nvar] = out_1[ind]
            if("hydro_scalar_03" in fields):
                out_1 = pscal3[jind::8].astype(np.float32)
                out_array[8*8+iind::8*nvar] = out_1[ind]
        f_hydro.seek(0,2)
        out_array.tofile(f_hydro)
