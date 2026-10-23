#!/usr/bin/env python3
"""Convert a RAMSES output directory output_NNNNN/ to a single HDF5 file.

Usage:
    ramses2hdf5.py output_00012 [-o output_00012.h5] [--compression gzip|zstd|none] [--cells]

HDF5 layout (format_version 1). By default the grid is stored oct by oct
(oct layout). Octs are stored level by level, and within a level file by
file, so that each file block is a contiguous spatial domain.
    /                        attrs: format, format_version and all info.txt entries
    /namelist, /compilation  text copies of namelist.txt and compilation.txt
    /amr/level_offsets       int64 (nlevel+1): oct range of each level
    /amr/file_offsets        int64 (nlevel, nfile+1): oct range of each file within a level
    /amr/ckey                int32 (noct, ndim): Cartesian index of the oct at its level
    /amr/refined             uint8 (noct): refinement bit mask, bit ind-1 for cell ind
    /<comp>/<var>            float32 (noct, 2**ndim) for comp = hydro, grav, rt, cr
    /particles/<prefix>/file_offsets   int64 (nfile+1)
    /particles/<prefix>/<field>        (npart,) or (npart, ncomp)
Cell ind within an oct runs with x fastest, then y, then z. The center of cell
ind of an oct at level ilevel is (2*ckey+1+offset(ind))*dx with
dx = boxlen/2**ilevel and offset = -1/2 or +1/2.
Refined cells are kept: their values are averages of their children, so sums
over all cells count refined regions twice. Use the refined bits to keep only
leaf cells, or ramses_io.leaf_cells().

With --cells, only leaf cells are stored (cell layout), ordered as in the
oct layout. /amr is replaced by /cells and the grid fields are 1D:
    /cells/level_offsets     int64 (nlevel+1): cell range of each level
    /cells/file_offsets      int64 (nlevel, nfile+1): cell range of each file within a level
    /cells/level             int8 (ncell): cell level
    /cells/icell             int32 (ncell, ndim): Cartesian index of the cell at its level
    /<comp>/<var>            float32 (ncell,)
The center of a cell is (icell+0.5)*dx with dx = boxlen/2**level.

Authors: Romain Teyssier (Princeton University, October 2026)
"""
import argparse
import os
import numpy as np
import h5py
from ramses_io import Output, cell_index, unpack_refined

format_version = 1

# Target chunk size in bytes
chunk_bytes = 2**20

def compression_opts(compression):
    """Return the h5py dataset keywords for a compression choice."""
    if compression=="gzip":
        return dict(compression="gzip",compression_opts=4,shuffle=True)
    if compression=="zstd":
        import hdf5plugin
        return dict(**hdf5plugin.Zstd(),shuffle=True)
    if compression=="none":
        return {}
    raise ValueError("Unknown compression "+compression)

def create_dataset(group, name, shape, dtype, opts):
    """Create a dataset chunked along its first axis, or contiguous if empty."""
    dtype = np.dtype(dtype)
    if shape[0]==0:
        return group.create_dataset(name,shape=shape,dtype=dtype)
    rowsize = dtype.itemsize*int(np.prod(shape[1:]))
    nrow = max(1,min(shape[0],chunk_bytes//rowsize))
    return group.create_dataset(name,shape=shape,dtype=dtype,
                                chunks=(nrow,)+tuple(shape[1:]),**opts)

def write_text(f, o, name):
    filename = os.path.join(o.path,name+".txt")
    if os.path.exists(filename):
        with open(filename) as t:
            f.create_dataset(name,data=t.read())

def write_attrs(f, o):
    f.attrs["format"] = "ramses"
    f.attrs["format_version"] = format_version
    for key, value in o.info.items():
        f.attrs[key] = value
    write_text(f,o,"namelist")
    write_text(f,o,"compilation")

def oct_offsets(o, comp):
    """Return level_offsets (nlevel+1) and file_offsets (nlevel, nfile+1)."""
    levels = o.levels()
    files = o.grid[comp]
    level_offsets = np.zeros(len(levels)+1,dtype=np.int64)
    file_offsets = np.zeros((len(levels),o.nfile+1),dtype=np.int64)
    for i, ilevel in enumerate(levels):
        file_offsets[i,0] = level_offsets[i]
        for ifile, g in enumerate(files):
            file_offsets[i,ifile+1] = file_offsets[i,ifile]+g.noct[ilevel]
        level_offsets[i+1] = file_offsets[i,-1]
    return level_offsets, file_offsets

def write_amr(f, o, opts):
    level_offsets, file_offsets = oct_offsets(o,"amr")
    noct = int(level_offsets[-1])
    grp = f.create_group("amr")
    grp.attrs["levelmin"] = o.levelmin
    grp.attrs["nlevelmax"] = o.nlevelmax
    grp.create_dataset("level_offsets",data=level_offsets)
    grp.create_dataset("file_offsets",data=file_offsets)
    ckey = create_dataset(grp,"ckey",(noct,o.ndim),np.int32,opts)
    refined = create_dataset(grp,"refined",(noct,),np.uint8,opts)
    for i, ilevel in enumerate(o.levels()):
        for ifile, g in enumerate(o.grid["amr"]):
            n = g.noct[ilevel]
            if n==0:
                continue
            start = file_offsets[i,ifile]
            k, r = g.read_amr(ilevel)
            ckey[start:start+n] = k
            refined[start:start+n] = r.astype(np.uint8)
    return level_offsets

def write_grid(f, o, comp, opts, amr_level_offsets):
    level_offsets, file_offsets = oct_offsets(o,comp)
    if not np.array_equal(level_offsets,amr_level_offsets):
        raise ValueError("Oct counts in "+comp+" files do not match amr files")
    noct = int(level_offsets[-1])
    ncell = 2**o.ndim
    grp = f.create_group(comp)
    names = o.var_names[comp]
    for ivar, name in enumerate(names):
        dset = create_dataset(grp,name,(noct,ncell),np.float32,opts)
        dset.attrs["index"] = ivar+1
    if comp=="rt" and o.rt_info is not None:
        for key, value in o.rt_info.items():
            grp.attrs[key] = value
    if comp=="cr" and o.cr_info is not None:
        for key, value in o.cr_info.items():
            grp.attrs[key] = value
    for i, ilevel in enumerate(o.levels()):
        for ifile, g in enumerate(o.grid[comp]):
            n = g.noct[ilevel]
            if n==0:
                continue
            start = file_offsets[i,ifile]
            block = g.read_level(ilevel)
            for ivar, name in enumerate(names):
                grp[name][start:start+n] = block[:,ivar,:]

def write_cells(f, o, opts):
    """Write the leaf cells (cell layout) of all grid components."""
    levels = o.levels()
    # Count leaf cells per level and file
    nleaf = np.zeros((len(levels),o.nfile),dtype=np.int64)
    for i, ilevel in enumerate(levels):
        for ifile, g in enumerate(o.grid["amr"]):
            _, bits = g.read_amr(ilevel)
            nleaf[i,ifile] = np.count_nonzero(~unpack_refined(bits,o.ndim))
    level_offsets = np.zeros(len(levels)+1,dtype=np.int64)
    file_offsets = np.zeros((len(levels),o.nfile+1),dtype=np.int64)
    for i in range(len(levels)):
        file_offsets[i,0] = level_offsets[i]
        file_offsets[i,1:] = level_offsets[i]+np.cumsum(nleaf[i])
        level_offsets[i+1] = file_offsets[i,-1]
    ncell = int(level_offsets[-1])
    grp = f.create_group("cells")
    grp.attrs["levelmin"] = o.levelmin
    grp.attrs["nlevelmax"] = o.nlevelmax
    grp.create_dataset("level_offsets",data=level_offsets)
    grp.create_dataset("file_offsets",data=file_offsets)
    level = create_dataset(grp,"level",(ncell,),np.int8,opts)
    icell = create_dataset(grp,"icell",(ncell,o.ndim),np.int32,opts)
    comps = [comp for comp in ["hydro","grav","rt","cr"] if comp in o.grid]
    for comp in comps:
        if not np.array_equal(oct_offsets(o,comp)[0],oct_offsets(o,"amr")[0]):
            raise ValueError("Oct counts in "+comp+" files do not match amr files")
        grp = f.create_group(comp)
        for ivar, name in enumerate(o.var_names[comp]):
            dset = create_dataset(grp,name,(ncell,),np.float32,opts)
            dset.attrs["index"] = ivar+1
        if comp=="rt" and o.rt_info is not None:
            for key, value in o.rt_info.items():
                grp.attrs[key] = value
        if comp=="cr" and o.cr_info is not None:
            for key, value in o.cr_info.items():
                grp.attrs[key] = value
    for i, ilevel in enumerate(levels):
        for ifile, g in enumerate(o.grid["amr"]):
            n = nleaf[i,ifile]
            if n==0:
                continue
            start = file_offsets[i,ifile]
            ckey, bits = g.read_amr(ilevel)
            leaf = ~unpack_refined(bits,o.ndim)
            level[start:start+n] = ilevel
            icell[start:start+n] = cell_index(ckey,o.ndim)[leaf]
            for comp in comps:
                block = o.grid[comp][ifile].read_level(ilevel)
                for ivar, name in enumerate(o.var_names[comp]):
                    f[comp][name][start:start+n] = block[:,ivar,:][leaf]

def write_part(f, o, prefix, opts):
    files = o.part[prefix]
    file_offsets = np.zeros(o.nfile+1,dtype=np.int64)
    for ifile, p in enumerate(files):
        file_offsets[ifile+1] = file_offsets[ifile]+p.npart
    npart = int(file_offsets[-1])
    grp = f.require_group("particles").create_group(prefix)
    grp.create_dataset("file_offsets",data=file_offsets)
    for name, dtype, ncomp in o.part_fields[prefix]:
        shape = (npart,) if ncomp==1 else (npart,ncomp)
        dset = create_dataset(grp,name,shape,dtype,opts)
        for ifile, p in enumerate(files):
            if p.npart==0:
                continue
            start = file_offsets[ifile]
            data = p.read_field(name)
            dset[start:start+p.npart] = data if ncomp==1 else data.T

def convert(path, h5name, compression="gzip", cells=False):
    """Convert the output directory path into the HDF5 file h5name,
    in oct layout, or in cell layout if cells is True."""
    o = Output(path)
    if "amr" not in o.grid:
        raise ValueError("No amr files found in "+path)
    opts = compression_opts(compression)
    with h5py.File(h5name,"w") as f:
        write_attrs(f,o)
        if cells:
            write_cells(f,o,opts)
        else:
            amr_level_offsets = write_amr(f,o,opts)
            for comp in ["hydro","grav","rt","cr"]:
                if comp in o.grid:
                    write_grid(f,o,comp,opts,amr_level_offsets)
        for prefix in o.part:
            write_part(f,o,prefix,opts)

if __name__=="__main__":
    parser = argparse.ArgumentParser(description="Convert a RAMSES output directory to HDF5")
    parser.add_argument("path", help="output directory, e.g. output_00012")
    parser.add_argument("-o", "--out", help="HDF5 file name (default: <path>.h5)")
    parser.add_argument("--compression", default="gzip", choices=["gzip","zstd","none"],
                        help="compression filter (default: gzip)")
    parser.add_argument("--cells", action="store_true",
                        help="store only leaf cells (cell layout) instead of octs")
    args = parser.parse_args()
    path = os.path.normpath(args.path)
    h5name = args.out if args.out else path+".h5"
    convert(path,h5name,args.compression,args.cells)
    print("Wrote "+h5name)
