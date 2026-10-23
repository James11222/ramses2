"""Low-level access to RAMSES output files.

This module parses the headers of an output directory output_NNNNN/ and
gives per-file, per-level access to the binary data, without loading
anything until asked. It is shared by ramses2hdf5.py and other readers.

File layouts (stream access, native endianness):
    amr.NNNNN       int32 ndim, levelmin, nlevelmax, noct(levelmin:nlevelmax)
                    then for each level and each oct: int32 ckey(ndim), int32 refined
    hydro.NNNNN     int32 ndim, nvar, levelmin, nlevelmax, noct(levelmin:nlevelmax)
    (grav, rt, cr)  then for each level: float32 (noct, nvar, 2**ndim)
    <prefix>.NNNNN  int32 ndim, npart, then each field column in header order
                    (prefix = part, star, sink, tree, trac, dust, core)

Authors: Romain Teyssier (Princeton University, October 2026)
"""
import numpy as np
import os

grid_components = ["amr","hydro","grav","rt","cr"]
part_components = ["part","star","sink","tree","trac","dust","core"]

# Data type of the oct Cartesian keys in amr files
ckey_dtype = np.dtype(np.int32)

def rd_info_txt(filename):
    """Parse a "key = value" text file (info.txt, rt_info.txt, cr_info.txt)
    into a dict. Values are converted to int, float, or a float array when
    several numbers are given; anything else is kept as a string."""
    info = {}
    with open(filename) as f:
        for line in f:
            if "=" not in line:
                continue
            key, value = line.split("=",1)
            key = key.strip()
            tokens = value.split()
            try:
                if len(tokens)==1:
                    try:
                        info[key] = int(tokens[0])
                    except ValueError:
                        info[key] = float(tokens[0])
                else:
                    info[key] = np.array([float(t) for t in tokens])
            except ValueError:
                info[key] = value.strip()
    return info

def rd_var_names(filename):
    """Return the variable names listed in a grid header file
    (hydro_header.txt, grav_header.txt, rt_header.txt or cr_header.txt),
    in file order."""
    with open(filename) as f:
        lines = [line.strip() for line in f]
    names = []
    for line in lines:
        if line.startswith("variable #"):
            names.append(line.split(":",1)[1].strip())
    return names

def rd_part_fields(filename):
    """Return the particle fields listed in a <prefix>_header.txt file,
    in file order, as a list of (name, numpy dtype, number of components)."""
    with open(filename) as f:
        lines = [line.strip() for line in f]
    fields = []
    for line in lines[lines.index("Particle fields")+1:]:
        if line=="":
            continue
        name, dtype, ncomp = line.split()
        fields.append((name,np.dtype(dtype),int(ncomp)))
    return fields

def unpack_refined(bits, ndim):
    """Convert the refined bit mask of each oct into a (noct, 2**ndim) bool array."""
    ind = np.arange(2**ndim)
    return ((bits[:,None] >> ind) & 1).astype(bool)

class GridFile:
    """One amr, hydro, grav, rt or cr file: header and byte offset of every level block."""

    def __init__(self, filename, component):
        self.filename = filename
        self.component = component
        head = np.fromfile(filename,dtype=np.int32,count=4)
        self.ndim = int(head[0])
        if component=="amr":
            self.nvar = 0
            self.levelmin = int(head[1])
            self.nlevelmax = int(head[2])
            hsize = 12
        else:
            self.nvar = int(head[1])
            self.levelmin = int(head[2])
            self.nlevelmax = int(head[3])
            hsize = 16
        nlevel = self.nlevelmax-self.levelmin+1
        noct = np.fromfile(filename,dtype=np.int32,count=nlevel,offset=hsize)
        self.noct = {}
        self.offset = {}
        offset = np.int64(hsize+4*nlevel)
        for ilevel in range(self.levelmin,self.nlevelmax+1):
            n = int(noct[ilevel-self.levelmin])
            self.noct[ilevel] = n
            self.offset[ilevel] = offset
            offset = offset + n*self.record_size()

    def record_size(self):
        """Number of bytes per oct."""
        if self.component=="amr":
            return self.ndim*ckey_dtype.itemsize+4
        else:
            return 4*self.nvar*2**self.ndim

    def read_amr(self, ilevel):
        """Return ckey (noct, ndim) and the refined bit mask (noct,) of a level."""
        n = self.noct[ilevel]
        rec = np.dtype([("ckey",ckey_dtype,(self.ndim,)),("refined",np.int32)])
        data = np.fromfile(self.filename,dtype=rec,count=n,offset=self.offset[ilevel])
        return data["ckey"], data["refined"]

    def read_level(self, ilevel):
        """Return the cell data of a level as a (noct, nvar, 2**ndim) float32 array."""
        n = self.noct[ilevel]
        if n==0:
            return np.zeros((0,self.nvar,2**self.ndim),dtype=np.float32)
        return np.memmap(self.filename,dtype=np.float32,mode="r",
                         offset=self.offset[ilevel],shape=(n,self.nvar,2**self.ndim))

class PartFile:
    """One particle file: number of particles and byte offset of every field column."""

    def __init__(self, filename, fields):
        self.filename = filename
        self.fields = fields
        head = np.fromfile(filename,dtype=np.int32,count=2)
        self.ndim = int(head[0])
        self.npart = int(head[1])
        self.offset = {}
        offset = np.int64(8)
        for name, dtype, ncomp in fields:
            self.offset[name] = offset
            offset = offset + ncomp*self.npart*dtype.itemsize

    def read_field(self, name):
        """Return a field as a (npart,) array, or (ncomp, npart) for vector fields."""
        for fname, dtype, ncomp in self.fields:
            if fname==name:
                break
        else:
            raise KeyError(name)
        data = np.fromfile(self.filename,dtype=dtype,count=ncomp*self.npart,
                           offset=self.offset[name])
        if ncomp==1:
            return data
        return data.reshape(ncomp,self.npart)

class Output:
    """A whole output directory: info, headers and the lists of data files.

    Attributes:
        info:        dict from info.txt
        ndim, nfile, levelmin, nlevelmax
        grid:        dict component -> list of GridFile (only components present)
        var_names:   dict component -> list of variable names (not for amr)
        part:        dict prefix -> list of PartFile (only particle types present)
        part_fields: dict prefix -> list of (name, dtype, ncomp)
        rt_info, cr_info: dicts from rt_info.txt and cr_info.txt when present
    """

    def __init__(self, path):
        self.path = path
        self.info = rd_info_txt(os.path.join(path,"info.txt"))
        self.ndim = self.info["ndim"]
        self.nfile = self.info["nfile"]
        self.levelmin = self.info["levelmin"]
        self.nlevelmax = self.info["levelmax"]
        self.rt_info = self._rd_optional_info("rt_info.txt")
        self.cr_info = self._rd_optional_info("cr_info.txt")

        self.grid = {}
        self.var_names = {}
        for comp in grid_components:
            if not os.path.exists(self._filename(comp,1)):
                continue
            files = [GridFile(self._filename(comp,ifile),comp) for ifile in range(1,self.nfile+1)]
            for g in files:
                if (g.ndim,g.levelmin,g.nlevelmax)!=(self.ndim,self.levelmin,self.nlevelmax):
                    raise ValueError("Header of "+g.filename+" does not match info.txt")
            self.grid[comp] = files
            if comp!="amr":
                names = rd_var_names(os.path.join(path,comp+"_header.txt"))
                if len(names)!=files[0].nvar:
                    raise ValueError("Found "+str(len(names))+" names in "+comp+"_header.txt"
                                     +" but nvar="+str(files[0].nvar))
                self.var_names[comp] = names

        self.part = {}
        self.part_fields = {}
        for prefix in part_components:
            header = os.path.join(path,prefix+"_header.txt")
            if not os.path.exists(header):
                continue
            fields = rd_part_fields(header)
            self.part_fields[prefix] = fields
            self.part[prefix] = [PartFile(self._filename(prefix,ifile),fields)
                                 for ifile in range(1,self.nfile+1)]

    def _filename(self, prefix, ifile):
        return os.path.join(self.path,prefix+"."+str(ifile).zfill(5))

    def _rd_optional_info(self, name):
        filename = os.path.join(self.path,name)
        if os.path.exists(filename):
            return rd_info_txt(filename)
        return None

    def levels(self):
        """List of the levels stored in the files."""
        return list(range(self.levelmin,self.nlevelmax+1))

    def noct_per_level(self, comp="amr"):
        """Dict level -> total number of octs over all files."""
        return {ilevel: sum(g.noct[ilevel] for g in self.grid[comp]) for ilevel in self.levels()}

    def npart(self, prefix):
        """Total number of particles of a given type over all files."""
        return sum(p.npart for p in self.part[prefix])
