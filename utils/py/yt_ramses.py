"""yt frontend for RAMSES output directories.

Usage:
    import yt
    import yt_ramses
    ds = yt.load("output_00012")

Importing this module registers RamsesDataset with yt.load. Data are read
through ramses_io.py. Only 3D outputs are supported for now.

Each amr.NNNNN file is a yt domain (domain_id = file number). Since the
domain decomposition is done separately at each level, the parent of an oct
or the child oct of a refined cell may be stored in another file. The oct
tree of each domain therefore also contains these octs, tagged with the
"foreign" domain id nfile+1, so that yt never counts their cells.

Authors: Romain Teyssier (Princeton University, October 2026)
"""
import os
import weakref
import numpy as np

from yt.data_objects.index_subobjects.octree_subset import OctreeSubset
from yt.data_objects.static_output import Dataset
from yt.fields.field_info_container import FieldInfoContainer
from yt.funcs import setdefaultattr
from yt.geometry.geometry_handler import YTDataChunk
from yt.geometry.oct_container import RAMSESOctreeContainer
from yt.geometry.oct_geometry_handler import OctreeIndex
from yt.utilities.io_handler import BaseIOHandler
from yt.utilities.physical_constants import kb, mp

from ramses_io import Output

# yt fluid type of each grid component
fluid_types = {"hydro": "ramses", "grav": "gravity", "rt": "rt", "cr": "cr"}

# yt name of particle fields, without the _x, _y, _z suffix for vector fields
particle_field_names = {
    "pos": "particle_position",
    "vel": "particle_velocity",
    "accel": "particle_acceleration",
    "angmom": "particle_angular_momentum",
    "mass": "particle_mass",
    "birth_id": "particle_identity",
    "level": "particle_level",
    "birth_date": "particle_birth_time",
    "merging_date": "particle_merging_time",
    "metallicity": "particle_metallicity",
}

def particle_field_list(fields):
    """Return the yt names of the particle fields of a header, with
    (field name, component) of the corresponding RAMSES field."""
    names = []
    for name, dtype, ncomp in fields:
        yname = particle_field_names.get(name,"particle_"+name)
        if ncomp==1:
            names.append((yname,name,0))
        else:
            for icomp in range(ncomp):
                names.append((yname+"_"+"xyz"[icomp],name,icomp))
    return names

def _rows(ckey):
    """View a (n, ndim) key array as a 1D array of rows, for set operations."""
    ckey = np.ascontiguousarray(ckey)
    return ckey.view(np.dtype((np.void,ckey.dtype.itemsize*ckey.shape[1]))).ravel()

def _unrows(rows, ckey_dtype, ndim):
    return rows.view(ckey_dtype).reshape(-1,ndim)

class RamsesDomainFile:
    """One amr.NNNNN file of the output and its oct tree."""

    _oct_handler_initialized = False

    def __init__(self, ds, domain_id):
        self.ds = ds
        self.domain_id = domain_id
        ifile = domain_id-1
        o = ds._output
        self.amr = o.grid["amr"][ifile]
        self.grid = {fluid_types[c]: o.grid[c][ifile] for c in fluid_types if c in o.grid}
        self.part = {p: o.part[p][ifile] for p in o.part}
        self.foreign_id = o.nfile+1

    def foreign_octs(self, own, refined):
        """Return, for each level, the keys of the octs that are missing from
        this file but are needed to build its tree: the ancestors of its octs
        and the children of its refined cells."""
        o = self.ds._output
        ndim = o.ndim
        levels = o.levels()
        ckey_dtype = own[levels[0]].dtype
        foreign = {ilevel: np.empty(0,dtype=_rows(np.zeros((1,ndim),ckey_dtype)).dtype)
                   for ilevel in levels}
        # Children of refined cells, cell ind-1 = i + 2*j + 4*k
        shift = np.array([[(ind>>idim)&1 for idim in range(ndim)] for ind in range(2**ndim)])
        for ilevel in levels[:-1]:
            ref = refined[ilevel]
            ioct, ind = np.nonzero(ref)
            if len(ioct)==0:
                continue
            kids = 2*own[ilevel][ioct].astype(np.int64)+shift[ind]
            kids = _rows(kids.astype(ckey_dtype))
            foreign[ilevel+1] = np.setdiff1d(kids,_rows(own[ilevel+1]))
        # Ancestors, from the finest level up
        for ilevel in levels[:0:-1]:
            keys = np.concatenate([_rows(own[ilevel]),foreign[ilevel]])
            if len(keys)==0:
                continue
            parents = _unrows(keys,ckey_dtype,ndim)//2
            missing = np.setdiff1d(_rows(parents),_rows(own[ilevel-1]))
            foreign[ilevel-1] = np.union1d(foreign[ilevel-1],missing)
        return {ilevel: _unrows(foreign[ilevel],ckey_dtype,ndim) for ilevel in levels}

    @property
    def oct_handler(self):
        if self._oct_handler_initialized:
            return self._oct_handler
        o = self.ds._output
        levels = o.levels()
        own = {}
        refined = {}
        for ilevel in levels:
            ckey, bits = self.amr.read_amr(ilevel)
            own[ilevel] = ckey
            refined[ilevel] = ((bits[:,None] >> np.arange(2**o.ndim)) & 1).astype(bool)
        foreign = self.foreign_octs(own,refined)

        counts = np.zeros(o.nfile+1,dtype=np.int64)
        counts[self.domain_id-1] = sum(len(own[ilevel]) for ilevel in levels)
        counts[self.foreign_id-1] = sum(len(foreign[ilevel]) for ilevel in levels)
        root_nodes = len(own[o.levelmin])+len(foreign[o.levelmin])

        oct_handler = RAMSESOctreeContainer(self.ds.domain_dimensions/2,
                                            self.ds.domain_left_edge,
                                            self.ds.domain_right_edge)
        oct_handler.allocate_domains(counts,root_nodes)
        # Add octs level by level, so that every parent already exists
        for ilevel in levels:
            for domain, ckey in [(self.foreign_id,foreign[ilevel]),(self.domain_id,own[ilevel])]:
                if len(ckey)==0:
                    continue
                pos = (ckey.astype(np.float64)+0.5)/2**(ilevel-1)
                oct_handler.add(domain,ilevel-o.levelmin,pos)
        oct_handler.finalize()
        self._oct_handler = oct_handler
        self._oct_handler_initialized = True
        return oct_handler

    def included(self, selector):
        domain_ids = self.oct_handler.domain_identify(selector)
        return self.domain_id in domain_ids

class RamsesDomainSubset(OctreeSubset):
    _domain_offset = 1
    _block_order = "F"

    @property
    def oct_handler(self):
        return self.domain.oct_handler

    def fill(self, ftype, fields, selector):
        """Return a dict field name -> cell values of the cells picked by selector."""
        oct_handler = self.oct_handler
        gfile = self.domain.grid[ftype]
        names = self.ds._var_names[ftype]
        cell_count = selector.count_oct_cells(oct_handler,self.domain_id)
        data = {f: np.zeros(cell_count,dtype="float64") for f in fields}
        if cell_count==0:
            return data
        levels, cell_inds, file_inds = oct_handler.file_index_octs(
            selector,self.domain_id,cell_count)
        for level in np.unique(levels):
            mask = levels==level
            block = gfile.read_level(int(level)+self.ds._output.levelmin)
            for f in fields:
                data[f][mask] = block[file_inds[mask],names.index(f),cell_inds[mask]]
        return data

class RamsesIndex(OctreeIndex):

    def __init__(self, ds, dataset_type="ramses_io"):
        self.dataset_type = dataset_type
        self.dataset = weakref.proxy(ds)
        self.index_filename = self.dataset.parameter_filename
        self.directory = os.path.dirname(self.index_filename)
        self.float_type = np.float64
        super().__init__(ds,dataset_type)

    def _initialize_oct_handler(self):
        self.domains = [RamsesDomainFile(self.dataset,ifile+1)
                        for ifile in range(self.dataset._output.nfile)]
        self.max_level = self.dataset.max_level

    def _detect_output_fields(self):
        o = self.dataset._output
        self.fluid_field_list = [(ftype,name) for ftype, names in self.dataset._var_names.items()
                                 for name in names]
        self.particle_field_list = [(ptype,yname) for ptype in o.part
                                    for yname, _, _ in particle_field_list(o.part_fields[ptype])]
        self.field_list = self.fluid_field_list+self.particle_field_list

    def _identify_base_chunk(self, dobj):
        if getattr(dobj,"_chunk_info",None) is None:
            domains = [dom for dom in self.domains if dom.included(dobj.selector)]
            base_region = getattr(dobj,"base_region",dobj)
            dobj._chunk_info = [RamsesDomainSubset(base_region,dom,self.dataset)
                                for dom in domains]
        dobj._current_chunk = list(self._chunk_all(dobj))[0]

    def _chunk_all(self, dobj):
        oobjs = getattr(dobj._current_chunk,"objs",dobj._chunk_info)
        yield YTDataChunk(dobj,"all",oobjs,None)

    def _chunk_spatial(self, dobj, ngz, sort=None, preload_fields=None):
        if ngz>0:
            raise NotImplementedError("Ghost zones are not supported yet")
        sobjs = getattr(dobj._current_chunk,"objs",dobj._chunk_info)
        for og in sobjs:
            yield YTDataChunk(dobj,"spatial",[og],None)

    def _chunk_io(self, dobj, cache=True, local_only=False):
        oobjs = getattr(dobj._current_chunk,"objs",dobj._chunk_info)
        for subset in oobjs:
            yield YTDataChunk(dobj,"io",[subset],None,cache=cache)

    def _get_particle_type_counts(self):
        o = self.dataset._output
        return {ptype: o.npart(ptype) for ptype in o.part}

rho_units = "code_density"
vel_units = "code_velocity"
pressure_units = "code_pressure"
b_units = "code_magnetic"
accel_units = "code_velocity / code_time"
specific_ener_units = "code_velocity**2"

class RamsesFieldInfo(FieldInfoContainer):
    known_other_fields = (
        ("density", (rho_units, ["density"], None)),
        ("velocity_x", (vel_units, ["velocity_x"], None)),
        ("velocity_y", (vel_units, ["velocity_y"], None)),
        ("velocity_z", (vel_units, ["velocity_z"], None)),
        ("thermal_pressure", (pressure_units, ["pressure"], None)),
        ("magnetic_field_x", (b_units, [], None)),
        ("magnetic_field_y", (b_units, [], None)),
        ("magnetic_field_z", (b_units, [], None)),
        ("metal_mass_fraction", ("", ["metallicity"], None)),
        ("xHI", ("", ["H_p0_fraction"], None)),
        ("xHII", ("", ["H_p1_fraction"], None)),
        ("xHeII", ("", ["He_p1_fraction"], None)),
        ("xHeIII", ("", ["He_p2_fraction"], None)),
        ("potential", (specific_ener_units, ["potential"], None)),
        ("accel_x", (accel_units, ["acceleration_x"], None)),
        ("accel_y", (accel_units, ["acceleration_y"], None)),
        ("accel_z", (accel_units, ["acceleration_z"], None)),
    )
    known_particle_fields = (
        ("particle_position_x", ("code_length", [], None)),
        ("particle_position_y", ("code_length", [], None)),
        ("particle_position_z", ("code_length", [], None)),
        ("particle_velocity_x", (vel_units, [], None)),
        ("particle_velocity_y", (vel_units, [], None)),
        ("particle_velocity_z", (vel_units, [], None)),
        ("particle_acceleration_x", (accel_units, [], None)),
        ("particle_acceleration_y", (accel_units, [], None)),
        ("particle_acceleration_z", (accel_units, [], None)),
        ("particle_mass", ("code_mass", [], None)),
        ("particle_identity", ("", ["particle_index"], None)),
        ("particle_level", ("", [], None)),
        ("particle_birth_time", ("code_time", [], None)),
        ("particle_merging_time", ("code_time", [], None)),
        ("particle_metallicity", ("", [], None)),
    )

class IOHandlerRamses(BaseIOHandler):
    _dataset_type = "ramses_io"

    def _read_fluid_selection(self, chunks, selector, fields, size):
        tr = {field: [] for field in fields}
        ftypes = {f[0] for f in fields}
        for chunk in chunks:
            for subset in chunk.objs:
                for ftype in ftypes:
                    names = [f for ft, f in fields if ft==ftype]
                    data = subset.fill(ftype,names,selector)
                    for f in names:
                        tr[ftype,f].append(data[f])
        return {field: np.concatenate(tr[field]) if tr[field] else np.empty(0,dtype="d")
                for field in fields}

    def _read_particle_subset(self, subset, ptype, ynames):
        """Read particle fields of one domain. Positions are divided by boxlen."""
        pfile = subset.domain.part[ptype]
        lookup = {yname: (name,icomp) for yname, name, icomp in
                  particle_field_list(subset.ds._output.part_fields[ptype])}
        tr = {}
        cache = {}
        for yname in ynames:
            name, icomp = lookup[yname]
            if name not in cache:
                cache[name] = pfile.read_field(name)
            data = cache[name]
            data = data if data.ndim==1 else data[icomp]
            if name=="pos":
                data = data/subset.ds.parameters["boxlen"]
            tr[yname] = data
        return tr

    def _read_particle_coords(self, chunks, ptf):
        pn = ["particle_position_"+ax for ax in "xyz"]
        for chunk in chunks:
            for subset in chunk.objs:
                for ptype in sorted(ptf):
                    rv = self._read_particle_subset(subset,ptype,pn)
                    yield ptype, tuple(np.asarray(rv[p],"=f8") for p in pn), 0.0

    def _read_particle_fields(self, chunks, ptf, selector):
        pn = ["particle_position_"+ax for ax in "xyz"]
        for chunk in chunks:
            for subset in chunk.objs:
                for ptype, field_list in sorted(ptf.items()):
                    rv = self._read_particle_subset(subset,ptype,set(field_list)|set(pn))
                    x, y, z = (np.asarray(rv[p],"=f8") for p in pn)
                    mask = selector.select_points(x,y,z,0.0)
                    if mask is None:
                        continue
                    for field in field_list:
                        yield (ptype,field), np.asarray(rv[field][mask],"=f8")

class RamsesDataset(Dataset):
    _index_class = RamsesIndex
    _field_info_class = RamsesFieldInfo

    def __init__(self, filename, dataset_type="ramses_io", units_override=None,
                 unit_system="cgs", default_species_fields=None):
        if os.path.isdir(filename):
            filename = os.path.join(filename,"info.txt")
        self.root_folder = os.path.dirname(os.path.abspath(filename))
        self._output = Output(self.root_folder)
        if self._output.ndim!=3:
            raise NotImplementedError("yt_ramses only supports 3D outputs")
        self._var_names = {fluid_types[c]: self._output.var_names[c]
                           for c in fluid_types if c in self._output.grid}
        super().__init__(filename,dataset_type,units_override=units_override,
                         unit_system=unit_system,default_species_fields=default_species_fields)
        self.fluid_types += tuple(self._var_names)
        self.particle_types = self.particle_types_raw = tuple(self._output.part)

    def __str__(self):
        return os.path.basename(self.root_folder)

    def _set_code_unit_attributes(self):
        boxlen = self.parameters["boxlen"]
        length_unit = self.parameters["unit_l"]
        density_unit = self.parameters["unit_d"]
        time_unit = self.parameters["unit_t"]
        mass_unit = density_unit*length_unit**3
        magnetic_unit = np.sqrt(4*np.pi*mass_unit/(time_unit**2*length_unit))
        pressure_unit = density_unit*(length_unit/time_unit)**2
        setdefaultattr(self,"density_unit",self.quan(density_unit,"g/cm**3"))
        setdefaultattr(self,"magnetic_unit",self.quan(magnetic_unit,"gauss"))
        setdefaultattr(self,"pressure_unit",self.quan(pressure_unit,"dyne/cm**2"))
        setdefaultattr(self,"time_unit",self.quan(time_unit,"s"))
        setdefaultattr(self,"mass_unit",self.quan(mass_unit,"g"))
        setdefaultattr(self,"velocity_unit",self.quan(length_unit,"cm")/self.time_unit)
        temperature_unit = self.velocity_unit**2*mp/kb
        setdefaultattr(self,"temperature_unit",temperature_unit.in_units("K"))
        # Only the length unit is scaled by boxlen, positions are in [0,1]
        setdefaultattr(self,"length_unit",self.quan(length_unit*boxlen,"cm"))

    def _parse_parameter_file(self):
        o = self._output
        self.parameters.update(o.info)
        self.unique_identifier = str(int(os.stat(self.parameter_filename).st_ctime))
        self.dimensionality = 3
        self.refine_by = 2
        self.min_level = o.levelmin-1
        self.max_level = o.nlevelmax-o.levelmin
        self.domain_left_edge = np.zeros(3,dtype="float64")
        self.domain_right_edge = np.ones(3,dtype="float64")
        self.domain_dimensions = np.ones(3,dtype="int32")*2**o.levelmin
        self._periodicity = (True,True,True)
        self.gamma = o.info["gamma"]
        self.current_time = o.info["time"]
        if o.info["aexp"]!=1 or o.info["H0"]!=1:
            self.cosmological_simulation = True
            self.current_redshift = 1/o.info["aexp"]-1
            self.omega_lambda = o.info["omega_l"]
            self.omega_matter = o.info["omega_m"]
            self.hubble_constant = o.info["H0"]/100
        else:
            self.cosmological_simulation = False
            self.current_redshift = 0
            self.omega_lambda = 0
            self.omega_matter = 0
            self.hubble_constant = 0

    @classmethod
    def _is_valid(cls, filename, *args, **kwargs):
        if os.path.isdir(filename):
            folder = filename
        elif os.path.basename(filename)=="info.txt":
            folder = os.path.dirname(filename)
        else:
            return False
        return (os.path.exists(os.path.join(folder,"info.txt"))
                and os.path.exists(os.path.join(folder,"amr.00001")))
