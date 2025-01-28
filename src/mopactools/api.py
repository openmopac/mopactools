# Python wrappers for the MOPAC API
import os
import copy
import math
import ctypes
import numpy as np
import scipy as sp
from . import binding

class mopac_system:
    model_dict = { "PM7":0, "PM6-D3H4":1, "PM6-ORG":2, "PM6":3, "AM1":4, "RM1":5 }
    def __init__(self, c_system=None):
        if c_system is None:
            self._as_parameter_ = None
            self.natom = 0
            self.natom_move = 0
            self.charge = 0
            self.spin = 0
            self.model = "PM7"
            self.epsilon = 0.0
            self.atom = np.array([], dtype=np.int32)
            self.coord = np.array([], dtype=np.float64)
            self.nlattice = 0
            self.nlattice_move = 0
            self.pressure = 0.0
            self.lattice = np.array([], dtype=np.float64)
            self.tolerance = 1.0
            self.max_time = 3600
        else:
            if not isinstance(c_system, binding.c_mopac_system):
                raise TypeError("mismatch between mopac_system and c_mopac_system")
            self._as_parameter_ = ctypes.pointer(c_system)
            self.natom = c_system.natom
            self.natom_move = c_system.natom_move
            self.charge = c_system.charge
            self.spin = c_system.spin
            self.model = c_system.model
            self.epsilon = c_system.epsilon
            if c_system.natom > 0:
                self.atom = np.ctypeslib.as_array(c_system.atom, (c_system.natom,))
                self.coord = np.ctypeslib.as_array(c_system.coord, (3*c_system.natom,))
            else:
                self.atom = np.array([], dtype=np.int32)
                self.coord = np.array([], dtype=np.float64)
            self.nlattice = c_system.nlattice
            self.nlattice_move = c_system.nlattice_move
            self.pressure = c_system.pressure
            if c_system.nlattice > 0:
                self.lattice = np.ctypeslib.as_array(c_system.lattice, (3*c_system.nlattice,))
            else:
                self.lattice = np.array([], dtype=np.float64)
            self.tolerance = c_system.tolerance
            self.max_time = c_system.max_time
    def attach(self):
            if self._as_parameter_ is None:
                self._as_parameter_ = ctypes.pointer(binding.c_mopac_system())
            self._as_parameter_[0].natom = self.natom
            self._as_parameter_[0].natom_move = self.natom_move
            self._as_parameter_[0].charge = self.charge
            self._as_parameter_[0].spin = self.spin
            if self.model in self.model_dict:
                self._as_parameter_[0].model = self.model_dict[self.model]
            else:
                self._as_parameter_[0].model = self.model
            self._as_parameter_[0].epsilon = self.epsilon
            self._as_parameter_[0].atom = self.atom.ctypes.data_as(ctypes.POINTER(ctypes.c_int))
            self._as_parameter_[0].coord = self.coord.ctypes.data_as(ctypes.POINTER(ctypes.c_double))
            self._as_parameter_[0].nlattice = self.nlattice
            self._as_parameter_[0].nlattice_move = self.nlattice_move
            self._as_parameter_[0].pressure = self.pressure
            self._as_parameter_[0].lattice = self.lattice.ctypes.data_as(ctypes.POINTER(ctypes.c_double))
            self._as_parameter_[0].tolerance = self.tolerance
            self._as_parameter_[0].max_time = self.max_time

class mopac_properties:
    def __init__(self, c_properties, system):
        if not isinstance(c_properties, binding.c_mopac_properties):
            raise TypeError("mismatch between mopac_properties and c_mopac_properties")
        if not isinstance(system, (binding.c_mopac_system, mopac_system)):
            raise TypeError("mopac_properties initialization needs a valid mopac_system for size information")
        self._as_parameter_ = ctypes.pointer(c_properties)
        self.heat = c_properties.heat
        self.dipole = np.ctypeslib.as_array(c_properties.dipole)
        if system.natom > 0:
            self.charge = np.ctypeslib.as_array(c_properties.charge, (system.natom,))
        else:
            self.charge = np.array([], dtype=np.float64)
        if system.natom_move > 0:
            self.coord_update = np.ctypeslib.as_array(c_properties.coord_update, (3*system.natom_move,))
            self.coord_deriv = np.ctypeslib.as_array(c_properties.coord_deriv, (3*system.natom_move,))
        else:
            self.coord_update = np.array([], dtype=np.float64)
            self.coord_deriv = np.array([], dtype=np.float64)
        if c_properties.freq:
            self.freq = np.ctypeslib.as_array(c_properties.freq, (3*system.natom_move,))
        else:
            self.freq = np.array([], dtype=np.float64)
        if c_properties.disp:
            self.disp = np.ctypeslib.as_array(c_properties.disp, (3*system.natom_move, 3*system.natom_move))
        else:
            self.disp = np.array([], dtype=np.float64)
        bo_indptr = np.ctypeslib.as_array(c_properties.bond_index, (system.natom+1,))
        if bo_indptr[-1] > 0:
            bo_indices = np.ctypeslib.as_array(c_properties.bond_atom, (bo_indptr[-1],))
            bo_data = np.ctypeslib.as_array(c_properties.bond_order, (bo_indptr[-1],))
        else:
            bo_indices = np.array([], dtype=np.int32)
            bo_data = np.array([], dtype=np.float64)
        self.bond_order = sp.sparse.csc_matrix((bo_data, bo_indices, bo_indptr), shape=(system.natom, system.natom))
        if system.nlattice_move > 0:
            self.lattice_update = np.ctypeslib.as_array(c_properties.lattice_update, (3*system.nlattice_move,))
            self.lattice_deriv = np.ctypeslib.as_array(c_properties.lattice_deriv, (3*system.nlattice_move,))
        else:
            self.lattice_update = np.array([], dtype=np.float64)
            self.lattice_deriv = np.array([], dtype=np.float64)
        self.stress = np.ctypeslib.as_array(c_properties.stress)
        self.error_msg = []
        for i in range(c_properties.nerror):
            self.error_msg.append(c_properties.error_msg[i].decode().split('\x00', 1)[0])
    def detach(self):
        if self._as_parameter_ is not None:
            self.dipole = self.dipole.copy()
            self.charge = self.charge.copy()
            self.coord_update = self.coord_update.copy()
            self.coord_deriv = self.coord_deriv.copy()
            self.freq = self.freq.copy()
            self.disp = self.disp.copy()
            self.bond_order = self.bond_order.copy()
            self.lattice_update = self.lattice_update.copy()
            self.lattice_deriv = self.lattice_deriv.copy()
            self.stress = self.stress.copy()
            binding.libmopac.destroy_mopac_properties(self)
            self._as_parameter_ = None
    def __del__(self):
        if self._as_parameter_ is not None:
            binding.libmopac.destroy_mopac_properties(self)

class mopac_state:
    def __init__(self, c_state=None):
        if c_state is None:
            self._as_parameter_ = None
            self.mpack = 0
            self.uhf = False
            self.pa = np.array([], dtype=np.float64)
            self.pb = np.array([], dtype=np.float64)
        else:
            if not isinstance(c_state, binding.c_mopac_state):
                raise TypeError("mismatch between mopac_state and c_mopac_state")
            self._as_parameter_ = ctypes.pointer(c_state)
            self.update()
    def attach(self):
        if self._as_parameter_ is None:
            self._as_parameter_ = ctypes.pointer(binding.c_mopac_state())
            self._as_parameter_[0].mpack = self.mpack
            binding.libmopac.create_mopac_state(self)
        elif self._as_parameter_[0].mpack != self.mpack:
            binding.libmopac.destroy_mopac_state(self)
            self._as_parameter_[0].mpack = self.mpack
            binding.libmopac.create_mopac_state(self)
        if (ctypes.addressof(self.pa.ctypes.data_as(ctypes.POINTER(ctypes.c_double))) == ctypes.addressof(self._as_parameter_[0].pa)
            and ctypes.addressof(self.pb.ctypes.data_as(ctypes.POINTER(ctypes.c_double))) == ctypes.addressof(self._as_parameter_[0].pb)):
            return
        self._as_parameter_[0].uhf = int(self.uhf)
        if self.mpack > 0:
            pa_ref = np.ctypeslib.as_array(self._as_parameter_[0].pa, (self.mpack,))
            np.copyto(pa_ref, self.pa)
            self.pa = pa_ref
            if self.uhf:
                pb_ref = np.ctypeslib.as_array(self._as_parameter_[0].pb, (self.mpack,))
                np.copyto(pb_ref, self.pb)
                self.pb = pb_ref
            else:
                self.pb = np.array([], dtype=np.float64)
        else:
            self.pa = np.array([], dtype=np.float64)
            self.pb = np.array([], dtype=np.float64)
    def update(self):
        self.mpack = self._as_parameter_[0].mpack
        self.uhf = bool(self._as_parameter_[0].uhf)
        if self._as_parameter_[0].mpack > 0:
            self.pa = np.ctypeslib.as_array(self._as_parameter_[0].pa, (self._as_parameter_[0].mpack,))
        else:
            self.pa = np.array([], dtype=np.float64)
        if self._as_parameter_[0].mpack > 0 and self.uhf:
            self.pb = np.ctypeslib.as_array(self._as_parameter_[0].pb, (self._as_parameter_[0].mpack,))
        else:
            self.pb = np.array([], dtype=np.float64)
    def detach(self):
        if self._as_parameter_ is not None:
            self.pa = self.pa.copy()
            self.pb = self.pb.copy()
            binding.libmopac.destroy_mopac_state(self)
            self._as_parameter_ = None
    def __del__(self):
        if self._as_parameter_ is not None:
            binding.libmopac.destroy_mopac_state(self)

class mozyme_state:
    def __init__(self, c_state=None):
        if c_state is None:
            self._as_parameter_ = None
            self.numat = 0
            self.nbonds = np.array([], dtype=np.int32)
            self.ibonds = np.array([], dtype=np.int32)
            self.iorbs = np.array([], dtype=np.int32)
            self.ncf = np.array([], dtype=np.int32)
            self.nce = np.array([], dtype=np.int32)
            self.icocc = np.array([], dtype=np.int32)
            self.icvir = np.array([], dtype=np.int32)
            self.cocc = np.array([], dtype=np.float64)
            self.cvir = np.array([], dtype=np.float64)
        else:
            if not isinstance(c_state, binding.c_mozyme_state):
                raise TypeError("mismatch between mozyme_state and c_mozyme_state")
            self._as_parameter_ = ctypes.pointer(c_state)
            self.update()
    def attach(self):
        if self._as_parameter_ is None:
            self._as_parameter_ = ctypes.pointer(binding.c_mopac_state())
            self._as_parameter_[0].numat = self.numat
            self._as_parameter_[0].noccupied = len(self.ncf)
            self._as_parameter_[0].nvirtual = len(self.nce)
            self._as_parameter_[0].icocc_dim = len(self.icocc)
            self._as_parameter_[0].icvir_dim = len(self.icvir)
            self._as_parameter_[0].cocc_dim = len(self.cocc)
            self._as_parameter_[0].cvir_dim = len(self.cvir)
            binding.libmopac.create_mozyme_state(self)
        elif (self._as_parameter_[0].numat != self.numat
              or self._as_parameter_[0].noccupied != len(self.ncf)
              or self._as_parameter_[0].nvirtual != len(self.nce)
              or self._as_parameter_[0].icocc_dim != len(self.icocc)
              or self._as_parameter_[0].icvir_dim != len(self.icvir)
              or self._as_parameter_[0].cocc_dim != len(self.cocc)
              or self._as_parameter_[0].cvir_dim != len(self.cvir)):
            binding.libmopac.destroy_mozyme_state(self)
            self._as_parameter_[0].numat = self.numat
            self._as_parameter_[0].noccupied = len(self.ncf)
            self._as_parameter_[0].nvirtual = len(self.nce)
            self._as_parameter_[0].icocc_dim = len(self.icocc)
            self._as_parameter_[0].icvir_dim = len(self.icvir)
            self._as_parameter_[0].cocc_dim = len(self.cocc)
            self._as_parameter_[0].cvir_dim = len(self.cvir)
            binding.libmopac.create_mozyme_state(self)
        if (ctypes.addressof(self.nbonds.ctypes.data_as(ctypes.POINTER(ctypes.c_int))) == ctypes.addressof(self._as_parameter_[0].nbonds)
            and ctypes.addressof(self.ibonds.ctypes.data_as(ctypes.POINTER(ctypes.c_int))) == ctypes.addressof(self._as_parameter_[0].ibonds)
            and ctypes.addressof(self.iorbs.ctypes.data_as(ctypes.POINTER(ctypes.c_int))) == ctypes.addressof(self._as_parameter_[0].iorbs)
            and ctypes.addressof(self.ncf.ctypes.data_as(ctypes.POINTER(ctypes.c_int))) == ctypes.addressof(self._as_parameter_[0].ncf)
            and ctypes.addressof(self.nce.ctypes.data_as(ctypes.POINTER(ctypes.c_int))) == ctypes.addressof(self._as_parameter_[0].nce)
            and ctypes.addressof(self.icocc.ctypes.data_as(ctypes.POINTER(ctypes.c_int))) == ctypes.addressof(self._as_parameter_[0].icocc)
            and ctypes.addressof(self.icvir.ctypes.data_as(ctypes.POINTER(ctypes.c_int))) == ctypes.addressof(self._as_parameter_[0].icvir)
            and ctypes.addressof(self.cocc.ctypes.data_as(ctypes.POINTER(ctypes.c_double))) == ctypes.addressof(self._as_parameter_[0].cocc)
            and ctypes.addressof(self.cvir.ctypes.data_as(ctypes.POINTER(ctypes.c_double))) == ctypes.addressof(self._as_parameter_[0].cvir)):
            return
        if self.numat > 0:
            nbonds_ref = np.ctypeslib.as_array(self._as_parameter_[0].nbonds, (self.numat,))
            np.copyto(nbonds_ref, self.nbonds)
            self.nbonds = nbonds_ref
            ibonds_ref = np.ctypeslib.as_array(self._as_parameter_[0].ibonds, (9*self.numat,))
            np.copyto(ibonds_ref, self.ibonds)
            self.ibonds = ibonds_ref
            iorbs_ref = np.ctypeslib.as_array(self._as_parameter_[0].iorbs, (self.numat,))
            np.copyto(iorbs_ref, self.iorbs)
            self.iorbs = iorbs_ref
        else:
            self.nbonds = np.array([], dtype=np.int32)
            self.ibonds = np.array([], dtype=np.int32)
            self.iorbs = np.array([], dtype=np.int32)
        if self._as_parameter_[0].noccupied > 0:
            ncf_ref = np.ctypeslib.as_array(self._as_parameter_[0].ncf, (self._as_parameter_[0].noccupied,))
            np.copyto(ncf_ref, self.ncf)
            self.ncf = ncf_ref
        else:
            self.ncf = np.array([], dtype=np.int32)
        if self._as_parameter_[0].nvirtual > 0:
            nce_ref = np.ctypeslib.as_array(self._as_parameter_[0].nce, (self._as_parameter_[0].nvirtual,))
            np.copyto(nce_ref, self.nce)
            self.nce = nce_ref
        else:
            self.nce = np.array([], dtype=np.int32)
        if self._as_parameter_[0].icocc_dim > 0:
            icocc_ref = np.ctypeslib.as_array(self._as_parameter_[0].icocc, (self._as_parameter_[0].icocc_dim,))
            np.copyto(icocc_ref, self.icocc)
            self.icocc = icocc_ref
        else:
            self.icocc = np.array([], dtype=np.int32)
        if self._as_parameter_[0].icvir_dim > 0:
            icvir_ref = np.ctypeslib.as_array(self._as_parameter_[0].icvir, (self._as_parameter_[0].icvir_dim,))
            np.copyto(icvir_ref, self.icvir)
            self.icvir = icvir_ref
        else:
            self.icvir = np.array([], dtype=np.int32)
        if self._as_parameter_[0].cocc_dim > 0:
            cocc_ref = np.ctypeslib.as_array(self._as_parameter_[0].cocc, (self._as_parameter_[0].cocc_dim,))
            np.copyto(cocc_ref, self.cocc)
            self.cocc = cocc_ref
        else:
            self.cocc = np.array([], dtype=np.float64)
        if self._as_parameter_[0].cvir_dim > 0:
            cvir_ref = np.ctypeslib.as_array(self._as_parameter_[0].cvir, (self._as_parameter_[0].cvir_dim,))
            np.copyto(cvir_ref, self.cvir)
            self.cvir = cvir_ref
        else:
            self.cvir = np.array([], dtype=np.float64)
    def update(self):
        self.numat = self._as_parameter_[0].numat
        if self._as_parameter_[0].numat > 0:
            self.nbonds = np.ctypeslib.as_array(self._as_parameter_[0].nbonds, (self._as_parameter_[0].numat,))
            self.ibonds = np.ctypeslib.as_array(self._as_parameter_[0].ibonds, (9*self._as_parameter_[0].numat,))
            self.iorbs = np.ctypeslib.as_array(self._as_parameter_[0].iorbs, (self._as_parameter_[0].numat,))
            if self._as_parameter_[0].noccupied > 0:
                self.ncf = np.ctypeslib.as_array(self._as_parameter_[0].ncf, (self._as_parameter_[0].noccupied,))
            else:
                self.ncf = np.array([], dtype=np.int32)
            if self._as_parameter_[0].nvirtual > 0:
                self.nce = np.ctypeslib.as_array(self._as_parameter_[0].nce, (self._as_parameter_[0].nvirtual,))
            else:
                self.nce = np.array([], dtype=np.int32)
            if self._as_parameter_[0].icocc_dim > 0:
                self.icocc = np.ctypeslib.as_array(self._as_parameter_[0].icocc, (self._as_parameter_[0].icocc_dim,))
            else:
                self.icocc = np.array([], dtype=np.int32)
            if self._as_parameter_[0].icvir_dim > 0:
                self.icvir = np.ctypeslib.as_array(self._as_parameter_[0].icvir, (self._as_parameter_[0].icvir_dim,))
            else:
                self.icvir = np.array([], dtype=np.int32)
            if self._as_parameter_[0].cocc_dim > 0:
                self.cocc = np.ctypeslib.as_array(self._as_parameter_[0].cocc, (self._as_parameter_[0].cocc_dim,))
            else:
                self.cocc = np.array([], dtype=np.float64)
            if self._as_parameter_[0].cvir_dim > 0:
                self.cvir = np.ctypeslib.as_array(self._as_parameter_[0].cvir, (self._as_parameter_[0].cvir_dim,))
            else:
                self.cvir = np.array([], dtype=np.float64)
        else:
            self.nbonds = np.array([], dtype=np.int32)
            self.ibonds = np.array([], dtype=np.int32)
            self.iorbs = np.array([], dtype=np.int32)
            self.ncf = np.array([], dtype=np.int32)
            self.nce = np.array([], dtype=np.int32)
            self.icocc = np.array([], dtype=np.int32)
            self.icvir = np.array([], dtype=np.int32)
            self.cocc = np.array([], dtype=np.float64)
            self.cvir = np.array([], dtype=np.float64)
    def detach(self):
        if self._as_parameter_ is not None:
            self.nbonds = self.nbonds.copy()
            self.ibonds = self.ibonds.copy()
            self.iorbs = self.iorbs.copy()
            self.ncf = self.ncf.copy()
            self.nce = self.nce.copy()
            self.icocc = self.icocc.copy()
            self.icvir = self.icvir.copy()
            self.cocc = self.cocc.copy()
            self.cvir = self.cvir.copy()
            binding.libmopac.destroy_mozyme_state(self)
            self._as_parameter_ = None
    def __del__(self):
        if self._as_parameter_ is not None:
            binding.libmopac.destroy_mozyme_state(self)

def scf(system, state):
    if not isinstance(system, mopac_system):
        raise TypeError("1st argument of scf must be a mopac_system")
    if not isinstance(state, (mopac_state, mozyme_state)):
        raise TypeError("2nd argument of scf must be a mopac_state or mozyme_state")
    c_properties = binding.c_mopac_properties()
    system.attach()
    state.attach()
    if isinstance(state, mopac_state):
        binding.libmopac.mopac_scf(system, state, ctypes.byref(c_properties))
    else:
        binding.libmopac.mozyme_scf(system, state, ctypes.byref(c_properties))
    state.update()
    return mopac_properties(c_properties, system)

def relax(system, state):
    if not isinstance(system, mopac_system):
        raise TypeError("1st argument of relax must be a mopac_system")
    if not isinstance(state, (mopac_state, mozyme_state)):
        raise TypeError("2nd argument of relax must be a mopac_state or mozyme_state")
    c_properties = binding.c_mopac_properties()
    system.attach()
    state.attach()
    if isinstance(state, mopac_state):
        binding.libmopac.mopac_relax(system, state, ctypes.byref(c_properties))
    else:
        binding.libmopac.mozyme_relax(system, state, ctypes.byref(c_properties))
    state.update()
    return mopac_properties(c_properties, system)

def vibe(system, state):
    if not isinstance(system, mopac_system):
        raise TypeError("1st argument of vibe must be a mopac_system")
    if not isinstance(state, (mopac_state, mozyme_state)):
        raise TypeError("2nd argument of vibe must be a mopac_state or mozyme_state")
    c_properties = binding.c_mopac_properties()
    system.attach()
    state.attach()
    if isinstance(state, mopac_state):
        binding.libmopac.mopac_vibe(system, state, ctypes.byref(c_properties))
    else:
        binding.libmopac.mozyme_vibe(system, state, ctypes.byref(c_properties))
    state.update()
    return mopac_properties(c_properties, system)

def from_input(filepath):
    if not os.path.isfile(filepath):
        raise ValueError("argument of from_input is not a valid file")
    status = binding.libmopac.run_mopac_from_input(ctypes.create_string_buffer(os.fsencode(filepath)))
    filename = os.fsdecode(filepath)
    if filename[-4:].upper() == ".MOP" or filename[-4:].upper() == ".DAT" or filename[-4:].upper() == ".ARC":
        outpath = filename[:-4] + ".out"
    else:
        outpath = filename + ".out"
    return (outpath, not bool(status))
