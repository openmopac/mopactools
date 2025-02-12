# Copyright 2025 Virginia Polytechnic Institute and State University
#
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
#    http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.

# Python wrappers for the MOPAC API
import os
import copy
import math
import ctypes
import numpy as np
import scipy
from . import binding

class mopac_system:
    """Information defining a MOPAC calculation of an atomistic system

    Attributes
    ----------
    natom : int
        The number of atoms in the system.
    natom_move : int
        The number of atoms in the system that are allowed to move.
        The first natom_move atoms in the list can move, and the last natom-natom_move atoms cannot.
    charge : int
        The total charge on the system.
    spin : int
        The excess electronic spin on the system.
        For example, spin=0 is a singlet and spin=1 is a triplet for an even number of electrons,
        and spin=0 is a doublet for an odd number of electrons.
    model : {"PM7", "PM6-D3H4", "PM6-ORG", "PM6", "AM1", "RM1"} or int
        The semiempirical model used for the calculation, or the corresponding int from model_dict.
    epsilon : float
        The dielectric constant of the implicit solvent that the atomistic system is embedded in.
        There is no implicit solvent if epsilon=1.0, the dielectric constant of empty vacuum.
    atom : ndarray of ints, or list of strs or ints
        An array containing the atomic numbers of the atoms being simulated.
        This can also be specified as a list of atomic numbers and/or element symbols,
        which will be converted into an array when mopac_system is attached to the MOPAC library.
    coord : ndarray
        An array containing the Cartesian coordinates of atoms in xyz order and Angstrom units.
    nlattice : int
        The number of lattice/translation vectors that define periodic boundary conditions.
    nlattice_move : int
        The number of lattice/translation vectors that are allowed to move.
        The first nlattice_move vectors can move, and the last nlattice-nlattice_move vectors cannot.
    pressure : float
        The external pressure applied to the system in units of Gigapascals (GPa).
    lattice : ndarray
        An array contains the Cartesian coordinates of lattice vectors in xyz order and Angstrom units.
    tolerance : float
        The relative numerical tolerances of the MOPAC calculation, with tolerance=1.0 as default.
        The default is appropriate for most calculations, and this is only intended to test convergence.
    max_time : int
        The maximum run time that a MOPAC calculation is allowed to use in seconds.

    Methods
    -------
    attach()
        Attach to a c_mopac_system instance that uses Python memory.
    """
    model_dict = { "PM7":0, "PM6-D3H4":1, "PM6-ORG":2, "PM6":3, "AM1":4, "RM1":5 }
    periodic_table = { "H":1, "He":2, "Li":3, "Be":4, "B":5, "C":6, "N":7, "O":8, "F":9, "Ne":10,
        "Na":11, "Mg":12, "Al":13, "Si":14, "P":15, "S":16, "Cl":17, "Ar":18, "K":19, "Ca":20,
        "Sc":21, "Ti":22, "V":23, "Cr":24, "Mn":25, "Fe":26, "Co":27, "Ni":28, "Cu":29, "Zn":30,
        "Ga":31, "Ge":32, "As":33, "Se":34, "Br":35, "Kr":36, "Rb":37, "Sr":38, "Y":39, "Zr":40,
        "Nb":41, "Mo":42, "Tc":43, "Ru":44, "Rh":45, "Pd":46, "Ag":47, "Cd":48, "In":49, "Sn":50,
        "Sb":51, "Te":52, "I":53, "Xe":54, "Cs":55, "Ba":56, "La":57, "Ce":58, "Pr":59, "Nd":60,
        "Pm":61, "Sm":62, "Eu":63, "Gd":64, "Tb":65, "Dy":66, "Ho":67, "Er":68, "Tm":69, "Yb":70,
        "Lu":71, "Hf":72, "Ta":73, "W":74, "Re":75, "Os":76, "Ir":77, "Pt":78, "Au":79, "Hg":80,
        "Tl":81, "Pb":82, "Bi":83, "Po":84, "At":85, "Rn":86, "Fr":87, "Ra":88, "Ac":89, "Th":90,
        "Pa":91, "U":92, "Np":93, "Pu":94, "Am":95, "Cm":96, "Bk":97, "Cf":98 }
    def __init__(self, c_system=None):
        """Initialize to default values, or the values in a c_mopac_system instance if provided"""
        if c_system is None:
            self._as_parameter_ = None
            self.natom = 0
            self.natom_move = 0
            self.charge = 0
            self.spin = 0
            self.model = "PM7"
            self.epsilon = 1.0
            self.atom = []
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
            if not isinstance(self.atom, np.ndarray):
                old_atom = self.atom
                self.atom = np.empty(self.natom, dtype=np.int32)
                for i in range(self.natom):
                    if isinstance(old_atom[i], int):
                        self.atom[i] = old_atom[i]
                    else:
                        if len(old_atom[i].strip()) == 2:
                            element = f"{old_atom[i][0].upper()}{old_atom[i][1].lower()}"
                        elif len(old_atom[i].strip()) == 1:
                            element = f"{old_atom[i][0].upper()}"
                        else:
                            raise TypeError(f"unknown element symbol, {old_atom[i]}, in mopac_system")
                        self.atom[i] = self.periodic_table[element]
            self._as_parameter_[0].atom = self.atom.ctypes.data_as(ctypes.POINTER(ctypes.c_int))
            self._as_parameter_[0].coord = self.coord.ctypes.data_as(ctypes.POINTER(ctypes.c_double))
            self._as_parameter_[0].nlattice = self.nlattice
            self._as_parameter_[0].nlattice_move = self.nlattice_move
            self._as_parameter_[0].pressure = self.pressure
            self._as_parameter_[0].lattice = self.lattice.ctypes.data_as(ctypes.POINTER(ctypes.c_double))
            self._as_parameter_[0].tolerance = self.tolerance
            self._as_parameter_[0].max_time = self.max_time

class mopac_properties:
    """Properties generated by a MOPAC calculation

    Attributes
    ----------
    heat : float
        The heat of formation of the system in units of kcal/mol.
    dipole : ndarray
        The dipole vector of the system in units of Debye.
    charge : ndarray
        The partial charge on each atom in the system.
    coord_update : ndarray
        The updated Cartesian coordinates of the moveable atoms in xyz format and Angstrom units.
    coord_deriv : ndarray
        The heat-of-formation gradients of the moveable atoms in xyz format and kcal/mol/Angstrom units.
    freq : ndarray
        The viibrational frequencies of the normal modes of the moveable atoms in units of 1/cm.
    disp : ndarray
        The displacement vectors of the normal modes of the moveable atoms in xyz format.
    bond_order : csc_matrix
        The bond-order matrix with diagonal entries corresponding to the valence of atoms.
    lattice_update : ndarray
        The updated Cartesian coordinates of the moveable lattice vectors in xyz format and Angstrom units.
    lattice_deriv : ndarray
        The heat-of-formation gradients of the moveable lattice vectors in xyz format and kcal/mol/Angstrom units.
    stress : ndarray
        The stress tensor in Voigt form (xx, yy, zz, yz, xz, xy) and Gigapascal units.

    Methods
    -------
    detach()
        Detach from a c_mopac_properties instance and move memory into Python.
    """
    def __init__(self, c_properties=None, system=None):
        """Initial to the values of a c_mopac_properties instance using size info from a (c_)mopac_system instance"""
        if c_properties is None:
            self._as_parameter_ = None
            self.heat = None
            self.dipole = None
            self.charge = None
            self.coord_update = None
            self.coord_deriv = None
            self.freq = None
            self.disp = None
            self.bond_order = None
            self.lattice_update = None
            self.lattice_deriv = None
            self.stress = None
        else:
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
            self.bond_order = scipy.sparse.csc_matrix((bo_data, bo_indices, bo_indptr), shape=(system.natom, system.natom))
            if system.nlattice_move > 0:
                self.lattice_update = np.ctypeslib.as_array(c_properties.lattice_update, (3*system.nlattice_move,))
                self.lattice_deriv = np.ctypeslib.as_array(c_properties.lattice_deriv, (3*system.nlattice_move,))
            else:
                self.lattice_update = np.array([], dtype=np.float64)
                self.lattice_deriv = np.array([], dtype=np.float64)
            self.stress = np.ctypeslib.as_array(c_properties.stress)
            self.error_msg = []
            for i in range(c_properties.nerror):
                self.error_msg.append(c_properties.error_msg[i].value.decode('utf-8'))
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
    """Electronic ground state in a conventional MOPAC calculation"""

    Attributes
    ----------
    mpack : int
        Size of the packed density matrix.
        An uninitialized MOPAC ground state is denoted by mpack=0.
    uhf : bool
        Flag specifying if the ground state uses spin-unrestricted electronic orbitals.
    pa : ndarray
        Electronic density matrix for the alpha electrons in a packed lower-triangle format.
    pb : ndarray
        If bool=True, electronic density matrix for the beta electrons in a packed lower-triangle format.

    Methods
    -------
    attach()
        Attach to a c_mopac_state instance that and move memory into MOPAC.
    detach()
        Detach from a c_mopac_state instance and move memory into Python.
    """
    def __init__(self, c_state=None):
        """Initialize to an uninitialized ground state, or the values in a c_mopac_state instance if provided"""
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
    """Electronic ground state in a linear-scaling MOZYME calculation"""
    
    Attributes
    ----------
    numat : int
        Number of "real" atoms with active electrons.
        An uninitialized MOZYME ground state is denoted by numat=0.
    nbonds : ndarray of ints
        Number of Lewis bonds for each real atom.
    ibonds : ndarray of ints
        List of Lewis-bonded real atoms for each real atom.
    iorbs : ndarray of ints
        Number of orbitals for each real atom.
    noccupied : int
        Number of occupied local molecular orbitals (LMOs).
    ncf : ndarray of ints
        Number of atoms in each occupied LMO.
    nvirtual : int
        Number of virtual LMOs.
    nce : ndarray of ints
        Number of atoms in each virtual LMO.
    icocc : ndarray of ints
        1-based index of each real atom in each occupied LMO.
    icvir : ndarray of ints
        1-based index of each real atom in each virtual LMO.
    cocc : ndarray
        orbital coefficients of each real atom in each occupied LMO.
    cvir : ndarray
        orbital coefficients of each real atom in each virtual LMO.

    Methods
    -------
    attach()
        Attach to a c_mozyme_state instance that and move memory into MOPAC.
    detach()
        Detach from a c_mozyme_state instance and move memory into Python.
    """
    def __init__(self, c_state=None):
        """Initialize to an uninitialized ground state, or the values in a c_mozyme_state instance if provided"""
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
    """Calculate an electronic ground state using a self-consistent field (SCF) cycle in MOPAC.

    Parameters
    ----------
    system : mopac_system
        The description of the atomistic system to be calculated.
    state : mopac_state or mozyme_state
        The electronic ground state that the calculation is initialized to.
        The type is used to determine whether to run a MOPAC or MOZYME calculation.

    Returns
    -------
    properties : mopac_properties
        The physical properties calculated by MOPAC,
        except vibrational frequencies, properties.freq, and displacement vectors, properties.disp,
        which are both set to None.
    """
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
    """Calculate an atomistic ground state by relaxing the geometry in MOPAC.

    Parameters
    ----------
    system : mopac_system
        The description of the atomistic system to be calculated.
    state : mopac_state or mozyme_state
        The electronic ground state that the calculation is initialized to.
        The type is used to determine whether to run a MOPAC or MOZYME calculation.

    Returns
    -------
    properties : mopac_properties
        The physical properties calculated by MOPAC,
        except vibrational frequencies, properties.freq, and displacement vectors, properties.disp,
        which are both set to None.
    """
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
    """Calculate the vibrational properties of a ground-state geometry in MOPAC.

    Parameters
    ----------
    system : mopac_system
        The description of the atomistic system to be calculated.
    state : mopac_state or mozyme_state
        The electronic ground state that the calculation is initialized to.
        The type is used to determine whether to run a MOPAC or MOZYME calculation.

    Returns
    -------
    properties : mopac_properties
        The physical properties calculated by MOPAC,
        including vibrational frequencies, properties.freq, and displacement vectors, properties.disp.
    """
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

def run_input(filepath):
    """Run a MOPAC calculation using the input file located at filepath.
    This function is equivalent to running MOPAC on the command line.
    For example, an input file at the path dir/name.mop should generate an ouput file at the path dir/name.out.
    """
    if not os.path.isfile(filepath):
        raise ValueError("argument of run_input is not a valid file")
    status = binding.libmopac.run_mopac_from_input(ctypes.create_string_buffer(os.fsencode(filepath)))
    return bool(status)

def version():
    """Return the version of MOPAC"""
    buffer = ctypes.create_string_buffer(21)
    binding.libmopac.get_mopac_version(buffer)
    return buffer.value.decode('utf-8')
