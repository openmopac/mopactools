import numpy as np
import mopactools.api as api

def mopac_water_in():
    system = api.mopac_system()
    system.natom = 3
    system.natom_move = 3
    system.atom = np.array([1, 1, 8])
    system.coord = np.array([0.76, 0.59, 0, -0.76, 0.59, 0, 0, 0, 0])
    return system

def mopac_water_out():
    properties = api.mopac_properties()
    properties.heat = -57.76975
    properties.coord_deriv = np.array([2.307865, 2.742432, 0, -2.307865, 2.711610, 0, 0, -5.454042, 0])
    properties.charge = np.array([0.322260, 0.322260, -0.644520])
    properties.dipole = np.array([0, 2.147, 0, 0])
    properties.stress = np.array([0, 0, 0, 0, 0, 0], dtype=float64)
    bond_index = np.array([1, 3, 5, 8])
    bond_atom = np.array([1, 3, 2, 3, 1, 2, 3])
    bond_order = np.array([0.896, 0.895, 0.896, 0.895, 0.895, 0.895, 1.791])
    properties.bond_order = sp.sparse.csc_matrix((bond_order, bond_atom, bond_index), shape=(3, 3))
    properties.error_msg = []
    return properties

def test_mopac_scf():
    system = mopac_water_in()
    state = api.mopac_state()
    properties = api.scf(system, state)
    print(properties)
