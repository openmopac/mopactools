import os
import numpy as np
import scipy
import pytest
import mopactools.api as api

@pytest.fixture
def setup_input():
    input_name = "test.mop"
    output_name = "test.out"
    input = open(input_name, "w")
    input.write('''1SCF
                   ---
                   ---
                   H     0.76 0.59 0.0
                   H    -0.76 0.59 0.0
                   O     0.0 0.0 0.0
                   ''')
    input.close()
    yield input_name, output_name
    for file in [input_name, output_name, "text.arc"]:
        try:
            os.remove(file)
        except FileNotFoundError:
            pass

def test_run_input(setup_input):
    input_name, output_name = setup_input
    error = api.run_input(input_name)
    assert not error, "run_input caused an error"
    output = open(output_name, "r")
    output_txt = output.read()
    assert "FINAL HEAT OF FORMATION =        -57.769" in output_txt, "run_input output file is missing results"

def mopac_water_in():
    system = api.mopac_system()
    system.natom = 3
    system.natom_move = 3
    system.atom = ["H", "H", "O"]
    system.coord = np.array([0.76, 0.59, 0, -0.76, 0.59, 0, 0, 0, 0])
    return system

def mopac_water_out():
    properties = api.mopac_properties()
    properties.heat = -57.76975
    properties.coord_update = np.array([0.76, 0.59, 0, -0.76, 0.59, 0, 0, 0, 0])
    properties.coord_deriv = np.array([2.307865, 2.742432, 0, -2.307865, 2.711610, 0, 0, -5.454042, 0])
    properties.charge = np.array([0.322260, 0.322260, -0.644520])
    properties.dipole = np.array([0, 2.147, 0])
    properties.stress = np.zeros(6)
    bond_index = np.array([0, 2, 4, 7])
    bond_atom = np.array([0, 2, 1, 2, 0, 1, 2])
    bond_order = np.array([0.896, 0.895, 0.896, 0.895, 0.895, 0.895, 1.791])
    properties.bond_order = scipy.sparse.csc_matrix((bond_order, bond_atom, bond_index), shape=(3, 3))
    properties.lattice_update = np.array([], dtype=np.float64)
    properties.lattice_deriv = np.array([], dtype=np.float64)
    properties.error_msg = []
    return properties

def property_distance(prop1, prop2):
    distance = np.zeros(12)
    distance[0] = np.abs(prop1.heat - prop2.heat)
    distance[1] = np.linalg.norm(prop1.dipole - prop2.dipole)
    distance[2] = np.linalg.norm(prop1.charge - prop2.charge)
    if prop1.coord_update is not None and len(prop1.coord_update) > 0:
        distance[3] = np.linalg.norm(prop1.coord_update - prop2.coord_update)
        distance[4] = np.linalg.norm(prop1.coord_deriv - prop2.coord_deriv)
    if prop1.freq is not None and len(prop1.freq) > 0:
        distance[5] = np.linalg.norm(prop1.freq - prop2.freq)
        distance[6] = np.linalg.norm(prop1.disp@prop1.freq@prop1.disp.T - prop2.disp@prop2.freq@prop2.disp.T)
    distance[7] = np.linalg.norm((prop1.bond_order - prop2.bond_order).todense())
    if prop1.lattice_update is not None and len(prop1.lattice_update) > 0:
        distance[8] = np.linalg.norm(prop1.lattice_update - prop2.lattice_update)
        distance[9] = np.linalg.norm(prop1.lattice_deriv - prop2.lattice_deriv)
    distance[10] = np.linalg.norm(prop1.stress - prop2.stress)
    distance[11] = np.abs(len(prop1.error_msg) - len(prop2.error_msg))
    return distance

def test_mopac_scf():
    system = mopac_water_in()
    state = api.mopac_state()
    properties = api.scf(system, state)
    ref = mopac_water_out()
    dist = property_distance(properties, ref)
    assert np.linalg.norm(dist, ord=np.inf) < 0.002, "mopac_scf calculation differs significantly from reference"
