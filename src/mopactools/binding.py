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

# Python bindings for the MOPAC API
import os
import itertools
from ctypes import *

# Define the API data structures
class c_mopac_system(Structure):
    """ ctypes binding of the mopac_system struct """
    _fields_ = [("natom", c_int),
                ("natom_move", c_int),
                ("charge", c_int),
                ("spin", c_int),
                ("model", c_int),
                ("epsilon", c_double),
                ("atom", POINTER(c_int)),
                ("coord", POINTER(c_double)),
                ("nlattice", c_int),
                ("nlattice_move", c_int),
                ("pressure", c_double),
                ("lattice", POINTER(c_double)),
                ("tolerance", c_double),
                ("max_time", c_int)]
class c_mopac_properties(Structure):
    """ ctypes binding of the mopac_properties struct """
    _fields_ = [("heat", c_double),
                ("dipole", c_double * 3),
                ("charge", POINTER(c_double)),
                ("coord_update", POINTER(c_double)),
                ("coord_deriv", POINTER(c_double)),
                ("freq", POINTER(c_double)),
                ("disp", POINTER(c_double)),
                ("bond_index", POINTER(c_int)),
                ("bond_atom", POINTER(c_int)),
                ("bond_order", POINTER(c_double)),
                ("lattice_update", POINTER(c_double)),
                ("lattice_deriv", POINTER(c_double)),
                ("stress", c_double * 6),
                ("nerror", c_int),
                ("error_msg", POINTER(c_char_p))]
class c_mopac_state(Structure):
    """ ctypes binding of the mopac_state struct """
    _fields_ = [("mpack", c_int),
                ("uhf", c_int),
                ("pa", POINTER(c_double)),
                ("pb", POINTER(c_double))]
class c_mozyme_state(Structure):
    """ ctypes binding of the mozyme_state struct """
    _fields_ = [("numat", c_int),
                ("nbonds", POINTER(c_int)),
                ("ibonds", POINTER(c_int)),
                ("iorbs", POINTER(c_int)),
                ("noccupied", c_int),
                ("ncf", POINTER(c_int)),
                ("nvirtual", c_int),
                ("nce", POINTER(c_int)),
                ("icocc_dim", c_int),
                ("icocc", POINTER(c_int)),
                ("icvir_dim", c_int),
                ("icvir", POINTER(c_int)),
                ("cocc_dim", c_int),
                ("cocc", POINTER(c_double)),
                ("cvir_dim", c_int),
                ("cvir", POINTER(c_double))]

# Load the MOPAC shared library, check all reasonable library names & paths
mopactools_path = os.path.dirname(__file__)
libmopac_name = ["libmopac.so", "libmopac.dylib", "libmopac.dll", "mopac.dll"]
libmopac_dir = [os.path.join(mopactools_path, "lib"), None]
for name, dir in itertools.product(libmopac_name, libmopac_dir):
    if dir is None:
        libmopac_path = name
    else:
        libmopac_path = os.path.join(dir, name)
    try:
        libmopac = CDLL(libmopac_path)
        break
    except OSError:
        pass
try:
    libmopac
except NameError:
    raise OSError("MOPAC library could not be found in the system path or mopactools/lib directory")

# Specify the argument lists of the API functions
libmopac.mopac_scf.argtypes = [POINTER(c_mopac_system), POINTER(c_mopac_state), POINTER(c_mopac_properties)]
libmopac.mopac_relax.argtypes = [POINTER(c_mopac_system), POINTER(c_mopac_state), POINTER(c_mopac_properties)]
libmopac.mopac_vibe.argtypes = [POINTER(c_mopac_system), POINTER(c_mopac_state), POINTER(c_mopac_properties)]
libmopac.mozyme_scf.argtypes = [POINTER(c_mopac_system), POINTER(c_mozyme_state), POINTER(c_mopac_properties)]
libmopac.mozyme_relax.argtypes = [POINTER(c_mopac_system), POINTER(c_mozyme_state), POINTER(c_mopac_properties)]
libmopac.mozyme_vibe.argtypes = [POINTER(c_mopac_system), POINTER(c_mozyme_state), POINTER(c_mopac_properties)]
libmopac.create_mopac_state.argtypes = [POINTER(c_mopac_state)]
libmopac.create_mozyme_state.argtypes = [POINTER(c_mozyme_state)]
libmopac.destroy_mopac_properties.argtypes = [POINTER(c_mopac_properties)]
libmopac.destroy_mopac_state.argtypes = [POINTER(c_mopac_state)]
libmopac.destroy_mozyme_state.argtypes = [POINTER(c_mozyme_state)]
libmopac.run_mopac_from_input.argtypes = [c_char_p]
libmopac.get_mopac_version.argtypes = [c_char_p]
