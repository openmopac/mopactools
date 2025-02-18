Example
=======

A minimal example using the MOPAC API to calculate the electronic ground state of water at a rough molecular geometry and print the heat of formation:

.. code::

    from mopactools import api
    import numpy as np

    system = api.MopacSystem()
    system.natom = 3
    system.atom = ["H", "H", "O"]
    system.coord = np.array([0.76, 0.59, 0, -0.76, 0.59, 0, 0, 0, 0])
    state = api.MopacState()
    properties = api.from_data(system, state)
    print("heat of formation = ", properties.heat)
