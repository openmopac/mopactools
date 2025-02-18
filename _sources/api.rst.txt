API
===

MOPAC's application programming interface (API) is intended to streamline access and improve performance for the most commonly used features,
while also providing access to MOPAC's conventional file-based operations.

The Python API is located in the ``api.py`` module of mopactools, which can be imported as::

    from mopactools import api

Conventional MOPAC calculations using an input file can then be run by using the function:

.. autofunction:: mopactools.api.from_file

The rest of the API functionality is presently limited to the most commonly used features in MOPAC, which is the calculation of an electronic ground
state, with an option to relax the atomic coordinates to their local energy minimum, and a further option to calculate vibrational properties at that
relaxed geometry. These features can be run using either the conventional MOPAC solver or the linear-scaling MOZYME solver. The MOZYME solver is faster
for larger systems, but only works under a narrower set of conditions. Specifically, the MOZYME solver only works for closed-shell systems
with a covalent and/or ionic bonding structure that can be inferred from the geometry. Otherwise, MOZYME will fail and report an error.

This functionality is accessed through the function:

.. autofunction:: mopactools.api.from_data

The Python data layout follows the underlying Fortran (with C bindings) data layout. There is a ``MopacSystem`` object for the atomistic system data
and computational options that define the calculation, a ``MopacState`` or ``MozymeState`` object that describes the electronic ground state at the beginning
and end of the calculation, and a ``MopacProperties`` object that stores the important physical properties that were evaluated during the calculation.
The ``MopacSystem`` object is strictly an input that is not altered by the calculation, while the ``MopacState`` and ``MozymeState`` objects act as both
inputs and outputs, so the input object is altered by the function call. The initial value of ``MopacState`` or ``MozymeState`` is used as the initial guess
for the self-consistent field (SCF) cycle, and the final value is the converged ground state at the end of the SCF cycle. An uninitialized ``MopacState`` or
``MozymeState`` instance correspond to the default guesses for the electronic ground state (an atomic guess for MOPAC and a Lewis-bonded guess for MOZYME).
The ``MopacProperties`` object is strictly an output, therefore it is returned by the API call rather than specified as a function argument.

.. note::
   The relaxed geometry of a system is stored in ``MopacProperties.coord_update`` rather than ``MopacSystem.coord`` because ``MopacSystem`` is not
   altered by calls to ``api.from_data``. The user has to update the geometry themselves if they want it to be used in subsequent calculations.

The specification of this API data structures is as follows:

.. autoclass:: mopactools.api.MopacSystem

.. autoclass:: mopactools.api.MopacState

.. autoclass:: mopactools.api.MozymeState

.. autoclass:: mopactools.api.MopacProperties

Finally, version information about the MOPAC shared library is provided by the module variable:

.. autodata:: mopactools.api.VERSION
   :annotation:

.. warning::
   None of the MOPAC API functionality is thread-safe. Typical usage of Python occurs in a single thread, so this is not relevant to normal operations.
   However, Python is slowly expanding its multi-threading support, so this may eventually be a problem. Making MOPAC itself thread-safe is not a realistic
   possibility without completely rewriting it, but the Python API wrapper could be made thread-safe if each thread loads a separate instance
   of the MOPAC shared library.
