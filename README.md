# mopactools

This Python package is primarily a Python wrapper for API calls to the [Molecular Orbital PACkage (MOPAC)](https://github.com/openmopac/mopac),
which includes launching conventional MOPAC calculations from an input file. For supported environments, it packages the MOPAC shared library
and can be used as a dependency for other Python packages that run MOPAC calculations. It is also intended to be
collection of pre-processing and post-processing tools to enhance the productivity of MOPAC users.

User contributions of new tools are welcome and encouraged. As an open-source project, the future of MOPAC now depends on contributions
from its large user base and the broader community of computational chemistry and materials science. While MOPAC's Fortran codebase is
not easily approachable to new developers, the intent is that an adjacent ecosystem of Python-based tools is more approachable. Also,
activity in this Python layer may motivate the development of more API-based connections to MOPAC's Fortran-based core functionality in the future.
