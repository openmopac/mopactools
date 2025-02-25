Installation
============

Mopactools is available on PyPI and can be installed using a standard ``pip install mopactools``. The PyPI distribution
includes a copy of the MOPAC shared library. Mopactools will also be available on the conda-forge channel of the Conda package
manager in the near future.

Otherwise, mopactools can be installed locally from a clone of its GitHub repository using a standard ``pip install .`` command
in the root of the source directory. For mopactools to function correctly, it must be able to find the MOPAC shared library in
the system path or in the ``lib`` subdirectory of the installed package directory. If a copy of the MOPAC shared library is
placed into the ``lib`` subdirectory of the source package, then the ``pip install .`` command will create another copy in the
installed package directory.
