Installation
============

Mopactools will soon be available on PyPI and conda-forge, which will handle all aspects of the installation. In particular,
they will package copies of the MOPAC shared library for convenient use within the Python software ecosystem.

Otherwise, mopactools can be installed locally using a standard ``pip install .`` command in the root of the source directory.
For mopactools to function correctly, it must be able to find the MOPAC shared library in the system path or in the ``lib``
subdirectory of the installed package directory. If a copy of the MOPAC shared library is placed into the ``lib`` subdirectory
of the source package, then the ``pip install .`` command will create another copy in the installed package directory.
