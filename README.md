PyGCGOpt
=========

This project provides an interface from Python to the [GCG Solver](https://gcg.or.rwth-aachen.de/).

Documentation
-------------

See [CHANGELOG.md](CHANGELOG.md) for added, removed or fixed functionality.

Installation
------------

**Using Conda**

[![Conda version](https://img.shields.io/conda/vn/conda-forge/pygcgopt?logo=conda-forge)](https://anaconda.org/conda-forge/pygcgopt)
[![Conda platforms](https://img.shields.io/conda/pn/conda-forge/pygcgopt?logo=conda-forge)](https://anaconda.org/conda-forge/pygcgopt)

Conda will install `GCG`, `SCIP` and `PySCIPOpt` automatically, hence everything can be installed in a single command:
```bash
conda install --channel conda-forge pygcgopt
```

**Using PyPI and from Source**

See [INSTALL.md](INSTALL.md) for instructions.
Please note that the latest PyGCGOpt version is usually only compatible with the latest major release of the SCIP Optimization Suite and the GCG Solver.
Information which version of PyGCGOpt is required for a given GCG version can also be found in [INSTALL.md](INSTALL.md).

