# Installation


## Using Conda

[![Conda version](https://img.shields.io/conda/vn/conda-forge/pygcgopt?logo=conda-forge)](https://anaconda.org/conda-forge/pygcgopt)
[![Conda platforms](https://img.shields.io/conda/pn/conda-forge/pygcgopt?logo=conda-forge)](https://anaconda.org/conda-forge/pygcgopt)

Conda will install `GCG`, `SCIP` and `PySCIPOpt` automatically, hence everything can be installed in a single command:
```bash
conda install --channel conda-forge pygcgopt
```

## Requirements

PyGCGOpt requires a working installation of the [SCIPOptSuite](https://scipopt.org) and the [GCG Solver](https://gcg.or.rwth-aachen.de/) which is usually included in the release. In addition, the Python interface for SCIP, [PySCIPOpt](https://github.com/scipopt/PySCIPOpt), has to be installed. Currently, GCG only runs on Linux and macOS, therefore you can use PyGCGOpt only on these operating systems.

Note that the latest PyGCGOpt version is usually only compatible with the latest major release of the SCIP Optimization Suite. PyGCGOpt requires at least SCIP 8.0 and GCG 3.5. The following table summarizes which versions of PyGCGOpt, GCG, PySCIPOpt, and SCIP are compatible:

|SCIP| PySCIPOpt | GCG | PyGCGOpt
|----|----|----|----|
9.0 | 5.0.0 | 3.6.0 | 0.4.0 |
8.1 | 4.4 | 3.5.5 | 0.3.0 |
8.0 | 4.0 | 3.5.0 | 0.1.0 |
7.0 | 3.x | - | - |
6.0 | 2.x | - | - |
5.0 | 1.4, 1.3 | - | - |
4.0 | 1.2, 1.1 | - | - |
3.2 | 1.0 | - | - |

We recommend installing PyGCGOpt and its Python dependencies in a [Python virtual environment](https://docs.python.org/3/tutorial/venv.html). In short, you can create a virtual environment and activate it with the following commands on Linux, macOS and Windows:
```
python3.9 -m venv venv
source venv/bin/activate
```

## Installing PyGCGOpt from PyPI

Simply run the following to install PyGCGOpt (and its requirement PySCIPOpt) with `pip`:
```
pip install PySCIPOpt
pip install PyGCGOpt
```

The following table summarizes which versions of PyGCGOpt and PySCIPOpt are compatible:

| PySCIPOpt | PyGCGOpt |
|----|----|
| 5.2.1 | 0.4.0 |


## Testing new installation

To test your brand-new installation of PyGCGOpt you need
[pytest](https://pytest.org/) on your system.

    pip install pytest

Here is the complete [installation
procedure](https://docs.pytest.org/en/latest/getting-started.html).

Tests can be run in the `PyGCGOpt` directory with: :

    py.test # all the available tests
    py.test tests/test_name.py # a specific tests/test_name.py (Unix)

Ideally, the status of your tests must be passed or skipped. Running
tests with pytest creates the `__pycache__` directory in `tests` and,
occasionally, a `model` file in the working directory. They can be
removed harmlessly.
