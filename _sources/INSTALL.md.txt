# Installation


## Using Conda

[![Conda version](https://img.shields.io/conda/vn/conda-forge/pygcgopt?logo=conda-forge)](https://anaconda.org/conda-forge/pygcgopt)
[![Conda platforms](https://img.shields.io/conda/pn/conda-forge/pygcgopt?logo=conda-forge)](https://anaconda.org/conda-forge/pygcgopt)

Conda will install `GCG`, `SCIP` and `PySCIPOpt` automatically, hence everything can be installed in a single command:
```bash
conda install --channel conda-forge pygcgopt
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
| 5.3.0 | 0.4.1 |


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
