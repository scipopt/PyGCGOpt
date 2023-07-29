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
8.0 | 4.0 | 3.5.x | 0.1.x |
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

## Installing SCIPOptSuite from the binary distribution

The SCIPOptSuite can either be installed from the binary distribution or it can be build and installed manually from source. If you want to compile the SCIPOptSuite from source, please skip to the next section.

To install a binary distribution of the SCIPOptSuite, navigate to the [download page](https://scipopt.org/index.php#download) and select at least version 8.0.0 and your operating system. After the installation, please check that you can run `scip` and `gcg` from the command line.

Next, skip to the instructions to install PySCIPOpt.

## Compiling SCIPOptSuite from source

To install the SCIPOptSuite from source, you need to [install SCIP and GCG using CMake](https://scipopt.org/doc/html/md_INSTALL.php#CMAKE). The Makefile system is not compatible with PyGCGOpt and PySCIPOpt!

As an example, from the root directory of the source distribution of the SCIPOptSuite, you can run the following commands to build and install locally:
```
mkdir build && cd build
cmake .. -DCMAKE_BUILD_TYPE=Debug -DSHARED=on -DCMAKE_INSTALL_PREFIX=./install -DGCG_DEV_BUILD=ON -DZIMPL=OFF -DIPOPT=OFF -DPAPILO=OFF
make -j4 && make install
```

If SCIP and GCG are not installed in the global path, you need to specify the install location using the environment variable `SCIPOPTDIR` before you can install PySCIPOpt and PyGCGOpt:

 - On Linux and macOS:\
   `export SCIPOPTDIR=<path_to_install_dir>`
 - On Windows:\
   `set SCIPOPTDIR=<path_to_install_dir>` (**cmd**, **Cmder**, **WSL**)\
   `$Env:SCIPOPTDIR = "<path_to_install_dir>"` (**powershell**)
  
`SCIPOPTDIR` needs to have a subdirectory `lib` that contains the
library, e.g. `libscip.so` and `libgcg.so` (for Linux) and a subdirectory `include` that
contains the corresponding header files:

    SCIPOPTDIR
      > lib
        > libscip.so ...
        > libgcg.so ...
      > include
        > scip
        > gcg
        > lpi
        > nlpi
        > ...

If GCG is installed in a different location than SCIP (this is *not* the default), you can specify the GCG install directory with the `GCGOPTDIR` environment variable.

To make the shared libraries available at runtime for PySCIPOpt and PyGCGOpt, you need to add them to the corresponding environment variables:
 - On Linux:
   `export LD_LIBRARY_PATH=<path_to_install_dir>/lib`
 - On macOS:
   `export DYLD_FALLBACK_LIBRARY_PATH=<path_to_install_dir>/lib`
 - On Windows: 
   `set PATH=%PATH%;<path_to_install_dir>\bin`


## Installing PySCIPOpt from PyPI

To install PySCIPOpt using pip, run the following:
```
pip install PySCIPOpt>=4.0.0
```


## Installing PySCIPOpt from source

For detailed instructions, please read the [installation instructions of PySCIPOpt](https://github.com/scipopt/PySCIPOpt/blob/master/INSTALL.md#building-everything-from-source).

In short, you can install PySCIPOpt from source with `pip`. For that, run the following command:
```
pip install <path/to/PySCIPOpt>
```

## Installing PyGCGOpt from PyPI

Simply run the folloing to install PyGCGOpt with `pip`:
```
pip install PyGCGOpt
```

## Install PyGCGOpt from source
-------------------------------

First, install the requirements for building PyGCGOpt:
```
pip install -r requirements.txt
```

Furthermore, you need to have the Python development files installed on your system (error message "Python.h not found"):
```
sudo apt-get install python3-dev  # for Python 3, on Linux
```

Then, compile and install PyGCGOpt:
```
pip install .
```

In case the compiler cannot find SCIP or GCG header files, please make sure that you installed the SCIPOptSuite globally or that you set the `SCIPOPTDIR` environment variable (see above).

### Building with debug information

To use debug information in PyGCGOpt you need to build it like this:

    python -m pip install --install-option="--debug" .

Be aware that you will need the **debug library** of the SCIPOptSuite for this to work
(`cmake .. -DCMAKE_BUILD_TYPE=Debug`).

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
