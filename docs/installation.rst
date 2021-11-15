Installation
============

*******************************
Building everything from source
*******************************

The interface is developed on Linux. The following instructions should, in theory, also work on macOS. It uses the CMake build system. The Makefile system should work in theory but has some unresolved dependency issues with shared libraries as of writing this.

First, we need to build GCG with custom parameters:

1. Navigate to your local gcg folder and make sure you have the current master branch checked out.
2. Create a subfolder for building: :code:`mkdir build && cd build`
3. Run CMake and setup building the shared library:

.. code:: bash

    cmake .. -DCMAKE_BUILD_TYPE=Release -DSHARED=on \
             -DCMAKE_INSTALL_PREFIX=./install -DZIMPL=OFF -DIPOPT=OFF -DPAPILO=OFF


4. Compile and install GCG: :code:`make && make install`

Now, we can build PyGCGOpt.

1. Navigate to the local PyGCGOpt folder.
2. Set environment variables to point to GCG:

.. code:: bash

    export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/path/to/gcg/build/install/lib
    export GCGOPTDIR=/path/to/gcg/build/install
    export SCIPOPTDIR=/path/to/gcg/build/install

3. Install cython: :code:`pip install cython`
4. (Optional) PyGCGOpt depends on PySCIPOpt during compile time and run time. Before you can compile PyGCGOpt you need to install PySCIPOpt from the package respository (see above) or have a local copy checked out. To install the local copy, execute :code:`pip install <path/to/PySCIPOpt>`.
5. Build the extension module: :code:`pip install .`

Now, the module can be imported in Python using :code:`import pygcgopt`.