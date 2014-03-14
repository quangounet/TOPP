TOPP
====

This is TOPP, the Time-Optimal Path Parameterization library by Quang-Cuong
Pham (cuong.pham@normalesup.org)

Requirements 
------------

The following software is required to run TOPP:

- Boost (1.47 or above) with Boost.Python
- Python (2.7 or above)

To integrate with OpenRAVE, you will also need:

- OpenRAVE (0.9 or above) with Python bindings (see "Notes on OpenRAVE integration" below for more details) 
- LAPACK (3.5.0 or above)

Installation
------------

Follow the standard installation procedure: from the TOPP directory,
  
    mkdir build
    cd build
    cmake ..
    make
    sudo make install

TOPP will be compiled with OpenRAVE support if it is found on your system.

See "Notes on OpenRAVE integration" below for more details.

Examples, Tutorials, Reference Manual
-------------------------------------

See the wiki https://github.com/quangounet/TOPP/wiki
