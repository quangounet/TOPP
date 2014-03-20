TOPP
====

This is TOPP, the Time-Optimal Path Parameterization library by Quang-Cuong
Pham (cuong.pham@normalesup.org)

Requirements 
------------

The following software is required to install TOPP:

- Boost (1.47 or above) with Boost.Python
- Python (2.7 or above)

If you need OpenRAVE support (for dynamics computations), the following software is also required:

- OpenRAVE (0.9 or above) with Python bindings
- LAPACK (3.5.0 or above)

Installation
------------

Follow the standard installation procedure: from the TOPP directory,
  
    mkdir build
    cd build
    cmake ..
    make
    sudo make install

TOPP will be compiled with OpenRAVE support if the latter is found on your system.

Examples, Tutorials, Reference Manual
-------------------------------------

See the wiki https://github.com/quangounet/TOPP/wiki
