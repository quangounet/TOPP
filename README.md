TOPP
====

This is TOPP, the Time-Optimal Path Parameterization library by Quang-Cuong
Pham (cuong.pham@normalesup.org)

If you use this library for your research, please reference the accompanying paper « A general, fast, and robust implementation of the time-optimal path parameterization algorithm » http://arxiv.org/abs/1312.6533 



Requirements 
------------

The following software is required to install TOPP:

- Boost (1.46 or above) with Boost.Python
- Python (2.7 or above)

If you need OpenRAVE support (for dynamics computations), the following software is also required:

- OpenRAVE (installed from source, "master" branch, see http://openrave.org/docs/latest_stable/coreapihtml/installation_linux.html **NB**: one must install from the "master" branch, not the "latest stable" branch, so change <code>git clone --branch latest_stable https://github.com/rdiankov/openrave.git</code> into <code>git clone --branch master https://github.com/rdiankov/openrave.git</code>)
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
