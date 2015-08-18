TOPP
====

This is TOPP, the Time-Optimal Path Parameterization library by Quang-Cuong Pham (cuong.pham@normalesup.org)

If you use this library for your research, please reference the accompanying paper « A general, fast, and robust implementation of the time-optimal path parameterization algorithm », IEEE Transactions on Robotics, vol. 30(6), pp. 1533-1540, 2014.

Many thanks to Stephane Caron, Rosen Diankov and Puttichai Lertkultanon for their contributions !

Requirements 
------------

The following software is required to install TOPP:

- Python (2.7 or above), with numpy, scipy and matplotlib
- Boost (1.46 or above), with Boost.Python

If you need OpenRAVE support (for dynamics computations), the following software is also required:

- OpenRAVE (installed from source, "master" branch, see http://openrave.org/docs/latest_stable/coreapihtml/installation_linux.html **NB**: one must install from the "master" branch, not the "latest stable" branch, so change <code>git clone --branch latest_stable https://github.com/rdiankov/openrave.git</code> into <code>git clone --branch master https://github.com/rdiankov/openrave.git</code>)

Installation
------------

Follow the standard installation procedure: from the TOPP directory,
  
    mkdir build
    cd build
    cmake ..
    make
    sudo make install

TOPP will be compiled with OpenRAVE support if the latter is found on your system.

Examples
--------

Please try the test files in the tests/ folder or copy-paste test cases from https://github.com/quangounet/TOPP/wiki/Quick-examples

Documentation, Tutorials...
---------------------------

See the wiki https://github.com/quangounet/TOPP/wiki
