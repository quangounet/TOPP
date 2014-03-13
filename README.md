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

Examples, tutorials, reference manual
-------------------------------------

See the wiki https://github.com/quangounet/TOPP/wiki


Notes on OpenRAVE integration
-----------------------------

We will suppose here that you are installing OpenRAVE from source. Let
OPENRAVE_DIR denote your OpenRAVE source folder, for instance:
    
    export OPENRAVE_DIR=~/openrave
    git clone https://github.com/rdiankov/openrave.git OPENRAVE_DIR

Install OpenRAVE (Linux instructions here:
http://openrave.org/docs/latest_stable/coreapihtml/installation_linux.html).
Supposing you kept the default installation path (i.e. /usr/local/), make
a symbolic link:

    /usr/local/include/openrave-0.9/openrave/python -> OPENRAVE_DIR/python

Add the Python bindings folder to your library path by exporting it to
LD_LIBRARY_PATH (you can put the following line in your .bashrc or .zshrc for
persistence):

    export LD_LIBRARY_PATH=$(openrave-config --python-dir)/openravepy/_openravepy_:$LD_LIBRARY_PATH
