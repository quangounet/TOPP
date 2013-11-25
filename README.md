TOPP
====

This is TOPP, the Time-Optimal Path Parameterization library by Quang-Cuong
Pham (cuong.pham@normalesup.org)

Requirements 
------------

- Boost library
- Python
- OpenRAVE (for TOPP with Torque Limits only)

Installation
------------

From the top folder:
  
    make release

will compile TOPPbindings.so in the current folder.

**OpenRAVE integration**

Let OPENRAVE_DIR denote your OpenRAVE source folder, for instance:
    
    export OPENRAVE_DIR=~/openrave
    git clone https://github.com/rdiankov/openrave.git OPENRAVE_DIR

Install OpenRAVE (Linux instructions here:
http://openrave.org/docs/latest_stable/coreapihtml/installation_linux.html).
Supposing you kept the default installation path (i.e. /usr/local/), make
a symbolic link

    /usr/local/include/openrave-0.9/openrave/python -> OPENRAVE_DIR/python

Add the Python bindings folder to your library path by exporting it to
LD_LIBRARY_PATH (you can put the following line in your .bashrc or .zshrc for
persistence):

    export LD_LIBRARY_PATH=$(openrave-config --python-dir)/openravepy/_openravepy_:$LD_LIBRARY_PATH

You will also need to patch file
/home/cuong/git/openrave/python/bindings/bindings.h (ask Stephane) to support
compilation with C++11 (-std=c++0x).

Testing
-------

From TOPP's top directory:
  
    python test_kinematic.py  # test with kinematic constraints
    python test_torques.py    # test with torque limits
