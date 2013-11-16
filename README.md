TOPP
====

This is TOPP, the Time-Optimal Path Parameterization library by Quang-Cuong Pham (cuong.pham@normalesup.org)

Requirements 
------------

- Boost library
- Python
- OpenRAVE (for TOPP with Torque Limits only)

Installation
------------

I) OpenRAVE stuffs
Assume that your openrave folder is /home/cuong/git/openrave/

1) In the TOPP folder, create a link to the openrave folder
ln -s /home/cuong/git/openrave .

2) In .bashrc add the following line
export LD_LIBRARY_PATH=/home/cuong/git/openrave/build/python/bindings:$LD_LIBRARY_PATH

3) In /home/cuong/git/openrave/python/bindings/openravepy_int.h, comment out the line
#define PY_ARRAY_UNIQUE_SYMBOL PyArrayHandle

4) Patch file /home/cuong/git/openrave/python/bindings/bindings.h (ask Stephane)




From the top folder:
  
    make

For TOPP with Kinematic Limits:
  
    python test_kinematic.py

For TOPP with Torque Limits:

    python test_torques.py
