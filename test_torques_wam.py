# -*- coding: utf-8 -*-
# Copyright (C) 2013 Quang-Cuong Pham <cuong.pham@normalesup.org>
#
# This file is part of the Time-Optimal Path Parameterization (TOPP) library.
# TOPP is free software: you can redistribute it and/or modify
# it under the terms of the GNU Lesser General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU Lesser General Public License for more details.
#
# You should have received a copy of the GNU Lesser General Public License
# along with this program.  If not, see <http://www.gnu.org/licenses/>.


import TOPPbindings
import TOPPpy
import TOPPopenravepy
import time
import string
import sys
from pylab import *
from numpy import *
from openravepy import *



ion()

########################### Robot ################################
env = Environment() # create openrave environment
#------------------------------------------#
robotfile = "robots/arm.robot.xml"
env.Load(robotfile)
robot=env.GetRobots()[0]
#------------------------------------------#
grav=[0,0,-9.8]
n=robot.GetDOF()
dof_lim=robot.GetDOFLimits()
vel_lim=robot.GetDOFVelocityLimits()
robot.SetDOFLimits(-10*ones(n),10*ones(n))
robot.SetDOFVelocityLimits(100*vel_lim)


############################ Tunings ############################
discrtimestep = 0.01;
integrationtimestep = 0.01;
bisectionprecision = 0.01;
passswitchpointnsteps = 5;
reparamtimestep = 0.01;
tuningsstring = "%f %f %f %d %f"%(discrtimestep,integrationtimestep,bisectionprecision,passswitchpointnsteps,reparamtimestep);


############################ Trajectory ############################
#------------------------------------------#
q0=[0,0,0,0]
q1=[ 2.32883,  1.61082,  0.97706,  1.94169]
v=1
qd0=[v,v,v,v]
qd1=[v,v,v,v]
T = 1.5
trajectorystring = "%f\n%d"%(T,4)
for i in range(4):
    a,b,c,d = TOPPpy.Interpolate3rdDegree(q0[i],q1[i],qd0[i],qd1[i],T)
    trajectorystring += "\n%f %f %f %f"%(d,c,b,a)
#------------------------------------------#
traj0 = TOPPpy.PiecewisePolynomialTrajectory.FromString(trajectorystring)


############################ Constraints ############################
#------------------------------------------#
taumin = array([-6,-15,-5,-4])
taumax = array([6,15,5,4])
vmax = [3,3,3,3]
t0 = time.time()
constraintstring = string.join([str(v) for v in vmax])
constraintstring += TOPPopenravepy.ComputeTorquesConstraints(robot,traj0,taumin,taumax,discrtimestep)
#------------------------------------------#


############################ Run TOPP ############################
t1 = time.time()
x = TOPPbindings.TOPPInstance("QuadraticConstraints",constraintstring,trajectorystring,tuningsstring);
t2 = time.time()
ret = x.RunComputeProfiles(0,0)
t3 = time.time()

if(ret == 1):
    x.ReparameterizeTrajectory()

t4 = time.time()

################ Plotting the MVC and the profiles #################
x.WriteProfilesList()
x.WriteSwitchPointsList()
profileslist = TOPPpy.ProfilesFromString(x.resprofilesliststring)
switchpointslist = TOPPpy.SwitchPointsFromString(x.switchpointsliststring)
TOPPpy.PlotProfiles(profileslist,switchpointslist,4)


##################### Plotting the trajectories #####################
if(ret == 1):
    x.WriteResultTrajectory()
    traj1 = TOPPpy.PiecewisePolynomialTrajectory.FromString(x.restrajectorystring)
    dtplot = 0.01
    TOPPpy.PlotKinematics(traj0,traj1,dtplot,vmax)


print "\n--------------"
print "Building TOPP Instance (including sampling dynamics in C++): ", t2-t1
print "Compute profiles (C++): ", t3-t2
print "Reparameterize trajectory (C++): ", t4-t3
print "Total: ", t4-t1 
print "Trajectory duration (estimate): ", x.resduration
print "Trajectory duration: ", traj1.duration

raw_input()
