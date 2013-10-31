# -*- coding: utf-8 -*-
# Copyright (C) 2013 Quang-Cuong Pham <cuong.pham@normalesup.org>
#
# This file is part of the Time-Optimal Path Parameterization (TOPP) library.
# TOPP is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <http://www.gnu.org/licenses/>.

import sys
sys.path.append('..')

import TOPPbindings
import TOPPpy
import TOPPopenravepy
import time
import string
from pylab import *
from numpy import *
from openravepy import *


ion()

########################### Robot ################################
env = Environment() # create openrave environment
#------------------------------------------#
robotfile = "../robots/twodof.robot.xml"
env.Load(robotfile)
robot=env.GetRobots()[0]
robot.SetTransform(array([[0,0,1,0],[0,1,0,0],[-1,0,0,0.3],[0,0,0,1]]))
#------------------------------------------#
grav=[0,0,-9.8]
n=robot.GetDOF()
dof_lim=robot.GetDOFLimits()
vel_lim=robot.GetDOFVelocityLimits()
robot.SetDOFLimits(-10*ones(n),10*ones(n))
robot.SetDOFVelocityLimits(100*vel_lim)


############################ Tunings ############################
discrtimestep = 0.005
integrationtimestep = discrtimestep
reparamtimestep = 0 #auto
passswitchpointnsteps = 5
tuningsstring = "%f %f %f %d"%(discrtimestep,integrationtimestep,reparamtimestep,passswitchpointnsteps)


############################ Trajectory ############################
#------------------------------------------#
T=1
[a1,b1,c1,a2,b2,c2] =  [3, -3, -3, 0, -2, -2] #[-3, 3, 3, -1, 0, -3]
trajectorystring = "%f\n%d\n%f %f %f\n%f %f %f"%(T,2,c1,b1,a1,c2,b2,a2)
trajectorystring =  """0.985516
2
0.0 -0.996909732549
0.0 -0.0785556181874"""

#------------------------------------------#
traj0 = TOPPpy.PiecewisePolynomialTrajectory.FromString(trajectorystring)


############################ Constraints ############################
#------------------------------------------#
#taumin = array([-15,-10])
#taumax = array([15,10])
vmax = array([0,0])
taumin = array([-11,-7])
taumax = array([11,7])
#vmax = array([0,0])
t0 = time.time()
constraintstring = string.join([str(x) for x in taumin]) + "\n" + string.join([str(a) for a in taumax]) + "\n" + string.join([str(a) for a in vmax])
print str(traj0)
constraintstring += TOPPopenravepy.ComputeTorquesConstraintsLegacy(robot,traj0,taumin,taumax,discrtimestep)
print "len(cstring) =", len(constraintstring)

#------------------------------------------#


############################ Run TOPP ############################
t1 = time.time()
x = TOPPbindings.TOPPInstance("TorqueLimits",constraintstring,trajectorystring,tuningsstring)
t2 = time.time()
#ret = x.RunComputeProfiles(0,0)
ret = x.RunVIP(0,1e-4)
print ret
t3 = time.time()
print x.sdendmin
print x.sdendmax

if(ret == 1) and False:
    x.ReparameterizeTrajectory()

t4 = time.time()

################ Plotting the MVC and the profiles #################
x.WriteProfilesList()
x.WriteSwitchPointsList()
profileslist = TOPPpy.ProfilesFromString(x.resprofilesliststring)
switchpointslist = TOPPpy.SwitchPointsFromString(x.switchpointsliststring)
TOPPpy.PlotProfiles(profileslist,switchpointslist,4)
axis([0,1,0,100])

##################### Plotting the trajectories #####################
if(ret == 1) and False:
    x.WriteResultTrajectory()
    traj1 = TOPPpy.PiecewisePolynomialTrajectory.FromString(x.restrajectorystring)
    dtplot = 0.01
    TOPPpy.PlotKinematics(traj0,traj1,dtplot,vmax)
    TOPPopenravepy.PlotTorques(robot,traj0,traj1,dtplot,taumin,taumax,3)


print "\n--------------"
print "Python preprocessing: ", t1-t0
print "Building TOPP Instance: ", t2-t1
print "Compute profiles: ", t3-t2
print "Reparameterize trajectory: ", t4-t3
print "Total: ", t4-t0
if(ret == 1) and False:
    print "Trajectory duration (estimate): ", x.resduration
    print "Trajectory duration: ", traj1.duration

raw_input()
