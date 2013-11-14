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
discrtimestep = 0.001
integrationtimestep = discrtimestep
reparamtimestep = 0 #auto
passswitchpointnsteps = 10
tuningsstring = "%f %f %f %d"%(discrtimestep,integrationtimestep,reparamtimestep,passswitchpointnsteps)


############################ Trajectory ############################
#------------------------------------------#

# TODO: integration autour du point singulier de gauche
sdbeg_min, sdbeg_max = 0.000000, 5.143234
trajectorystring = """0.727992
2
-0.0539850761315 0.744141496569 -1.85622830904
-0.300968741716 0.668022030388 -2.0090717382"""

# TODO: alpha au-dessus de beta sous la MVC
sdbeg_min, sdbeg_max = 0.000000, 8.808698
trajectorystring = """1.000000
2
-0.420982931773 -0.111291845962 1.77275676137 -1.33551682399
-0.630425778805 0.0894067072732 0.0852189728776 0.290595287026"""



#------------------------------------------#
traj0 = TOPPpy.PiecewisePolynomialTrajectory.FromString(trajectorystring)


############################ Constraints ############################
#------------------------------------------#
vmax = array([0,0])
taumin = array([-11,-7])
taumax = array([11,7])
t0 = time.time()
constraintstring = string.join([str(x) for x in taumin]) + "\n" + string.join([str(a) for a in taumax]) + "\n" + string.join([str(a) for a in vmax])
constraintstring += TOPPopenravepy.ComputeTorquesConstraintsLegacy(robot,traj0,taumin,taumax,discrtimestep)
#------------------------------------------#


############################ Run TOPP ############################
t1 = time.time()
x = TOPPbindings.TOPPInstance("TorqueLimits",constraintstring,trajectorystring,tuningsstring)
t2 = time.time()
#ret = x.RunComputeProfiles(0,0)
ret = x.RunVIP(sdbeg_min, sdbeg_max)
print ret
t3 = time.time()

#print x.resduration
print "sdendmin =", x.sdendmin
print "sdendmax =", x.sdendmax

# if(ret == 1):
#     x.ReparameterizeTrajectory()

t4 = time.time()

################ Plotting the MVC and the profiles #################
x.WriteProfilesList()
x.WriteSwitchPointsList()
profileslist = TOPPpy.ProfilesFromString(x.resprofilesliststring)
switchpointslist = TOPPpy.SwitchPointsFromString(x.switchpointsliststring)

def replot():
    clf()
    TOPPpy.PlotProfiles(profileslist,switchpointslist,1)
    axis([0, traj0.duration, 0, 100])

replot()
#TOPPpy.PlotAlphaBeta(x, axis())



##################### Plotting the trajectories #####################
# if(ret == 1):
#     x.WriteResultTrajectory()
#     traj1 = TOPPpy.PiecewisePolynomialTrajectory.FromString(x.restrajectorystring)
#     dtplot = 0.01
#     TOPPpy.PlotKinematics(traj0,traj1,dtplot,vmax)
#     TOPPopenravepy.PlotTorques(robot,traj0,traj1,dtplot,taumin,taumax,3)


# print "\n--------------"
# print "Python preprocessing: ", t1-t0
# print "Building TOPP Instance: ", t2-t1
# print "Compute profiles: ", t3-t2
# print "Reparameterize trajectory: ", t4-t3
# print "Total: ", t4-t0
# if(ret == 1):
#     print "Trajectory duration (estimate): ", x.resduration
#     print "Trajectory duration: ", traj1.duration

raw_input()
