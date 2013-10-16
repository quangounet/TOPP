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
discrtimestep = 0.01
integrationtimestep = 0.01
reparamtimestep = 0.01
passswitchpointnsteps = 5
tuningsstring = "%f %f %f %d"%(discrtimestep,integrationtimestep,reparamtimestep,passswitchpointnsteps)


############################ Trajectory ############################
#------------------------------------------#
T=1
[a1,b1,c1,a2,b2,c2] =  [3, -3, -3, 0, -2, -2] #[-3, 3, 3, -1, 0, -3]
trajectorystring = "%f\n%d\n%f %f %f\n%f %f %f"%(T,2,c1,b1,a1,c2,b2,a2)
#------------------------------------------#
traj0 = TOPPpy.PiecewisePolynomialTrajectory.FromString(trajectorystring)


############################ Constraints ############################
#------------------------------------------#
taumin = array([-13,-7])
taumax = array([13,7])
vmax = array([0,0])
t0 = time.time()
constraintstring = string.join([str(x) for x in taumin]) + "\n" + string.join([str(a) for a in taumax]) + "\n" + string.join([str(a) for a in vmax])
constraintstring += TOPPopenravepy.ComputeTorquesConstraintsLegacy(robot,traj0,taumin,taumax,discrtimestep)
#------------------------------------------#


############################ Run TOPP ############################
t1 = time.time()
x = TOPPbindings.TOPPInstance("TorqueLimits",constraintstring,trajectorystring,tuningsstring)
t2 = time.time()
x.RunVIP(1,4)
t3 = time.time()


################ Plotting the MVC and the profiles #################
x.WriteProfilesList()
x.WriteSwitchPointsList()
profileslist = TOPPpy.ProfilesFromString(x.resprofilesliststring)
switchpointslist = TOPPpy.SwitchPointsFromString(x.switchpointsliststring)
TOPPpy.PlotProfiles(profileslist,switchpointslist,4)


print "\n--------------"
print "Python preprocessing: ", t1-t0
print "Building TOPP Instance: ", t2-t1
print "Compute profiles: ", t3-t2
print "Total: ", t3-t0
print "(sdendmin,sdendmax) = (", x.sdendmin, ",", x.sdendmax, ")"

raw_input()
