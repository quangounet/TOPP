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
robotfile = "robots/twodof.robot.xml"
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
discrtimestep = 0.01;
integrationtimestep = 0.01;
bisectionprecision = 0.01;
passswitchpointnsteps = 5;
reparamtimestep = 0.01;
tuningsstring = "%f %f %f %d %f"%(discrtimestep,integrationtimestep,bisectionprecision,passswitchpointnsteps,reparamtimestep);


############################ Trajectory ############################
#------------------------------------------#
T=1
[a1,b1,c1,a2,b2,c2] =  [3, -3, -3, 0, -2, -2] #[-3, 3, 3, -1, 0, -3]
trajectorystring = "%f\n%d\n%f %f %f\n%f %f %f"%(T,2,c1,b1,a1,c2,b2,a2)
#------------------------------------------#
traj0 = TOPPpy.PiecewisePolynomialTrajectory(trajectorystring)


############################ Constraints ############################
#------------------------------------------#
taumin = array([-15,-10])
taumax = array([15,10])
vmax = array([2.5,3])
constraintstring = string.join([str(x) for x in taumin]) + "\n" + string.join([str(a) for a in taumax]) + "\n" + string.join([str(a) for a in vmax])
constraintstring += TOPPopenravepy.ComputeConstraintsLegacy(robot,traj0,taumin,taumax,discrtimestep)
#------------------------------------------#


############################ Run TOPP ############################
x = TOPPbindings.TOPPInstance("TorqueLimits",constraintstring,trajectorystring,tuningsstring)
ret = x.RunPP(1e-4,1e-4)


################ Plotting the MVC and the profiles #################
x.WriteProfilesList()
profileslist = TOPPpy.ProfilesFromString(x.resprofilesliststring)
TOPPpy.PlotProfiles(profileslist,4)


##################### Plotting the trajectories #####################
x.WriteResultTrajectory()
traj1 = TOPPpy.PiecewisePolynomialTrajectory(x.restrajectorystring)
dtplot = 0.01
TOPPpy.PlotKinematics(traj0,traj1,dtplot,vmax)
TOPPopenravepy.PlotTorques(robot,traj0,traj1,dtplot,taumin,taumax,3)



raw_input()
