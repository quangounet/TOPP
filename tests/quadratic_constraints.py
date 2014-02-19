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

import string
from pylab import *
from numpy import *
from openravepy import *
import TOPPbindings
import TOPPpy
import TOPPopenravepy

# Robot (OpenRAVE)
env = Environment()
env.Load("robots/barrettwam.robot.xml")
env.SetViewer('qtcoin')
env.GetViewer().SetCamera(array([[ 0.92038107,  0.00847738, -0.39093071,  0.69997793],
       [ 0.39101295, -0.02698477,  0.91998951, -1.71402919],
       [-0.00275007, -0.9995999 , -0.02815103,  0.40470174],
       [ 0.        ,  0.        ,  0.        ,  1.        ]]))
robot=env.GetRobots()[0]
grav=[0,0,-9.8]
ndof=robot.GetDOF()
dof_lim=robot.GetDOFLimits()
vel_lim=robot.GetDOFVelocityLimits()
robot.SetDOFLimits(-10*ones(ndof),10*ones(ndof)) # Overrides robot joint limits
robot.SetDOFVelocityLimits(100*vel_lim) # Override robot velocity limits

# Trajectory
q0 = zeros(ndof)
q1 = zeros(ndof)
qd0 = zeros(ndof)
qd1 = zeros(ndof)
q1[0:4] = [2.32883,  1.61082,  0.97706,  1.94169] #Target DOF values for the shoulder and elbow joints
T = 1.5
trajectorystring = "%f\n%d"%(T,ndof)
for i in range(ndof):
    a,b,c,d = TOPPpy.Interpolate3rdDegree(q0[i],q1[i],qd0[i],qd1[i],T)
    trajectorystring += "\n%f %f %f %f"%(d,c,b,a)
traj0 = TOPPpy.PiecewisePolynomialTrajectory.FromString(trajectorystring)

# Constraints
discrtimestep = 0.005
vmax = zeros(ndof)
taumin = zeros(ndof)
taumax = zeros(ndof)
vmax[0:4] = [2,2,2,2]  # Velocity limits, only for the shoulder and elbow joints
taumin[0:4] = [-30,-50,-25,-15] # Torque limits, only for the shoulder and elbow joints
taumax[0:4] = [30,50,25,15]
constraintstring = str(discrtimestep) + "\n";
constraintstring += string.join([str(v) for v in vmax])
constraintstring += TOPPopenravepy.ComputeTorquesConstraints(robot,traj0,taumin,taumax,discrtimestep)

# Run TOPP
x = TOPPbindings.TOPPInstance(None,"QuadraticConstraints",constraintstring,trajectorystring);
ret = x.RunComputeProfiles(0,0)
if(ret == 1):
    x.ReparameterizeTrajectory()

# Display results
ion()    
x.WriteProfilesList()
x.WriteSwitchPointsList()
profileslist = TOPPpy.ProfilesFromString(x.resprofilesliststring)
switchpointslist = TOPPpy.SwitchPointsFromString(x.switchpointsliststring)
TOPPpy.PlotProfiles(profileslist,switchpointslist,4)
if(ret == 1):
    x.WriteResultTrajectory()
    traj1 = TOPPpy.PiecewisePolynomialTrajectory.FromString(x.restrajectorystring)
    dtplot = 0.01
    TOPPpy.PlotKinematics(traj0,traj1,dtplot,vmax)
    TOPPopenravepy.PlotTorques(robot,traj0,traj1,dtplot,taumin,taumax,3)
    print "Trajectory duration before TOPP: ", traj0.duration
    print "Trajectory duration after TOPP: ", traj1.duration
else:
    print "Trajectory is not time-parameterizable"

# Execute trajectory
if(ret == 1):
    TOPPopenravepy.Execute(robot,traj1)
