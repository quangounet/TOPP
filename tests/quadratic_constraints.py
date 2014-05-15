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

print "\n************************************\nNB: This test file requires OpenRAVE\n************************************\n"

import string,time
from pylab import *
from numpy import *
from openravepy import *
from TOPP import TOPPbindings
from TOPP import TOPPpy
from TOPP import TOPPopenravepy
from TOPP import Trajectory
from TOPP import Utilities

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
q0[0:7] = [-2,-1,-1.5,-2.5,0,0,0]
q1[0:7] = [2,-3,1.5,2,1,1,1]
T = 1.5
trajectorystring = "%f\n%d"%(T,ndof)
for i in range(ndof):
    a,b,c,d = Utilities.Interpolate3rdDegree(q0[i],q1[i],qd0[i],qd1[i],T)
    trajectorystring += "\n%f %f %f %f"%(d,c,b,a)
traj0 = Trajectory.PiecewisePolynomialTrajectory.FromString(trajectorystring)

t0 = time.time()
# Constraints
discrtimestep = 0.001
vmax = zeros(ndof)
taumin = zeros(ndof)
taumax = zeros(ndof)
vmax[0:7] = 1.5*vel_lim[0:7]  # Velocity limits
taumin[0:7] = -robot.GetDOFMaxTorque()[0:7] # Torque limits
taumax[0:7] = robot.GetDOFMaxTorque()[0:7] # Torque limits
constraintstring = str(discrtimestep) + "\n";
constraintstring += string.join([str(v) for v in vmax])
constraintstring += TOPPopenravepy.ComputeTorquesConstraints(robot,traj0,taumin,taumax,discrtimestep)

t1 = time.time()

# Run TOPP
x = TOPPbindings.TOPPInstance(None,"QuadraticConstraints",constraintstring,trajectorystring);

# Testing constraints serialization
x.WriteConstraints()
x = TOPPbindings.TOPPInstance(None,"QuadraticConstraints",x.outconstraintstring,trajectorystring);

ret = x.RunComputeProfiles(0,0)
if(ret == 1):
    x.ReparameterizeTrajectory()

print  t1-t0, time.time()-t1

# Display results
ion()
x.WriteProfilesList()
x.WriteSwitchPointsList()
profileslist = TOPPpy.ProfilesFromString(x.resprofilesliststring)
switchpointslist = TOPPpy.SwitchPointsFromString(x.switchpointsliststring)
TOPPpy.PlotProfiles(profileslist,switchpointslist,4)
if(ret == 1):
    x.WriteResultTrajectory()
    traj1 = Trajectory.PiecewisePolynomialTrajectory.FromString(x.restrajectorystring)
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

raw_input()
