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

import string
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
env.Load("robots/twodof.robot.xml")
env.SetViewer('qtcoin')
env.GetViewer().SetCamera(array([[ 0.00846067,  0.4334184 , -0.9011531 ,  0.84555054],
       [ 0.99938498,  0.0270039 ,  0.02237072,  0.01155015],
       [ 0.03403053, -0.90078814, -0.43292337,  0.65048862],
       [ 0.        ,  0.        ,  0.        ,  1.        ]]))
robot = env.GetRobots()[0]
robot.SetTransform(array([[0, 0, 1, 0],
                          [0, 1, 0, 0],
                          [-1, 0, 0, 0.3],
                          [0, 0, 0, 1]]))
grav = [0, 0, -9.8]
ndof = robot.GetDOF()
dof_lim = robot.GetDOFLimits()
vel_lim = robot.GetDOFVelocityLimits()
robot.SetDOFLimits(-10 * ones(ndof), 10 * ones(ndof)) # Overrides robot joint limits
robot.SetDOFVelocityLimits(100 * vel_lim) # Overrides robot velocity limits

# Trajectory
q0 = [0,0]
q1 = [5*pi/4,-pi/2]
qd0 = [1,1]
qd1 = [1,1]
T = 2
trajectorystring = "%f\n%d"%(T,ndof)
for i in range(ndof):
    a,b,c,d = Utilities.Interpolate3rdDegree(q0[i],q1[i],qd0[i],qd1[i],T)
    trajectorystring += "\n%f %f %f %f"%(d,c,b,a)
traj0 = Trajectory.PiecewisePolynomialTrajectory.FromString(trajectorystring)

# Constraints
discrtimestep = 0.002
vmax = array([5, 5])
taumin = array([-25, -10])
taumax = array([25, 10])
constraintstring = str(discrtimestep) + "\n";
constraintstring += string.join([str(v) for v in vmax]) + "\n"
constraintstring += string.join([str(t) for t in taumin]) + "\n" + string.join([str(t) for t in taumax])

#Run TOPP
x = TOPPbindings.TOPPInstance(robot,"TorqueLimitsRave", constraintstring, trajectorystring)
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
