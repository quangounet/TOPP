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

import string
from pylab import *
from numpy import *
from openravepy import *
from TOPP import TOPPbindings
from TOPP import TOPPpy
from TOPP import TOPPopenravepy
from TOPP import Trajectory

# Robot
env = Environment()
env.SetViewer('qtcoin')
env.Load('data/lab1.env.xml')
robot = env.GetRobots()[0]
RaveSetDebugLevel(DebugLevel.Debug)
ndof = robot.GetDOF()

# Finding a trajectory using OpenRAVE RRT with constraintparabolicsmoothing
vmax = 2 * ones(ndof)
amax = 10* ones(ndof)
robot.SetDOFVelocityLimits(vmax)
robot.SetDOFAccelerationLimits(amax)
robot.SetDOFValues(zeros(ndof))
robot.SetActiveDOFs(range(4)) # set the shoulder and elbow joints to be the only active joints
params = Planner.PlannerParameters()
params.SetRobotActiveJoints(robot)
params.SetGoalConfig([0.1,pi/2,pi/3,pi/2.1])
params.SetExtraParameters("""<_postprocessing planner="constraintparabolicsmoother">
    <_fStepLength>0</_fStepLength>
    <minswitchtime>0.5</minswitchtime>
    <_nmaxiterations>40</_nmaxiterations>
</_postprocessing>""")
planner=RaveCreatePlanner(env,'birrt')
planner.InitPlan(robot, params)
ravetraj0 = RaveCreateTrajectory(env,'')
planner.PlanPath(ravetraj0)
topptraj0 = TOPPopenravepy.FromRaveTraj(robot,ravetraj0)

# Constraints
discrtimestep = 0.005
taumin = zeros(ndof)
taumax = zeros(ndof)
taumin[0:4] = [-30,-50,-25,-15] # Torque limits, only for the shoulder and elbow joints
taumax[0:4] = [30,50,25,15]
constraintstring = str(discrtimestep) + "\n";
constraintstring += string.join([str(v) for v in vmax[0:4]])
constraintstring += TOPPopenravepy.ComputeTorquesConstraints(robot,topptraj0,taumin,taumax,discrtimestep)

# Run TOPP
x = TOPPbindings.TOPPInstance(None,"QuadraticConstraints",constraintstring,str(topptraj0));
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
    topptraj1 = Trajectory.PiecewisePolynomialTrajectory.FromString(x.restrajectorystring)
    #topptraj1 = x.GetOpenRAVEResultTrajectory(env)
    dtplot = 0.01
    TOPPpy.PlotKinematics(topptraj0,topptraj1,dtplot,vmax)
    TOPPopenravepy.PlotTorques(robot,topptraj0,topptraj1,dtplot,taumin,taumax,3)
    print "Trajectory duration before TOPP: ", topptraj0.duration
    print "Trajectory duration after TOPP: ", topptraj1.duration
else:
    print "Trajectory is not time-parameterizable"

# Execute trajectory
if(ret == 1):
    spec = ravetraj0.GetConfigurationSpecification()
    ravetraj1 = TOPPopenravepy.ToRaveTraj(robot,spec,topptraj1)
    robot.GetController().SetPath(ravetraj1)


raw_input()
