# -*- coding: utf-8 -*-
# Copyright (C) 2014 Quang-Cuong Pham <cuong.pham@normalesup.org>
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

print("\n***********************************************\nNB: This test file requires OpenRAVE and cvxopt\n***********************************************\n")


import time
from pylab import *
from numpy import *
from openravepy import *
from TOPP import Trajectory
from TOPP import TOPPpy
from TOPP import TOPPopenravepy
from TOPP import TOPPbindings
from TOPP import ClosedChain
from TOPP import Bimanual

try:
    input = raw_input
except NameError:
    pass

ion()

############# Load environment and robot ############

env = Environment()
env.SetViewer('qtcoin')
env.Load('../robots/ArmFull.xml')
env.Load('../robots/ArmCut.xml')
robot1 = env.GetRobots()[0]
robot2 = env.GetRobots()[1]

t = RaveCreateKinBody(env,'')
t.InitFromBoxes(array([array([0.3,0,-0.06,0.5,0.2,0.02])]),True)
t.SetName('Box')
g=t.GetLinks()[0].GetGeometries()[0]
g.SetAmbientColor([0.5,0.5,0.5])
g.SetDiffuseColor([0.5,0.5,0.5])
env.Add(t,True)


M = array([[ 0.99981043, -0.01138393,  0.01579602,  0.23296329],
       [-0.01788201, -0.21589215,  0.97625346, -0.95132697],
       [-0.00770337, -0.97635085, -0.21605479,  0.3695007 ],
       [ 0.        ,  0.        ,  0.        ,  1.        ]])
env.GetViewer().SetCamera(M)

T1 = array([[1.,0,0,0],
            [0,0,-1,0],
            [0,1,0,0],
            [0,0,0,1]])
T2 = array(T1)
T2[0,3] = 0.6

constrainedlink = 9
obj = robot1.GetLinks()[constrainedlink]

taumin = -20*ones(6)
taumax = 20*ones(6)

robot = Bimanual.Robot()
robot.robot1 = robot1
robot.robot2 = robot2
robot.T1 = T1
robot.T2 = T2
robot.freedofs = [0,1,2]
robot.dependentdofs = [3,4,5]
robot.constrainedlink = 9
robot.Gdofs = [0,1,2]
robot.Sdofs = [3,4,5]
robot.actuated = [True]*6
robot.taumin = taumin
robot.taumax = taumax
robot.vmax = [3]*6
robot.q_range = [[-pi,pi],[-pi,pi],[-pi,pi]]    

tunings = Bimanual.Tunings()
tunings.duration = 1
tunings.nchunks = 10
tunings.chunksubdiv = 20
tunings.tol_jacobian = 1e-2
tunings.discrtimestep = 5e-3

############# Kinematics to find start and goal configurations ############

# Start configuration
# Robot 1
q_start1 = array([ 1.5, -0.91191789, -0.52655972])
robot1.SetTransform(T1)
robot1.SetDOFValues(q_start1)
desiredpose = Bimanual.Getxztheta(obj.GetTransform())
# Robot 2
desiredpose2 = desiredpose + [0,0,pi]
robot1.SetTransform(T2)
q_start2, error = Bimanual.IK3(robot1,desiredpose2,q_start=array([pi/2,pi/2,0]))
robot1.SetTransform(T1)
robot1.SetDOFValues(q_start1)
robot2.SetDOFValues(q_start2[0:2])
q_start = zeros(6)
q_start[0:3] = q_start1
q_start[3:6] = q_start2
robot2.SetTransform(T2)
#print "Initial config"
#input()

# Goal configuration
# Robot 1
q_goal1 = q_start1 + array([0.5,-0.5,-0.2])
robot1.SetTransform(T1)
robot1.SetDOFValues(q_goal1)
desiredpose = Bimanual.Getxztheta(obj.GetTransform())
# Robot 2
desiredpose2 = desiredpose + [0,0,pi]
robot1.SetTransform(T2)
q_goal2, error = Bimanual.IK3(robot1,desiredpose2,q_start=array([pi,0,0])-q_goal1)
robot1.SetTransform(T1)
robot1.SetDOFValues(q_goal1)
robot2.SetDOFValues(q_goal2[0:2])
q_goal = zeros(6)
q_goal[0:3] = q_goal1
q_goal[3:6] = q_goal2
#print "Goal config"
#input()


########### Interpolate between the two configurations ##############

freestartvelocities = array([0.1,0.1,0.1])
freegoalvelocities = array([0.5,0,0])
trajtotal = Bimanual.Interpolate(robot, tunings, q_start, q_goal, freestartvelocities, freegoalvelocities)
trajectorystring = str(trajtotal)


################ Bobrow with actuation redundancy ##################

constraintstring = str(tunings.discrtimestep)
constraintstring += "\n" + " ".join([str(v) for v in robot.vmax])
scaledowncoef = 0.99
robot.taumin = taumin * scaledowncoef # safety bound
robot.taumax = taumax * scaledowncoef
t0 = time.time()
constraintstring += Bimanual.ComputeConstraints(robot,tunings,trajtotal)
robot.taumin = taumin
robot.taumax = taumax
t1 = time.time()

x = TOPPbindings.TOPPInstance(None,"PolygonConstraints",constraintstring,trajectorystring)
x.integrationtimestep = 1e-3
ret = x.RunComputeProfiles(1,1)
t2 = time.time()

print("Compute Polygon constraints:", t1-t0, "seconds")
print("Note : Compute Polygon is faster with the C++ version, checkout branch tomas-develop")
print("Run TOPP:", t2-t1, "seconds")

x.WriteProfilesList()
x.WriteSwitchPointsList()
profileslist = TOPPpy.ProfilesFromString(x.resprofilesliststring)
switchpointslist = TOPPpy.SwitchPointsFromString(x.switchpointsliststring)

if(ret == 1):
    x.ReparameterizeTrajectory()
    x.WriteResultTrajectory()
    trajtotal2 = Trajectory.PiecewisePolynomialTrajectory.FromString(x.restrajectorystring)

    dt=0.01
    robot1.SetTransform(T1)
    robot2.SetTransform(T2)
    traj = trajtotal
    for t in arange(0,traj.duration,dt):
        q = traj.Eval(t)
        robot1.SetDOFValues(q[0:3])
        robot2.SetDOFValues(q[3:5])
        time.sleep(dt)

    TOPPpy.PlotProfiles(profileslist,switchpointslist,4)
    TOPPpy.PlotKinematics(trajtotal,trajtotal2,dt,robot.vmax)
    Bimanual.PlotTorques(robot,trajtotal,trajtotal2,dt)

input()

