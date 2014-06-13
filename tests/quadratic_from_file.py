# -*- coding: utf-8 -*-

from pylab import *
from numpy import *


from TOPP import TOPPbindings
from TOPP import TOPPpy
from TOPP import Trajectory

with open("mujin16.topp.traj", "r") as f:
    trajectorystring = f.read()

with open("mujin16.topp.constraints", "r") as f:
    constraintstring = f.read()

x = TOPPbindings.TOPPInstance(
    None, "QuadraticConstraints", constraintstring, trajectorystring)

x.integrationtimestep = 1e-5
ret = x.RunComputeProfiles(0, 0)
if(ret == 1):
    x.ReparameterizeTrajectory()

ion()
x.WriteProfilesList()
x.WriteSwitchPointsList()
profileslist = TOPPpy.ProfilesFromString(x.resprofilesliststring)
switchpointslist = TOPPpy.SwitchPointsFromString(x.switchpointsliststring)
TOPPpy.PlotProfiles(profileslist, switchpointslist, 4)
#TOPPpy.PlotAlphaBeta(x)
traj0 = Trajectory.PiecewisePolynomialTrajectory.FromString(trajectorystring)
TOPPpy.PlotKinematics(traj0,traj0)


if(ret == 1):
    x.WriteResultTrajectory()
    traj0 = Trajectory.PiecewisePolynomialTrajectory.FromString(trajectorystring)
    traj1 = Trajectory.PiecewisePolynomialTrajectory.FromString(x.restrajectorystring)
    print "Trajectory duration before TOPP: ",  traj0.duration
    print "Trajectory duration after TOPP: ",  traj1.duration
else:
    print "Trajectory is not time-parameterizable"

raw_input()
