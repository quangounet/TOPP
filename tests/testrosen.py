import string
from pylab import *
from numpy import *
from openravepy import *
from TOPP import TOPPbindings
from TOPP import TOPPpy
from TOPP import TOPPopenravepy
from TOPP import Trajectory
from TOPP import Utilities

ion()

trajfile = "data/mujin2.topp.traj"
constraintsfile = "data/mujin2.topp.constraints"

handle = open(trajfile,"r")
trajectorystring = handle.read()
traj0 = Trajectory.PiecewisePolynomialTrajectory.FromString(trajectorystring)

handle = open(constraintsfile,"r")
constraintsstring = handle.read()


x = TOPPbindings.TOPPInstance(None,"QuadraticConstraints",constraintsstring,trajectorystring);

ret = x.RunComputeProfiles(0,0)

if(ret == 1):
    x.ReparameterizeTrajectory()


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
    TOPPpy.PlotKinematics(traj0,traj1,dtplot)
    print "Trajectory duration before TOPP: ", traj0.duration
    print "Trajectory duration after TOPP: ", traj1.duration
else:
    print "Trajectory is not time-parameterizable"


raw_input()
