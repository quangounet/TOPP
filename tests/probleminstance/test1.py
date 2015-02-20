from openravepy import *
from TOPP import *
import matplotlib.pyplot as plt
import pylab

pylab.ion()

env = Environment()
env.Load('tomas-develop_hrp2_14.env.xml')

collisionChecker = RaveCreateCollisionChecker(env, 'ode')
env.SetCollisionChecker(collisionChecker)

robot = env.GetRobots()[0]

problem = int(raw_input('input the problem number: '))

with open('constraints' + str(problem), 'r') as f:
    constraintsstring = f.read()

with open('trajectory' + str(problem), 'r') as f:
    trajectorystring = f.read()

integrationtimestep = 1e-3
reparamtimestep = 1e-2
passswitchpointnsteps = 5
discrtimestep = 1e-3

x = TOPPbindings.TOPPInstance(robot, 'PolygonConstraints', constraintsstring, trajectorystring)
x.integrationtimestep = integrationtimestep
x.reparamtimestep = reparamtimestep
x.passswitchpointnsteps = passswitchpointnsteps

retvip = x.RunVIP(0, 0)
print "VIP returns", retvip
print "(sdendmin, sdendmax) = ({0}, {1})".format(x.sdendmin, x.sdendmax)

"""

ret = x.RunComputeProfiles(0, 0)

x.WriteProfilesList()
x.WriteSwitchPointsList()

fignum = 1
profileslist = TOPPpy.ProfilesFromString(x.resprofilesliststring)
switchpointslist = TOPPpy.SwitchPointsFromString(x.switchpointsliststring)
TOPPpy.PlotProfiles(profileslist, switchpointslist, fignum)

"""
