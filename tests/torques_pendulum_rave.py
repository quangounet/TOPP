# -*- coding: utf-8 -*-
# Copyright (C) 2013 Quang-Cuong Pham <cuong.pham@normalesup.org>
#
# This file is part of the Time-Optimal Path Parameterization (TOPP) library.
# TOPP is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <http://www.gnu.org/licenses/>.

import sys
sys.path.append('..')

import TOPPbindings
import TOPPpy
import time
from pylab import ion, array, ones, axis, clf
from openravepy import Environment

ion()

########################### Robot ################################
env = Environment()  # create openrave environment
#------------------------------------------#
robotfile = "../robots/twodof.robot.xml"
env.Load(robotfile)
robot = env.GetRobots()[0]
robot.SetTransform(array([[0, 0, 1, 0],
                          [0, 1, 0, 0],
                          [-1, 0, 0, 0.3],
                          [0, 0, 0, 1]]))
#------------------------------------------#
grav = [0, 0, -9.8]
n = robot.GetDOF()
dof_lim = robot.GetDOFLimits()
vel_lim = robot.GetDOFVelocityLimits()
robot.SetDOFLimits(-10 * ones(n), 10 * ones(n))
robot.SetDOFVelocityLimits(100 * vel_lim)


############################ Tunings ############################
discrtimestep = 0.0001
integrationtimestep = discrtimestep
reparamtimestep = 0  # auto
passswitchpointnsteps = 10
tuningsstring = "%f %f %f %d" % (discrtimestep, integrationtimestep,
                                 reparamtimestep, passswitchpointnsteps)


############################ Trajectory ############################

# TODO: no solution for discrtimestep=1e-3, solutions when =1e-4
sdbeg_min, sdbeg_max = 0.0, 0.0001
trajectorystring = """1.000000
2
0.0 0.0 -1.04335890092 0.547353298841
0.0 0.0 -2.24856384056 1.36915743417"""

# TODO: x.resduration is too small (~1e-324)
sdbeg_min, sdbeg_max = 0.0, 0.0001
trajectorystring = """1.000000
2
0.0 -0.316841821141 -0.517973306853 0.413832196233
0.0 -1.83079350169 1.68090295997 -0.480535237087"""


#------------------------------------------#
traj0 = TOPPpy.PiecewisePolynomialTrajectory.FromString(trajectorystring)

############################ Constraints ############################
vmax = array([0, 0])
taumin = array([-11, -7])
taumax = array([11, 7])
t0 = time.time()
constraintstring = ' '.join([str(x) for x in taumin])
constraintstring += "\n" + ' '.join([str(a) for a in taumax])
constraintstring += "\n" + ' '.join([str(a) for a in vmax])
#------------------------------------------#


############################ Run TOPP ############################
t1 = time.time()
x = TOPPbindings.TOPPInstance("TorqueLimitsRave", constraintstring,
                              trajectorystring, tuningsstring, robot)
t2 = time.time()
#ret = x.RunComputeProfiles(0,0)
ret = x.RunVIP(sdbeg_min, sdbeg_max)
print ret
t3 = time.time()

#print x.resduration
print "sdendmin =", x.sdendmin
print "sdendmax =", x.sdendmax

# if ret == 1:
#     x.ReparameterizeTrajectory()

t4 = time.time()

################ Plotting the MVC and the profiles #################
x.WriteProfilesList()
x.WriteSwitchPointsList()
profileslist = TOPPpy.ProfilesFromString(x.resprofilesliststring)
switchpointslist = TOPPpy.SwitchPointsFromString(x.switchpointsliststring)


def replot(prec=None):
    cur_axis = axis()
    clf()
    TOPPpy.PlotProfiles(profileslist, switchpointslist, 1)
    axis(cur_axis)

axis([0, traj0.duration, 0, 100])
replot()
#TOPPpy.PlotAlphaBeta(x)

# validate with RRT in the (s, sd) plane, in case no solution was found
if x.sdendmin < 0 and raw_input("Check with RRT? [y/N] ") == 'y':
    rrt = TOPPpy.TryRRT(x, traj0, sdbeg_min, sdbeg_max, discrtimestep)
    if rrt.found_solution():
        rrt.plot_solution()


##################### Plotting the trajectories #####################
# if(ret == 1):
#     x.WriteResultTrajectory()
#     traj1 = PiecewisePolynomialTrajectory.FromString(x.restrajectorystring)
#     dtplot = 0.01
#     TOPPpy.PlotKinematics(traj0,traj1,dtplot,vmax)
#     TOPPopenravepy.PlotTorques(robot,traj0,traj1,dtplot,taumin,taumax,3)


print "\n--------------"
print "Python preprocessing: ", t1 - t0
print "Building TOPP Instance: ", t2 - t1
print "Compute profiles: ", t3 - t2
print "Reparameterize trajectory: ", t4 - t3
print "Total: ", t4 - t0
if ret == 1:
    print "Trajectory duration (estimate): ", x.resduration

raw_input()
