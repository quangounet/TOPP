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

sdbeg_min, sdbeg_max = 0.0, 4.11542810233
trajectorystring = """1.000000
2
-0.0950348403575 0.00709188111616 -0.735957024034 0.46461228402
-0.165204811611 0.164744178028 -4.02513169685 2.5711050706"""

# TODO: forward profile crosses the MVC
#discrtimestep = 0.005
#sdbeg_min, sdbeg_max = 0.0, 0.0001
#trajectorystring = """1.000000
#2
#0.0 -0.000441411097554 -0.16346942929 0.109925764182
#0.0 -0.0902570804637 -0.811688277595 0.600976616276"""

# TODO: minimum velocity curve in the beginning
sdbeg_min, sdbeg_max = 0.000000, 8.808698
trajectorystring = """1.000000
2
-0.420982931773 -0.111291845962 1.77275676137 -1.33551682399
-0.630425778805 0.0894067072732 0.0852189728776 0.290595287026"""

# TODO: lift degree limitation
trajectorystring = """1.000000
2
-0.537652888623 0.234321204132 0.0
0.183981719702 0.488818919925 3.33066907388e-16
1.000000
2
-0.303331684468 0.234321204132 0.0315593379862 0.507751607858
0.672800639677 0.488818919925 -2.25421681376 1.50275249907"""




# Trajectory 510
discrtimestep = 0.005000
integrationtimestep = discrtimestep
sdbeg_min, sdbeg_max = 0.000000, 2.963261
tuningsstring = """0.005000 0.005000 0.000000 10"""
constraintstring = """-11.0 -7.0
11.0 7.0
0 0"""
trajectorystring = """1.000000
2
0.0598839777986 -0.767543731454 17.0567694943 -13.207517087
0.43299721927 1.80131573114 -5.85537754605 3.62106459564"""

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
#ret = x.RunComputeProfiles(0, 1e-4)
#print "RunComputeProfiles:", ret
ret = x.RunVIP(sdbeg_min, sdbeg_max)
print "RunVIP:", ret
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
msg = "No solution found. Check with RRT? [y/N] "
if x.sdendmin < 0 and raw_input(msg) == 'y':
    rrt = TOPPpy.TryRRT(x, traj0, sdbeg_min, sdbeg_max, discrtimestep, max_nodes=200)
    rrt.plot_tree()
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
