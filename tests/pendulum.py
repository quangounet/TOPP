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

import TOPPpy
import TOPPopenravepy
import time
import unittest
from pylab import array, ones, ion, axis
from openravepy import Environment

vmax = [0, 0]
discrtimestep = 0.001
integrationtimestep = discrtimestep
reparamtimestep = 0  # auto
passswitchpointnsteps = 5
robotfile = "../robots/twodof.robot.xml"
dtplot = 0.01

# "QuadraticConstraints" or "TorqueLimits"
constraints_type = "TorqueLimits"


def append_traj(traj_list, traj_str, tauref, sd_min=0., sd_max=1e-4):
    traj_list.append((traj_str, sd_min, sd_max, +tauref, -tauref))

tau_8_4 = array([8.000, 4.000])
tau_11_7 = array([11.000, 7.000])

############################## Test Trajectories ##############################

#
# Traversable trajectories
#

traversable_trajs = []

append_traj(traversable_trajs, """1.000000
2
0.0 0.280641438732 -1.38946387942 0.84562395633
0.0 3.46940504119 -10.3267357647 5.62007589178""", tau_11_7)

append_traj(traversable_trajs, """1.000000
2
0.0 -0.420982931761 -1.66533453694e-16 1.11022302463e-16
0.0 -0.630425778814 1.66533453694e-16 -2.22044604925e-16""", tau_11_7)

append_traj(traversable_trajs, """1.000000
2
0.0 -0.626567074367 -1.66533453694e-16
0.0 -0.458558522489 0.0""", tau_11_7)

append_traj(traversable_trajs, """1.000000
2
0.0 0.0935471462441 -1.10530824184 0.748562611233
0.0 1.15646834706 -4.94523847985 2.55151530108""", tau_11_7)

append_traj(traversable_trajs, """1.000000
2
0.0 0.0280641438732 -0.974342237093 0.683079608862
0.0 0.346940504119 -3.32618279396 1.74198745813""", tau_11_7)

append_traj(traversable_trajs, """1.000000
2
-0.496005602127 -0.496005602078 -7.64552243845 5.49594098906
-0.879406406487 -0.879406406399 4.27317393732 -2.51436112444""", tau_11_7)

append_traj(traversable_trajs, """1.000000
2
-0.0539850762113 -0.0539850762059 0.131202643008 -0.0710088455766
-0.300968741813 -0.300968741783 0.480216768411 -0.231499819896""", tau_11_7)

#
# Non-traversable trajectories
#

impossible_trajs = []

append_traj(impossible_trajs, """1.000000
2
0.0 0.0 9.6622280813 -6.6326750309
0.0 0.0 0.289254321757 -0.184727934151""", tau_8_4)

append_traj(impossible_trajs, """1.000000
2
0.0 0.0 11.0005555849 -7.9710025345
0.0 0.0 0.232496359281 -0.127969971675""", tau_8_4)

append_traj(impossible_trajs, """1.000000
2
-0.0477763549985 -0.129943080005 -9.55936508575 6.59549186717
-0.353220535126 -0.447152621088 1.86760269163 -1.06722953542""", tau_11_7)

append_traj(impossible_trajs, """1.000000
2
-0.0539850762113 -0.0539850762059 -0.609921814791 -0.446413715041
-0.300968741813 -0.300968741783 -0.511789732519 0.583618460003""", tau_11_7)

#
# Test cases
#

class TorquePendulumExec(unittest.TestCase):
    def setUp(self):
        self.discrtimestep = discrtimestep
        self.dtplot = dtplot
        self.constraints_type = constraints_type
        self.env = Environment()  # create openrave environment
        self.env.Load(robotfile)
        self.robot = self.env.GetRobots()[0]
        self.robot.SetTransform(array([
            [0, 0, 1, 0],
            [0, 1, 0, 0],
            [-1, 0, 0, 0.3],
            [0, 0, 0, 1]]))

        # Robot
        n = self.robot.GetDOF()
        vel_lim = self.robot.GetDOFVelocityLimits()
        self.robot.SetDOFLimits(-10 * ones(n), 10 * ones(n))
        self.robot.SetDOFVelocityLimits(100 * vel_lim)
        self.robot.GetEnv().GetPhysicsEngine().SetGravity([0, 0, -9.81])

        # Tunings
        self.tuningsstring = "%f %f %f %d" % (discrtimestep,
                                              integrationtimestep,
                                              reparamtimestep,
                                              passswitchpointnsteps)

        self.ret = None
        self.ret_vip = None

    def run_topp(self, traj_str, sd_min, sd_max, tau_min, tau_max):
        print traj_str
        from TOPPopenravepy import ComputeTorquesConstraintsLegacy
        from TOPPbindings import TOPPInstance
        self.cur_tau_min = tau_min
        self.cur_tau_max = tau_max
        self.traj0 = TOPPpy.PiecewisePolynomialTrajectory.FromString(traj_str)
        self.t0 = time.time()
        self.constraintstring = ' '.join([str(v) for v in tau_min]) + "\n"
        self.constraintstring += ' '.join([str(v) for v in tau_max]) + "\n"
        self.constraintstring += ' '.join([str(v) for v in vmax])
        self.constraintstring += ComputeTorquesConstraintsLegacy(self.robot,
                                                                 self.traj0,
                                                                 tau_min,
                                                                 tau_max,
                                                                 discrtimestep)
        self.t1 = time.time()
        self.topp = TOPPInstance(self.constraints_type,
                                 self.constraintstring,
                                 traj_str,
                                 self.tuningsstring)
        self.t2 = time.time()
        self.ret = self.topp.RunComputeProfiles(0, 0)
        self.t3 = time.time()
        if self.ret == 1:
            self.topp.ReparameterizeTrajectory()
            print "---- Time duration =", self.topp.resduration
        else:
            print "---- Impossible from (0,0)"
        self.t4 = time.time()

        # run VIP as well
        self.topp = TOPPInstance(self.constraints_type,
                                 self.constraintstring,
                                 traj_str,
                                 self.tuningsstring)
        self.ret_vip = self.topp.RunVIP(sd_min, sd_max)
        self.sd_end_min = self.topp.sdendmin
        self.sd_end_max = self.topp.sdendmax
        print "---- Ret_vip =", self.ret_vip,
        print "and sd:", (self.sd_end_min, self.sd_end_max)

    def print_comp_times(self):
        print "Python preprocessing: ", (self.t1 - self.t0)
        print "Building TOPP Instance: ", (self.t2 - self.t1)
        print "Compute profiles: ", (self.t3 - self.t2)
        print "Reparameterize trajectory: ", (self.t4 - self.t3)
        print "Total: ", (self.t4 - self.t0)
        if self.ret == 1:
            print "Trajectory duration (estimate): ", self.topp.resduration
            print "Trajectory duration: ", self.traj1.duration

    def plot_result(self):
        from TOPPpy import ProfilesFromString, SwitchPointsFromString
        from TOPPpy import PiecewisePolynomialTrajectory
        ion()
        profileslist = ProfilesFromString(self.topp.resprofilesliststring)
        splist_str = self.topp.switchpointsliststring
        switchpointslist = SwitchPointsFromString(splist_str)
        TOPPpy.PlotProfiles(profileslist, switchpointslist, 4)
        if self.ret == 1:
            self.topp.WriteResultTrajectory()
            restraj = self.topp.restrajectorystring
            traj1 = PiecewisePolynomialTrajectory.FromString(restraj)
            TOPPpy.PlotKinematics(self.traj0, traj1, self.dtplot, vmax)
            TOPPopenravepy.PlotTorques(self.robot, self.traj0, traj1,
                                       self.dtplot, self.cur_tau_min,
                                       self.cur_tau_max, 3)

    def test_traversable(self):
        for i, trajuple in enumerate(traversable_trajs):
            print "\n\n\nTest traversable"
            self.run_topp(*trajuple)
            self.assertEqual(self.ret, 1)
            self.assertNotEqual(self.ret_vip, 0)

    def test_impossible(self):
        for i, trajuple in enumerate(impossible_trajs):
            print "\n\n\nTest impossible"
            self.run_topp(*trajuple)
            self.assertNotEqual(self.ret, 1)
            self.assertEqual(self.ret_vip, 0)


if __name__ == '__main__':
    unittest.main()
