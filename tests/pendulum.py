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
from pylab import array, ones, ion
from openravepy import Environment

VERBOSE = False  # set to True if -v is a command-line argument

vmax = [0, 0]
discrtimestep = 0.001
integrationtimestep = discrtimestep
reparamtimestep = 0  # auto
passswitchpointnsteps = 10
robotfile = "../robots/twodof.robot.xml"
dtplot = 0.01

def append_traj(traj_list, traj_str, tauref, sd_min=0., sd_max=1e-4):
    traj_list.append((traj_str, sd_min, sd_max, -tauref, +tauref))

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

append_traj(impossible_trajs, """1.000000
2
-0.496005602127 -0.496005602078 -7.64552243845 5.49594098906
-0.879406406487 -0.879406406399 4.27317393732 -2.51436112444""", tau_11_7)


class TorquePendulumExec(unittest.TestCase):
    def setUp(self):
        self.discrtimestep = discrtimestep
        self.dtplot = dtplot
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
        from TOPPopenravepy import ComputeTorquesConstraintsLegacy,ComputeTorquesConstraints
        from TOPPbindings import TOPPInstance
        self.cur_traj_str = traj_str
        self.cur_tau_min = tau_min
        self.cur_tau_max = tau_max
        self.traj0 = TOPPpy.PiecewisePolynomialTrajectory.FromString(traj_str)
        self._t0 = time.time()
        
        # QuadraticConstraints (topp0)
        self.constraintstring0 = ' '.join([str(v) for v in vmax])
        self.constraintstring0 += TOPPopenravepy.ComputeTorquesConstraints(self.robot,self.traj0,tau_min,tau_max,discrtimestep)
        self.topp0 = TOPPInstance("QuadraticConstraints",
                                 self.constraintstring0,
                                 traj_str,
                                 self.tuningsstring,False)
        self.ret0 = self.topp0.RunComputeProfiles(0, 0)


        # TorqueLimits (topp)
        self.constraintstring = ' '.join([str(v) for v in tau_min]) + "\n"
        self.constraintstring += ' '.join([str(v) for v in tau_max]) + "\n"
        self.constraintstring += ' '.join([str(v) for v in vmax])
        self.constraintstring += ComputeTorquesConstraintsLegacy(self.robot,
                                                                 self.traj0,
                                                                 tau_min,
                                                                 tau_max,
                                                                 discrtimestep)
        self._t1 = time.time()
        self.topp = TOPPInstance("TorqueLimits",
                                 self.constraintstring,
                                 traj_str,
                                 self.tuningsstring,False)
        self._t2 = time.time()
        self.ret = self.topp.RunComputeProfiles(0, 0)
        self._t3 = time.time()
        if self.ret == 1:
            self.topp.ReparameterizeTrajectory()
        self._t4 = time.time()

        # TorqueLimits AVP (topp1)
        self.topp1 = TOPPInstance("TorqueLimits",
                                 self.constraintstring,
                                 traj_str,
                                 self.tuningsstring,False)
        self._t5 = time.time()
        self.ret_vip = self.topp1.RunVIP(sd_min, sd_max)
        self._t6 = time.time()
        self.sd_beg_min = sd_min
        self.sd_beg_max = sd_max
        self.sd_end_min = self.topp1.sdendmin
        self.sd_end_max = self.topp1.sdendmax

    def print_info(self):
        print "Torque limits:"
        print "- min:", self.cur_tau_min
        print "- max:", self.cur_tau_max
        print "Trajectory string:"
        print self.cur_traj_str

        print "\nComputation times:"
        print "- Python preprocessing:      %.2f s" % (self._t1 - self._t0)
        print "- Building TOPP Instance:    %.2f s" % (self._t2 - self._t1)
        print "- Compute profiles:          %.2f s" % (self._t3 - self._t2)
        print "- Reparameterize trajectory: %.2f s" % (self._t4 - self._t3)
        print "- Building TOPP Instance:    %.2f s" % (self._t5 - self._t4)
        print "- AV propagation:            %.2f s" % (self._t6 - self._t5)
        print "Total: %.2f s" % (self._t6 - self._t0)

        print "\nResult:"
        if self.ret == 1:
            print "- Old traj. duration: ", self.traj0.duration
            print "- New traj. duration (estimate): ", self.topp.resduration
        print "- AVP return code:", self.ret_vip
        print "- AVP sd_beg:", (self.sd_beg_min, self.sd_beg_max)
        print "- AVP sd_end:", (self.sd_end_min, self.sd_end_max)

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
            self.run_topp(*trajuple)
            if VERBOSE:
                print "\n\n---------------------------------------------------"
                print "Traversable test results"
                self.print_info()
            self.assertEqual(self.ret, 1)
            self.assertTrue(abs(self.topp0.resduration-self.topp.resduration)<1e-2)
            self.assertNotEqual(self.ret_vip, 0)

    def test_impossible(self):
        for i, trajuple in enumerate(impossible_trajs):
            self.run_topp(*trajuple)
            if VERBOSE:
                print "\n\n---------------------------------------------------"
                print "Impossible test results"
                self.print_info()
            self.assertNotEqual(self.ret, 1)
            self.assertEqual(self.ret_vip, 0)


if __name__ == '__main__':
    VERBOSE = '-v' in sys.argv
    unittest.main()
