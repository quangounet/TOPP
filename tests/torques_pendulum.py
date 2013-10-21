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
import string
from pylab import array, ones, ion
from openravepy import Environment


class TorquePendulumExec(object):
    def __init__(self,
                 taumin=array([-8, -4]),
                 taumax=array([+8, +4]),
                 vmax=[0, 0],
                 discrtimestep=0.01,
                 integrationtimestep=0.01,
                 reparamtimestep=0.01,
                 passswitchpointnsteps=2,
                 dtplot=0.01):
        self.taumin = taumin
        self.taumax = taumax
        self.discrtimestep = discrtimestep
        self.dtplot = dtplot
        self.vmax = vmax
        self.env = Environment()  # create openrave environment
        robotfile = "../robots/twodof.robot.xml"
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

    def run(self, traj_str):
        from TOPPopenravepy import ComputeTorquesConstraints
        from TOPPbindings import TOPPInstance
        self.traj0 = TOPPpy.PiecewisePolynomialTrajectory.FromString(traj_str)
        self.t0 = time.time()
        self.constraintstring = string.join([str(v) for v in self.vmax])
        self.constraintstring += ComputeTorquesConstraints(self.robot,
                                                           self.traj0,
                                                           self.taumin,
                                                           self.taumax,
                                                           self.discrtimestep)
        self.t1 = time.time()
        self.topp = TOPPInstance("QuadraticConstraints",
                                 self.constraintstring,
                                 traj_str,
                                 self.tuningsstring)
        self.t2 = time.time()
        self.ret = self.topp.RunComputeProfiles(0, 0)
        self.t3 = time.time()
        if self.ret == 1:
            self.topp.ReparameterizeTrajectory()
        self.t4 = time.time()
        self.topp.WriteProfilesList()
        self.topp.WriteSwitchPointsList()

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
            TOPPpy.PlotKinematics(self.traj0, traj1, self.dtplot, self.vmax)
            TOPPopenravepy.PlotTorques(self.robot, self.traj0, traj1,
                                       self.dtplot, self.taumin, self.taumax,
                                       3)


trajectories = []

trajectories.append("""1.000000
2
0.0 -1.4059993022 -6.56388799235 4.82829464097
0.0 0.0 0.0 0.0""")

trajectories.append("""1.000000
2
0.0 0.0 0.200440827515 -0.132913533868
0.0 0.0 -1.71157060946 1.19264863791""")

if __name__ == '__main__':
    topp_exec = TorquePendulumExec()
    try:
        for i, traj_str in enumerate(trajectories):
            print "Running on trajectory %d... " % i
            topp_exec.run(traj_str)
            print "\n----- Trajectory %d: OK -----" % i
            topp_exec.print_comp_times()
            print ""
    except:
        print ""
        print "-------------------- TEST FAILURE --------------------"
        print "Error for trajectory:"
        print traj_str
        raise
