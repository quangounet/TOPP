# -*- coding: utf-8 -*-
# Copyright (C) 2013 Quang-Cuong Pham <cuong.pham@normalesup.org >
#
# This file is part of the Time-Optimal Path Parameterization (TOPP) library.
# TOPP is free software: you can redistribute it and/or modify
# it under the terms of the GNU Lesser General Public License as published by
# the Free Software Foundation,  either version 3 of the License,  or
# at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU Lesser General Public License for more details.
#
# You should have received a copy of the GNU Lesser General Public License
# along with this program.  If not,  see <http://www.gnu.org/licenses/ > .

import sys
sys.path.append('../..')

import TOPPbindings
import TOPPopenravepy
import pylab
import time

from openravepy import Environment
from Trajectory import PiecewisePolynomialTrajectory

robotfile = "../../../robots/hrp4r.dae"
instancevar = {}  # read from test file


def play_trajectory(traj):
    dt = float(instancevar['tuningsstring'].split(' ')[0])
    s = instancevar['constraintstring'].split('\n')[0]
    activedofs_vect = pylab.array(map(float, s.split(' ')))
    active_dofs = [i for i, v in enumerate(activedofs_vect) if v > 0]
    robot.SetActiveDOFs(active_dofs)
    for t in pylab.arange(0, traj.duration, dt):
        robot.SetActiveDOFValues(traj.Eval(t))
        time.sleep(dt)


view_transform = pylab.array(
    [[-0.81072723, -0.02540599, -0.58487254,  1.37787211],
     [+0.58541559, -0.04056306, -0.80971799,  1.72345638],
     [-0.00315253, -0.99885393,  0.04775864,  0.70553762],
     [+0.,          0.,          0.,          1.]])


if __name__ == "__main__":
    assert len(sys.argv) > 1, "Please supply a test file"
    testfile = sys.argv[1]
    execfile(testfile, instancevar)
    env = Environment()
    baselinkname = "BODY"
    robot = TOPPopenravepy.LoadFloat(env, robotfile, baselinkname)
    n = robot.GetDOF()
    dof_lim = robot.GetDOFLimits()
    vel_lim = robot.GetDOFVelocityLimits()
    dof_lim[1][18] = 0.5
    dof_lim[1][28] = 1.2
    dof_lim[1][29] = 0.8
    robot.SetDOFLimits(dof_lim[0], dof_lim[1])
    robot.SetDOFVelocityLimits(100 * vel_lim)
    env.SetViewer('qtcoin')
    env.GetViewer().SetCamera(view_transform)

    print "Computing TOPP..."
    x = TOPPbindings.TOPPInstance("ZMPTorqueLimits",
                                  instancevar['constraintstring'],
                                  instancevar['trajectorystring'],
                                  instancevar['tuningsstring'], robot)
    ret = x.RunComputeProfiles(0, 0)
    x.WriteResultTrajectory()

    if len(x.restrajectorystring) > 0:
        print "Showing reparameterized trajectory..."
        traj = PiecewisePolynomialTrajectory.FromString(x.restrajectorystring)
        play_trajectory(traj)
    else:
        print "TOPP found no feasible reparameterization."
