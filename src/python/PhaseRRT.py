# -*- coding: utf-8 -*-
# Copyright (C) 2013 Stéphane Caron <caron@ynl.t.u-tokyo.ac.jp>
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

"""Standard RRT in the (s, sd) phase space."""

import numpy
import pylab
import time


class __PhaseRRT(object):
    class Node(object):
        """Node for phase-space RRT."""

        def __init__(self, s, sd, parent=None):
            self.s = s
            self.sd = sd
            self.parent = parent

    def __init__(self, topp_inst, traj, sdbegmin, sdbegmax, ds):
        sd_start = numpy.linspace(sdbegmin, sdbegmax, 42)
        self.topp_inst = topp_inst
        self.traj = traj
        self.ds = ds
        self.sdbegmax = sdbegmax
        self.nodes = [self.Node(0., sd) for sd in sd_start]
        self.end_node = None
        self.max_reached_s = 0.
        self.max_reached_sd = self.sdbegmax

    def found_solution(self):
        return self.end_node is not None

    def plot_tree(self):
        cur_axis = pylab.axis()
        for node in self.nodes:
            s, sd = node.s, node.sd
            pylab.plot([s], [sd], 'bo')
            if node.parent:
                ps, psd = node.parent.s, node.parent.sd
                pylab.plot([ps, s], [psd, sd], 'g-', linewidth=1)
        pylab.axis(cur_axis)

    def plot_path(self, node):
        if node.parent:
            s, sd, ps, psd = node.s, node.sd, node.parent.s, node.parent.sd
            pylab.plot([ps, s], [psd, sd], 'r-', linewidth=3)
            return self.plot_path(node.parent)

    def plot_solution(self):
        cur_axis = pylab.axis()
        if self.end_node:
            self.plot_path(self.end_node)
        pylab.axis(cur_axis)

    def steer(self, node, target):
        """Returns True iff the steering reached the target."""
        interp_step = (target.sd - node.sd) / (target.s - node.s)
        interp_sd = lambda s: (s - node.s) * interp_step + node.sd
        for s in numpy.arange(node.s, target.s, self.ds):
            sd = interp_sd(s)
            alpha = self.topp_inst.GetAlpha(s, sd)
            beta = self.topp_inst.GetBeta(s, sd)
            # stepping condition is: alpha / sd <= sdd <= beta / sd
            if sd <= 0 or not (alpha <= interp_step * sd <= beta):
                return False
        return True

    def extend(self, target, k=10):
        from random import sample
        candidates = [node for node in self.nodes if node.s < target.s]
        if len(candidates) > k:
            candidates = sample(candidates, k)
        for candidate in candidates:
            if not self.steer(candidate, target):
                continue
            new_node = self.Node(target.s, target.sd, candidate)
            self.nodes.append(new_node)
            if target.s >= self.traj.duration:
                self.end_node = new_node
            if target.s > self.max_reached_s:
                self.max_reached_s = target.s
            if target.sd > 0.75 * self.max_reached_sd:
                self.max_reached_sd *= 1.25

    def run(self, max_nodes, time_budget):
        """Runs until the time budget is exhausted."""
        smax = self.traj.duration
        svar = smax / 10.
        start_time = time.time()
        while not self.found_solution():
            if len(self.nodes) > max_nodes \
               or time.time() - start_time > time_budget:
                break
            if pylab.random() < 0.1:
                s = pylab.random() * self.traj.duration
            else:
                s = pylab.normal(.5 * (smax + self.max_reached_s), svar)
                s = max(0., min(smax, s))
            sd = pylab.random() * self.max_reached_sd
            self.extend(self.Node(s, sd))
            if sd < self.max_reached_sd / 10:  # happens 1/10 times
                sdend = pylab.random() * self.max_reached_sd
                self.extend(self.Node(smax, sdend))
        print "RRT run time: %d s" % int(time.time() - start_time)


def TryRRT(topp_inst, traj, sdbegmin, sdbegmax, ds=1e-3, max_nodes=500,
           time_budget=360):
    rrt = __PhaseRRT(topp_inst, traj, sdbegmin, sdbegmax, ds)
    rrt.run(max_nodes, time_budget)
    return rrt
