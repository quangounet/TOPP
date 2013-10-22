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


from pylab import array, arange, double, gca, plot, figure
from pylab import clf, hold, axis, title, zeros

import TOPPbindings

from Trajectory import PiecewisePolynomialTrajectory
from Trajectory import NoTrajectoryFound


################### Public interface ######################

class Tunings(object):
    def __init__(self, dt, mvc_dt=None, integ_dt=None, sd_prec=None,
                 switchpoint_steps=20, reparam_dt=None):
        """

        mvc_dt -- time step for discretizing the MVC
        integ_dt -- time step for integrating velocity profiles
        sd_prec -- precision for sdot search around switch points
        switchpoint_steps -- number of time steps to integrate around dynamic
                             singularities
        reparam_dt -- time step for reparametrization

        """
        self.mvc_tstep = mvc_dt if mvc_dt else dt
        self.integ_tstep = integ_dt if integ_dt else dt
        self.reparam_tstep = reparam_dt if reparam_dt else dt
        self.sd_prec = sd_prec if sd_prec else dt
        self.switchpoint_steps = switchpoint_steps

    def __str__(self):
        return "%f %f %f %d %f" % (
            self.mvc_tstep, self.integ_tstep, self.sd_prec,
            self.switchpoint_steps, self.reparam_tstep)


class TorqueConstraints(object):
    def __init__(self, tau_min, tau_max, v_max):
        self.tau_min = tau_min
        self.tau_max = tau_max
        self.v_max = v_max

    def __str__(self):
        tau_min_str = vect2str(self.tau_min)
        tau_max_str = vect2str(self.tau_max)
        #v_max_str = vect2str(self.v_max)
        # zlap* does not work for now:
        v_max_str = vect2str(zeros(len(self.tau_min)))
        return "%s\n%s\n%s" % (tau_min_str, tau_max_str, v_max_str)


class RaveTorqueInstance(object):
    def __init__(self, constraints, openrave_robot, traj, tunings):
        assert isinstance(constraints, TorqueConstraints)
        assert isinstance(traj, PiecewisePolynomialTrajectory)

        self.constraints = constraints
        self.robot = openrave_robot
        self.tunings = tunings
        self.traj = traj

        buffsize = 100000
        input_str = str(constraints) + self.get_dynamics_str()

        assert len(input_str) < buffsize, \
            "%d is bigger than buffer size" % len(input_str)
        assert len(str(self.traj)) < buffsize
        assert len(str(self.tunings)) < buffsize

        # print "Trajectory:"
        # print str(self.traj)
        # print "q(0.0)  =", self.traj.Eval(0)
        # print "qd(0.0) =", self.traj.Evald(0)
        # print ""
        # print "q(0.5)  =", self.traj.Eval(0.5)
        # print "qd(0.5) =", self.traj.Evald(0.5)
        # print ""
        # print "q(1.0)  =", self.traj.Eval(1)
        # print "qd(1.0) =", self.traj.Evald(1)
        # print "--"
        # print ""
        raw_input()

        self.solver = TOPPbindings.TOPPInstance(
            "TorqueLimits", input_str, str(self.traj), str(self.tunings))

    def get_dynamics_str(self):
        dt = self.tunings.mvc_tstep
        nb_steps = 1 + int((self.traj.duration + 1e-10) / dt)
        trange = arange(nb_steps) * dt
        invdyn_str = ''
        for i, t in enumerate(trange):
            q = self.traj.Eval(t)
            qd = self.traj.Evald(t)
            qdd = self.traj.Evaldd(t)
            with self.robot:
                self.robot.SetDOFValues(q)
                self.robot.SetDOFVelocities(qd)
                args, kwargs = (qdd, None), {'returncomponents': True}
                tm, tc, tg = self.robot.ComputeInverseDynamics(*args, **kwargs)
                to = self.robot.ComputeInverseDynamics(qd) - tc - tg
                invdyn_str += '\n' + vect2str(to)       # a vector
                invdyn_str += '\n' + vect2str(tm + tc)  # b vector
                invdyn_str += '\n' + vect2str(tg)       # c vector
        return invdyn_str

    def parametrize_path(self):
        return_code = self.solver.ReparameterizeTrajectory()
        if return_code < 0:
            raise NoTrajectoryFound

        self.solver.WriteResultTrajectory()
        traj_str = self.solver.restrajectorystring
        return PiecewisePolynomialTrajectory.FromString(traj_str)

    def propagate_velocity_interval(self, sd_min, sd_max):
        return_code = self.solver.RunVIP(sd_min, sd_max)
        if return_code == 0:
            raise NoTrajectoryFound

        sd_end_min = self.solver.sdendmin
        sd_end_max = self.solver.sdendmax
        return (sd_end_min, sd_end_max)


###################### Utilities #########################

def vect2str(v):
    return ' '.join(map(str, v))


def Interpolate3rdDegree(q0, q1, qd0, qd1, T):
    a = ((qd1 - qd0) * T - 2 * (q1 - q0 - qd0 * T)) / T ** 3
    b = (3 * (q1 - q0 - qd0 * T) - (qd1 - qd0) * T) / T ** 2
    c = qd0
    d = q0
    return a, b, c, d


def BezierToPolynomial(T, p0, p1, p2, p3):
    a = -p0 + 3 * p1 - 3 * p2 + p3
    b = 3 * p0 - 6 * p1 + 3 * p2
    c = -3 * p0 + 3 * p1
    d = 1
    return a / (T * T * T), b / (T * T), c / T, d


def BezierToTrajectoryString(Tv, p0v, p1v, p2v, p3v):
    nchunks = len(Tv)
    dimension = len(p0v[0])
    trajectorystring = ""
    for i in range(nchunks):
        if i > 0:
            trajectorystring += "\n"
        trajectorystring += str(Tv[i]) + "\n" + str(dimension)
        for j in range(dimension):
            a, b, c, d = BezierToPolynomial(Tv[i], p0v[i][j], p1v[i][j],
                                            p2v[i][j], p3v[i][j])
            trajectorystring += "\n%f %f %f %f" % (d, c, b, a)
    return trajectorystring


################# Reading from string #####################

def ProfileFromLines(lines):
    l = lines[0]
    [duration, dt] = [double(x) for x in l.split(' ')]
    l = lines[1]
    sarray = array([double(x) for x in l.split(' ')])
    l = lines[2]
    sdarray = array([double(x) for x in l.split(' ')])
    return [duration, dt, sarray, sdarray]


def ProfilesFromString(s):
    s = s.strip(" \n")
    profileslist = []
    lines = [l.strip(" \n") for l in s.split('\n')]
    n = len(lines) / 3
    for i in range(n):
        profileslist.append(ProfileFromLines(lines[3 * i:3 * i + 3]))
    return profileslist


def SwitchPointsFromString(s):
    if(len(s))==0:
        return []
    s = s.strip(" \n")
    switchpointslist = []
    lines = [l.strip(" \n") for l in s.split('\n')]
    for l in lines:
        switchpointslist.append(VectorFromString(l))
    return switchpointslist


def VectorFromString(s):
    # left for compatibility TODO: remove?
    s = s.strip(" \n")
    return array([double(x) for x in s.split(' ')])


################# Compute constraints #####################

def ComputeKinematicConstraints(traj, amax, discrtimestep):
    # Sample the dynamics constraints
    ndiscrsteps = int((traj.duration + 1e-10) / discrtimestep) + 1
    constraintstring = ""
    for i in range(ndiscrsteps):
        t = i * discrtimestep
        qd = traj.Evald(t)
        qdd = traj.Evaldd(t)
        constraintstring += "\n" + vect2str(+qd) + " " + vect2str(-qd)
        constraintstring += "\n" + vect2str(+qdd) + " " + vect2str(-qdd)
        constraintstring += "\n" + vect2str(-amax) + " " + vect2str(-amax)
    return constraintstring


######################## Plots ############################

def PlotProfiles(profileslist0, switchpointslist=[], figstart=0):
    profileslist = list(profileslist0)
    figure(figstart)
    clf()
    hold('on')
    mvcbobrow = profileslist.pop(0)
    plot(mvcbobrow[2], mvcbobrow[3], 'm--', linewidth=4)
    mvcdirect = profileslist.pop(0)
    plot(mvcdirect[2], mvcdirect[3], 'c--', linewidth=4)
    for p in profileslist:
        plot(p[2], p[3], linewidth=2)
    if len(profileslist) > 0:
        M = 2 * max([max(p[3]) for p in profileslist])
    else:
        M = 20
        bobrow = filter((lambda x: x < M), mvcbobrow[3])
        direct = filter((lambda x: x < M), mvcdirect[3])
        if len(bobrow) > 0:
            M = max(M, max(bobrow))
        if len(direct) > 0:
            M = max(M, max(direct))
    for sw in switchpointslist:
        if sw[2] == 0:
            plot(sw[0], sw[1], 'ro', markersize=8)
        if sw[2] == 1:
            plot(sw[0], sw[1], 'go', markersize=8)
    axis([0, mvcbobrow[0], 0, M])
    title('MVCs and profiles')


def PlotKinematics(traj0, traj1, dt=0.01, vmax=[], amax=[], figstart=0):
    colorcycle = ['r', 'g', 'b', 'm', 'c', 'y']
    colorcycle = colorcycle[0:traj0.dimension]
    Tmax = max(traj0.duration, traj1.duration)
    # Joint angles
    figure(figstart)
    clf()
    hold('on')
    ax = gca()
    ax.set_color_cycle(colorcycle)
    traj0.Plot(dt, '--')
    traj1.Plot(dt)
    title('Joint values')
    # Velocity
    figure(figstart + 1)
    clf()
    hold('on')
    ax = gca()
    ax.set_color_cycle(colorcycle)
    traj0.Plotd(dt, '--')
    traj1.Plotd(dt)
    for v in vmax:
        plot([0, Tmax], [v, v], '-.')
    for v in vmax:
        plot([0, Tmax], [-v, -v], '-.')
    if len(vmax) > 0:
        Vmax = 1.2 * max(vmax)
        axis([0, Tmax, -Vmax, Vmax])
    title('Joint velocities')
    # Acceleration
    figure(figstart + 2)
    clf()
    ax = gca()
    ax.set_color_cycle(colorcycle)
    hold('on')
    traj0.Plotdd(dt, '--')
    traj1.Plotdd(dt)
    for a in amax:
        plot([0, Tmax], [a, a], '-.')
    for a in amax:
        plot([0, Tmax], [-a, -a], '-.')
    if len(amax) > 0:
        Amax = 1.2 * max(amax)
        axis([0, Tmax, -Amax, Amax])
    title('Joint accelerations')
