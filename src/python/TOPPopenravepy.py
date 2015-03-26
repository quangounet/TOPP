# - * - coding: utf - 8 - * -
# Copyright (C) 2013 Quang - Cuong Pham <cuong.pham@normalesup.org>
#
# This file is part of the Time-Optimal Path Parameterization (TOPP) library.
# TOPP is free software: you can redistribute it and/or modify
# it under the terms of the GNU Lesser General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# at your option, any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU Lesser General Public License for more details.
#
# You should have received a copy of the GNU Lesser General Public License
# along with this program.  If not, see <http://www.gnu.org/licenses/>.

import openravepy
import string
import time

import Trajectory
import TOPPbindings

from QuadraticConstraints import QuadraticConstraints
from pylab import arange, array, cross, dot, inv, norm, random, zeros
from pylab import arcsin, arctan2, cos, pi, sin


class RAVEBindings(QuadraticConstraints):
    """Wrapper around TOPPbindings for OpenRAVE robot models.

    robot -- OpenRAVE robot object
    pbname -- string identifier of the underlying TOPPbindings problem
    constring -- constraints string (problem - dependant)
    trajstring -- trajectory string

    """

    def __init__(self, robot, pbname, constring, trajstring,
                 discrtimestep=None, integrationtimestep=None):
        if discrtimestep is None:
            discrtimestep = 1e-2
        if integrationtimestep is None:
            integrationtimestep = 1e-3
        self.discrtimestep = discrtimestep
        self.integrationtimestep = integrationtimestep
        self.solver = TOPPbindings.TOPPInstance(robot, pbname, constring,
                                                trajstring)


def ToRaveTraj(robot, spec, traj):
    nchunks = len(traj.chunkslist)
    if nchunks < 1:
        print "TOPP trajectory has less than 1 chunk"
        return None
    timespec = openravepy.ConfigurationSpecification()
    timespec.AddGroup('deltatime', 1, 'linear')
    posspec = robot.GetActiveConfigurationSpecification()
    velspec = posspec.ConvertToVelocitySpecification()
    ravetraj = openravepy.RaveCreateTrajectory(robot.GetEnv(), '')
    ravetraj.Init(spec)
    ravetraj.Insert(0, array([0]), timespec)
    ravetraj.Insert(0, traj.Eval(0), posspec, True)
    ravetraj.Insert(0, traj.Evald(0), velspec, True)
    for i in range(nchunks):
        c = traj.chunkslist[i]
        ravetraj.Insert(i + 1, array([c.duration]), timespec)
        ravetraj.Insert(i + 1, c.Eval(c.duration), posspec, True)
        ravetraj.Insert(i + 1, c.Evald(c.duration), velspec, True)
    return ravetraj


def FromRaveTraj(robot, traj):
    N = traj.GetNumWaypoints()
    if N < 2:
        print "RAVE trajectory has less than 2 waypoints"
        return None
    timespec = openravepy.ConfigurationSpecification()
    timespec.AddGroup('deltatime', 1, 'linear')
    posspec = robot.GetActiveConfigurationSpecification()
    velspec = posspec.ConvertToVelocitySpecification()
    chunkslist = []
    vel = traj.GetWaypoint(0, velspec)
    ndof = len(vel)
    for i in range(N - 1):
        pos = traj.GetWaypoint(i, posspec)
        deltatime = traj.GetWaypoint(i + 1, timespec)[0]
        if deltatime < 1e-5:
            continue
        nextvel = traj.GetWaypoint(i + 1, velspec)
        polynomialslist = []
        for j in range(ndof):
            x = pos[j]
            v = vel[j]
            a = (nextvel[j] - vel[j]) / deltatime
            polynomialslist.append(Trajectory.Polynomial([x, v, a / 2]))
        chunkslist.append(Trajectory.Chunk(deltatime, polynomialslist))
        vel = nextvel
    return Trajectory.PiecewisePolynomialTrajectory(chunkslist)


def ComputeTorquesConstraints(robot, traj, taumin, taumax, discrtimestep):
    """Sample the dynamics constraints."""
    ndiscrsteps = int((traj.duration + 1e-10) / discrtimestep) + 1
    constraintstring = ""
    for i in range(ndiscrsteps):
        t = i * discrtimestep
        q = traj.Eval(t)
        qd = traj.Evald(t)
        qd2 = zeros(robot.GetDOF())
        qd2[robot.GetActiveDOFIndices()] = qd
        qdd = traj.Evaldd(t)
        qdd2 = zeros(robot.GetDOF())
        qdd2[robot.GetActiveDOFIndices()] = qdd
        with robot:
            robot.SetActiveDOFValues(q)
            robot.SetActiveDOFVelocities(qd)
            tm, tc, tg = robot.ComputeInverseDynamics(qdd2, None,
                                                      returncomponents=True)
            to = robot.ComputeInverseDynamics(qd2) - tc - tg
            avect = []
            bvect = []
            cvect = []
            for i in range(len(taumax)):
                if abs(taumax[i]) > 1e-10:
                    avect.append(to[i])
                    bvect.append(tm[i] + tc[i])
                    cvect.append(tg[i] - taumax[i])
            for i in range(len(taumin)):
                if abs(taumin[i]) > 1e-10:
                    avect.append(-to[i])
                    bvect.append(-tm[i] - tc[i])
                    cvect.append(-tg[i] + taumin[i])
            constraintstring += "\n" + string.join([str(a) for a in avect])
            constraintstring += "\n" + string.join([str(b) for b in bvect])
            constraintstring += "\n" + string.join([str(c) for c in cvect])
    return constraintstring


def PlotTorques(robot, traj0, traj1, dt=0.001, taumin=[], taumax=[],
                figstart=0):
    from pylab import figure, clf, hold, gca, plot, axis, title, xlabel, ylabel
    colorcycle = ['r', 'g', 'b', 'm', 'c', 'y', 'k']
    colorcycle = colorcycle[0:traj0.dimension]
    Tmax = max(traj0.duration, traj1.duration)
    tvect0, tauvect0 = ComputeTorques(traj0, robot, dt)
    tvect1, tauvect1 = ComputeTorques(traj1, robot, dt)
    figure(figstart)
    clf()
    hold('on')
    ax = gca()
    ax.set_color_cycle(colorcycle)
    plot(tvect0, tauvect0, '--', linewidth=2)
    ax.set_color_cycle(colorcycle)
    plot(tvect1, tauvect1, linewidth=2)
    ax.set_color_cycle(colorcycle)
    for a in taumax:
        plot([0, Tmax], [a, a], '-.')
    ax.set_color_cycle(colorcycle)
    for a in taumin:
        plot([0, Tmax], [a, a], '-.')
    if len(taumax) > 0:
        axis([0, Tmax, 1.2 * min(taumin), 1.2 * max(taumax)])
    title('Joint torques', fontsize=20)
    xlabel('Time (s)', fontsize=18)
    ylabel('Joint torques (Nm)', fontsize=18)


def PlotTorques2(tvect0, torques0, tvect1, torques1, taumin=[], taumax=[],
                 figstart=0):
    from pylab import figure, clf, hold, gca, plot, axis, title, xlabel, ylabel
    colorcycle = ['r', 'g', 'b', 'm', 'c', 'y', 'k']
    Tmax = max(max(tvect0), max(tvect1))
    figure(figstart)
    clf()
    hold('on')
    ax = gca()
    ax.set_color_cycle(colorcycle)
    plot(tvect0, torques0, '--', linewidth=2)
    ax.set_color_cycle(colorcycle)
    plot(tvect1, torques1, linewidth=2)
    for a in taumax:
        plot([0, Tmax], [a, a], '-.')
    for a in taumin:
        plot([0, Tmax], [a, a], '-.')
    if len(taumax) > 0:
        axis([0, Tmax, 1.2 * min(taumin), 1.2 * max(taumax)])
    title('Joint torques', fontsize=20)
    xlabel('Time (s)', fontsize=18)
    ylabel('Joint torques (Nm)', fontsize=18)


def PlotZMP(robot, traj0, traj1, zmplimits, dt=0.01, figstart=0, border=0):
    from pylab import figure, clf, hold, plot, axis, title, xlabel, ylabel
    xmin, xmax, ymin, ymax = zmplimits
    xminf, xmaxf = xmin - border, xmax + border
    yminf, ymaxf = ymin - border, ymax + border
    tvect0, xzmp0, yzmp0, com0 = ComputeZMP(traj0, robot, dt)
    tvect1, xzmp1, yzmp1, com1 = ComputeZMP(traj1, robot, dt)
    com0, com1 = array(com0), array(com1)
    figure(figstart)
    clf()
    hold('on')
    plot(tvect0, xzmp0, 'r--', linewidth=2)
    plot(tvect0, yzmp0, 'g--', linewidth=2)
    plot(tvect1, xzmp1, 'r', linewidth=2)
    plot(tvect1, yzmp1, 'g', linewidth=2)
    plot(tvect0, com0[:, 0], 'm', linewidth=2)
    plot(tvect0, com0[:, 1], 'c', linewidth=2)
    Tmax = max(traj0.duration, traj1.duration)
    plot([0, Tmax], [xmin, xmin], 'r-.')
    plot([0, Tmax], [xmax, xmax], 'r-.')
    plot([0, Tmax], [ymin, ymin], 'g-.')
    plot([0, Tmax], [ymax, ymax], 'g-.')
    axis([0, Tmax, 1.2 * min(xmin, ymin), 1.2 * max(xmax, ymax)])
    title('Coordinates of the ZMP', fontsize=20)
    xlabel('Time (s)', fontsize=18)
    ylabel('ZMP (m)', fontsize=18)
    figure(figstart + 1)
    clf()
    plot([xminf, xminf, xmaxf, xmaxf, xminf], [yminf, ymaxf, ymaxf, yminf,
                                               yminf], 'k', linewidth=2)
    plot([xmin, xmin, xmax, xmax, xmin], [ymin, ymax, ymax, ymin, ymin], 'k--',
         linewidth=2)
    plot(xzmp0, yzmp0, 'r--', linewidth=3)
    plot(xzmp1, yzmp1, 'r', linewidth=3)
    plot(xzmp0[0], yzmp0[0], 'rs', markersize=14)
    plot(xzmp1[0], yzmp1[0], 'rs', markersize=15)
    plot(com1[0, 0], com1[0, 1], 'bs', markersize=8)
    plot(xzmp0[-1], yzmp0[-1], 'r*', markersize=25)
    plot(xzmp1[-1], yzmp1[-1], 'r*', markersize=25)
    plot(com1[:, 0], com1[:, 1], 'b', linewidth=3)
    plot(com1[-1, 0], com1[-1, 1], 'b*', markersize=15)
    axis([xminf - 0.01, xmax + 0.01, ymin - 0.01, ymax + 0.01])
    axis('equal')
    title('Spatial trajectory of the ZMP under the left foot', fontsize=20)
    xlabel('Anteroposterior axis (m)', fontsize=18)
    ylabel('Mediolateral axis (m)', fontsize=18)


def Fill(robot, q):
    if hasattr(robot, 'activedofs'):
        n = robot.GetDOF()
        qfilled = robot.qdefault
        counter = 0
        for i in range(n):
            if robot.activedofs[i] > 0.1:
                qfilled[i] = q[counter]
                counter += 1
        return qfilled
    else:
        return q


def Trim(robot, q):
    if hasattr(robot, 'activedofs'):
        n = robot.GetDOF()
        qtrimmed = []
        for i in range(n):
            if robot.activedofs[i] > 0.1:
                qtrimmed.append(q[i])
        return qtrimmed
    else:
        return q


def Execute(robot, traj, dt=0.01):
    tvect = arange(0, traj.duration + dt, dt)
    for t in tvect:
        q = traj.Eval(t)
        robot.SetDOFValues(Fill(robot, q))
        time.sleep(dt)


def ComputeTorques(traj, robot, dt):
    tvect = arange(0, traj.duration + dt, dt)
    tauvect = []
    for t in tvect:
        with robot:
            q = traj.Eval(t)
            qd = traj.Evald(t)
            qdd = traj.Evaldd(t)
            qdd2 = zeros(robot.GetDOF())
            qdd2[robot.GetActiveDOFIndices()] = qdd
            robot.SetActiveDOFValues(Fill(robot, q))
            robot.SetActiveDOFVelocities(Fill(robot, qd))
            tau = robot.ComputeInverseDynamics(Fill(robot, qdd2), None,
                                               returncomponents=False)
            tauvect.append(Trim(robot, tau))
    return tvect, array(tauvect)


def CheckCollisionStaticStabilityConfig(robot, q, baselimits):
    """Return true if in collision or not quasi - statically stable."""
    with robot:
        robot.SetDOFValues(Fill(robot, q))
        n = len(robot.GetLinks())
        totalmass = 0
        com = zeros(3)
        for i in range(n):
            if hasattr(robot, 'activelinks') and robot.activelinks[i] < 0.1:
                continue
            ci = robot.GetLinks()[i].GetGlobalCOM()
            mass = robot.GetLinks()[i].GetMass()
            com += mass * ci
            totalmass += mass
        com = com / totalmass
        if com[0] < baselimits[0] or com[0] > baselimits[1] \
                or com[1] < baselimits[2] or com[1] > baselimits[3]:
            return True
        return robot.GetEnv().CheckCollision(robot) \
            or robot.CheckSelfCollision()


def CheckCollisionStaticStabilitySegment(robot, qstart, qend, steplength,
                                         baselimits):
    """
    Return true if segment is in collision or not quasi - statically stable.
    Checks discretized configurations starting from qstart.
    """
    n = int(norm(qend - qstart) / steplength) + 1
    dq = (qend - qstart) / n
    for i in range(n + 1):
        q = qstart + i * dq
        if CheckCollisionStaticStabilityConfig(robot, q, baselimits):
            return True
    return False


def ComputeZMPConfig(robot, q, qd, qdd):
    n = len(robot.GetLinks())
    g = robot.GetEnv().GetPhysicsEngine().GetGravity()
    with robot:
        robot.SetDOFValues(Fill(robot, q))
        robot.SetDOFVelocities(Fill(robot, qd))
        com_pos = array([k.GetGlobalCOM() for k in robot.GetLinks()])
        vel = robot.GetLinkVelocities()
        acc = robot.GetLinkAccelerations(Fill(robot, qdd))
        transforms = [k.GetTransform()[0:3, 0:3] for k in robot.GetLinks()]
        masses = [k.GetMass() for k in robot.GetLinks()]
        localCOM = [k.GetLocalCOM() for k in robot.GetLinks()]
    tau0 = array([0., 0., 0.])
    totalmass = sum(masses)
    if hasattr(robot, 'activelinks'):
        totalmass = 0
        for k in range(n):
            if robot.activelinks[k] > 0.1:
                totalmass += masses[k]
    f02 = totalmass * g[2]
    com = zeros(3)
    for i in range(n):
        if hasattr(robot, 'activelinks') and robot.activelinks[i] < 0.1:
            continue
        # Compute the inertia matrix in the global frame
        R = transforms[i]
        ri = dot(R, localCOM[i])
        omegai = vel[i, 3:6]
        omegadi = acc[i, 3:6]
        #com_vel = vel[i, 0:3] + cross(omegai, ri)
        ci = com_pos[i]
        cidd = acc[i, 0:3] + cross(omegai, cross(omegai, ri)) + cross(omegadi,
                                                                      ri)
        tau0 += masses[i] * cross(ci, g - cidd)
        f02 -= masses[i] * cidd[2]
        com += masses[i] * ci

    return -tau0[1] / f02, tau0[0] / f02, com / totalmass


def ComputeZMP(traj, robot, dt):
    tvect = arange(0, traj.duration + dt, dt)
    xzmp = []
    yzmp = []
    com = []
    for t in tvect:
        q = traj.Eval(t)
        qd = traj.Evald(t)
        qdd = traj.Evaldd(t)
        x, y, comq = ComputeZMPConfig(robot, q, qd, qdd)
        xzmp.append(x)
        yzmp.append(y)
        com.append(comq)
    return tvect, xzmp, yzmp, com


class HRP4Robot():
    def __init__(self, robot, baselimits, collisioncheckstep):
        self.robot = robot
        self.baselimits = baselimits
        self.collisioncheckstep = collisioncheckstep
        self.qlowerlimit = []
        self.qupperlimit = []
        n = robot.GetDOF()
        for i in range(n):
            if hasattr(robot, 'activedofs'):
                if robot.activedofs[i] > 0.1:
                    self.qlowerlimit.append(robot.GetDOFLimits()[0][i])
                    self.qupperlimit.append(robot.GetDOFLimits()[1][i])
            else:
                self.qlowerlimit.append(robot.GetDOFLimits()[0][i])
                self.qupperlimit.append(robot.GetDOFLimits()[1][i])
                self.weights.append(1)
        self.weights = abs(array(self.qupperlimit) - array(self.qlowerlimit))
        self.weights = self.weights / norm(self.weights)

    def CheckCollisionSegment(self, qstart, qend):
        return CheckCollisionStaticStabilitySegment(self.robot, qstart, qend,
                                                    self.collisioncheckstep,
                                                    self.baselimits)

    def CheckCollisionConfig(self, q):
        return CheckCollisionStaticStabilityConfig(self.robot, q,
                                                   self.baselimits)

    def RandomConfig(self):
        q = zeros(len(self.qlowerlimit))
        for i in range(len(self.qlowerlimit)):
            deltaq = self.qupperlimit[i] - self.qlowerlimit[i]
            q[i] = self.qlowerlimit[i] + random() * deltaq
        return q


def RotFromAngles(a):
    """See Tait - Bryan X1 - Y2 - Z3 on Wikipedia."""
    h1 = a[0]
    h2 = a[1]
    h3 = a[2]
    c1 = cos(h1)
    s1 = sin(h1)
    c2 = cos(h2)
    s2 = sin(h2)
    c3 = cos(h3)
    s3 = sin(h3)
    return array([
        [c2 * c3, - c2 * s3, s2],
        [c1 * s3 + c3 * s1 * s2, c1 * c3 - s1 * s2 * s3, - c2 * s1],
        [s1 * s3 - c1 * c3 * s2, c3 * s1 + c1 * s2 * s3, c1 * c2]])


def AnglesFromRot(R):
    # Hypothesis 1
    h2a = arcsin(R[0, 2])
    cos2 = cos(h2a)
    if abs(cos2) < 1e-8:
        cos2 += 1e-8
    h1a = arctan2(-R[1, 2] / cos2, R[2, 2] / cos2)
    h3a = arctan2(-R[0, 1] / cos2, R[0, 0] / cos2)
    # Hypothesis 2
    h2b = pi - h2a
    cos2 = cos(h2b)
    if abs(cos2) < 1e-8:
        cos2 += 1e-8
    h1b = arctan2(-R[1, 2] / cos2, R[2, 2] / cos2)
    h3b = arctan2(-R[0, 1] / cos2, R[0, 0] / cos2)
    # Testing
    resa = abs(sum(R - RotFromAngles([h1a, h2a, h3a])))
    resb = abs(sum(R - RotFromAngles([h1b, h2b, h3b])))
    if resa < resb:
        return [h1a, h2a, h3a]
    else:
        return [h1b, h2b, h3b]


def JointValuesFromTransform(robot, Tdesired):
    """
    Find the dummy joint values such that the Transform of the baselink is
    Tdesired. T0 is the initial Transform of the baselink.
    """
    R = dot(inv(robot.baselinkinittransform[0:3, 0:3]), Tdesired[0:3, 0:3])
    [h1, h2, h3] = AnglesFromRot(R)
    offset = dot(R, robot.baselinkinittransform[0:3, 3])
    [s1, s2, s3] = Tdesired[0:3, 3] - offset
    return [s1, s2, s3, h1, h2, h3]


def LoadFloat(env, robotfile, baselinkname):
    xml = """
<robot>
  <kinbody>
    <body name="root">
    <mass type="mimicgeom">
      <total>0</total>
    </mass>
    </body>
  </kinbody>
  <kinbody>
    <body name="dummy1">
    <mass type="mimicgeom">
      <total>0</total>
    </mass>
    </body>
  </kinbody>
  <kinbody>
    <body name="dummy2">
    <mass type="mimicgeom">
      <total>0</total>
    </mass>
    </body>
  </kinbody>
  <kinbody>
    <body name="dummy3">
    <mass type="mimicgeom">
      <total>0</total>
    </mass>
    </body>
  </kinbody>
  <kinbody>
    <body name="dummy4">
    <mass type="mimicgeom">
      <total>0</total>
    </mass>
    </body>
  </kinbody>
  <kinbody>
    <body name="dummy5">
    <mass type="mimicgeom">
      <total>0</total>
    </mass>
    </body>
  </kinbody>
  <robot file="%s">
    <kinbody>
      <body name="%s">
      </body>
      <joint name="slider1" type="slider" circular="true">
        <body>root</body>
        <body>dummy1</body>
        <axis>1 0 0</axis>
      </joint>
      <joint name="slider2" type="slider" circular="true">
        <body>dummy1</body>
        <body>dummy2</body>
        <axis>0 1 0</axis>
      </joint>
      <joint name="slider3" type="slider" circular="true">
        <body>dummy2</body>
        <body>dummy3</body>
        <axis>0 0 1</axis>
      </joint>
      <joint name="hinge1" type="hinge" circular="true">
        <body>dummy3</body>
        <body>dummy4</body>
        <axis>1 0  0</axis>
      </joint>
      <joint name="hinge2" type="hinge" circular="true">
        <body>dummy4</body>
        <body>dummy5</body>
        <axis>0 1 0</axis>
      </joint>
      <joint name="hinge3" type="hinge" circular="true">
        <body>dummy5</body>
        <body>%s</body>
        <axis>0 0 1</axis>
      </joint>
    </kinbody>
  </robot>
</robot>
""" % (robotfile, baselinkname, baselinkname)
    robot = env.ReadRobotData(xml)
    env.Add(robot)
    baselink = robot.GetLink(baselinkname)
    robot.baselinkinittransform = baselink.GetTransform()
    return robot
