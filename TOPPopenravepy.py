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

import string
from numpy import *
from pylab import *
from openravepy import *


def ComputeTorquesConstraints(robot,traj,taumin,taumax,discrtimestep):
    # Sample the dynamics constraints
    ndiscrsteps = int((traj.duration+1e-10)/discrtimestep)+1;
    constraintstring = ""
    for i in range(ndiscrsteps):
        t = i*discrtimestep
        q = traj.Eval(t)
        qd = traj.Evald(t)
        qdd = traj.Evaldd(t)
        with robot:
            robot.SetDOFValues(q)
            robot.SetDOFVelocities(qd)
            tm,tc,tg = robot.ComputeInverseDynamics(qdd,None,returncomponents=True)
            to = robot.ComputeInverseDynamics(qd) - tc - tg
            constraintstring += "\n" + string.join([str(x) for x in to]) + " " + string.join([str(x) for x in -to])
            constraintstring += "\n" + string.join([str(x) for x in tm+tc]) + " " + string.join([str(x) for x in -tm-tc]) 
            constraintstring += "\n" + string.join([str(x) for x in tg-taumax]) + " " + string.join([str(x) for x in -tg+taumin])
    return constraintstring


def ComputeTorquesConstraintsLegacy(robot,traj,taumin,taumax,discrtimestep):
    # Sample the dynamics constraints
    ndiscrsteps = int((traj.duration+1e-10)/discrtimestep)+1;
    constraintstring = ""
    for i in range(ndiscrsteps):
        t = i*discrtimestep
        q = traj.Eval(t)
        qd = traj.Evald(t)
        qdd = traj.Evaldd(t)
        with robot:
            robot.SetDOFValues(q)
            robot.SetDOFVelocities(qd)
            tm,tc,tg = robot.ComputeInverseDynamics(qdd,None,returncomponents=True)
            to = robot.ComputeInverseDynamics(qd) - tc - tg
            constraintstring += "\n" + string.join([str(x) for x in to])
            constraintstring += "\n" + string.join([str(x) for x in tm+tc])
            constraintstring += "\n" + string.join([str(x) for x in tg])
    return constraintstring


def PlotTorques(robot,traj0,traj1,dt=0.001,taumin=[],taumax=[],figstart=0):
    colorcycle = ['r','g','b','m','c','y']
    colorcycle = colorcycle[0:traj0.dimension]
    Tmax = max(traj0.duration,traj1.duration)
    tvect0,tauvect0 = ComputeTorques(traj0,robot,dt)
    tvect1,tauvect1 = ComputeTorques(traj1,robot,dt)
    figure(figstart)
    clf()
    hold('on')
    ax=gca()        
    ax.set_color_cycle(colorcycle)
    plot(tvect0,tauvect0,'--',linewidth=2)
    plot(tvect1,tauvect1,linewidth=2)
    for a in taumax:
        plot([0,Tmax],[a,a],'-.')
    for a in taumin:
        plot([0,Tmax],[a,a],'-.')
    if(len(taumax)>0):
        axis([0,Tmax,1.2*min(taumin),1.2*max(taumax)])
    title('Joint torques')


def ComputeTorques(traj,robot,dt):
    tvect = arange(0,traj.duration+dt,dt)
    tauvect = []
    for t in tvect:
        with robot:
            q = traj.Eval(t)
            qd = traj.Evald(t)
            qdd = traj.Evaldd(t)
            robot.SetDOFValues(q)
            robot.SetDOFVelocities(qd)
            tau = robot.ComputeInverseDynamics(qdd,None,returncomponents=False)
            tauvect.append(tau)
    return tvect,array(tauvect)

