# -*- coding: utf-8 -*-
# Copyright (C) 2013 Quang-Cuong Pham <cuong.pham@normalesup.org>
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

print "\n************************************\nNB: This test file requires OpenRAVE\n************************************\n"

import string,time,pickle
from pylab import *
from numpy import *
from openravepy import *
from TOPP import TOPPbindings
from TOPP import TOPPpy
from TOPP import TOPPopenravepy
from TOPP import Trajectory
from TOPP import Utilities

ion()


# Robot (OpenRAVE)
env = Environment()
env.Load("robots/barrettwam.robot.xml")
robot=env.GetRobots()[0]
grav=[0,0,-9.8]
ndof=robot.GetDOF()
dof_lim=robot.GetDOFLimits()
vel_lim=robot.GetDOFVelocityLimits()
robot.SetDOFLimits(-100*ones(ndof),100*ones(ndof)) # Overrides robot joint limits for TOPP computations
robot.SetDOFVelocityLimits(100*vel_lim) # Override robot velocity limits for TOPP computations
RaveSetDebugLevel(0)



# Set-up
ntraj = 1000
ncurve = 1
discrtimestep = 0.005
vmax = zeros(ndof)
taumin = zeros(ndof)
taumax = zeros(ndof)
vmax[0:7] = vel_lim[0:7]  # Velocity limits
taumin[0:7] = -robot.GetDOFMaxTorque()[0:7] # Torque limits
taumax[0:7] = robot.GetDOFMaxTorque()[0:7] # Torque limits


# Test begins here
nfail = 0
nsingulartreateds = 0
ntangenttreateds = 0
for j in range(ntraj):
    print j
    p0a = Utilities.vect2str_mintos(rand(ndof)*2*pi-pi)
    p0b = Utilities.vect2str_mintos(rand(ndof)*2*pi-pi)
    p1a = Utilities.vect2str_mintos(rand(ndof)*2*pi-pi)
    p1b = Utilities.vect2str_mintos(rand(ndof)*2*pi-pi)
    s = '%d'%ncurve
    s+= '\n1.0 ' + p0a + ' ' + p0b
    for k in range(ncurve-1):    
        a = rand(ndof)*2*pi-pi
        b = rand(ndof)*2*pi-pi
        c = 2*b-a
        pa = Utilities.vect2str(a)
        pb = Utilities.vect2str(b)
        pc = Utilities.vect2str(c)
        s+= ' ' + pa + ' ' + pb + '\n1.0 ' + pb + ' ' + pc
    s+= ' ' + p1a + ' ' + p1b
    Tv,p0v,p1v,p2v,p3v = TOPPpy.string2p(s)
    trajectorystring = TOPPpy.BezierToTrajectoryString(Tv,p0v,p1v,p2v,p3v)
    traj0 = Trajectory.PiecewisePolynomialTrajectory.FromString(trajectorystring)
    constraintstring = str(discrtimestep)
    constraintstring += "\n" + string.join([str(v) for v in vmax])
    constraintstring += TOPPopenravepy.ComputeTorquesConstraints(robot,traj0,taumin,taumax,discrtimestep)
    x = TOPPbindings.TOPPInstance(None,"QuadraticConstraints",constraintstring,trajectorystring)
    x.extrareps = 5
    ret = x.RunComputeProfiles(0,0)
    x.WriteProfilesList()
    x.WriteSwitchPointsList()
    profileslist = TOPPpy.ProfilesFromString(x.resprofilesliststring)
    switchpointslist = TOPPpy.SwitchPointsFromString(x.switchpointsliststring)
    #TOPPpy.PlotProfiles(profileslist,switchpointslist,4)
    #raw_input()
    if(ret == 1):
        x.ReparameterizeTrajectory()
        x.WriteResultTrajectory()
        traj1 = Trajectory.PiecewisePolynomialTrajectory.FromString(x.restrajectorystring)
        #tvect1, tauvect1 = TOPPopenravepy.ComputeTorques(traj1, robot, 0.01)
        nsingulartreateds += x.nsingulartreated
        ntangenttreateds += x.ntangenttreated
        print x.ntangenttreated, x.nsingulartreated, ntangenttreateds, nsingulartreateds
    else:
        nfail += 1
        print ">>>>>>>>>>>>>>>>>>> TOPP could not retime ", nfail
        TOPPpy.PlotProfiles(profileslist,switchpointslist,4)
        TOPPopenravepy.PlotTorques(robot,traj0,traj1,0.01,taumin,taumax,3)
        raw_input()

      
print "\nNumber of failures:", nfail
print "Number of singularities:", nsingulartreateds      
