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

import sys
sys.path.append('..')

import os

import TOPPbindings
import TOPPpy
import time
import string
from pylab import *
from numpy import *


plotting = True

ion()

############################ Tunings ############################

integrationtimestep = 0 # auto
reparamtimestep = 0 # auto
passswitchpointnsteps = 0
discrtimestep = 1/100.
tuningsstring = "%f %f %f %d"%(discrtimestep,integrationtimestep,reparamtimestep,passswitchpointnsteps)

############################ Constraints ############################
ndof = 30
v = 1.5
a = 1
j = 28
vmax = v*ones(ndof)
amax = a*ones(ndof)
constraintstring = string.join([str(v) for v in amax]) + "\n"
constraintstring += string.join([str(v) for v in vmax])


############################ Trajectory ############################

trajfile = 'testfiles/traj-%d-%d'%(ndof,j)
h = open(trajfile,'r')
s = h.read()
h.close()      
Tv,p0v,p1v,p2v,p3v = TOPPpy.string2p(s)
trajectorystring = TOPPpy.BezierToTrajectoryString(Tv,p0v,p1v,p2v,p3v)

# ############################ Mintos ############################
# limitfile = 'testfiles/limits-%d-%f-%f'%(ndof,v,a)
# command = "./timeopt %s %s %d 1 > /tmp/res-%d-%f"%(trajfile,limitfile,gridres2,ndof,v)
# os.system(command)
# res = open("/tmp/res-%d-%f"%(ndof,v),"r").read()
# lines = [l.strip(" \n") for l in res.split('\n')]
# resline = ""
# for l in lines:
#     if(len(l)>0 and l[0]=="O"):
#         resline = l
#         break
# print "Mintos traj duration: ", (float(resline.split(" ")[2].strip(",")))
# print "Mintos comput time: ", (float(resline.split(" ")[4]))


############################ Run TOPP ############################
t0 = time.time()
x = TOPPbindings.TOPPInstance("KinematicLimits",constraintstring,trajectorystring,tuningsstring,False);
t1 = time.time()
ret = x.RunComputeProfiles(0,0)
t2 = time.time()

if(ret == 1):
    x.ReparameterizeTrajectory()

t3 = time.time()


################ Plotting the MVC and the profiles #################
x.WriteProfilesList()
x.WriteSwitchPointsList()
profileslist = TOPPpy.ProfilesFromString(x.resprofilesliststring)
switchpointslist = TOPPpy.SwitchPointsFromString(x.switchpointsliststring)
if plotting:
    TOPPpy.PlotProfiles(profileslist,switchpointslist,4)
    axis([0,sum(Tv),0,1])


##################### Plotting the trajectories #####################
if(ret == 1):
    x.WriteResultTrajectory()
    traj1 = TOPPpy.PiecewisePolynomialTrajectory.FromString(x.restrajectorystring)
    if plotting:
        TOPPpy.PlotKinematics(traj1,traj1,0.01,vmax,amax)

print "TOPP traj duration: ", x.resduration
print "TOPP comput time: ", t2-t0


raw_input()
