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


import TOPPbindings
import TOPPpy
import time
import string
import sys
from pylab import *
from numpy import *


ion()

############################ Tunings ############################
discrtimestep = 0.01;
integrationtimestep = 0.01;
bisectionprecision = 0.01;
passswitchpointnsteps = 20;
reparamtimestep = 0.01;
tuningsstring = "%f %f %f %d %f"%(discrtimestep,integrationtimestep,bisectionprecision,passswitchpointnsteps,reparamtimestep);


############################ Trajectory ############################
#------------------------------------------#
trajectorystring = "2 \n 2\n 1 1 0 1\n 0 2 0 -1\n 3\n 2\n 11 13 6 0.1666666666666\n -4 -10 -6 0.5";
#------------------------------------------#
traj0 = TOPPpy.PiecewisePolynomialTrajectory.FromString(trajectorystring)


############################ Constraints ############################
#------------------------------------------#
amax = array([15,10])
vmax = array([20,10])
t0 = time.time()
constraintstring = string.join([str(v) for v in vmax])
constraintstring += TOPPpy.ComputeKinematicConstraints(traj0,amax,discrtimestep)
t1 = time.time()
#------------------------------------------#


############################ Run TOPP ############################
x = TOPPbindings.TOPPInstance("QuadraticConstraints",constraintstring,trajectorystring,tuningsstring);
t2 = time.time()
ret = x.RunPP(0,0)
t3 = time.time()

################ Plotting the MVC and the profiles #################
x.WriteProfilesList()
profileslist = TOPPpy.ProfilesFromString(x.resprofilesliststring)
TOPPpy.PlotProfiles(profileslist,4)


##################### Plotting the trajectories #####################
x.WriteResultTrajectory()
traj1 = TOPPpy.PiecewisePolynomialTrajectory.FromString(x.restrajectorystring)
dtplot = 0.01
TOPPpy.PlotKinematics(traj0,traj1,dtplot,vmax,amax)


print "\n--------------"
print "Sampling dynamics: ", t1-t0
print "Building TOPP Instance: ", t2-t1
print "TOPP proper (C++): ", t3-t2
print "Total: ", t3-t0 




raw_input()
