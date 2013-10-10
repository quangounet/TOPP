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
import TOPPopenravepy
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
traj0 = TOPPpy.PiecewisePolynomialTrajectory(trajectorystring)


############################ Constraints ############################
#------------------------------------------#
amax = array([15,10])
vmax = array([20,10])
constraintstring = string.join([str(v) for v in amax]) + "\n"
constraintstring += string.join([str(v) for v in vmax])
#------------------------------------------#


############################ Run TOPP ############################
x = TOPPbindings.TOPPInstance("KinematicLimits",constraintstring,trajectorystring,tuningsstring);
ret = x.RunPP(0,0)


################ Plotting the MVC and the profiles #################
x.WriteProfilesList()
profileslist = TOPPpy.ProfilesFromString(x.resprofilesliststring)
TOPPpy.PlotProfiles(profileslist,4)


##################### Plotting the trajectories #####################
x.WriteResultTrajectory()
traj1 = TOPPpy.PiecewisePolynomialTrajectory(x.restrajectorystring)
dtplot = 0.01
TOPPpy.PlotKinematics(traj0,traj1,dtplot,vmax,amax)



raw_input()
