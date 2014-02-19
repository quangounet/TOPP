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

import string
from pylab import *
from numpy import *
import TOPPbindings
import TOPPpy

# Trajectory
trajectorystring = "2 \n 2\n 1 1 0 1\n 0 2 0 -1\n 3\n 2\n 11 13 6 0.1666666666666\n -4 -10 -6 0.5"

# Constraints
discrtimestep = 0.01
vmax = array([10,10]) # Velocity limits
amax = array([15,10]) # Acceleration limits
constraintstring = str(discrtimestep) + "\n";  # Discretization time step
constraintstring += string.join([str(v) for v in vmax]) + "\n" 
constraintstring += string.join([str(a) for a in amax])

# Run TOPP
x = TOPPbindings.TOPPInstance(None,"KinematicLimits",constraintstring,trajectorystring);
ret = x.RunComputeProfiles(0,0)
x.ReparameterizeTrajectory()

# Display results
ion()
x.WriteProfilesList()
x.WriteSwitchPointsList()
profileslist = TOPPpy.ProfilesFromString(x.resprofilesliststring)
switchpointslist = TOPPpy.SwitchPointsFromString(x.switchpointsliststring)
TOPPpy.PlotProfiles(profileslist,switchpointslist,4)
x.WriteResultTrajectory()
traj1 = TOPPpy.PiecewisePolynomialTrajectory.FromString(x.restrajectorystring)
traj0 = TOPPpy.PiecewisePolynomialTrajectory.FromString(trajectorystring)
dtplot = 0.01
TOPPpy.PlotKinematics(traj0,traj1,dtplot,vmax,amax)
print "Trajectory duration before TOPP: ", traj0.duration
print "Trajectory duration after TOPP: ", traj1.duration
