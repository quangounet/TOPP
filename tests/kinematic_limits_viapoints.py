# -*- coding: utf-8 -*-
# Copyright (C) 2016 Quang-Cuong Pham <cuong.pham@normalesup.org>
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

import string, time
from pylab import *
from numpy import *
from TOPP import TOPPbindings
from TOPP import TOPPpy
from TOPP import Trajectory
from TOPP import Utilities

# A two-dof path going through 5 viapoints (0,1) - (1,1) - (5,1) - (3,2) - (5,4)
path = array([[0,1,5,3,5],[1,1,1,2,4]])
traj0 = Utilities.InterpolateViapoints(path) # Interpolate using splines

# Constraints
vmax = 2*ones(traj0.dimension)  # Velocity limits
amax = 10*ones(traj0.dimension) # Acceleration limits

# Set up the TOPP instance
trajectorystring = str(traj0)
discrtimestep = 0.005
uselegacy = True
t0 = time.time()
if uselegacy: #Using the legacy KinematicLimits (a bit faster but not fully supported)
    constraintstring = str(discrtimestep)
    constraintstring += "\n" + string.join([str(v) for v in vmax])
    constraintstring += "\n" + string.join([str(a) for a in amax])
    x = TOPPbindings.TOPPInstance(None,"KinematicLimits",constraintstring,trajectorystring);
else: #Using the general QuadraticConstraints (fully supported)
    constraintstring = str(discrtimestep)
    constraintstring += "\n" + string.join([str(v) for v in vmax])
    constraintstring += TOPPpy.ComputeKinematicConstraints(traj0, amax, discrtimestep) 
    x = TOPPbindings.TOPPInstance(None,"QuadraticConstraints",constraintstring,trajectorystring);

# Run TOPP
t1 = time.time()
ret = x.RunComputeProfiles(0,0)
x.ReparameterizeTrajectory()
t2 = time.time()

print "Using legacy:", uselegacy
print "Discretization step:", discrtimestep
print "Setup TOPP:", t1-t0
print "Run TOPP:", t2-t1
print "Total:", t2-t0

# Display results
ion()
x.WriteProfilesList()
x.WriteSwitchPointsList()
profileslist = TOPPpy.ProfilesFromString(x.resprofilesliststring)
switchpointslist = TOPPpy.SwitchPointsFromString(x.switchpointsliststring)
TOPPpy.PlotProfiles(profileslist,switchpointslist,4)
x.WriteResultTrajectory()
traj1 = Trajectory.PiecewisePolynomialTrajectory.FromString(x.restrajectorystring)
dtplot = 0.01
TOPPpy.PlotKinematics(traj0,traj1,dtplot,vmax,amax)

raw_input()
