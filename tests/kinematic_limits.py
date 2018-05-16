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

import time
from pylab import *
from numpy import *
from TOPP import TOPPbindings
from TOPP import TOPPpy
from TOPP import Trajectory

# Trajectory
ndof = 5
trajectorystring = """1.0
5
-0.495010 1.748820 -2.857899 1.965396
0.008319 0.004494 1.357524 -1.289918
-0.354081 1.801074 -1.820616 0.560944
0.221734 -1.320792 3.297177 -2.669786
-0.137741 0.105246 0.118968 -0.051712
1.0
5
0.361307 1.929207 -4.349490 2.275776
0.080419 -1.150212 2.511645 -1.835906
0.187321 -0.157326 -0.355785 0.111770
-0.471667 -2.735793 7.490559 -4.501124
0.034761 0.188049 -1.298730 1.553443"""
traj0 = Trajectory.PiecewisePolynomialTrajectory.FromString(trajectorystring)

# Constraints
vmax = 2*ones(ndof)  # Velocity limits
amax = 10*ones(ndof) # Acceleration limits

# Set up the TOPP instance
discrtimestep = 0.01
uselegacy = True
t0 = time.time()
if uselegacy: #Using the legacy KinematicLimits (a bit faster but not fully supported)
    constraintstring = str(discrtimestep)
    constraintstring += "\n" + " ".join([str(v) for v in vmax])
    constraintstring += "\n" + " ".join([str(a) for a in amax])
    x = TOPPbindings.TOPPInstance(None,"KinematicLimits",constraintstring,trajectorystring);
else: #Using the general QuadraticConstraints (fully supported)
    constraintstring = str(discrtimestep)
    constraintstring += "\n" + " ".join([str(v) for v in vmax])
    constraintstring += TOPPpy.ComputeKinematicConstraints(traj0, amax, discrtimestep) 
    x = TOPPbindings.TOPPInstance(None,"QuadraticConstraints",constraintstring,trajectorystring);

# Run TOPP
t1 = time.time()
ret = x.RunComputeProfiles(0,0)
x.ReparameterizeTrajectory()
t2 = time.time()

print("Using legacy:", uselegacy)
print("Discretization step:", discrtimestep)
print("Setup TOPP:", t1-t0)
print("Run TOPP:", t2-t1)
print("Total:", t2-t0)

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

input()
