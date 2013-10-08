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


import TOPPbindings
import TOPPpy
import time
import sys

from pylab import arange, array, axis, clf, figure, hold, ion, plot

# Constraints :
# amax0, amax1 \n vmax0, vmax1 (vmax = 0 means no velocity constraints)
amax0 = 15
amax1 = 10
vmax0 = 20
vmax1 = 10
constraintstring = "%f %f\n %f %f" % (amax0, amax1, vmax0, vmax1)


# Tunings :
# - time step for discretizing the MVC
# - time step for integrating the profiles
# - precision for sdot search around switch points
# - number of time steps to integrate around dynamic singularities
# - time step for reparameterization
tuningsstring = "0.01 0.01 0.01 20 0.01"


# Trajectory :
# Chunk 1 time duration
# Chunk 1 number of dofs
# Coefficients of the polynomial of dof 1 (from lower to higher)
# Coefficients of the polynomial of dof 2 (from lower to higher) ...
# Chunk 2 ...
trajectorystring = \
    "2 \n 2\n 1 1 0 1\n 0 2 0 -1\n " + \
    "3\n 2\n 11 13 6 0.1666666666666\n -4 -10 -6 0.5"


# Run TOPP
start = time.time()
x = TOPPbindings.TOPPInstance("KinematicLimits", constraintstring,
                              trajectorystring, tuningsstring)
ret = x.RunPP(1e-4, 1e-4)

if ret == 0:
    sys.exit()

print "Computation time: ", (time.time() - start)
print "Duration reparameterized trajectory: ", x.resduration

# Run Velocity Interval Propagation
#x.RunVIP(1,2)
#print "[sdendmin,sdendmax] = [", x.sdendmin ,",", x.sdendmax, "]"


# Display results
x.WriteResultTrajectory()
traj0 = TOPPpy.PiecewisePolynomialTrajectory.FromString(trajectorystring)
traj1 = TOPPpy.PiecewisePolynomialTrajectory.FromString(x.restrajectorystring)


# Verification
ion()
dt = 0.1
tvect = arange(0, traj1.duration + dt, dt)
qdd = array([traj1.Evaldd(t) for t in tvect])
print "Max acceleration: ", max(abs(qdd[:, 0])), ",", max(abs(qdd[:, 1]))
Tmax = max(traj0.duration, traj1.duration)
Vmax = 1.2 * max(vmax0, vmax1)
Amax = 1.2 * max(amax0, amax1)

# Plotting
figure(0)
clf()
hold('on')
traj0.Plot(dt)
traj1.Plot(dt, '--')

figure(1)
clf()
hold('on')
traj0.Plotd(dt)
traj1.Plotd(dt, '--')
plot([0, Tmax], [vmax0, vmax0], 'r-.')
plot([0, Tmax], [-vmax0, -vmax0], 'r-.')
plot([0, Tmax], [vmax1, vmax1], 'g-.')
plot([0, Tmax], [-vmax1, -vmax1], 'g-.')
axis([0, Tmax, -Vmax, Vmax])

figure(2)
clf()
hold('on')
traj0.Plotdd(dt)
traj1.Plotdd(dt, '--')
plot([0, Tmax], [amax0, amax0], 'r-.')
plot([0, Tmax], [-amax0, -amax0], 'r-.')
plot([0, Tmax], [amax1, amax1], 'g-.')
plot([0, Tmax], [-amax1, -amax1], 'g-.')
axis([0, Tmax, -Amax, Amax])

x.WriteProfilesList()
profileslist = TOPPpy.ProfilesFromString(x.resprofilesliststring)
figure(3)
clf()
hold('on')
mvcbobrow = profileslist.pop(0)
plot(mvcbobrow[2], mvcbobrow[3], 'm--', linewidth=4)
mvcdirect = profileslist.pop(0)
plot(mvcdirect[2], mvcdirect[3], 'c--', linewidth=4)
for p in profileslist:
    plot(p[2], p[3])

axis([0, mvcbobrow[0], 0, 2 * max([max(p[3]) for p in profileslist])])

raw_input()
