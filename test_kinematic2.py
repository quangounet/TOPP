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
import string
import time
import sys
from pylab import *


# Constraints :
# amax0 , amax1 \n vmax0, vmax1 (vmax = 0 means no velocity constraints)
amax0 = 15
amax1 = 10
amax = array([amax0,amax1])
vmax0 = 20
vmax1 = 10
constraintstring = "%f %f"%(vmax0,vmax1);


# Tunings :
# - time step for discretizing the MVC
# - time step for integrating the profiles
# - precision for sdot search around switch points
# - number of time steps to integrate around dynamic singularities
# - time step for reparameterization
discrtimestep = 0.01;
integrationtimestep = 0.01;
bisectionprecision = 0.1;
passswitchpointnsteps = 20;
reparamtimestep = 0.01;
tuningsstring = "%f %f %f %d %f"%(discrtimestep,integrationtimestep,bisectionprecision,passswitchpointnsteps,reparamtimestep);



# Trajectory :
# Chunk 1 time duration
# Chunk 1 number of dofs
# Coefficients of the polynomial of dof 1 (from lower to higher)
# Coefficients of the polynomial of dof 2 (from lower to higher) ...
# Chunk 2 ...
trajectorystring = "2 \n 2\n 1 1 0 1\n 0 2 0 -1\n 3\n 2\n 11 13 6 0.1666666666666\n -4 -10 -6 0.5";


# Sampling the dynamics of the trajectory in python
traj0 = TOPPpy.PiecewisePolynomialTrajectory.FromString(trajectorystring)
ndiscrsteps = int((traj0.duration+1e-10)/discrtimestep)+1;

start = time.time()

for i in range(ndiscrsteps):
    t = i*discrtimestep
    q=traj0.Eval(t)
    qd=traj0.Evald(t)
    qdd=traj0.Evaldd(t)
    constraintstring += "\n" + string.join([str(x) for x in qd]) + " " + string.join([str(x) for x in -qd])
    constraintstring += "\n" + string.join([str(x) for x in qdd]) + " " + string.join([str(x) for x in -qdd])
    constraintstring += "\n" + string.join([str(x) for x in -amax]) + " " + string.join([str(x) for x in -amax])


t1 = time.time()



# Run TOPP

x = TOPPbindings.TOPPInstance("QuadraticConstraints",constraintstring,trajectorystring,tuningsstring);
ret = x.RunPP(1e-4,1e-4)

t2 = time.time()

if(ret==0):
    sys.exit()

print "Sampling time (Python): ", t1 - start
print "Parameterization time (C++): ", t2 - t1
print "Total time: ", t2 -start
print "Duration reparameterized trajectory: ", x.resduration

# Run Velocity Interval Propagation
#x.RunVIP(1,2)
#print "[sdendmin,sdendmax] = [", x.sdendmin ,",", x.sdendmax, "]"


# Display results
x.WriteResultTrajectory()
traj1 = TOPPpy.PiecewisePolynomialTrajectory.FromString(x.restrajectorystring)


# Verification
ion()
dt = 0.1
tvect = arange(0,traj1.duration+dt,dt)
qdd = array([traj1.Evaldd(t) for t in tvect])
print "Max acceleration: ", max(abs(qdd[:,0])) ,"," , max(abs(qdd[:,1]))
Tmax = max(traj0.duration,traj1.duration)
Vmax = 1.2*max(vmax0,vmax1)
Amax = 1.2*max(amax0,amax1)

# Plotting
figure(0)
clf()
hold('on')
traj0.Plot(dt)
traj1.Plot(dt,'--')

figure(1)
clf()
hold('on')
traj0.Plotd(dt)
traj1.Plotd(dt,'--')
plot([0,Tmax],[vmax0,vmax0],'r-.')
plot([0,Tmax],[-vmax0,-vmax0],'r-.')
plot([0,Tmax],[vmax1,vmax1],'g-.')
plot([0,Tmax],[-vmax1,-vmax1],'g-.')
axis([0,Tmax,-Vmax,Vmax])

figure(2)
clf()
hold('on')
traj0.Plotdd(dt)
traj1.Plotdd(dt,'--')
plot([0,Tmax],[amax0,amax0],'r-.')
plot([0,Tmax],[-amax0,-amax0],'r-.')
plot([0,Tmax],[amax1,amax1],'g-.')
plot([0,Tmax],[-amax1,-amax1],'g-.')
axis([0,Tmax,-Amax,Amax])

x.WriteProfilesList()
profileslist = TOPPpy.ProfilesFromString(x.resprofilesliststring)
figure(3)
clf()
hold('on')
mvcbobrow = profileslist.pop(0)
plot(mvcbobrow[2],mvcbobrow[3],'m--',linewidth=4)
mvcdirect = profileslist.pop(0)
plot(mvcdirect[2],mvcdirect[3],'c--',linewidth=4)
for p in profileslist:
    plot(p[2],p[3])

axis([0,mvcbobrow[0],0,2*max([max(p[3]) for p in profileslist])])



raw_input()

