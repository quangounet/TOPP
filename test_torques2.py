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
from openravepy import *


ion()

######################## Useful functions ##############

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



# Load robot
env = Environment() # create openrave environment
env.Load('robots/twodof.robot.xml')
robot=env.GetRobots()[0]
grav=[0,0,-9.8]
n=robot.GetDOF()
dof_lim=robot.GetDOFLimits()
vel_lim=robot.GetDOFVelocityLimits()
robot.SetDOFLimits(-10*ones(n),10*ones(n))
robot.SetDOFVelocityLimits(100*vel_lim)
robot.SetTransform(array([[0,0,1,0],[0,1,0,0],[-1,0,0,0.3],[0,0,0,1]]))


# Parameters
taumin = array([-15,-10])
taumax = array([15,10])
vmax = [2.5,3]

discrtimestep = 0.01;
integrationtimestep = 0.01;
bisectionprecision = 0.01;
passswitchpointnsteps = 5;
reparamtimestep = 0.01;

T=1
[a1,b1,c1,a2,b2,c2] = [3, -3, -3, 0, -2, -2]

#Good example:[-3, 3, 3, -1, 0, -3]

tuningsstring = "%f %f %f %d %f"%(discrtimestep,integrationtimestep,bisectionprecision,passswitchpointnsteps,reparamtimestep);
trajectorystring = "%f\n%d\n%f %f %f\n%f %f %f"%(T,2,c1,b1,a1,c2,b2,a2)
constraintstring = string.join([str(a) for a in vmax])


t0 = time.time()

# Sampling the dynamics of the trajectory in python
traj0 = TOPPpy.PiecewisePolynomialTrajectory(trajectorystring)
ndiscrsteps = int((traj0.duration+1e-10)/discrtimestep)+1;

start = time.time()

for i in range(ndiscrsteps):
    t = i*discrtimestep
    q=traj0.Eval(t)
    qd=traj0.Evald(t)
    qdd=traj0.Evaldd(t)
    with robot:
        robot.SetDOFValues(q)
        robot.SetDOFVelocities(qd)
        tm,tc,tg = robot.ComputeInverseDynamics(qdd,None,returncomponents=True)
        to = robot.ComputeInverseDynamics(qd) - tc - tg
        constraintstring += "\n" + string.join([str(x) for x in to]) + " " + string.join([str(x) for x in -to])
        constraintstring += "\n" + string.join([str(x) for x in tm+tc]) + " " + string.join([str(x) for x in -tm-tc]) 
        constraintstring += "\n" + string.join([str(x) for x in tg-taumax]) + " " + string.join([str(x) for x in -tg+taumin])


t1 = time.time()

# Solve in C++
x = TOPPbindings.TOPPInstance("QuadraticConstraints",constraintstring,trajectorystring,tuningsstring);
ret = x.RunPP(1e-4,1e-4)

t2 = time.time()

print "Sampling time (Python): ", t1 - start
print "Parameterization time (C++): ", t2 - t1
print "Total time: ", t2 -start
print "Duration reparameterized trajectory: ", x.resduration


# if(ret==0):
#     print "Trajectory not time-parameterizable"
#     sys.exit()



# Computations
x.WriteProfilesList()
profileslist = TOPPpy.ProfilesFromString(x.resprofilesliststring)



##### BEGIN TRAJ COMPUTATIONS #####

x.WriteResultTrajectory()
traj1 = TOPPpy.PiecewisePolynomialTrajectory(x.restrajectorystring)
dt = 0.001
tvect0,tauvect0 = ComputeTorques(traj0,robot,dt)
tvect1,tauvect1 = ComputeTorques(traj1,robot,dt)
print "Max torques: ", max(abs(tauvect1[:,0])) ,"," , max(abs(tauvect1[:,1]))


Tmax = max(traj0.duration,traj1.duration)
Vmax = 1.2*max(vmax)
Taumax = 1.2*max(max(taumax),abs(min(taumin)))


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
plot([0,Tmax],[vmax[0],vmax[0]],'r-.')
plot([0,Tmax],[-vmax[0],-vmax[0]],'r-.')
plot([0,Tmax],[vmax[1],vmax[1]],'g-.')
plot([0,Tmax],[-vmax[1],-vmax[1]],'g-.')
axis([0,Tmax,-Vmax,Vmax])

figure(2)
clf()
hold('on')
traj0.Plotdd(dt)
traj1.Plotdd(dt,'--')

figure(3)
clf()
hold('on')
ax=gca()
ax.set_color_cycle(['r','g','b'])
plot(tvect0,tauvect0)
ax.set_color_cycle(['r','g','b'])
plot(tvect1,tauvect1,'--')
plot([0,Tmax],[taumax[0],taumax[0]],'r-.')
plot([0,Tmax],[taumin[0],taumin[0]],'r-.')
plot([0,Tmax],[taumax[1],taumax[1]],'g-.')
plot([0,Tmax],[taumin[1],taumin[1]],'g-.')
axis([0,Tmax,-Taumax,Taumax])


##### END TRAJ COMPUTATIONS #####

figure(4)
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
