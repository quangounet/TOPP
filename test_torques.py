import TOPPbindings
import TOPPpy
import time
import string
from pylab import *
from openravepy import *


# Load robot
env = Environment() # create openrave environment
#env.SetViewer('qtcoin') # attach viewer (optional)
env.Load('robots/twodof.robot.xml')
robot=env.GetRobots()[0]
grav=[0,0,-9.8]
n=robot.GetDOF()
dof_lim=robot.GetDOFLimits()
vel_lim=robot.GetDOFVelocityLimits()
robot.SetDOFLimits(-10*ones(n),10*ones(n))
robot.SetDOFVelocityLimits(100*vel_lim)
robot.SetTransform(array([[0,0,1,0],[0,1,0,0],[-1,0,0,0.3],[0,0,0,1]]))
k=robot.GetLinks()[0]
g=k.GetGeometries()[0]
g.SetAmbientColor([0.0,0,0])
g.SetDiffuseColor([0.0,0,0])
k=robot.GetLinks()[1]
g=k.GetGeometries()[0]
g.SetAmbientColor([0.6,0,0])
g.SetDiffuseColor([0.6,0,0])
k=robot.GetLinks()[2]
g=k.GetGeometries()[0]
g.SetAmbientColor([0,0,0.6])
g.SetDiffuseColor([0,0,0.6])
#V=env.GetViewer()
M=array([[ 0,   0,  -1,   1.1],
         [  1,   0,  0,   0],
         [  0,  -1,  0,   0.3],
         [  0,   0,   0,   1]])
#V.SetCamera(M)


# Parameters
taumin = [-1010000,-10000000]
taumax = [1000000,1000000]
discrtimestep = 0.01;
integrationtimestep = 0.01;
sdprecision = 0.01;
passswitchpointnsteps = 20;
reparamtimestep = 0.01;

tuningsstring = "%f %f %f %d %f"%(discrtimestep,integrationtimestep,sdprecision,passswitchpointnsteps,reparamtimestep);
trajectorystring = "2 \n 2\n 1 1 0 1\n 0 2 0 -1\n 3\n 2\n 11 13 6 0.1666666666666\n -4 -10 -6 0.5";
constraintstring = string.join([str(x) for x in taumin]) + "\n" + string.join([str(a) for a in taumax])


# Sampling trajectory
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
        constraintstring += "\n" + string.join([str(x) for x in to])
        constraintstring += "\n" + string.join([str(x) for x in tm+tc])
        constraintstring += "\n" + string.join([str(x) for x in tc])

duration0 = time.time()-start

x = TOPPbindings.TOPPProblem("TorqueLimits",constraintstring,trajectorystring,tuningsstring);
x.RunPP(1e-4,1e-4)

duration1 = time.time()-start

print "Computation time: ", duration0, " ", duration1-duration0
print "Trajectory duration: ", x.resduration


ion()
dt = 0.1

x.WriteProfilesList()
profileslist = TOPPpy.ProfilesFromString(x.resprofilesliststring)
figure(3)
clf()
hold('on')
mvc = profileslist.pop(0)
plot(mvc[2],mvc[3],'k',linewidth=2)
for p in profileslist:
    plot(p[2],p[3])

axis([0,mvc[0],0,2*max([max(p[3]) for p in profileslist])])

raw_input()
