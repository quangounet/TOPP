import TOPPbindings
import TOPPpy
import time
from pylab import *

constraintstring = "15 10\n  0 0";
tuningsstring = "0.01 0.01 0.01 20 0.01";
trajectorystring = "2 \n 2\n 1 1 0 1\n 0 2 0 -1\n 3\n 2\n 11 13 6 0.1666666666666\n -4 -10 -6 0.5";


traj0 = TOPPpy.PiecewisePolynomialTrajectory(trajectorystring)

start = time.time()

reps = 10
for i in range(reps):
    x = TOPPbindings.TOPPProblem("KinematicLimits",constraintstring,trajectorystring,tuningsstring);
    x.RunPP(1e-4,1e-4)
    #x.RunVIP(1,2)

print "Computation time: ", (time.time()-start)/reps
print "Trajectory duration: ", x.resduration
#print "[sdendmin,sdendmax] = [", x.sdendmin ,",", x.sdendmax, "]"


ion()
dt = 0.1


x.WriteResultTrajectory()
traj1 = TOPPpy.PiecewisePolynomialTrajectory(x.restrajectorystring)

dt = 0.1
tvect = arange(0,traj1.duration+dt,dt)
qdd = array([traj1.Evaldd(t) for t in tvect])
print "Max qdd: ", max(abs(qdd[:,0])) ,"," , max(abs(qdd[:,1])) 

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

figure(2)
clf()
hold('on')
traj0.Plotdd(dt)
traj1.Plotdd(dt,'--')


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


#x.WriteResult()
#print x.restrajectorystring;
