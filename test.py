import TOPP
import TOPPpy
import time
from pylab import *

constraintstring = "4 4\n  0 0";
tuningsstring = "0.01 0.01 0.01 10 0.01";
trajectorystring = "2 \n 2\n 1 1 0 1\n 0 2 0 -1\n 3\n 2\n 11 13 6 0.1666666666666\n -4 -10 -6 0.5";


traj0 = TOPPpy.PieceWisePolyTrajectory(trajectorystring)

start = time.time()

reps = 10
for i in range(reps):
    x = TOPP.TOPPProblem(constraintstring,trajectorystring,tuningsstring);
    x.Solve()

print "Computation time: ", (time.time()-start)/reps
print "Trajectory duration: ", x.resduration


x.WriteResultTrajectory()
traj1 = TOPPpy.PieceWisePolyTrajectory(x.restrajectorystring)
x.WriteProfilesList()
profileslist = TOPPpy.ProfilesFromString(x.resprofilesliststring)

dt = 0.1
tvect = arange(0,traj1.duration+dt,dt)
qdd = array([traj1.Evaldd(t) for t in tvect])
print "Max qdd: ", max(abs(qdd[:,0])) ,"," , max(abs(qdd[:,1])) 




ion()

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

figure(3)
clf()
hold('on')
for p in profileslist:
    plot(p[2],p[3])



show()

raw_input()


#x.WriteResult()
#print x.restrajectorystring;
