import TOPP
import time

constraintstring = "1.5 4\n  0 0";
tuningsstring = "0.01 0.05 0.01 10 0.01";
trajectorystring = "2 \n 2\n 1 1 0 1\n 0 2 0 -1\n 3\n 2\n 11 13 6 0.1666666666666\n -4 -10 -6 0.5";


start = time.time()

for i in range(10):
    x = TOPP.TOPPProblem(constraintstring,trajectorystring,tuningsstring);
    x.Solve()
    x.WriteResultTrajectory()

print "Computation time: ", time.time()-start
print "Trajectory duration: ", x.resduration


#x.WriteResult()
#print x.restrajectorystring;
