import TOPP
x = TOPP.TOPPProblem("con","traj","tun");
print x.trajectorystring;
x.Solve()
print x.restrajectorystring;
