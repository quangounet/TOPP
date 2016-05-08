import string
import time
import numpy as np
import matplotlib.pyplot as plt
import pylab
import random
from TOPP import TOPPbindings
from TOPP import TOPPpy
from TOPP import Trajectory
from TOPP import Utilities

pylab.ion()
rng = random.SystemRandom()

generate_new_traj = True

def TrajString5thDegree(q0, q1, qs0, qs1, qss0, qss1, T):
    trajectorystring = ''
    ndof = len(q0)
    trajectorystring += '{0:f}\n{1:d}'.format(T, ndof)
    for k in range(ndof):
        a, b, c, d, e, f = Utilities.Interpolate5thDegree(q0[k], q1[k], qs0[k], 
                                                          qs1[k], qss0[k], qss1[k], T)
        trajectorystring += '\n{0:f} {1:f} {2:f} {3:f} {4:f} {5:f}'.\
        format(f, e, d, c, b, a)
    return trajectorystring

# Trajectory
ndof = 2
if generate_new_traj:
    q0 = [random.uniform(-1, 1) for i in xrange(ndof)]
    qd0 = [random.uniform(-1, 1) for i in xrange(ndof)]
    qdd0 = [random.uniform(-1, 1) for i in xrange(ndof)]

    q1 = [random.uniform(-1, 1) for i in xrange(ndof)]
    qd1 = [random.uniform(-1, 1) for i in xrange(ndof)]
    qdd1 = [random.uniform(-1, 1) for i in xrange(ndof)]

    # q2 = [random.uniform(-1, 1) for i in xrange(ndof)]
    # qd2 = [random.uniform(-1, 1) for i in xrange(ndof)]
    # qdd2 = [random.uniform(-1, 1) for i in xrange(ndof)]

    # trajectorystring0 = TrajString5thDegree(q0, q1, qd0, qd1, qdd0, qdd1, 1.0)
    # trajectorystring1 = TrajString5thDegree(q1, q2, qd1, qd2, qdd1, qdd2, 1.0)
    # trajectorystring = trajectorystring0 + '\n' + trajectorystring1
    trajectorystring = TrajString5thDegree(q0, q1, qd0, qd1, qdd0, qdd1, 1.0)
    
else:
    trajectorystring = """1.000000
2
0.870017 0.717759 -0.046980 -16.327319 24.723340 -9.925695
0.836242 0.042483 0.192225 -11.994327 18.208341 -7.333639
"""
#     trajectorystring = """1.000000
# 2
# 0.891655 0.466621 -0.332083 -7.440755 10.133411 -3.791841
# -0.261058 -0.570469 0.412685 6.795503 -11.005398 4.478107
# 1.000000
# 2
# -0.072991 -0.945367 0.227714 14.873039 -21.686929 8.533979
# -0.150630 -0.989646 -0.452124 20.878187 -30.454072 12.060103"""
        
# Constraints
vmax = 10*np.ones(ndof)  # Velocity limits
amax = 10*np.ones(ndof) # Acceleration limits
jmax = 1000
discrtimestep = 0.005

traj0 = Trajectory.PiecewisePolynomialTrajectory.FromString(trajectorystring)
constraintstring = str(discrtimestep)
constraintstring += "\n" + string.join([str(v) for v in vmax])
constraintstring += "\n" + string.join([str(a) for a in amax])
x = TOPPbindings.TOPPInstance(None, "KinematicLimits", constraintstring, trajectorystring)
x.integrationtimestep = 1e-4

sdstart = 0
sdend = 0
ret = x.RunComputeProfiles(sdstart, sdend)
x.ReparameterizeTrajectory()
x.WriteProfilesList()
x.WriteSwitchPointsList()
x.WriteResultTrajectory()
traj1 = Trajectory.PiecewisePolynomialTrajectory.FromString(x.restrajectorystring)

fignum = 1
profileslist = TOPPpy.ProfilesFromString(x.resprofilesliststring)
switchpointslist = TOPPpy.SwitchPointsFromString(x.switchpointsliststring)
TOPPpy.PlotProfiles(profileslist, switchpointslist, fignum)

################################################################################

import Integration
