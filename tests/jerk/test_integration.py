import string
import time
import numpy as np
import pylab
import matplotlib.pyplot as plt
from TOPP import TOPPbindings
from TOPP import TOPPpy
from TOPP import Trajectory

pylab.ion()

# Sample trajectory
ndof = 1
trajectorystring = """1.0\n1\n0.0 1.0"""
    
# Constraints
vmax = 5*np.ones(ndof)  # Velocity limits
amax = 10*np.ones(ndof) # Acceleration limits
discrtimestep = 0.005

traj0 = Trajectory.PiecewisePolynomialTrajectory.FromString(trajectorystring)
constraintstring = str(discrtimestep)
constraintstring += "\n" + string.join([str(v) for v in vmax])
constraintstring += "\n" + string.join([str(a) for a in amax])
x = TOPPbindings.TOPPInstance(None, "KinematicLimits", constraintstring, trajectorystring)

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

jmax = 50
plt.axis([0, 1, 0, 3.5])
res0 = Integration.IntegrateFWFollowingDelta(x, traj0, 0, 0.5, jerkmax=jmax)
res1 = Integration.IntegrateFWFollowingBeta(x, res0[0][-1], 0.55, sdstart=res0[1][-1])
res2 = Integration.IntegrateBWFollowingDelta(x, traj0, 1.0, 0.5, jerkmax=jmax)
res3 = Integration.IntegrateBWFollowingAlpha(x, res2[0][0], 0.45, sdstart=res2[1][0])

istart = 117 # manually found
sstart = res1[0][istart]
sdstart = res1[1][istart]
res4 = Integration.IntegrateFWFollowingGamma(x, traj0, sstart, 1.0 - sstart, 
                                             sdstart=sdstart, 
                                             sddstart=x.GetBeta(sstart, sdstart))


plt.plot(res0[0], res0[1], 'y', linewidth=2, label='delta') # delta
plt.plot(res1[0], res1[1], 'g', linewidth=2, label='beta') # beta
plt.plot(res2[0], res2[1], 'y', linewidth=2) # delta
plt.plot(res3[0], res3[1], 'r', linewidth=2, label='alpha') # alpha
plt.plot(res4[0], res4[1], 'c', linewidth=2, label='gamma') # gamma
plt.legend()

################################################################################

from mpl_toolkits.mplot3d import Axes3D
fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')
ax.set_xlabel('s')
ax.set_ylabel('sd')
ax.set_zlabel('sdd')

ax.plot(res0[0], res0[1], res0[2], color='y', linewidth=2)

beta_res1 = [x.GetBeta(s, sd) for (s, sd) in zip(res1[0], res1[1])]
ax.plot(res1[0], res1[1], beta_res1, color='g', linewidth=2)

ax.plot(res2[0], res2[1], res2[2], color='y', linewidth=2)

alpha_res3 = [x.GetAlpha(s, sd) for (s, sd) in zip(res3[0], res3[1])]
ax.plot(res3[0], res3[1], alpha_res3, color='r', linewidth=2)

ax.plot(res4[0], res4[1], res4[2], color='c', linewidth=2)

Integration.PlotVectorFields(x, traj0, ax, jmax, prec=5, length=0.5)
