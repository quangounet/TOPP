from numpy import *
from pylab import *
import time
import TOPP
from TOPP import TOPPpy
from TOPP import TOPPbindings
from TOPP import Trajectory

def Extractabc(abc):
    lista = [float(x) for x in abc[0].split()]
    listb = [float(x) for x in abc[1].split()]
    listc = [float(x) for x in abc[2].split()]
    n= len(lista)/6
    a = zeros((n,6))
    b = zeros((n,6))
    c = zeros((n,6))
    for i in range(n):
        a[i,:] = lista[i*6:i*6+6]
        b[i,:] = listb[i*6:i*6+6]
        c[i,:] = listc[i*6:i*6+6]
    return a, b, c

trajstr = """0.500000
3
0.0 0.0 -35.9066153846 47.930797594
0.0 0.0 -0.645566686001 1.11351913336
0.0 0.0 8.3609538376 -11.3580450529
0.500000
3
0.0 0.1 -0.159247917034 0.0789972227119
0.0 0.1 -32.7720979649 43.5627972865
0.0 0.1 0.958473557774 -1.41129807703"""
traj = Trajectory.PiecewisePolynomialTrajectory.FromString(trajstr)
inertia = eye(3)
vmax = ones(3)
accelmax = ones(3)
discrtimestep= 1e-3
constraintsstr = str(discrtimestep)
constraintsstr += "\n" + " ".join([str(a) for a in accelmax]) 
for v in inertia:
    constraintsstr += "\n" + " ".join([str(i) for i in v])
#When Inertia is an Identity matrix, angular accelerations are the same as torques
t0 = time.time()
abc = TOPPbindings.RunComputeSO3Constraints(trajstr,constraintsstr)
a,b,c = Extractabc(abc)

t1 = time.time()

topp_inst = TOPP.QuadraticConstraints(traj, discrtimestep, vmax, list(a), list(b), list(c))

x = topp_inst.solver
ret = x.RunComputeProfiles(0,0)
    
x.ReparameterizeTrajectory()
t2 = time.time()
print("Compute a,b,c:", t1-t0)
print("Run TOPP:", t2-t1)
print("Total:", t2-t0)

# Display results
ion()
x.WriteProfilesList()
x.WriteSwitchPointsList()
profileslist = TOPPpy.ProfilesFromString(x.resprofilesliststring)
switchpointslist = TOPPpy.SwitchPointsFromString(x.switchpointsliststring)
TOPPpy.PlotProfiles(profileslist,switchpointslist,4)
    
x.WriteResultTrajectory()

traj1 = Trajectory.PiecewisePolynomialTrajectory.FromString(x.restrajectorystring)
dtplot = 0.01
TOPPpy.PlotKinematics(traj,traj1,dtplot,vmax,accelmax)


input()
