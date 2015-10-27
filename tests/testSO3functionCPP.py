from numpy import *

import TOPP
from TOPP import TOPPpy
from TOPP import TOPPbindings
from TOPP import Trajectory
from pylab import *

def Extractabc(abc): #This function is currently in lie.py
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

with open ("data/lie.traj", "r") as myfile0:
    trajstr=myfile0.read()

with open ("data/lie.constraints", "r") as myfile1:
    constraintsstr=myfile1.read()

vmax = ones(3)
discrtimestep= 1e-2

abc = TOPPbindings.RunComputeSO3Constraints(trajstr,constraintsstr)
a,b,c = Extractabc(abc)

traj = Trajectory.PiecewisePolynomialTrajectory.FromString(trajstr)
topp_inst = TOPP.QuadraticConstraints(traj, discrtimestep, vmax, list(a), list(b), list(c))

x = topp_inst.solver

ret = x.RunComputeProfiles(0,0)
if ret == 1:
    x.ReparameterizeTrajectory()
    x.WriteResultTrajectory()

traj1 = Trajectory.PiecewisePolynomialTrajectory.FromString(x.restrajectorystring)


