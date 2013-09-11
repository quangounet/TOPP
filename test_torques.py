import TOPPbindings
import TOPPpy
import time
import string
from pylab import *
from numpy import *
from openravepy import *


# Old version
import MintimeTrajectory
import MintimeProblemGeneric
import MintimeProblemTorque
import MintimeProfileIntegrator
import Kinodynamic

ion()




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
taumin = [-13,-5]
taumax = [13,5]
discrtimestep = 0.01;
integrationtimestep = 0.01;
sdprecision = 0.01;
passswitchpointnsteps = 20;
reparamtimestep = 0.01;

T=1
[a1,b1,c1,a2,b2,c2]=[2, -1, 0, 1, 3, -3]

tuningsstring = "%f %f %f %d %f"%(discrtimestep,integrationtimestep,sdprecision,passswitchpointnsteps,reparamtimestep);
#trajectorystring = "2 \n 2\n 1 1 0 1\n 0 2 0 -1\n 3\n 2\n 11 13 6 0.1666666666666\n -4 -10 -6 0.5";
trajectorystring = "%f\n%d\n%f %f %f\n%f %f %f"%(T,2,c1,b1,a1,c2,b2,a2)
constraintstring = string.join([str(x) for x in taumin]) + "\n" + string.join([str(a) for a in taumax])


t0 = time.time()

# Old version
tunings=Kinodynamic.Tunings()
tunings.grav=grav
tunings.t_step=discrtimestep
tunings.dt_integ=integrationtimestep # time step to integrate the limiting curves
tunings.disc_thr=1e15
tunings.width=10
tunings.palier=10
tunings.tolerance_ends=1e-2
tunings.smart_mode=False
tunings.zi_palier=10

tunings.dicho_steps=10
tunings.n_close=5
tunings.n_rand=20

pwp_traj=MintimeTrajectory.PieceWisePolyTrajectory([[poly1d([a1,b1,c1]),poly1d([a2,b2,c2])]],[T])
traj=pwp_traj.GetSampleTraj(T,tunings.t_step)
pb=MintimeProblemTorque.MintimeProblemTorque(robot,traj)
pb.set_dynamics_limits([taumin,taumax])
pb.disc_thr=tunings.disc_thr
pb.preprocess()

algo=MintimeProfileIntegrator.MintimeProfileIntegrator(pb)
algo.dt_integ=tunings.dt_integ
algo.width=tunings.width
algo.palier=tunings.palier
algo.tolerance_ends=tunings.tolerance_ends
algo.smart_mode=tunings.smart_mode
algo.zi_palier=tunings.zi_palier
status=algo.compute_limiting_curves()


algo.sdot_init=1e-4 # initial value of sdot
algo.sdot_final=1e-4 # final value of sdot
algo.smooth_threshold=5
algo.smooth_window=100
algo.integrate_all_profiles()
algo.integrate_final()
s_res=algo.s_res
sdot_res=algo.sdot_res
undersample_coef=int(round(pb.t_step/algo.dt_integ))
s_res_u=s_res[range(1,len(s_res),undersample_coef)]
sdot_res_u=sdot_res[range(1,len(s_res),undersample_coef)]
traj2=pwp_traj.ResampleTraj(s_res_u,sdot_res_u,pb.t_step)


figure(5)
clf()
algo.plot_profiles()
axis([0,1,0,5])
grid('on')

t1 = time.time()


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
        constraintstring += "\n" + string.join([str(x) for x in tg])


x = TOPPbindings.TOPPProblem("TorqueLimits",constraintstring,trajectorystring,tuningsstring);
x.RunPP(1e-4,1e-4)


t2 = time.time()

#print "Computation time: ", duration0, " ", duration1-duration0

print "Trajectory duration old : ", traj2.duration , " computed in : ", t1-t0
print "Trajectory duration new: ", x.resduration , " computed in : ", t2-t1


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
