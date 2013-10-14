import TOPPbindings
import TOPPpy
import time
import string
import sys
from pylab import *
from numpy import *
from openravepy import *
import unittest



class TOPPTestCase(unittest.TestCase):

    def setUp(self):
        self.durationtolerance = 1e-2
        


class Torques2DTestCase(TOPPTestCase):
    def runTest(self):
        # Tunings
        discrtimestep = 0.01;
        integrationtimestep = 0.01;
        bisectionprecision = 0.01;
        passswitchpointnsteps = 5;
        reparamtimestep = 0.01;
        tuningsstring = "%f %f %f %d %f"%(discrtimestep,integrationtimestep,bisectionprecision,passswitchpointnsteps,reparamtimestep);

        # Trajectory
        T=1
        [a1,b1,c1,a2,b2,c2] = [-3, 3, 3, -1, 0, -3]
        trajectorystring = "%f\n%d\n%f %f %f\n%f %f %f"%(T,2,c1,b1,a1,c2,b2,a2)
        traj0 = TOPPpy.PiecewisePolynomialTrajectory(trajectorystring)
        ndiscrsteps = int((traj0.duration+1e-10)/discrtimestep)+1;

        # Constraints
        taumin = array([-15,-10])
        taumax = array([15,10])
        vmax = [2.5,3]
        constraintstring = string.join([str(a) for a in vmax])
        # Load robot
        env = Environment() # create openrave environment
        env.Load('robots/twodof.robot.xml')
        robot=env.GetRobots()[0]
        grav=[0,0,-9.8]
        n=robot.GetDOF()
        dof_lim=robot.GetDOFLimits()
        vel_lim=robot.GetDOFVelocityLimits()
        robot.SetDOFLimits(-10*ones(n),10*ones(n))
        robot.SetDOFVelocityLimits(100*vel_lim)
        robot.SetTransform(array([[0,0,1,0],[0,1,0,0],[-1,0,0,0.3],[0,0,0,1]]))
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
                constraintstring += "\n" + string.join([str(x) for x in to]) + " " + string.join([str(x) for x in -to])
                constraintstring += "\n" + string.join([str(x) for x in tm+tc]) + " " + string.join([str(x) for x in -tm-tc]) 
                constraintstring += "\n" + string.join([str(x) for x in tg-taumax]) + " " + string.join([str(x) for x in -tg+taumin])

        # Run TOPP
        x = TOPPbindings.TOPPInstance("QuadraticConstraints",constraintstring,trajectorystring,tuningsstring);
        ret = x.RunPP(1e-4,1e-4)

        referenceduration = 0.73795011802
        self.assertTrue(abs(x.resduration-referenceduration) <= self.durationtolerance)


class Torques2DZlajpahTestCase(TOPPTestCase):
    def runTest(self):
        # Tunings
        discrtimestep = 0.01;
        integrationtimestep = 0.01;
        bisectionprecision = 0.01;
        passswitchpointnsteps = 5;
        reparamtimestep = 0.01;
        tuningsstring = "%f %f %f %d %f"%(discrtimestep,integrationtimestep,bisectionprecision,passswitchpointnsteps,reparamtimestep);

        # Trajectory
        T=1
        [a1,b1,c1,a2,b2,c2] =  [3, -3, -3, 0, -2, -2]
        trajectorystring = "%f\n%d\n%f %f %f\n%f %f %f"%(T,2,c1,b1,a1,c2,b2,a2)
        traj0 = TOPPpy.PiecewisePolynomialTrajectory(trajectorystring)
        ndiscrsteps = int((traj0.duration+1e-10)/discrtimestep)+1;

        # Constraints
        taumin = array([-15,-10])
        taumax = array([15,10])
        vmax = [2.5,3]
        constraintstring = string.join([str(a) for a in vmax])
        # Load robot
        env = Environment() # create openrave environment
        env.Load('robots/twodof.robot.xml')
        robot=env.GetRobots()[0]
        grav=[0,0,-9.8]
        n=robot.GetDOF()
        dof_lim=robot.GetDOFLimits()
        vel_lim=robot.GetDOFVelocityLimits()
        robot.SetDOFLimits(-10*ones(n),10*ones(n))
        robot.SetDOFVelocityLimits(100*vel_lim)
        robot.SetTransform(array([[0,0,1,0],[0,1,0,0],[-1,0,0,0.3],[0,0,0,1]]))
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
                constraintstring += "\n" + string.join([str(x) for x in to]) + " " + string.join([str(x) for x in -to])
                constraintstring += "\n" + string.join([str(x) for x in tm+tc]) + " " + string.join([str(x) for x in -tm-tc]) 
                constraintstring += "\n" + string.join([str(x) for x in tg-taumax]) + " " + string.join([str(x) for x in -tg+taumin])

        # Run TOPP
        x = TOPPbindings.TOPPInstance("QuadraticConstraints",constraintstring,trajectorystring,tuningsstring);
        ret = x.RunPP(1e-4,1e-4)

        referenceduration = 0.879688710063
        self.assertTrue(abs(x.resduration-referenceduration) <= self.durationtolerance)
