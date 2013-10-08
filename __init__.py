import pylab

import TOPPbindings

from TOPPpy import PiecewisePolynomialTrajectory


class NoTrajectoryFound(Exception):
    pass


def __vect_to_str(v):
    return ' '.join(map(str, v))


class Tunings(object):
    def __init__(self, dt, mvc_dt=None, integ_dt=None, sd_prec=None,
                 switchpoint_steps=20, reparam_dt=None):
        """

        mvc_dt -- time step for discretizing the MVC
        integ_dt -- time step for integrating velocity profiles
        sd_prec -- precision for sdot search around switch points
        switchpoint_steps -- number of time steps to integrate around dynamic
                             singularities
        reparam_dt -- time step for reparametrization

        """
        self.mvc_tstep = mvc_dt if mvc_dt else dt
        self.integ_tstep = integ_dt if integ_dt else dt
        self.reparam_tstep = reparam_dt if reparam_dt else dt
        self.sd_prec = sd_prec if sd_prec else dt
        self.switchpoint_steps = switchpoint_steps

    def __str__(self):
        kron = [self.mvc_tstep, self.integ_tstep, self.sd_prec]
        kron.extend([self.switchpoint_steps, self.reparam_tstep])
        return "%f %f %f %d %f" % kron


class TorqueConstraints(object):
    def __init__(self, tau1, tau2=None):
        self.tau_min = tau1 if tau2 else -tau1
        self.tau_max = tau2 if tau2 else +tau1

    def __str__(self):
        tau_min_str = __vect_to_str(self.tau_min)
        tau_max_str = __vect_to_str(self.tau_max)
        return "%s\n%s" % (tau_min_str, tau_max_str)


class RaveTorqueInstance(object):
    def __init__(self, constraints, openrave_robot, traj, tunings):
        self.constraints = constraints
        self.robot = openrave_robot
        self.tunings = tunings
        self.traj = traj

        self.solver = TOPPbindings.TOPPProblem(
            "TorqueLimits", self.get_dynamics_str(self.robot), str(self.traj),
            str(self.tunings))

    def get_dynamics_str(self):
        poly_traj = PiecewisePolynomialTrajectory(str(self.traj))
        assert abs(poly_traj.duration - self.traj.duration) < 1e-6

        nb_steps = 1 + int((self.duration + 1e-10) / self.tunings.disc_tstep)
        trange = pylab.arange(nb_steps) * self.tunings.disc_tstep
        invdyn_str = ''

        for i, t in enumerate(trange):
            q = poly_traj.Eval(t)
            qd = poly_traj.Evald(t)
            qdd = poly_traj.Evaldd(t)
            with self.robot:
                self.robot.SetDOFValues(q)
                self.robot.SetDOFVelocities(qd)
                args, kwargs = (qdd, None), {'returncomponents': True}
                tm, tc, tg = self.robot.ComputeInverseDynamics(*args, **kwargs)
                to = self.robot.ComputeInverseDynamics(qd) - tc - tg
                invdyn_str += "\n" + __vect_to_str(to)       # a vector
                invdyn_str += "\n" + __vect_to_str(tm + tc)  # b vector
                invdyn_str += "\n" + __vect_to_str(tg)       # c vector

        return invdyn_str

    def parametrize_path(self):
        return_code = self.solver.RunPP(1e-4, 1e-4)
        if return_code == 0:
            raise NoTrajectoryFound

        self.solver.WriteResultTrajectory()
        traj_str = self.solver.restrajectorystring
        return PiecewisePolynomialTrajectory.FromString(traj_str)

    def propagate_velocity_interval(self):
        return_code = self.solver.RunVIP(1e-4, 1e-4)
        if return_code == 0:
            raise NoTrajectoryFound

        sd_end_min = self.solver.sdendmin
        sd_end_max = self.solver.sdendmax
        return (sd_end_min, sd_end_max)
