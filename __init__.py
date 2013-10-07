import pylab

import TOPPpy
import TOPPbindings


class NoTrajectoryFound(Exception):
    pass


def __vect_to_str(v):
    return ' '.join(map(str, v))


class Tunings(object):
    def __init__(self, mvc_tstep, integ_tstep, sd_prec, switchpoint_steps,
                 reparam_tstep):
        """

        mvc_tstep -- time step for discretizing the MVC
        integ_tstep -- time step for integrating velocity profiles
        sd_prec -- precision for sdot search around switch points
        switchpoint_steps -- number of time steps to integrate around dynamic
                             singularities
        reparam_tstep -- time step for reparametrization

        """
        self.mvc_tstep = mvc_tstep
        self.integ_tstep = integ_tstep
        self.reparam_tstep = reparam_tstep
        self.sd_prec = sd_prec
        self.switchpoint_steps = switchpoint_steps

    def __str__(self):
        kron = [self.mvc_tstep, self.integ_tstep, self.sd_prec]
        kron.extend([self.switchpoint_steps, self.reparam_tstep])
        return "%f %f %f %d %f" % kron


class Trajectory(object):
    def __init__(self, duration, degree, coefs):
        self.duration = duration
        self.degree = degree
        self.coefs = coefs

    def __str__(self):
        kron = [self.duration, self.degree]
        kron.extend(self.coefs)
        return "%f\n%d\n%f %f %f\n%f %f %f" % kron


class TorqueConstraints(object):
    def __init__(self, tau_min, tau_max):
        self.tau_min = tau_min
        self.tau_max = tau_max

    def __str__(self):
        tau_min_str = __vect_to_str(self.tau_min)
        tau_max_str = __vect_to_str(self.tau_max)
        return "%s\n%s" % (tau_min_str, tau_max_str)


class RaveTOPPInstance(object):
    def __init__(self, constraints, openrave_robot, traj, tunings):
        self.constraints = constraints
        self.robot = openrave_robot
        self.tunings = tunings
        self.traj = traj

        self.solver = TOPPbindings.TOPPProblem(
            "TorqueLimits", self.get_dynamics_str(self.robot), str(self.traj),
            str(self.tunings))

    def get_dynamics_str(self):
        poly_traj = TOPPpy.PiecewisePolynomialTrajectory(str(self.traj))
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
        return self.solver.restrajectory

    def propagate_velocity_interval(self):
        return_code = self.solver.RunVIP(1e-4, 1e-4)
        if return_code == 0:
            raise NoTrajectoryFound

        sd_end_min = self.solver.sdendmin
        sd_end_max = self.solver.sdendmax
        return (sd_end_min, sd_end_max)
