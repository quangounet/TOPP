import pylab

import TOPPbindings

from TOPPpy import Polynomial, Chunk, PiecewisePolynomialTrajectory


class NoTrajectoryFound(Exception):
    pass


def _vect_to_str(v):
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
        return "%f %f %f %d %f" % (
            self.mvc_tstep, self.integ_tstep, self.sd_prec,
            self.switchpoint_steps, self.reparam_tstep)


class TorqueConstraints(object):
    def __init__(self, tau_min, tau_max, v_max):
        self.tau_min = tau_min
        self.tau_max = tau_max
        self.v_max = v_max

    def __str__(self):
        tau_min_str = _vect_to_str(self.tau_min)
        tau_max_str = _vect_to_str(self.tau_max)
        v_max_str = _vect_to_str(self.v_max)
        return "%s\n%s\n%s" % (tau_min_str, tau_max_str, v_max_str)


class RaveTorqueInstance(object):
    def __init__(self, constraints, openrave_robot, traj, tunings):
        assert isinstance(constraints, TorqueConstraints)
        assert isinstance(traj, PiecewisePolynomialTrajectory)

        self.constraints = constraints
        self.robot = openrave_robot
        self.tunings = tunings
        self.traj = traj

        input_str = str(constraints) + self.get_dynamics_str()
        self.solver = TOPPbindings.TOPPInstance(
            "TorqueLimits", input_str, str(self.traj), str(self.tunings))

    def get_dynamics_str(self):
        dt = self.tunings.mvc_tstep
        nb_steps = 1 + int((self.traj.duration + 1e-10) / dt)
        trange = pylab.arange(nb_steps) * dt
        invdyn_str = ''
        for i, t in enumerate(trange):
            q = self.traj.Eval(t)
            qd = self.traj.Evald(t)
            qdd = self.traj.Evaldd(t)
            with self.robot:
                self.robot.SetDOFValues(q)
                self.robot.SetDOFVelocities(qd)
                args, kwargs = (qdd, None), {'returncomponents': True}
                tm, tc, tg = self.robot.ComputeInverseDynamics(*args, **kwargs)
                to = self.robot.ComputeInverseDynamics(qd) - tc - tg
                invdyn_str += '\n' + _vect_to_str(to)       # a vector
                invdyn_str += '\n' + _vect_to_str(tm + tc)  # b vector
                invdyn_str += '\n' + _vect_to_str(tg)       # c vector
        return invdyn_str

    def parametrize_path(self):
        return_code = self.solver.RunPP(1e-4, 1e-4)
        if return_code == 0:
            raise NoTrajectoryFound

        self.solver.WriteResultTrajectory()
        traj_str = self.solver.restrajectorystring
        return PiecewisePolynomialTrajectory.FromString(traj_str)

    def propagate_velocity_interval(self, sd_min, sd_max):
        return_code = self.solver.RunVIP(sd_min, sd_max)
        if return_code == 0:
            raise NoTrajectoryFound

        sd_end_min = self.solver.sdendmin
        sd_end_max = self.solver.sdendmax
        return (sd_end_min, sd_end_max)
