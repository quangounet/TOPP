import TOPPpy
import TOPPbindings

from TOPPpy import vect2str


class RaveInstance(TOPPpy.RaveInstance):
    def __init__(self, robot, traj, activedofs, activelinks, taumin,
                 taumax, zmplimits, vmax, qdefault, **kwargs):
        super(RaveInstance, self).__init__(robot, traj, taumin, taumax, vmax,
                                           **kwargs)

        buffsize = 200000
        tuningsstring = "%f %f %f %d" % (
            self.discrtimestep, self.integrationtimestep,
            self.reparamtimestep, self.passswitchpointnsteps)

        constraintstring = vect2str(activedofs)
        constraintstring += "\n" + vect2str(activelinks)
        constraintstring += "\n" + vect2str(taumin)
        constraintstring += "\n" + vect2str(taumax)
        constraintstring += "\n" + vect2str(zmplimits)
        constraintstring += "\n" + vect2str(vmax)
        constraintstring += "\n" + vect2str(qdefault)

        print "tuningsstring = \"\"\"" + tuningsstring + "\"\"\"\n"
        print "constraintstring = \"\"\"" + constraintstring + "\"\"\"\n"
        print "trajectorystring = \"\"\"" + str(traj) + "\"\"\"\n"

        assert len(constraintstring) < buffsize, \
            "%d is bigger than buffer size" % len(constraintstring)
        assert len(str(traj)) < buffsize
        assert len(tuningsstring) < buffsize

        self.solver = TOPPbindings.TOPPInstance(
            "ZMPTorqueLimits", constraintstring, str(traj), tuningsstring,
            robot)


def AVP(robot, traj, sdbegmin, sdbegmax, activedofs, activelinks, taumin,
        taumax, zmplimits, vmax, qdefault, **kwargs):
    rave_instance = RaveInstance(
        robot, traj, activedofs, activelinks, taumin, taumax, zmplimits, vmax,
        qdefault, **kwargs)
    return rave_instance.GetAVP(sdbegmin, sdbegmax)


def Reparameterize(robot, traj, sdbegmin, sdbegmax, activedofs, activelinks,
                   taumin, taumax, zmplimits, vmax, qdefault, **kwargs):
    rave_instance = RaveInstance(
        robot, traj, activedofs, activelinks, taumin, taumax, zmplimits, vmax,
        qdefault, **kwargs)
    return rave_instance.GetTrajectory(sdbegmin, sdbegmax)
