import TOPPpy
import TOPPbindings

from TOPPpy import vect2str


class RaveInstance(TOPPpy.RaveInstance):
    def __init__(self, robot, traj, activedofs, activelinks, taumin,
                 taumax, zmplimits, vmax, qdefault, support_foot, **kwargs):
        super(RaveInstance, self).__init__(robot, traj, taumin, taumax, vmax,
                                           **kwargs)
        constraintstring = "%f" % self.discrtimestep
        constraintstring += "\n" + vect2str(activedofs)
        constraintstring += "\n" + vect2str(activelinks)
        constraintstring += "\n" + vect2str(vmax)
        constraintstring += "\n" + vect2str(taumin)
        constraintstring += "\n" + vect2str(taumax)
        constraintstring += "\n" + vect2str(zmplimits)
        constraintstring += "\n" + vect2str(qdefault)
        constraintstring += "\n" + support_foot

        print "constraintstring = \"\"\"" + constraintstring + "\"\"\"\n"
        print "trajectorystring = \"\"\"" + str(traj) + "\"\"\"\n"

        self.solver = TOPPbindings.TOPPInstance(
            robot, "ZMPTorqueLimits", constraintstring, str(traj))


def AVP(robot, traj, sdbegmin, sdbegmax, activedofs, activelinks, taumin,
        taumax, zmplimits, vmax, qdefault, support_foot, **kwargs):
    rave_instance = RaveInstance(
        robot, traj, activedofs, activelinks, taumin, taumax, zmplimits, vmax,
        qdefault, support_foot, **kwargs)
    return rave_instance.GetAVP(sdbegmin, sdbegmax)


def Reparameterize(robot, traj, sdbegmin, sdbegmax, activedofs, activelinks,
                   taumin, taumax, zmplimits, vmax, qdefault, support_foot,
                   **kwargs):
    rave_instance = RaveInstance(
        robot, traj, activedofs, activelinks, taumin, taumax, zmplimits, vmax,
        qdefault, support_foot, **kwargs)
    return rave_instance.GetTrajectory(sdbegmin, sdbegmax)
