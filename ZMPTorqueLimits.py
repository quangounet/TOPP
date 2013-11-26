import TOPPpy
import TOPPbindings

from TOPPpy import vect2str


class RaveInstance(TOPPpy.RaveInstance):
    def __init__(self, robot, traj, tunings, activedofs, activelinks, taumin,
                 taumax, zmplimits, vmax, **kwargs):
        super(RaveInstance, self).__init__(robot, traj, taumin, taumax, vmax,
                                           **kwargs)

        buffsize = 200000
        tunstring = "%f %f %f %d" % (self.discrtimestep,
                                     self.integrationtimestep,
                                     self.reparamtimestep,
                                     self.passswitchpointnsteps)

        trajstring = str(traj)

        constraintstring = vect2str(activedofs)
        constraintstring += "\n" + vect2str(activelinks)
        constraintstring += "\n" + vect2str(taumin)
        constraintstring += "\n" + vect2str(taumax)
        constraintstring += "\n" + vect2str(zmplimits)
        constraintstring += "\n" + vect2str(vmax)

        print "tuningsstring = \"\"\"" + tunstring + "\"\"\"\n"
        print "constraintstring = \"\"\"" + constraintstring + "\"\"\"\n"
        print "trajectorystring = \"\"\"" + trajstring + "\"\"\"\n"

        assert len(constraintstring) < buffsize, \
            "%d is bigger than buffer size" % len(constraintstring)
        assert len(trajstring) < buffsize
        assert len(tunstring) < buffsize

        self.solver = TOPPbindings.TOPPInstance(
            "ZMPTorqueLimits", constraintstring, trajstring, tunstring, robot)
