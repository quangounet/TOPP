from TOPPpy import RAVEBindings
from Utilities import vect2str

class Bindings(RAVEBindings):
    """Bindings for the 'TorqueLimitsRave' problem."""

    def __init__(self, robot, traj, activedofs, activelinks, taumin, taumax,
                 zmplimits, vmax, qdefault, support_foot, discrtimestep=None,
                 integrationtimestep=None):
        constring = "%f" % discrtimestep
        constring += "\n" + vect2str(activedofs)
        constring += "\n" + vect2str(activelinks)
        constring += "\n" + vect2str(vmax)
        constring += "\n" + vect2str(taumin)
        constring += "\n" + vect2str(taumax)
        constring += "\n" + vect2str(zmplimits)
        constring += "\n" + vect2str(qdefault)
        constring += "\n" + support_foot
        trajstring = str(traj)
        super(Bindings, self).__init__(robot, "ZMPTorqueLimits", constring,
                                       trajstring)
