# -*- coding: utf-8 -*-
# Copyright (C) 2013 Stéphane Caron <caron@ynl.t.u-tokyo.ac.jp>
#
# This file is part of the Time-Optimal Path Parameterization (TOPP) library.
# TOPP is free software: you can redistribute it and/or modify
# it under the terms of the GNU Lesser General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU Lesser General Public License for more details.
#
# You should have received a copy of the GNU Lesser General Public License
# along with this program.  If not, see <http://www.gnu.org/licenses/>.


from TOPPopenravepy import RAVEBindings
from Utilities import vect2str


class ZMPTorqueLimits(RAVEBindings):
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
        super(ZMPTorqueLimits, self).__init__(
            robot, "ZMPTorqueLimits", constring, trajstring)
