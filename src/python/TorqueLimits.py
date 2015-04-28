# -*- coding: utf-8 -*-
# Copyright (C) 2013 Stéphane Caron <caron@ynl.t.u-tokyo.ac.jp>
#
# This file is part of the Time-Optimal Path Parameterization (TOPP) library.
# TOPP is free software: you can redistribute it and/or modify
# it under the terms of the GNU Lesser General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# at your option, any later version.
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


class TorqueLimits(RAVEBindings):
    """Bindings for the 'TorqueLimitsRave' problem."""

    def __init__(self, robot, traj, taumin, taumax, vmax, discrtimestep=None,
                 integrationtimestep=None):
        constring = str(discrtimestep) + "\n"
        constring += vect2str(taumin) + "\n"
        constring += vect2str(taumax) + "\n"
        constring += vect2str([0, 0])  # TODO: non-zero vmax
        trajstring = str(traj)
        super(TorqueLimits, self).__init__(
            robot, "TorqueLimitsRave", constring, trajstring)
