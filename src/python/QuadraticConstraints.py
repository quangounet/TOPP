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


import TOPPpy

from Errors import NoTrajectoryFound, TOPP_OK
from Trajectory import PiecewisePolynomialTrajectory
from TOPPbindings import TOPPInstance
from TOPPpy import ProfilesFromString, SwitchPointsFromString, PlotProfiles


class QuadraticConstraints(object):
    def __init__(self, traj, discrtimestep, vmax, a, b, c):
        constraintstring = str(discrtimestep)
        constraintstring += "\n" + ' '.join(map(str, vmax))
        for i, _ in enumerate(a):
            constraintstring += "\n" + ' '.join(map(str, a[i]))
            constraintstring += "\n" + ' '.join(map(str, b[i]))
            constraintstring += "\n" + ' '.join(map(str, c[i]))
        self.discrtimestep = discrtimestep
        self.solver = TOPPInstance(
            None, "QuadraticConstraints", constraintstring, str(traj))

    def AVP(self, sdmin, sdmax):
        return_code = self.solver.RunAVP(sdmin, sdmax)
        if return_code != TOPP_OK:
            raise NoTrajectoryFound(return_code)
        sdendmin = self.solver.sdendmin
        sdendmax = self.solver.sdendmax
        return (sdendmin, sdendmax)

    def Reparameterize(self, sdbeg=0., sdend=0.):
        return_code = self.solver.RunComputeProfiles(sdbeg, sdend)
        if return_code != TOPP_OK:
            raise NoTrajectoryFound(return_code)

        return_code = self.solver.ReparameterizeTrajectory()
        if return_code < 0:
            raise NoTrajectoryFound(return_code)

        self.solver.WriteResultTrajectory()
        traj_str = self.solver.restrajectorystring
        return PiecewisePolynomialTrajectory.FromString(traj_str)

    def PlotProfiles(self):
        self.solver.WriteProfilesList()
        self.solver.WriteSwitchPointsList()
        profileslist = ProfilesFromString(self.solver.resprofilesliststring)
        switchpointslist = SwitchPointsFromString(
            self.solver.switchpointsliststring)
        PlotProfiles(profileslist, switchpointslist)

    def PlotAlphaBeta(self):
        return TOPPpy.PlotAlphaBeta(self.solver)
