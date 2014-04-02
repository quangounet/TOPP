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

error_msg = [  # needs to be in sync with Errors.h
    "unspecified error",
    "everything OK",
    "trajectory too short",
    "MVC hit the sd=0 axis",
    "some CLC error",
    "sdbegmin is too high",
    "sdendmin is too high",
    "forward integration hit the sd=0 axis",
    "backward integration hit the sd=0 axis",
    "forward integration failed",
    "backward integration failed",
]


class NoTrajectoryFound(Exception):
    def __init__(self, error_code=None, wrapper=None):
        """Report when no trajectory is found.

        error_code -- return code from ComputeProfiles
        wrapper -- TOPP python wrapper (e.g. QuadraticConstraints object)

        """
        self.error_code = error_code
        self.wrapper = wrapper

    def __str__(self):
        return error_msg[self.error_code]
