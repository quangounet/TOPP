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

TOPP_UNSPEC = 0
TOPP_OK = 1
TOPP_CANNOT_PREPROCESS = 2
TOPP_SHORT_TRAJ = 3
TOPP_MVC_HIT_ZERO = 4
TOPP_CLC_ERROR = 5
TOPP_SDBEGMIN_TOO_HIGH = 6
TOPP_SDENDMIN_TOO_HIGH = 7
TOPP_FWD_HIT_ZERO = 8
TOPP_BWD_HIT_ZERO = 9
TOPP_FWD_FAIL = 10
TOPP_BWD_FAIL = 11

MESSAGES = {
    TOPP_UNSPEC: "unspecified error",
    TOPP_OK: "everything OK",
    TOPP_CANNOT_PREPROCESS: "cannot preprocess trajectory",
    TOPP_SHORT_TRAJ: "trajectory too short",
    TOPP_MVC_HIT_ZERO: "MVC hit the sd=0 axis",
    TOPP_CLC_ERROR: "some CLC error",
    TOPP_SDBEGMIN_TOO_HIGH: "sdbegmin is too high",
    TOPP_SDENDMIN_TOO_HIGH: "sdendmin is too high",
    TOPP_FWD_HIT_ZERO: "forward integration hit the sd=0 axis",
    TOPP_BWD_HIT_ZERO: "backward integration hit the sd=0 axis",
    TOPP_FWD_FAIL: "forward integration failed",
    TOPP_BWD_FAIL: "backward integration failed"
}


class NoTrajectoryFound(Exception):
    def __init__(self, error_code):
        self.error_code = error_code

    def __str__(self):
        return MESSAGES[self.error_code]
