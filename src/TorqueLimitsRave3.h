// -*- coding: utf-8 -*-
// Copyright (C) 2013 Quang-Cuong Pham <cuong.pham@normalesup.org> & Rosen Diankov <rosen.diankov@gmail.com>
//
// This file is part of the Time-Optimal Path Parameterization (TOPP) library.
// TOPP is free software: you can redistribute it and/or modify
// it under the terms of the GNU Lesser General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// at your option) any later version.
//
// This program is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU Lesser General Public License for more details.
//
// You should have received a copy of the GNU Lesser General Public License
// along with this program.  If not, see <http://www.gnu.org/licenses/>.
#ifndef TOPP_TORQUELIMITSRAVE3_H
#define TOPP_TORQUELIMITSRAVE3_H

#ifdef WITH_OPENRAVE

#include <openrave/openrave.h>

#include "TOPP.h"

namespace TOPP {

/// \brief gets all the information from OpenRAVE structures. Use the active DOFs of the robot
///
/// Uses torque limits that are a function of the current speed
/// Bounds torque/inertia limits
class TorqueLimitsRave3 : public QuadraticConstraints
{
public:
    TorqueLimitsRave3(OpenRAVE::RobotBasePtr probot, OpenRAVE::TrajectoryBaseConstPtr ptraj, dReal discrtimestep);

protected:
    //////////////// Specific members and methods //////////////////////
    std::vector<dReal> taumin, taumax; // Torque limits
    OpenRAVE::RobotBasePtr _probot;
};

}

#endif
#endif
