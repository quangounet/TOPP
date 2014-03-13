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
#ifndef TOPP_TORQUELIMITSRAVE_H
#define TOPP_TORQUELIMITSRAVE_H

#ifdef WITH_OPENRAVE

#include <openrave/openrave.h>

#include "TOPP.h"
#include "TorqueLimits.h"

namespace TOPP {

/// \brief convers an openrave trajectory to a TOPP::Trajectory
void ConvertToTOPPTrajectory(OpenRAVE::TrajectoryBaseConstPtr pintraj, const OpenRAVE::ConfigurationSpecification& posspec, TOPP::Trajectory& outtraj);

/// \brief convers a TOPP:trajectory to an OpenRAVE trajectory using posspec.
///
/// posspec.GetDOF() has to equal intraj.dimension
void ConvertToOpenRAVETrajectory(const TOPP::Trajectory& intraj, OpenRAVE::TrajectoryBasePtr pouttraj, const OpenRAVE::ConfigurationSpecification& posspec);

/// \brief simple version of torque limits
class TorqueLimitsRave : public TorqueLimits {
public:
    TorqueLimitsRave(OpenRAVE::RobotBasePtr probot, std::string& constraintsstring, Trajectory* ptraj);

};

/// \brief gets all the information from OpenRAVE structures. Use the active DOFs of the robot
///
/// Uses torque limits that are a function of the current speed
/// Bounds torque/inertia limits
class TorqueLimitsRave2 : public Constraints
{
public:
    TorqueLimitsRave2(OpenRAVE::RobotBasePtr probot, OpenRAVE::TrajectoryBaseConstPtr ptraj, dReal discrtimestep);

    //////////////// Overloaded methods //////////////////////
    virtual std::pair<dReal,dReal> SddLimits(dReal s, dReal sd);
    virtual dReal SdLimitBobrowInit(dReal s);
    virtual void FindSingularSwitchPoints();
    virtual void ComputeSlopeDynamicSingularity(dReal s, dReal sd, std::vector<dReal>& slopesvector);

protected:
    //////////////// Specific members and methods //////////////////////
    std::vector<dReal> taumin, taumax; // Torque limits
    std::vector<std::vector<dReal> > avect, bvect, cvect; // Dynamics coefficients: avect = M(q) * dq/ds; bvect = M(q) * d^2q/ds^2 + dq^T*C(q)*dq; cvect = G(q)
    void InterpolateDynamics(dReal s, std::vector<dReal>& a, std::vector<dReal>& b, std::vector<dReal>& c); // Linearly interpolate dynamics coefficients

    OpenRAVE::RobotBasePtr _probot;
};

}

#endif
#endif
