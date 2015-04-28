// -*- coding: utf-8 -*-
// Copyright (C) 2013 Quang-Cuong Pham <cuong.pham@normalesup.org>
//
// This file is part of the Time-Optimal Path Parameterization (TOPP) library.
// TOPP is free software: you can redistribute it and/or modify
// it under the terms of the GNU Lesser General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// at your option, any later version.
//
// This program is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU Lesser General Public License for more details.
//
// You should have received a copy of the GNU Lesser General Public License
// along with this program.  If not, see <http://www.gnu.org/licenses/>.


#ifndef TorqueLimits_H
#define TorqueLimits_H


#include "TOPP.h"


namespace TOPP {

class TorqueLimits : public Constraints {
public:
    TorqueLimits() : Constraints(){
    }
    TorqueLimits(const std::string& constraintsstring);

    //////////////// Overloaded methods //////////////////////
    std::pair<dReal,dReal> SddLimits(dReal s, dReal sd);
    dReal SdLimitBobrowInit(dReal s);
    void FindSingularSwitchPoints();
    void ComputeSlopeDynamicSingularity(dReal s, dReal sd, std::vector<dReal>& slopesvector);

    //////////////// Specific members and methods //////////////////////
    std::vector<dReal> taumin, taumax; // Torque limits
    std::vector<std::vector<dReal> > avect, bvect, cvect; // Dynamics coefficients (torque case, different from general case of QuadraticConstraints)
    void InterpolateDynamics(dReal s, std::vector<dReal>& a, std::vector<dReal>& b, std::vector<dReal>& c); // Linearly interpolate dynamics coefficients
};
}

#endif
