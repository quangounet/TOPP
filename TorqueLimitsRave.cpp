// -*- coding: utf-8 -*-
// Copyright (C) 2013 Quang-Cuong Pham <cuong.pham@normalesup.org>
//
// This file is part of the Time-Optimal Path Parameterization (TOPP) library.
// TOPP is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// at your option) any later version.
//
// This program is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with this program.  If not, see <http://www.gnu.org/licenses/>.


#include "TorqueLimitsRave.h"

#define CLA_Nothing 0

using namespace OpenRAVE;

namespace TOPP {

TorqueLimitsRave::TorqueLimitsRave(const std::string& constraintsstring, Trajectory* ptraj, const Tunings& tunings, RobotBasePtr probot){
    int buffsize = BUFFSIZE;  // TODO: remove this dirty string interface!
    std::vector<dReal> tmpvect;
    char buff[buffsize];
    std::istringstream iss(constraintsstring);
    iss.getline(buff,buffsize);
    VectorFromString(std::string(buff),taumin);
    iss.getline(buff,buffsize);
    VectorFromString(std::string(buff),taumax);
    iss.getline(buff,buffsize);
    VectorFromString(std::string(buff),vmax);
    hasvelocitylimits = VectorMax(vmax) > TINY;
    maxrep = 1;

    // Define the avect, bvect, cvect
    int ndiscrsteps = int((ptraj->duration+1e-10)/tunings.discrtimestep)+1;
    std::vector<dReal> q(ptraj->dimension), qd(ptraj->dimension), qdd(ptraj->dimension), torquesimple;
    boost::array< std::vector< dReal >, 3 > torquecomponents;
    for(int i = 0; i<ndiscrsteps; i++) {
        dReal s = i*tunings.discrtimestep;
        ptraj->Eval(s,q);
        ptraj->Evald(s,qd);
        ptraj->Evaldd(s,qdd);
        probot->SetDOFValues(q,CLA_Nothing);
        probot->SetDOFVelocities(qd,CLA_Nothing);
        probot->ComputeInverseDynamics(torquecomponents,qdd);
        probot->ComputeInverseDynamics(torquesimple,qd);
        avect.push_back(VectorAdd(VectorAdd(torquesimple,torquecomponents[1],1,-1),torquecomponents[2],1,-1));
        bvect.push_back(VectorAdd(torquecomponents[0],torquecomponents[1]));
        cvect.push_back(torquecomponents[2]);
    }
}


}
