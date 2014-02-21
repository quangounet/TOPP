// -*- coding: utf-8 -*-
// Copyright (C) 2013 Quang-Cuong Pham <cuong.pham@normalesup.org>
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


#include "TorqueLimitsRave.h"

#define CLA_Nothing 0

using namespace OpenRAVE;

namespace TOPP {

TorqueLimitsRave::TorqueLimitsRave(RobotBasePtr probot, std::string& constraintsstring, Trajectory* ptraj){

    trajectory = *ptraj;

    int buffsize = BUFFSIZE;  // TODO: remove this dirty string interface!
    std::vector<dReal> tmpvect;
    char buff[buffsize];
    std::istringstream iss(constraintsstring);
    iss.getline(buff,buffsize);
    discrtimestep = atof(buff);
    iss.getline(buff,buffsize);
    VectorFromString(std::string(buff),vmax);
    iss.getline(buff,buffsize);
    VectorFromString(std::string(buff),taumin);
    iss.getline(buff,buffsize);
    VectorFromString(std::string(buff),taumax);
    hasvelocitylimits = VectorMax(vmax) > TINY;

    trajectory = *ptraj;

    int ndof = trajectory.dimension;

    // Define the avect, bvect, cvect
    int ndiscrsteps = int((trajectory.duration+1e-10)/discrtimestep)+1;
    std::vector<dReal> q(ndof), qd(ndof), qdd(ndof), tmp0(ndof), tmp1(ndof), torquesimple;
    boost::array< std::vector< dReal >, 3 > torquecomponents;
    {
        EnvironmentMutex::scoped_lock lock(probot->GetEnv()->GetMutex()); // lock environment
        for(int i = 0; i<ndiscrsteps; i++) {
            dReal s = i*discrtimestep;
            trajectory.Eval(s,q);
            trajectory.Evald(s,qd);
            trajectory.Evaldd(s,qdd);
            probot->SetDOFValues(q,CLA_Nothing);
            probot->SetDOFVelocities(qd,CLA_Nothing);
            probot->ComputeInverseDynamics(torquesimple,qd);
            probot->ComputeInverseDynamics(torquecomponents,qdd);
            VectorAdd(torquesimple,torquecomponents[1],tmp0,1,-1);
            VectorAdd(tmp0,torquecomponents[2],tmp1,1,-1);
            avect.push_back(tmp1);
            VectorAdd(torquecomponents[0],torquecomponents[1],tmp0);
            bvect.push_back(tmp0);
            cvect.push_back(torquecomponents[2]);
        }
    }
}


}
