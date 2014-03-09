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
#ifdef WITH_OPENRAVE

#include "TorqueLimitsRave.h"

using namespace OpenRAVE;

namespace TOPP {


TorqueLimitsRave::TorqueLimitsRave(RobotBasePtr probot, std::string& constraintsstring, Trajectory* ptraj){
    trajectory = *ptraj;
    int ndof = trajectory.dimension;
    std::istringstream iss(constraintsstring);
    iss >> discrtimestep;
    ReadVectorFromStream(iss, ndof, vmax);
    ReadVectorFromStream(iss, ndof, taumin);
    ReadVectorFromStream(iss, ndof, taumax);
    hasvelocitylimits = false;
    FOREACH(itv, vmax) {
        if( std::abs(*itv) > TINY ) {
            hasvelocitylimits = true;
            break;
        }
    }

    // Define the avect, bvect, cvect
    int ndiscrsteps = int((trajectory.duration+1e-10)/discrtimestep)+1;
    std::vector<dReal> q(ndof), qd(ndof), qdd(ndof), tmp0(ndof), tmp1(ndof), torquesimple;
    boost::array< std::vector< dReal >, 3 > torquecomponents;
    {
        avect.resize(ndiscrsteps);
        bvect.resize(ndiscrsteps);
        cvect.resize(ndiscrsteps);
        EnvironmentMutex::scoped_lock lock(probot->GetEnv()->GetMutex()); // lock environment
        for(int i = 0; i<ndiscrsteps; i++) {
            dReal s = i*discrtimestep;
            trajectory.Eval(s,q);
            trajectory.Evald(s,qd);
            trajectory.Evaldd(s,qdd);
            probot->SetDOFValues(q,KinBody::CLA_Nothing);
            probot->SetDOFVelocities(qd,KinBody::CLA_Nothing);
            probot->ComputeInverseDynamics(torquesimple,qd);
            probot->ComputeInverseDynamics(torquecomponents,qdd);
            VectorAdd(torquesimple,torquecomponents[1], bvect[i], 1, -1);
            VectorAdd(bvect[i], torquecomponents[2], avect[i], 1, -1);
            VectorAdd(torquecomponents[0],torquecomponents[1], bvect[i]);
            cvect[i] = torquecomponents[2];
        }
    }
}


}

#endif
