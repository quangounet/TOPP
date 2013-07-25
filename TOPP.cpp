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

#include "TOPP.h"

namespace TOPP {


// Trajectory



// Constraints


Constraints::Constraints(Trajectory trajectory0, Tunings tunings0){
    trajectory = trajectory0;
    tunings = tunings0;

    //Preprocessing
    SampleDynamics();
    ComputeMVC();
    ComputeSwitchPoints();
}


dReal Constraints::SdLimitMVC(dReal s){
    return 0;
}

dReal Constraints::SdLimitDirect(dReal s){
    return 0;
}

std::pair<dReal,dReal> Constraints::SddLimits(dReal s, dReal sd){
    std::pair<dReal,dReal> res;
    res.first=0;
    res.second=0;
}

void GetSwitchPoints(std::list<SwitchPoint>& switchpointslist){
    int x=0;
}


// Profile

Profile::Profile(std::list<dReal>& slist, std::list<dReal>& sdlist, std::list<dReal>&  sddlist, dReal integrationtimestep0){
    svect = std::vector<dReal>(slist.begin(), slist.end());
    sdvect = std::vector<dReal>(sdlist.begin(), slist.end());
    sddvect = std::vector<dReal>(sddlist.begin(), slist.end());
    // TODO: handle the case of last step with variable width
    integrationtimestep = integrationtimestep0;
    nsteps = svect.size();
    duration = integrationtimestep * nsteps;
}

bool Profile::Evalall(dReal t, dReal& s, dReal& sd, dReal& sdd){
    int istep;
    if (t < 0 || t > duration) {
        return false;
    }
    if(duration-t <= 1e-10) {
        istep = nsteps-1;
    }
    else{
        istep = (int) floor(t/integrationtimestep);
    }
    dReal tremain = t - istep * integrationtimestep;
    s = svect.at(istep) + tremain*sdvect.at(istep) + 0.5*tremain*tremain*sddvect.at(istep);
    sd = sdvect.at(istep) + tremain*sddvect.at(istep);
    sdd = sddvect.at(istep);
    return true;
}



// Integration

int IntegrateForward(Constraints& constraints, dReal sstart, dReal sdstart, Profile& profile, dReal maxsteps, std::list<Profile> testprofileslist){
    dReal dt = constraints.tunings.integrationtimestep;
    dReal scur = sstart, sdcur = sdstart;
    std::list<dReal> slist, sdlist, sddlist;
    bool cont = true;
    int returntype;

    while(cont) {
        if(slist.size()>maxsteps) {
            returntype = RT_MAXSTEPS;
            break;
        }
        else if(scur > constraints.trajectory.duration) {
            //TODO: change the time step of previous step to reach the end
            returntype = RT_END;
            break;
        }
        else if(sdcur < 0) {
            //TODO: double check whether alpha>beta
            returntype = RT_BOTTOM;
            break;
        }
        else if(sdcur > constraints.SdLimitMVC(scur)) {
            //TODO: Zlajpah
            returntype = RT_MVC;
            break;
        }
        else{
            slist.push_back(scur);
            sdlist.push_back(sdcur);
            std ::pair<dReal,dReal> sddlimits = constraints.SddLimits(scur,sdcur);
            dReal alpha = sddlimits.first;
            dReal beta = sddlimits.second;
            sddlist.push_back(beta);
            dReal sdnext = sdcur + dt * beta;
            dReal snext = scur + dt * sdcur + 0.5*dt*dt*beta;
            scur = snext;
            sdcur = sdnext;
        }
    }

    profile = Profile(slist,sdlist,sddlist,dt);
    return returntype;
}


} // end namespace TOPP










