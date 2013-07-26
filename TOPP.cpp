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

void Trajectory::Eval(dReal s, std::vector<dReal>& q){

}

void Trajectory::Evald(dReal s, std::vector<dReal>& qd){

}

void Trajectory::Evaldd(dReal s, std::vector<dReal>& qdd){

}





// Constraints


void Constraints::Preprocess(Trajectory& trajectory0, Tunings& tunings0){
    trajectory = trajectory0;
    tunings = tunings0;
    Discretize();
    ComputeMVC();
    FindSwitchPoints();
}


void Constraints::Discretize(){
    ndiscrsteps = int((trajectory.duration+TINY)/tunings.discrtimestep);
    tunings.discrtimestep = trajectory.duration/ndiscrsteps;
    ndiscrsteps++;
    discrsvect.resize(0);
    for(int i=0; i<ndiscrsteps; i++) {
        discrsvect.push_back(i*tunings.discrtimestep);
    }
}


void Constraints::ComputeMVC(){
    for(int i=0; i<ndiscrsteps; i++) {
        mvc.push_back(SdLimitMVC(discrsvect[i]));
    }
}


dReal Constraints::SdLimitDirect(dReal s){
    return 0;
}



void Constraints::FindSwitchPoints(){
    FindTangentSwitchPoints();
    FindSingularSwitchPoints();
    FindDiscontinuousSwitchPoints();
}

void Constraints::AddSwitchPoint(int i, int switchpointtype){
    int iadd = i+1;
    if(mvc[i-1]<std::min(mvc[i],mvc[i+1])) {
        iadd = i-1;
    }
    else{
        if(mvc[i]<mvc[i+1]) {
            iadd = i;
        }
    }
}



void Constraints::FindTangentSwitchPoints(){
    if(ndiscrsteps<3) {
        return;
    }
    int i = 1;
    dReal s,sd,snext,sdnext,alpha,beta,diff,diffprev;
    std::pair<dReal,dReal> sddlimits;

    s = mvc[i];
    snext = mvc[i+1];
    sd = SdLimitMVC(s);
    sdnext = SdLimitMVC(snext);
    sddlimits = SddLimits(s,sd);
    alpha = sddlimits.first;
    beta = sddlimits.second;
    diffprev = alpha/sd - (sdnext-sd)/tunings.discrtimestep;

    for(int i=2; i<ndiscrsteps-1; i++) {
        s = mvc[i];
        snext = mvc[i+1];
        sd = SdLimitMVC(s);
        sdnext = SdLimitMVC(snext);
        sddlimits = SddLimits(s,sd);
        alpha = sddlimits.first;
        beta = sddlimits.second;
        diff = alpha/sd - (sdnext-sd)/tunings.discrtimestep;
        if(diffprev>0 && diff<0) {
            AddSwitchPoint(i,SP_TANGENT);
        }
        diffprev = diff;
    }
}

void Constraints::FindDiscontinuousSwitchPoints(){

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

bool Profile::FindTimestepIndex(dReal t, int& index, dReal& remainder){
    if (t < 0 || t > duration) {
        return false;
    }
    if(duration-t <= 1e-10) {
        index = nsteps-1;
    }
    else{
        index = (int) floor(t/integrationtimestep);
    }
    remainder = t - index * integrationtimestep;
    return true;
}

dReal Profile::Eval(dReal t){
    int index;
    dReal remainder;
    if(FindTimestepIndex(t, index, remainder)) {
        return svect[index] + remainder*sdvect[index] + 0.5*remainder*remainder*sddvect[index];
    }
    else{
        return INF;
    }
}

dReal Profile::Evald(dReal t){
    int index;
    dReal remainder;
    if(FindTimestepIndex(t, index, remainder)) {
        return sdvect[index] + remainder*sddvect[index];
    }
    else{
        return INF;
    }
}


dReal Profile::Evaldd(dReal t){
    int index;
    dReal remainder;
    if(FindTimestepIndex(t, index, remainder)) {
        return sddvect[index];
    }
    else{
        return INF;
    }
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










