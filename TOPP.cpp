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



////////////////////////////////////////////////////////////////////
/////////////////////////// Constraints ////////////////////////////
////////////////////////////////////////////////////////////////////

void Constraints::Preprocess(Trajectory& trajectory0, const Tunings& tunings0){
    ptrajectory = &trajectory0;
    tunings = tunings0;
    Discretize();
    ComputeMVC();
    FindSwitchPoints();
}


void Constraints::Discretize(){
    ndiscrsteps = int((ptrajectory->duration+TINY)/tunings.discrtimestep);
    tunings.discrtimestep = ptrajectory->duration/ndiscrsteps;
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
    switchpointslist.push_back(SwitchPoint(discrsvect[iadd],switchpointtype));
}



void Constraints::FindTangentSwitchPoints(){
    if(ndiscrsteps<3) {
        return;
    }
    int i = 1;
    dReal s,sd,snext,sdnext,alpha,beta,diff,diffprev;
    std::pair<dReal,dReal> sddlimits;

    s = discrsvect[i];
    snext = discrsvect[i+1];
    sd = SdLimitMVC(s);
    sdnext = SdLimitMVC(snext);
    sddlimits = SddLimits(s,sd);
    alpha = sddlimits.first;
    beta = sddlimits.second;
    diffprev = alpha/sd - (sdnext-sd)/tunings.discrtimestep;

    for(int i=2; i<ndiscrsteps-1; i++) {
        s = discrsvect[i];
        snext = discrsvect[i+1];
        sd = SdLimitMVC(s);
        sdnext = SdLimitMVC(snext);
        sddlimits = SddLimits(s,sd);
        alpha = sddlimits.first;
        beta = sddlimits.second;
        diff = alpha/sd - (sdnext-sd)/tunings.discrtimestep;
        if(diffprev*diff<0) {
            AddSwitchPoint(i,SPT_TANGENT);
        }
        diffprev = diff;
    }
}

void Constraints::FindDiscontinuousSwitchPoints(){

}




////////////////////////////////////////////////////////////////////
//////////////////////////// Profile ///////////////////////////////
////////////////////////////////////////////////////////////////////

Profile::Profile(std::list<dReal>& slist, std::list<dReal>& sdlist, std::list<dReal>&  sddlist, dReal integrationtimestep0){
    svect = std::vector<dReal>(slist.begin(), slist.end());
    sdvect = std::vector<dReal>(sdlist.begin(), sdlist.end());
    sddvect = std::vector<dReal>(sddlist.begin(), sddlist.end());
    // TODO: handle the case of last step with variable width
    integrationtimestep = integrationtimestep0;
    nsteps = svect.size();
    duration = integrationtimestep * (nsteps-1);
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


//Find t from s, starting search from startindex
bool Profile::Invert(dReal s,  dReal& t, bool searchbackward){
    if(currentindex<0 || currentindex>nsteps-1) {
        return false;
    }
    if(s<svect[0]-TINY || s>svect[nsteps-1]+TINY) {
        return false;
    }
    if(not searchbackward) {
        if(svect[currentindex]>s) {
            return false;
        }
        while(currentindex+1<=nsteps-1 && svect[currentindex+1]<s) {
            currentindex++;
        }
        if(currentindex+1>nsteps-1) {
            return false;
        }
        dReal tres;
        assert(SolveQuadraticEquation(svect[currentindex]-s,sdvect[currentindex],0.5*sddvect[currentindex],0,integrationtimestep,tres));
        t = currentindex*integrationtimestep + tres;
        return true;
    }
    else{
        if(currentindex==0 || svect[currentindex]<s) {
            return false;
        }
        while(currentindex-1>=0 && svect[currentindex-1]>s) {
            currentindex--;
        }
        if(currentindex-1<0) {
            return false;
        }
        dReal tres;
        assert(SolveQuadraticEquation(svect[currentindex-1]-s,sdvect[currentindex-1],0.5*sddvect[currentindex-1],0,integrationtimestep,tres));
        t = (currentindex-1)*integrationtimestep + tres;
        return true;
    }
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





////////////////////////////////////////////////////////////////////
//////////////////////// Integration ///////////////////////////////
////////////////////////////////////////////////////////////////////


int IntegrateForward(Constraints& constraints, dReal sstart, dReal sdstart, Profile& profile, int maxsteps, std::list<Profile>& testprofileslist){
    dReal dt = constraints.tunings.integrationtimestep;
    dReal scur = sstart, sdcur = sdstart;
    std::list<dReal> slist, sdlist, sddlist;
    bool cont = true;
    int returntype;

    // Initialize the currentindex of the profiles for search purpose
    if(testprofileslist.size()>0) {
        std::list<Profile>::iterator it = testprofileslist.begin();
        while(it != testprofileslist.end()) {
            it->currentindex = 0;
            it++;
        }
    }

    // Integrate forward
    while(cont) {
        if(slist.size()>maxsteps) {
            returntype = IRT_MAXSTEPS;
            break;
        }
        else if(scur > constraints.ptrajectory->duration) {
            //TODO: change the time step of previous step to reach the end
            returntype = IRT_END;
            break;
        }
        else if(sdcur < 0) {
            //TODO: double check whether alpha>beta
            returntype = IRT_BOTTOM;
            break;
        }
        else if(sdcur > constraints.SdLimitMVC(scur)) {
            //TODO: Zlajpah
            returntype = IRT_MVC;
            break;
        }
        else if(IsAboveProfilesList(scur,sdcur,testprofileslist)) {
            slist.push_back(scur);
            sdlist.push_back(sdcur);
            sddlist.push_back(0);
            returntype = IRT_ABOVEPROFILES;
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



int IntegrateBackward(Constraints& constraints, dReal sstart, dReal sdstart, Profile& profile, int maxsteps, std::list<Profile>& testprofileslist){
    dReal dt = constraints.tunings.integrationtimestep;
    dReal scur = sstart, sdcur = sdstart;
    std::list<dReal> slist, sdlist, sddlist;
    bool cont = true;
    int returntype;
    bool searchbackward = true;

    // Initialize the currentindex of the profiles for search purpose
    if(testprofileslist.size()>0) {
        std::list<Profile>::iterator it = testprofileslist.begin();
        while(it != testprofileslist.end()) {
            it->currentindex = it->nsteps-1;
            it++;
        }
    }

    // Integrate backward
    while(cont) {
        if(slist.size()>maxsteps) {
            returntype = IRT_MAXSTEPS;
            break;
        }
        else if(scur < 0) {
            //TODO: change the time step of previous step to reach the end
            returntype = IRT_END;
            break;
        }
        else if(sdcur < 0) {
            //TODO: double check whether alpha>beta
            returntype = IRT_BOTTOM;
            break;
        }
        else if(IsAboveProfilesList(scur,sdcur,testprofileslist,searchbackward)) {
            slist.push_back(scur);
            sdlist.push_back(sdcur);
            sddlist.push_back(0);
            returntype = IRT_ABOVEPROFILES;
            break;
        }
        else{
            slist.push_back(scur);
            sdlist.push_back(sdcur);
            std ::pair<dReal,dReal> sddlimits = constraints.SddLimits(scur,sdcur);
            dReal alpha = sddlimits.first;
            dReal beta = sddlimits.second;
            sddlist.push_back(alpha);
            dReal sdnext = sdcur - dt * alpha;
            dReal snext = scur - dt * sdcur + 0.5*dt*dt*alpha;
            scur = snext;
            sdcur = sdnext;
        }
    }
    slist.reverse();
    sdlist.reverse();
    sddlist.reverse();
    profile = Profile(slist,sdlist,sddlist,dt);
    return returntype;
}


int ComputeLimitingCurves(){

}



int IntegrateAllProfiles(){

}




////////////////////////////////////////////////////////////////////
///////////////////////// Utilities ////////////////////////////////
////////////////////////////////////////////////////////////////////


bool SolveQuadraticEquation(dReal a0, dReal a1, dReal a2, dReal lowerbound, dReal upperbound, dReal& sol){
    dReal delta = a1*a1- 4*a0*a2;
    if(delta<0) {
        return false;
    }
    if(a2<0) {
        a2=-a2;
        a1=-a1;
        a0=-a0;
    }
    if(a2<=TINY) {
        if(abs(a1)<=TINY) {
            return false; // Must be in a weird case
        }
        sol = -a0/a1;
        return sol>=lowerbound && sol<=upperbound;
    }
    sol = (-a1-sqrt(delta))/(2*a2);
    if(sol>=lowerbound-TINY && sol<=upperbound+TINY) {
        return true;
    }
    sol = (-a1+sqrt(delta))/(2*a2);
    return sol>=lowerbound-TINY && sol<=upperbound+TINY;
}


bool IsAboveProfilesList(dReal s, dReal sd, std::list<Profile>& testprofileslist, bool searchbackward){
    dReal t;
    std::list<Profile>::iterator it = testprofileslist.begin();
    while(it != testprofileslist.end()) {
        if(it->Invert(s,t,searchbackward)) {
            if(sd > it->Eval(t)) {
                return true;
            }
        }
        it++;
    }
    return false;
}



} // end namespace TOPP

