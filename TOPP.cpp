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

Tunings::Tunings(const std::string& tuningsstring){
    std::istringstream iss(tuningsstring);
    std::string sub;
    iss >> sub;
    discrtimestep = atof(sub.c_str());
    iss >> sub;
    integrationtimestep = atof(sub.c_str());
    iss >> sub;
    threshold = atof(sub.c_str());
    iss >> sub;
    passswitchpointnsteps = atof(sub.c_str());
    iss >> sub;
    reparamtimestep = atof(sub.c_str());
}





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
    std::list<SwitchPoint>::iterator it = switchpointslist.begin();
    dReal s = discrsvect[iadd];
    dReal sd = mvc[iadd];
    if(sd>=INF) {
        return;
    }
    while(it!=switchpointslist.end()) {
        if(s == it->s) {
            return;
        }
        if(s<=it->s) {
            break;
        }
        it++;
    }
    switchpointslist.insert(it,SwitchPoint(s,sd,switchpointtype));
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


Profile::Profile(std::list<Profile>& profileslist, dReal integrationtimestep0){
    integrationtimestep = integrationtimestep0;

    // First, reinitalize the profiles for search
    std::list<Profile>::iterator it = profileslist.begin();
    while(it != profileslist.end()) {
        it->currentindex = 0;
        it++;
    }

    // Now integrate
    dReal scur = 0, sdcur, sd, sdd, t;
    dReal dt = 0.01;
    dReal dtsq = dt*dt;

    dReal ks = 0;
    dReal ksd = 0;

    Profile profile;
    dReal tres;

    //if(ComputeLowestSd(scur,sdcur,sdd,profileslist)) {
    scur = 0.01;

    FindLowestProfile(scur,profile,tres,profileslist);
    sdcur = profile.Evald(tres);
    while(true) {
        //std::cout << scur << " " << sdcur << "   " << "\n";
        if(!FindLowestProfile(scur,profile,tres,profileslist)) {
            break;
        }
        sd = profile.Evald(tres);
        sdd = profile.Evaldd(tres);
        if(profile.Eval(tres+dt)<INF) {
            sdd += ks*(profile.Eval(tres+dt)-scur) + ksd*(profile.Evald(tres+dt)-sdcur);
        }
        svect.push_back(scur);
        sdvect.push_back(sdcur);
        sddvect.push_back(sdd);
        scur += sdcur*dt + 0.5*sdd*dtsq;
        sdcur += sdd*dt;
        //std::cout << scur << " " << sdcur << "   " << "\n\n";

    }
    //}

    nsteps = svect.size();
    duration = dt * (nsteps-1);
}


bool Profile::FindTimestepIndex(dReal t, int &index, dReal& remainder){
    if (t < -TINY || t > duration+TINY) {
        return false;
    }
    if(t<0) {
        t = 0;
    }
    if(t>duration) {
        t = duration;
    }
    if(duration-t <= TINY) {
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
        assert(SolveQuadraticEquation(svect[currentindex]-s,sdvect[currentindex],0.5*sddvect[currentindex],tres,0,integrationtimestep));
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
        assert(SolveQuadraticEquation(svect[currentindex-1]-s,sdvect[currentindex-1],0.5*sddvect[currentindex-1],tres,0,integrationtimestep));
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

void Profile::Print(){
    for(dReal t=0; t<=duration; t+=integrationtimestep) {
        std::cout<< Eval(t) << " ; " << Evald(t) << " ; " << Evaldd(t) <<  "\n";
    }
}



void Profile::Write(std::stringstream& ss, dReal dt){
    ss << duration << " " << dt << "\n";
    for(dReal t=0; t<=duration; t+=dt) {
        ss << Eval(t) << " ";
    }
    ss << "\n";
    for(dReal t=0; t<=duration; t+=dt) {
        ss << Evald(t) << " ";
    }
}




////////////////////////////////////////////////////////////////////
//////////////////////// Integration ///////////////////////////////
////////////////////////////////////////////////////////////////////


int IntegrateForward(Constraints& constraints, dReal sstart, dReal sdstart, Profile& profile, int maxsteps, std::list<Profile>&testprofileslist){
    dReal dt = constraints.tunings.integrationtimestep;
    dReal dtsq = dt*dt;
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
            dReal snext = scur + dt * sdcur + 0.5*dtsq*beta;
            scur = snext;
            sdcur = sdnext;
        }
    }

    profile = Profile(slist,sdlist,sddlist,dt);
    profile.forward = true;
    return returntype;
}



int IntegrateBackward(Constraints& constraints, dReal sstart, dReal sdstart, Profile& profile, int maxsteps, std::list<Profile>&testprofileslist){
    dReal dt = constraints.tunings.integrationtimestep;
    dReal dtsq = dt*dt;
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
            dReal snext = scur - dt * sdcur + 0.5*dtsq*alpha;
            scur = snext;
            sdcur = sdnext;
        }
    }
    slist.reverse();
    sdlist.reverse();
    sddlist.reverse();
    sddlist.pop_front();
    sddlist.push_back(0);
    profile = Profile(slist,sdlist,sddlist,dt);
    profile.forward = false;
    return returntype;
}


bool PassSwitchPoint(Constraints& constraints, dReal s, dReal sd){
    int ret;
    Profile profile;
    ret = IntegrateBackward(constraints,s,sd,profile,constraints.tunings.passswitchpointnsteps);
    if(ret==IRT_MAXSTEPS||ret==IRT_END) {
        ret = IntegrateForward(constraints,s,sd,profile,constraints.tunings.passswitchpointnsteps);
        if(ret==IRT_MAXSTEPS||ret==IRT_END) {
            return true;
        }
    }
    return false;
}


dReal BisectionSearch(Constraints& constraints, dReal s, dReal sdbottom, dReal sdtop, int position){
    if(position!=1 && PassSwitchPoint(constraints,s,sdtop)) {
        return sdtop;
    }
    if(sdtop-sdbottom<constraints.tunings.threshold) {
        if(position!=-1 && PassSwitchPoint(constraints,s,sdbottom)) {
            return sdbottom;
        }
        return -1;
    }
    dReal sdmid = (sdbottom+sdtop)*0.5;
    return std::max(BisectionSearch(constraints,s,sdbottom,sdmid,-1),BisectionSearch(constraints,s,sdmid,sdtop,1));
}



bool AddressSwitchPoint(Constraints& constraints, const SwitchPoint &switchpoint, dReal& sbackward, dReal& sdbackward, dReal& sforward, dReal& sdforward){
    dReal s = switchpoint.s;
    dReal sd = switchpoint.sd;

    //    if(switchpoint.switchpointtype == SPT_TANGENT || switchpoint.switchpointtype == SPT_DISCONTINUOUS) {
    if(true) {
        dReal sdtop = BisectionSearch(constraints,s,0,sd,0);
        if(sdtop<=0) {
            return false;
        }
        sbackward = s;
        sforward = s;
        sdbackward = sdtop;
        sdforward = sdtop;
        return true;
    }
    else{
        // here switchpointtype == SP_SINGULAR

    }
}



int ComputeLimitingCurves(Constraints& constraints, std::list<Profile>&resprofileslist){
    std::list<SwitchPoint> switchpointslist0(constraints.switchpointslist);
    Profile tmpprofile;
    dReal sswitch, sdswitch, sbackward, sdbackward, sforward, sdforward;
    int status = CLC_OK;
    int integratestatus;

    while(switchpointslist0.size()>0) {
        SwitchPoint switchpoint = switchpointslist0.front();
        switchpointslist0.pop_front();
        sswitch = switchpoint.s;
        sdswitch = switchpoint.sd;
        if(IsAboveProfilesList(sswitch,sdswitch,resprofileslist,false,true)) {
            continue;
        }

        // Address Switch Point
        if(!AddressSwitchPoint(constraints,switchpoint,sbackward,sdbackward,sforward,sdforward)) {
            status = CLC_SP;
            break;
        }

        // Integrate backward
        integratestatus = IntegrateBackward(constraints,sbackward,sdbackward,tmpprofile);
        if(tmpprofile.nsteps>1) {
            resprofileslist.push_back(tmpprofile);
        }

        // Integrate forward
        integratestatus = IntegrateForward(constraints,sforward,sdforward,tmpprofile);
        if(tmpprofile.nsteps>1) {
            resprofileslist.push_back(tmpprofile);
        }

    }
}



int IntegrateAllProfiles(Constraints& constraints, std::list<Profile>&resprofileslist){

}




////////////////////////////////////////////////////////////////////
///////////////////////// Utilities ////////////////////////////////
////////////////////////////////////////////////////////////////////


void VectorFromString(const std::string& s,std::vector<dReal>& resvect){
    std::istringstream iss(s);
    std::string sub;
    dReal value;
    resvect.resize(0);
    while(iss.good()) {
        iss >> sub;
        value = atof(sub.c_str());
        resvect.push_back(value);
    }
}


bool SolveQuadraticEquation(dReal a0, dReal a1, dReal a2, dReal& sol, dReal lowerbound, dReal upperbound){
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
        if(std::abs(a1)<=TINY) {
            return false; // Must be in a weird case
        }
        sol = -a0/a1;
        return sol>=lowerbound-TINY && sol<=upperbound+TINY;
    }
    sol = (-a1-sqrt(delta))/(2*a2);
    if(sol>=lowerbound-TINY && sol<=upperbound+TINY) {
        return true;
    }
    sol = (-a1+sqrt(delta))/(2*a2);
    return sol>=lowerbound-TINY && sol<=upperbound+TINY;
}


bool IsAboveProfilesList(dReal s, dReal sd, std::list<Profile>&testprofileslist, bool searchbackward, bool reinitialize){
    dReal t;
    std::list<Profile>::iterator it = testprofileslist.begin();
    while(it != testprofileslist.end()) {
        if(reinitialize) {
            if(searchbackward) {
                it->currentindex = it->nsteps-1;
            }
            else{
                it->currentindex = 0;
            }
        }
        if(it->Invert(s,t,searchbackward)) {
            if(sd > it->Evald(t)) {
                return true;
            }
        }
        it++;
    }
    return false;
}


bool FindLowestProfile(dReal s, Profile& profile, dReal& tres, std::list<Profile>&testprofileslist){
    dReal t;
    dReal sdmin;
    dReal sdtmp,sdd;
    int i = 0;
    int finali = 0;
    sdmin = INF;
    dReal snext,sdnext,snext2,sdnext2;

    dReal dt = 0.01;

    std::list<Profile>::iterator it = testprofileslist.begin();
    while(it != testprofileslist.end()) {
        if(it->Invert(s,t)) {
            sdtmp = it->Evald(t);
            if(sdtmp < sdmin) {
                sdmin = sdtmp;
                profile = *it;
                tres = t;
                sdd = it->Evaldd(t);
                snext = it->Eval(t+dt);
                sdnext = it->Evald(t+dt);
                snext2 = s + dt*sdmin + 0.5*sdd*dt*dt;
                sdnext2 = sdmin + dt*sdd;
                finali = i;
            }
        }
        it++;
        i++;
    }
    //std::cout << s << "/" << snext << "/" << snext2 << "  " << sdmin << "/" << sdnext << "/" << sdnext2 << " index: " << finali << "\n";
    return sdmin < INF;
}



} // end namespace TOPP

