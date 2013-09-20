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
/////////////////////////// Tunings ////////////////////////////////
////////////////////////////////////////////////////////////////////

Tunings::Tunings(const std::string& tuningsstring){
    std::istringstream iss(tuningsstring);
    std::string sub;
    iss >> sub;
    discrtimestep = atof(sub.c_str());
    iss >> sub;
    integrationtimestep = atof(sub.c_str());
    iss >> sub;
    sdprecision = atof(sub.c_str());
    iss >> sub;
    passswitchpointnsteps = atof(sub.c_str());
    iss >> sub;
    reparamtimestep = atof(sub.c_str());
}





////////////////////////////////////////////////////////////////////
/////////////////////////// Constraints ////////////////////////////
////////////////////////////////////////////////////////////////////

void Constraints::Preprocess(Trajectory& trajectory0, const Tunings& tunings0){
    trajectory = trajectory0;
    tunings = tunings0;
    Discretize();
    DiscretizeDynamics();
    ComputeMVCBobrow();
    ComputeMVCDirect();
    FindSwitchPoints();
}


void Constraints::Discretize(){
    ndiscrsteps = int((trajectory.duration+TINY)/tunings.discrtimestep);
    ndiscrsteps++;
    discrsvect.resize(0);
    for(int i=0; i<ndiscrsteps; i++) {
        discrsvect.push_back(i*tunings.discrtimestep);
    }
}


void Constraints::ComputeMVCBobrow(){
    for(int i=0; i<ndiscrsteps; i++) {
        mvcdirect.push_back(SdLimitBobrow(discrsvect[i]));
    }
}

void Constraints::ComputeMVCDirect(){
    for(int i=0; i<ndiscrsteps; i++) {
        mvcdirect.push_back(SdLimitDirect(discrsvect[i]));
    }
}



void Constraints::WriteMVCBobrow(std::stringstream& ss, dReal dt){
    dReal duration = trajectory.duration;
    ss << duration << " " << dt << "\n";
    for(dReal t=0; t<=duration; t+=dt) {
        ss << t << " ";
    }
    ss << "\n";
    for(dReal t=0; t<=duration; t+=dt) {
        ss << SdLimitBobrow(t) << " ";
    }
}


void Constraints::FindSwitchPoints(){
    FindSingularSwitchPoints();
    FindTangentSwitchPoints();
    FindDiscontinuousSwitchPoints();
}

void Constraints::AddSwitchPoint(int i, int switchpointtype){
    int iadd = i+1;
    if(mvcbobrow[i-1]<std::min(mvcbobrow[i],mvcbobrow[i+1])) {
        iadd = i-1;
    }
    else{
        if(mvcbobrow[i]<mvcbobrow[i+1]) {
            iadd = i;
        }
    }
    std::list<SwitchPoint>::iterator it = switchpointslist.begin();
    dReal s = discrsvect[iadd];
    dReal sd = mvcbobrow[iadd];
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
    sd = SdLimitBobrow(s);
    sdnext = SdLimitBobrow(snext);
    sddlimits = SddLimits(s,sd);
    alpha = sddlimits.first;
    beta = sddlimits.second;
    diffprev = alpha/sd - (sdnext-sd)/tunings.discrtimestep;

    for(int i=2; i<ndiscrsteps-1; i++) {
        s = discrsvect[i];
        snext = discrsvect[i+1];
        sd = SdLimitBobrow(s);
        sdnext = SdLimitBobrow(snext);
        sddlimits = SddLimits(s,sd);
        alpha = sddlimits.first;
        beta = sddlimits.second;
        diff = alpha/sd - (sdnext-sd)/tunings.discrtimestep;
        if(diffprev*diff<0) {
            AddSwitchPoint(i,SP_TANGENT);
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
        index = int(t/integrationtimestep);
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


int IntegrateForward(Constraints& constraints, dReal sstart, dReal sdstart, dReal dt,  Profile& resprofile, int maxsteps, std::list<Profile>&testprofileslist, bool testmvc){
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
            returntype = INT_MAXSTEPS;
            break;
        }
        else if(scur > constraints.trajectory.duration) {
            //TODO: change the time step of previous step to reach the end
            returntype = INT_END;
            break;
        }
        else if(sdcur < 0) {
            //TODO: double check whether alpha>beta
            returntype = INT_BOTTOM;
            break;
        }
        else if(testmvc && sdcur > std::min(constraints.SdLimitDirect(scur),constraints.SdLimitBobrow(scur))) {
            if(constraints.SdLimitDirect(scur)>=constraints.SdLimitBobrow(scur)) {
                returntype = INT_MVCBOBROW;
                break;
            }
            //Here we have sddirect < sdcur < sdbobrow
            bool cont2 = true;
            while(cont2) {
                // Last step before going above MVCDirect
                if(slist.size()>=1) {
                    scur = slist.back();
                    sdcur = sdlist.back();
                }
                std ::pair<dReal,dReal> sddlimits = constraints.SddLimits(scur,sdcur);
                dReal alpha = sddlimits.first, beta = sddlimits.second;
                dReal snextalpha, sdnextalpha, sdnextdirectalpha, snextbeta, sdnextbeta, sdnextdirectbeta;
                sdnextalpha = sdcur + dt * alpha;
                snextalpha = scur + dt * sdcur + 0.5*dtsq*alpha;
                sdnextdirectalpha = constraints.SdLimitDirect(snextalpha);
                sdnextbeta = sdcur + dt * beta;
                snextbeta = scur + dt * sdcur + 0.5*dtsq*beta;
                sdnextdirectbeta = constraints.SdLimitDirect(snextbeta);
                if(sdnextbeta<=sdnextdirectbeta) {
                    // Beta points below MVCDirect
                    // Follow beta
                    scur = snextbeta;
                    sdcur = snextbeta;
                    slist.push_back(scur);
                    sdlist.push_back(sdcur);
                    sddlist.push_back(beta);
                    break;
                }
                else if(sdnextalpha>=sdnextdirectalpha) {
                    // Alpha points above MVCDirect
                    // Follow MVCDirect until alpha points below MVCDirect
                    // Then add this point to the switchpointslist
                    while(scur+dt<=constraints.trajectory.duration) {
                        scur += dt;
                        sdcur = constraints.SdLimitDirect(scur);
                        std ::pair<dReal,dReal> sddlimits2 = constraints.SddLimits(scur,sdcur);
                        dReal alpha2 = sddlimits2.first;
                        dReal snextalpha2, sdnextalpha2, sdnextdirectalpha2;
                        sdnextalpha2 = sdcur + dt * alpha;
                        snextalpha2 = scur + dt * sdcur + 0.5*dtsq*alpha2;
                        sdnextdirectalpha2 = constraints.SdLimitDirect(snextalpha2);
                        if(sdnextalpha<sdnextdirectalpha) {
                            // Add point to switchpointlist
                            constraints.zlajpahlist.push_back(SwitchPoint(scur,sdcur,SP_ZLAJPAH));
                            cont = false;
                            cont2 = false;
                            break;
                        }
                    }
                }
                else{
                    // MVCDirect is between beta and alpha
                    // Slide along MVCDirect
                }
            }
        }
        else if(IsAboveProfilesList(scur,sdcur,testprofileslist)) {
            slist.push_back(scur);
            sdlist.push_back(sdcur);
            sddlist.push_back(0);
            returntype = INT_PROFILE;
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

    resprofile = Profile(slist,sdlist,sddlist,dt);
    resprofile.forward = true;
    return returntype;
}



int IntegrateBackward(Constraints& constraints, dReal sstart, dReal sdstart, dReal dt, Profile& resprofile,  int maxsteps, std::list<Profile>&testprofileslist, bool testmvc){
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
            returntype = INT_MAXSTEPS;
            break;
        }
        else if(scur < 0) {
            //TODO: change the time step of previous step to reach the end
            returntype = INT_END;
            break;
        }
        else if(sdcur < 0) {
            //TODO: double check whether alpha>beta
            returntype = INT_BOTTOM;
            break;
        }
        else if(testmvc && sdcur > constraints.SdLimitBobrow(scur)) {
            returntype = INT_MVCBOBROW;
            break;
        }
        else if(IsAboveProfilesList(scur,sdcur,testprofileslist,searchbackward)) {
            slist.push_back(scur);
            sdlist.push_back(sdcur);
            sddlist.push_back(0);
            returntype = INT_PROFILE;
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
    if(slist.size()>0) {
        slist.reverse();
        sdlist.reverse();
        sddlist.reverse();
        sddlist.pop_front();
        sddlist.push_back(0);
    }
    resprofile = Profile(slist,sdlist,sddlist,dt);
    resprofile.forward = false;
    return returntype;
}


bool PassSwitchPoint(Constraints& constraints, dReal s, dReal sd, dReal dt){
    int ret;
    Profile resprofile;
    ret = IntegrateBackward(constraints,s,sd,dt,resprofile,constraints.tunings.passswitchpointnsteps);
    if(ret==INT_MAXSTEPS||ret==INT_END) {
        ret = IntegrateForward(constraints,s,sd,dt,resprofile,constraints.tunings.passswitchpointnsteps);
        if(ret==INT_MAXSTEPS||ret==INT_END) {
            return true;
        }
    }
    return false;
}


dReal BisectionSearch(Constraints& constraints, dReal s, dReal sdbottom, dReal sdtop, dReal dt, int position){
    if(position!=1 && PassSwitchPoint(constraints,s,sdtop,dt)) {
        return sdtop;
    }
    if(sdtop-sdbottom<constraints.tunings.sdprecision) {
        if(position!=-1 && PassSwitchPoint(constraints,s,sdbottom,dt)) {
            return sdbottom;
        }
        return -1;
    }
    dReal sdmid = (sdbottom+sdtop)*0.5;
    return std::max(BisectionSearch(constraints,s,sdbottom,sdmid,dt,-1),BisectionSearch(constraints,s,sdmid,sdtop,dt,1));
}



bool AddressSwitchPoint(Constraints& constraints, const SwitchPoint &switchpoint, dReal& sbackward, dReal& sdbackward, dReal& sforward, dReal& sdforward){
    dReal s = switchpoint.s;
    dReal sd = switchpoint.sd;
    dReal discr = constraints.tunings.discrtimestep;
    dReal dt;
    int ret;
    Profile resprofile;

    if(switchpoint.switchpointtype == SP_TANGENT || switchpoint.switchpointtype == SP_DISCONTINUOUS || switchpoint.switchpointtype == SP_ZLAJPAH) {
        dReal sdtop = BisectionSearch(constraints,s,0,sd,dt,0);
        dt = constraints.tunings.discrtimestep/2;
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
        dt = discr/2;
        ret = IntegrateBackward(constraints,s-discr,sd*0.9,dt,resprofile,constraints.tunings.passswitchpointnsteps,voidprofileslist,false);
        if(ret==INT_MAXSTEPS||ret==INT_END) {
            sbackward = resprofile.Eval(0);
            sdbackward = resprofile.Evald(0);
            if(sdbackward>constraints.SdLimitBobrow(sbackward)) {
                return false;
            }
            ret = IntegrateForward(constraints,s+discr,sd*0.9,dt,resprofile,constraints.tunings.passswitchpointnsteps,voidprofileslist,false);
            if(ret==INT_MAXSTEPS||ret==INT_END) {
                sforward = resprofile.Eval(resprofile.duration);
                sdforward = resprofile.Evald(resprofile.duration);
                if(sdforward>constraints.SdLimitBobrow(sforward)) {
                    return false;
                }
                return true;
            }
        }
    }
    return false;
}



int ComputeLimitingCurves(Constraints& constraints, std::list<Profile>&resprofileslist, bool zlajpah){
    std::list<SwitchPoint>* switchpointslist0;
    if(zlajpah) {
        *switchpointslist0 = std::list<SwitchPoint>(constraints.switchpointslist);
    }
    else{
        switchpointslist0 = &(constraints.zlajpahlist);
    }

    Profile tmpprofile;
    dReal sswitch, sdswitch, sbackward, sdbackward, sforward, sdforward;
    int integratestatus;

    while(switchpointslist0->size()>0) {
        SwitchPoint switchpoint = switchpointslist0->front();
        switchpointslist0->pop_front();
        sswitch = switchpoint.s;
        sdswitch = switchpoint.sd;
        if(IsAboveProfilesList(sswitch,sdswitch,resprofileslist,false,true)) {
            continue;
        }

        // Address Switch Point
        if(!AddressSwitchPoint(constraints,switchpoint,sbackward,sdbackward,sforward,sdforward)) {
            continue;
        }

        // Add middle part
        if(sforward-sbackward>TINY) {
            std::list<dReal> slist, sdlist, sddlist;
            dReal dtmod = 2*(sforward-sbackward)/(sdforward+sdbackward);
            dReal sdd = (sdforward-sdbackward)/dtmod;
            slist.push_back(sbackward);
            slist.push_back(sforward);
            sdlist.push_back(sdbackward);
            sdlist.push_back(sdforward);
            sddlist.push_back(sdd);
            sddlist.push_back(0);
            resprofileslist.push_back(Profile(slist,sdlist,sddlist,dtmod));
        }

        // Integrate backward
        integratestatus = IntegrateBackward(constraints,sbackward,sdbackward,constraints.tunings.integrationtimestep,tmpprofile);
        if(tmpprofile.nsteps>1) {
            resprofileslist.push_back(tmpprofile);
        }
        if(integratestatus == INT_BOTTOM) {
            return CLC_BOTTOM;
        }

        // Integrate forward
        integratestatus = IntegrateForward(constraints,sforward,sdforward,constraints.tunings.integrationtimestep,tmpprofile);
        if(tmpprofile.nsteps>1) {
            resprofileslist.push_back(tmpprofile);
        }
        if(integratestatus == INT_BOTTOM) {
            return CLC_BOTTOM;
        }
    }
    return CLC_OK;
}



int PP(Constraints& constraints, Trajectory& trajectory, Tunings& tunings, dReal sdbeg, dReal sdend, Trajectory& restrajectory, std::list<Profile>& resprofileslist){
    constraints.Preprocess(trajectory,tunings);
    if(VectorMin(constraints.mvcbobrow) <= TINY) {
        std::cout << "MVCBobrow hit 0\n";
        return 0;
    }
    Profile resprofile;
    int ret;
    ret = ComputeLimitingCurves(constraints,resprofileslist);
    if(ret!=CLC_OK) {
        std::cout << "CLC failed\n";
        return 0;
    }
    // Zlajpah
    while(constraints.zlajpahlist.size()>0) {
        int ret = ComputeLimitingCurves(constraints,resprofileslist,true);
        if(ret!=CLC_OK) {
            std::cout << "CLC failed\n";
            return 0;
        }
    }


    ret = IntegrateForward(constraints,0,sdbeg,constraints.tunings.integrationtimestep,resprofile,1e5,resprofileslist);
    if(ret==INT_BOTTOM) {
        std::cout << "FW reached 0\n";
        return 0;
    }
    resprofileslist.push_back(resprofile);
    // Zlajpah
    while(constraints.zlajpahlist.size()>0) {
        int ret = ComputeLimitingCurves(constraints,resprofileslist,true);
        if(ret!=CLC_OK) {
            std::cout << "CLC failed\n";
            return 0;
        }
    }

    ret = IntegrateBackward(constraints,trajectory.duration,sdend,constraints.tunings.integrationtimestep,resprofile,1e5,resprofileslist);
    if(ret==INT_BOTTOM) {
        std::cout << "BW reached 0\n";
        return 0;
    }
    resprofileslist.push_back(resprofile);
    trajectory.Reparameterize(resprofileslist,tunings.reparamtimestep,restrajectory);
    return 1;
}



int VIP(Constraints& constraints, Trajectory& trajectory, Tunings& tunings, dReal sdbegmin, dReal sdbegmax, dReal& sdendmin, dReal& sdendmax, std::list<Profile>&resprofileslist){
    constraints.Preprocess(trajectory,tunings);
    Profile tmpprofile;
    dReal tres;
    dReal smallincrement = constraints.tunings.integrationtimestep*2;

    // Compute the limiting curves
    int resclc = ComputeLimitingCurves(constraints,resprofileslist);
    if(resclc == CLC_SWITCH || resclc == CLC_BOTTOM) {
        return 0;
    }

    // Determine the lowest profile at t=0
    FindLowestProfile(smallincrement,tmpprofile,tres,resprofileslist);
    dReal bound = std::min(tmpprofile.Evald(tres),constraints.mvcbobrow[0]);
    if(sdbegmin>bound) {
        return 0;
    }
    sdbegmax = std::min(sdbegmax,bound);


    // Compute sdendmax by integrating forward from (0,sdbegmax)
    int resintfw = IntegrateForward(constraints,0,sdbegmax,constraints.tunings.integrationtimestep,tmpprofile,1e5,resprofileslist);
    resprofileslist.push_back(tmpprofile);
    if(resintfw == INT_BOTTOM || resintfw == INT_MVCBOBROW) {
        return 0;
    }
    if(resintfw == INT_END) {
        sdendmax = tmpprofile.Evald(tmpprofile.duration);
    }
    else if(resintfw == INT_PROFILE) {
        // Look for the lowest profile at the end
        if(FindLowestProfile(trajectory.duration-smallincrement,tmpprofile,tres,resprofileslist)) {
            sdendmax = tmpprofile.Evald(tres);
        }
        else{
            // No profile reaches the end, consider the MVC BOBROW instead
            sdendmax = constraints.mvcbobrow[constraints.mvcbobrow.size()-1];
            int resintbw = IntegrateBackward(constraints,trajectory.duration,sdendmax,constraints.tunings.integrationtimestep,tmpprofile,1e5,resprofileslist);
            if(resintbw == INT_BOTTOM || resintbw == INT_MVCBOBROW) {
                return 0;
            }
        }
    }

    // Integrate from (send,0)
    int resintbw = IntegrateBackward(constraints,trajectory.duration,0,constraints.tunings.integrationtimestep,tmpprofile,1e5,resprofileslist);
    if((resintbw == INT_END && tmpprofile.Evald(0)>=sdbegmin) || resintbw == INT_PROFILE) {
        resprofileslist.push_back(tmpprofile);
        sdendmin = 0;
        return 1;
    }

    dReal sdupper = sdendmax, sdlower = 0;
    Profile bestprofile;
    while(sdupper-sdlower > tunings.sdprecision) {
        dReal sdtest = (sdupper + sdlower)/2;
        int resintbw2 = IntegrateBackward(constraints,trajectory.duration,sdtest,constraints.tunings.integrationtimestep,tmpprofile,1e5,resprofileslist);
        if((resintbw2 == INT_END && tmpprofile.Evald(0)>=sdbegmin) || resintbw2 == INT_PROFILE) {
            sdupper = sdtest;
            bestprofile = tmpprofile;
        }
        else{
            sdlower = sdtest;
        }
    }
    sdendmin = sdupper;
    resprofileslist.push_back(bestprofile);

    return 1;
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


dReal VectorMin(const std::vector<dReal>& v){
    std::vector<dReal>::const_iterator it = v.begin();
    dReal res = INF;
    while(it!=v.end()) {
        res = std::min(res,*it);
        it++;
    }
    return res;
}


dReal VectorMax(const std::vector<dReal>& v){
    std::vector<dReal>::const_iterator it = v.begin();
    dReal res = -INF;
    while(it!=v.end()) {
        res = std::max(res,*it);
        it++;
    }
    return res;
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
    return sdmin < INF;
}



} // end namespace TOPP

