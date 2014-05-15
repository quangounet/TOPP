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


bool Constraints::Preprocess() {
    resprofileslist.resize(0);
    // Change discrtimestep so as it becomes a divisor of trajectory duration
    int ndiscrsteps = int((trajectory.duration+1e-10)/discrtimestep);
    if(ndiscrsteps<1) {
        return false;
    }

    discrtimestep = trajectory.duration/ndiscrsteps;
    Discretize();

    std::chrono::time_point<std::chrono::system_clock> t0,t1,t2,t3;
    std::chrono::duration<double> d1,d2,d3;

    t0 = std::chrono::system_clock::now();
    ComputeMVCBobrow();
    t1 = std::chrono::system_clock::now();
    ComputeMVCCombined();
    t2 = std::chrono::system_clock::now();
    FindSwitchPoints();
    t3 = std::chrono::system_clock::now();

    //d1 = t1-t0;
    //d2 = t2-t1;
    //d3 = t3-t2;

    //std::cout << d1.count() <<  " " << d2.count() << " " << d3.count() << "\n";

    if(passswitchpointnsteps == 0) {
        passswitchpointnsteps = 5;
    }

    // Set integration timestep automatically if it is initially set to 0
    dReal meanmvc = 0;
    if(integrationtimestep == 0) {
        for(int i=0; i<ndiscrsteps; i++) {
            meanmvc += std::min(mvccombined[i],10.);
        }
        meanmvc /= ndiscrsteps;
        meanmvc = std::min(1.,meanmvc);
        integrationtimestep = discrtimestep/meanmvc;
        //std::cout << "\n--------------\nIntegration timestep: " << integrationtimestep << "\n";
    }

    return true;
}


void Constraints::Discretize() {
    ndiscrsteps = int((trajectory.duration+TINY)/discrtimestep);
    ndiscrsteps++;
    discrsvect.resize(0);
    for(int i=0; i<ndiscrsteps; i++) {
        discrsvect.push_back(i*discrtimestep);
    }
}


void Constraints::ComputeMVCBobrow() {
    for(int i=0; i<ndiscrsteps; i++) {
        mvcbobrow.push_back(SdLimitBobrowInit(discrsvect[i]));
    }
}


void Constraints::ComputeMVCCombined(){
    for(int i=0; i<ndiscrsteps; i++) {
        mvccombined.push_back(SdLimitCombinedInit(discrsvect[i]));
    }
}


dReal Constraints::Interpolate1D(dReal s, const std::vector<dReal>& v) {
    assert(s>=-TINY && s<=trajectory.duration+TINY);
    if(s<0) {
        s=0;
    }
    if(s>=trajectory.duration) {
        int n = ndiscrsteps-1;
        return v[n];
    }
    int n = int(s/discrtimestep);
    dReal coef = (s-n*discrtimestep)/discrtimestep;
    return (1-coef)*v[n] + coef*v[n+1];
}


dReal Constraints::SdLimitCombinedInit(dReal s){
    dReal res = SdLimitBobrow(s);
    std::vector<dReal> qd(trajectory.dimension);
    if(hasvelocitylimits) {
        trajectory.Evald(s, qd);
        for(int i=0; i<trajectory.dimension; i++) {
            if(std::abs(qd[i])>TINY) {
                res = std::min(res,vmax[i]/std::abs(qd[i]));
            }
        }
    }
    return res;
}


dReal Constraints::SdLimitCombined(dReal s) {
    return Interpolate1D(s,mvccombined);
}


dReal Constraints::SdLimitBobrow(dReal s) {
    return Interpolate1D(s,mvcbobrow);
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


void Constraints::WriteMVCDirect(std::stringstream& ss, dReal dt){
    std::vector<dReal> qd(trajectory.dimension);
    dReal duration = trajectory.duration;
    ss << duration << " " << dt << "\n";
    for(dReal t=0; t<=duration; t+=dt) {
        ss << t << " ";
    }
    ss << "\n";
    for(dReal t=0; t<=duration; t+=dt) {
        dReal res = INF;
        trajectory.Evald(t, qd);
        for(int i=0; i<trajectory.dimension; i++) {
            if(std::abs(qd[i])>TINY) {
                res = std::min(res,vmax[i]/std::abs(qd[i]));
            }
        }
        ss << res << " ";
    }
}


void Constraints::FindSwitchPoints() {
    FindSingularSwitchPoints();
    FindTangentSwitchPoints();
    FindDiscontinuousSwitchPoints();
    TrimSwitchPoints();
}


void Constraints::AddSwitchPoint(int i, int switchpointtype, dReal sd){
    dReal s = discrsvect[i];
    // If no sd is specified, then take the value of the mvc
    // (The case when sd is specified corresponds to a singular switchpoint in some cases)
    if(sd<0) {
        sd = mvcbobrow[i];
    }
    if(sd > MAXSD) {
        return;
    }
    std::list<SwitchPoint>::iterator it = switchpointslist.begin();
    while(it!=switchpointslist.end()) {
        if(s == it->s) {
            return;
        }
        if(s<=it->s) {
            break;
        }
        it++;
    }
    SwitchPoint sw(s,sd,switchpointtype);
    if(switchpointtype == SP_SINGULAR) {
        sw.slopesvector = std::vector<dReal>();
        ComputeSlopeDynamicSingularity(s,sd,sw.slopesvector);
    }
    switchpointslist.insert(it,sw);
}

void Constraints::FindTangentSwitchPoints(){
    if(ndiscrsteps<3)
        return;
    int i = 1;
    dReal s,sd,snext,sdnext,alpha,diff,diffprev,tangent,prevtangent;
    std::pair<dReal,dReal> sddlimits;

    s = discrsvect[i];
    snext = discrsvect[i+1];
    sd = SdLimitBobrow(s);
    sdnext = SdLimitBobrow(snext);
    tangent = (sdnext-sd)/discrtimestep;
    prevtangent = (sd - SdLimitBobrow(discrsvect[i-1]))/discrtimestep;
    sddlimits = SddLimits(s,sd);
    alpha = sddlimits.first;
    //beta = sddlimits.second;
    diffprev = alpha/sd - tangent;

    for(int i=2; i<ndiscrsteps-1; i++) {
        s = discrsvect[i];
        snext = discrsvect[i+1];
        sd = SdLimitBobrow(s);
        sdnext = SdLimitBobrow(snext);
        sddlimits = SddLimits(s,sd);
        alpha = sddlimits.first;
        if(std::abs(prevtangent-tangent)>2 && prevtangent < 0 && tangent >0) {
            AddSwitchPoint(i,SP_DISCONTINUOUS);
        }
        prevtangent = tangent;
        tangent = (sdnext-sd)/discrtimestep;
        //if(std::abs(tangent-prevtangent)>1.) {
        //    continue;
        //}
        //beta = sddlimits.second;
        diff = alpha/sd - tangent;
        if(diffprev*diff<0 && std::abs(diff)<1) {
            AddSwitchPoint(i,SP_TANGENT);
        }
        diffprev = diff;
    }
}

void Constraints::FindDiscontinuousSwitchPoints() {
    if(ndiscrsteps<3)
        return;
    int i = 0;
    dReal sd, sdn, sdnn;
    sd = SdLimitBobrow(discrsvect[i]);
    sdn = SdLimitBobrow(discrsvect[i+1]);
    for(int i=0; i<ndiscrsteps-3; i++) {
        sdnn = SdLimitBobrow(discrsvect[i+2]);
        if(std::abs(sdnn-sdn)>100*std::abs(sdn-sd)) {
            if(sdn<sdnn) {
                AddSwitchPoint(i+1,SP_DISCONTINUOUS);
            }
            else{
                AddSwitchPoint(i+2,SP_DISCONTINUOUS);
            }
        }
        sd = sdn;
        sdn = sdnn;
    }
}


void InsertInSpList(std::list<SwitchPoint>& splist, SwitchPoint sp){
    if(splist.size()==0) {
        splist.push_back(sp);
        return;
    }
    std::list<SwitchPoint>::iterator it = splist.begin();
    while(it!=splist.end()) {
        if(sp.switchpointtype == SP_SINGULAR) {
            if((it->switchpointtype == SP_SINGULAR && it->sd >= sp.sd) ||  it->switchpointtype!=SP_SINGULAR) {
                splist.insert(it,sp);
                return;
            }
        }
        else if(it->switchpointtype!=SP_SINGULAR && it->sd >= sp.sd ) {
            splist.insert(it,sp);
            return;
        }
        it++;

    }
    splist.push_back(sp);
}


void Constraints::TrimSwitchPoints() {
    dReal radius = discrtimestep*2.1;
    dReal scur = -1, snext, sdcur = -1, sdnext;
    int stypenext;
    std::list<SwitchPoint>::iterator itcur;
    std::vector<dReal> slopesvector;
    std::list<SwitchPoint>::iterator it = switchpointslist.begin();

    // Merge singular points
    while(it!=switchpointslist.end()) {
        if(it->switchpointtype == SP_SINGULAR) {
            scur = it->s;
            sdcur = it->sd;
            itcur = it;
            break;
        }
        it++;
    }
    if(scur>=0) {
        it++;
        while(it!=switchpointslist.end()) {
            snext = it->s;
            sdnext = it->sd;
            stypenext = it->switchpointtype;
            if(stypenext == SP_SINGULAR) {
                if(snext-scur<radius) {
                    if(sdcur < sdnext) {
                        slopesvector = it->slopesvector;
                        it =  switchpointslist.erase(it);
                    }
                    else{
                        slopesvector = itcur->slopesvector;
                        switchpointslist.erase(itcur);
                        scur = snext;
                        sdcur = sdnext;
                        itcur = it;
                        it++;
                    }
                    for(int j=0; j<int(slopesvector.size()); j++) {
                        itcur->slopesvector.push_back(slopesvector[j]);
                    }
                }
                else{
                    scur = snext;
                    sdcur = sdnext;
                    itcur = it;
                    it++;
                }
            }
            else{
                it++;
            }

        }
    }

    // Merge non singular points
    it = switchpointslist.begin();
    while(it!=switchpointslist.end()) {
        if(it->switchpointtype != SP_SINGULAR) {
            scur = it->s;
            sdcur = it->sd;
            itcur = it;
            break;
        }
        it++;
    }
    if(scur>=0) {
        it++;
        while(it!=switchpointslist.end()) {
            snext = it->s;
            sdnext = it->sd;
            stypenext = it->switchpointtype;
            if(stypenext != SP_SINGULAR) {
                if(snext-scur<radius) {
                    if(sdcur < sdnext) {
                        slopesvector = it->slopesvector;
                        it =  switchpointslist.erase(it);
                    }
                    else{
                        slopesvector = itcur->slopesvector;
                        switchpointslist.erase(itcur);
                        scur = snext;
                        sdcur = sdnext;
                        itcur = it;
                        it++;
                    }
                    for(int j=0; j<int(slopesvector.size()); j++) {
                        itcur->slopesvector.push_back(slopesvector[j]);
                    }
                }
                else{
                    scur = snext;
                    sdcur = sdnext;
                    itcur = it;
                    it++;
                }
            }
            else{
                it++;
            }

        }
    }

    // Sort by sd
    std::list<SwitchPoint> splist;
    it = switchpointslist.begin();
    while(it!=switchpointslist.end()) {
        InsertInSpList(splist,*it);
        it++;
    }
    switchpointslist = splist;
}



////////////////////////////////////////////////////////////////////
/////////////////// Quadratic Constraints //////////////////////////
////////////////////////////////////////////////////////////////////


QuadraticConstraints::QuadraticConstraints(const std::string& constraintsstring) {
    std::vector<dReal> tmpvect;
    std::string buff;
    std::istringstream iss(constraintsstring);
    getline(iss, buff, '\n');
    discrtimestep = atof(buff.c_str());
    getline(iss, buff, '\n');
    VectorFromString(buff, vmax);
    while(iss.good()) {
        getline(iss, buff, '\n');
        VectorFromString(buff, tmpvect);
        avect.push_back(tmpvect);
        getline(iss, buff, '\n');
        VectorFromString(buff,tmpvect);
        bvect.push_back(tmpvect);
        getline(iss, buff, '\n');
        VectorFromString(buff,tmpvect);
        cvect.push_back(tmpvect);
    }
    nconstraints = int(avect.front().size());
    hasvelocitylimits =  VectorMax(vmax) > TINY;
}


void QuadraticConstraints::CheckInput() {
    if ((int)vmax.size() != trajectory.dimension) {
        std::ostringstream msg;
        msg << "vmax has dimension " << vmax.size()
            << " but trajectory has dimension " << trajectory.dimension << ".";
        std::cout << "[TOPP] " << msg.str() << std::endl;
        throw TOPPException(msg.str());
    }
}


void QuadraticConstraints::InterpolateDynamics(dReal s, std::vector<dReal>& a, std::vector<dReal>& b, std::vector<dReal>& c) {
    a.resize(nconstraints);
    b.resize(nconstraints);
    c.resize(nconstraints);
    assert(s>=-TINY && s<=trajectory.duration+TINY);
    if(s < 0)
        s = 0;
    if(s >= trajectory.duration-TINY) {
        int n = ndiscrsteps-1;
        for(int i = 0; i < nconstraints; i++) {
            a[i]= avect[n][i];
            b[i]= bvect[n][i];
            c[i]= cvect[n][i];
        }
        return;
    }

    int n = int(s/discrtimestep);
    dReal coef = (s-n*discrtimestep)/discrtimestep;
    for (int i = 0; i < nconstraints; i++) {
        a[i] = (1-coef)*avect[n][i] + coef*avect[n+1][i];
        b[i] = (1-coef)*bvect[n][i] + coef*bvect[n+1][i];
        c[i] = (1-coef)*cvect[n][i] + coef*cvect[n+1][i];
    }
}

void QuadraticConstraints::FixStart(dReal& sstart,dReal& sdstart){
    std::vector<dReal> a, b, c, a2, b2, c2;
    sstart = 0;
    dReal delta=TINY2, ap, bp, cp, slope;
    InterpolateDynamics(0,a,b,c);
    InterpolateDynamics(delta,a2,b2,c2);
    dReal sdcurrent = INF;
    int indexcurrent = 0;
    for(int j=0; j<int(a.size()); j++) {
        if(std::abs(a[j])<TINY2) {
            if(c[j]*b[j]<0) {
                if(sqrt(-c[j]/b[j])<sdcurrent) {
                    sdcurrent = sqrt(-c[j]/b[j]);
                    indexcurrent = j;
                }
            }
        }
    }
    if(sdcurrent<1e3) {
        sstart = integrationtimestep;
        ap = (a2[indexcurrent]-a[indexcurrent])/delta;
        bp = (b2[indexcurrent]-b[indexcurrent])/delta;
        cp = (c2[indexcurrent]-c[indexcurrent])/delta;
        slope = (-bp*sdcurrent*sdcurrent-cp)/((2*b[indexcurrent]+ap)*sdcurrent);
        sdstart = sdcurrent+slope*sstart;
        resprofileslist.push_back(StraightProfile(0,sstart,sdcurrent,sdstart));
    }
}

void QuadraticConstraints::FixEnd(dReal& sendnew,dReal& sdendnew){
    std::vector<dReal> a, b, c, a2, b2, c2;
    dReal delta=TINY2, ap, bp, cp, slope;
    dReal send = trajectory.duration;
    sendnew = send;
    InterpolateDynamics(send,a,b,c);
    InterpolateDynamics(send-delta,a2,b2,c2);
    dReal sdcurrent = INF;
    int indexcurrent = 0;
    for(int j=0; j<int(a.size()); j++) {
        if(std::abs(a[j])<TINY2) {
            if(-c[j]*b[j]>0) {
                if(sqrt(-c[j]/b[j])<sdcurrent) {
                    sdcurrent = sqrt(-c[j]/b[j]);
                    indexcurrent = j;
                }
            }
        }
    }
    if(sdcurrent<1e3) {
        dReal stub = integrationtimestep;
        sendnew = send-stub;
        ap = (a[indexcurrent]-a2[indexcurrent])/delta;
        bp = (b[indexcurrent]-b2[indexcurrent])/delta;
        cp = (c[indexcurrent]-c2[indexcurrent])/delta;
        slope = (-bp*sdcurrent*sdcurrent-cp)/((2*b[indexcurrent]+ap)*sdcurrent);
        sdendnew = sdcurrent-slope*stub;
        resprofileslist.push_back(StraightProfile(sendnew,send,sdendnew,sdcurrent));
    }
}



void QuadraticConstraints::ComputeSlopeDynamicSingularity(dReal s, dReal sd, std::vector<dReal>& slopesvector) {
    dReal delta = TINY2, s2, ap, bp, cp, slope;
    std::vector<dReal> a, b, c, a2, b2, c2;
    if(s>delta) {
        delta = -delta;
    }
    s2 = s + delta;
    InterpolateDynamics(s,a,b,c);
    InterpolateDynamics(s2,a2,b2,c2);

    slopesvector.resize(0);
    for(int i=0; i<int(a.size()); i++) {
        ap = (a2[i]-a[i])/delta;
        bp = (b2[i]-b[i])/delta;
        cp = (c2[i]-c[i])/delta;
        slope = (-bp*sd*sd-cp)/((2*b[i]+ap)*sd);
        slopesvector.push_back(slope);
    }
}



std::pair<dReal,dReal> QuadraticConstraints::SddLimits(dReal s, dReal sd){
    dReal dtsq = integrationtimestep;
    dtsq = dtsq*dtsq;
    dReal safetybound = discrtimestep/dtsq;
    dReal alpha = -safetybound;
    dReal beta = safetybound;
    dReal sdsq = sd*sd;
    std::vector<dReal> a, b, c;

    dReal alpha_i, beta_i;
    InterpolateDynamics(s,a,b,c);

    for(int i=0; i<nconstraints; i++) {
        if(std::abs(a[i])<TINY) {
            if(b[i]*sdsq+c[i]>0) {
                // Constraint not satisfied
                beta = -INF;
                alpha = INF;
            }
            continue;
        }
        if(a[i]>0) {
            beta_i = (-sdsq*b[i]-c[i])/a[i];
            beta = std::min(beta,beta_i);
        }
        else{
            alpha_i = (-sdsq*b[i]-c[i])/a[i];
            alpha = std::max(alpha,alpha_i);
        }
    }
    std::pair<dReal,dReal> result(alpha,beta);
    return result;
}


dReal QuadraticConstraints::SdLimitBobrowInit(dReal s){
    std::vector<dReal> a, b, c;
    InterpolateDynamics(s,a,b,c);
    if(VectorNorm(a)<TINY) {
        if(s<1e-2) {
            s+=1e-3;
        }
        else{
            s-=1e-3;
        }
        InterpolateDynamics(s,a,b,c);
    }
    std::pair<dReal,dReal> sddlimits = SddLimits(s,0);
    if(sddlimits.first > sddlimits.second) {
        return 0;
    }

    dReal sdmin = INF;
    for(int k=0; k<nconstraints; k++) {
        for(int m=k+1; m<nconstraints; m++) {
            dReal num, denum, r;
            // If we have a pair of alpha and beta bounds, then determine the sdot for which the two bounds are equal
            if(a[k]*a[m]<0) {
                num = a[k]*c[m]-a[m]*c[k];
                denum = -a[k]*b[m]+a[m]*b[k];
                if(std::abs(denum)>TINY) {
                    r = num/denum;
                    if(r>=0) {
                        sdmin = std::min(sdmin,sqrt(r));
                    }
                }
            }
        }
    }
    return sdmin;
}


void QuadraticConstraints::FindSingularSwitchPoints() {
    if(ndiscrsteps<3) {
        return;
    }
    int i = 0;
    std::vector<dReal> a,aprev,b,c;

    InterpolateDynamics(discrsvect[i],aprev,b,c);

    for(int i=1; i<ndiscrsteps-1; i++) {
        InterpolateDynamics(discrsvect[i],a,b,c);
        dReal minsd = mvcbobrow[i];
        bool found = false;
        for(int j=0; j< int(a.size()); j++) {
            if(a[j]*aprev[j]<0) {
                if(c[j]/b[j]<0) {
                    found = true;
                    minsd = std::min(minsd,sqrt(-c[j]/b[j]));
                }
            }
        }
        if(found) {
            AddSwitchPoint(i,SP_SINGULAR,minsd);
        }
        aprev = a;
    }
}


////////////////////////////////////////////////////////////////////
//////////////////////////// Profile ///////////////////////////////
////////////////////////////////////////////////////////////////////


Profile::Profile(std::list<dReal>& slist, std::list<dReal>& sdlist, std::list<dReal>&  sddlist, dReal integrationtimestep0) {
    svect = std::vector<dReal>(slist.begin(), slist.end());
    sdvect = std::vector<dReal>(sdlist.begin(), sdlist.end());
    sddvect = std::vector<dReal>(sddlist.begin(), sddlist.end());
    // TODO: handle the case of last step with variable width
    integrationtimestep = integrationtimestep0;
    nsteps = svect.size();
    duration = integrationtimestep * (nsteps-1);
}


bool Profile::FindTimestepIndex(dReal t, int &index, dReal& remainder){
    if (t < -TINY || t > duration+TINY)
        return false;
    else if (t < 0)
        t = 0;
    else if (t > duration)
        t = duration;

    if (duration - t <= TINY)
        index = nsteps - 1;
    else
        index = int(t / integrationtimestep);

    remainder = t - index * integrationtimestep;
    return true;
}


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
        if(!SolveQuadraticEquation(svect[currentindex]-s,sdvect[currentindex],0.5*sddvect[currentindex],tres,0,integrationtimestep)) {
            //std::cout << "***************** Inversion warning: tres=" << tres << " while integrationtimestep= "<< integrationtimestep << "****************\n";
            SolveQuadraticEquation(svect[currentindex]-s,sdvect[currentindex],0.5*sddvect[currentindex],tres,0,integrationtimestep);
        }
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
        if(!SolveQuadraticEquation(svect[currentindex-1]-s,sdvect[currentindex-1],0.5*sddvect[currentindex-1],tres,0,integrationtimestep)) {
            //std::cout << "***************** Inversion warning: tres=" << tres << " while integrationtimestep= "<< integrationtimestep << "****************\n";
            SolveQuadraticEquation(svect[currentindex-1]-s,sdvect[currentindex-1],0.5*sddvect[currentindex-1],tres,0,integrationtimestep);
        }
        t = (currentindex-1)*integrationtimestep + tres;
        return true;
    }
}


dReal Profile::Eval(dReal t){
    int index;
    dReal remainder;
    if(FindTimestepIndex(t, index, remainder))
        return svect[index] + remainder*sdvect[index]
               + 0.5*remainder*remainder*sddvect[index];
    return INF;
}


dReal Profile::Evald(dReal t) {
    int index;
    dReal remainder;
    if (FindTimestepIndex(t, index, remainder))
        return sdvect[index] + remainder * sddvect[index];
    return INF;
}


dReal Profile::Evaldd(dReal t) {
    int index;
    dReal remainder;
    if (FindTimestepIndex(t, index, remainder))
        return sddvect[index];
    return INF;
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
////////////////// Address switch points ///////////////////////////
////////////////////////////////////////////////////////////////////

// Determine whether one can go more than passswitchpointnsteps steps away from (s,sd)
bool PassSwitchPoint(Constraints& constraints, dReal s, dReal sd, dReal dt){
    int ret;
    Profile resprofile;
    bool testaboveexistingprofiles = false, testmvc = true, zlajpah = false;
    ret = IntegrateBackward(constraints,s,sd,dt,resprofile,constraints.passswitchpointnsteps);
    if(ret==INT_MAXSTEPS||ret==INT_END) {
        ret = IntegrateForward(constraints,s,sd,dt,resprofile,constraints.passswitchpointnsteps,testaboveexistingprofiles, testmvc, zlajpah);
        if(ret==INT_MAXSTEPS||ret==INT_END) {
            return true;
        }
    }
    return false;
}

// Bisection search to find the highest sd that allows going through a tangent or discontinuous switch point
dReal BisectionSearch(Constraints& constraints, dReal s, dReal sdbottom, dReal sdtop, dReal dt, int position){
    if(position!=1 && PassSwitchPoint(constraints,s,sdtop,dt)) {
        return sdtop;
    }
    if(sdtop-sdbottom<constraints.bisectionprecision) {
        if(position!=-1 && PassSwitchPoint(constraints,s,sdbottom,dt)) {
            return sdbottom;
        }
        return -1;
    }
    dReal sdmid = (sdbottom+sdtop)*0.5;
    return std::max(BisectionSearch(constraints,s,sdbottom,sdmid,dt,-1),BisectionSearch(constraints,s,sdmid,sdtop,dt,1));
}

// Return false if cannot integrate backward or forward from the switchpoint. This does not mean that the algorithm must fail.
// Return true otherwise. In this case, (sbackward,sdbackward) is the point where the backward integration should start and (sforward,sdforward) is the point where the forward integration should start.
bool AddressSwitchPoint(Constraints& constraints, const SwitchPoint &switchpoint, dReal& sbackward, dReal& sdbackward, dReal& sforward, dReal& sdforward){
    dReal s = switchpoint.s;
    dReal sd = switchpoint.sd;
    dReal dt = constraints.integrationtimestep;
    bool testaboveexistingprofiles = false, testmvc = true, zlajpah = false;
    Profile resprofile;

    // Singular switchpoints
    if(switchpoint.switchpointtype == SP_SINGULAR) {
        dReal bestsstep = 0.1;
        dReal bestslope = 0;
        dReal bestscore = INF;
        dReal sstep = 1e-3/3.1;
        while(sstep <= 0.01) {
            sstep = sstep * 3.1;
            for(int i = 0; i<int(switchpoint.slopesvector.size()); i++) {
                sforward = std::min(s + sstep,constraints.trajectory.duration);
                sbackward = std::max(s - sstep,0.);
                dReal slope = switchpoint.slopesvector[i];
                //if(std::abs(slope*sstep)>0.1) {
                //    continue;
                //}
                sdforward = sd + (sforward-s)*slope;
                sdbackward = sd - (s-sbackward)*slope;
                bool canintegrate = false;
                int ret = IntegrateBackward(constraints,sbackward,sdbackward,dt,resprofile,constraints.passswitchpointnsteps);
                if(ret==INT_MAXSTEPS||ret==INT_END) {
                    ret = IntegrateForward(constraints,sforward,sdforward,dt,resprofile,constraints.passswitchpointnsteps,testaboveexistingprofiles, testmvc, zlajpah);
                    if(ret==INT_MAXSTEPS||ret==INT_END) {
                        canintegrate = true;
                    }
                }
                if(!canintegrate) {
                    continue;
                }
                dReal alphabackward = constraints.SddLimits(sbackward,sdbackward).first/sd;
                dReal betaforward = constraints.SddLimits(sforward,sdforward).second/sd;
                dReal slopediff1 = std::abs(alphabackward-slope);
                dReal slopediff2 = std::abs(betaforward-slope);
                if(sdbackward > 0.01 && sdforward > 0.01) {
                    //dReal score = (slopediff1+slopediff2)/std::abs(std::log(sstep));
                    dReal score = slopediff1+slopediff2;
                    if(score<bestscore) {
                        bestsstep = sstep;
                        bestslope = slope;
                        bestscore = score;
                    }
                }
            }
        }
        sforward = std::min(s + bestsstep,constraints.trajectory.duration);
        sbackward = std::max(s - bestsstep,0.);
        sdforward = sd + (sforward-s)*bestslope;
        sdbackward = sd - (s-sbackward)*bestslope;
        return bestscore<INF;
    }
    // Tangent, Discontinuous and Zlajpah switchpoints
    else{
        dReal sdtop = BisectionSearch(constraints,s,sd*constraints.loweringcoef,sd,constraints.integrationtimestep,0);
        if(sdtop<=0) {
            return false;
        }
        sbackward = s;
        sforward = s;
        sdbackward = sdtop;
        sdforward = sdtop;
        return true;
    }
    return false;
}


////////////////////////////////////////////////////////////////////
////////////// Utilities for integration ///////////////////////////
////////////////////////////////////////////////////////////////////


// Compute the sdd that allows sliding along the curve
dReal ComputeSlidesdd(Constraints& constraints, dReal s, dReal sd, dReal dt) {
    dReal sddtest, snext, sdnext_int, sdnext_mvc, dtsq;
    std ::pair<dReal,dReal> sddlimits = constraints.SddLimits(s,sd);
    dReal alpha = sddlimits.first;
    dReal beta = sddlimits.second;
    dtsq = dt*dt;

    //Check alpha
    snext = s + dt*sd + 0.5*dtsq*alpha;
    sdnext_int = sd + dt*alpha;
    if(snext > constraints.trajectory.duration || snext<0) {
        //std::cout << "Compute slide fin traj\n";
        return 0;
    }
    sdnext_mvc = constraints.SdLimitCombined(snext);
    if(sdnext_mvc < sdnext_int) {
        //std::cout << "Cannot slide alpha \n";
        return 0;
    }

    //Check beta
    snext = s + dt*sd + 0.5*dtsq*beta;
    sdnext_int = sd + dt*beta;
    if(snext > constraints.trajectory.duration || snext<0) {
        //std::cout << "Compute slide fin traj\n";
        return 0;
    }
    sdnext_mvc = constraints.SdLimitCombined(snext);
    if(sdnext_mvc > sdnext_int) {
        //std::cout << "Cannot slide beta \n";
        return 0;
    }

    //Determine the optimal acceleration by bisection
    while(beta-alpha>constraints.bisectionprecision) {
        sddtest = (beta+alpha)/2;
        snext = s + dt*sd + 0.5*dtsq*sddtest;
        sdnext_int = sd + dt*sddtest;
        sdnext_mvc = constraints.SdLimitCombined(snext);
        if(sdnext_int > sdnext_mvc) {
            beta = sddtest;
        }
        else{
            alpha = sddtest;
        }
    }
    return alpha;
}


// Compute the sdd that allows sliding along the curve
dReal ComputeSlidesddBackward(Constraints& constraints, dReal s, dReal sd, dReal dt){
    dReal sddtest, sprev, sdprev_int, sdprev_mvc, dtsq;
    std ::pair<dReal,dReal> sddlimits = constraints.SddLimits(s,sd);
    dReal alpha = sddlimits.first;
    dReal beta = sddlimits.second;
    dtsq = dt*dt;

    //Check alpha
    sprev = s - dt*sd + 0.5*dtsq*alpha;
    sdprev_int = sd - dt*alpha;
    if(sprev > constraints.trajectory.duration || sprev<0) {
        //std::cout << "Compute slide fin traj\n";
        return 0;
    }
    sdprev_mvc = constraints.SdLimitCombined(sprev);
    if(sdprev_mvc > sdprev_int) {
        //std::cout << "Cannot slide alpha \n";
        return 0;
    }

    //Check beta
    sprev = s - dt*sd + 0.5*dtsq*beta;
    sdprev_int = sd - dt*beta;
    if(sprev > constraints.trajectory.duration || sprev<0) {
        //std::cout << "Compute slide fin traj\n";
        return 0;
    }
    sdprev_mvc = constraints.SdLimitCombined(sprev);
    if(sdprev_mvc < sdprev_int) {
        //std::cout << "Cannot slide beta \n";
        return 0;
    }

    //Determine the optimal acceleration by bisection
    while(beta-alpha>constraints.bisectionprecision) {
        sddtest = (beta+alpha)/2;
        sprev = s - dt*sd + 0.5*dtsq*sddtest;
        sdprev_int = sd - dt*sddtest;
        sdprev_mvc = constraints.SdLimitCombined(sprev);
        if(sdprev_int < sdprev_mvc) {
            beta = sddtest;
        }
        else{
            alpha = sddtest;
        }
    }
    return beta;
}


// Determine the relative positions of the flow and the slope of the MVC
int FlowVsMVC(Constraints& constraints, dReal s, dReal sd, int flag, dReal dt) {
    //flag = 1 : alpha flow
    //flag = 2 : beta flow
    //return = -1 if flow points below MVC
    //return = +1 if flow points above MVC
    //return = 0 if end of trajectory
    dReal flow, snext, sdnextflow, sdnextmvc;
    if(s > constraints.trajectory.duration || s<0) {
        return 0;
    }
    std ::pair<dReal,dReal> sddlimits = constraints.SddLimits(s,sd);
    if(flag==1) {
        flow = sddlimits.first;
    }
    else{
        flow = sddlimits.second;
    }
    snext = s + dt * sd + 0.5*dt*dt*flow;
    sdnextflow = sd + dt * flow;
    if(snext > constraints.trajectory.duration || snext<0) {
        return 0;
    }
    sdnextmvc = constraints.SdLimitCombined(snext);
    if(sdnextflow <= sdnextmvc) {
        return -1;
    }
    else{
        return 1;
    }
}

// Determine the relative positions of the alpha flow and the slope of the MVC
int FlowVsMVCBackward(Constraints& constraints, dReal s, dReal sd, dReal dt){
    //return = -1 if alpha points below MVC
    //return = +1 if alpha points above MVC
    //return = 0 if end of trajectory
    dReal alpha, sprev, sdprevalpha, sdprevmvc;
    if(s > constraints.trajectory.duration || s<0) {
        return 0;
    }
    std ::pair<dReal,dReal> sddlimits = constraints.SddLimits(s,sd);
    alpha = sddlimits.first;
    sprev = s - dt * sd + 0.5*dt*dt*alpha;
    sdprevalpha = sd - dt * alpha;
    if(sprev > constraints.trajectory.duration || sprev<0) {
        return 0;
    }
    sdprevmvc = constraints.SdLimitCombined(sprev);
    if(sdprevalpha <= sdprevmvc) {
        return -1;
    }
    else{
        return 1;
    }
}


////////////////////////////////////////////////////////////////////
////////////////////// Integration /////////////////////////////////
////////////////////////////////////////////////////////////////////


int IntegrateForward(Constraints& constraints, dReal sstart, dReal sdstart, dReal dt,  Profile& resprofile, int maxsteps, bool testaboveexistingprofiles, bool testmvc, bool zlajpah){
    dReal dtsq = dt*dt;
    dReal scur = sstart, sdcur = sdstart, snext, sdnext;
    std::list<dReal> slist, sdlist, sddlist;
    bool cont = true;
    int returntype = -1; // should be changed

    // Initialize the currentindex of the profiles for search purpose
    if(testaboveexistingprofiles && constraints.resprofileslist.size()>0) {
        std::list<Profile>::iterator it = constraints.resprofileslist.begin();
        while(it != constraints.resprofileslist.end()) {
            it->currentindex = 0;
            it++;
        }
    }

    // Integrate forward
    while(cont) {
        if(int(slist.size()) > maxsteps) {
            returntype = INT_MAXSTEPS;
            break;
        }
        else if(scur > constraints.trajectory.duration) {
            //TODO: change the time step of previous step to reach the end
            slist.push_back(scur);
            sdlist.push_back(sdcur);
            sddlist.push_back(0);
            returntype = INT_END;
            break;
        }
        else if(sdcur < 0) {
            //TODO: double check whether alpha>beta
            slist.push_back(scur);
            sdlist.push_back(sdcur);
            sddlist.push_back(0);
            returntype = INT_BOTTOM;
            break;
        }
        else if(IsAboveProfilesList(scur,sdcur,constraints.resprofileslist)) {
            slist.push_back(scur);
            sdlist.push_back(sdcur);
            sddlist.push_back(0);
            returntype = INT_PROFILE;
            break;
        }
        else if(zlajpah && testmvc && sdcur >= constraints.SdLimitCombined(scur)-TINY2) {
            if(constraints.SdLimitBobrow(scur)-constraints.SdLimitCombined(scur)<TINY2) {
                slist.push_back(scur);
                sdlist.push_back(sdcur);
                sddlist.push_back(0);
                returntype = INT_MVC;
                break;
            }

            // Lower the sd to the MVC and integrate backward
            if(slist.size()==0) {
                sdcur = constraints.SdLimitCombined(scur);
                Profile tmpprofile;
                //std::cout << "Integrate backward from (" <<scur << "," << sdcur  <<  ")\n";
                int res3 = IntegrateBackward(constraints,scur,sdcur,constraints.integrationtimestep,tmpprofile,1e5,true,true,true);
                if(res3 == INT_BOTTOM) {
                    //std::cout << "BW reached 0 (From Zlajpah)\n";
                    return INT_BOTTOM;
                }
                //std::cout << "BW size " << tmpprofile.nsteps << "\n";
                if(tmpprofile.nsteps>1) {
                    // Add the backward profile to resprofilelist
                    constraints.resprofileslist.push_back(tmpprofile);
                }
            }

            //Now we have sdcombined < sdcur <= sdbobrow
            //We know for sure that beta points above the MVCCombined because we must have reached the MVCCombined following beta
            //There are 2 cases:
            //a) alpha points above MVCCombined (trap case); then step along MVCCombined until alpha points below MVCCombined and integrate backward from that point. Cf. Zlajpah ICRA 1996.
            //b) alpha points below MVCCombined (slide case); then slide along MVCCombined, which is admissible, until either
            // b1) alpha points above MVCCombined (trapped)
            // b2) beta points below MVCCombined (exit slide)

            int res = FlowVsMVC(constraints,scur,sdcur,1,dt);
            if(res == 0) {
                // Most probably we arrived at the end
                slist.push_back(scur);
                sdlist.push_back(sdcur);
                sddlist.push_back(0);
                returntype = INT_END;
                break;
            }
            else if(res == 1) {
                // Case a
                //std::cout <<"\nZlajpah trap ("<< scur << "," << sdcur << ") \n";
                // Insert the profile calculated so far into the resprofileslist
                // And reinitialize the profile
                Profile profile1 = Profile(slist,sdlist,sddlist,dt);
                profile1.forward = true;
                if(profile1.nsteps>1) {
                    constraints.resprofileslist.push_back(profile1);
                }
                slist.resize(0);
                sdlist.resize(0);
                sddlist.resize(0);
                // Now step along the MVCCombined
                while(true) {
                    snext = scur + dt;
                    if(snext>constraints.trajectory.duration) {
                        cont = false;
                        //std::cout << "End traj\n";
                        slist.push_back(scur);
                        sdlist.push_back(sdcur);
                        sddlist.push_back(0);
                        returntype = INT_END; // ZLAJPAH END
                        break;
                    }
                    sdnext = constraints.SdLimitCombined(snext);
                    //std::cout <<"Next ("<< snext << "," << sdnext << ") \n";
                    int res2 = FlowVsMVC(constraints,snext,sdnext,1,dt);
                    if(res2 == 0) {
                        cont = false;
                        //std::cout << "End traj\n";
                        slist.push_back(scur);
                        sdlist.push_back(sdcur);
                        sddlist.push_back(0);
                        returntype = INT_END; // ZLAJPAH END
                        break;
                    }
                    else if(res2 == -1) {
                        // Alpha points below the MVC
                        //std::cout << "End Zlajpah trap (" <<snext << "," << sdnext  <<  ")\n";
                        scur = snext;
                        sdcur = sdnext;
                        break;
                    }
                    scur = snext;
                    sdcur = sdnext;
                }
            }
            else {
                // Case b

                //std::cout <<"\nSlide ("<< scur << "," << sdcur << ") \n";
                while(true) {
                    if(int(slist.size()) > maxsteps) {
                        cont = false;
                        returntype = INT_MAXSTEPS;
                        break;
                    }
                    else if(scur > constraints.trajectory.duration) {
                        cont = false;
                        slist.push_back(scur);
                        sdlist.push_back(sdcur);
                        sddlist.push_back(0);
                        returntype = INT_END;
                        break;
                    }
                    else if(sdcur < 0) {
                        cont = false;
                        slist.push_back(scur);
                        sdlist.push_back(sdcur);
                        sddlist.push_back(0);
                        returntype = INT_BOTTOM;
                        break;
                    }
                    else if(IsAboveProfilesList(scur,sdcur,constraints.resprofileslist)) {
                        cont = false;
                        slist.push_back(scur);
                        sdlist.push_back(sdcur);
                        sddlist.push_back(0);
                        returntype = INT_PROFILE;
                        break;
                    }
                    dReal slidesdd = ComputeSlidesdd(constraints,scur,sdcur,dt);
                    snext = scur + dt * sdcur + 0.5*dtsq*slidesdd;
                    sdnext = sdcur + dt * slidesdd;
                    slist.push_back(scur);
                    sdlist.push_back(sdcur);
                    sddlist.push_back(slidesdd);
                    scur = snext;
                    sdcur = sdnext;
                    //std::cout <<"Next ("<< snext << "," << sdnext << ") \n";
                    int res1 = FlowVsMVC(constraints,snext,sdnext,1,dt);
                    if(res1 == 0) {
                        cont = false;
                        //std::cout << "End traj\n";
                        break;
                    }
                    else if(res1 == 1) {
                        // Case b1
                        //std::cout << "End slide with trap (" <<snext << "," << sdnext  <<  ")\n";
                        break;
                    }
                    int res2 = FlowVsMVC(constraints,snext,sdnext,2,dt);
                    if(res2 == 0) {
                        cont = false;
                        //std::cout << "End traj\n";
                        break;
                    }
                    else if(res2 == -1) {
                        // Case b2
                        //std::cout << "End slide with exit (" <<snext << "," << sdnext  <<  ")\n";
                        break;
                    }
                }
            }
        }
        else if((!zlajpah) && testmvc && sdcur > constraints.SdLimitCombined(scur)+TINY2) {
            slist.push_back(scur);
            sdlist.push_back(sdcur);
            constraints.zlajpahlist.push_back(std::pair<dReal,dReal>(scur,sdcur));
            returntype = INT_MVC;
            break;
        }
        else{
            slist.push_back(scur);
            sdlist.push_back(sdcur);
            dReal beta = constraints.SddLimitBeta(scur,sdcur);
            sddlist.push_back(beta);
            dReal snext = scur + dt * sdcur + 0.5*dtsq*beta;
            dReal sdnext = sdcur + dt * beta;
            scur = snext;
            sdcur = sdnext;
        }
    }
    resprofile = Profile(slist,sdlist,sddlist,dt);
    resprofile.forward = true;
    return returntype;
}


int IntegrateBackward(Constraints& constraints, dReal sstart, dReal sdstart, dReal dt, Profile& resprofile,  int maxsteps, bool testaboveexistingprofiles, bool testmvc,bool zlajpah){
    dReal dtsq = dt*dt;
    dReal scur = sstart, sdcur = sdstart;
    std::list<dReal> slist, sdlist, sddlist;
    bool cont = true;
    int returntype = INT_END;
    bool searchbackward = true;
    dReal alphabk = INF;

    // Initialize the currentindex of the profiles for search purpose
    if(testaboveexistingprofiles && constraints.resprofileslist.size()>0) {
        std::list<Profile>::iterator it = constraints.resprofileslist.begin();
        while(it != constraints.resprofileslist.end()) {
            it->currentindex = it->nsteps-1;
            it++;
        }
    }

    // Integrate backward
    while(cont) {
        if(int(slist.size()) > maxsteps) {
            returntype = INT_MAXSTEPS;
            break;
        }
        else if(scur < 0) {
            //TODO: change the time step of previous step to reach the end
            slist.push_back(scur);
            sdlist.push_back(sdcur);
            sddlist.push_back(0);
            returntype = INT_END;
            break;
        }
        else if(sdcur < 0) {
            //TODO: double check whether alpha>beta
            slist.push_back(scur);
            sdlist.push_back(sdcur);
            sddlist.push_back(0);
            returntype = INT_BOTTOM;
            break;
        }
        else if(IsAboveProfilesList(scur,sdcur,constraints.resprofileslist,searchbackward)) {
            slist.push_back(scur);
            sdlist.push_back(sdcur);
            sddlist.push_back(0);
            returntype = INT_PROFILE;
            break;
        }
        else if(zlajpah && testmvc && sdcur >= constraints.SdLimitCombined(scur)-TINY2 && FlowVsMVCBackward(constraints,scur,sdcur,dt) != -1) {
            if(sdcur > constraints.SdLimitBobrow(scur)) {
                slist.push_back(scur);
                sdlist.push_back(sdcur);
                sddlist.push_back(0);
                returntype = INT_MVC;
                break;
            }
            while(true) {
                if(scur <0) {
                    cont = false;
                    slist.push_back(scur);
                    sdlist.push_back(sdcur);
                    sddlist.push_back(0);
                    returntype = INT_END;
                    break;
                }
                else if(sdcur < 0) {
                    cont = false;
                    slist.push_back(scur);
                    sdlist.push_back(sdcur);
                    sddlist.push_back(0);
                    returntype = INT_BOTTOM;
                    break;
                }
                if(sdcur > constraints.SdLimitBobrow(scur)) {
                    slist.push_back(scur);
                    sdlist.push_back(sdcur);
                    sddlist.push_back(0);
                    returntype = INT_MVC;
                    break;
                }
                else if(IsAboveProfilesList(scur,sdcur,constraints.resprofileslist)) {
                    cont = false;
                    slist.push_back(scur);
                    sdlist.push_back(sdcur);
                    sddlist.push_back(0);
                    returntype = INT_PROFILE;
                    break;
                }
                //std::cout <<"Slide from ("<< scur << "," << sdcur << ") \n";
                dReal slidesdd = ComputeSlidesddBackward(constraints,scur,sdcur,dt);
                dReal sprev = scur - dt * sdcur + 0.5*dtsq*slidesdd;
                dReal sdprev = sdcur - dt * slidesdd;
                slist.push_back(scur);
                sdlist.push_back(sdcur);
                sddlist.push_back(slidesdd);
                scur = sprev;
                sdcur = sdprev;
                int res1 = FlowVsMVCBackward(constraints,sprev,sdprev,dt);
                if(res1 == 0) {
                    cont = false;
                    //std::cout << "End traj\n";
                    returntype = INT_END;
                    break;
                }
                else if(res1 == -1) {
                    // Case b1
                    //std::cout << "End slide with exit (" <<sprev << "," << sdprev  <<  ")\n";
                    break;
                }
            }
        }
        else if(!zlajpah && testmvc && sdcur > constraints.SdLimitCombined(scur)+TINY2) {
            slist.push_back(scur);
            sdlist.push_back(sdcur);
            sddlist.push_back(0);
            returntype = INT_MVC;
            break;
        }
        else{
            slist.push_back(scur);
            sdlist.push_back(sdcur);
            //std ::pair<dReal,dReal> sddlimits = constraints.SddLimits(scur,sdcur);
            //dReal alpha = sddlimits.first;
            dReal alpha = constraints.SddLimitAlpha(scur,sdcur);
            if(sddlist.size()==0) {
                alphabk = alpha;
            }
            sddlist.push_back(alpha);
            //std::cout << scur << " " << sdcur << " " << alpha << "\n";
            dReal sprev = scur - dt * sdcur + 0.5*dtsq*alpha;
            dReal sdprev = sdcur - dt * alpha;
            scur = sprev;
            sdcur = sdprev;
            //alpha = constraints.SddLimitAlpha(std::min(constraints.trajectory.duration,std::max(0.,scur)),std::max(0.,sdcur));
        }
    }
    if(sddlist.size()>0) {
        slist.reverse();
        sdlist.reverse();
        sddlist.reverse();
        sddlist.pop_front();
        if(alphabk == INF) {
            sddlist.push_back(sddlist.back());
        }
        else{
            sddlist.push_back(alphabk);
        }
    }
    resprofile = Profile(slist,sdlist,sddlist,dt);
    resprofile.forward = false;
    return returntype;
}


////////////////////////////////////////////////////////////////////
//////////////////////// Limiting curves ///////////////////////////
////////////////////////////////////////////////////////////////////

Profile StraightProfile(dReal sbackward,dReal sforward,dReal sdbackward,dReal sdforward){
    std::list<dReal> slist, sdlist, sddlist;
    dReal dtmod = 2 * (sforward - sbackward) / (sdforward + sdbackward);
    dReal sdd = (sdforward - sdbackward) / dtmod;
    slist.push_back(sbackward);
    slist.push_back(sforward);
    sdlist.push_back(sdbackward);
    sdlist.push_back(sdforward);
    sddlist.push_back(sdd);
    sddlist.push_back(0);
    return Profile(slist,sdlist,sddlist,dtmod);
}


int ComputeLimitingCurves(Constraints& constraints){
    std::list<SwitchPoint> switchpointslist0(constraints.switchpointslist);

    Profile tmpprofile;
    dReal sswitch, sdswitch, sbackward, sdbackward, sforward, sdforward;
    int integratestatus;
    bool testaboveexistingprofiles = true, testmvc = true, zlajpah = false;

    while(switchpointslist0.size() > 0) {
        SwitchPoint switchpoint = switchpointslist0.front();
        switchpointslist0.pop_front();
        sswitch = switchpoint.s;
        sdswitch = switchpoint.sd;
        if(IsAboveProfilesList(sswitch,sdswitch,constraints.resprofileslist,false,0.1*sdswitch))
            continue;
        if(sdswitch > constraints.SdLimitCombined(sswitch)+TINY2)
            continue;

        //std::cout << sswitch << "\n";

        // Address Switch Point
        if (!AddressSwitchPoint(constraints, switchpoint, sbackward,
                                sdbackward, sforward, sdforward))
            continue;

        // Add middle part
        if (sforward - sbackward > TINY) {
            constraints.resprofileslist.push_back(StraightProfile(sbackward,sforward,sdbackward,sdforward));
            //std::cout << "Main stub : " << constraints.resprofileslist.back().duration << "\n";

        }

        // Integrate backward
        integratestatus = IntegrateBackward(constraints, sbackward, sdbackward,
                                            constraints.integrationtimestep, tmpprofile);
        if(tmpprofile.nsteps>2) {
            constraints.resprofileslist.push_back(tmpprofile);
            //std::cout << "Backward : "  << tmpprofile.nsteps << " " << tmpprofile.duration << "\n";
        }

        if(integratestatus == INT_BOTTOM)
            return CLC_BOTTOM;

        // Integrate forward
        integratestatus = IntegrateForward(constraints, sforward, sdforward,
                                           constraints.integrationtimestep, tmpprofile, 1e5,
                                           testaboveexistingprofiles, testmvc, zlajpah);
        if(tmpprofile.nsteps>2) {
            constraints.resprofileslist.push_back(tmpprofile);
            //std::cout << "Forward : " << tmpprofile.nsteps << " " << tmpprofile.duration << "\n";
        }

        if(integratestatus == INT_BOTTOM)
            return CLC_BOTTOM;
    }
    return CLC_OK;
}


int ComputeProfiles(Constraints& constraints, dReal sdbeg, dReal sdend){
    std::chrono::time_point<std::chrono::system_clock> t0,t1,t2,t3;
    std::chrono::duration<double> d1,d2,d3,dtot;

    t0 = std::chrono::system_clock::now();

    bool retprocess = constraints.Preprocess();
    if(!retprocess) {
        std::cout << "Trajectory duration too short\n";
        return TOPP_SHORT_TRAJ;
    }

    if(VectorMin(constraints.mvcbobrow) <= TINY) {
        //std::cout << "[TOPP] MVCBobrow hit 0\n";
        return TOPP_MVC_HIT_ZERO;
    }

    Profile resprofile;
    int ret;

    t1 = std::chrono::system_clock::now();

    bool integrateprofilesstatus = true;

    std::string message;
    for(int rep=0; rep<constraints.extrareps+1; rep++) {

        // Lower integrationtimestep
        if(rep>0) {
            constraints.integrationtimestep /= 2;
            //constraints.passswitchpointnsteps *= 2;
            std::cout << rep << "!!!!!!!!!!!!!!!!!!!!!!!!!! Try lower integration timestep: " << constraints.integrationtimestep << "!!!!!!!!!!!!!!!!!!!!!!!!\n";

        }
        constraints.resprofileslist.resize(0);
        constraints.zlajpahlist.resize(0);

        // flags
        bool testaboveexistingprofiles = true, testmvc = true, zlajpah = false;


        /////////////////  Compute the CLC /////////////////////
        ret = ComputeLimitingCurves(constraints);
        if(ret!=CLC_OK) {
            message = "CLC failed\n";
            integrateprofilesstatus = false;
            continue;
        }


        Profile tmpprofile;
        dReal smallincrement = constraints.integrationtimestep*2;
        dReal bound,tres;


        /////////////////  Integrate from start /////////////////////
        // Fix start
        dReal sstartnew = 0, sdstartnew = sdbeg;
        constraints.FixStart(sstartnew,sdstartnew);
        if(sstartnew<=TINY2 || sdstartnew > constraints.SdLimitCombined(0) - TINY2) {
            sdstartnew = sdbeg;
        }
        // Check whether sdend > CLC or MVC
        if(FindLowestProfile(smallincrement,tmpprofile,tres,constraints.resprofileslist)) {
            bound = std::min(tmpprofile.Evald(tres),constraints.mvccombined[0]);
        }
        else{
            bound = constraints.mvccombined[0];
        }
        if(sdstartnew > bound) {
            message = "Sdbeg > CLC or MVC\n";
            integrateprofilesstatus = false;
            continue;
        }
        // Integrate from s = 0
        ret = IntegrateForward(constraints,sstartnew,sdstartnew,constraints.integrationtimestep,resprofile,1e5,testaboveexistingprofiles,testmvc,zlajpah);
        if(resprofile.nsteps>1) {
            constraints.resprofileslist.push_back(resprofile);
        }
        if(ret==INT_BOTTOM) {
            message = "FW reached 0";
            integrateprofilesstatus = false;
            continue;
        }


        /////////////////  Integrate from end /////////////////////
        // Fix end
        dReal sendnew = constraints.trajectory.duration, sdendnew = sdend;
        constraints.FixEnd(sendnew,sdendnew);
        if(constraints.trajectory.duration-sendnew<=TINY2 || sdendnew > constraints.SdLimitCombined(sendnew)- TINY2) {
            sdendnew = sdend;
        }
        // Check whether sdend > CLC or MVC
        if(FindLowestProfile(constraints.trajectory.duration-smallincrement,tmpprofile,tres,constraints.resprofileslist)) {
            bound = std::min(tmpprofile.Evald(tres),constraints.mvccombined[constraints.mvccombined.size()-1]);
        }
        else{
            bound = constraints.mvccombined[constraints.mvccombined.size()-1];
        }
        if(sdendnew > bound) {
            message = "Sdend > CLC or MVC\n";
            integrateprofilesstatus = false;
            continue;
        }
        // Integrate back from s = send
        ret = IntegrateBackward(constraints,sendnew,sdendnew,constraints.integrationtimestep,resprofile,1e5,testaboveexistingprofiles,testmvc);
        if(resprofile.nsteps>1) {
            constraints.resprofileslist.push_back(resprofile);
        }
        if(ret==INT_BOTTOM) {
            message = "BW reached 0";
            integrateprofilesstatus = false;
            continue;
        }


        /////////////////////// Zlajpah //////////////////////////
        // Add separation points between MVCBobrow and MVCCombined to zlajpahlist
        bool active = false;
        for(int i = 0; i<constraints.ndiscrsteps; i++) {
            if(std::abs(constraints.mvcbobrow[i]-constraints.mvccombined[i])<TINY2) {
                if(!active) {
                    active = true;
                }
            }
            else{
                if(active) {
                    active = false;
                    constraints.zlajpahlist.push_back(std::pair<dReal,dReal>(constraints.discrsvect[i],constraints.mvccombined[i]));
                }
            }
        }

        // Integrate forward from Zlajpah points
        zlajpah = true;
        std::list<std::pair<dReal,dReal> >::iterator zit = constraints.zlajpahlist.begin();
        bool zlajpaherror = false;
        while(zit!=constraints.zlajpahlist.end()) {
            dReal zs = zit->first;
            dReal zsd = zit->second;
            ret = IntegrateForward(constraints,zs,zsd,constraints.integrationtimestep,resprofile,1e5,testaboveexistingprofiles,testmvc,zlajpah);
            if(resprofile.nsteps>1) {
                constraints.resprofileslist.push_back(resprofile);
            }
            if(ret==INT_BOTTOM) {
                zlajpaherror = true;
                break;
            }
            zit++;
        }
        if(zlajpaherror) {
            message = "Zlajpah error";
            integrateprofilesstatus = false;
            continue;
        }


        /////////////////////  Final checks /////////////////////////
        // Reset currentindex
        std::list<Profile>::iterator it = constraints.resprofileslist.begin();
        while(it != constraints.resprofileslist.end()) {
            it->currentindex = 0;
            it++;
        }

        // Estimate resulting trajectory duration
        constraints.resduration = 0;
        Profile profile;
        //dReal ds = constraints.discrtimestep;
        dReal ds = 1e-2;
        int nsamples = int((constraints.trajectory.duration-TINY)/ds);
        dReal s,sdcur,sdnext;
        if(FindLowestProfile(0,profile,tres,constraints.resprofileslist)) {
            sdcur = profile.Evald(tres);
        }
        else{
            message = "CLC discontinuous at 0";
            integrateprofilesstatus = false;
            continue;
        }
        bool clcdiscontinuous = false;
        for(int i=1; i<=nsamples; i++) {
            s = i*ds;
            if(FindLowestProfile(s,profile,tres,constraints.resprofileslist)) {
                sdnext = profile.Evald(tres);
                constraints.resduration += 2*ds/(sdcur+sdnext);
                sdcur = sdnext;
            }
            else{
                clcdiscontinuous = true;
                break;
                std::cout << "CLC discontinuous: " << s << "\n";
            }
        }
        if(clcdiscontinuous) {
            message = std::string("CLC discontinuous : ") + std::to_string((long double)s);
            integrateprofilesstatus = false;
            continue;
        }

        integrateprofilesstatus = true;
        break;

    }


    if(!integrateprofilesstatus) {
        std::cout << "[TOPP::ComputeProfiles] " << message << "\n";
        return TOPP_UNSPEC;
    }

    t2 = std::chrono::system_clock::now();

    d1 = t1-t0;
    d2 = t2-t1;
    d3 = t3-t2;
    dtot = t3-t0;

    //std::cout << "Constraints preprocessing : " << d1.count() << " \n ";
    //std::cout << "Profiles calculation : " <<  d2.count() << " \n ";
    //std::cout << "Duration calculation : " <<  d3.count() << " \n ";
    //std::cout << "Total PP : " <<  dtot.count() << " \n ";

    return TOPP_OK;
}


int VIP(Constraints& constraints, dReal sdbegmin, dReal sdbegmax, dReal& sdendmin, dReal& sdendmax){
    if (constraints.trajectory.duration <= 0) {
        std::cout << "[TOPP::VIP] Warning : trajectory duration is <= 0 \n";
        return TOPP_SHORT_TRAJ;
    }

    bool retprocess = constraints.Preprocess();
    if(!retprocess) {
        std::cout << "Trajectory duration too short\n";
        return TOPP_SHORT_TRAJ;
    }

    if(VectorMin(constraints.mvcbobrow) <= TINY) {
        std::cout << "[TOPP::VIP] MVCBobrow hit 0 \n";
        return TOPP_MVC_HIT_ZERO;
    }

    Profile tmpprofile;
    dReal tres;
    dReal smallincrement = constraints.integrationtimestep*2;

    // Compute the limiting curves
    int resclc = ComputeLimitingCurves(constraints);
    if(resclc == CLC_SWITCH || resclc == CLC_BOTTOM) {
        std::cout << "[TOPP::VIP] resclc == CLC_SWITCH or CLC_BOTTOM \n";
        return TOPP_CLC_ERROR;
    }

    // Determine the lowest profile at t=0
    dReal bound;
    if(FindLowestProfile(smallincrement,tmpprofile,tres,constraints.resprofileslist))
        bound = std::min(tmpprofile.Evald(tres),constraints.mvccombined[0]);
    else // just to make sure the profile is below mvccombined
        bound = constraints.mvccombined[0];

    if(sdbegmin>bound) {
        std::cout << "[TOPP::VIP] sdbegmin is above the CLC or the combined MVC \n";
        return TOPP_SDBEGMIN_TOO_HIGH;
    }

    sdbegmax = std::min(sdbegmax,bound-constraints.bisectionprecision);

    bool testaboveexistingprofiles = true, testmvc = true, zlajpah = true;

    // Compute sdendmax by integrating forward from (0,sdbegmax)
    int resintfw = IntegrateForward(constraints,0,sdbegmax,constraints.integrationtimestep,tmpprofile,1e5,testaboveexistingprofiles,testmvc,zlajpah);
    constraints.resprofileslist.push_back(tmpprofile);
    if(resintfw == INT_BOTTOM) {
        std::cout << "[TOPP::VIP] Forward integration hit sd = 0 \n";
        return TOPP_FWD_HIT_ZERO;
    }
    else if (resintfw == INT_END && tmpprofile.Evald(tmpprofile.duration) <= constraints.mvccombined[constraints.mvccombined.size() - 1])
        sdendmax = tmpprofile.Evald(tmpprofile.duration);
    else if (resintfw == INT_MVC || resintfw == INT_PROFILE) {
        // Look for the lowest profile at the end
        if(FindLowestProfile(constraints.trajectory.duration-smallincrement,tmpprofile,tres,constraints.resprofileslist) && tmpprofile.Evald(tres) <= constraints.mvccombined[constraints.mvccombined.size()-1])
            sdendmax = tmpprofile.Evald(tres);
        else {
            // No profile reaches the end, consider the MVC instead
            sdendmax = constraints.mvccombined[constraints.mvccombined.size()-1]-smallincrement;
            int count = 0;
            int resintbw;
            dReal dtint = constraints.integrationtimestep;
            // If integrating from sdendmax fails with INT_BOTTOM or INT_MVC,
            // then the trajectory is not traversable. However, since
            // integrating backward from a high sdendmax can be risky, we give
            // three chances by decreasing the value of the integration step
            while(count<3) {
                count++;
                resintbw = IntegrateBackward(constraints,constraints.trajectory.duration,sdendmax,dtint,tmpprofile,1e5);
                if(resintbw != INT_BOTTOM && resintbw != INT_MVC) {
                    break;
                }
                dtint /= 3.3;
            }
            if(resintbw == INT_BOTTOM || resintbw == INT_MVC || (resintbw == INT_END && tmpprofile.Evald(0)<sdbegmin))
                return TOPP_BWD_FAIL;
        }
    }

    // Integrate from (send,0). If succeeds, then sdendmin=0 and exits
    int resintbw = IntegrateBackward(constraints,constraints.trajectory.duration,0,constraints.integrationtimestep,tmpprofile,1e5);
    if((resintbw == INT_END && tmpprofile.Evald(0)>=sdbegmin) || resintbw == INT_PROFILE) {
        constraints.resprofileslist.push_back(tmpprofile);
        sdendmin = 0;
        return TOPP_OK;
    }

    // Determine sdendmin by bisection
    dReal sdupper = sdendmax, sdlower = 0;
    Profile bestprofile;
    while(sdupper-sdlower > constraints.bisectionprecision) {
        dReal sdtest = (sdupper + sdlower)/2;
        int resintbw2 = IntegrateBackward(constraints,constraints.trajectory.duration,sdtest,constraints.integrationtimestep,tmpprofile,1e5);
        if((resintbw2 == INT_END && tmpprofile.Evald(0)>=sdbegmin) || resintbw2 == INT_PROFILE) {
            sdupper = sdtest;
            bestprofile = tmpprofile;
        }
        else
            sdlower = sdtest;
    }
    sdendmin = sdupper;
    constraints.resprofileslist.push_back(bestprofile);

    return TOPP_OK;
}


int VIPBackward(Constraints& constraints, dReal& sdbegmin, dReal& sdbegmax, dReal sdendmin, dReal sdendmax) {
    if (constraints.trajectory.duration <= 0) {
        std::cout << "[TOPP::VIPBackward] Warning : trajectory duration is <= 0 \n";
        return TOPP_SHORT_TRAJ;
    }

    bool retprocess = constraints.Preprocess();
    if(!retprocess) {
        std::cout << "Trajectory duration too short\n";
        return TOPP_SHORT_TRAJ;
    }

    if(VectorMin(constraints.mvcbobrow) <= TINY) {
        std::cout << "[TOPP::VIPBackward] MVCBobrow hit 0 \n";
        return TOPP_MVC_HIT_ZERO;
    }

    Profile tmpprofile;
    dReal tres;
    dReal smallincrement = constraints.integrationtimestep*2;

    // Compute the limiting curves
    int resclc = ComputeLimitingCurves(constraints);
    if (resclc == CLC_SWITCH || resclc == CLC_BOTTOM) {
        std::cout << "[TOPP::VIPBackward] resclc == CLC_SWITCH or CLC_BOTTOM \n";
        return TOPP_CLC_ERROR;
    }

    // Determine the lowest profile at send
    dReal bound;
    if (FindLowestProfile(constraints.trajectory.duration, tmpprofile, tres, constraints.resprofileslist))
        bound = std::min(tmpprofile.Evald(tres), constraints.mvccombined[constraints.mvccombined.size() - 1]);
    else // just to make sure the profile is below mvccombined
        bound = constraints.mvccombined[constraints.mvccombined.size() - 1];

    if (sdendmin > bound) {
        std::cout << "[TOPP::VIPBackward] sdendmin is above the combined MVC \n";
        return TOPP_SDENDMIN_TOO_HIGH;
    }

    sdendmax = std::min(sdendmax, bound-constraints.bisectionprecision);

    bool testaboveexistingprofiles = true, testmvc = true, zlajpah = true;

    // Compute sdbegmax by integrating backward from (send, sdendmax)
    int resintbw = IntegrateBackward(constraints, constraints.trajectory.duration, sdendmax, constraints.integrationtimestep, tmpprofile, 1e5, testaboveexistingprofiles, testmvc, zlajpah);
    constraints.resprofileslist.push_back(tmpprofile);
    if (resintbw == INT_BOTTOM) {
        std::cout << "[TOPP::VIPBackward] Backward integration hits sd = 0 \n";
        return TOPP_BWD_HIT_ZERO;
    }
    // Phi profile reaches s = 0
    else if (resintbw == INT_END && tmpprofile.Evald(0) <= constraints.mvccombined[0]) {
        // Phi profile reaches s = 0 and is below mvccombined
        sdbegmax = tmpprofile.Evald(0);
    }
    // Phi profile hits MVC or CLC
    else if (resintbw == INT_MVC || resintbw == INT_PROFILE) {
        // Look for the lowest profile at s = 0
        if (FindLowestProfile(smallincrement, tmpprofile, tres, constraints.resprofileslist) && tmpprofile.Evald(tres) <= constraints.mvccombined[0])
            sdbegmax = tmpprofile.Evald(tres);
        else {
            // No profile reaches s = 0, consider the MVC instead
            sdbegmax = constraints.mvccombined[0] - smallincrement;
            int count = 0;
            int resintfw;
            dReal dint = constraints.integrationtimestep;
            // If integrating from sdbegmax fails with INT_BOTTOM or INT_MVC, then the trajectory is not traversable. However, since integrating forward from a high sdbegmax can be risky, we give three chances by decreasing the value of the integration step.
            while (count < 3) {
                count++;
                resintfw = IntegrateForward(constraints, 0, sdbegmax, dint, tmpprofile, 1e5);
                if (resintfw != INT_BOTTOM && resintfw != INT_MVC) {
                    break;
                }
                dint /= 3.3;
            }
            if (resintfw == INT_BOTTOM || resintfw == INT_MVC || (resintfw == INT_END && tmpprofile.Evald(tmpprofile.duration) < sdendmin))
                return TOPP_FWD_FAIL;
        }
    }

    // Integrate from (0, 0). If succeeds, then sdbegmin = 0 and exit.
    int resintfw = IntegrateForward(constraints, 0, 0, constraints.integrationtimestep, tmpprofile, 1e5);
    if ((resintfw == INT_END && tmpprofile.Evald(tmpprofile.duration) >= sdendmin) || resintfw == INT_PROFILE) {
        constraints.resprofileslist.push_back(tmpprofile);
        sdbegmin = 0;
        return TOPP_OK;
    }

    // If sdbegmin != 0, determine sdbegmin by bisection
    dReal sdupper = sdbegmax, sdlower = 0;
    Profile bestprofile;
    while (sdupper - sdlower > constraints.bisectionprecision) {
        dReal sdtest = (sdupper + sdlower)/2;
        int resintfw2 = IntegrateForward(constraints, 0, sdtest, constraints.integrationtimestep, tmpprofile, 1e5);
        if ((resintfw2 == INT_END && tmpprofile.Evald(tmpprofile.duration) >= sdendmin) || resintfw2 == INT_PROFILE) {
            sdupper = sdtest;
            bestprofile = tmpprofile;
        }
        else
            sdlower = sdtest;
    }
    sdbegmin = sdupper;
    constraints.resprofileslist.push_back(bestprofile);

    return TOPP_OK;
}


////////////////////////////////////////////////////////////////////
///////////////////////// Utilities ////////////////////////////////
////////////////////////////////////////////////////////////////////


void VectorFromString(const std::string& s,std::vector<dReal>&resvect){
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

void ReadVectorFromStream(std::istream& s, size_t N, std::vector<dReal>& resvect)
{
    resvect.resize(N);
    for(size_t i = 0; i < N; ++i) {
        s >> resvect[i];
    }
    if( !s ) {
        throw TOPP_EXCEPTION_FORMAT("failed to read %d items from stream", N, 0);
    }
}

dReal VectorMin(const std::vector<dReal>&v){
    dReal res = INF;
    for(int i=0; i<int(v.size()); i++) {
        res = std::min(res,v[i]);
    }
    return res;
}


dReal VectorMax(const std::vector<dReal>&v){
    dReal res = -INF;
    for(int i=0; i<int(v.size()); i++) {
        res = std::max(res,v[i]);
    }
    return res;
}


void PrintVector(const std::vector<dReal>& v){
    std::cout << "[";
    for(int i=0; i<int(v.size()); i++) {
        std::cout<< v[i] << ", ";
    }
    std::cout << "] \n ";
}

void VectorAdd(const std::vector<dReal>&a, const std::vector<dReal>&b,  std::vector<dReal>&res, dReal coefa, dReal coefb){
    res.resize(a.size());
    assert(a.size() == b.size());
    for(int i=0; i<int(a.size()); i++) {
        res[i] = coefa*a[i]+coefb*b[i];
    }
}

void VectorMultScalar(const std::vector<dReal>&a, std::vector<dReal>&res, dReal scalar){
    res.resize(a.size());
    for(int i=0; i<int(a.size()); i++) {
        res[i] = scalar*a[i];
    }
}


dReal VectorNorm(const std::vector<dReal>&v){
    dReal norm = 0;
    for(int i=0; i<int(v.size()); i++) {
        norm += v[i]*v[i];
    }
    return sqrt(norm);
}


bool SolveQuadraticEquation(dReal a0, dReal a1, dReal a2, dReal& sol, dReal lowerbound, dReal upperbound) {
    dReal delta = a1*a1- 4*a0*a2;
    if(delta<-TINY2) {
        return false;
    }
    delta = std::max(delta,0.);
    if(a2<0) {
        a2=-a2;
        a1=-a1;
        a0=-a0;
    }
    if(a2<=TINY2) {
        if(std::abs(a1)<=TINY2) {
            return false; // Must be in a weird case
        }
        sol = -a0/a1;
        return sol>=lowerbound-TINY2 && sol<=upperbound+TINY2;
    }
    sol = (-a1-sqrt(delta))/(2*a2);
    if(sol>=lowerbound-TINY2 && sol<=upperbound+TINY2) {
        return true;
    }
    if(sol>upperbound+TINY2) {
        return false;
    }
    sol = (-a1+sqrt(delta))/(2*a2);
    return sol>=lowerbound-TINY2 && sol<=upperbound+TINY2;
}


bool IsAboveProfilesList(dReal s, dReal sd, std::list<Profile>&resprofileslist, bool searchbackward, dReal softborder){
    dReal t;
    std::list<Profile>::iterator it = resprofileslist.begin();
    while(it != resprofileslist.end()) {
        if(searchbackward) {
            it->currentindex = it->nsteps-1;
        }
        else{
            it->currentindex = 0;
        }

        if(it->Invert(s,t,searchbackward)) {
            if(sd > it->Evald(t) + softborder) {
                return true;
            }
        }
        it++;
    }
    return false;
}


bool FindLowestProfile(dReal s, Profile& profile, dReal& tres, std::list<Profile>&resprofileslist){
    dReal t;
    dReal sdmin;
    dReal sdtmp;
    int i = 0;
    sdmin = INF;

    std::list<Profile>::iterator it = resprofileslist.begin();
    while(it != resprofileslist.end()) {
        if(it->Invert(s,t)) {
            sdtmp = it->Evald(t);
            if(sdtmp < sdmin) {
                sdmin = sdtmp;
                profile = *it;
                tres = t;
            }
        }
        it++;
        i++;
    }
    return (sdmin < INF);
}


void CheckInsert(std::list<std::pair<dReal,dReal> >& reslist, std::pair<dReal,dReal> e, bool inverse){
    std::list<std::pair<dReal,dReal> >::iterator it = reslist.begin();
    while(it != reslist.end()) {
        if(inverse) {
            if(e.first >= it->first && e.second >= it->second) {
                return;
            }
            if(e.first <= it->first && e.second <= it->second) {
                it = reslist.erase(it);
            }
            else{
                it++;
            }
        }
        else{
            if(e.first <= it->first && e.second <= it->second) {
                return;
            }
            if(e.first >= it->first && e.second >= it->second) {
                it = reslist.erase(it);
            }
            else{
                it++;
            }
        }
    }
    reslist.push_back(e);
}


void FindMaxima(const std::list<std::pair<dReal,dReal> >& origlist, std::list<std::pair<dReal,dReal> >& reslist, bool inverse){
    reslist.resize(0);
    if(origlist.size()==0) {
        return;
    }
    std::list<std::pair<dReal,dReal> >::const_iterator it = origlist.begin();
    while(it != origlist.end()) {
        CheckInsert(reslist,*it,inverse);
        it++;
    }
}



} // end namespace TOPP
