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
    switchpointslist.clear();
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
        for(size_t i=0; i< mvccombined.size(); i++) {
            meanmvc += std::min(mvccombined[i],10.);
        }
        meanmvc /= mvccombined.size();
        meanmvc = std::min(1.,meanmvc);
        integrationtimestep = discrtimestep/meanmvc;
        //std::cout << "\n--------------\nIntegration timestep: " << integrationtimestep << "\n";
    }

    return true;
}


void Constraints::Discretize() {
    ndiscrsteps = int((trajectory.duration+TINY)/discrtimestep);
    ndiscrsteps++;
    discrsvect.resize(ndiscrsteps);
    for(int i=0; i<ndiscrsteps; i++) {
        discrsvect[i] = i*discrtimestep;
    }
}


void Constraints::ComputeMVCBobrow() {
    mvcbobrow.resize(ndiscrsteps);
    for(int i=0; i<ndiscrsteps; i++) {
        mvcbobrow[i] = SdLimitBobrowInit(discrsvect[i]);
    }
}


void Constraints::ComputeMVCCombined()
{
    mvccombined.resize(ndiscrsteps);
    for(int i=0; i<ndiscrsteps; i++) {
        mvccombined[i] = SdLimitCombinedInit(discrsvect[i]);
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
            if(std::abs(qd[i])>TINY && std::abs(vmax[i])>0) {
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
            if(std::abs(qd[i])>TINY && std::abs(vmax[i])>0) {
                res = std::min(res,vmax[i]/std::abs(qd[i]));
            }
        }
        ss << res << " ";
    }
}


void Constraints::FindSwitchPoints()
{
    switchpointslist.clear();
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
        sw.slopesvector.resize(0);
        ComputeSlopeDynamicSingularity(s,sd,sw.slopesvector);
    }
    switchpointslist.insert(it,sw);
}

bool CompareSwitchPoint(const SwitchPoint& sw0, const SwitchPoint& sw1)
{
    return sw0.s < sw1.s;
}

void Constraints::AddSwitchPoint2(dReal s, dReal sd, int switchpointtype)
{
    // If no sd is specified, then take the value of the mvc
    // (The case when sd is specified corresponds to a singular switchpoint in some cases)
    if(sd > MAXSD) {
        return;
    }
    SwitchPoint sw(s,sd,switchpointtype);
    std::list<SwitchPoint>::iterator it = std::lower_bound(switchpointslist.begin(), switchpointslist.end(), sw, CompareSwitchPoint);
    if( it != switchpointslist.end() ) {
        if( s >= it->s+TINY ) {
            std::cout << "switch point already exists, type=" << it->switchpointtype;
            return;
        }
    }
    if(switchpointtype == SP_SINGULAR) {
        sw.slopesvector.resize(0);
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
            AddSwitchPoint2(s,sd,SP_DISCONTINUOUS);
        }
        prevtangent = tangent;
        tangent = (sdnext-sd)/discrtimestep;
        //if(std::abs(tangent-prevtangent)>1.) {
        //    continue;
        //}
        //beta = sddlimits.second;
        diff = alpha/sd - tangent;
        if(diffprev*diff<0 && std::abs(diff)<1) {
            AddSwitchPoint2(s,sd,SP_TANGENT);
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
    // also look for the start of the chunks for the trajectory
    //std::list<dReal>::const_iterator itchunkstart = trajectory.chunkcumulateddurationslist.begin();
    //int nLastAddedSwitchIndex = -1;
    for(int i=0; i<ndiscrsteps-2; i++) {
        sdnn = SdLimitBobrow(discrsvect[i+2]);
        if(std::abs(sdnn-sdn)>100*std::abs(sdn-sd)) {
            if(sdn<sdnn) {
                AddSwitchPoint2(discrsvect[i+1],mvcbobrow[i+1],SP_DISCONTINUOUS);
            }
            else{
                AddSwitchPoint2(discrsvect[i+2],mvcbobrow[i+2],SP_DISCONTINUOUS);
            }
        }
        // if( trajectory.degree <= 3 ) {
        //     // if the trajectory degree is <= 3, then the accelerations will not be differentiable at the trajectory chunk edges.
        //     // therefore add those discontinuity points.
        //     // perhaps there's a better way to compute this, but the above threshold doesn't catch it.
        //     if( itchunkstart != trajectory.chunkcumulateddurationslist.end() && *itchunkstart <= discrsvect[i+2]+TINY ) {
        //         if( nLastAddedSwitchIndex < i+1 ) {
        //             AddSwitchPoint2(discrsvect[i+1],mvcbobrow[i+1],SP_DISCONTINUOUS);
        //             nLastAddedSwitchIndex = i+1;
        //         }
        //         ++itchunkstart;
        //     }
        // }
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
    // Find the first singular point, if any
    while(it!=switchpointslist.end()) {
        if(it->switchpointtype == SP_SINGULAR) {
            scur = it->s;
            sdcur = it->sd;
            itcur = it;
            break;
        }
        it++;
    }
    // Merge consecutive singular points that are in a small radius
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

    // Merge non-singular switchpoints
    // Find the first non singular switchpoint, if any
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
    // Merge consecutive non-singular switchpoints that are in a small radius
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
                        // don't set scur/sdcur since then radius will be sliding
                        //scur = snext;
                        //sdcur = sdnext;
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


void QuadraticConstraints::WriteConstraints(std::stringstream& ss){
    ss << discrtimestep << "\n";
    for(int i=0; i<int(vmax.size()); i++) {
        ss << vmax[i] << " ";
    }
    ss << "\n";
    for(int i=0; i<int(avect.size()); i++) {
        for(int j=0; j<int(avect[0].size()); j++) {
            ss << avect[i][j] << " ";
        }
        ss << "\n";
        for(int j=0; j<int(avect[0].size()); j++) {
            ss << bvect[i][j] << " ";
        }
        ss << "\n";
        for(int j=0; j<int(avect[0].size()); j++) {
            ss << cvect[i][j] << " ";
        }
        if(i<int(avect.size())-1) {
            ss << "\n";
        }
    }
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
    BOOST_ASSERT(s>=-TINY && s<=trajectory.duration+TINY);
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
    dReal coef = (s-n*discrtimestep);
    if( std::abs(coef) <= TINY ) {
        a = avect[n];
        b = bvect[n];
        c = cvect[n];
    }
    else {
        coef /= discrtimestep;
        for (int i = 0; i < nconstraints; i++) {
            a[i] = (1-coef)*avect[n][i] + coef*avect[n+1][i];
            b[i] = (1-coef)*bvect[n][i] + coef*bvect[n+1][i];
            c[i] = (1-coef)*cvect[n][i] + coef*cvect[n+1][i];
        }
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
    dReal idelta=1/delta;
    slopesvector.resize(a.size());
    for(size_t i=0; i< a.size(); i++) {
        ap = (a2[i]-a[i])*idelta;
        bp = (b2[i]-b[i])*idelta;
        cp = (c2[i]-c[i])*idelta;
        if(std::abs((2*b[i]+ap)*sd)>TINY) {
            slope = (-bp*sd*sd-cp)/((2*b[i]+ap)*sd);
        }
        else{
            slope = 0;
        }
        slopesvector[i] = slope;
    }
}

std::pair<dReal,dReal> QuadraticConstraints::SddLimits(dReal s, dReal sd){
    dReal dtsq = integrationtimestep;
    dtsq = dtsq*dtsq;
    dReal alpha = -INF;
    dReal beta = INF;
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

// Compute the SdLimitBobrow after removing one inequality (sdot^\dag in Pham 2014)
dReal QuadraticConstraints::SdLimitBobrowExclude(dReal s, int iexclude){
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
            if (k==iexclude || m==iexclude) {
                continue;
            }
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
    std::vector<dReal> a,b,c, aprev, bprev, cprev;

    InterpolateDynamics(discrsvect[i],aprev,bprev,cprev);

    for(int i=1; i<ndiscrsteps-1; i++) {
        InterpolateDynamics(discrsvect[i],a,b,c);
        dReal minsd = INF;
        dReal mins = INF;
        bool found = false;
        for(int j=0; j< int(a.size()); j++) {
            if(a[j]*aprev[j]<=0) {
                dReal adiff = a[j] - aprev[j];
                dReal ccur=c[j], bcur=b[j], scur=discrsvect[i];
                if( fabs(adiff) > TINY ) {
                    // compute the zero-crossing and linearly interpolate dynamics
                    dReal interp=-aprev[j]/adiff;
                    scur = discrsvect[i-1] + interp*(discrsvect[i]-discrsvect[i-1]);
                    bcur = bprev[j] + interp*(b[j]-bprev[j]);
                    ccur = cprev[j] + interp*(c[j]-cprev[j]);
                }

                dReal f = ccur/bcur;
                if(f<0) {
                    dReal sdstar = sqrt(-f);
                    dReal sdplus = SdLimitBobrowExclude(scur,j);
                    if(sdplus >0 && sdplus < sdstar) {
                        continue;
                    }
                    if( !found || sdstar < minsd ) {
                        found = true;
                        minsd = sdstar;
                        mins = scur;
                    }
                }
            }
        }
        if(found) {
            AddSwitchPoint2(mins, minsd, SP_SINGULAR);
        }
        aprev.swap(a);
        bprev.swap(b);
        cprev.swap(c);
    }
}



////////////////////////////////////////////////////////////////////
//////////////////////////// Profile ///////////////////////////////
////////////////////////////////////////////////////////////////////


Profile::Profile(const std::list<dReal>& slist, const std::list<dReal>& sdlist, const std::list<dReal>&  sddlist, dReal integrationtimestep, bool forward) : integrationtimestep(integrationtimestep), forward(forward) {
    if( forward ) {
        svect.insert(svect.end(), slist.begin(), slist.end());
        sdvect.insert(sdvect.end(), sdlist.begin(), sdlist.end());
        sddvect.insert(sddvect.end(), sddlist.begin(), sddlist.end());
    }
    else {
        svect.insert(svect.end(), slist.rbegin(), slist.rend());
        sdvect.insert(sdvect.end(), sdlist.rbegin(), sdlist.rend());
        sddvect.insert(sddvect.end(), sddlist.rbegin(), sddlist.rend());
    }
    BOOST_ASSERT(svect.size()>0);
    BOOST_ASSERT(svect.size()==sdvect.size());
    BOOST_ASSERT(svect.size()==sddvect.size());
    // TODO: handle the case of last step with variable width
    nsteps = svect.size();
    duration = integrationtimestep * (nsteps-1);
}

Profile::Profile(const std::vector<dReal>& svectin, const std::vector<dReal>& sdvectin, const std::vector<dReal>&  sddvectin, dReal integrationtimestep, bool forward) : integrationtimestep(integrationtimestep), forward(forward) {
    BOOST_ASSERT(svectin.size()>0);
    BOOST_ASSERT(svectin.size()==sdvectin.size());
    BOOST_ASSERT(svectin.size()==sddvectin.size());
    if( forward ) {
        svect = svectin;
        sdvect = sdvectin;
        sddvect = sddvectin;
    }
    else {
        svect.resize(svectin.size());
        std::copy(svectin.rbegin(), svectin.rend(), svect.begin());
        sdvect.resize(sdvectin.size());
        std::copy(sdvectin.rbegin(), sdvectin.rend(), sdvect.begin());
        sddvect.resize(sddvectin.size());
        std::copy(sddvectin.rbegin(), sddvectin.rend(), sddvect.begin());
    }
    // TODO: handle the case of last step with variable width
    nsteps = svect.size();
    duration = integrationtimestep * (nsteps-1);
}

bool Profile::FindTimestepIndex(dReal t, int &index, dReal& remainder) const {
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


bool Profile::Invert(dReal s,  dReal& t) const {
    if( svect.size() == 0 ) {
        return false;
    }
    if(s<svect[0]-TINY || s>svect.back()+TINY) {
        return false;
    }
    std::vector<dReal>::const_iterator it = std::lower_bound(svect.begin(), svect.end(), s);
    size_t index = it - svect.begin();
    if( index == 0 ) {
        t = 0;
    }
    else {
        index -= 1;
        dReal tres;
        if(!SolveQuadraticEquation(svect[index]-s,sdvect[index],0.5*sddvect[index],tres,0,integrationtimestep)) {
            //std::cout << "***************** Inversion warning: tres=" << tres << " while integrationtimestep= "<< integrationtimestep << "****************\n";
            SolveQuadraticEquation(svect[index]-s,sdvect[index],0.5*sddvect[index],tres,0,integrationtimestep);
        }
        t = index*integrationtimestep + tres;
    }
    return true;
}

dReal Profile::Eval(dReal t) const {
    int index;
    dReal remainder;
    if(FindTimestepIndex(t, index, remainder))
        return svect[index] + remainder*sdvect[index]
               + 0.5*remainder*remainder*sddvect[index];
    return INF;
}


dReal Profile::Evald(dReal t) const {
    int index;
    dReal remainder;
    if (FindTimestepIndex(t, index, remainder))
        return sdvect[index] + remainder * sddvect[index];
    return INF;
}


dReal Profile::Evaldd(dReal t) const {
    int index;
    dReal remainder;
    if (FindTimestepIndex(t, index, remainder))
        return sddvect[index];
    return INF;
}


void Profile::Print() const {
    for(dReal t=0; t<=duration; t+=integrationtimestep) {
        std::cout<< Eval(t) << " ; " << Evald(t) << " ; " << Evaldd(t) <<  "\n";
    }
}


void Profile::Write(std::stringstream& ss, dReal dt) const {
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
        const dReal magicconst = 3.1; // what is this??
        const dReal slopethresh = 0.1;
        const dReal stepthresh = 0.04; /// slopes within 0.02 of the singularity are also pretty big, so need longer to stabilize them.
        dReal sstep = 1e-3/magicconst;
        while(sstep <= stepthresh) {
            sstep = sstep * magicconst;
            for(int i = 0; i<int(switchpoint.slopesvector.size()); i++) {
                sforward = std::min(s + sstep,constraints.trajectory.duration);
                sbackward = std::max(s - sstep,0.);
                dReal slope = switchpoint.slopesvector[i];
                // Avoid too big stubs
                if(std::abs(slope*sstep)>slopethresh) {
                    continue;
                }
                sdforward = sd + (sforward-s)*slope;
                sdbackward = sd - (s-sbackward)*slope;
                bool canintegrate = false;
                int ret = IntegrateBackward(constraints,sbackward,sdbackward,dt,resprofile,constraints.passswitchpointnsteps);
                if(!(ret==INT_MAXSTEPS||ret==INT_END)) {
                    // need to try smaller discretization in case slope has noise
                    ret = IntegrateBackward(constraints,sbackward,sdbackward,dt*0.1,resprofile,constraints.passswitchpointnsteps);
                }
                if(ret==INT_MAXSTEPS||ret==INT_END) {
                    ret = IntegrateForward(constraints,sforward,sdforward,dt,resprofile,constraints.passswitchpointnsteps,testaboveexistingprofiles, testmvc, zlajpah);
                    if(!(ret==INT_MAXSTEPS||ret==INT_END)) {
                        // need to try smaller discretization in case slope has noise
                        ret = IntegrateForward(constraints,sforward,sdforward,dt*0.1,resprofile,constraints.passswitchpointnsteps,testaboveexistingprofiles, testmvc, zlajpah);
                    }
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
                    if(score<bestscore && score<10) {
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
        constraints.nsingulartreated++;
        return bestscore<INF;
    }
    // Tangent, Discontinuous and Zlajpah switchpoints
    else{
        dReal sdtop = -1;
        dReal testsd = sd;
        for(int itry = 0; itry < 5; ++itry) {
            testsd *= constraints.loweringcoef;
            sdtop = BisectionSearch(constraints,s,testsd,sd,constraints.integrationtimestep,0);
            if(sdtop>0) {
                break;
            }
        }
        if(sdtop<=0) {
            return false;
        }
        sbackward = s;
        sforward = s;
        sdbackward = sdtop;
        sdforward = sdtop;
        constraints.ntangenttreated++;
        return true;
    }
    return false;
}



////////////////////////////////////////////////////////////////////
////////////////// Recovery mechanism //////////////////////////////
////////////////////////////////////////////////////////////////////

// Determine whether one fill the gap
bool IntegrateRecover(Constraints& constraints, dReal s, dReal sd, dReal dt, Profile& profilebw, Profile& profilefw){
    int ret;
    ret = IntegrateBackward(constraints,s,sd,dt,profilebw);
    if(ret==INT_PROFILE) {
        ret = IntegrateForward(constraints,s,sd,dt,profilefw);
        if(ret==INT_PROFILE) {
            return true;
        }
    }
    return false;
}

bool RecoverGap(Constraints& constraints, dReal s, dReal dt){
    int ntries = 100;
    dReal sdtop =  constraints.SdLimitCombined(s);
    dReal dsd = sdtop / ntries;
    Profile profilebw, profilefw;
    for(dReal sd = sdtop; sd > 0; sd -= dsd) {
        if(IntegrateRecover(constraints, s, sd, dt, profilebw, profilefw)) {
            constraints.resprofileslist.push_back(profilebw);
            constraints.resprofileslist.push_back(profilefw);
            return true;
        }
    }
    return false;
}


dReal Recover(Constraints& constraints, dReal ds){
    dReal s = 0, gapstart = 0, gapend = constraints.trajectory.duration;
    while(s<=constraints.trajectory.duration) {
        ProfileSample lowestsample = FindLowestProfileFast(s, INF, constraints.resprofileslist);
        if( lowestsample.itprofile == constraints.resprofileslist.end() ) {
            gapstart = s;
            s += ds;
            while(s<=constraints.trajectory.duration) {
                ProfileSample lowestsample = FindLowestProfileFast(s, INF, constraints.resprofileslist);
                if( lowestsample.itprofile == constraints.resprofileslist.end() ) {
                    s += ds;
                }
                else{
                    gapend = s - ds;
                    break;
                }
            }
            if(!RecoverGap(constraints,(gapstart+gapend)/2, ds)) {
                return gapstart;
            }
            else{
                std::cout << "Recovered gap [" << gapstart << "," << gapend << "]\n";
            }
        }
        s += ds;
    }
    return -1;
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

    // //Check alpha
    // snext = s + dt*sd + 0.5*dtsq*alpha;
    // sdnext_int = sd + dt*alpha;
    // sdnext_mvc = constraints.SdLimitCombined(snext);
    // if(snext > constraints.trajectory.duration || snext<0) {
    //     //std::cout << "Compute slide fin traj\n";
    // 	if(sdnext_mvc < sdnext_int) {
    // 	    return INF; // hits the MVC right before the end of the trajectory
    // 	}
    //     return -INF; // really reach the end of the trajectory
    // }    
    // if(sdnext_mvc < sdnext_int) {
    //     //std::cout << "Cannot slide alpha \n";
    //     return INF;
    // }

    // //Check beta
    // snext = s + dt*sd + 0.5*dtsq*beta;
    // sdnext_int = sd + dt*beta;
    // if(snext > constraints.trajectory.duration || snext<0) {
    //     //std::cout << "Compute slide fin traj\n";
    //     return -INF;
    // }
    // sdnext_mvc = constraints.SdLimitCombined(snext);
    // if(sdnext_mvc > sdnext_int) {
    //     //std::cout << "Cannot slide beta \n";
    //     return beta;
    // }
    
    //Actually no need here to determine whether trajectory.duration is reached or not
    //Check alpha
    snext = s + dt*sd + 0.5*dtsq*alpha;
    sdnext_int = sd + dt*alpha;
    sdnext_mvc = constraints.SdLimitCombined(snext);
    if(sdnext_mvc < sdnext_int) {
	return INF;
    }
    
    //Check beta
    snext = s + dt*sd + 0.5*dtsq*beta;
    sdnext_int = sd + dt*beta;
    sdnext_mvc = constraints.SdLimitCombined(snext);
    if(sdnext_mvc > sdnext_int) {
	return INF;
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
        return alpha;
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
        return INF;
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
    // used the cached values for memory performance
    BOOST_ASSERT(!constraints._busingcache);
    std::vector<dReal>& svect=constraints._svectcache; svect.resize(0);
    std::vector<dReal>& sdvect=constraints._sdvectcache; sdvect.resize(0);
    std::vector<dReal>& sddvect=constraints._sddvectcache; sddvect.resize(0);
    constraints._busingcache = true;
    resprofile.Reset();

    bool cont = true;
    int returntype = -1; // should be changed

//    // Initialize the currentindex of the profiles for search purpose
//    if(testaboveexistingprofiles && constraints.resprofileslist.size()>0) {
//        std::list<Profile>::iterator it = constraints.resprofileslist.begin();
//        while(it != constraints.resprofileslist.end()) {
//            it++;
//        }
//    }

    // Integrate forward
    while(cont) {
        if(int(svect.size()) > maxsteps) {
            returntype = INT_MAXSTEPS;
            break;
        }
        else if(scur > constraints.trajectory.duration) {
            //TODO: change the time step of previous step to reach the end
	    // CHANGED
            svect.push_back(scur);
            sdvect.push_back(sdcur);
            sddvect.push_back(0);
            returntype = INT_END;
            break;
        }
        else if(sdcur < 0) {
            //TODO: double check whether alpha>beta
            svect.push_back(scur);
            sdvect.push_back(sdcur);
            sddvect.push_back(0);
            returntype = INT_BOTTOM;
            break;
        }
        else if(IsAboveProfilesList(scur,sdcur,constraints.resprofileslist)) {
            svect.push_back(scur);
            sdvect.push_back(sdcur);
            sddvect.push_back(0);
            returntype = INT_PROFILE;
            break;
        }
        else if(zlajpah && testmvc && sdcur >= constraints.SdLimitCombined(scur)-TINY2) {
            if(constraints.SdLimitBobrow(scur)-constraints.SdLimitCombined(scur)<TINY2) {
                svect.push_back(scur);
                sdvect.push_back(sdcur);
                sddvect.push_back(0);
                returntype = INT_MVC;
                break;
            }

            // Lower the sd to the MVC and integrate backward
            if(svect.size()==0) {
                sdcur = constraints.SdLimitCombined(scur);
                Profile tmpprofile;
                //std::cout << "Integrate backward from (" <<scur << "," << sdcur  <<  ")\n";
                std::vector<dReal> svecttemp = svect, sdvecttemp = sdvect, sddvecttemp = sddvect;
                constraints._busingcache = false;
                int res3 = IntegrateBackward(constraints,scur,sdcur,constraints.integrationtimestep,tmpprofile,1e5,true,true,true);
                svect = svecttemp;
                sdvect = sdvecttemp;
                sddvect = sddvecttemp;
                constraints._busingcache = true;

                if(res3 == INT_BOTTOM) {
                    //std::cout << "BW reached 0 (From Zlajpah)\n";
                    constraints._busingcache = false;
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
            //a) alpha points above MVCCombined (trap case); then step along MVCCombined until alpha points below MVCCombined 
	    //   and integrate backward from that point. Cf. Zlajpah ICRA 1996.
            //b) alpha points below MVCCombined (slide case); then slide along MVCCombined, which is admissible, until either
            // b1) alpha points above MVCCombined (trapped)
            // b2) beta points below MVCCombined (exit slide)

            int res = FlowVsMVC(constraints,scur,sdcur,1,dt);
	    //FlowVsMVC
	    //flag = 1 : alpha flow
	    //return = -1 if flow points below MVC
	    //return = +1 if flow points above MVC
	    //return = 0 if end of trajectory

            if(res == 0) {
                // Most probably we arrived at the end
                svect.push_back(scur);
                sdvect.push_back(sdcur);
                sddvect.push_back(0);
                returntype = INT_END;
                break;
            }
            else if(res == 1) {
                // Case a
                // Insert the profile calculated so far into the resprofileslist
                if( svect.size() > 1 ) {
                    // And reinitialize the profile
                    constraints.resprofileslist.push_back(Profile(svect,sdvect,sddvect,dt, true));
                }
                svect.resize(0);
                sdvect.resize(0);
                sddvect.resize(0);
                // Now step along the MVCCombined
                while(true) {
                    snext = scur + dt;
                    if(snext>constraints.trajectory.duration) {
                        cont = false;
                        svect.push_back(scur);
                        sdvect.push_back(sdcur);
                        sddvect.push_back(0);
                        returntype = INT_END; // ZLAJPAH END
			// should the return type be INT_MVC instead?
			// because what reaches the end is not the integrated profile
                        break;
                    }
                    sdnext = constraints.SdLimitCombined(snext);
                    int res2 = FlowVsMVC(constraints,snext,sdnext,1,dt);
                    if(res2 == 0) {
			// how come? we have already check if snext > traj.duration
                        cont = false;
                        svect.push_back(scur);
                        sdvect.push_back(sdcur);
                        sddvect.push_back(0);
                        returntype = INT_END; // ZLAJPAH END
			// should the return type be INT_MVC instead?
			// because what reaches the end is not the integrated profile
                        break;
                    }
                    else if(res2 == -1) {
                        // Alpha points below the MVC
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
		// slide (integrate) along the MVC until the profile hits something
                while(true) {
                    if(int(svect.size()) > maxsteps) {
                        cont = false;
                        returntype = INT_MAXSTEPS;
                        break;
                    }
                    else if(scur >= constraints.trajectory.duration) {
                        cont = false;
                        svect.push_back(scur);
                        sdvect.push_back(sdcur);
                        sddvect.push_back(0);
                        returntype = INT_END;
                        break;
                    }
                    else if(sdcur < 0) {
			// maybe some numerical error happens here
			// because if the MVC does not hits zero, how would this happen?
                        cont = false;
                        svect.push_back(scur);
                        sdvect.push_back(sdcur);
                        sddvect.push_back(0);
                        returntype = INT_BOTTOM;
                        break;
                    }
                    else if(IsAboveProfilesList(scur,sdcur,constraints.resprofileslist)) {
                        cont = false;
                        svect.push_back(scur);
                        sdvect.push_back(sdcur);
                        sddvect.push_back(0);
                        returntype = INT_PROFILE;
                        break;
                    }
                    dReal slidesdd = ComputeSlidesdd(constraints,scur,sdcur,dt);
		    // if beta points below the MVC, slidesdd already returns beta
                    if(slidesdd == INF) {
			// slidesdd == INF in case alpha is pointing above the MVC
                        cont = false;
                        svect.push_back(scur);
                        sdvect.push_back(sdcur);
                        sddvect.push_back(0);
                        returntype = INT_MVC;
                        break;
                    }
                    snext = scur + dt * sdcur + 0.5*dtsq*slidesdd;
                    sdnext = sdcur + dt * slidesdd;
                    svect.push_back(scur);
                    sdvect.push_back(sdcur);
                    sddvect.push_back(slidesdd);
                    scur = snext;
                    sdcur = sdnext;
		    
		    // check beta flow
		    int res1 = FlowVsMVC(constraints,snext,sdnext,2,dt);
		    if (res1 == -1) {
			// exit sliding (and continue with normal integration)
			break;
		    }

		    // it seems like the following if-elses are not neccessary
		    
                    // //std::cout <<"Next ("<< snext << "," << sdnext << ") \n";
                    // int res1 = FlowVsMVC(constraints,snext,sdnext,1,dt);
                    // if(res1 == 0) {
                    //     cont = false;
                    //     //std::cout << "End traj\n";
		    // 	returntype = INT_END;
                    //     break;
                    // }
                    // else if(res1 == 1) {
                    //     // Case b1
                    //     //std::cout << "End slide with trap (" <<snext << "," << sdnext  <<  ")\n";
                    //     break;
                    // }
                    // int res2 = FlowVsMVC(constraints,snext,sdnext,2,dt);
                    // if(res2 == 0) {
                    //     cont = false;
                    //     //std::cout << "End traj\n";
                    //     break;
                    // }
                    // else if(res2 == -1) {
                    //     // Case b2
                    //     //std::cout << "End slide with exit (" <<snext << "," << sdnext  <<  ")\n";
                    //     break;
                    // }
                }
            }
        }
        else if((!zlajpah) && testmvc && sdcur > constraints.SdLimitCombined(scur)+TINY2) {
            svect.push_back(scur);
            sdvect.push_back(sdcur);
            sddvect.push_back(0);
            constraints.zlajpahlist.push_back(std::pair<dReal,dReal>(scur,sdcur));
            returntype = INT_MVC;
            break;
        }
        else{
            svect.push_back(scur);
            sdvect.push_back(sdcur);
            dReal beta = constraints.SddLimitBeta(scur,sdcur);
            sddvect.push_back(beta);
            dReal snext = scur + dt * sdcur + 0.5*dtsq*beta;
            dReal sdnext = sdcur + dt * beta;
	    
	    // adjust the last step size such that snext is constraints.trajectory.duration
	    if (snext > constraints.trajectory.duration) {
		dReal send = constraints.trajectory.duration;
		dReal dtnew;
		if (!SolveQuadraticEquation(scur - send, scur, 0.5*beta, dtnew, 0.0)) {
		    std::cout << "[TOPP::IntegrateForward] Solving for dtnew failed.\n";
		    scur = snext;
		    sdcur = sdnext;
		}
		else {
		    sdnext = sdcur + dtnew * beta;
		    scur = send;
		    sdcur = sdnext;
		}
	    }
	    else {
		scur = snext;
		sdcur = sdnext;
	    }
        }
    }
    resprofile = Profile(svect,sdvect,sddvect,dt, true);
    constraints._busingcache = false;
    return returntype;
}


int IntegrateBackward(Constraints& constraints, dReal sstart, dReal sdstart, dReal dt, Profile& resprofile,  int maxsteps, bool testaboveexistingprofiles, bool testmvc,bool zlajpah){
    dReal dtsq = dt*dt;
    dReal scur = sstart, sdcur = sdstart;
    //std::list<dReal> slist, sdlist, sddlist;
    BOOST_ASSERT(!constraints._busingcache);
    std::vector<dReal>& svect=constraints._svectcache; svect.resize(0);
    std::vector<dReal>& sdvect=constraints._sdvectcache; sdvect.resize(0);
    std::vector<dReal>& sddvect=constraints._sddvectcache; sddvect.resize(0);
    constraints._busingcache = true;
    bool cont = true;
    int returntype = INT_END;
    bool searchbackward = true;
    dReal alphabk = INF;
    resprofile.Reset();

//    // Initialize the currentindex of the profiles for search purpose
//    if(testaboveexistingprofiles && constraints.resprofileslist.size()>0) {
//        std::list<Profile>::iterator it = constraints.resprofileslist.begin();
//        while(it != constraints.resprofileslist.end()) {
//            it++;
//        }
//    }

    // Integrate backward
    while(cont) {
        if(int(svect.size()) > maxsteps) {
            returntype = INT_MAXSTEPS;
            break;
        }
        else if(scur < 0) {
            //TODO: change the time step of previous step to reach the end
	    // CHANGED
            svect.push_back(scur);
            sdvect.push_back(sdcur);
            sddvect.push_back(0);
            returntype = INT_END;
            break;
        }
        else if(sdcur < 0) {
            //TODO: double check whether alpha>beta
            svect.push_back(scur);
            sdvect.push_back(sdcur);
            sddvect.push_back(0);
            returntype = INT_BOTTOM;
            break;
        }
        else if(IsAboveProfilesList(scur,sdcur,constraints.resprofileslist,searchbackward)) {
            svect.push_back(scur);
            sdvect.push_back(sdcur);
            sddvect.push_back(0);
            returntype = INT_PROFILE;
            break;
        }
        else if(zlajpah && testmvc && sdcur >= constraints.SdLimitCombined(scur)-TINY2 && FlowVsMVCBackward(constraints,scur,sdcur,dt) != -1) {
            if(sdcur > constraints.SdLimitBobrow(scur)) {
                svect.push_back(scur);
                sdvect.push_back(sdcur);
                sddvect.push_back(0);
                returntype = INT_MVC;
                break;
            }
            while(true) {
                if(scur <0) {
                    cont = false;
                    svect.push_back(scur);
                    sdvect.push_back(sdcur);
                    sddvect.push_back(0);
                    returntype = INT_END;
                    break;
                }
                else if(sdcur < 0) {
                    cont = false;
                    svect.push_back(scur);
                    sdvect.push_back(sdcur);
                    sddvect.push_back(0);
                    returntype = INT_BOTTOM;
                    break;
                }
                if(sdcur > constraints.SdLimitBobrow(scur)) {
                    svect.push_back(scur);
                    sdvect.push_back(sdcur);
                    sddvect.push_back(0);
                    returntype = INT_MVC;
                    break;
                }
                else if(IsAboveProfilesList(scur,sdcur,constraints.resprofileslist)) {
                    cont = false;
                    svect.push_back(scur);
                    sdvect.push_back(sdcur);
                    sddvect.push_back(0);
                    returntype = INT_PROFILE;
                    break;
                }
                //std::cout <<"Slide from ("<< scur << "," << sdcur << ") \n";
                dReal slidesdd = ComputeSlidesddBackward(constraints,scur,sdcur,dt);
                if(slidesdd == INF) {
                    cont = false;
                    svect.push_back(scur);
                    sdvect.push_back(sdcur);
                    sddvect.push_back(0);
                    returntype = INT_MVC;
                    break;
                }
                dReal sprev = scur - dt * sdcur + 0.5*dtsq*slidesdd;
                dReal sdprev = sdcur - dt * slidesdd;
                svect.push_back(scur);
                sdvect.push_back(sdcur);
                sddvect.push_back(slidesdd);
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
            svect.push_back(scur);
            sdvect.push_back(sdcur);
            sddvect.push_back(0);
            returntype = INT_MVC;
            break;
        }
        else{
            svect.push_back(scur);
            sdvect.push_back(sdcur);
            //std ::pair<dReal,dReal> sddlimits = constraints.SddLimits(scur,sdcur);
            //dReal alpha = sddlimits.first;
            dReal alpha = constraints.SddLimitAlpha(scur,sdcur);
            if(sddvect.size()==0) {
                alphabk = alpha;
            }
            sddvect.push_back(alpha);
            //std::cout << scur << " " << sdcur << " " << alpha << "\n";
            dReal sprev = scur - dt * sdcur + 0.5*dtsq*alpha;
            dReal sdprev = sdcur - dt * alpha;

	    //alpha = constraints.SddLimitAlpha(std::min(constraints.trajectory.duration,std::max(0.,scur)),std::max(0.,sdcur));

	    // adjust the last step size such that sprev is 0
	    if (sprev < 0) {
		dReal sbeg = 0;
		dReal dtnew;
		if (!SolveQuadraticEquation(scur - sbeg, -scur, 0.5*alpha, dtnew, 0.0)) {
		    std::cout << "[TOPP::IntegrateBackward] Solving for dtnew failed.\n";
		    scur = sprev;
		    sdcur = sdprev;
		}
		else{
		    sdprev = sdcur - dtnew * alpha;
		    scur = sbeg;
		    sdcur = sdprev;
		}
	    }
	    else {
		scur = sprev;
		sdcur = sdprev;
	    }
        }
    }
    if(sddvect.size()>0) {
        dReal sddfront = alphabk == INF ? sddvect.front() : alphabk;
        sddvect.pop_back();
        sddvect.insert(sddvect.begin(), sddfront);
        resprofile = Profile(svect,sdvect,sddvect,dt, false);
    }

    constraints._busingcache = false;
    return returntype;
}


////////////////////////////////////////////////////////////////////
//////////////////////// Limiting curves ///////////////////////////
////////////////////////////////////////////////////////////////////

Profile StraightProfile(dReal sbackward,dReal sforward, dReal sdbackward, dReal sdforward)
{
    std::vector<dReal> slist, sdlist, sddlist;
    dReal dtmod = 2 * (sforward - sbackward) / (sdforward + sdbackward);
    dReal sdd = (sdforward - sdbackward) / dtmod;
    slist.push_back(sbackward);
    slist.push_back(sforward);
    sdlist.push_back(sdbackward);
    sdlist.push_back(sdforward);
    sddlist.push_back(sdd);
    sddlist.push_back(0);
    return Profile(slist,sdlist,sddlist,dtmod,true);
}


int ComputeLimitingCurves(Constraints& constraints){
    //std::list<SwitchPoint> switchpointslist0(constraints.switchpointslist);

    Profile tmpprofile;
    dReal sswitch, sdswitch, sbackward, sdbackward, sforward, sdforward;
    int integratestatus;
    bool testaboveexistingprofiles = true, testmvc = true, zlajpah = false;
    constraints.ntangenttreated = 0;
    constraints.nsingulartreated = 0;
    int switchpointindex = 0; // for debugging purposes;
    for(std::list<SwitchPoint>::const_iterator itswitchpoint = constraints.switchpointslist.begin(); itswitchpoint != constraints.switchpointslist.end(); ++itswitchpoint, ++switchpointindex) {
        //while(switchpointslist0.size() > 0) {
        //SwitchPoint switchpoint = switchpointslist0.front();
        //switchpointslist0.pop_front();
        const SwitchPoint& switchpoint = *itswitchpoint;
        sswitch = switchpoint.s;
        sdswitch = switchpoint.sd;
        if(IsAboveProfilesList(sswitch,sdswitch,constraints.resprofileslist,false,0))
            continue;
        if(sdswitch > constraints.SdLimitCombined(sswitch)+TINY2)
            continue;

        // Address Switch Point
        if (!AddressSwitchPoint(constraints, switchpoint, sbackward,
                                sdbackward, sforward, sdforward))
            continue;

        bool shiller = false;
        bool rescale = false;

        if(!shiller) {
            // Add middle part
            if(rescale) {
                dReal slope = (sdforward-sdbackward)/(sforward-sbackward);
                sbackward = sswitch - constraints.integrationtimestep;
                sforward = sswitch + constraints.integrationtimestep;
                sdbackward = sdswitch - slope*constraints.integrationtimestep;
                sdforward = sdswitch + slope*constraints.integrationtimestep;
                std::cout << "toto " << constraints.integrationtimestep << "\n";
            }
            if (sforward - sbackward > TINY) {
                constraints.resprofileslist.push_back(StraightProfile(sbackward,sforward,sdbackward,sdforward));
            }
        }
        else{
            sbackward = sswitch - constraints.integrationtimestep;
            sforward = sswitch + constraints.integrationtimestep;
            sdbackward = constraints.SdLimitCombined(sbackward);
            sdforward = constraints.SdLimitCombined(sforward);
            constraints.resprofileslist.push_back(StraightProfile(sswitch,sforward,sdswitch,sdforward));
            constraints.resprofileslist.push_back(StraightProfile(sbackward,sswitch,sdbackward,sdswitch));
        }

        // Integrate backward
        integratestatus = IntegrateBackward(constraints, sbackward, sdbackward,
                                            constraints.integrationtimestep, tmpprofile);
        if(tmpprofile.nsteps>2) {
            constraints.resprofileslist.push_back(tmpprofile);
        }

        if(integratestatus == INT_BOTTOM) {
            std::cerr << str(boost::format("IntegrateBackward INT_BOTTOM, s=%.15e, sd=%.15e\n")%sbackward%sdbackward);
            return CLC_BOTTOM;
        }

        // Integrate forward
        integratestatus = IntegrateForward(constraints, sforward, sdforward,
                                           constraints.integrationtimestep, tmpprofile, 1e5,
                                           testaboveexistingprofiles, testmvc, zlajpah);
        if(tmpprofile.nsteps>2) {
            constraints.resprofileslist.push_back(tmpprofile);
            //std::cout << "Forward : " << tmpprofile.nsteps << " " << tmpprofile.duration << "\n";
        }

        if(integratestatus == INT_BOTTOM) {
            std::cerr << str(boost::format("IntegrateForward INT_BOTTOM, s=%.15e, sd=%.15e\n")%sforward%sdforward);
            return CLC_BOTTOM;
        }
    }
    return CLC_OK;
}


int ComputeProfiles(Constraints& constraints, dReal sdbeg, dReal sdend){
    std::chrono::time_point<std::chrono::system_clock> t0,t1,t2,t3;
    std::chrono::duration<double> d1,d2,d3,dtot;

    t0 = std::chrono::system_clock::now();

    bool retprocess = constraints.Preprocess();
    if(!retprocess) {
        std::cout << "Cannot preprocess\n";
        return TOPP_CANNOT_PREPROCESS;
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
            constraints.stepthresh *= 2;
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
            message = "CLC failed";
            std::cout << message << std::endl;
            integrateprofilesstatus = false;
            continue;
        }


        Profile tmpprofile;

        /////////////////  Integrate from start /////////////////////
        ret = IntegrateForward(constraints,0,sdbeg,constraints.integrationtimestep,resprofile,1e5,testaboveexistingprofiles,testmvc,zlajpah);
        if(ret != INT_BOTTOM && resprofile.nsteps>2) {
            constraints.resprofileslist.push_back(resprofile);
        }
        if(ret == INT_BOTTOM || resprofile.nsteps < 5) {
            // Integration failed. However, if qd(0) = 0,  there was probably a singularity, and one can try different values for sdot
            std::vector<dReal> qd(constraints.trajectory.dimension);
            constraints.trajectory.Evald(0,qd);
            if(VectorNorm(qd) <= 5e-2) {
                dReal s = 1e-2;
                dReal ntrials = 1000.;
                dReal mvcsd = std::min(constraints.SdLimitCombined(s),10.);
                dReal incr = mvcsd / (ntrials+0.001);
                for(dReal trialsd = incr; trialsd < mvcsd; trialsd += incr) {
                    dReal betas = constraints.SddLimitBeta(s,trialsd);
                    dReal betas2 = constraints.SddLimitBeta(s*2,trialsd+betas*s);
                    if(betas2 > betas) {
                        constraints.resprofileslist.push_back(StraightProfile(0,s,trialsd-s*betas,trialsd));
                        ret = IntegrateForward(constraints,s,trialsd,constraints.integrationtimestep,resprofile,1e5,testaboveexistingprofiles,testmvc,zlajpah);
                        if(resprofile.nsteps>1) {
                            constraints.resprofileslist.push_back(resprofile);
                        }
                        break;
                    }
                }
            }
        }
        // Now if it still fails, shouganai
        if(ret==INT_BOTTOM) {
            message = "Start profile hit 0";
            std::cout << message << std::endl;
            integrateprofilesstatus = false;
            continue;
        }


        /////////////////  Integrate from end /////////////////////
        ret = IntegrateBackward(constraints,constraints.trajectory.duration,sdend,constraints.integrationtimestep,resprofile,1e5,testaboveexistingprofiles,testmvc);
        if(ret != INT_BOTTOM && resprofile.nsteps>2) {
            constraints.resprofileslist.push_back(resprofile);
        }
        if(ret == INT_BOTTOM || resprofile.nsteps < 5) {
            // Integration failed. However, if qd(send) = 0, there was probably a singularity, and one can try different values for sdot
            dReal send = constraints.trajectory.duration;
            std::vector<dReal> qd(constraints.trajectory.dimension);
            constraints.trajectory.Evald(send,qd);
            //std::cout << send << " " << qd[0] << " " << qd[1] << " " << qd[2]  << " VN\n";
            if(VectorNorm(qd) <= 5e-2) {
                dReal ds = 1e-2;
                dReal s = send - ds;
                dReal ntrials = 1000.;
                dReal mvcsd = std::min(constraints.SdLimitCombined(s),10.);
                dReal incr = mvcsd / (ntrials+0.001);
                for(dReal trialsd = incr; trialsd < mvcsd; trialsd += incr) {
                    dReal alphas = constraints.SddLimitAlpha(s,trialsd);
                    dReal alphas2 = constraints.SddLimitAlpha(s-ds,trialsd-alphas*ds);
                    //std::cout << trialsd << " " << alphas << " " << alphas2 << "\n";
                    if(alphas2 < alphas) {
                        constraints.resprofileslist.push_back(StraightProfile(s,send,trialsd,trialsd+alphas*ds));
                        ret = IntegrateBackward(constraints,s,trialsd,constraints.integrationtimestep,resprofile,1e5,testaboveexistingprofiles,testmvc);
                        if(resprofile.nsteps>1) {
                            constraints.resprofileslist.push_back(resprofile);
                        }
                        break;
                    }
                }

            }
        }
        if(ret==INT_BOTTOM) {
        }
        // Now if it still fails, shouganai
        if(ret==INT_BOTTOM) {
            message = "End profile hit 0";
            std::cout << message << std::endl;
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
            std::cout << message << std::endl;
            integrateprofilesstatus = false;
            continue;
        }


        /////////////////////  Final checks /////////////////////////
        // Check whether CLC is continuous and recover if not
        dReal sdiscontinuous = Recover(constraints,constraints.integrationtimestep);
        if(sdiscontinuous>-0.5) {
            message = str(boost::format("Could not recover from CLC discontinuous s=%.15e")%sdiscontinuous);
        }


        // Estimate resulting trajectory duration
        constraints.resduration = 0;
        //Profile profile;
        dReal ds = constraints.integrationtimestep;
        int nsamples = int((constraints.trajectory.duration-TINY)/ds);
        dReal s,sdcur,sdnext;
        ProfileSample lowestsample = FindLowestProfileFast(0, INF, constraints.resprofileslist);
        if( lowestsample.itprofile != constraints.resprofileslist.end() ) {
            sdcur = lowestsample.sd;
        }
        else {
            message = "CLC discontinuous at 0";
            std::cout << message << std::endl;
            integrateprofilesstatus = false;
            continue;
        }
        bool clcdiscontinuous = false;
        for(int i=1; i<=nsamples; i++) {
            s = i*ds;
            ProfileSample lowestsample = FindLowestProfileFast(s, INF, constraints.resprofileslist);
            if( lowestsample.itprofile != constraints.resprofileslist.end() ) {
                sdnext = lowestsample.sd;
                constraints.resduration += 2*ds/(sdcur+sdnext);
                sdcur = sdnext;
            }
            else{
                clcdiscontinuous = true;
                lowestsample = FindLowestProfileFast(s, INF, constraints.resprofileslist);
                break;
            }
        }
        if(clcdiscontinuous) {
            message = str(boost::format("CLC discontinuous s=%.15e")%s);
            std::cout << message << std::endl;
            integrateprofilesstatus = false;
            continue;
        }

        integrateprofilesstatus = true;
        break;

    }


    if(!integrateprofilesstatus) {
        std::cout << "[TOPP::ComputeProfiles] failed: " << message << std::endl;
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
        std::cout << "Cannot preprocess\n";
        return TOPP_CANNOT_PREPROCESS;
    }

    if(VectorMin(constraints.mvcbobrow) <= TINY) {
        std::cout << "[TOPP::VIP] MVCBobrow hit 0 \n";
        return TOPP_MVC_HIT_ZERO;
    }

    Profile tmpprofile;
    //dReal tres;
    dReal smallincrement = constraints.integrationtimestep*2;

    // Compute the limiting curves
    int resclc = ComputeLimitingCurves(constraints);
    if(resclc == CLC_SWITCH || resclc == CLC_BOTTOM) {
        std::cout << "[TOPP::VIP] resclc == CLC_SWITCH or CLC_BOTTOM \n";
        return TOPP_CLC_ERROR;
    }

    // dReal sdiscontinuous = Recover(constraints,constraints.integrationtimestep);
    // if(sdiscontinuous>-0.5) {
    //     std::cout <<  "Could not recover from CLC discontinuous\n";
    //     return TOPP_CLC_ERROR;
    // }


    // Determine the lowest profile at t=0
    dReal bound;
    ProfileSample lowestsample = FindLowestProfileFast(smallincrement, INF, constraints.resprofileslist);
    if( lowestsample.itprofile != constraints.resprofileslist.end() ) {
        bound = std::min(lowestsample.sd,constraints.mvccombined[0]);
    }
    else { // just to make sure the profile is below mvccombined
        bound = constraints.mvccombined[0];
    }

//    if(FindLowestProfile(smallincrement,tmpprofile,tres,constraints.resprofileslist))
//        bound = std::min(tmpprofile.Evald(tres),constraints.mvccombined[0]);
//    else // just to make sure the profile is below mvccombined
//        bound = constraints.mvccombined[0];

    if(sdbegmin>bound) {
        std::cout << "[TOPP::VIP] sdbegmin is above the CLC or the combined MVC \n";
        return TOPP_SDBEGMIN_TOO_HIGH;
    }

    sdbegmax = std::min(sdbegmax,bound-constraints.bisectionprecision);

    bool testaboveexistingprofiles = true, testmvc = true, zlajpah = true;

    // Compute sdendmax by integrating forward from (0,sdbegmax)
    int resintfw = IntegrateForward(constraints,0,sdbegmax,constraints.integrationtimestep,tmpprofile,1e5,testaboveexistingprofiles,testmvc,zlajpah);

    if (resintfw == INT_BOTTOM) std::cout << "[VIP] INT_BOTTOM\n";
    else if (resintfw == INT_END) std::cout << "[VIP] INT_END\n";
    else if (resintfw == INT_MVC) std::cout << "[VIP] INT_MVC\n";
    else if (resintfw == INT_PROFILE) std::cout << "[VIP] INT_PROFILE\n";
    else if (resintfw == INT_MAXSTEPS) std::cout << "[VIP] INT_MAXSTEPS\n";
    else std::cout << "[VIP] resintfw = " << resintfw << "\n";
    if (tmpprofile.Evald(tmpprofile.duration) <= constraints.mvccombined[constraints.mvccombined.size() - 1]) {
	std::cout << "tmpprofile is not above the MVC\n";
    }
    else {
	std::cout << "tmpprofile is above the MVC\n";
    }

    constraints.resprofileslist.push_back(tmpprofile);
    if(resintfw == INT_BOTTOM) {
        std::cout << "[TOPP::VIP] Forward integration hit sd = 0 \n";
        return TOPP_FWD_HIT_ZERO;
    }
    else if (resintfw == INT_END && tmpprofile.Evald(tmpprofile.duration) <= constraints.mvccombined[constraints.mvccombined.size() - 1]) {
        sdendmax = tmpprofile.Evald(tmpprofile.duration);
	std::cout << "[VIP] sdendmax " << sdendmax << "\n";
	std::cout << "[VIP] tmpprofile.duration " << tmpprofile.duration << "\n";
    }
    else if (resintfw == INT_MVC || resintfw == INT_PROFILE) {
        // Look for the lowest profile at the end
//        if(FindLowestProfile(constraints.trajectory.duration-smallincrement,tmpprofile,tres,constraints.resprofileslist) && tmpprofile.Evald(tres) <= constraints.mvccombined[constraints.mvccombined.size()-1])
//            sdendmax = tmpprofile.Evald(tres);

        ProfileSample lowestsample = FindLowestProfileFast(constraints.trajectory.duration-smallincrement, constraints.mvccombined[constraints.mvccombined.size()-1], constraints.resprofileslist);
        if( lowestsample.itprofile != constraints.resprofileslist.end() ) {
            sdendmax = lowestsample.sd;
        }
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
        std::cout << "Cannot preprocess\n";
        return TOPP_CANNOT_PREPROCESS;
    }

    if(VectorMin(constraints.mvcbobrow) <= TINY) {
        std::cout << "[TOPP::VIPBackward] MVCBobrow hit 0 \n";
        return TOPP_MVC_HIT_ZERO;
    }

    Profile tmpprofile;
    //dReal tres;
    dReal smallincrement = constraints.integrationtimestep*2;

    // Compute the limiting curves
    int resclc = ComputeLimitingCurves(constraints);
    if (resclc == CLC_SWITCH || resclc == CLC_BOTTOM) {
        std::cout << "[TOPP::VIPBackward] resclc == CLC_SWITCH or CLC_BOTTOM \n";
        return TOPP_CLC_ERROR;
    }

    // Determine the lowest profile at send
    dReal bound;
    ProfileSample lowestsample = FindLowestProfileFast(constraints.trajectory.duration-smallincrement, INF, constraints.resprofileslist);
    if( lowestsample.itprofile != constraints.resprofileslist.end() ) {
        bound = std::min(lowestsample.sd, constraints.mvccombined[constraints.mvccombined.size() - 1]);
    }
    else {
        bound = constraints.mvccombined[constraints.mvccombined.size() - 1];
    }
//    if (FindLowestProfile(constraints.trajectory.duration-smallincrement, tmpprofile, tres, constraints.resprofileslist))
//        bound = std::min(tmpprofile.Evald(tres), constraints.mvccombined[constraints.mvccombined.size() - 1]);
//    else // just to make sure the profile is below mvccombined
//        bound = constraints.mvccombined[constraints.mvccombined.size() - 1];

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
//        if (FindLowestProfile(smallincrement, tmpprofile, tres, constraints.resprofileslist) && tmpprofile.Evald(tres) <= constraints.mvccombined[0])
//            sdbegmax = tmpprofile.Evald(tres);

        ProfileSample lowestsample = FindLowestProfileFast(smallincrement, constraints.mvccombined[0], constraints.resprofileslist);
        if( lowestsample.itprofile != constraints.resprofileslist.end() ) {
            sdbegmax = lowestsample.sd;
        }
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


dReal EmergencyStop(Constraints& constraints, dReal sdbeg, Trajectory& restrajectory) {

    dReal dt = constraints.integrationtimestep;
    if(dt == 0) {
        dt = constraints.discrtimestep;
    }

    int ndiscrsteps = int((constraints.trajectory.duration+1e-10)/constraints.discrtimestep);
    if(ndiscrsteps<1) {
        return false;
    }

    constraints.discrtimestep = constraints.trajectory.duration/ndiscrsteps;
    constraints.Discretize();

    dReal dtsq = dt*dt;

    dReal scur = 0, sdcur = sdbeg, snext, sdnext, alpha, beta, sprev = 0;
    std::list<dReal> slist, sdlist, sddlist;
    std::pair<dReal,dReal> sddlimits;
    std::vector<dReal> qd(constraints.trajectory.dimension);

    bool cont = true;
    int returntype = -1;

    while(cont) {

        if(sdcur<0) {
            // Could achieve emergency stop
            slist.push_back(scur);
            sdlist.push_back(sdcur);
            sddlist.push_back(0);
            returntype = INT_BOTTOM;
            break;
        }

        if(scur>constraints.trajectory.duration) {
            // Reached the end of the trajectory
            slist.push_back(scur);
            sdlist.push_back(sdcur);
            sddlist.push_back(0);
            std::cout << "[ES] Reached end\n";
            returntype = INT_END;
            break;
        }

        if(constraints.hasvelocitylimits) {
            constraints.trajectory.Evald(scur, qd);
            for(int i=0; i<constraints.trajectory.dimension; i++) {
                // Violated velocity constraint
                if(std::abs(qd[i])>TINY && sdcur > 1e-2 + constraints.vmax[i]/std::abs(qd[i])) {
                    slist.push_back(scur);
                    sdlist.push_back(sdcur);
                    sddlist.push_back(0);
                    std::cout << sdcur << ">" << constraints.vmax[i]/std::abs(qd[i]) << " " << "[ES] Violated velocity constraint\n";
                    returntype = INT_MVC;
                    break;
                }
            }
            if(returntype == INT_MVC) {
                break;
            }
        }

        sddlimits = constraints.SddLimits(scur,sdcur);
        alpha = sddlimits.first;
        beta = sddlimits.second;
        if(alpha>beta + 1e-2) {
            // Violated the acceleration constraints
            slist.push_back(scur);
            sdlist.push_back(sdcur);
            sddlist.push_back(0);
            std::cout << "[ES] Violated acceleration constraint\n";
            returntype = INT_MVC;
            break;
        }

        // Else, integrate forward following alpha
        slist.push_back(scur);
        sdlist.push_back(sdcur);
        sddlist.push_back(alpha);
        snext = scur + dt * sdcur + 0.5*dtsq*alpha;
        sdnext = sdcur + dt * alpha;
        sprev = scur;
        scur = snext;
        sdcur = sdnext;
    }

    if(returntype == INT_BOTTOM) {
        constraints.resprofileslist.resize(0);
        constraints.resprofileslist.push_back(Profile(slist,sdlist,sddlist,dt,true));
        int ret = constraints.trajectory.Reparameterize(constraints, restrajectory, sprev);
        if(ret > 0) {
            return sprev;
        }
        else{
            return 0;
        }
    }
    else{
        return 0;
    }


}

////////////////////////////////////////////////////////////////////
///////////////////////// Utilities ////////////////////////////////
////////////////////////////////////////////////////////////////////


void VectorFromString(std::string& s,std::vector<dReal>&resvect){
    s.erase(std::find_if(s.rbegin(), s.rend(), std::bind1st(std::not_equal_to<char>(), ' ')).base(), s.end()); //remove trailing spaces
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


void ReadVectorFromStream(std::istream& s, size_t N, std::vector<dReal>&resvect)
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


void PrintVector(const std::vector<dReal>&v){
    std::cout << "[";
    for(int i=0; i<int(v.size()); i++) {
        std::cout<< v[i] << ", ";
    }
    std::cout << "] \n ";
}

void VectorAdd(const std::vector<dReal>&a, const std::vector<dReal>&b,  std::vector<dReal>&res, dReal coefa, dReal coefb){
    res.resize(a.size());
    BOOST_ASSERT(a.size() == b.size());
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


bool SolveQuadraticEquation(dReal a0, dReal a1, dReal a2, dReal& sol, dReal lowerbound, dReal upperbound, dReal epsilon) {
    dReal delta2 = a1*a1- 4*a0*a2;
    if(delta2 < -epsilon) {
        return false;
    }
    delta2 = std::max(delta2,0.);
    if(a2<0) {
        a2=-a2;
        a1=-a1;
        a0=-a0;
    }
    if(a2<=epsilon) {
        if(std::abs(a1)<=epsilon) {
            return false; // Must be in a weird case
        }
        sol = -a0/a1;
        return sol>=lowerbound-epsilon && sol<=upperbound+epsilon;
    }
    dReal delta = sqrt(delta2);
    sol = (-a1-delta)/(2*a2);
    if(sol>=lowerbound-epsilon && sol<=upperbound+epsilon) {
        return true;
    }
    if(sol>upperbound+epsilon) {
        return false;
    }
    sol = (-a1+delta)/(2*a2);
    return sol>=lowerbound-epsilon && sol<=upperbound+epsilon;
}

bool IsAboveProfilesList(dReal s, dReal sd, const std::list<Profile>&resprofileslist, bool searchbackward, dReal softborder){
    ProfileSample sample = FindLowestProfileFast(s, INF, resprofileslist);
    if( sample.itprofile != resprofileslist.end() ) {
        if( sd > sample.sd+softborder ) {
            return true;
        }
    }
    return false;
//    dReal t;
//    for(std::list<Profile>::const_iterator it = resprofileslist.begin(); it != resprofileslist.end(); ++it) {
//        if(it->Invert(s,t)) {
//            if(sd > it->Evald(t) + softborder) {
//                return true;
//            }
//        }
//    }
//    return false;
}

bool FindLowestProfile(dReal s, Profile& profile, dReal& tres, const std::list<Profile>& resprofileslist){
    dReal t;
    dReal sdmin;
    dReal sdtmp;
    int i = 0;
    sdmin = INF;

    std::list<Profile>::const_iterator itbestprofile = resprofileslist.end();
    std::list<Profile>::const_iterator it = resprofileslist.begin();
    while(it != resprofileslist.end()) {
        if(it->Invert(s,t)) {
            sdtmp = it->Evald(t);
            if(sdtmp < sdmin) {
                sdmin = sdtmp;
                itbestprofile = it;
                tres = t;
            }
        }
        it++;
        i++;
    }
    if( itbestprofile != resprofileslist.end() ) {
        profile = *itbestprofile;
        return true;
    }
    return false;
}

ProfileSample FindLowestProfileFast(dReal scur, dReal sdmax, const std::list<Profile>& resprofileslist)
{
    ProfileSample bestprofile;
    bestprofile.itprofile = resprofileslist.end();
    bestprofile.sd = sdmax;
    for(std::list<Profile>::const_iterator itprofile = resprofileslist.begin(); itprofile != resprofileslist.end(); ++itprofile) {
        if( scur >= itprofile->svect.at(0) && scur <= itprofile->svect.back() ) {
            // find the sdcur
            dReal sdcur = 0, sddcur = 0, tcur=0;
            std::vector<dReal>::const_iterator its = std::lower_bound(itprofile->svect.begin(), itprofile->svect.end(), scur);
            size_t index = its - itprofile->svect.begin();
            if( index == 0 ) {
                sdcur = itprofile->sdvect.at(0);
                sddcur = itprofile->sddvect.at(0);
            }
            else {
                // have to interpolate quadratically
                // know that s + sd*t + 0.5*sdd*t*t = scur
                // solve for t: t = (-sd +- sqrt(sd*sd - 2*(s-scur)*sdd))/(sdd)
                // sdcur = sd + sdd*t = +- sqrt(sd*sd - 2*(s-scur)*sdd)
                dReal sd = itprofile->sdvect.at(index-1);
                sddcur = itprofile->sddvect.at(index-1);
                dReal discr = sd*sd - 2*sddcur*(itprofile->svect.at(index-1) - scur);
                if( discr < 0 ) {
                    if( discr > -TINY ) {
                        discr = 0;
                    }
                    else {
                        continue;
                    }
                }
                else if( discr > 0 ) {
                    discr = sqrt(discr);
                }

                // have to take the sign such that 0 <= t <= integrationtimestep
                // test lowest sd first
                if( sddcur > TINY ) {
                    tcur = (-sd + discr)/sddcur;
                    if( tcur < 0 || tcur > itprofile->integrationtimestep ) {
                        tcur =(-sd - discr)/sddcur;
                        if( tcur < 0 || tcur > itprofile->integrationtimestep ) {
                            // could not solve
                            continue;
                        }
                        sdcur = -discr;
                    }
                    else {
                        sdcur = discr;
                    }
                }
                else if( sddcur < -TINY ) {
                    tcur = (-sd - discr)/sddcur;
                    if( tcur < 0 || tcur > itprofile->integrationtimestep ) {
                        tcur =(-sd + discr)/sddcur;
                        if( tcur < 0 || tcur > itprofile->integrationtimestep ) {
                            // could not solve
                            continue;
                        }
                        sdcur = discr;
                    }
                    else {
                        sdcur = -discr;
                    }
                }
                else {
                    // sdd is very close to 0, so ignore
                    tcur = (scur - itprofile->svect.at(index-1))/sd;
                    if( tcur < 0 || tcur > itprofile->integrationtimestep ) {
                        continue;
                    }
                    sdcur = sd;
                }

                index -= 1;
            }
            if( sdcur < bestprofile.sd ) {
                bestprofile.s = scur;
                bestprofile.sd = sdcur;
                bestprofile.sdd = sddcur;
                bestprofile.sindex = index;
                bestprofile.t = tcur;
                bestprofile.itprofile = itprofile;
            }
        }
    }
    return bestprofile;
}

ProfileSample FindEarliestProfileIntersection(dReal sstart, dReal sdstart, dReal sddstart, dReal tmax, const std::list<Profile>& resprofileslist, std::list<Profile>::const_iterator itprofileexclude, dReal& tintersect)
{
    dReal smax = sstart + tmax*(sdstart + 0.5*tmax*sddstart);
    //dReal sdmax = sdstart + tmax*sddstart;

    ProfileSample bestprofile;
    bestprofile.itprofile = resprofileslist.end();
    for(std::list<Profile>::const_iterator itprofile = resprofileslist.begin(); itprofile != resprofileslist.end(); ++itprofile) {
        if( itprofile == itprofileexclude ) {
            continue;
        }
        if( smax >= itprofile->svect.at(0) && sstart < itprofile->svect.back()-TINY ) {
            std::vector<dReal>::const_iterator its0 = std::lower_bound(itprofile->svect.begin(), itprofile->svect.end(), sstart);
            size_t index0 = its0 - itprofile->svect.begin();
            if( index0 == 0 ) {
                // only reason this would happen is if s and itprofile->svect[0] are close
                if( fabs(sstart-itprofile->svect.at(0)) <= TINY && fabs(sdstart-itprofile->sdvect.at(0)) <= TINY ) {
                    bestprofile.s = itprofile->svect.at(0);
                    bestprofile.sd = itprofile->sdvect.at(0);
                    bestprofile.sdd = itprofile->sddvect.at(0);
                    bestprofile.sindex = 0;
                    bestprofile.t = 0;
                    bestprofile.itprofile = itprofile;
                    return bestprofile;
                }
            }
            else {
                --index0;
            }

            std::vector<dReal>::const_iterator its1 = std::lower_bound(itprofile->svect.begin(), itprofile->svect.end(), smax);
            size_t index1 = its1 - itprofile->svect.begin();

            // need to check all segments in [index0, index1) for intersection with [sstart, sdstart, sddstart]
            for(size_t index = index0; index < index1; ++index) {
                dReal snext = itprofile->svect.at(index);
                dReal sdnext = itprofile->sdvect.at(index);
                dReal sddnext = itprofile->sddvect.at(index);
                tintersect=INF; // time from sstart
                if( fabs(sddstart) <= TINY ) {
                    if( fabs(sddnext) <= TINY ) {
                        if( fabs(sdstart-sdnext) > TINY ) {
                            //RAVELOG_ERROR_FORMAT("sddstart and sddnext are both close to 0 at s=%.15e, sd diff=%.15e, don't know that to do", sstart%fabs(sdstart-sdnext));
                            continue;
                        }
                        else {
                            dReal t = (snext - sstart)/sdstart;
                            if( t >= 0 && t <= tmax ) {
                                bestprofile.s = snext;
                                bestprofile.sd = sdnext;
                                bestprofile.sdd = sddnext;
                                bestprofile.sindex = index;
                                bestprofile.t = 0;
                                bestprofile.itprofile = itprofile;
                                return bestprofile;
                            }
                            else {
                                continue;
                            }
                        }
                    }
                    bestprofile.t = (sdstart - sdnext)/sddnext;
                    bestprofile.s = snext + bestprofile.t * (sdnext + bestprofile.t*sddnext*0.5);
                    bestprofile.sd = sdstart;
                    bestprofile.sdd = sddnext;
                    tintersect = (bestprofile.s - sstart)/sdstart;
                }
                else if( fabs(sddnext) <= TINY ) {
                    tintersect = (sdnext - sdstart)/sddstart;
                    bestprofile.s = sstart + tintersect * (sdstart + tintersect*sddstart*0.5);
                    bestprofile.sd = sdnext;
                    bestprofile.sdd = sddnext;
                    bestprofile.t = (bestprofile.s - snext)/sdnext;
                }
                else if( fabs(sddstart-sddnext) <= TINY ) {
                    // two profiles will never intersect
                    continue;
                }
                else {
                    // s(sd) = sstart + (sd*sd - sdstart*sdstart)/(2*sddstart) => astart * sd*sd + cstart
                    dReal astart = 1/(2*sddstart), cstart = (sstart - sdstart*sdstart/(2*sddstart));
                    dReal anext = 1/(2*sddnext), cnext = (snext - sdnext*sdnext/(2*sddnext));
                    dReal ad = astart - anext;
                    dReal sd2 = (cnext-cstart)/ad;
                    if( sd2 < 0 ) {
                        // two profiles never intersect!
                        continue;
                    }
                    bestprofile.sd = sqrt(sd2);
                    bestprofile.s = astart*sd2 + cstart;
                    bestprofile.sdd = sddnext;
                    tintersect = (bestprofile.sd-sdstart)/sddstart;
                    bestprofile.t = (bestprofile.sd - sdnext)/sddnext;
                }

                if( tintersect >= 0 && tintersect <= tmax && bestprofile.t >= 0 && bestprofile.t <= itprofile->integrationtimestep+TINY && bestprofile.s > sstart ) {
                    if( bestprofile.t < 0 ) {
                        //std::cerr << "best profile is negative" << std::endl;
                        continue;
                    }
                    bestprofile.itprofile = itprofile;
                    bestprofile.sindex = index;
                    return bestprofile;
                }
            }
        }
    }

    // reached end, so didn't find any intersections
    bestprofile.itprofile = resprofileslist.end();
    return bestprofile;
}

void CheckInsert(std::list<std::pair<dReal,dReal> >&reslist, std::pair<dReal,dReal> e, bool inverse){
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


void FindMaxima(const std::list<std::pair<dReal,dReal> >&origlist, std::list<std::pair<dReal,dReal> >&reslist, bool inverse){
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
