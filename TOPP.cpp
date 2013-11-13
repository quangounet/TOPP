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


#include "TOPP.h"


namespace TOPP {


////////////////////////////////////////////////////////////////////
/////////////////////////// Tunings ////////////////////////////////
////////////////////////////////////////////////////////////////////


Tunings::Tunings(const std::string& tuningsstring) {
    std::istringstream iss(tuningsstring);
    std::string sub;

    iss >> sub;
    discrtimestep = atof(sub.c_str());

    iss >> sub;
    integrationtimestep = atof(sub.c_str());

    iss >> sub;
    reparamtimestep = atof(sub.c_str());

    iss >> sub;
    passswitchpointnsteps = atof(sub.c_str());

    bisectionprecision = 0.01;
    loweringcoef = 0.95;
}


////////////////////////////////////////////////////////////////////
/////////////////////////// Constraints ////////////////////////////
////////////////////////////////////////////////////////////////////


void Constraints::Preprocess(Trajectory& trajectory0, Tunings& tunings0) {
    trajectory = trajectory0;
    // Change discrtimestep so as it becomes a divisor of trajectory duration
    int ndiscrsteps = int((trajectory.duration+1e-10)/tunings0.discrtimestep);
    tunings0.discrtimestep = trajectory.duration/ndiscrsteps;
    //std::cout << tunings0.discrtimestep << "\n";
    tunings = tunings0;
    Discretize();
    ComputeMVCBobrow();
    ComputeMVCCombined();
    FindSwitchPoints();

    // Set the reparam timestep automatically if it is initially set to 0
    dReal meanmvc = 0;
    if(tunings.integrationtimestep == 0) {
        for(int i=0; i<ndiscrsteps; i++) {
            meanmvc += std::min(mvccombined[i],10.);
        }
        meanmvc /= ndiscrsteps;
        meanmvc = std::min(1.,meanmvc);
        tunings.integrationtimestep = tunings.discrtimestep/meanmvc;
        //std::cout << "\n--------------\nIntegration timestep: " << tunings.integrationtimestep << "\n";
    }
}


void Constraints::Discretize() {
    ndiscrsteps = int((trajectory.duration+TINY)/tunings.discrtimestep);
    ndiscrsteps++;
    discrsvect.resize(0);
    for(int i=0; i<ndiscrsteps; i++) {
        discrsvect.push_back(i*tunings.discrtimestep);
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
    int n = int(s/tunings.discrtimestep);
    dReal coef = (s-n*tunings.discrtimestep)/tunings.discrtimestep;
    return (1-coef)*v[n] + coef*v[n+1];
}


dReal Constraints::SdLimitCombinedInit(dReal s){
    dReal res = SdLimitBobrow(s);
    if(hasvelocitylimits) {
        std::vector<dReal> qd(trajectory.dimension);
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


void Constraints::WriteMVCCombined(std::stringstream& ss, dReal dt){
    dReal duration = trajectory.duration;
    ss << duration << " " << dt << "\n";
    for(dReal t=0; t<=duration; t+=dt) {
        ss << t << " ";
    }
    ss << "\n";
    for(dReal t=0; t<=duration; t+=dt) {
        ss << SdLimitCombined(t) << " ";
    }
}


void Constraints::FindSwitchPoints() {
    FindSingularSwitchPoints();
    FindTangentSwitchPoints();
    FindDiscontinuousSwitchPoints();
}


void Constraints::AddSwitchPoint(int i, int switchpointtype, dReal sd){
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
    // If no sd is specified, then take the value of the mvc
    // (The case when sd is specified corresponds to a singular switchpoint in some cases)
    if(sd<0) {
        sd = mvcbobrow[iadd];
    }
    if(sd > MAXSD) {
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
    if(ndiscrsteps<3)
        return;
    int i = 1;
    dReal s,sd,snext,sdnext,alpha,diff,diffprev;
    std::pair<dReal,dReal> sddlimits;

    s = discrsvect[i];
    snext = discrsvect[i+1];
    sd = SdLimitBobrow(s);
    sdnext = SdLimitBobrow(snext);
    sddlimits = SddLimits(s,sd);
    alpha = sddlimits.first;
    //beta = sddlimits.second;
    diffprev = alpha/sd - (sdnext-sd)/tunings.discrtimestep;

    for(int i=2; i<ndiscrsteps-1; i++) {
        s = discrsvect[i];
        snext = discrsvect[i+1];
        sd = SdLimitBobrow(s);
        sdnext = SdLimitBobrow(snext);
        sddlimits = SddLimits(s,sd);
        alpha = sddlimits.first;
        //beta = sddlimits.second;
        diff = alpha/sd - (sdnext-sd)/tunings.discrtimestep;
        if(diffprev*diff<0) {
            AddSwitchPoint(i,SP_TANGENT);
        }
        diffprev = diff;
    }
}


void Constraints::FindDiscontinuousSwitchPoints() {

}



////////////////////////////////////////////////////////////////////
/////////////////// Quadratic Constraints //////////////////////////
////////////////////////////////////////////////////////////////////


QuadraticConstraints::QuadraticConstraints(const std::string& constraintsstring) {
    int buffsize = BUFFSIZE;
    std::vector<dReal> tmpvect;
    char buff[buffsize];
    std::istringstream iss(constraintsstring);
    iss.getline(buff,buffsize);
    VectorFromString(std::string(buff),vmax);
    while(iss.good()) {
        iss.getline(buff,buffsize);
        VectorFromString(std::string(buff),tmpvect);
        avect.push_back(tmpvect);
        iss.getline(buff,buffsize);
        VectorFromString(std::string(buff),tmpvect);
        bvect.push_back(tmpvect);
        iss.getline(buff,buffsize);
        VectorFromString(std::string(buff),tmpvect);
        cvect.push_back(tmpvect);
    }
    nconstraints = int(avect.front().size());
    hasvelocitylimits =  VectorMax(vmax) > TINY;
    maxrep = 1;
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

    int n = int(s/tunings.discrtimestep);
    dReal coef = (s-n*tunings.discrtimestep)/tunings.discrtimestep;
    for (int i = 0; i < nconstraints; i++) {
        a[i] = (1-coef)*avect[n][i] + coef*avect[n+1][i];
        b[i] = (1-coef)*bvect[n][i] + coef*bvect[n+1][i];
        c[i] = (1-coef)*cvect[n][i] + coef*cvect[n+1][i];
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

    std::vector<std::pair<dReal,int> > vp;
    for(int i=0; i<trajectory.dimension; i++) {
        vp.push_back(std::pair<dReal,int>(std::abs(a[i]),i));
    }
    std::sort(vp.begin(),vp.end());
    slopesvector.resize(0);
    for(int i=0; i<trajectory.dimension; i++) {
        ap = (a[vp[i].second]-a2[vp[i].second])/delta;
        bp = (b[vp[i].second]-b2[vp[i].second])/delta;
        cp = (c[vp[i].second]-c2[vp[i].second])/delta;
        slope = (-bp*sd*sd-cp)/((2*b[vp[i].second]+ap)*sd);
        slopesvector.push_back(slope);
    }
}



std::pair<dReal,dReal> QuadraticConstraints::SddLimits(dReal s, dReal sd){
    dReal dtsq = tunings.integrationtimestep;
    dtsq = dtsq*dtsq;
    dReal safetybound = tunings.discrtimestep/dtsq;
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
    std::pair<dReal,dReal> sddlimits = SddLimits(s,0);
    if(sddlimits.first > sddlimits.second) {
        return 0;
    }
    std::vector<dReal> a, b, c;
    InterpolateDynamics(s,a,b,c);

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
        for(int j=0; j<trajectory.dimension; j++) {
            if(a[j]*aprev[j]<0) {
                AddSwitchPoint(i,SP_SINGULAR);
                continue;
            }
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
    ret = IntegrateBackward(constraints,s,sd,dt,resprofile,constraints.tunings.passswitchpointnsteps);
    if(ret==INT_MAXSTEPS||ret==INT_END) {
        ret = IntegrateForward(constraints,s,sd,dt,resprofile,constraints.tunings.passswitchpointnsteps,testaboveexistingprofiles, testmvc, zlajpah);
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
    if(sdtop-sdbottom<constraints.tunings.bisectionprecision) {
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
    dReal discr = constraints.tunings.discrtimestep;
    dReal dt = 0; // should be changed
    Profile resprofile;

    // Tangent, Discontinuous and Zlajpah switchpoints
    if(switchpoint.switchpointtype == SP_TANGENT || switchpoint.switchpointtype == SP_DISCONTINUOUS || switchpoint.switchpointtype == SP_ZLAJPAH) {
        dt = discr/2;
        dReal sdtop = BisectionSearch(constraints,s,sd*constraints.tunings.loweringcoef,sd,dt,0);
        if(sdtop<=0) {
            return false;
        }
        sbackward = s;
        sforward = s;
        sdbackward = sdtop;
        sdforward = sdtop;
        return true;
    }
    // Singular switchpoints
    else{
        dReal sstep = constraints.tunings.passswitchpointnsteps*constraints.tunings.integrationtimestep;
        std::vector<dReal> slopesvector;
        dReal bestslope = 0;
        dReal bestscore = INF;
        constraints.ComputeSlopeDynamicSingularity(s,sd,slopesvector);
        sforward = s + sstep;
        sbackward = s - sstep;
        if(sforward>constraints.trajectory.duration || sbackward<0) {
            return false;
        }
        for(int i = 0; i<int(slopesvector.size()); i++) {
            dReal slope = slopesvector[i];
            sdforward = sd + sstep*slope;
            sdbackward = sd - sstep*slope;
            dReal alphabackward = constraints.SddLimits(sbackward,sdbackward).first/sd;
            dReal betaforward = constraints.SddLimits(sforward,sdforward).second/sd;
            dReal bob1 = constraints.SdLimitBobrow(sbackward)-sdbackward;
            dReal bob2 = constraints.SdLimitBobrow(sforward)-sdforward;
            dReal slopediff1 = std::abs(alphabackward-slope);
            dReal slopediff2 = std::abs(betaforward-slope);
            //std::cout << s << " " << bob1 << " " << bob2  << " " << slopediff1 << " " << slopediff2 << "\n";
            if(bob1>=0 && bob2>=0) {
                dReal score = slopediff1+slopediff2;
                if(score<bestscore) {
                    bestslope = slope;
                    bestscore = score;
                }
            }
        }
        sdforward = sd + sstep*bestslope;
        sdbackward = sd - sstep*bestslope;
        return bestscore<INF;

        // here switchpointtype == SP_SINGULAR
        // dt = discr/2;
        // ret = IntegrateBackward(constraints,s-discr,sd*constraints.tunings.loweringcoef,dt,resprofile,constraints.tunings.passswitchpointnsteps,testaboveexistingprofiles,testmvc);
        // if(ret==INT_MAXSTEPS || ret==INT_END || ret==INT_PROFILE) {
        //     sbackward = resprofile.Eval(0);
        //     sdbackward = resprofile.Evald(0);
        //     if(sdbackward>constraints.SdLimitBobrow(sbackward)) {
        //         return false;
        //     }
        //     ret = IntegrateForward(constraints,s+discr,sd*constraints.tunings.loweringcoef,dt,resprofile,constraints.tunings.passswitchpointnsteps,testaboveexistingprofiles,testmvc,zlajpah);
        //     if(ret==INT_MAXSTEPS || ret==INT_END || ret==INT_PROFILE) {
        //         sforward = resprofile.Eval(resprofile.duration);
        //         sdforward = resprofile.Evald(resprofile.duration);
        //         if(sdforward>constraints.SdLimitBobrow(sforward)) {
        //             return false;
        //         }
        //         return true;
        //     }
        // }
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
    while(beta-alpha>constraints.tunings.bisectionprecision) {
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
    while(beta-alpha>constraints.tunings.bisectionprecision) {
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
                int res3 = IntegrateBackward(constraints,scur,sdcur,constraints.tunings.integrationtimestep,tmpprofile,1e5,true,true,true);
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
        else if((!zlajpah) && testmvc && sdcur > constraints.SdLimitCombined(scur)) {
            slist.push_back(scur);
            sdlist.push_back(sdcur);
            constraints.zlajpahlist.push_back(std::pair<dReal,dReal>(scur,sdcur));
            returntype = INT_MVC;
            break;
        }
        else{
            slist.push_back(scur);
            sdlist.push_back(sdcur);
            std ::pair<dReal,dReal> sddlimits = constraints.SddLimits(scur,sdcur);
            dReal beta = sddlimits.second;
            sddlist.push_back(beta);
            //std::cout << scur << "," << sdcur << "," << beta << "\n";
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
            // // Lower the sd to the MVC
            // if(slist.size()>0) {
            //     dReal sprev = slist.back(), sdprev = sdlist.back();
            //     dReal slidesdd = ComputeSlidesddBackward(constraints,sprev,sdprev,-dt);
            //     scur = sprev - sdprev*dt + 0.5*dtsq*slidesdd;
            //     sdcur = sdprev - dt*slidesdd;
            //     sddlist.pop_back();
            //     sddlist.push_back(slidesdd);
            // }
            //Now we have sdcombined < sdcur <= sdbobrow
            // Slide along MVCCombined, which is admissible, until either
            // alpha points above MVCCombined (trapped)

            //std::cout <<"\nSlide backward ("<< scur << "," << sdcur << ") \n";
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
        else if(!zlajpah && testmvc && sdcur > constraints.SdLimitCombined(scur)) {
            slist.push_back(scur);
            sdlist.push_back(sdcur);
            returntype = INT_MVC;
            break;
        }
        else{
            slist.push_back(scur);
            sdlist.push_back(sdcur);
            std ::pair<dReal,dReal> sddlimits = constraints.SddLimits(scur,sdcur);
            dReal alpha = sddlimits.first;
            //std::cout << scur << " " << sdcur << " " << alpha << "\n";
            sddlist.push_back(alpha);
            dReal sprev = scur - dt * sdcur + 0.5*dtsq*alpha;
            dReal sdprev = sdcur - dt * alpha;
            scur = sprev;
            sdcur = sdprev;
        }
    }
    if(sddlist.size()>0) {
        slist.reverse();
        sdlist.reverse();
        sddlist.reverse();
        sddlist.pop_front();
        sddlist.push_back(sddlist.back());
    }
    resprofile = Profile(slist,sdlist,sddlist,dt);
    resprofile.forward = false;
    return returntype;
}


////////////////////////////////////////////////////////////////////
//////////////////////// Limiting curves ///////////////////////////
////////////////////////////////////////////////////////////////////


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
        if(IsAboveProfilesList(sswitch,sdswitch,constraints.resprofileslist))
            continue;
        if(sdswitch > constraints.SdLimitCombined(sswitch)+TINY2)
            continue;

        // Address Switch Point
        if (!AddressSwitchPoint(constraints, switchpoint, sbackward,
                                sdbackward, sforward, sdforward))
            continue;

        // Add middle part
        if (sforward - sbackward > TINY) {
            std::list<dReal> slist, sdlist, sddlist;
            dReal dtmod = 2 * (sforward - sbackward) / (sdforward + sdbackward);
            dReal sdd = (sdforward - sdbackward) / dtmod;
            slist.push_back(sbackward);
            slist.push_back(sforward);
            sdlist.push_back(sdbackward);
            sdlist.push_back(sdforward);
            sddlist.push_back(sdd);
            sddlist.push_back(0);
            constraints.resprofileslist.push_back(Profile(slist,sdlist,sddlist,dtmod));
        }

        // Integrate backward
        integratestatus = IntegrateBackward(constraints, sbackward, sdbackward,
                                            constraints.tunings.integrationtimestep, tmpprofile);
        if(tmpprofile.nsteps>1)
            constraints.resprofileslist.push_back(tmpprofile);
        if(integratestatus == INT_BOTTOM)
            return CLC_BOTTOM;

        // Integrate forward
        integratestatus = IntegrateForward(constraints, sforward, sdforward,
                                           constraints.tunings.integrationtimestep, tmpprofile, 1e5,
                                           testaboveexistingprofiles, testmvc, zlajpah);
        if(tmpprofile.nsteps>1)
            constraints.resprofileslist.push_back(tmpprofile);
        if(integratestatus == INT_BOTTOM)
            return CLC_BOTTOM;
    }
    return CLC_OK;
}


int ComputeProfiles(Constraints& constraints, Trajectory& trajectory, Tunings& tunings, dReal sdbeg, dReal sdend){
    std::chrono::time_point<std::chrono::system_clock> t0,t1,t2,t3;
    std::chrono::duration<double> d1,d2,d3,dtot;

    t0 = std::chrono::system_clock::now();

    constraints.Preprocess(trajectory,tunings);
    if(VectorMin(constraints.mvcbobrow) <= TINY) {
        //std::cout << "MVCBobrow hit 0\n";
        return 0;
    }
    Profile resprofile;
    int ret;

    t1 = std::chrono::system_clock::now();

    bool integrateprofilesstatus = true;

    std::string message;
    for(int rep=0; rep<constraints.maxrep; rep++) {

        // Lower integrationtimestep
        if(rep>0) {
            constraints.tunings.integrationtimestep /= 3.3;
            //std::cout << "!!!!!!!!!!!!!!!!!!!!!!!!!! Try lower integration timestep: " << constraints.tunings.integrationtimestep << "!!!!!!!!!!!!!!!!!!!!!!!!\n";
        }
        constraints.resprofileslist.resize(0);
        constraints.zlajpahlist.resize(0);

        // Compute the CLC
        ret = ComputeLimitingCurves(constraints);
        if(ret!=CLC_OK) {
            message = "CLC failed\n";
            integrateprofilesstatus = false;
            continue;
        }

        // flags
        bool testaboveexistingprofiles = true, testmvc = true, zlajpah = false;

        // Integrate forward from s = 0
        ret = IntegrateForward(constraints,0,sdbeg,constraints.tunings.integrationtimestep,resprofile,1e5,testaboveexistingprofiles,testmvc,zlajpah);
        if(resprofile.nsteps>1) {
            constraints.resprofileslist.push_back(resprofile);
        }
        if(ret==INT_BOTTOM) {
            message = "FW reached 0\n";
            integrateprofilesstatus = false;
            continue;
        }

        // Integrate backward from s = s_end
        ret = IntegrateBackward(constraints,trajectory.duration,sdend,constraints.tunings.integrationtimestep,resprofile,1e5,testaboveexistingprofiles,testmvc);
        if(resprofile.nsteps>1) {
            constraints.resprofileslist.push_back(resprofile);
        }
        if(ret==INT_BOTTOM) {
            message = "BW reached 0\n";
            integrateprofilesstatus = false;
            continue;
        }

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
            ret = IntegrateForward(constraints,zs,zsd,constraints.tunings.integrationtimestep,resprofile,1e5,testaboveexistingprofiles,testmvc,zlajpah);
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

        // Reset currentindex
        std::list<Profile>::iterator it = constraints.resprofileslist.begin();
        while(it != constraints.resprofileslist.end()) {
            it->currentindex = 0;
            it++;
        }

        // Estimate resulting trajectory duration
        constraints.resduration = 0;
        Profile profile;
        dReal ds = tunings.discrtimestep/10;
        dReal sdprev,sdnext;
        dReal tres;
        if(FindLowestProfile(0,profile,tres,constraints.resprofileslist)) {
            sdprev = profile.Evald(tres);
        }
        else{
            message = "CLC discontinuous";
            integrateprofilesstatus = false;
            continue;
        }
        bool clcdiscontinuous = false;
        for(dReal s=ds; s<=trajectory.duration; s+=ds) {
            if(FindLowestProfile(s,profile,tres,constraints.resprofileslist)) {
                sdnext = profile.Evald(tres);
                constraints.resduration += 2*ds/(sdprev+sdnext);
                sdprev = sdnext;
            }
            else{
                clcdiscontinuous = true;
                break;
            }
        }
        if(clcdiscontinuous) {
            message = "CLC discontinuous";
            integrateprofilesstatus = false;
            continue;
        }

        integrateprofilesstatus = true;
        break;

    }


    if(!integrateprofilesstatus) {
        std::cout << message << "\n";
        return 0;
    }

    t2 = std::chrono::system_clock::now();

    // d1 = t1-t0;
    // d2 = t2-t1;
    // d3 = t3-t2;
    // dtot = t3-t0;

    //std::cout << "Constraints preprocessing: " << d1.count() << "\n";
    //std::cout << "Profiles calculation: " <<  d2.count() << "\n";
    //std::cout << "Duration calculation: " <<  d3.count() << "\n";
    //std::cout << "Total PP: " <<  dtot.count() << "\n";

    return 1;
}


int VIP(Constraints& constraints, Trajectory& trajectory, Tunings& tunings, dReal sdbegmin, dReal sdbegmax, dReal& sdendmin, dReal& sdendmax){
    if (trajectory.duration <= 0) {
        std::cout << "[TOPP] Warning: trajectory duration is <= 0\n";
        return 0;
    }

    constraints.Preprocess(trajectory,tunings);
    if(VectorMin(constraints.mvcbobrow) <= TINY) {
        //std::cout << "MVCBobrow hit 0\n";
        return 0;
    }

    Profile tmpprofile;
    dReal tres;
    dReal smallincrement = constraints.tunings.integrationtimestep*2;

    // Compute the limiting curves
    int resclc = ComputeLimitingCurves(constraints);
    if(resclc == CLC_SWITCH || resclc == CLC_BOTTOM) {
        return 0;
    }
    // Determine the lowest profile at t=0
    dReal bound;
    if(FindLowestProfile(smallincrement,tmpprofile,tres,constraints.resprofileslist)) {
        bound = std::min(tmpprofile.Evald(tres),constraints.mvccombined[0]);
    }
    // just to make sure the profile is below mvccombined
    else{
        bound = constraints.mvccombined[0];
    }

    if(sdbegmin>bound) {
        return 0;
    }

    sdbegmax = std::min(sdbegmax,bound-tunings.bisectionprecision);

    bool testaboveexistingprofiles = true, testmvc = true, zlajpah = true;

    // Compute sdendmax by integrating forward from (0,sdbegmax)
    int resintfw = IntegrateForward(constraints,0,sdbegmax,constraints.tunings.integrationtimestep,tmpprofile,1e5,testaboveexistingprofiles,testmvc,zlajpah);
    constraints.resprofileslist.push_back(tmpprofile);
    if(resintfw == INT_BOTTOM || resintfw == INT_MVC)
        return 0;
    else if (resintfw == INT_END && tmpprofile.Eval(tmpprofile.duration)<= constraints.mvccombined[constraints.mvccombined.size()-1]) {
        sdendmax = tmpprofile.Evald(tmpprofile.duration);
    }
    else if (resintfw == INT_PROFILE) {
        // Look for the lowest profile at the end
        if(FindLowestProfile(trajectory.duration-smallincrement,tmpprofile,tres,constraints.resprofileslist) && tmpprofile.Evald(tres) <= constraints.mvccombined[constraints.mvccombined.size()-1]) {
            sdendmax = tmpprofile.Evald(tres);
        }
        else {
            // No profile reaches the end, consider the MVC instead
            sdendmax = constraints.mvccombined[constraints.mvccombined.size()-1];
            int count = 0;
            int resintbw;
            dReal dtint = constraints.tunings.integrationtimestep;
            // If integrating from sdendmax fails with INT_BOTTOM or INT_MVC, then the trajectory is not traversable. However, since integrating backward from a high sdendmax can be risky, we give three chances by decreasing the value of the integration step
            while(count<3) {
                count++;
                resintbw = IntegrateBackward(constraints,trajectory.duration,sdendmax,dtint,tmpprofile,1e5);
                if(resintbw != INT_BOTTOM && resintbw != INT_MVC) {
                    break;
                }
                dtint /= 3.3;
            }
            if(resintbw == INT_BOTTOM || resintbw == INT_MVC)
                return 0;
        }
    }

    // Integrate from (send,0). If succeeds, then sdendmin=0 and exits
    int resintbw = IntegrateBackward(constraints,trajectory.duration,0,constraints.tunings.integrationtimestep,tmpprofile,1e5);
    if((resintbw == INT_END && tmpprofile.Evald(0)>=sdbegmin) || resintbw == INT_PROFILE) {
        constraints.resprofileslist.push_back(tmpprofile);
        sdendmin = 0;
        return 1;
    }

    // Determine sdendmin by bisection
    dReal sdupper = sdendmax, sdlower = 0;
    Profile bestprofile;
    while(sdupper-sdlower > tunings.bisectionprecision) {
        dReal sdtest = (sdupper + sdlower)/2;
        int resintbw2 = IntegrateBackward(constraints,trajectory.duration,sdtest,constraints.tunings.integrationtimestep,tmpprofile,1e5);
        if((resintbw2 == INT_END && tmpprofile.Evald(0)>=sdbegmin) || resintbw2 == INT_PROFILE) {
            sdupper = sdtest;
            bestprofile = tmpprofile;
        }
        else
            sdlower = sdtest;
    }
    sdendmin = sdupper;
    constraints.resprofileslist.push_back(bestprofile);

    return 1;
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


dReal VectorMin(const std::vector<dReal>&v){
    std::vector<dReal>::const_iterator it = v.begin();
    dReal res = INF;
    while(it!=v.end()) {
        res = std::min(res,*it);
        it++;
    }
    return res;
}


dReal VectorMax(const std::vector<dReal>&v){
    std::vector<dReal>::const_iterator it = v.begin();
    dReal res = -INF;
    while(it!=v.end()) {
        res = std::max(res,*it);
        it++;
    }
    return res;
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


bool IsAboveProfilesList(dReal s, dReal sd, std::list<Profile>&resprofileslist, bool searchbackward){
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
            if(sd > it->Evald(t) + TINY2) {
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

} // end namespace TOPP
