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



#include "KinematicLimits.h"


namespace TOPP {

KinematicLimits::KinematicLimits(const std::string& constraintsstring){
    int buffsize = BUFFSIZE;
    char buff[buffsize];
    std::istringstream iss(constraintsstring);
    iss.getline(buff,buffsize);
    discrtimestep = atof(buff);
    iss.getline(buff,buffsize);
    VectorFromString(std::string(buff),vmax);
    iss.getline(buff,buffsize);
    VectorFromString(std::string(buff),amax);
    hasvelocitylimits =  VectorMax(vmax) > TINY;
}


void KinematicLimits::Discretize(){
    dReal s;
    std::vector<dReal> qd(trajectory.dimension), qdd(trajectory.dimension);
    ndiscrsteps = int((trajectory.duration+TINY)/discrtimestep);
    ndiscrsteps++;
    discrsvect.resize(0);
    for(int i=0; i<ndiscrsteps; i++) {
        s = i*discrtimestep;
        trajectory.Evald(s,qd);
        trajectory.Evaldd(s,qdd);
        discrsvect.push_back(s);
        qdvect.push_back(qd);
        qddvect.push_back(qdd);
    }
}


void KinematicLimits::InterpolateDynamics(dReal s, std::vector<dReal>& qd, std::vector<dReal>& qdd){
    qd.resize(trajectory.dimension);
    qdd.resize(trajectory.dimension);
    assert(s>=-TINY && s<=trajectory.duration+TINY);
    if(s<0) {
        s=0;
    }
    if(s>=trajectory.duration-TINY) {
        int n = ndiscrsteps-1;
        for(int i=0; i<trajectory.dimension; i++) {
            qd[i]= qdvect[n][i];
            qdd[i]= qddvect[n][i];
        }
        return;
    }
    int n = int(s/discrtimestep);
    dReal coef = (s-n*discrtimestep)/discrtimestep;
    for(int i=0; i<trajectory.dimension; i++) {
        qd[i] = (1-coef)*qdvect[n][i] + coef*qdvect[n+1][i];
        qdd[i] = (1-coef)*qddvect[n][i] + coef*qddvect[n+1][i];
    }
}


void KinematicLimits::ComputeSlopeDynamicSingularity(dReal s, dReal sd, std::vector<dReal>& slopesvector) {
    dReal delta = TINY2, s2, qdp, qddp, slope;
    std::vector<dReal> qd, qdd, qd2, qdd2;
    if(s>delta) {
        delta = -delta;
    }
    s2 = s + delta;
    InterpolateDynamics(s,qd,qdd);
    InterpolateDynamics(s2,qd2,qdd2);

    slopesvector.resize(0);
    for(int i=0; i<trajectory.dimension; i++) {
        qdp = (qd2[i]-qd[i])/delta;
        qddp = (qdd2[i]-qdd[i])/delta;
        slope = (-qddp*sd*sd)/((2*qdd[i]+qdp)*sd);
        //std::cout << i << " " << slope << "***\n";
        slopesvector.push_back(slope);
    }
}


std::pair<dReal,dReal> KinematicLimits::SddLimits(dReal s, dReal sd){
    dReal dtsq = integrationtimestep;
    dtsq = dtsq*dtsq;
    dReal safetybound = discrtimestep/dtsq;
    dReal alpha = -safetybound;
    dReal beta = safetybound;
    dReal sdsq = sd*sd;
    dReal alpha_i, beta_i;
    std::vector<dReal> qd, qdd;
    InterpolateDynamics(s,qd,qdd);
    for(int i=0; i<trajectory.dimension; i++) {
        if(std::abs(qd[i])<=TINY) {
            continue;
        }
        if(qd[i]>0) {
            alpha_i = (-amax[i]-sdsq*qdd[i])/qd[i];
            beta_i = (amax[i]-sdsq*qdd[i])/qd[i];
            alpha = std::max(alpha,alpha_i);
            beta = std::min(beta,beta_i);
        }
        else{
            alpha_i = (amax[i]-sdsq*qdd[i])/qd[i];
            beta_i = (-amax[i]-sdsq*qdd[i])/qd[i];
            alpha = std::max(alpha,alpha_i);
            beta = std::min(beta,beta_i);
        }
    }
    std::pair<dReal,dReal> res(alpha,beta);
    return res;
}


dReal KinematicLimits::SddLimitAlpha(dReal s, dReal sd){
    dReal dtsq = integrationtimestep;
    dtsq = dtsq*dtsq;
    dReal safetybound = discrtimestep/dtsq;
    dReal alpha = -safetybound;
    dReal sdsq = sd*sd;
    dReal alpha_i;
    std::vector<dReal> qd, qdd;
    InterpolateDynamics(s,qd,qdd);
    for(int i=0; i<trajectory.dimension; i++) {
        if(std::abs(qd[i])<=TINY) {
            continue;
        }
        if(qd[i]>0) {
            alpha_i = (-amax[i]-sdsq*qdd[i])/qd[i];
            alpha = std::max(alpha,alpha_i);
        }
        else{
            alpha_i = (amax[i]-sdsq*qdd[i])/qd[i];
            alpha = std::max(alpha,alpha_i);
        }
    }
    return alpha;
}

dReal KinematicLimits::SddLimitBeta(dReal s, dReal sd){
    dReal dtsq = integrationtimestep;
    dtsq = dtsq*dtsq;
    dReal safetybound = discrtimestep/dtsq;
    dReal beta = safetybound;
    dReal sdsq = sd*sd;
    dReal beta_i;
    std::vector<dReal> qd, qdd;
    InterpolateDynamics(s,qd,qdd);
    for(int i=0; i<trajectory.dimension; i++) {
        if(std::abs(qd[i])<=TINY) {
            continue;
        }
        if(qd[i]>0) {
            beta_i = (amax[i]-sdsq*qdd[i])/qd[i];
            beta = std::min(beta,beta_i);
        }
        else{
            beta_i = (-amax[i]-sdsq*qdd[i])/qd[i];
            beta = std::min(beta,beta_i);
        }
    }
    return beta;
}


dReal KinematicLimits::SdLimitBobrowInit(dReal s){
    std::pair<dReal,dReal> sddlimits = KinematicLimits::SddLimits(s,0);
    if(sddlimits.first > sddlimits.second) {
        return 0;
    }
    std::list<std::pair<dReal,dReal> > alpha_list, beta_list;
    std::vector<dReal> qd(trajectory.dimension), qdd(trajectory.dimension);
    //std::chrono::time_point<std::chrono::system_clock> t0,t1,t2,t3;
    //std::chrono::duration<double> d1,d2,d3;

    dReal a,b;

    //t0 = std::chrono::system_clock::now();
    trajectory.Evald(s, qd);
    trajectory.Evaldd(s, qdd);
    //t1 = std::chrono::system_clock::now();
    for(int i=0; i<trajectory.dimension; i++) {
        a = amax[i]/qd[i];
        b = -qdd[i]/qd[i];
        if(qd[i] > 0) {
            CheckInsert(alpha_list,std::pair<dReal,dReal>(-a,b));
            CheckInsert(beta_list,std::pair<dReal,dReal>(a,b),true);
        }
        else{
            CheckInsert(alpha_list,std::pair<dReal,dReal>(a,b));
            CheckInsert(beta_list,std::pair<dReal,dReal>(-a,b),true);
        }
    }
    //t2 = std::chrono::system_clock::now();
    //std::cout << alpha_list.size() << " " << beta_list.size() << "\n";
    dReal sdmin = INF;
    std::list<std::pair<dReal,dReal> >::iterator italpha = alpha_list.begin();
    while(italpha != alpha_list.end()) {
        std::list<std::pair<dReal,dReal> >::iterator itbeta = beta_list.begin();
        while(itbeta != beta_list.end()) {
            dReal num, denum, r;
            num = itbeta->first - italpha->first;
            denum = italpha->second - itbeta->second;
            if(std::abs(denum) > TINY) {
                r = num/denum;
                if(r>=0) {
                    sdmin = std::min(sdmin,sqrt(r));
                }
            }
            itbeta++;
        }
        italpha++;
    }
    //t3 = std::chrono::system_clock::now();

    //d1 = t1-t0;
    //d2 = t2-t1;
    //d3 = t3-t2;
    //std::cout << d1.count() <<  " " << d2.count() << " " << d3.count() << "\n";



    return sdmin;
}


void KinematicLimits::FindSingularSwitchPoints(){
    if(ndiscrsteps<3) {
        return;
    }
    int i = 0;
    std::vector<dReal> qd(trajectory.dimension),qdprev(trajectory.dimension),qdd(trajectory.dimension);
    std::vector<dReal> a_alpha(trajectory.dimension), a_beta(trajectory.dimension);
    trajectory.Evald(discrsvect[i],qdprev);

    for(int i=1; i<ndiscrsteps-1; i++) {
        trajectory.Evald(discrsvect[i],qd);
        trajectory.Evaldd(discrsvect[i],qdd);
        for(int j=0; j<trajectory.dimension; j++) {
            if(qd[j] > 0) {
                a_alpha[j] = -amax[j];
                a_beta[j] = amax[j];
            }
            else{
                a_alpha[j] = amax[j];
                a_beta[j] = -amax[j];
            }
        }
        dReal minsd = mvcbobrow[i];
        bool found = false;
        for(int j=0; j<trajectory.dimension; j++) {
            if(qd[j]*qdprev[j]<0) {
                found = true;
                if(a_beta[j]/qdd[j]<0) {
                    minsd = std::min(minsd,sqrt(-a_beta[j]/qdd[j]));
                }
                else{
                    minsd = std::min(minsd,sqrt(a_beta[j]/qdd[j]));
                }
            }
        }
        if(found) {
            AddSwitchPoint(i,SP_SINGULAR,minsd);
        }
        qdprev = qd;
    }
}
}
