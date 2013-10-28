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



#include "KinematicLimits.h"


namespace TOPP {

KinematicLimits::KinematicLimits(const std::string& constraintsstring){
    int buffsize = BUFFSIZE;
    char buff[buffsize];
    std::istringstream iss(constraintsstring);
    iss.getline(buff,buffsize);
    VectorFromString(std::string(buff),amax);
    iss.getline(buff,buffsize);
    VectorFromString(std::string(buff),vmax);
    hasvelocitylimits =  VectorMax(vmax) > TINY;
}


void KinematicLimits::Discretize(){
    dReal s;
    std::vector<dReal> qd(trajectory.dimension), qdd(trajectory.dimension);
    ndiscrsteps = int((trajectory.duration+TINY)/tunings.discrtimestep);
    ndiscrsteps++;
    discrsvect.resize(0);
    for(int i=0; i<ndiscrsteps; i++) {
        s = i*tunings.discrtimestep;
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
    if(s>=trajectory.duration) {
        int n = ndiscrsteps-1;
        for(int i=0; i<trajectory.dimension; i++) {
            qd[i]= qdvect[n][i];
            qdd[i]= qddvect[n][i];
        }
        return;
    }
    int n = int(s/tunings.discrtimestep);
    dReal coef = (s-n*tunings.discrtimestep)/tunings.discrtimestep;
    for(int i=0; i<trajectory.dimension; i++) {
        qd[i] = (1-coef)*qdvect[n][i] + coef*qdvect[n+1][i];
        qdd[i] = (1-coef)*qddvect[n][i] + coef*qddvect[n+1][i];
    }
}


void KinematicLimits::ComputeSlopeDynamicSingularity(dReal s, dReal sd, std::vector<dReal>& slopesvector) {
    dReal delta = 0.001, s2, qdp, qddp, slope;
    std::vector<dReal> qd, qdd, qd2, qdd2;
    if(s>delta) {
        delta = -delta;
    }
    s2 = s + delta;
    InterpolateDynamics(s,qd,qdd);
    InterpolateDynamics(s2,qd2,qdd2);

    std::vector<std::pair<dReal,int> > vp;
    for(int i=0; i<trajectory.dimension; i++) {
        vp.push_back(std::pair<dReal,int>(std::abs(qd[i]),i));
    }
    std::sort(vp.begin(),vp.end());
    slopesvector.resize(0);
    for(int i=0; i<trajectory.dimension; i++) {
        qdp = (qd2[vp[i].second]-qd[vp[i].second])/delta;
        qddp = (qdd2[vp[i].second]-qdd[vp[i].second])/delta;
        slope = (-qddp*sd*sd)/((2*qdd[vp[i].second]+qdp)*sd);
        //std::cout << vp[i].second << " " << slope << "***\n";
        slopesvector.push_back(slope);
    }
}


std::pair<dReal,dReal> KinematicLimits::SddLimits(dReal s, dReal sd){
    dReal dtsq = tunings.integrationtimestep;
    dtsq = dtsq*dtsq;
    dReal safetybound = tunings.discrtimestep/dtsq;
    dReal alpha = -safetybound;
    dReal beta = safetybound;
    dReal sdsq = sd*sd;
    dReal a_alpha_i, a_beta_i, alpha_i, beta_i;
    std::vector<dReal> qd, qdd;
    InterpolateDynamics(s,qd,qdd);
    for(int i=0; i<trajectory.dimension; i++) {
        if(std::abs(qd[i])<=TINY) {
            continue;
        }
        if(qd[i]>0) {
            a_alpha_i = -amax[i];
            a_beta_i = amax[i];
        }
        else{
            a_alpha_i = amax[i];
            a_beta_i = -amax[i];
        }
        alpha_i = (a_alpha_i-sdsq*qdd[i])/qd[i];
        beta_i = (a_beta_i-sdsq*qdd[i])/qd[i];
        alpha = std::max(alpha,alpha_i);
        beta = std::min(beta,beta_i);
    }
    std::pair<dReal,dReal> res(alpha,beta);
    return res;
}


dReal KinematicLimits::SdLimitBobrowInit(dReal s){
    std::pair<dReal,dReal> sddlimits = KinematicLimits::SddLimits(s,0);
    if(sddlimits.first > sddlimits.second) {
        return 0;
    }
    std::vector<dReal> a_alpha(trajectory.dimension), a_beta(trajectory.dimension);
    std::vector<dReal> qd(trajectory.dimension), qdd(trajectory.dimension);
    trajectory.Evald(s, qd);
    trajectory.Evaldd(s, qdd);
    for(int i=0; i<trajectory.dimension; i++) {
        if(qd[i] > 0) {
            a_alpha[i] = -amax[i];
            a_beta[i] = amax[i];
        }
        else{
            a_alpha[i] = amax[i];
            a_beta[i] = -amax[i];
        }
    }
    dReal sdmin = INF;
    for(int k=0; k<trajectory.dimension; k++) {
        for(int m=k+1; m<trajectory.dimension; m++) {
            dReal num, denum, r;
            num = qd[m]*a_alpha[k]-qd[k]*a_beta[m];
            denum = qd[m]*qdd[k]-qd[k]*qdd[m];
            if(std::abs(denum) > TINY) {
                r = num/denum;
                if(r>=0) {
                    sdmin = std::min(sdmin,sqrt(r));
                }
            }
            num = qd[k]*a_alpha[m]-qd[m]*a_beta[k];
            denum = -denum;
            if(std::abs(denum) > TINY) {
                r = num/denum;
                if(r>=0) {
                    sdmin = std::min(sdmin,sqrt(r));
                }
            }
        }
    }
    return sdmin;
}


void KinematicLimits::FindSingularSwitchPoints(){
    if(ndiscrsteps<3) {
        return;
    }
    int i = 0;
    std::vector<dReal> qd(trajectory.dimension),qdprev(trajectory.dimension);
    trajectory.Evald(discrsvect[i],qdprev);

    for(int i=1; i<ndiscrsteps-1; i++) {
        trajectory.Evald(discrsvect[i],qd);
        for(int j=0; j<trajectory.dimension; j++) {
            if(qd[j]*qdprev[j]<0) {
                AddSwitchPoint(i,SP_SINGULAR);
                continue;
            }
        }
        qdprev = qd;
    }
}
}
