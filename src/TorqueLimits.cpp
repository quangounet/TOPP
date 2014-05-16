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


#include "TorqueLimits.h"

namespace TOPP {

TorqueLimits::TorqueLimits(const std::string& constraintsstring){
    std::vector<dReal> tmpvect;
    std::string buff;
    std::istringstream iss(constraintsstring);
    getline(iss, buff,'\n');
    discrtimestep = atof(buff.c_str());
    getline(iss, buff,'\n');
    VectorFromString(buff,taumin);
    getline(iss, buff, '\n');
    VectorFromString(buff,taumax);
    getline(iss, buff, '\n');
    VectorFromString(buff, vmax);
    while(iss.good()) {
        getline(iss, buff,'\n');
        VectorFromString(buff,tmpvect);
        avect.push_back(tmpvect);
        getline(iss, buff,'\n');
        VectorFromString(buff,tmpvect);
        bvect.push_back(tmpvect);
        getline(iss, buff,'\n');
        VectorFromString(buff,tmpvect);
        cvect.push_back(tmpvect);
    }
    hasvelocitylimits = VectorMax(vmax) > TINY;
}


void TorqueLimits::InterpolateDynamics(dReal s, std::vector<dReal>& a, std::vector<dReal>& b, std::vector<dReal>& c){
    a.resize(trajectory.dimension);
    b.resize(trajectory.dimension);
    c.resize(trajectory.dimension);
    assert(s >= -TINY && s <= trajectory.duration + TINY);
    if(s < 0)
        s = 0;
    if(s >= trajectory.duration - TINY) {
        int n = ndiscrsteps - 1;
        for(int i = 0; i < trajectory.dimension; i++) {
            a[i] = avect[n][i];
            b[i] = bvect[n][i];
            c[i] = cvect[n][i];
        }
        return;
    }
    int n = int(s / discrtimestep);
    dReal coef = (s - n * discrtimestep) / discrtimestep;
    for(int i = 0; i < trajectory.dimension; i++) {
        a[i] = (1-coef)*avect[n][i] + coef*avect[n+1][i];
        b[i] = (1-coef)*bvect[n][i] + coef*bvect[n+1][i];
        c[i] = (1-coef)*cvect[n][i] + coef*cvect[n+1][i];
    }
}

void TorqueLimits::ComputeSlopeDynamicSingularity(dReal s, dReal sd, std::vector<dReal>& slopesvector) {
    dReal delta = TINY2, s2, ap, bp, cp, slope;
    std::vector<dReal> a, b, c, a2, b2, c2;
    if(s>delta) {
        delta = -delta;
    }
    s2 = s + delta;
    InterpolateDynamics(s,a,b,c);
    InterpolateDynamics(s2,a2,b2,c2);

    slopesvector.resize(0);
    for(int i=0; i<trajectory.dimension; i++) {
        if(std::abs(taumax[i]-taumin[i])<TINY) {
            continue;
        }
        ap = (a2[i]-a[i])/delta;
        bp = (b2[i]-b[i])/delta;
        cp = (c2[i]-c[i])/delta;
        slope = (-bp*sd*sd-cp)/((2*b[i]+ap)*sd);
        slopesvector.push_back(slope);
    }
}

std::pair<dReal,dReal> TorqueLimits::SddLimits(dReal s, dReal sd){
    dReal dtsq = integrationtimestep;
    dtsq = dtsq*dtsq;
    dReal safetybound = discrtimestep/dtsq;
    dReal alpha = -safetybound;
    dReal beta = safetybound;
    dReal sdsq = sd*sd;
    std::vector<dReal> a, b, c;

    dReal taumin_i, taumax_i, alpha_i, beta_i;
    InterpolateDynamics(s,a,b,c);

    for(int i=0; i<trajectory.dimension; i++) {
        if(std::abs(a[i])<TINY || std::abs(taumax[i]-taumin[i])<TINY) {
            continue;
        }
        if(a[i]>0) {
            taumin_i = taumin[i];
            taumax_i = taumax[i];
        }
        else{
            taumin_i = taumax[i];
            taumax_i = taumin[i];
        }
        alpha_i = (taumin_i-sdsq*b[i]-c[i])/a[i];
        beta_i = (taumax_i-sdsq*b[i]-c[i])/a[i];
        alpha = std::max(alpha,alpha_i);
        beta = std::min(beta,beta_i);
    }
    return std::make_pair(alpha,beta);
}


dReal TorqueLimits::SdLimitBobrowInit(dReal s){
    std::pair<dReal,dReal> sddlimits = TorqueLimits::SddLimits(s,0);
    if(sddlimits.first > sddlimits.second) {
        return 0;
    }
    std::vector<dReal> tau_alpha(trajectory.dimension), tau_beta(trajectory.dimension);
    std::vector<dReal> a, b, c;
    InterpolateDynamics(s,a,b,c);

    for(int i=0; i<trajectory.dimension; i++) {
        if(std::abs(taumax[i]-taumin[i])<TINY) {
            continue;
        }
        if(a[i] > 0) {
            tau_alpha[i] = taumin[i];
            tau_beta[i] = taumax[i];
        }
        else{
            tau_alpha[i] = taumax[i];
            tau_beta[i] = taumin[i];
        }
    }
    dReal sdmin = INF;
    for(int k=0; k<trajectory.dimension; k++) {
        for(int m=k+1; m<trajectory.dimension; m++) {
            if(std::abs(taumax[k]-taumin[k])<TINY || std::abs(taumax[m]-taumin[m])<TINY) {
                continue;
            }
            dReal num, denum, r;
            num = a[k]*(tau_alpha[m]-c[m])-a[m]*(tau_beta[k]-c[k]);
            denum = a[k]*b[m]-a[m]*b[k];
            if(std::abs(denum) > TINY) {
                r = num/denum;
                if(r>=0) {
                    sdmin = std::min(sdmin,sqrt(r));
                }
            }
            num = a[m]*(tau_alpha[k]-c[k])-a[k]*(tau_beta[m]-c[m]);
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


void TorqueLimits::FindSingularSwitchPoints(){
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
        for(int j=0; j<trajectory.dimension; j++) {
            if(std::abs(taumax[j]-taumin[j])<TINY) {
                continue;
            }
            if(a[j]*aprev[j]<0) {
                dReal r = (taumin[j]-c[j])/b[j];
                if(r<0) {
                    r = (taumax[j]-c[j])/b[j];
                }
                if(r>0) {
                    found = true;
                    minsd = std::min(minsd,sqrt(r));
                }
            }
        }
        if(found) {
            //std::cout << discrsvect[i] << "," << minsd << "\n";
            AddSwitchPoint(i,SP_SINGULAR,minsd);
        }
        aprev = a;
    }
}
}
