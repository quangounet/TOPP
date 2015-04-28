// -*- coding: utf-8 -*-
// Copyright (C) 2013 Quang-Cuong Pham <cuong.pham@normalesup.org>
//
// This file is part of the Time-Optimal Path Parameterization (TOPP) library.
// TOPP is free software: you can redistribute it and/or modify
// it under the terms of the GNU Lesser General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// at your option, any later version.
//
// This program is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU Lesser General Public License for more details.
//
// You should have received a copy of the GNU Lesser General Public License
// along with this program.  If not, see <http://www.gnu.org/licenses/>.



#include "PolygonConstraints.h"


namespace TOPP {

PolygonConstraints::PolygonConstraints(const std::string& constraintsstring){
    std::string buff;
    std::istringstream iss(constraintsstring);
    std::vector<dReal> tmpvect;
    std::vector<std::pair<dReal,dReal> > polygon;
    polygonvector.resize(0);
    getline(iss, buff, '\n');
    discrtimestep = atof(buff.c_str());
    getline(iss, buff, '\n');
    VectorFromString(buff,vmax);
    while(iss.good()) {
        getline(iss, buff, '\n');
        VectorFromString(buff, tmpvect);
        polygon.resize(0);
        for(int i=0; i <(int) tmpvect.size() /2; i++) {
            polygon.push_back(std::pair<dReal,dReal>(tmpvect[2*i],tmpvect[2*i+1]));
        }
        polygonvector.push_back(polygon);
    }
    hasvelocitylimits =  VectorMax(vmax) > TINY;
}


void PolygonConstraints::Discretize(){
    dReal s;
    ndiscrsteps = int((trajectory.duration+TINY)/discrtimestep);
    ndiscrsteps++;
    discrsvect.resize(0);
    for(int i=0; i<ndiscrsteps; i++) {
        s = i*discrtimestep;
        discrsvect.push_back(s);
    }
}


std::pair<dReal,dReal> PolygonConstraints::SdLimitBobrowInitDiscrete(int i){
    std::vector<std::pair<dReal,dReal> > polygon = polygonvector[i];
    dReal sdmax = 0;
    dReal a = 0;
    for(int j=0; j<(int) polygon.size(); j++) {
        if(sqrt(polygon[j].second) > sdmax) {
            a = polygon[j].first;
            sdmax = sqrt(polygon[j].second);
        }
    }
    return std::pair<dReal,dReal>(a,sdmax);
}

void PolygonConstraints::ComputeMVCBobrow() {
    mvcbobrow.resize(ndiscrsteps);
    accelonMVC.resize(ndiscrsteps);
    for(int i=0; i<ndiscrsteps; i++) {
        std::pair<dReal,dReal> p = SdLimitBobrowInitDiscrete(i);
        accelonMVC[i] = p.first;
        mvcbobrow[i] = p.second;
    }
}

// Dummy since we reimplement ComputeMVCBobrow
dReal PolygonConstraints::SdLimitBobrowInit(dReal s){
    BOOST_ASSERT(false);
    return 0;
}

// Need to take the square root before returning since the polygon is in the sd^2 space
std::pair<dReal,dReal> PolygonConstraints::SddLimitsDiscrete(int i, dReal sd){
    dReal alpha, beta;
    std::vector<std::pair<dReal,dReal> > polygon = polygonvector[i];
    if(sd <= 0) {
        alpha = polygon[0].first;
        beta = polygon[polygon.size()-1].first;
        return std::pair<dReal,dReal>(alpha,beta);
    }
    if(sd >= mvcbobrow[i]) {
        alpha = accelonMVC[i];
        beta = accelonMVC[i];
        return std::pair<dReal,dReal>(alpha,beta);
    }
    int stage1 = 0;
    dReal sdsq = sd*sd;
    for(int j=0; j<(int) polygon.size(); j++) {
        if(polygon[j].second > sdsq) {
            stage1 = j;
            dReal coef = (sdsq-polygon[j-1].second)/(polygon[j].second-polygon[j-1].second);
            alpha = (1-coef)*polygon[j-1].first + coef*polygon[j].first;
            break;
        }
    }
    for(int j=stage1; j<(int) polygon.size(); j++) {
        if(polygon[j].second < sdsq) {
            dReal coef = (sdsq-polygon[j].second)/(polygon[j-1].second-polygon[j].second);
            beta = (1-coef)*polygon[j].first + coef*polygon[j-1].first;
            break;
        }
    }
    //std::cout << i*discrtimestep << " " << mvcbobrow[i] << " " << sd << " " << alpha << " " << beta << "\n";
    return std::pair<dReal,dReal>(alpha,beta);
}

std::pair<dReal,dReal> PolygonConstraints::SddLimits(dReal s, dReal sd){
    BOOST_ASSERT(s>=-TINY && s<=trajectory.duration+TINY);
    if(s<0) {
        s=0;
    }
    if(s>=trajectory.duration) {
        int n = ndiscrsteps-1;
        return SddLimitsDiscrete(n, sd);
    }
    int n = int(s/discrtimestep);
    dReal coef = (s-n*discrtimestep)/discrtimestep;
    std::pair<dReal,dReal> a, b;
    a = SddLimitsDiscrete(n, sd);
    b = SddLimitsDiscrete(n+1, sd);
    return std::pair<dReal,dReal>((1-coef)*a.first + coef*b.first,(1-coef)*a.second + coef*b.second);
}

void PolygonConstraints::ComputeSlopeDynamicSingularity(dReal s, dReal sd, std::vector<dReal>& slopesvector) {
    // Simply returns a large number of discretized slopes
    dReal maxslope = 10;
    int nslopes = 100;
    slopesvector.resize(2*nslopes);
    for(int i=0; i<nslopes; i++) {
        dReal ee = i*maxslope/nslopes;
        slopesvector[2*i] = ee*ee;
        slopesvector[2*i+1] = -ee*ee;
    }
}


void PolygonConstraints::FindSingularSwitchPoints(){
    // Simply walk on the MVC and detect discontinuity of the tangent vector
    if(ndiscrsteps<3)
        return;
    int i = 1;
    dReal sd,sdnext,tangent,prevtangent;
    sd = mvcbobrow[i];
    sdnext = mvcbobrow[i+1];
    prevtangent = (sdnext-sd)/discrtimestep;
    for(int i=2; i<ndiscrsteps-1; i++) {
        sd = mvcbobrow[i];
        sdnext = mvcbobrow[i+1];
        tangent = (sdnext-sd)/discrtimestep;
        if(std::abs(tangent-prevtangent)>1) {
            AddSwitchPoint2(discrsvect[i],sd,SP_SINGULAR);
        }
        prevtangent = tangent;
    }
}

}
