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

#ifndef TOPP_H
#define TOPP_H

#include <iostream>
#include <vector>
#include <list>
#include <math.h>
#include <stdlib.h>
#include <assert.h>


typedef double dReal;

#define TINY 1e-10
#define INF 1e15

namespace TOPP {


enum IntegrationReturnType {
    RT_MAXSTEPS = 0,
    RT_END = 1,
    RT_BOTTOM = 2,
    RT_MVC = 3
};

enum SwitchPointType {
    SP_TANGENT = 0,
    SP_SINGULAR = 1,
    SP_DISCONTINUOUS = 2
};


class Trajectory {
public:
    Trajectory(){
    };
    int dimension;
    dReal duration;
    void Eval(dReal s, std::vector<dReal>& q);
    void Evald(dReal s, std::vector<dReal>& qd);
    void Evaldd(dReal s, std::vector<dReal>& qdd);
};


class Tunings {
public:
    Tunings(){
    }
    dReal discrtimestep;
    dReal integrationtimestep;
};


class SwitchPoint {
    dReal s;
    int switchpointtype;
};


class Profile {
public:
    Profile(std::list<dReal>& slist, std::list<dReal>& sdlist, std::list<dReal>&  sddlist, dReal integrationtimestep);
    std::vector<dReal> svect, sdvect, sddvect;
    dReal integrationtimestep;
    dReal duration;
    int nsteps;
    bool FindTimestepIndex(dReal t, int& index, dReal& remainder);
    dReal Eval(dReal t);
    dReal Evald(dReal t);
    dReal Evaldd(dReal t);
};


class Constraints {
public:

    Trajectory trajectory;
    Tunings tunings;
    int ndiscrsteps;
    std::vector<dReal> discrsvect;
    std::vector<dReal> mvc;
    std::list<SwitchPoint> switchpointslist;



    //////////////////////// General ///////////////////////////

    Constraints(){
    }
    virtual void Preprocess(const Trajectory& trajectory, const Tunings& tunings);
    void Discretize();
    void ComputeMVC();



    //////////////////////// Limits ///////////////////////////

    // upper limit on sd given by the MVC
    virtual dReal SdLimitMVC(dReal s){
        throw "Virtual method not implemented";
    }

    // upper limit on sd given by the direct velocity constraints
    dReal SdLimitDirect(dReal s);

    // lower and upper limits on sdd
    virtual std::pair<dReal,dReal> SddLimits(dReal s, dReal sd){
        throw "Virtual method not implemented";
    }



    ///////////////////////// Switch Points ///////////////////////
    void FindSwitchPoints();
    void AddSwitchPoint(int i, int switchpointtype);
    void FindTangentSwitchPoints();
    virtual void FindSingularSwitchPoints(){
        throw "Virtual method not implemented";
    };
    void FindDiscontinuousSwitchPoints();
    virtual void HandleSingularSwitchPoint(SwitchPoint sp){
        throw "Virtual method not implemented";
    }
};


// Useful functions

//Solves a0 + a1*x + a2*x^2 = 0 in the interval [lowerbound,upperbound]
//If more than one solution, returns the smallest
bool SolveQuadraticEquation(dReal a0, dReal a1, dReal a2, dReal lowerbound, dReal upperbound, dReal& sol);

}

#endif // TOPP_H
