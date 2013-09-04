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
#include <cmath>
//#include <armadillo>
#include <stdlib.h>
#include <assert.h>


typedef double dReal;

#define TINY 1e-10
#define TINY2 1e-5
#define INF 1e15

namespace TOPP {



////////////////////////////////////////////////////////////////////
/////////////////////////// Tunings ////////////////////////////////
////////////////////////////////////////////////////////////////////

class Tunings {
public:
    Tunings(){
    }
    dReal discrtimestep;
    dReal integrationtimestep;
    dReal threshold;
    int passswitchpointnsteps;
};




////////////////////////////////////////////////////////////////////
/////////////////////////// Trajectory /////////////////////////////
////////////////////////////////////////////////////////////////////

class Trajectory {
public:
    Trajectory(){
    };
    int dimension;
    dReal duration;
    virtual void Eval(dReal s, std::vector<dReal>& q){
    }
    virtual void Evald(dReal s, std::vector<dReal>& qd){
    }
    virtual void Evaldd(dReal s, std::vector<dReal>& qdd){
    }
};



////////////////////////////////////////////////////////////////////
/////////////////////////// Switch Point ///////////////////////////
////////////////////////////////////////////////////////////////////

enum SwitchPointType {
    SPT_TANGENT = 0,
    SPT_SINGULAR = 1,
    SPT_DISCONTINUOUS = 2
};

class SwitchPoint {
public:
    SwitchPoint(dReal s0, dReal sd0, int switchpointtype0){
        s = s0;
        sd = sd0;
        switchpointtype = switchpointtype0;
    }
    dReal s, sd;
    int switchpointtype;
};




////////////////////////////////////////////////////////////////////
/////////////////////////// Profile ////////////////////////////////
////////////////////////////////////////////////////////////////////

class Profile {
public:
    Profile(std::list<dReal>& slist, std::list<dReal>& sdlist, std::list<dReal>&  sddlist, dReal integrationtimestep);
    Profile(std::list<Profile>& profileslist, dReal integrationtimestep0);
    Profile(){
    }
    std::vector<dReal> svect, sdvect, sddvect;
    dReal integrationtimestep;
    dReal duration;
    bool forward;
    int nsteps;
    int currentindex;
    bool FindTimestepIndex(dReal t, int& index, dReal& remainder);
    bool Invert(dReal s,  dReal& t, bool searchbackward=false);
    dReal Eval(dReal t);
    dReal Evald(dReal t);
    dReal Evaldd(dReal t);
    void Print();
};




////////////////////////////////////////////////////////////////////
/////////////////////////// Constraints ////////////////////////////
////////////////////////////////////////////////////////////////////

class Constraints {
public:

    Trajectory* ptrajectory;
    Tunings tunings;
    int ndiscrsteps;
    std::vector<dReal> discrsvect;
    std::vector<dReal> mvc;
    std::list<SwitchPoint> switchpointslist;

    //////////////////////// General ///////////////////////////

    Constraints(){
    }
    virtual void Preprocess(Trajectory& trajectory, const Tunings& tunings);
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




////////////////////////////////////////////////////////////////////
//////////////////////// Integration ///////////////////////////////
////////////////////////////////////////////////////////////////////

enum IntegrationReturnType {
    IRT_MAXSTEPS = 0,
    IRT_END = 1,
    IRT_BOTTOM = 2,
    IRT_MVC = 3,
    IRT_ABOVEPROFILES = 4
};

enum CLCReturnType {
    CLC_OK = 0,   // no problem
    CLC_SP = 1    // cannot cross a switchpoint
};


static std::list<Profile> voidprofileslist;

int IntegrateForward(Constraints& constraints, dReal sstart, dReal sdstart, Profile& profile, int maxsteps=1e5, std::list<Profile>& testprofileslist = voidprofileslist);

int IntegrateBackward(Constraints& constraints, dReal sstart, dReal sdstart, Profile& profile, int maxsteps=1e5, std::list<Profile>& testprofileslist = voidprofileslist);

int ComputeLimitingCurves(Constraints& constraints, std::list<Profile>& resprofileslist);

int IntegrateAllProfiles();




////////////////////////////////////////////////////////////////////
///////////////////////// Utilities ////////////////////////////////
////////////////////////////////////////////////////////////////////

//Solves a0 + a1*x + a2*x^2 = 0 in the interval [lowerbound,upperbound]
//If more than one solution, returns the smallest
bool SolveQuadraticEquation(dReal a0, dReal a1, dReal a2, dReal lowerbound, dReal upperbound, dReal& sol);

bool IsAboveProfilesList(dReal s, dReal sd, std::list<Profile>& testprofileslist, bool searchbackward=false, bool reinitialize=false);

bool FindLowestProfile(dReal s, Profile& profile, dReal& tres, std::list<Profile>&testprofileslist);

}

#endif // TOPP_H
