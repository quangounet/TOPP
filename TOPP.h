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


#ifndef TOPP_H
#define TOPP_H

#include <iostream>
#include <sstream>
#include <string>
#include <iomanip>

#include <vector>
#include <list>
#include <cmath>
//#include <armadillo>
#include <stdlib.h>
#include <assert.h>

typedef double dReal;

#define TINY 1e-10
#define TINY2 1e-5
#define INF 1.0e15



namespace TOPP {



////////////////////////////////////////////////////////////////////
/////////////////////////// Tunings ////////////////////////////////
////////////////////////////////////////////////////////////////////

class Tunings {
public:
    Tunings(){
    }
    Tunings(const std::string& tuningsstring);
    dReal discrtimestep;
    dReal integrationtimestep;
    dReal sdprecision;
    int passswitchpointnsteps;
    dReal reparamtimestep;
};




////////////////////////////////////////////////////////////////////
/////////////////////////// Switch Point ///////////////////////////
////////////////////////////////////////////////////////////////////

enum SwitchPointType {
    SP_TANGENT = 0,
    SP_SINGULAR = 1,
    SP_DISCONTINUOUS = 2,
    SP_ZLAJPAH = 3
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
    Profile(std::list<dReal>&slist, std::list<dReal>&sdlist, std::list<dReal>&sddlist, dReal integrationtimestep);
    Profile(){
    }
    std::vector<dReal> svect, sdvect, sddvect;
    dReal integrationtimestep;
    dReal duration;
    bool forward;
    int nsteps;
    int currentindex;
    bool FindTimestepIndex(dReal t, int &index, dReal& remainder);
    bool Invert(dReal s,  dReal& t, bool searchbackward=false);
    dReal Eval(dReal t);
    dReal Evald(dReal t);
    dReal Evaldd(dReal t);
    void Print();
    void Write(std::stringstream& ss, dReal dt=0.001);

};



////////////////////////////////////////////////////////////////////
/////////////////////////// Trajectory /////////////////////////////
////////////////////////////////////////////////////////////////////


class Polynomial {
public:
    Polynomial(const std::vector<dReal>& coefficientsvector);
    Polynomial(const std::string& s);
    Polynomial(){
    }
    void InitFromCoefficientsVector(const std::vector<dReal>&coefficientsvector);
    int degree;
    std::vector<dReal> coefficientsvector;
    std::vector<dReal> coefficientsvectord;
    std::vector<dReal> coefficientsvectordd;
    dReal Eval(dReal s);
    dReal Evald(dReal s);
    dReal Evaldd(dReal s);
    void Write(std::stringstream& ss);
};


class Chunk {
public:
    Chunk(dReal duration, const std::vector<Polynomial>& polynomialsvector);
    Chunk(){
    };
    int dimension;
    int degree;
    dReal duration;
    dReal sbegin, send;
    std::vector<Polynomial> polynomialsvector;
    void Eval(dReal s, std::vector<dReal>&q);
    void Evald(dReal s, std::vector<dReal>&qd);
    void Evaldd(dReal s, std::vector<dReal>&qdd);
    void Write(std::stringstream& ss);
};


class Trajectory {
public:
    Trajectory(const std::list<Chunk>& chunkslist);
    Trajectory(const std::string& trajectorystring);
    Trajectory(){
    }
    void InitFromChunksList(const std::list<Chunk>&chunkslist);

    int dimension;
    dReal duration;
    int degree;

    std::list<Chunk> chunkslist;
    std::list<dReal> chunkdurationslist;
    std::list<dReal> chunkcumulateddurationslist;

    void FindChunkIndex(dReal s, int& index, dReal& remainder);
    void Eval(dReal s, std::vector<dReal>&q);
    void Evald(dReal s, std::vector<dReal>&qd);
    void Evaldd(dReal s, std::vector<dReal>&qdd);
    void ComputeChunk(dReal t, dReal t0, dReal s, dReal sd, dReal sdd, const Chunk& currentchunk, Chunk& newchunk);
    void SPieceToChunks(dReal s, dReal sd, dReal sdd, dReal T, int& currentchunkindex, dReal& processedcursor, std::list<Chunk>::iterator& itcurrentchunk, std::list<Chunk>& chunkslist);
    void Reparameterize(std::list<Profile>& profileslist, dReal integrationtimestep, Trajectory& restrajectory);
    void Write(std::stringstream& ss);
};





////////////////////////////////////////////////////////////////////
/////////////////////////// Constraints ////////////////////////////
////////////////////////////////////////////////////////////////////

class Constraints {
public:

    Trajectory trajectory;
    Tunings tunings;
    int ndiscrsteps;
    std::vector<dReal> discrsvect;
    std::vector<dReal> mvcbobrow;
    std::vector<dReal> mvccombined;
    bool hasvelocitylimits;
    std::vector<dReal> vmax;

    std::list<SwitchPoint> switchpointslist;
    std::list<SwitchPoint> zlajpahlist;

    //////////////////////// General ///////////////////////////

    Constraints(){
    }
    virtual void Preprocess(Trajectory& trajectory, const Tunings &tunings);
    void Discretize();
    void ComputeMVCBobrow();
    void ComputeMVCCombined();
    void WriteMVCBobrow(std::stringstream& ss, dReal dt=0.01);
    void WriteMVCCombined(std::stringstream& ss, dReal dt=0.01);

    // discretize dynamics
    virtual void DiscretizeDynamics(){
        std::cout << "Virtual method not implemented\n";
        throw "Virtual method not implemented";
    }


    dReal Interpolate1D(dReal s, const std::vector<dReal>& v);


    //////////////////////// Limits ///////////////////////////

    // upper limit on sd given by acceleration constraints (MVC)

    dReal SdLimitBobrow(dReal s);
    virtual dReal SdLimitBobrowInit(dReal s){
        std::cout << "Virtual method not implemented\n";
        throw "Virtual method not implemented";
    }

    // upper limit on sd given by the combined velocity constraints
    dReal SdLimitCombined(dReal s);
    dReal SdLimitCombinedInit(dReal s);

    // lower and upper limits on sdd
    virtual std::pair<dReal,dReal> SddLimits(dReal s, dReal sd){
        std::cout << "Virtual method not implemented\n";
        throw "Virtual method not implemented";
    }

    ///////////////////////// Switch Points ///////////////////////
    void FindSwitchPoints();
    void AddSwitchPoint(int i, int switchpointtype);
    void FindTangentSwitchPoints();
    virtual void FindSingularSwitchPoints(){
        std::cout << "Virtual method not implemented\n";
        throw "Virtual method not implemented";
    };
    void FindDiscontinuousSwitchPoints();
};





////////////////////////////////////////////////////////////////////
//////////////////////// Integration ///////////////////////////////
////////////////////////////////////////////////////////////////////

enum IntegrationReturnType {
    INT_MAXSTEPS = 0,     // reached maximum number of steps
    INT_END = 1,          // reached s = send (FW) or s = 0 (BW)
    INT_BOTTOM = 2,       // crossed sd = 0
    INT_MVC = 3,          // crossed the MVC
    INT_PROFILE = 4       // crossed a profile
};

enum CLCReturnType {
    CLC_OK = 0,        // no problem
    CLC_SWITCH = 1,    // could not pass a switchpoint
    CLC_BOTTOM = 2     // crossed sd = 0
};


static std::list<Profile> voidprofileslist;

int IntegrateForward(Constraints& constraints, dReal sstart, dReal sdstart, dReal dt, Profile& resprofile, int maxsteps=1e5, std::list<Profile>& testprofileslist = voidprofileslist, bool testmvc=true, bool zlajpah=true);

int IntegrateBackward(Constraints& constraints, dReal sstart, dReal sdstart, dReal dt, Profile& resprofile, int maxsteps=1e5, std::list<Profile>& testprofileslist = voidprofileslist, bool testmvc=true);

int ComputeLimitingCurves(Constraints& constraints, std::list<Profile>&resprofileslist, bool zlajpah = false);

// Path parameterization
int PP(Constraints& constraints, Trajectory& trajectory, Tunings& tunings, dReal sdbeg, dReal sdend, Trajectory& restrajectory, std::list<Profile>&resprofileslist);

// Velocity Interval Propagation
int VIP(Constraints& constraints, Trajectory& trajectory, Tunings& tunings, dReal sdbegmin, dReal sdbegmax, dReal& sdendmin, dReal& sdendmax, std::list<Profile>&resprofileslist);





//////// ////////////////////////////////////////////////////////////
///////////////////////// Utilities ////////////////////////////////
////////////////////////////////////////////////////////////////////

// Find the smallest and largest element of a vector
dReal VectorMin(const std::vector<dReal>& v);
dReal VectorMax(const std::vector<dReal>& v);

// Read a vector of dReal from a space separated string
void VectorFromString(const std::string& s,std::vector<dReal>&resvect);

//Solve a0 + a1*x + a2*x^2 = 0 in the interval [lowerbound,upperbound]
//If more than one solution, returns the smallest
bool SolveQuadraticEquation(dReal a0, dReal a1, dReal a2, dReal& sol, dReal lowerbound=-INF, dReal upperbound=INF);

// Check whether the point (s,sd) is above at least one profile in the list
bool IsAboveProfilesList(dReal s, dReal sd, std::list<Profile>&testprofileslist, bool searchbackward=false, bool reinitialize=false);

// Find the lowest profile at s
bool FindLowestProfile(dReal s, Profile& profile, dReal& tres, std::list<Profile>&testprofileslist);

}

#endif // TOPP_H
