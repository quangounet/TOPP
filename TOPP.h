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
#include <ctime>
#include <chrono>
#include <algorithm>

#include <vector>
#include <list>
#include <cmath>
//#include <armadillo>
#include <cstdio>
#include <cstdlib>
#include <cassert>

typedef double dReal;

#define TINY 1e-10
#define TINY2 1e-5
#define INF 1.0e15
#define MAXSD 200

#include "Trajectory.h"


namespace TOPP {


const int BUFFSIZE = 300000;


////////////////////////////////////////////////////////////////////
/////////////////////////// Tunings ////////////////////////////////
////////////////////////////////////////////////////////////////////

// Some tunings parameters
class Tunings {
public:
    Tunings(){
    }
    Tunings(const std::string& tuningsstring);
    dReal discrtimestep; // Time step to discretize the trajectory, usually 0.01
    dReal integrationtimestep; // Time step to integrate the profiles, usually 0.01
    dReal reparamtimestep; // Time step for the reparameterization, usually 0.01. If 0, set reparamtimestep so as to keep the number of discretization points unchanged
    int passswitchpointnsteps; // Number of steps to integrate from the switch point when assessing whether it's addressable
    dReal bisectionprecision; //Precision for the sd search, set to 0.01 by default
    dReal loweringcoef; //While addressing switchpoints, lower sd by loweringcoef. Set to 0.9
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

// Velocity profile
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
    // Find t such that Eval(t) = s
    // Return false if no such t
    bool Invert(dReal s,  dReal& t, bool searchbackward=false);
    dReal Eval(dReal t);
    dReal Evald(dReal t);
    dReal Evaldd(dReal t);
    void Print();
    void Write(std::stringstream& ss, dReal dt=0.001);

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
    int maxrep; // Max number of reps to try integrating profiles (reduce integrationtimestep at each rep). Set to 1.
    std::vector<dReal> vmax;

    std::list<SwitchPoint> switchpointslist; // list of switch points
    std::list<std::pair<dReal,dReal> > zlajpahlist; // list of zlajpah points
    std::list<Profile> resprofileslist; // resulting profiles

    dReal resduration;


    //////////////////////// General ///////////////////////////

    Constraints(){
    }
    virtual void Preprocess(Trajectory& trajectory, Tunings &tunings);

    // Discretize the time interval
    virtual void Discretize();

    // Compute the MVC given by acceleration constraints
    void ComputeMVCBobrow();

    // Compute the combined MVC (incorporating pure velocity constraints)
    void ComputeMVCCombined();

    // Write the MVC to stringstreams
    void WriteMVCBobrow(std::stringstream& ss, dReal dt=0.01);
    void WriteMVCCombined(std::stringstream& ss, dReal dt=0.01);

    // Linear interpolation
    dReal Interpolate1D(dReal s, const std::vector<dReal>& v);


    //////////////////////// Limits ///////////////////////////

    // Upper limit on sd given by acceleration constraints (Bobrow)
    dReal SdLimitBobrow(dReal s);

    // Compute the maximum velocity curve due to dynamics at s
    // Called at initialization
    virtual dReal SdLimitBobrowInit(dReal s){
        std::cout << "Virtual method not implemented\n";
        throw "Virtual method not implemented";
    }

    // Upper limit on sd after incorporating pure velocity constraints
    dReal SdLimitCombined(dReal s);
    dReal SdLimitCombinedInit(dReal s);

    // Pair of (lower,upper) limits on sdd
    virtual std::pair<dReal,dReal> SddLimits(dReal s, dReal sd){
        std::cout << "Virtual method not implemented\n";
        throw "Virtual method not implemented";
    }


    ///////////////////////// Switch Points ///////////////////////

    // Find all switch points, add them to switchpointslist
    void FindSwitchPoints();
    void FindTangentSwitchPoints();
    void FindDiscontinuousSwitchPoints();
    // Trim nearby switch points, priority is given to singular switch points
    void TrimSwitchPoints();

    // Compute the slope of the profiles near a dynamic singularity
    virtual void ComputeSlopeDynamicSingularity(dReal s, dReal sd, std::vector<dReal>& slopesvector){
        std::cout << "Virtual method not implemented\n";
        throw "Virtual method not implemented";
    }
    virtual void FindSingularSwitchPoints(){
        std::cout << "Virtual method not implemented\n";
        throw "Virtual method not implemented";
    };

    // Add a switch point to switchpointslist
    void AddSwitchPoint(int i, int switchpointtype, dReal sd = -1);

};




////////////////////////////////////////////////////////////////////
/////////////////// Quadratic Constraints //////////////////////////
////////////////////////////////////////////////////////////////////

// Constraints of the form a(s)sdd + b(s)sd^2 + c(s) <= 0
class QuadraticConstraints : public Constraints {
public:
    QuadraticConstraints() : Constraints(){
    }
    QuadraticConstraints(const std::string& constraintsstring);

    //////////////// Overloaded methods //////////////////////
    std::pair<dReal,dReal> SddLimits(dReal s, dReal sd);
    dReal SdLimitBobrowInit(dReal s);
    void FindSingularSwitchPoints();
    void ComputeSlopeDynamicSingularity(dReal s, dReal sd, std::vector<dReal>& slopesvector);

    //////////////// Specific members and methods //////////////////////
    int nconstraints;  // Number of constraints
    std::vector<std::vector<dReal> > avect, bvect, cvect;  // Dynamics coefficients
    void InterpolateDynamics(dReal s, std::vector<dReal>& a, std::vector<dReal>& b, std::vector<dReal>& c);   // Linearly interpolate the dynamics coefficients a,b,c

};




////////////////////////////////////////////////////////////////////
//////////////////////// Integration ///////////////////////////////
////////////////////////////////////////////////////////////////////

// Return type for the integration methods
enum IntegrationReturnType {
    INT_MAXSTEPS = 0,     // reached maximum number of steps
    INT_END = 1,          // reached s = send (FW) or s = 0 (BW)
    INT_BOTTOM = 2,       // crossed sd = 0
    INT_MVC = 3,          // crossed the MVC
    INT_PROFILE = 4       // crossed a profile
};

// Return type for the CLC computation
enum CLCReturnType {
    CLC_OK = 0,        // no problem
    CLC_SWITCH = 1,    // could not pass a switchpoint
    CLC_BOTTOM = 2     // crossed sd = 0
};

// Integrate forward from (sstart,sdstart)
int IntegrateForward(Constraints& constraints, dReal sstart, dReal sdstart, dReal dt, Profile& resprofile, int maxsteps=1e5, bool testaboveexistingprofiles=true, bool testmvc=true, bool zlajpah=false);

// Integrate backward from (sstart,sdstart)
int IntegrateBackward(Constraints& constraints, dReal sstart, dReal sdstart, dReal dt, Profile& resprofile, int maxsteps=1e5, bool testaboveexistingprofiles=true, bool testmvc=true, bool zlajpah=false);

// Compute the CLC
int ComputeLimitingCurves(Constraints& constraints);

// Compute all the profiles (CLC, forward from 0, backward from send)
int ComputeProfiles(Constraints& constraints, Trajectory& trajectory, Tunings& tunings, dReal sdbeg, dReal sdend);

// Velocity Interval Propagation
int VIP(Constraints& constraints, Trajectory& trajectory, Tunings& tunings, dReal sdbegmin, dReal sdbegmax, dReal& sdendmin, dReal& sdendmax);




//////// ////////////////////////////////////////////////////////////
///////////////////////// Utilities ////////////////////////////////
////////////////////////////////////////////////////////////////////

// Find the smallest and largest element of a vector
dReal VectorMin(const std::vector<dReal>& v);
dReal VectorMax(const std::vector<dReal>& v);
std::vector<dReal> VectorAdd(const std::vector<dReal>& a, const std::vector<dReal>& b, dReal coefa=1, dReal coefb=1);

// Read a vector of dReal from a space-separated string
void VectorFromString(const std::string& s,std::vector<dReal>&resvect);

// Solve a0 + a1*x + a2*x^2 = 0 in the interval [lowerbound,upperbound]
// Return false if no solution in [lowerbound,upperbound], true otherwise
// If more than one solution, choose the smallest
bool SolveQuadraticEquation(dReal a0, dReal a1, dReal a2, dReal& sol, dReal lowerbound=-INF, dReal upperbound=INF);

// Check whether the point (s,sd) is above at least one profile in the list
bool IsAboveProfilesList(dReal s, dReal sd, std::list<Profile>& resprofileslist, bool searchbackward=false);

// Find the lowest profile at s
// Return false if no profile covers s, true otherwise
bool FindLowestProfile(dReal s, Profile& profile, dReal& tres, std::list<Profile>& resprofileslist);

}

#endif // TOPP_H
