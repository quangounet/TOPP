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

#include "Trajectory.h"

#include <boost/lexical_cast.hpp>
#include <boost/format.hpp>

#ifdef WITH_OPENRAVE
#include <openrave/openrave.h>
#endif

#ifdef _MSC_VER
#include <boost/typeof/std/string.hpp>
#include <boost/typeof/std/vector.hpp>
#include <boost/typeof/std/list.hpp>
#include <boost/typeof/std/map.hpp>
#include <boost/typeof/std/string.hpp>
#include <boost/typeof/typeof.hpp>

#define FOREACH(it, v) for(BOOST_TYPEOF(v) ::iterator it = (v).begin(); it != (v).end(); (it)++)
#define FOREACH_NOINC(it, v) for(BOOST_TYPEOF(v) ::iterator it = (v).begin(); it != (v).end(); )

#define FOREACHC(it, v) for(BOOST_TYPEOF(v) ::const_iterator it = (v).begin(); it != (v).end(); (it)++)
#define FOREACHC_NOINC(it, v) for(BOOST_TYPEOF(v) ::const_iterator it = (v).begin(); it != (v).end(); )

#else

#if __cplusplus > 199711L || defined(__GXX_EXPERIMENTAL_CXX0X__)
#define FOREACH(it, v) for(decltype((v).begin()) it = (v).begin(); it != (v).end(); (it)++)
#define FOREACH_NOINC(it, v) for(decltype((v).begin()) it = (v).begin(); it != (v).end(); )
#define FOREACHC FOREACH
#define FOREACHC_NOINC FOREACH_NOINC
#else
#define FOREACH(it, v) for(typeof((v).begin())it = (v).begin(); it != (v).end(); (it)++)
#define FOREACH_NOINC(it, v) for(typeof((v).begin())it = (v).begin(); it != (v).end(); )
#define FOREACHC FOREACH
#define FOREACHC_NOINC FOREACH_NOINC
#endif

#endif

namespace TOPP {

#ifdef WITH_OPENRAVE
typedef OpenRAVE::dReal dReal;
#else
typedef double dReal;
#endif

#define TINY 1e-10
#define TINY2 1e-5
#define INF 1.0e15
#define MAXSD 200

/// \brief Exception that all OpenRAVE internal methods throw; the error codes are held in \ref OpenRAVEErrorCode.
class TOPPException : public std::exception
{
public:
    TOPPException() : std::exception(), _s("unknown exception"), _errorcode(0) {
    }
    TOPPException(const std::string& s, int errorcode=0) : std::exception() {
        _errorcode = errorcode;
        _s = "openrave (";
        _s += boost::lexical_cast<std::string>(_errorcode);
        _s += "): ";
        _s += s;
    }
    virtual ~TOPPException() throw() {
    }
    char const* what() const throw() {
        return _s.c_str();
    }
    const std::string& message() const {
        return _s;
    }
    int GetCode() const {
        return _errorcode;
    }
private:
    std::string _s;
    int _errorcode;
};

#define TOPP_EXCEPTION_FORMAT0(s, errorcode) TOPP::TOPPException(boost::str(boost::format("[%s:%d] " s)%(__PRETTY_FUNCTION__)%(__LINE__)),errorcode)

/// adds the function name and line number to an TOPP exception
#define TOPP_EXCEPTION_FORMAT(s, args,errorcode) TOPP::TOPPException(boost::str(boost::format("[%s:%d] " s)%(__PRETTY_FUNCTION__)%(__LINE__)%args),errorcode)

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
    std::vector<dReal> slopesvector;
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

    // Tuning parameters
    dReal discrtimestep; // Time step to discretize the trajectory when computing the MVC
    dReal integrationtimestep; // Time step to integrate the profiles
    dReal reparamtimestep; // Time step to reparameterize the trajectory based on the optimal profile
    int passswitchpointnsteps; // Number of steps to integrate around a switch point
    int extrareps; // Number of reps to lower the integrationtimestep when integration fails
    dReal bisectionprecision; // Precision when determining the sd for tangent switch point
    dReal loweringcoef; // Lowerbound that multiplies sd when searching sd for tangent switch point

    // Maximum Velocity Curves
    int ndiscrsteps; // Number of discretization steps, around trajectory.duration/discrtimestep
    std::vector<dReal> discrsvect; // Discretization points on the s axis
    std::vector<dReal> mvcbobrow;
    std::vector<dReal> mvccombined;
    bool hasvelocitylimits;
    std::vector<dReal> vmax;

    std::list<SwitchPoint> switchpointslist; // list of switch points
    std::list<std::pair<dReal,dReal> > zlajpahlist; // list of zlajpah points
    std::list<Profile> resprofileslist; // resulting profiles

    dReal resduration;


    //////////////////////// General ///////////////////////////

    Constraints(){
    }

    // Compute the MVCs and the switchpoints and other initializations
    virtual bool Preprocess();

    // Discretize the time interval
    virtual void Discretize();

    // Compute the MVC given by acceleration constraints
    virtual void ComputeMVCBobrow();

    // Compute the combined MVC (incorporating pure velocity constraints)
    virtual void ComputeMVCCombined();

    // Write the MVC to stringstreams
    virtual void WriteMVCBobrow(std::stringstream& ss, dReal dt=0.01);
    virtual void WriteMVCDirect(std::stringstream& ss, dReal dt=0.01);
    virtual void WriteExtra(std::stringstream& ss){
        return;
    }
    virtual void WriteConstraints(std::stringstream& ss){
        return;
    }

    // Linear interpolation
    virtual dReal Interpolate1D(dReal s, const std::vector<dReal>& v);


    //////////////////////// Limits ///////////////////////////

    // Upper limit on sd given by acceleration constraints (Bobrow)
    virtual dReal SdLimitBobrow(dReal s);

    // Compute the maximum velocity curve due to dynamics at s
    // Called at initialization
    virtual dReal SdLimitBobrowInit(dReal s){
        std::cout << "Virtual method not implemented\n";
        throw TOPPException("Virtual method not implemented");
    }

    // Upper limit on sd after incorporating pure velocity constraints
    virtual dReal SdLimitCombined(dReal s);
    virtual dReal SdLimitCombinedInit(dReal s);

    // Pair of (lower,upper) limits on sdd
    virtual std::pair<dReal,dReal> SddLimits(dReal s, dReal sd){
        std::cout << "Virtual method not implemented\n";
        throw TOPPException("Virtual method not implemented");
    }
    virtual dReal SddLimitAlpha(dReal s, dReal sd){
        return SddLimits(s, sd).first;
    }
    virtual dReal SddLimitBeta(dReal s, dReal sd){
        return SddLimits(s, sd).second;
    }


    ///////////////////////// Switch Points ///////////////////////

    // Find all switch points, add them to switchpointslist
    virtual void FindSwitchPoints();
    virtual void FindTangentSwitchPoints();
    virtual void FindDiscontinuousSwitchPoints();

    // Switch points that are close to each other will be replaced by a single swtich point
    // Modifies switchpointslist
    virtual void TrimSwitchPoints();

    // Compute the slope of the profiles near a dynamic singularity
    virtual void ComputeSlopeDynamicSingularity(dReal s, dReal sd, std::vector<dReal>& slopesvector){
        std::cout << "Virtual method not implemented\n";
        throw TOPPException("Virtual method not implemented");
    }

    // Fix the integration at s=0 when there is a singularity there
    // If there's nothing to do then sstartnew = 0
    // Else sstartnew > 0 and sdstartnew will be the value that allows going through the singularity
    virtual void FixStart(dReal& sstartnew,dReal& sdstartnew){
        sstartnew = 0;
    }

    // Fix the integration at s=send when there is a singularity there
    // If there's nothing to do then sendnew = send
    // Else sendnew < send and sdendnew will be the value that allows going through the singularity
    virtual void FixEnd(dReal& sendnew,dReal& sdendnew){
        sendnew = trajectory.duration;
    }

    // Find all the singular switch points
    // Add them to switchpointslist
    virtual void FindSingularSwitchPoints(){
        std::cout << "Virtual method not implemented\n";
        throw TOPPException("Virtual method not implemented");
    };

    // Add a switch point to switchpointslist
    virtual void AddSwitchPoint(int i, int switchpointtype, dReal sd = -1);

};




////////////////////////////////////////////////////////////////////
/////////////////// Quadratic Constraints //////////////////////////
////////////////////////////////////////////////////////////////////

// Constraints of the form a(s)sdd + b(s)sd^2 + c(s) <= 0
// See our article http://arxiv.org/abs/1312.6533 for more details
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
    void WriteConstraints(std::stringstream& ss);

    //////////////// Specific members and methods //////////////////////
    int nconstraints;  // Number of constraints
    std::vector<std::vector<dReal> > avect, bvect, cvect;  // Dynamics coefficients. avect[i], bvect[i], cvect[i] are vectors of length 2*ndof where the first ndof are the upper limit, the next ndof are for the lower limit. These incorporate any upper/lower limits.

    void InterpolateDynamics(dReal s, std::vector<dReal>& a, std::vector<dReal>& b, std::vector<dReal>& c);   // Linearly interpolate the dynamics coefficients a,b,c
    void FixStart(dReal& sstartnew,dReal& sdstartnew);
    void FixEnd(dReal& sendnew,dReal& sdendnew);

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
int ComputeProfiles(Constraints& constraints, dReal sdbeg, dReal sdend);

// Velocity Interval Propagation
int VIP(Constraints& constraints, dReal sdbegmin, dReal sdbegmax, dReal& sdendmin, dReal& sdendmax);

// Velocity Interval Propagation from sdend backwards
int VIPBackward(Constraints& constraints, dReal& sdbegmin, dReal& sdbegmax, dReal sdendmin, dReal sdendmax);

// Emergency stopping (integrate forward the minimum acceleration)
dReal  EmergencyStop(Constraints& constraints, dReal sdbeg, Trajectory& restrajectory);

//////// ////////////////////////////////////////////////////////////
///////////////////////// Utilities ////////////////////////////////
////////////////////////////////////////////////////////////////////

// Various operations on Vectors
dReal VectorMin(const std::vector<dReal>& v);
dReal VectorMax(const std::vector<dReal>& v);
void VectorAdd(const std::vector<dReal>&a, const std::vector<dReal>&b,  std::vector<dReal>& res, dReal coefa=1, dReal coefb=1);
void VectorMultScalar(const std::vector<dReal>&a, std::vector<dReal>& res, dReal scalar);
dReal VectorNorm(const std::vector<dReal>&v);
void PrintVector(const std::vector<dReal>& v);

// Read a vector of dReal from a space-separated string
void VectorFromString(const std::string& s,std::vector<dReal>&resvect);

/// \brief read N items from a stream and put them into vector
void ReadVectorFromStream(std::istream& s, size_t N, std::vector<dReal>& resvect);

// Solve a0 + a1*x + a2*x^2 = 0 in the interval [lowerbound,upperbound]
// Return false if no solution in [lowerbound,upperbound], true otherwise
// If more than one solution, choose the smallest
bool SolveQuadraticEquation(dReal a0, dReal a1, dReal a2, dReal& sol, dReal lowerbound=-INF, dReal upperbound=INF);

// Check whether the point (s,sd) is above at least one profile in the list
bool IsAboveProfilesList(dReal s, dReal sd, std::list<Profile>&resprofileslist, bool searchbackward=false, dReal softborder=TINY2);

// Find the lowest profile at s
// Return false if no profile covers s, true otherwise
bool FindLowestProfile(dReal s, Profile& profile, dReal& tres, std::list<Profile>& resprofileslist);

void CheckInsert(std::list<std::pair<dReal,dReal> >& reslist, std::pair<dReal,dReal> e, bool inverse = false);
void FindMaxima(const std::list<std::pair<dReal,dReal> >& origlist, std::list<std::pair<dReal,dReal> >& reslist, bool inverse = false);

Profile StraightProfile(dReal sbackward,dReal sforward,dReal sdbackward,dReal sdforward);



}

#endif // TOPP_H
