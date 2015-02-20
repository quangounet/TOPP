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
#include <cstdio>
#include <cstdlib>
#include <cassert>

#include "Errors.h"
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
#define BOBROWEXCLUDENOTDEFINED -1
    
////////////////////////////////////////////////////////////////////
/////////////////////////// Exception //////////////////////////////
////////////////////////////////////////////////////////////////////
    
class TOPPException : public std::exception {
    /// \brief Exception that all OpenRAVE internal methods throw; the error codes are held in \ref OpenRAVEErrorCode.
    public:
    TOPPException() : std::exception(), _s("unknown exception"), _errorcode(0) {
	}
    TOPPException(const std::string& s, int errorcode=0) : std::exception() {
	    _errorcode = errorcode;
	    _s = "topp (";
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
#define TOPP_EXCEPTION_FORMAT(s, args, errorcode) TOPP::TOPPException(boost::str(boost::format("[%s:%d] " s)%(__PRETTY_FUNCTION__)%(__LINE__)%args),errorcode)
    
    
////////////////////////////////////////////////////////////////////
//////////////////////// Switch Point //////////////////////////////
////////////////////////////////////////////////////////////////////
    
enum SwitchPointType {
    SP_TANGENT = 0,
    SP_SINGULAR = 1,
    SP_DISCONTINUOUS =2,
    SP_ZLAJPAH = 3
};

class SwitchPoint {
 public:
    SwitchPoint(dReal s0, dReal sd0, int switchpointtype0) {
	s = s0;
	sd = sd0;
	switchpointtype = switchpointtype0;
    }
    dReal s, sd;                     ///< position of the switchpoint
    int switchpointtype;             ///< type
    std::vector<dReal> slopesvector; ///< slope of the profiles near the singular switchpoint
};

////////////////////////////////////////////////////////////////////
/////////////////////////// Profile ////////////////////////////////
////////////////////////////////////////////////////////////////////

class Profile {
    /// \brief Velocity profile obtained from integrations
 public:
    Profile(const std::list<dReal>& slist, const std::list<dReal>& sdlist, const std::list<dReal>& sddlist, dReal integrationtimestep, bool forward);
    Profile(const std::vector<dReal>& svect, const std::vector<dReal>& sdvect, const std::vector<dReal>& sddvect, dReal integrationtimestep, bool forward);
    Profile() {
    }

    inline void Reset() {
	svect.resize(0);
	sdvect.resize(0);
	sddvect.resize(0);
	nsteps = 0;
    }

    std::vector<dReal> svect, sdvect, sddvect; ///< points on the integrated profile. 
    /// svect is in increasing order if it is integrated forward, decreasing order otherwise
    
    dReal integrationtimestep;
    int nsteps;     ///< nsteps = svect.size() -1
    dReal duration; ///< duration = nsteps*integrationtimestep
    bool forward;   ///< if true, the profile is from forward integration, otherwise, it is from backward integration
    
    bool FindTimestepIndex(dReal t, int &index, dReal &remainder) const;
    bool Invert(dReal s, dReal &t) const; 
    /// Finds t such that Eval(t) = s, return false if no such t
    
    dReal Eval(dReal t) const;
    dReal Evald(dReal t) const;
    dReal Evaldd(dReal t) const;
    void Print() const;
    void Write(std::stringstream &ss, dReal dt = 0.01) const;

};

////////////////////////////////////////////////////////////////////
///////////////////////// Constraints //////////////////////////////
////////////////////////////////////////////////////////////////////

class Constraints {
 public:
    
    /// \param Tuning parameters
    dReal discrtimestep;       ///< timestep to discretize the trajectory when computing the MVC
    dReal integrationtimestep; ///< timestep to integrate
    dReal reparamtimestep;     ///< timestep to reparamterize the trajectory based on the optimal profile
    int passswitchpointnsteps; ///< number of steps to integrate around a switch point
    int extrareps;             ///< number of times to lower integrationtimestep when integration fails
    dReal bisectionprecision;  ///< precision when determining the sd for a tangent switch point
    dReal loweringcoef;        ///< lower bound that multiplies sd when searching an sd for a tangent switch point

    /// \param MVC
    /// now modified to handle the case which sd has one valid interval [sdmin, sdmax] where sdmin > 0 (having some island)
    std::vector<dReal> discrsvect;   ///< discretized points on the s-axis
    int ndiscrsteps;                 ///< ndiscrsteps = discrvect.size() - 1
    std::vector<dReal> mvcbobrow;    ///< normal MVC above which has an infeasible area
    std::vector<dReal> mvcbobrowlower;   ///< special MVC below which has an infeasible area
    std::vector<dReal> mvccombined;
    bool hasvelocitylimits;
    std::vector<dReal> vmax;

    std::list<SwitchPoint> switchpointslist;         ///< list of switch points, ordered by s (ascending)
    std::list<SwitchPoint> switchpointslistlower;    ///< list of switch points on the lower MVC, ordered by s (ascending)
    std::list<std::pair<dReal, dReal> > zlajpahlist; ///< list of zlajpah points
    
    std::list<Profile> resprofileslist; ///< resulting profiles from integrations
    dReal resduration;                  ///< duration of the resulting profile    
    Trajectory trajectory;
    
    dReal stepthresh; ///< threshold for amount of s to step around a singularity.
    /// larger values stabilizes the graph, but can make trajectory slower.
    
    ////////////////////////////// General //////////////////////////////
    
    Constraints() {
	_busingcache = false;
	stepthresh = 0.01;
    }

    virtual void CheckInput() {
	/// Checks input after this->trajectory has been set (from TOPPbindings)
    }  
    
    virtual bool Preprocess();
    /// Computes the MVCs and the switchpoints and other initializations
    
    virtual void Discretize();
    /// Discretizes the time interval
    
    virtual void ComputeMVCBobrow();
    /// Computes the MVC given by acceleration constraints

    virtual void ComputeMVCBobrow2();
    /// Computes the lower MVC and the upper MVC given by acceleration constraints
        
    virtual void ComputeMVCCombined();
    /// Computes the combined MVC (incorporating pure velocity constraints)
    
    virtual void WriteMVCBobrow(std::stringstream& ss, dReal dt = 0.01);
    virtual void WriteMVCDirect(std::stringstream& ss, dReal dt = 0.01);
    virtual void WriteExtra(std::stringstream& ss) {
        return;
    }
    virtual void WriteConstraints(std::stringstream& ss) {
        return;
    }

    virtual std::vector<std::vector<dReal> > GetABCConstraints(dReal s) {
	std::cout << "Virtual method not implemented\n";
        throw TOPPException("Virtual method not implemented");
    }
    
    virtual dReal Interpolate1D(dReal s, const std::vector<dReal>& v);
    /// Linearly interpolates between two points in 1D
    
    ////////////////////////////// Limits //////////////////////////////
    
    virtual dReal SdLimitBobrowInit(dReal s) {
        /// Computes an upper limit on sd given by acceleration constraints (Bobrow)
	/// Called at initialization
	std::cout << "Virtual method not implemented\n";
        throw TOPPException("Virtual method not implemented");
    }

    virtual dReal SdLimitBobrow(dReal s);
    /// Computes an upper limit on sd given by acceleration constraints (Bobrow)
    /// Called after mvcbobrow exists
    
    dReal SdLimitBobrowExclude(dReal s, int iexclude) {
        /// Computes an upper limit on sd given by acceleration constraints (except the iexclude-th constraint)
	return BOBROWEXCLUDENOTDEFINED;
    }
    
    virtual dReal SdLimitCombinedInit(dReal s);
    /// Computes an upper limit on sd after incorporating pure velocity constraints

    virtual dReal SdLimitCombined(dReal s);
    
    virtual std::pair<dReal,dReal> SddLimits(dReal s, dReal sd) {
	/// Computes lower and upper limits on sdd
	/// Returns (\alpha, \beta)
        std::cout << "Virtual method not implemented\n";
        throw TOPPException("Virtual method not implemented");
    }
    
    ///
    virtual dReal SdLimitBobrowInitUpper(dReal s) {
        /// Computes an upper limit on sd given by acceleration constraints (Bobrow)
	/// Called at initialization
	std::cout << "Virtual method not implemented\n";
        throw TOPPException("Virtual method not implemented");
    }
    virtual dReal SdLimitBobrowInitLower(dReal s) {
        /// Computes a lower limit on sd given by acceleration constraints (Bobrow)
	/// Called at initialization
	std::cout << "Virtual method not implemented\n";
        throw TOPPException("Virtual method not implemented");
    }
    virtual dReal SdLimitBobrowExcludeUpper(dReal s, int iexclude) {
        /// Computes an upper limit on sd given by acceleration constraints (except the iexclude-th constraint)
	return BOBROWEXCLUDENOTDEFINED;
    }
    virtual dReal SdLimitBobrowExcludeLower(dReal s, int iexclude) {
        /// Computes a lower limit on sd given by acceleration constraints (except the iexclude-th constraint)
	return BOBROWEXCLUDENOTDEFINED;
    }
    virtual dReal SdLimitBobrowUpper(dReal s);                      ///< Computes sd at s on mvcbobrow
    virtual dReal SdLimitBobrowLower(dReal s);                      ///< Comptues sd at s on mvcbobrowlower
    ///
    
    virtual dReal SddLimitAlpha(dReal s, dReal sd) {
        return SddLimits(s, sd).first;
    }
    
    virtual dReal SddLimitBeta(dReal s, dReal sd) {
        return SddLimits(s, sd).second;
    }

    ////////////////////////////// Switch Points //////////////////////////////
    
    int ntangenttreated;  ///< number of tangent switchpoints that are treated (for stats purpose)
    int nsingulartreated; ///< number of singular switchpoints that are treated (for stats purpose)
    
    virtual void FindSwitchPoints();
    virtual void FindTangentSwitchPoints();
    virtual void FindDiscontinuousSwitchPoints();
    virtual void FindSingularSwitchPoints() {
        std::cout << "Virtual method not implemented\n";
        throw TOPPException("Virtual method not implemented");
    }

    virtual void TrimSwitchPoints();
    /// Replaces switch points that are close to each other by a single switch point
    /// Modifies switchpointslist
    
    virtual void ComputeSlopeDynamicSingularity(dReal s, dReal sd, std::vector<dReal>& slopesvector) {
	/// Computes the slope of the profiles near a dynamic singularity at (s, sd)
        std::cout << "Virtual method not implemented\n";
        throw TOPPException("Virtual method not implemented");
    }

    virtual void AddSwitchPoint(int i, int switchpointtype, dReal sd = -1);
    virtual void AddSwitchPoint2(dReal s, dReal sd, int switchpointtype);
    virtual void AddSwitchPointLower(dReal s, dReal sd, int switchpointtype);

    std::vector<dReal> _svectcache, _sdvectcache, _sddvectcache; ///< caches
    bool _busingcache;
    
};

////////////////////////////////////////////////////////////////////
///////////////////// Quadratic Constraints ////////////////////////
////////////////////////////////////////////////////////////////////

class QuadraticConstraints : public Constraints {
    /// \brief Constraints of the form a(s)sdd + b(s)sd^2 + c(s) <= 0
    /// See our article http://arxiv.org/abs/1312.6533 for more details
 public:
    QuadraticConstraints() : Constraints() {
    }
    QuadraticConstraints(const std::string& constraintsstring);
    
    ////////////////////////////// Overloaded Methods //////////////////////////////
    dReal SdLimitBobrowInit(dReal s);
    dReal SdLimitBobrowExclude(dReal s, int iexclude);
    std::pair<dReal, dReal> SddLimits(dReal s, dReal sd);
    void CheckInput();
    void ComputeSlopeDynamicSingularity(dReal s, dReal sd, std::vector<dReal>& slopesvector);
    void WriteConstraints(std::stringstream& ss);
    
    dReal SdLimitBobrowInitUpper(dReal s);                  ///< Computes mvcbobrow
    dReal SdLimitBobrowInitLower(dReal s);                  ///< Computes mvcbobrowlower
    dReal SdLimitBobrowExcludeUpper(dReal s, int iexclude);
    dReal SdLimitBobrowExcludeLower(dReal s, int iexclude);
    
    void FindSingularSwitchPoints();                        ///< Finds all singular switch points
    void FindTangentSwitchPoints();                         ///< Finds all tangent switch points
    void FindDiscontinuousSwitchPoints();                   ///< Finds all discontinuous switch points
    /// Finds all switch points on both mvcbobrow & mvcbobrowlower

    ////////////////////////////// Specific Members & Methods //////////////////////////////
    int nconstraints;  ///< number of constraints
    bool hasislands;
    // int nislands;      ///< number of islands
    std::vector<std::vector<dReal> > avect, bvect, cvect;
    
    void InterpolateDynamics(dReal s, std::vector<dReal>& a, std::vector<dReal>& b, std::vector<dReal>& c);
    /// Linearly interpolates the dynamics coefficients a, b, c
    
    virtual void FixStart(dReal& sstartnew, dReal& sdstartnew, dReal timestep);
    void FixEnd(dReal& sendnew, dReal& sdendnew);

};

////////////////////////////////////////////////////////////////////
////////////////////////// Integration /////////////////////////////
////////////////////////////////////////////////////////////////////

enum IntegrationReturnType {
    INT_MAXSTEPS = 0,  ///< reached maximum number of steps
    INT_END = 1,       ///< reached s = send (FW) or s = 0 (BW)
    INT_BOTTOM = 2,    ///< crossed sd = 0
    INT_MVC = 3,       ///< crossed the MVC
    INT_PROFILE = 4    ///< crossed a profile
};

enum CLCReturnType {
    CLC_OK = 0,        ///< no problem
    CLC_SWITCH = 1,    ///< could not pass a switchpoint
    CLC_BOTTOM = 2     ///< crossed sd = 0
};
 
int IntegrateForward(Constraints& constraints, dReal sstart, dReal sdstart, dReal dt, Profile& resprofile, 
		     int maxsteps = 1e7, bool testaboveexistingprofiles = true, bool testmvc = true, bool zlajpah = false);
/// Integrates forward from (sstart, sdstart)
/// \param[in] dt the integration timestep
/// \param[in] maxsteps the maximum steps to integrate for

int IntegrateBackward(Constraints& constraints, dReal sstart, dReal sdstart, dReal dt, Profile& resprofile, 
		      int maxsteps = 1e7, bool testaboveexistingprofiles = true, bool testmvc = true, bool zlajpah = false);
/// Integrate backward from (sstart,sdstart)
/// \param[in] dt the integration timestep
/// \param[in] maxsteps the maximum steps to integrate for

int ComputeLimitingCurves(Constraints& constraints);
/// Computes the CLC

int ComputeProfiles(Constraints& constraints, dReal sdbeg, dReal sdend);
/// Compute all the profiles (CLC, forward from 0, backward from send)

int VIP(Constraints& constraints, dReal sdbegmin, dReal sdbegmax, dReal& sdendmin, dReal& sdendmax);
int VIPBackward(Constraints& constraints, dReal& sdbegmin, dReal& sdbegmax, dReal sdendmin, dReal sdendmax);
///? TO DO:: change the names to AVPFW and AVPBW

dReal  EmergencyStop(Constraints& constraints, dReal sdbeg, Trajectory& restrajectory);
/// Does emergency stopping (integrate forward the minimum acceleration)

////////////////////////////////////////////////////////////////////
/////////////////////////// Utilities //////////////////////////////
////////////////////////////////////////////////////////////////////

/// \brief Various operations on std::vector<dReal>
dReal VectorMin(const std::vector<dReal>& v);
dReal VectorMax(const std::vector<dReal>& v);
void VectorAdd(const std::vector<dReal>& a, const std::vector<dReal>& b,  std::vector<dReal>& res, dReal coefa = 1, dReal coefb = 1);
void VectorMultScalar(const std::vector<dReal>& a, std::vector<dReal>& res, dReal scalar);
dReal VectorNorm(const std::vector<dReal>& v);
void PrintVector(const std::vector<dReal>& v);

void VectorFromString(std::string& s,std::vector<dReal>&resvect);
/// Reads an std::vector<dReal> from a space-separated string

void ReadVectorFromStream(std::istream& s, size_t N, std::vector<dReal>& resvect);
/// Reads N items from a stream istream and put them into resvect

bool SolveQuadraticEquation(dReal a0, dReal a1, dReal a2, dReal& sol, dReal lowerbound = -INF, dReal upperbound = INF, dReal epsilon = TINY);
/// Solves a0 + a1*x + a2*x^2 = 0 in the interval [lowerbound,upperbound]
/// Returns false if no solution in [lowerbound,upperbound], true otherwise
/// If more than one solution, the smallest is chosen
/// epsilon is used to threshold the coefficients and results to avoid dividing by 0

Profile StraightProfile(dReal sbackward, dReal sforward, dReal sdbackward, dReal sdforward);

bool FindLowestProfile(dReal s, Profile& profile, dReal& tres, const std::list<Profile>& resprofileslist);
/// Find the lowest profile (in terms of sd) at s
/// Return false if no profile covers s, true otherwise

class ProfileSample {
    /// \brief A location in the profile
 public:
 ProfileSample() : s(0), sd(0), sdd(0), sindex(0), t(0) {
    }
    std::list<Profile>::const_iterator itprofile;
    dReal s, sd, sdd; ///< samples of svect, sdvect, and sddvector
    int sindex;       ///< index into svect such that svect[sindex] <= s <= svect[sindex+1]
    dReal t;          ///< time from svect[index] to s
};

ProfileSample FindLowestProfileFast(dReal scur, dReal sdmax, const std::list<Profile>& resprofileslist);
/// Finds the lowest profile (in terms of sd) at s and returns a ProfileSample
/// If there's no profile that covers s, profilesample.itprofile = resprofileslist.end() (which means invalid)
/// \param scur the s where to look for a low sd
/// \param sdmax need to look for sd that are <= this value

ProfileSample FindEarliestProfileIntersection(dReal sstart, dReal sdstart, dReal sddstart, dReal tmax, const std::list<Profile>& resprofileslist, 
					      std::list<Profile>::const_iterator itprofileexclude, dReal& tintersect);
/// Finds if any profiles in resprofileslist intersects with a ramp [s, sd, sdd] during times [0, tmax].
/// \param tintersect The time to the intersection point computed as sstart + tintersect*(sdstart + tintersect*0.5*sddstart)
/// \param itprofileexclude profile to exclude for intersection, usually the profile that [sstart, sdstart, sddstart] come from.
/// \return the intersection sample point

bool IsAboveProfilesList(dReal s, dReal sd, const std::list<Profile>& resprofileslist, bool searchbackward = false, dReal softborder = TINY2);
/// Checks whether the point (s, sd) is above any profile in resprofileslist

void CheckInsert(std::list<std::pair<dReal,dReal> >& reslist, std::pair<dReal,dReal> e, bool inverse = false);
void FindMaxima(const std::list<std::pair<dReal,dReal> >& origlist, std::list<std::pair<dReal,dReal> >& reslist, bool inverse = false);

}

#endif /// TOPP_H
