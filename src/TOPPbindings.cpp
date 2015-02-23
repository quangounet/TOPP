// -*- coding: utf-8 -*-
// Copyright (C) 2013 Quang-Cuong Pham <cuong.pham@normalesup.org> & Rosen Diankov <rosen.diankov@gmail.com>
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
#include "TOPP.h"
#include "KinematicLimits.h"
#include "TorqueLimits.h"
#include "PolygonConstraints.h"

#include <boost/format.hpp>
#include <boost/python.hpp>
#include <boost/python/module.hpp>
#include <boost/python/def.hpp>

#ifdef WITH_OPENRAVE
//#include "ComputePolygon.h"
#include "TorqueLimitsRave.h"
#include "TorqueLimitsRave3.h"
//#include "ZMPTorqueLimits.h"
#include "FrictionLimits.h"
#include "openravepy.h"
#endif

using namespace boost::python;
using namespace TOPP;


class TOPPInstance {
public:
    TOPPInstance(object o, std::string problemtype, std::string constraintsstring, std::string trajectorystring) {

        TOPP::Trajectory* ptrajectory = new TOPP::Trajectory(trajectorystring);

        if (problemtype.compare("KinematicLimits") == 0) {
            pconstraints.reset(new KinematicLimits(constraintsstring));
            pconstraints->trajectory = *ptrajectory;
        }
        else if (problemtype.compare("TorqueLimits") == 0) {
            pconstraints.reset(new TorqueLimits(constraintsstring));
            pconstraints->trajectory = *ptrajectory;
        }
        else if (problemtype.compare("PolygonConstraints") == 0) {
            pconstraints.reset(new PolygonConstraints(constraintsstring));
            pconstraints->trajectory = *ptrajectory;
        }
        else if (problemtype.compare("QuadraticConstraints") == 0) {
            pconstraints.reset(new QuadraticConstraints(constraintsstring));
            pconstraints->trajectory = *ptrajectory;
            pconstraints->CheckInput();
        }

#ifdef WITH_OPENRAVE
        else if (problemtype.compare("TorqueLimitsRave") == 0) {
            _probot = openravepy::GetRobot(o);
            pconstraints.reset(new TorqueLimitsRave(_probot, constraintsstring, ptrajectory));
        }
        else if (problemtype.compare("FrictionLimits") == 0) {
            _probot = openravepy::GetRobot(o);
            pconstraints.reset(new FrictionLimits(_probot, constraintsstring, ptrajectory));
        }
        // else if (problemtype.compare("ZMPTorqueLimits")==0) {
        //     _probot = openravepy::GetRobot(o);
        //     pconstraints.reset(new ZMPTorqueLimits(_probot,constraintsstring,ptrajectory));
        // }
#endif
        else {
            throw TOPP_EXCEPTION_FORMAT("cannot create %s problem type", problemtype, 0);
        }

        // Set default public tuning parameters
        integrationtimestep = 0;
        reparamtimestep = 0;
        passswitchpointnsteps = 5;
        extrareps = 0;

        // Set default private tuning parameters
        pconstraints->bisectionprecision = 0.01;
        pconstraints->loweringcoef = 0.95;
    }

    TOPPInstance(object orobot, std::string problemtype, object otrajectory, TOPP::dReal discrtimestep)
    {
#ifdef WITH_OPENRAVE
        _probot = openravepy::GetRobot(orobot);
        std::list<OpenRAVE::TrajectoryBasePtr> listtrajectories;
        OpenRAVE::TrajectoryBasePtr ptrajectory = openravepy::GetTrajectory(otrajectory);
        if( !!ptrajectory ) {
            listtrajectories.push_back(ptrajectory);
        }
        else {
            // might be a list of trajectories
            size_t num = len(otrajectory);
            for (size_t i = 0; i < num; ++i) {
                OpenRAVE::TrajectoryBasePtr ptrajectory = openravepy::GetTrajectory(otrajectory[i]);
                BOOST_ASSERT(!!ptrajectory);
                listtrajectories.push_back(ptrajectory);
            }
        }
        BOOST_ASSERT(listtrajectories.size()>0);
        if (problemtype.compare("TorqueLimitsRave2") == 0) {
            pconstraints.reset(new TorqueLimitsRave2(_probot,listtrajectories.front(),discrtimestep));
        }
        else if (problemtype.compare("TorqueLimitsRave3") == 0) {
            pconstraints.reset(new TorqueLimitsRave3(_probot,listtrajectories.front(),discrtimestep));
        }
        else {
            throw TOPP_EXCEPTION_FORMAT("cannot create %s problem type", problemtype, 0);
        }

        // Set default public tuning parameters
        integrationtimestep = 0;
        reparamtimestep = 0;
        passswitchpointnsteps = 5;
        extrareps = 0;

        // Set default private tuning parameters
        pconstraints->bisectionprecision = 0.01;
        pconstraints->loweringcoef = 0.95;
#else
        throw TOPPException("not compiled with openrave");
#endif
    }

    boost::shared_ptr<Constraints> pconstraints;
    TOPP::Trajectory restrajectory;

    std::string outconstraintstring;
    std::string restrajectorystring;
    std::string resextrastring;
    std::string resprofilesliststring;
    std::string resfinalprofilestring;
    std::string switchpointsliststring;
    int ntangenttreated;
    int nsingulartreated;
    TOPP::dReal resduration;
    TOPP::dReal sdendmin, sdendmax;
    TOPP::dReal sdbegmin, sdbegmax;

    // Tuning parameters
    TOPP::dReal integrationtimestep, reparamtimestep;
    int passswitchpointnsteps, extrareps;
    
#ifdef WITH_OPENRAVE
    OpenRAVE::RobotBasePtr _probot;
#endif
    
    TOPP::dReal GetAlpha(TOPP::dReal s, TOPP::dReal sd) {
        std::pair<TOPP::dReal, TOPP::dReal> sdd_lim = pconstraints->SddLimits(s, sd);
        return sdd_lim.first;
    }
    
    
    TOPP::dReal GetBeta(TOPP::dReal s, TOPP::dReal sd) {
        std::pair<TOPP::dReal, TOPP::dReal> sdd_lim = pconstraints->SddLimits(s, sd);
        return sdd_lim.second;
    }
    
    void GetABC(TOPP::dReal s) {
	std::vector<std::vector<TOPP::dReal> > res = pconstraints->GetABCConstraints(s);
        std::cout << res[0][0] << " " << res[0][1] << "\n";
        std::cout << res[1][0] << " " << res[1][1] << "\n";
        std::cout << res[2][0] << " " << res[2][1] << "\n";
    }
    
    int RunComputeProfiles(TOPP::dReal sdbeg, TOPP::dReal sdend) {
        // Set tuning parameters
        pconstraints->integrationtimestep = integrationtimestep;
        pconstraints->passswitchpointnsteps = passswitchpointnsteps;
        pconstraints->extrareps = extrareps;
        pconstraints->stepthresh = 0.01;

        int res = ComputeProfiles(*pconstraints,sdbeg,sdend);
        resduration = pconstraints->resduration;
        return res;
    }


    int ReparameterizeTrajectory(TOPP::dReal reparamtimestep=0)
    {
        // Set tuning parameters
        pconstraints->reparamtimestep = reparamtimestep;

        int ret = pconstraints->trajectory.Reparameterize(*pconstraints, restrajectory);
        return ret;
    }


    int RunVIP(TOPP::dReal sdbeg1, TOPP::dReal sdbeg2) {
        // Set tuning parameters
        pconstraints->integrationtimestep = integrationtimestep;
        pconstraints->passswitchpointnsteps = passswitchpointnsteps;
        pconstraints->extrareps = extrareps;

        sdbegmin = sdbeg1;
        sdbegmax = sdbeg2;
        int ret = VIP(*pconstraints, sdbegmin, sdbegmax,
                      sdendmin, sdendmax);
        if(ret == 0) {
            sdendmin = -1;
            sdendmax = -1;
        }
        return ret;
    }

    int RunVIPBackward(TOPP::dReal sdend1, TOPP::dReal sdend2) {
        // Set tuning parameters
        pconstraints->integrationtimestep = integrationtimestep;
        pconstraints->passswitchpointnsteps = passswitchpointnsteps;
        pconstraints->extrareps = extrareps;

        sdendmin = sdend1;
        sdendmax = sdend2;
        int ret = VIPBackward(*pconstraints, sdbegmin, sdbegmax,
                              sdendmin, sdendmax);
        if(ret == 0) {
            sdbegmin = -1;
            sdbegmax = -1;
        }
        return ret;
    }

#ifdef WITH_OPENRAVE
    object ExtractOpenRAVETrajectoryFromProfiles(object opyenv)
    {
        OpenRAVE::TrajectoryBasePtr ptraj = OpenRAVE::RaveCreateTrajectory(_probot->GetEnv());
        bool bsuccess = TOPP::ExtractOpenRAVETrajectoryFromProfiles(*pconstraints, 0, _probot->GetActiveConfigurationSpecification(), ptraj);
        if( !bsuccess ) {
            return object();
        }
        return openravepy::toPyTrajectory(ptraj, opyenv);
    }

    object GetOpenRAVEResultTrajectory(object opyenv)
    {
        OpenRAVE::TrajectoryBasePtr ptraj = OpenRAVE::RaveCreateTrajectory(_probot->GetEnv());
        ConvertToOpenRAVETrajectory(restrajectory, ptraj, _probot->GetActiveConfigurationSpecification());
        return openravepy::toPyTrajectory(ptraj, opyenv);
    }
#endif

    TOPP::dReal RunEmergencyStop(TOPP::dReal sdbeg) {
        // Set tuning parameters
        pconstraints->integrationtimestep = integrationtimestep;
        pconstraints->passswitchpointnsteps = passswitchpointnsteps;
        pconstraints->reparamtimestep = reparamtimestep;
        pconstraints->extrareps = extrareps;

        TOPP::dReal res = EmergencyStop(*pconstraints, sdbeg, restrajectory);
        return res;
    }

    void WriteResultTrajectory() {
        // std::stringstream ss; ss << std::setprecision(std::numeric_limits<OpenRAVE::dReal>::digits10+1);
        //std::stringstream ss; ss << std::setprecision(std::numeric_limits<OpenRAVE::dReal>::digits10+1);
        // printf("WriteResultTrajectory: %d %f %d blah\n",
        //        restrajectory.dimension, restrajectory.duration,
        //        restrajectory.degree);
        std::stringstream ss;
        ss << std::setprecision(17);
        restrajectory.Write(ss);
        restrajectorystring = ss.str();
    }

    void WriteProfilesList() {
        std::list<Profile>::iterator itprofile = pconstraints->resprofileslist.begin();
        //std::stringstream ss; ss << std::setprecision(std::numeric_limits<OpenRAVE::dReal>::digits10+1);
        //std::stringstream ss; ss << std::setprecision(std::numeric_limits<OpenRAVE::dReal>::digits10+1);
        std::stringstream ss;
        ss << std::setprecision(17);

        TOPP::dReal dt = 1e-4;
        pconstraints->WriteMVCBobrow(ss, dt);
        ss << "\n";
        pconstraints->WriteMVCDirect(ss, dt);
        ss << "\n";
        while (itprofile != pconstraints->resprofileslist.end()) {
            itprofile->Write(ss, dt);
            ss << "\n";
            itprofile++;
        }
        resprofilesliststring = ss.str();
    }
    
    void WriteFinalProfile() {
	TOPP::dReal ds = 1e-4;
	TOPP::dReal s = 0;
	TOPP::dReal dur = pconstraints->trajectory.duration;
	std::stringstream ss; ss << std::setprecision(std::numeric_limits<OpenRAVE::dReal>::digits10+1);
	std::string separator = "";
	
	ss << dur << " " << ds << "\n";
	while (s < dur) {
	    ss << separator << s;
	    separator = " ";
	    s += ds;
	}
	ss << separator << dur << "\n";
	s = 0;
	separator = "";
	while (s < pconstraints->trajectory.duration) {
	    TOPP::ProfileSample sample = TOPP::FindLowestProfileFast(s, 1e30, pconstraints->resprofileslist);
	    ss << separator << sample.sd;
	    separator = " ";
	    s += ds;
	}
	TOPP::ProfileSample sample = TOPP::FindLowestProfileFast(dur, 1e30, pconstraints->resprofileslist);
	ss << separator << sample.sd;
	resfinalprofilestring = ss.str();
    }

    std::string SerializeInputTrajectory() {
        //std::stringstream ss; ss << std::setprecision(std::numeric_limits<OpenRAVE::dReal>::digits10+1);
        std::stringstream ss;
        ss << std::setprecision(17);

        pconstraints->trajectory.Write(ss);
        return ss.str();
    }

    void WriteSwitchPointsList(){
        //std::stringstream ss; ss << std::setprecision(std::numeric_limits<OpenRAVE::dReal>::digits10+1);
        //std::stringstream ss; ss << std::setprecision(std::numeric_limits<OpenRAVE::dReal>::digits10+1);
        std::stringstream ss;
        ss << std::setprecision(17);
        std::list<SwitchPoint>::iterator itsw = pconstraints->switchpointslist.begin();
        while (itsw != pconstraints->switchpointslist.end()) {
            ss << itsw->s << " " << itsw->sd << " " << itsw->switchpointtype << "\n";
            itsw++;
        }
	// /// switch points on mvcbobrowlower
	// itsw = pconstraints->switchpointslistlower.begin();
	// while (itsw != pconstraints->switchpointslistlower.end()) {
        //     ss << itsw->s << " " << itsw->sd << " " << itsw->switchpointtype << "\n";
        //     itsw++;
        // }
        switchpointsliststring = ss.str();
        ntangenttreated = pconstraints->ntangenttreated;
        nsingulartreated = pconstraints->nsingulartreated;
    }

    // Write Constraints (currently works only for QuadraticConstraints)
    void WriteConstraints(){
        //std::stringstream ss; ss << std::setprecision(std::numeric_limits<OpenRAVE::dReal>::digits10+1);
        std::stringstream ss;
        ss << std::setprecision(17);

        pconstraints->WriteConstraints(ss);
        outconstraintstring = ss.str();
    }

    // Extra string, such as the coordinates of the ZMP (depending on the application)

    void WriteExtra(){
        //std::stringstream ss; ss << std::setprecision(std::numeric_limits<OpenRAVE::dReal>::digits10+1);
        std::stringstream ss;
        ss << std::setprecision(17);
        pconstraints->WriteExtra(ss);
        resextrastring = ss.str();
    }

};


BOOST_PYTHON_MEMBER_FUNCTION_OVERLOADS(ReparameterizeTrajectory_overloads, ReparameterizeTrajectory, 0, 1)

BOOST_PYTHON_MODULE(TOPPbindings) {
    using namespace boost::python;
    class_<TOPPInstance>("TOPPInstance", init<object,std::string,std::string,std::string>())
	.def(init<object,std::string,object,TOPP::dReal>(args("openraverobot","problemtype","openravetrajectory","discrtimestep")))
	.def_readwrite("integrationtimestep", &TOPPInstance::integrationtimestep)
	.def_readwrite("reparamtimestep", &TOPPInstance::reparamtimestep)
	.def_readwrite("passswitchpointnsteps", &TOPPInstance::passswitchpointnsteps)
	.def_readwrite("extrareps", &TOPPInstance::extrareps)
	.def_readonly("restrajectorystring", &TOPPInstance::restrajectorystring)
	.def_readonly("outconstraintstring", &TOPPInstance::outconstraintstring)
	.def_readonly("resprofilesliststring", &TOPPInstance::resprofilesliststring)
	.def_readonly("resfinalprofilestring", &TOPPInstance::resfinalprofilestring)
	.def_readonly("resextrastring", &TOPPInstance::resextrastring)
	.def_readonly("switchpointsliststring", &TOPPInstance::switchpointsliststring)
	.def_readonly("ntangenttreated", &TOPPInstance::ntangenttreated)
	.def_readonly("nsingulartreated", &TOPPInstance::nsingulartreated)
	.def_readonly("resduration", &TOPPInstance::resduration)
	.def_readonly("sdbegmin", &TOPPInstance::sdbegmin)
	.def_readonly("sdbegmax", &TOPPInstance::sdbegmax)
	.def_readonly("sdendmin", &TOPPInstance::sdendmin)
	.def_readonly("sdendmax", &TOPPInstance::sdendmax)
	.def_readonly("pconstraints", &TOPPInstance::pconstraints)
	.def("GetAlpha",&TOPPInstance::GetAlpha)
	.def("GetBeta",&TOPPInstance::GetBeta)
	.def("GetABC", &TOPPInstance::GetABC)
	.def("RunComputeProfiles",&TOPPInstance::RunComputeProfiles)
	.def("ReparameterizeTrajectory",&TOPPInstance::ReparameterizeTrajectory, ReparameterizeTrajectory_overloads(args("reparamtimestep")))
	.def("RunVIP",&TOPPInstance::RunVIP)
	.def("RunVIPBackward",&TOPPInstance::RunVIPBackward)
	.def("RunEmergencyStop",&TOPPInstance::RunEmergencyStop)
	.def("WriteResultTrajectory",&TOPPInstance::WriteResultTrajectory)
	.def("WriteProfilesList",&TOPPInstance::WriteProfilesList)
	.def("WriteFinalProfile", &TOPPInstance::WriteFinalProfile)
	.def("WriteExtra",&TOPPInstance::WriteExtra)
	.def("WriteConstraints",&TOPPInstance::WriteConstraints)
	.def("WriteSwitchPointsList",&TOPPInstance::WriteSwitchPointsList)
	.def("SerializeInputTrajectory", &TOPPInstance::SerializeInputTrajectory)
#ifdef WITH_OPENRAVE
	.def("ExtractOpenRAVETrajectoryFromProfiles", &TOPPInstance::ExtractOpenRAVETrajectoryFromProfiles, args("pyenv"), "extract openrave trajectory directly from profiles without reparameterizing")
	.def("GetOpenRAVEResultTrajectory", &TOPPInstance::GetOpenRAVEResultTrajectory, args("pyenv"), "resulting re-parameterized trajectory")
#endif
	;
}
