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

#include <boost/format.hpp>
#include <boost/python.hpp>
#include <boost/python/module.hpp>
#include <boost/python/def.hpp>

#ifdef WITH_OPENRAVE
#include "TorqueLimitsRave.h"
#include "openravepy.h"
#endif

using namespace boost::python;
using namespace TOPP;


class TOPPInstance {
public:
    TOPPInstance(object o, std::string problemtype, std::string constraintsstring, std::string trajectorystring){

        TOPP::Trajectory* ptrajectory = new TOPP::Trajectory(trajectorystring);

        if (problemtype.compare("KinematicLimits")==0) {
            pconstraints.reset(new KinematicLimits(constraintsstring));
            pconstraints->trajectory = *ptrajectory;
        }
        else if (problemtype.compare("TorqueLimits")==0) {
            pconstraints.reset(new TorqueLimits(constraintsstring));
            pconstraints->trajectory = *ptrajectory;
        }
        else if (problemtype.compare("QuadraticConstraints")==0) {
            pconstraints.reset(new QuadraticConstraints(constraintsstring));
            pconstraints->trajectory = *ptrajectory;
            pconstraints->CheckInput();
        }

#ifdef WITH_OPENRAVE
        else if (problemtype.compare("TorqueLimitsRave")==0) {
            _probot = openravepy::GetRobot(o);
            pconstraints.reset(new TorqueLimitsRave(_probot,constraintsstring,ptrajectory));
        }
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

    TOPPInstance(object orobot, std::string problemtype, object otrajectory, dReal discrtimestep)
    {
#ifdef WITH_OPENRAVE
        _probot = openravepy::GetRobot(orobot);
        OpenRAVE::TrajectoryBasePtr ptrajectory = openravepy::GetTrajectory(otrajectory);
        if (problemtype.compare("TorqueLimitsRave2")==0) {
            pconstraints.reset(new TorqueLimitsRave2(_probot,ptrajectory,discrtimestep));
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

    std::string restrajectorystring;
    std::string resextrastring;
    std::string resprofilesliststring;
    std::string switchpointsliststring;
    int ntangenttreated;
    int nsingulartreated;
    TOPP::dReal resduration;
    TOPP::dReal sdendmin, sdendmax;
    TOPP::dReal sdbegmin, sdbegmax;

    // Tuning parameters
    dReal integrationtimestep, reparamtimestep;
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


    int RunComputeProfiles(TOPP::dReal sdbeg, TOPP::dReal sdend){
        // Set tuning parameters
        pconstraints->integrationtimestep = integrationtimestep;
        pconstraints->passswitchpointnsteps = passswitchpointnsteps;
        pconstraints->extrareps = extrareps;

        int res = ComputeProfiles(*pconstraints,sdbeg,sdend);
        resduration = pconstraints->resduration;
        return res;
    }


    int ReparameterizeTrajectory(){
        // Set tuning parameters
        pconstraints->reparamtimestep = reparamtimestep;

        int ret = pconstraints->trajectory.Reparameterize(*pconstraints, restrajectory);
        return ret;
    }


    int RunVIP(TOPP::dReal sdbeg1, TOPP::dReal sdbeg2){
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

    int RunVIPBackward(TOPP::dReal sdend1, TOPP::dReal sdend2){
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
    object GetOpenRAVEResultTrajectory(object opyenv)
    {
        OpenRAVE::TrajectoryBasePtr ptraj = OpenRAVE::RaveCreateTrajectory(_probot->GetEnv());
        ConvertToOpenRAVETrajectory(restrajectory, ptraj, _probot->GetActiveConfigurationSpecification());
        return openravepy::toPyTrajectory(ptraj, opyenv);
    }

#endif

    void WriteResultTrajectory(){
        // std::stringstream ss; ss << std::setprecision(std::numeric_limits<OpenRAVE::dReal>::digits10+1);
        std::stringstream ss; ss << std::setprecision(11);
        // printf("WriteResultTrajectory: %d %f %d blah\n",
        //        restrajectory.dimension, restrajectory.duration,
        //        restrajectory.degree);
        restrajectory.Write(ss);
        restrajectorystring = ss.str();
    }

    void WriteProfilesList(){
        std::list<Profile>::iterator itprofile = pconstraints->resprofileslist.begin();
        //std::stringstream ss; ss << std::setprecision(std::numeric_limits<OpenRAVE::dReal>::digits10+1);
        std::stringstream ss; ss << std::setprecision(11);
        TOPP::dReal dt = 1e-4;
        pconstraints->WriteMVCBobrow(ss,dt);
        ss << "\n";
        pconstraints->WriteMVCDirect(ss,dt);
        ss << "\n";
        while(itprofile!=pconstraints->resprofileslist.end()) {
            itprofile->Write(ss,dt);
            ss << "\n";
            itprofile++;
        }
        resprofilesliststring = ss.str();
    }

    void WriteSwitchPointsList(){
        //std::stringstream ss; ss << std::setprecision(std::numeric_limits<OpenRAVE::dReal>::digits10+1);
        std::stringstream ss; ss << std::setprecision(11);
        std::list<SwitchPoint>::iterator itsw = pconstraints->switchpointslist.begin();
        while(itsw != pconstraints->switchpointslist.end()) {
            ss << itsw->s << " " << itsw->sd << " " << itsw->switchpointtype << "\n";
            itsw++;
        }
        switchpointsliststring = ss.str();
        ntangenttreated = pconstraints->ntangenttreated;
        nsingulartreated = pconstraints->nsingulartreated;
    }

    // Extra string, such as the coordinates of the ZMP (depending on the application)
    void WriteExtra(){
        //std::stringstream ss; ss << std::setprecision(std::numeric_limits<OpenRAVE::dReal>::digits10+1);
        std::stringstream ss; ss << std::setprecision(11);
        pconstraints->WriteExtra(ss);
        resextrastring = ss.str();

    }

};


BOOST_PYTHON_MODULE(TOPPbindings) {
    using namespace boost::python;
    class_<TOPPInstance>("TOPPInstance", init<object,std::string,std::string,std::string>())
    .def(init<object,std::string,object,dReal>(args("openraverobot","problemtype","openravetrajectory","discrtimestep")))
    .def_readwrite("integrationtimestep", &TOPPInstance::integrationtimestep)
    .def_readwrite("reparamtimestep", &TOPPInstance::reparamtimestep)
    .def_readwrite("passswitchpointnsteps", &TOPPInstance::passswitchpointnsteps)
    .def_readwrite("extrareps", &TOPPInstance::extrareps)
    .def_readonly("restrajectorystring", &TOPPInstance::restrajectorystring)
    .def_readonly("resprofilesliststring", &TOPPInstance::resprofilesliststring)
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
    .def("RunComputeProfiles",&TOPPInstance::RunComputeProfiles)
    .def("ReparameterizeTrajectory",&TOPPInstance::ReparameterizeTrajectory)
    .def("RunVIP",&TOPPInstance::RunVIP)
    .def("RunVIPBackward",&TOPPInstance::RunVIPBackward)
    .def("WriteResultTrajectory",&TOPPInstance::WriteResultTrajectory)
    .def("WriteProfilesList",&TOPPInstance::WriteProfilesList)
    .def("WriteExtra",&TOPPInstance::WriteExtra)
    .def("WriteSwitchPointsList",&TOPPInstance::WriteSwitchPointsList)
#ifdef WITH_OPENRAVE
    .def("GetOpenRAVEResultTrajectory", &TOPPInstance::GetOpenRAVEResultTrajectory, args("pyenv"))
#endif
    ;
}
