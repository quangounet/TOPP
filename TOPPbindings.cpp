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


#include "TOPP.h"
#include "KinematicLimits.h"
#include "TorqueLimits.h"


#include <boost/python.hpp>
#include <boost/python/module.hpp>
#include <boost/python/def.hpp>

using namespace TOPP;

class TOPPInstance {
public:
    TOPPInstance(std::string problemtype, std::string constraintsstring,std::string trajectorystring,std::string tuningsstring){
        if(problemtype.compare("KinematicLimits")==0) {
            pconstraints = new KinematicLimits(constraintsstring);
        }
        else if(problemtype.compare("TorqueLimits")==0) {
            pconstraints = new TorqueLimits(constraintsstring);
        }
        ptrajectory = new Trajectory(trajectorystring);
        tunings = Tunings(tuningsstring);
    }

    Constraints* pconstraints;
    Trajectory* ptrajectory;
    Trajectory restrajectory;
    std::list<Profile> resprofileslist;

    Tunings tunings;
    std::string restrajectorystring;
    std::string resprofilesliststring;
    dReal resduration;
    dReal sdendmin,sdendmax;



    int RunPP(dReal sdbeg, dReal sdend){
        int ret = PP(*pconstraints,*ptrajectory,tunings,sdbeg,sdend,restrajectory,resprofileslist);
        if(ret) {
            resduration = restrajectory.duration;
        }
        return ret;
    }

    void RunVIP(dReal sdbegmin, dReal sdbegmax){
        int ret = VIP(*pconstraints,*ptrajectory,tunings,sdbegmin,sdbegmax,sdendmin,sdendmax,resprofileslist);
        if(ret == 0) {
            sdendmin = -1;
            sdendmax = -1;
        }
    }


    void WriteResultTrajectory(){
        std::stringstream ss;
        restrajectory.Write(ss);
        restrajectorystring = ss.str();
    }

    void WriteProfilesList(){
        std::list<Profile>::iterator itprofile = resprofileslist.begin();
        std::stringstream ss;
        pconstraints->WriteMVCBobrow(ss);
        ss << "\n";
        pconstraints->WriteMVCCombined(ss);
        ss << "\n";
        while(itprofile!=resprofileslist.end()) {
            itprofile->Write(ss);
            ss << "\n";
            itprofile++;
        }
        resprofilesliststring = ss.str();
    }
};





BOOST_PYTHON_MODULE(TOPPbindings)
{
    using namespace boost::python;
    class_<TOPPInstance>("TOPPInstance", init<std::string,std::string,std::string,std::string>())
    .def_readonly("restrajectorystring", &TOPPInstance::restrajectorystring)
    .def_readonly("resprofilesliststring", &TOPPInstance::resprofilesliststring)
    .def_readonly("resduration", &TOPPInstance::resduration)
    .def_readonly("sdendmin", &TOPPInstance::sdendmin)
    .def_readonly("sdendmax", &TOPPInstance::sdendmax)
    .def("RunPP",&TOPPInstance::RunPP)
    .def("RunVIP",&TOPPInstance::RunVIP)
    .def("WriteResultTrajectory",&TOPPInstance::WriteResultTrajectory)
    .def("WriteProfilesList",&TOPPInstance::WriteProfilesList);

}
