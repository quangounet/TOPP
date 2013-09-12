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


#include "TOPP.h"
#include "KinematicLimits.h"
#include "TorqueLimits.h"


#include <boost/python.hpp>
#include <boost/python/module.hpp>
#include <boost/python/def.hpp>

using namespace TOPP;

class TOPPProblem {
public:
    TOPPProblem(std::string problemtype, std::string constraintsstring,std::string trajectorystring,std::string tuningsstring){
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
        pconstraints->WriteMVC(ss);
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
    class_<TOPPProblem>("TOPPProblem", init<std::string,std::string,std::string,std::string>())
    .def_readonly("restrajectorystring", &TOPPProblem::restrajectorystring)
    .def_readonly("resprofilesliststring", &TOPPProblem::resprofilesliststring)
    .def_readonly("resduration", &TOPPProblem::resduration)
    .def_readonly("sdendmin", &TOPPProblem::sdendmin)
    .def_readonly("sdendmax", &TOPPProblem::sdendmax)
    .def("RunPP",&TOPPProblem::RunPP)
    .def("RunVIP",&TOPPProblem::RunVIP)
    .def("WriteResultTrajectory",&TOPPProblem::WriteResultTrajectory)
    .def("WriteProfilesList",&TOPPProblem::WriteProfilesList);

}
