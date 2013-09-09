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


    void RunPP(dReal sdbeg, dReal sdend){
        int res = PP(*pconstraints,*ptrajectory,tunings,sdbeg,sdend,restrajectory,resprofileslist);
        resduration = restrajectory.duration;
    }

    void RunVIP(dReal sdbegmin, dReal sdbegmax){
        int res = VIP(*pconstraints,*ptrajectory,tunings,sdbegmin,sdbegmax,sdendmin,sdendmax,resprofileslist);
        if(res == 0) {
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
