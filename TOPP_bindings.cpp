#include "PiecewisePolynomialTrajectory.h"
#include "KinematicLimits.h"


#include <boost/python.hpp>
#include <boost/python/module.hpp>
#include <boost/python/def.hpp>

using namespace TOPP;

class TOPPProblem {
public:
    TOPPProblem(std::string constraintsstring,std::string trajectorystring,std::string tuningsstring){
        constraints = KinematicLimits(constraintsstring);
        ptrajectory = new PiecewisePolynomialTrajectory(trajectorystring);
        tunings = Tunings(tuningsstring);
    }

    KinematicLimits constraints;
    PiecewisePolynomialTrajectory* ptrajectory;
    PiecewisePolynomialTrajectory restrajectory;
    std::list<Profile> resprofileslist;

    Tunings tunings;
    std::string restrajectorystring;
    std::string resprofilesliststring;
    dReal resduration;

    void Solve(){
        constraints.Preprocess(*ptrajectory,tunings);
        Profile profile;
        ComputeLimitingCurves(constraints,resprofileslist);
        IntegrateForward(constraints,0,1e-4,profile,1e5,resprofileslist);
        resprofileslist.push_back(profile);
        IntegrateBackward(constraints,ptrajectory->duration,1e-4,profile,1e5,resprofileslist);
        resprofileslist.push_back(profile);
        ptrajectory->Reparameterize(resprofileslist,tunings.reparamtimestep,restrajectory);
        resduration = restrajectory.duration;
    }

    void WriteResultTrajectory(){
        std::stringstream ss;
        restrajectory.Write(ss);
        restrajectorystring = ss.str();
    }

    void WriteProfilesList(){
        std::list<Profile>::iterator itprofile = resprofileslist.begin();
        std::stringstream ss;
        while(itprofile!=resprofileslist.end()) {
            itprofile->Write(ss);
            ss << "\n";
            itprofile++;
        }
        resprofilesliststring = ss.str();
    }


};





BOOST_PYTHON_MODULE(TOPP)
{
    using namespace boost::python;
    class_<TOPPProblem>("TOPPProblem", init<std::string,std::string,std::string>())
    .def_readonly("restrajectorystring", &TOPPProblem::restrajectorystring)
    .def_readonly("resprofilesliststring", &TOPPProblem::resprofilesliststring)
    .def_readonly("resduration", &TOPPProblem::resduration)
    .def("Solve",&TOPPProblem::Solve)
    .def("WriteResultTrajectory",&TOPPProblem::WriteResultTrajectory)
    .def("WriteProfilesList",&TOPPProblem::WriteProfilesList);

}



