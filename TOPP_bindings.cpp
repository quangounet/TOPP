#include "PiecewisePolynomialTrajectory.h"
#include <boost/python.hpp>
#include <boost/python/module.hpp>
#include <boost/python/def.hpp>

class TOPPProblem {
public:
    TOPPProblem(std::string constraintsstring0,std::string trajectorystring0,std::string tuningsstring0){
        constraintsstring = constraintsstring0;
        trajectorystring = trajectorystring0;
        tuningsstring = tuningsstring0;
    }
    std::string constraintsstring;
    std::string trajectorystring;
    std::string tuningsstring;
    std::string restrajectorystring;
    void Solve(){
        restrajectorystring = trajectorystring + ".parameterized";
    }
};





BOOST_PYTHON_MODULE(TOPP)
{
    using namespace boost::python;
    class_<TOPPProblem>("TOPPProblem", init<std::string,std::string,std::string>())
    .def_readonly("trajectorystring", &TOPPProblem::trajectorystring)
    .def_readonly("restrajectorystring", &TOPPProblem::restrajectorystring)
    .def("Solve",&TOPPProblem::Solve);
}

