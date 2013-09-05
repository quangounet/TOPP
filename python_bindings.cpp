#include "PiecewisePolynomialTrajectory.h"
#include <boost/python/module.hpp>
#include <boost/python/def.hpp>



BOOST_PYTHON_MODULE(TOPP_bindings)
{
    using namespace boost::python;
    def("greet", greet);
}
