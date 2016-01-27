#include "TOPP.h"
#include <boost/python.hpp>

namespace TOPP {
void ComputeSO3Constraints(const std::string& SO3trajstring, const std::string& constraintsstring, boost::python::list& resstringlist);

boost::multi_array<dReal, 2> SkewFromVect(const std::vector<dReal>& vect);

std::vector<dReal> Cross(const std::vector<dReal>& a, const std::vector<dReal>& b );

dReal DotVect(const std::vector<dReal>& a, const std::vector<dReal>& b );

std::vector<dReal> AddVect(const std::vector<dReal>& a, const std::vector<dReal>& b,dReal coefa, dReal coefb  );

boost::multi_array<dReal, 2> MatrixAdd(const boost::multi_array<dReal, 2>& A, const boost::multi_array<dReal, 2>& B, dReal coefA, dReal coefB);

boost::multi_array<dReal, 2> MatricesMult3(const boost::multi_array<dReal, 2>& A, const boost::multi_array<dReal, 2>& B);

std::vector<dReal> MatrixMultVector(const boost::multi_array<dReal, 2>& M, const std::vector<dReal>& v);

boost::multi_array<dReal, 2> Eye();

std::string VectToString(const std::vector<dReal>& vect, const bool& IsCvect);
}
