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


#include <openrave-core.h>

#include "TOPP.h"

/* Undefine macros that Python.h defines to avoid redefinition warning.  */
// imported from:
// https://github.com/Slicer/VTK/commit/58e382c061ec870de18ef1516db3009c30ae6577
#undef _POSIX_C_SOURCE
#undef _POSIX_THREADS
#undef _XOPEN_SOURCE // Only defined by Python.h in python 2.3 and up

#define NO_IMPORT_ARRAY
#include "openrave/python/bindings/openravepy_int.h"

#include <boost/python.hpp>
#include <boost/python/module.hpp>
#include <boost/python/def.hpp>


using namespace OpenRAVE;

namespace TOPP {

class FrictionLimits : public QuadraticConstraints {
public:
    FrictionLimits(RobotBasePtr probot, std::string& constraintsstring, Trajectory* ptraj);

    int ndof;
    int nlink;
    dReal mb; // bottle's mass
    dReal mt; // tray's mass
    dReal eps;
    std::vector<dReal> objspecs; // dx(half width), dy(hale depth), half height (of the bottle)
    std::vector<dReal> activelinks;
    std::vector<KinBody::LinkPtr> linksvector; // Vector of pointers to the links
    std::vector<int> dofsvector;
    std::vector<dReal> mass;

    Vector ZMP(std::vector<dReal>& qfilled, std::vector<dReal>& qdfilled, std::vector<dReal>& qddfilled, bool withangularmomentum=false);
    // C = Mv
    Vector MatrixMultVector(const boost::multi_array<dReal, 2>& M, const std::vector<dReal>& v);
    Vector MatrixMultVector(const boost::multi_array<dReal, 2>& M, const Vector& v);

    boost::multi_array<dReal, 2> ExtractI(const RaveTransformMatrix<dReal>& H);
    boost::multi_array<dReal, 2> ExtractR(const RaveTransform<dReal>& H);
    Vector ExtractT(const RaveTransform<dReal>& H);

    // C = aA + bB
    void MatrixAdd(const boost::multi_array<dReal, 2>& A, const boost::multi_array<dReal, 2>& B, boost::multi_array<dReal, 2>& C, dReal coefA, dReal coefB);
    boost::multi_array<dReal, 2> MatrixTrans(const boost::multi_array<dReal, 2>& A);
    boost::multi_array<dReal, 2> MatricesMult3(const boost::multi_array<dReal, 2>& A, const boost::multi_array<dReal, 2>& B);
};



}
