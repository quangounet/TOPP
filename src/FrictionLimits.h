// -*- coding: utf-8 -*-
// Copyright (C) 2013 Quang-Cuong Pham <cuong.pham@normalesup.org>
//
// This file is part of the Time-Optimal Path Parameterization (TOPP) library.
// TOPP is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// at your option, any later version.
//
// This program is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with this program.  If not, see <http://www.gnu.org/licenses/>.


#include <openrave/openrave.h>
#include "openravepy.h"

#include "TOPP.h"

using namespace OpenRAVE;

namespace TOPP {

class FrictionLimits : public QuadraticConstraints {
public:
    FrictionLimits(RobotBasePtr probot, std::string& constraintsstring, Trajectory* ptraj);

    int ndof;
    int nlink;
    int nbottle;
    dReal mt; // tray's mass
    dReal mu;
    dReal eps;
    
    std::vector<dReal> mbvect; // bottles' masses
    std::vector<dReal> dxvect, dyvect, bottlehvect;
    std::vector<dReal> objspecs; // dx(half width), dy(hale depth), half height (of the bottle)
    std::vector<KinBody::LinkPtr> linksvector; // Vector of pointers to the links

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
