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

class ZMPTorqueLimits : public QuadraticConstraints {
public:
    ZMPTorqueLimits(RobotBasePtr probot0, std::string& constraintsstring, Trajectory* ptraj);

    RobotBasePtr probot;
    int ndof; // number of active DOFs
    int nlink0;
    dReal totalmass;
    std::vector<dReal> taumin, taumax; // Torque limits
    std::vector<dReal> zmplimits; // xmin,xmax,ymin,ymax
    std::vector<dReal> activelinks;
    std::vector<KinBody::LinkPtr> linksvector; // Vector of pointers to the links
    std::vector<int> dofsvector; // Vector of indices of active dofs
    std::vector<dReal> mass;
    std::vector<dReal> qdefault;
    std::string supportfootlinkname;

    Vector COM(std::vector<dReal>& qfilled);
    Vector ZMP(std::vector<dReal>& qfilled, std::vector<dReal>& qdfilled, std::vector<dReal>& qddfilled, bool withangularmomentum=false);
    void Fill(const std::vector<dReal>&q, std::vector<dReal>&qfilled);
    void Trim(const std::vector<dReal>&q, std::vector<dReal>&qtrimmed);

    void WriteExtra(std::stringstream& ss);

    // Multiply a matrix and a vector taking into account activedofs
    Vector MatrixMultVector(const boost::multi_array<dReal,2>& M, const std::vector<dReal>& v);

    /// \brief Get the transpose of the contact (rotational and translational)
    /// Jacobian for a given foot.
    ///
    /// \param[in] foot LinkPtr to the targeted foot
    /// \param[out] mjacobiantrans Transpose of the contact Jacobian matrix.
    void GetFootJacobianTranspose(const KinBody::LinkPtr foot,
                                  boost::multi_array<dReal,2>& mjacobiantrans);

    /// \brief Compute the inverse dynamics (actuated joint torques + contact
    /// wrench) for the humanoid robot in a single-support setting. The support
    /// foot is assumed to be the one with lowest z-coordinate.
    ///
    /// \param[out] doftorques The output actuated joint torques (base link
    /// coordinates will be zero).
    /// \param[in] dofaccelerations The dof accelerations of the current robot state.
    void ComputeInverseDynamicsSingleSupport(std::vector<dReal>& doftorques,
                                             const std::vector<dReal>& dofaccelerations);

    /// \brief Compute the separated inverse dynamics
    /// for the humanoid robot in a single-support setting. The support
    /// foot is assumed to be the one with lowest z-coordinate. Actuated
    /// torques are such that:
    ///
    /// torques = M(dofvalues) * dofaccel + C(dofvalues,dofvel) * dofvel +
    ///           G(dofvalues) - J(dofvalues)^T * Wext
    ///
    /// with:
    ///
    /// torques -- generalized forces associated with dofvalues
    /// M -- manipulator inertia tensor
    /// C -- coriolis and centripetal effects
    /// G -- gravity loading + base link angular acceleration contribution
    /// J -- contact Jacobian at the support foot
    /// Wext -- contact wrench exerted at the support foot
    ///
    /// \param[out] doftorquecomponents A set of 3 torques:
    ///    [M(dofvalues) * dofaccel,
    ///     C(dofvalues,dofvel) * dofvel,
    ///     G(dofvalues) - J(dofvalues)^T * Wext]
    /// coordinates will be zero).
    /// \param[in] dofaccelerations The dof accelerations of the current robot state.
    void ComputeInverseDynamicsSingleSupport(boost::array<std::vector<dReal>, 3>&
                                             doftorquecomponents, const std::vector<dReal>& dofaccelerations);
};


void MatrixAdd(const boost::multi_array<dReal,2>& A, const boost::multi_array<dReal,2>& B, boost::multi_array<dReal,2>& C, dReal coefA=1, dReal coefB=1);

}
