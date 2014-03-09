// -*- coding: utf-8 --*
// Copyright (C) 2014 Rosen Diankov <rosen.diankov@mujin.co.jp> & Cuong Pham
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
// \author Rosen Diankov
#include <openrave/openrave.h>

#define PY_ARRAY_UNIQUE_SYMBOL PyArrayHandle
#include <boost/python.hpp>
//#include <pyconfig.h>
//#include <numpy/arrayobject.h>
//#include <openrave/xmlreaders.h>

// declared from openravepy_int
namespace openravepy {

OpenRAVE::Transform ExtractTransform(const boost::python::object& oraw);
OpenRAVE::TransformMatrix ExtractTransformMatrix(const boost::python::object& oraw);
boost::python::object toPyArray(const OpenRAVE::TransformMatrix& t);
boost::python::object toPyArray(const OpenRAVE::Transform& t);

OpenRAVE::XMLReadablePtr ExtractXMLReadable(boost::python::object o);
boost::python::object toPyXMLReadable(OpenRAVE::XMLReadablePtr p);
bool ExtractIkParameterization(boost::python::object o, OpenRAVE::IkParameterization& ikparam);
boost::python::object toPyIkParameterization(const OpenRAVE::IkParameterization& ikparam);
boost::python::object toPyIkParameterization(const std::string& serializeddata);

boost::python::object toPyPlannerParameters(OpenRAVE::PlannerBase::PlannerParametersPtr params);

boost::python::object toPyEnvironment(boost::python::object);
boost::python::object toPyKinBody(OpenRAVE::KinBodyPtr, boost::python::object opyenv);
boost::python::object toPyKinBodyLink(OpenRAVE::KinBody::LinkPtr plink, boost::python::object opyenv);

//EnvironmentBasePtr GetEnvironment(boost::python::boost::python::object);
OpenRAVE::TrajectoryBasePtr GetTrajectory(boost::python::object);
OpenRAVE::KinBodyPtr GetKinBody(boost::python::object);
OpenRAVE::KinBody::LinkPtr GetKinBodyLink(boost::python::object);
OpenRAVE::RobotBasePtr GetRobot(boost::python::object o);

OpenRAVE::EnvironmentBasePtr GetEnvironment(boost::python::object o);
}
