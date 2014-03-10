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
#ifdef WITH_OPENRAVE

#include "TorqueLimitsRave.h"

using namespace OpenRAVE;

namespace TOPP {

TorqueLimitsRave::TorqueLimitsRave(RobotBasePtr probot, std::string& constraintsstring, Trajectory* ptraj){
    trajectory = *ptraj;
    int ndof = trajectory.dimension;
    std::istringstream iss(constraintsstring);
    iss >> discrtimestep;
    ReadVectorFromStream(iss, ndof, vmax);
    ReadVectorFromStream(iss, ndof, taumin);
    ReadVectorFromStream(iss, ndof, taumax);
    hasvelocitylimits = false;
    FOREACH(itv, vmax) {
        if( std::abs(*itv) > TINY ) {
            hasvelocitylimits = true;
            break;
        }
    }

    // Define the avect, bvect, cvect
    int ndiscrsteps = int((trajectory.duration+1e-10)/discrtimestep)+1;
    std::vector<dReal> q(ndof), qd(ndof), qdd(ndof), tmp0(ndof), tmp1(ndof), torquesimple;
    boost::array< std::vector< dReal >, 3 > torquecomponents;
    {
        avect.resize(ndiscrsteps);
        bvect.resize(ndiscrsteps);
        cvect.resize(ndiscrsteps);
        EnvironmentMutex::scoped_lock lock(probot->GetEnv()->GetMutex()); // lock environment
        for(int i = 0; i<ndiscrsteps; i++) {
            dReal s = i*discrtimestep;
            trajectory.Eval(s,q);
            trajectory.Evald(s,qd);
            trajectory.Evaldd(s,qdd);
            probot->SetDOFValues(q,KinBody::CLA_Nothing);
            probot->SetDOFVelocities(qd,KinBody::CLA_Nothing);
            probot->ComputeInverseDynamics(torquesimple,qd);
            probot->ComputeInverseDynamics(torquecomponents,qdd);
            VectorAdd(torquesimple,torquecomponents[1], bvect[i], 1, -1);
            VectorAdd(bvect[i], torquecomponents[2], avect[i], 1, -1);
            VectorAdd(torquecomponents[0],torquecomponents[1], bvect[i]);
            cvect[i] = torquecomponents[2];
        }
    }
}

void ConvertToTrajectory(OpenRAVE::TrajectoryBaseConstPtr pintraj, const OpenRAVE::ConfigurationSpecification& posspec, Trajectory& outtraj)
{
    int N = pintraj->GetNumWaypoints();
    if( N < 2 ) {
        throw TOPPException("openrave trajectory has less than 2 waypoints");
    }
    if( pintraj->GetDuration() <= 0 ) {
        throw TOPPException("openrave trajectory is not retimed");
    }
    OpenRAVE::ConfigurationSpecification timespec;
    timespec.AddDeltaTimeGroup();
    std::vector<OpenRAVE::ConfigurationSpecification::Group>::const_iterator itposgroup = pintraj->GetConfigurationSpecification().FindCompatibleGroup(posspec._vgroups.at(0).name, false);
    BOOST_ASSERT(itposgroup!= pintraj->GetConfigurationSpecification()._vgroups.end());

    int degree = 0;
    if( itposgroup->interpolation == "linear" ) {
        degree = 1;
    }
    else if( itposgroup->interpolation == "quadratic" ) {
        degree = 2;
    }
    else if( itposgroup->interpolation == "cubic" ) {
        degree = 3;
    }
    else if( itposgroup->interpolation == "quadric" ) {
        degree = 4;
    }
    else if( itposgroup->interpolation == "quintic" ) {
        degree = 5;
    }
    else {
        throw TOPP_EXCEPTION_FORMAT("unknown interpolation method '%s'", itposgroup->interpolation, 0);
    }

    // go up to accelerations
    std::vector<ConfigurationSpecification> derivspecs(std::min(degree,3));
    derivspecs[0] = posspec;
    for(size_t i = 1; i < derivspecs.size(); ++i) {
        derivspecs[i] = posspec.ConvertToDerivativeSpecification(i);
    }

    std::list<Chunk> listchunks;
    std::vector< std::vector<dReal> > vprevpoint(derivspecs.size()), vnewpoint(derivspecs.size());
    std::vector<dReal> vdeltatime;
    std::vector<Polynomial> polynomialsvector(posspec.GetDOF());
    std::vector<dReal> coefficientsvector(degree+1);
    for(size_t ideriv = 0; ideriv < derivspecs.size(); ++ideriv) {
        pintraj->GetWaypoint(0, vprevpoint[ideriv], derivspecs[ideriv]);
    }

    for(size_t iwaypoint = 1; iwaypoint < pintraj->GetNumWaypoints(); ++iwaypoint) {
        for(size_t ideriv = 0; ideriv < derivspecs.size(); ++ideriv) {
            pintraj->GetWaypoint(iwaypoint, vnewpoint[ideriv], derivspecs[ideriv]);
        }
        pintraj->GetWaypoint(iwaypoint, vdeltatime, timespec);
        dReal deltatime = vdeltatime.at(0);
        dReal ideltatime = 1/deltatime;
        for(int idof = 0; idof < posspec.GetDOF(); ++idof) {
            // try not to use the acceleration
            dReal px = vnewpoint.at(0).at(idof) - vprevpoint.at(0).at(idof);
            if( degree == 2 ) {
                coefficientsvector[0] = vprevpoint.at(0).at(idof);
                coefficientsvector[1] = vprevpoint.at(1).at(idof);
                coefficientsvector[2] = (px*ideltatime - vprevpoint.at(1).at(idof))*ideltatime;
            }
            else if( degree == 3 ) {
                coefficientsvector[0] = vprevpoint.at(0).at(idof);
                coefficientsvector[1] = vprevpoint.at(1).at(idof);
                coefficientsvector[2] = (3*px*ideltatime - 2*vprevpoint.at(1).at(idof) - vnewpoint.at(1).at(idof))*ideltatime;
                coefficientsvector[3] = (-2*px*ideltatime + vprevpoint.at(1).at(idof) + vnewpoint.at(1).at(idof))*ideltatime*ideltatime;
            }
            else {
                throw TOPP_EXCEPTION_FORMAT("do not degree %d yet", degree, 0);
            }
            polynomialsvector[idof].InitFromCoefficientsVector(coefficientsvector);
        }
        listchunks.push_back(Chunk(deltatime, polynomialsvector));
        vprevpoint.swap(vnewpoint);
    }
    outtraj.InitFromChunksList(listchunks);
}

TorqueLimitsRave2::TorqueLimitsRave2(RobotBasePtr probot, OpenRAVE::TrajectoryBaseConstPtr ptraj, dReal discrtimestep)
{
    EnvironmentMutex::scoped_lock lock(probot->GetEnv()->GetMutex()); // lock environment
    RobotBase::RobotStateSaver robotsaver(probot, KinBody::Save_LinkTransformation|KinBody::Save_LinkVelocities);
    OPENRAVE_ASSERT_OP((int)probot->GetActiveDOFIndices().size(),==,probot->GetActiveDOF()); // don't allow affine dofs
    this->discrtimestep = discrtimestep;
    ConvertToTrajectory(ptraj, probot->GetActiveConfigurationSpecification(), trajectory);
    int ndof = trajectory.dimension;
    probot->GetActiveDOFVelocityLimits(vmax);
    hasvelocitylimits = false;
    FOREACH(itv, vmax) {
        if( std::abs(*itv) > TINY ) {
            hasvelocitylimits = true;
            break;
        }
    }

    std::vector<dReal> vgearratios(probot->GetActiveDOF());
    taumin.resize(probot->GetActiveDOF());
    taumax.resize(probot->GetActiveDOF());
    for(int i = 0; i < probot->GetActiveDOF(); ++i) {
        KinBody::JointPtr pjoint = probot->GetJointFromDOFIndex(probot->GetActiveDOFIndices()[i]);
        dReal maxtorque = pjoint->GetMaxTorque(probot->GetActiveDOFIndices()[i] - pjoint->GetDOFIndex());
        // ElectricMotorActuatorInfo is a new spec in openrave
        ElectricMotorActuatorInfoPtr infoElectricMotor = pjoint->GetInfo()._infoElectricMotor;
        if( !!infoElectricMotor ) {
            dReal gear_ratio = infoElectricMotor->gear_ratio;
            //maxtorque =
        }
        else {
            std::cout << "could not find electric motor definition for joint " << pjoint->GetName() << std::endl;
        }
        taumax[i] = maxtorque;
        taumin[i] = -taumax[i];
    }

    // Define the avect, bvect, cvect
    int ndiscrsteps = int((trajectory.duration+1e-10)/discrtimestep)+1;
    std::vector<dReal> q(ndof), qd(ndof), qdd(ndof), vfullvalues(probot->GetDOF()), torquesimple;
    boost::array< std::vector< dReal >, 3 > torquecomponents;
    {
        avect.resize(ndiscrsteps);
        bvect.resize(ndiscrsteps);
        cvect.resize(ndiscrsteps);
        for(int i = 0; i<ndiscrsteps; i++) {
            dReal s = i*discrtimestep;
            trajectory.Eval(s,q);
            trajectory.Evald(s,qd);
            trajectory.Evaldd(s,qdd);
            probot->SetActiveDOFValues(q,KinBody::CLA_Nothing);
            probot->SetActiveDOFVelocities(qd,KinBody::CLA_Nothing);
            for(int idof = 0; idof < ndof; ++idof) {
                vfullvalues[probot->GetActiveDOFIndices()[idof]] = qd[idof];
            }
            probot->ComputeInverseDynamics(torquesimple,vfullvalues);
            for(int idof = 0; idof < ndof; ++idof) {
                vfullvalues[probot->GetActiveDOFIndices()[idof]] = qdd[idof];
            }
            probot->ComputeInverseDynamics(torquecomponents,vfullvalues);
            VectorAdd(torquesimple,torquecomponents[1], bvect[i], 1, -1);
            VectorAdd(bvect[i], torquecomponents[2], avect[i], 1, -1); // avect = torquesimple - torquecomponents[1] - torquecomponents[2]
            VectorAdd(torquecomponents[0],torquecomponents[1], bvect[i]); // bvect = torquecomponents[0] + torquecomponents[1]
            cvect[i] = torquecomponents[2];
        }
    }
}

}

#endif
