// -*- coding: utf-8 -*-
// Copyright (C) 2013 Quang-Cuong Pham <cuong.pham@normalesup.org> & Rosen Diankov <rosen.diankov@gmail.com>
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

#include "TorqueLimitsRave3.h"
#include "TorqueLimitsRave.h"

using namespace OpenRAVE;

namespace TOPP {

TorqueLimitsRave3::TorqueLimitsRave3(RobotBasePtr probot, OpenRAVE::TrajectoryBaseConstPtr ptraj, dReal discrtimestep)
{
    EnvironmentMutex::scoped_lock lock(probot->GetEnv()->GetMutex()); // lock environment
    _probot = probot;
    RobotBase::RobotStateSaver robotsaver(probot, KinBody::Save_LinkTransformation|KinBody::Save_LinkVelocities);
    OPENRAVE_ASSERT_OP((int)probot->GetActiveDOFIndices().size(),==,probot->GetActiveDOF()); // don't allow affine dofs
    this->discrtimestep = discrtimestep;
    ConvertToTOPPTrajectory(ptraj, probot->GetActiveConfigurationSpecification(), trajectory);
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
        int iaxis = probot->GetActiveDOFIndices()[i] - pjoint->GetDOFIndex();
        dReal maxtorque = pjoint->GetMaxTorque(iaxis);
        dReal maxinertia = pjoint->GetMaxInertia(iaxis);
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
    if( ndiscrsteps > 1 ) {
        discrtimestep = trajectory.duration/(ndiscrsteps-1);
    }
    std::vector<dReal> q(ndof), qd(ndof), qdd(ndof), vfullvalues(probot->GetDOF()), torquesimple;
    probot->GetDOFValues(vfullvalues);
    boost::array< std::vector< dReal >, 3 > torquecomponents;
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
        // Torque components with qdd
        for(int idof = 0; idof < ndof; ++idof) {
            vfullvalues[probot->GetActiveDOFIndices()[idof]] = qdd[idof];
        }
        probot->ComputeInverseDynamics(torquecomponents,vfullvalues);
        // Torque simple with qd
        for(int idof = 0; idof < ndof; ++idof) {
            vfullvalues[probot->GetActiveDOFIndices()[idof]] = qd[idof];
        }
        probot->ComputeInverseDynamics(torquesimple,vfullvalues);
        // Note that the a,b,c in QuadraticConstraints are not the same as those in TorqueLimits
        // See TOPPopenravepy.ComputeTorquesConstraints for details
        avect[i].resize(ndof*2);
        bvect[i].resize(ndof*2);
        cvect[i].resize(ndof*2);
        // Constraints tau < taumax
        for(int idof = 0; idof < ndof; ++idof) {
            dReal a = torquesimple[idof] - torquecomponents[1][idof] - torquecomponents[2][idof];
            dReal b = torquecomponents[0][idof] + torquecomponents[1][idof];
            dReal c = torquecomponents[2][idof];
            if( avect[i][idof] >= 0 ) {
                // upper limits
                avect[i][idof] = a;
                bvect[i][idof] = b;
                cvect[i][idof] = c - taumax[idof];
                // lower limits
                avect[i][ndof+idof] = -a;
                bvect[i][ndof+idof] = -b;
                cvect[i][ndof+idof] = -c + taumin[idof];
            }
            else {
                // upper limits
                avect[i][idof] = a;
                bvect[i][idof] = b;
                cvect[i][idof] = c - taumin[idof];
                // lower limits
                avect[i][ndof+idof] = -a;
                bvect[i][ndof+idof] = -b;
                cvect[i][ndof+idof] = -c + taumax[idof];
            }
        }
    }
    nconstraints = int(avect.front().size());
}

} // end namespace TOPP

#endif
