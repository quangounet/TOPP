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


#include "ZMPTorqueLimits.h"

#include <boost/numeric/ublas/lu.hpp>

#define CLA_Nothing 0 // no warning when setting DOF values or velocities

// \brief Different in z-coordinate used to detect the support foot in the
// single-support case.
const dReal SINGLE_SUPPORT_Z_DIFF = 1e-2; 

using namespace OpenRAVE;

namespace TOPP {

ZMPTorqueLimits::ZMPTorqueLimits(const std::string& constraintsstring, Trajectory* ptraj, const Tunings& tunings, RobotBasePtr probot0){
    int buffsize = BUFFSIZE;  // TODO: remove this dirty string interface!
    std::vector<dReal> tmpvect, activedofs, activelinks0;
    char buff[buffsize];
    std::istringstream iss(constraintsstring);
    iss.getline(buff,buffsize);
    VectorFromString(std::string(buff),activedofs);
    iss.getline(buff,buffsize);
    VectorFromString(std::string(buff),activelinks0);
    iss.getline(buff,buffsize);
    VectorFromString(std::string(buff),taumin);
    iss.getline(buff,buffsize);
    VectorFromString(std::string(buff),taumax);
    iss.getline(buff,buffsize);
    VectorFromString(std::string(buff),zmplimits);
    iss.getline(buff,buffsize);
    VectorFromString(std::string(buff),vmax);
    hasvelocitylimits = VectorMax(vmax) > TINY;
    maxrep = 5;

    probot = probot0;
    activelinks = activelinks0;

    //Check soundness
    assert(int(activedofs.size()) == probot->GetDOF());
    assert(activelinks.size() == probot->GetLinks().size());
    assert(zmplimits.size() == 4);

    // DOFs
    for(int i=0; i<int(activedofs.size()); i++) {
        if(activedofs[i]>TINY) {
            dofsvector.push_back(i);
        }
    }
    ndof = int(dofsvector.size());
    assert(ndof == ptraj->dimension);

    // Links
    nlink0 = int(activelinks.size());
    linksvector = probot->GetLinks();
    for(int i=0; i < nlink0; i++) {
        mass.push_back(linksvector[i]->GetMass());
        totalmass += mass[i];
    }

    // ZMP limits
    dReal xmin = zmplimits[0];
    dReal xmax = zmplimits[1];
    dReal ymin = zmplimits[2];
    dReal ymax = zmplimits[3];

    cout << "ZMP limits: ";
    cout << xmin << " " << xmax << " " <<  ymin << " " << ymax << "\n";

    Vector g = probot->GetEnv()->GetPhysicsEngine()->GetGravity();
    std::vector<dReal> q(ndof), qd(ndof), qdd(ndof), qfilled, qdfilled, qddfilled;
    probot->GetDOFValues(qfilled);
    probot->GetDOFVelocities(qdfilled);
    qddfilled.resize(probot->GetDOF());

    // Torque intermediate variables
    std::vector<dReal> a,b,c, torquesimple;
    boost::array<std::vector<dReal>,3> torquecomponents;
    std::vector<dReal> torquesimpletrimmed(ndof), torquecomponents0trimmed(ndof), torquecomponents1trimmed(ndof), torquecomponents2trimmed(ndof);

    // ZMP intermediate variables
    std::vector<std::pair<Vector,Vector> > linkaccelerations;
    std::vector<dReal> qplusdeltaqdfilled(qfilled.size());
    boost::multi_array< dReal, 2 > jacobian, jacobiandelta, jacobiandiff(boost::extents[3][probot->GetDOF()]);

    // Initialize jacobian diff
    probot->CalculateJacobian(0,probot->GetLinks()[0]->GetGlobalCOM(),jacobian);

    dReal delta = TINY2;
    Vector ci, ciVg, q1, ciVq1, q2, ciVq2, q3, ciVq3;

    int ndiscrsteps = int((ptraj->duration+1e-15)/tunings.discrtimestep)+1;
    dReal dt = ptraj->duration/(ndiscrsteps-1);


    // for(int t = 0; t<ndiscrsteps; t++) {
    //     dReal s = t*tunings.discrtimestep;
    //     ptraj->Eval(s,q);
    //     ptraj->Evald(s,qd);
    //     ptraj->Evaldd(s,qdd);
    //     Fill(q,qfilled);
    //     Fill(qd,qdfilled);
    //     Fill(qdd,qddfilled);
    //     std::cout << ZMP(qfilled,qdfilled,qddfilled) << "\n";
    // }

    {
        EnvironmentMutex::scoped_lock lock(probot->GetEnv()->GetMutex());
        for(int t = 0; t<ndiscrsteps; t++) {
            dReal s = t*dt;
            ptraj->Eval(s,q);
            ptraj->Evald(s,qd);
            ptraj->Evaldd(s,qdd);
            Fill(q,qfilled);
            Fill(qd,qdfilled);
            Fill(qdd,qddfilled);
            probot->SetDOFValues(qfilled,CLA_Nothing);
            probot->SetDOFVelocities(qdfilled,CLA_Nothing);
            a.resize(0);
            b.resize(0);
            c.resize(0);

            // Torque limits
            //probot->ComputeInverseDynamics(torquesimple,qdfilled);
            //probot->ComputeInverseDynamics(torquecomponents,qddfilled);
            this->ComputeInverseDynamicsSingleSupport(torquesimple, qdfilled);
            this->ComputeInverseDynamicsSingleSupport(torquecomponents, qddfilled);
            Trim(torquesimple,torquesimpletrimmed);
            Trim(torquecomponents[0],torquecomponents0trimmed);
            Trim(torquecomponents[1],torquecomponents1trimmed);
            Trim(torquecomponents[2],torquecomponents2trimmed);
            for(int j=0; j<ndof; j++) {
                // Add inequalities only when taumax != 0
                if(std::abs(taumax[j])>TINY2) {
                    a.push_back(torquesimpletrimmed[j] - torquecomponents1trimmed[j] - torquecomponents2trimmed[j]);
                    a.push_back(-a.back());
                    b.push_back(torquecomponents0trimmed[j] + torquecomponents1trimmed[j]);
                    b.push_back(-b.back());
                    c.push_back(torquecomponents2trimmed[j] - taumax[j]);
                    c.push_back(-torquecomponents2trimmed[j] + taumin[j]);
                }
            }

            // ZMP limits (only processed when xmax>xmin)
            if(xmax-xmin>=TINY2) {
                dReal norm_qd = VectorNorm(qd);
                bool qdiszero = norm_qd<TINY;
                if(!qdiszero) {
                    VectorAdd(qfilled,qdfilled,qplusdeltaqdfilled,1,delta/norm_qd);
                }
                Vector tau,h;
                Vector Atau, Btau, Ctau, Cisum, Ah, Bh, Ch;
                for(int i=0; i < int(nlink0); i++) {
                    if(activelinks[i]<=TINY) {
                        continue;
                    }
                    // Set DOFValues to q and extract jacobian
                    probot->SetDOFValues(qfilled,CLA_Nothing);
                    ci = linksvector[i]->GetGlobalCOM();
                    probot->CalculateJacobian(i,ci,jacobian);

                    // Set DOFValues to qplusdeltaqdfilled and extract jacobian
                    if(!qdiszero) {
                        probot->SetDOFValues(qplusdeltaqdfilled,CLA_Nothing);
                        probot->CalculateJacobian(i,linksvector[i]->GetGlobalCOM(),jacobiandelta);
                        //Calculate the derivative of the jacobian
                        MatrixAdd(jacobiandelta,jacobian,jacobiandiff,norm_qd/delta,-norm_qd/delta);
                    }

                    // Compute the components
                    q1 = MatrixMultVector(jacobian,qd);
                    q2 = MatrixMultVector(jacobian,qdd);
                    if(!qdiszero) {
                        q3 = MatrixMultVector(jacobiandiff,qd);
                    }
                    else{
                        q3 = 0.*q1;
                    }
                    ciVq1 = ci.cross(q1);
                    ciVq2 = ci.cross(q2);
                    ciVq3 = ci.cross(q3);
                    // h =  mass[i]*(Ah*sdd + Bh*sd^2 + Ch)
                    Ah += mass[i]*(-q1);
                    Bh += mass[i]*(-q2-q3);
                    // tau = mass[i]*(Atau*sdd + Btau*sd^2 + Ctau)
                    Atau += mass[i]*(-ciVq1);
                    Btau += mass[i]*(-ciVq2-ciVq3);
                    Cisum += mass[i]*ci;
                }
                dReal Ch2 = totalmass*g[2];
                Ctau = Cisum.cross(g);
                //std::cout << -(Btau[1]+Ctau[1])/(Bh[2]+Ch2) << " " << (Btau[0]+Ctau[0])/(Bh[2]+Ch2) <<"\n";
                // x_zmp <= x_max
                dReal a_xmax = (Atau[1]+xmax*Ah[2]);
                dReal b_xmax = (Btau[1]+xmax*Bh[2]);
                dReal c_xmax = (Ctau[1]+xmax*Ch2);
                // x_zmp >= x_min
                dReal a_xmin = (-Atau[1]-xmin*Ah[2]);
                dReal b_xmin = (-Btau[1]-xmin*Bh[2]);
                dReal c_xmin = (-Ctau[1]-xmin*Ch2);
                // y_zmp <= y_max
                dReal a_ymax = (-Atau[0]+ymax*Ah[2]);
                dReal b_ymax = (-Btau[0]+ymax*Bh[2]);
                dReal c_ymax = (-Ctau[0]+ymax*Ch2);
                // y_zmp >= y_min
                dReal a_ymin = (Atau[0]-ymin*Ah[2]);
                dReal b_ymin = (Btau[0]-ymin*Bh[2]);
                dReal c_ymin = (Ctau[0]-ymin*Ch2);

                a.push_back(a_xmax);
                b.push_back(b_xmax);
                c.push_back(c_xmax);

                a.push_back(a_xmin);
                b.push_back(b_xmin);
                c.push_back(c_xmin);

                a.push_back(a_ymax);
                b.push_back(b_ymax);
                c.push_back(c_ymax);

                a.push_back(a_ymin);
                b.push_back(b_ymin);
                c.push_back(c_ymin);
            }
            avect.push_back(a);
            bvect.push_back(b);
            cvect.push_back(c);

        }
    }

    nconstraints = int(avect.front().size());

}


Vector ZMPTorqueLimits::COM(std::vector<dReal>& qfilled){
    Vector com;
    {
        EnvironmentMutex::scoped_lock lock(probot->GetEnv()->GetMutex());
        probot->SetDOFValues(qfilled,CLA_Nothing);
        for(int i=0; i < int(nlink0); i++) {
            if(activelinks[i]<=TINY) {
                continue;
            }
            com += linksvector[i]->GetMass() * linksvector[i]->GetGlobalCOM();
        }
    }
    return 1/totalmass * com;
}


Vector ZMPTorqueLimits::ZMP(std::vector<dReal>& qfilled, std::vector<dReal>& qdfilled, std::vector<dReal>& qddfilled, bool withangularmomentum){
    Vector tau0;
    Vector g = probot->GetEnv()->GetPhysicsEngine()->GetGravity();
    dReal f02 = totalmass * g[2];
    std::vector<std::pair<Vector,Vector> > linkvelocities, linkaccelerations;
    std::vector<dReal> zeros(4);
    Vector ci, cid, cidd, localcom, linvel, angvel, linacc, angacc, ri;
    Transform T;
    {
        EnvironmentMutex::scoped_lock lock(probot->GetEnv()->GetMutex());
        probot->SetDOFValues(qfilled,CLA_Nothing);
        probot->SetDOFVelocities(qdfilled,CLA_Nothing);
        probot->GetLinkVelocities(linkvelocities);
        probot->GetLinkAccelerations(qddfilled,linkaccelerations);
        for(int i=0; i < int(nlink0); i++) {
            if(activelinks[i]<=TINY) {
                continue;
            }
            ri = linksvector[i]->GetTransform().rotate(linksvector[i]->GetLocalCOM());
            ci = linksvector[i]->GetGlobalCOM();
            linvel = linkvelocities[i].first;
            angvel = linkvelocities[i].second;
            linacc = linkaccelerations[i].first;
            angacc = linkaccelerations[i].second;
            cid = linvel + angvel.cross(ri);
            cidd = linacc + angvel.cross(angvel.cross(ri))+angacc.cross(ri);
            tau0 +=  mass[i]*ci.cross(g-cidd);
            f02 -= mass[i]*cidd[2];
        }
    }
    dReal temp = tau0[0];
    tau0[0] = -tau0[1];
    tau0[1] = temp;
    tau0[2] = 0;
    //std::cout << 1/f02 * tau0 << "\n";
    return 1/f02 * tau0;
}


void ZMPTorqueLimits::Fill(const std::vector<dReal>&q, std::vector<dReal>&qfilled){
    for(int i=0; i<ndof; i++) {
        qfilled[dofsvector[i]]=q[i];
    }
}


void ZMPTorqueLimits::Trim(const std::vector<dReal>&q, std::vector<dReal>&qtrimmed){
    for(int i=0; i<ndof; i++) {
        qtrimmed[i]=q[dofsvector[i]];
    }
}


Vector ZMPTorqueLimits::MatrixMultVector(const boost::multi_array<dReal,2>& M, const std::vector<dReal>& v){
    Vector res;
    assert(M.shape()[0] == 3);
    for(int i=0; i<3; i++) {
        res[i] = 0;
        for(int j=0; j<int(v.size()); j++) {
            res[i] += M[i][dofsvector[j]]*v[j];
        }
    }
    return res;
}


void MatrixAdd(const boost::multi_array<dReal,2>& A, const boost::multi_array<dReal,2>& B, boost::multi_array<dReal,2>& C, dReal coefA, dReal coefB){
    for(int i=0; i<int(A.shape()[0]); i++) {
        for(int j=0; j<int(A.shape()[1]); j++) {
            C[i][j] = coefA*A[i][j] + coefB*B[i][j];
        }
    }
}


void ZMPTorqueLimits::GetFootJacobianTranspose(const KinBody::LinkPtr foot, 
        boost::multi_array<dReal,2>& mjacobiantrans) {
    int linkindex = foot->GetIndex();
    Transform t = foot->GetTransform();

    std::vector<dReal> transjacobian;
    std::vector<dReal> rotjacobian;
    probot->CalculateJacobian(linkindex, t.trans, transjacobian);
    probot->CalculateRotationJacobian(linkindex, t.rot, rotjacobian);

    int nbdof = probot->GetDOF();
    mjacobiantrans.resize(boost::extents[nbdof][7]);
    for (int dof = 0; dof < nbdof; dof++) {
        mjacobiantrans[dof][0] = rotjacobian[0 * nbdof + dof];
        mjacobiantrans[dof][1] = rotjacobian[1 * nbdof + dof];
        mjacobiantrans[dof][2] = rotjacobian[2 * nbdof + dof];
        mjacobiantrans[dof][3] = rotjacobian[3 * nbdof + dof];
        mjacobiantrans[dof][4] = transjacobian[0 * nbdof + dof];
        mjacobiantrans[dof][5] = transjacobian[1 * nbdof + dof];
        mjacobiantrans[dof][6] = transjacobian[2 * nbdof + dof];
    }
}


void ZMPTorqueLimits::ComputeInverseDynamicsSingleSupport(
        boost::array<std::vector<dReal>, 3>& doftorquecomponents, 
        const std::vector<dReal>& dofaccelerations) {
    int nbdof = probot->GetDOF();
    int nbactivedof = nbdof - 6;
    int baselinkstart = nbactivedof;

    // Detect the support foot
    KinBody::LinkPtr leftfoot = probot->GetLink("L_FOOT_LINK");
    KinBody::LinkPtr rightfoot = probot->GetLink("R_FOOT_LINK");
    Transform lefttrans = leftfoot->GetTransform();
    Transform righttrans = rightfoot->GetTransform();
    dReal zdiff = lefttrans.trans.z - righttrans.trans.z;
    assert(zdiff < -SINGLE_SUPPORT_Z_DIFF || zdiff > +SINGLE_SUPPORT_Z_DIFF);
    KinBody::LinkPtr foot = (zdiff < 0) ? leftfoot : rightfoot;

    // Compute the flying-base inverse dynamics
    boost::multi_array<dReal,2> mjacobiantrans;
    std::vector<dReal> taubaselink(6);
    
    this->GetFootJacobianTranspose(foot, mjacobiantrans);
    probot->ComputeInverseDynamics(doftorquecomponents, dofaccelerations);
    for (int i = 0; i < 6; i++) {
        int j = baselinkstart + i;
        taubaselink[i] = doftorquecomponents[0][j] + doftorquecomponents[1][j] + doftorquecomponents[2][j];
    }

    // Compute the contact wrench from the (virtual) base link torques
    boost::multi_array<dReal,2> mactuatedjacobiantrans(boost::extents[nbactivedof][7]);
    boost::numeric::ublas::matrix<dReal> U(6, 7);
    for (int j = 0; j < 7; j++) {
        for (int i = 0; i < nbactivedof; i++)
            mactuatedjacobiantrans[i][j] = mjacobiantrans[i][j];
        for (int i = 0; i < 6; i++)
            U(i, j) = mjacobiantrans[baselinkstart + i][j];
    }

    boost::numeric::ublas::permutation_matrix<size_t> P(6); // U has 6 rows
    boost::numeric::ublas::vector<dReal> x(6);
    for (int i = 0; i < 6; i++)
        x(i) = taubaselink[i];

    boost::numeric::ublas::lu_factorize(U, P);
    boost::numeric::ublas::lu_substitute(U, P, x);

    std::vector<dReal> externalwrench(7);
    for (int i = 0; i < 7; i++)
        externalwrench[i] = x[i];

    // Feed the contact wrench back into the free-flying torque components
    Vector backtorques = MatrixMultVector(mactuatedjacobiantrans, externalwrench);
    for (int i = 0; i < nbactivedof; i++)
        doftorquecomponents[2][i] -= backtorques[i];
    for (int i = nbactivedof; i < nbdof; i++)
        doftorquecomponents[2][baselinkstart + i] -= taubaselink[i];
}


void ZMPTorqueLimits::ComputeInverseDynamicsSingleSupport(
        std::vector<dReal>& doftorques, 
        const std::vector<dReal>& dofaccelerations) {
    boost::array<std::vector<dReal>, 3> doftorquecomponents;
    this->ComputeInverseDynamicsSingleSupport(doftorquecomponents, dofaccelerations);
    int nbdof = probot->GetDOF();
    doftorques.resize(nbdof);
    for (int i = 0; i < nbdof; i++)
        doftorques[i] = doftorquecomponents[0][i] + doftorquecomponents[1][i] + doftorquecomponents[2][i];
}


}
