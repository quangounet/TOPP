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

#define CLA_Nothing 0

using namespace OpenRAVE;

namespace TOPP {

ZMPTorqueLimits::ZMPTorqueLimits(const std::string& constraintsstring, Trajectory* ptraj, const Tunings& tunings, RobotBasePtr probot0){
    int buffsize = BUFFSIZE;  // TODO: remove this dirty string interface!
    std::vector<dReal> tmpvect;
    char buff[buffsize];
    std::istringstream iss(constraintsstring);
    iss.getline(buff,buffsize);
    VectorFromString(std::string(buff),taumin);
    iss.getline(buff,buffsize);
    VectorFromString(std::string(buff),taumax);
    iss.getline(buff,buffsize);
    VectorFromString(std::string(buff),zmplimits);
    iss.getline(buff,buffsize);
    VectorFromString(std::string(buff),vmax);
    hasvelocitylimits = VectorMax(vmax) > TINY;
    maxrep = 1;

    probot = probot0;

    dReal xmin = zmplimits[0];
    dReal xmax = zmplimits[1];
    dReal ymin = zmplimits[2];
    dReal ymax = zmplimits[3];

    // General
    linksvector = probot->GetLinks();
    nlinks = int(linksvector.size());

    for(int i=0; i < int(nlinks); i++) {
        mass.push_back(linksvector[i]->GetMass());
        totalmass += mass[i];
    }
    int ndof = ptraj->dimension;
    assert(ndof == probot->GetDOF());
    Vector g = probot->GetEnv()->GetPhysicsEngine()->GetGravity();
    std::vector<dReal> q(ndof), qd(ndof), qdd(ndof);

    // Torque intermediate variables
    std::vector<dReal> a,b,c, torquesimple;
    boost::array<std::vector<dReal>,3> torquecomponents;

    // ZMP intermediate variables
    std::vector<std::pair<Vector,Vector> > linkaccelerations;
    std::vector<dReal> qplusdeltaqd(ndof);
    boost::multi_array< dReal, 2 > jacobian, jacobiandelta, jacobiandiff(boost::extents[3][ndof]);

    // Initialize jacobian diff
    probot->CalculateJacobian(0,linksvector[0]->GetGlobalCOM(),jacobian);

    dReal delta = TINY2;
    Vector ci, ciVg, ci_qXqdd, ciVci_qXqdd, qdXci_qqXqd, ciVqdXci_qqXqd, Atau, Btau, Ctau, Ah, Bh, Ch;
    dReal a_xmax=0, b_xmax=0, c_xmax=0;
    dReal a_xmin=0, b_xmin=0, c_xmin=0;
    dReal a_ymax=0, b_ymax=0, c_ymax=0;
    dReal a_ymin=0, b_ymin=0, c_ymin=0;

    std::cout << "\n\n\n\n";

    int ndiscrsteps = int((ptraj->duration+1e-10)/tunings.discrtimestep)+1;

    for(int t = 0; t<ndiscrsteps; t++) {
        dReal s = t*tunings.discrtimestep;
        ptraj->Eval(s,q);
        ptraj->Evald(s,qd);
        ptraj->Evaldd(s,qdd);
        //std::cout << COM(q) << "\n";
        ZMP(q,qd,qdd);
    }

    std::cout << "\n\n\n\n";

    {
        EnvironmentMutex::scoped_lock lock(probot->GetEnv()->GetMutex());
        for(int t = 0; t<ndiscrsteps; t++) {
            dReal s = t*tunings.discrtimestep;
            ptraj->Eval(s,q);
            ptraj->Evald(s,qd);
            ptraj->Evaldd(s,qdd);
            probot->SetDOFValues(q,CLA_Nothing);
            probot->SetDOFVelocities(qd,CLA_Nothing);
            a.resize(0);
            b.resize(0);
            c.resize(0);

            // Torque limits
            probot->ComputeInverseDynamics(torquesimple,qd);
            probot->ComputeInverseDynamics(torquecomponents,qdd);
            for(int j=0; j<ndof; j++) {
                // Add inequalities only when taumax != 0
                if(std::abs(taumax[j])>TINY2) {
                    a.push_back(torquesimple[j] - torquecomponents[1][j] - torquecomponents[2][j]);
                    a.push_back(-a.back());
                    b.push_back(torquecomponents[0][j] + torquecomponents[1][j]);
                    b.push_back(-b.back());
                    c.push_back(torquecomponents[2][j] - taumax[j]);
                    c.push_back(-torquecomponents[2][j] + taumin[j]);
                }
            }

            // ZMP limits (only processed when xmax>xmin)
            if(xmax-xmin>=TINY2) {
                dReal norm_qd = VectorNorm(qd);
                VectorAdd(q,qd,qplusdeltaqd,1,delta/norm_qd);

                dReal tau0 = 0, tau1 = 0, h2 = 0;

                for(int i=0; i < int(nlinks); i++) {
                    // Set DOFValues to q and extract jacobian
                    probot->SetDOFValues(q,CLA_Nothing);
                    ci = linksvector[i]->GetGlobalCOM();
                    probot->CalculateJacobian(i,ci,jacobian);

                    // Set DOFValues to qplusdeltaqd and extract jacobian
                    probot->SetDOFValues(qplusdeltaqd,CLA_Nothing);
                    probot->CalculateJacobian(i,linksvector[i]->GetGlobalCOM(),jacobiandelta);
                    //Calculate the derivative of the jacobian
                    MatrixAdd(jacobiandelta,jacobian,jacobiandiff,norm_qd/delta,-norm_qd/delta);
                    // Compute the components
                    ci_qXqdd = MatrixMultVector(jacobian,qdd);
                    ciVci_qXqdd = ci.cross(ci_qXqdd);
                    qdXci_qqXqd = MatrixMultVector(jacobiandiff,qd);
                    ciVqdXci_qqXqd = ci.cross(qdXci_qqXqd);
                    ciVg = ci.cross(g);

                    // tau = mass[i]*(Atau*sdd + Btau*sd^2 + Ctau)
                    Atau = -ciVci_qXqdd;
                    Btau = -ciVci_qXqdd - ciVqdXci_qqXqd;
                    Ctau = ciVg;
                    // h = Ah*sdd + Bh*sd^2 + Ch
                    Ah = -ci_qXqdd;
                    Bh = -ci_qXqdd - qdXci_qqXqd;
                    Ch = g;


                    // // Testing
                    // dReal sdotsq = 1;
                    // tau0 += mass[i]*(-sdotsq*ciVci_qXqdd[0]-sdotsq*ciVqdXci_qqXqd[0]+ciVg[0]);
                    // tau1 += mass[i]*(-sdotsq*ciVci_qXqdd[1]-sdotsq*ciVqdXci_qqXqd[1]+ciVg[1]);
                    // h2 += mass[i]*(-sdotsq*ci_qXqdd[2]-sdotsq*qdXci_qqXqd[2]+g[2]);

                    // x_zmp <= x_max
                    a_xmax += mass[i] * (Atau[1]+xmax*Ah[2]);
                    b_xmax += mass[i] * (Btau[1]+xmax*Bh[2]);
                    c_xmax += mass[i] * (Ctau[1]+xmax*Ch[2]);

                    // x_zmp >= x_min
                    a_xmin += mass[i] * (-Atau[1]-xmin*Ah[2]);
                    b_xmin += mass[i] * (-Btau[1]-xmin*Bh[2]);
                    c_xmin += mass[i] * (-Ctau[1]-xmin*Ch[2]);

                    // y_zmp <= y_max
                    a_ymax += mass[i] * (-Atau[0]+ymax*Ah[2]);
                    b_ymax += mass[i] * (-Btau[0]+ymax*Bh[2]);
                    c_ymax += mass[i] * (-Ctau[0]+ymax*Ch[2]);

                    // y_zmp <= y_min
                    a_ymin += mass[i] * (Atau[0]-ymin*Ah[2]);
                    b_ymin += mass[i] * (Btau[0]-ymin*Bh[2]);
                    c_ymin += mass[i] * (Ctau[0]-ymin*Ch[2]);
                }

                a.push_back(a_xmax);
                a.push_back(a_xmin);
                a.push_back(a_ymax);
                a.push_back(a_ymin);
                b.push_back(b_xmax);
                b.push_back(b_xmin);
                b.push_back(b_ymax);
                b.push_back(b_ymin);
                c.push_back(c_xmax);
                c.push_back(c_xmin);
                c.push_back(c_ymax);
                c.push_back(c_ymin);

            }
            avect.push_back(a);
            bvect.push_back(b);
            cvect.push_back(c);
        }
    }

    nconstraints = int(avect.front().size());

}


Vector ZMPTorqueLimits::COM(std::vector<dReal>& q){
    Vector com;
    {
        EnvironmentMutex::scoped_lock lock(probot->GetEnv()->GetMutex());
        probot->SetDOFValues(q,CLA_Nothing);
        for(int i=0; i < int(nlinks); i++) {
            com += linksvector[i]->GetMass() * linksvector[i]->GetGlobalCOM();
        }
    }
    return 1/totalmass * com;
}


Vector ZMPTorqueLimits::ZMP(std::vector<dReal>& q, std::vector<dReal>& qd, std::vector<dReal>& qdd, bool withangularmomentum){
    Vector tau0;
    Vector g = probot->GetEnv()->GetPhysicsEngine()->GetGravity();
    dReal f02 = totalmass * g[2];
    std::vector<std::pair<Vector,Vector> > linkvelocities, linkaccelerations;
    std::vector<dReal> zeros(4);
    Vector ci, cid, cidd, localcom, linvel, angvel, linacc, angacc, ri;
    Transform T;
    {
        EnvironmentMutex::scoped_lock lock(probot->GetEnv()->GetMutex());
        probot->SetDOFValues(q,CLA_Nothing);
        probot->SetDOFVelocities(qd,CLA_Nothing);
        probot->GetLinkVelocities(linkvelocities);
        probot->GetLinkAccelerations(qdd,linkaccelerations);
        for(int i=0; i < int(nlinks); i++) {
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
    std::cout << 1/f02 * tau0 << "\n";
    return 1/f02 * tau0;
}


Vector MatrixMultVector(const boost::multi_array<dReal,2>& M, const std::vector<dReal>& v){
    Vector res;
    assert(M.shape()[0] == 3);
    assert(M.shape()[1] == v.size());
    for(int i=0; i<3; i++) {
        res[i] = 0;
        for(int j=0; j<int(v.size()); j++) {
            res[i] += M[i][j]*v[j];
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


}
