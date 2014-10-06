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


#include "FrictionLimits.h"
#include <math.h>
#include <cmath>

#define CLA_Nothing 0

using namespace OpenRAVE;

namespace TOPP {

FrictionLimits::FrictionLimits(RobotBasePtr probot, std::string& constraintsstring, Trajectory* ptraj){

    trajectory = *ptraj;
    std::string buff;
    std::istringstream iss(constraintsstring);

    getline(iss, buff, '\n');
    discrtimestep = atof(buff.c_str());
    getline(iss, buff, '\n');
    VectorFromString(buff, vmax);
    ndof = probot->GetDOF();
    getline(iss, buff, '\n');
    nbottle = atoi(buff.c_str());

    for (int i = 0; i < nbottle; i++) {
        objspecs.resize(0);
        getline(iss, buff, '\n');
        VectorFromString(buff, objspecs);
        dxvect.push_back(objspecs[0]);
        dyvect.push_back(objspecs[1]);
        bottlehvect.push_back(objspecs[2]);
    }

    getline(iss, buff, '\n');
    mu = atof(buff.c_str());

    hasvelocitylimits = VectorMax(vmax) > TINY;

    std::vector<dReal> amax;
    probot->GetDOFAccelerationLimits(amax);

    //Check Soundness
    BOOST_ASSERT(ndof == ptraj->dimension);

    //Links
    linksvector = probot->GetLinks();
    nlink = linksvector.size();

    int traylinkindex = nlink - (nbottle + 1);

    for (int i = 0; i < nbottle; i++) {
        mbvect.push_back(linksvector[nlink - (i + 1)]->GetMass());
    }

    Vector g = probot->GetEnv()->GetPhysicsEngine()->GetGravity();
    Vector worldx(1, 0, 0), worldy(0, 1, 0), worldz(0, 0, 1);

    std::vector<dReal> q(ndof), qd(ndof), qdd(ndof);

    //vectors storing inequalities coefficients
    std::vector<dReal> a, b, c;

    //intermediate variables
    std::vector<dReal> qplusdeltaqd(ndof);
    boost::multi_array<dReal, 2> jacobianb_p, jacobianb_o, jacobianb_pdelta, jacobianb_odelta,
    jacobianb_pdiff(boost::extents[3][ndof]), jacobianb_odiff(boost::extents[3][ndof]);
    boost::multi_array<dReal, 2> jacobiant_p, jacobiant_o, jacobiant_pdelta, jacobiant_odelta,
    jacobiant_pdiff(boost::extents[3][ndof]), jacobiant_odiff(boost::extents[3][ndof]);

    eps = 1e-10;
    dReal delta = TINY2;
    int ndiscrsteps = int((ptraj->duration + 1e-10)/discrtimestep) + 1;

    {
        EnvironmentMutex::scoped_lock lock(probot->GetEnv()->GetMutex());
        for (int t = 0; t < ndiscrsteps; t++) {
            dReal s = t*discrtimestep;

            ptraj->Eval(s, q);
            ptraj->Evald(s, qd);
            ptraj->Evaldd(s, qdd);

            probot->SetDOFValues(q, CLA_Nothing);
            probot->SetDOFVelocities(qd, CLA_Nothing);

            //intermediate variables for Jacobian derivatives
            dReal norm_qd = VectorNorm(qd);
            bool qdiszero = norm_qd < TINY;

            if(!qdiszero) {
                VectorAdd(q, qd, qplusdeltaqd, 1, delta/norm_qd);
		// VectorAdd(q, qd, qplusdeltaqd, 1, delta);
            }

            RaveTransform<dReal> Htray = probot->GetLinks()[traylinkindex]->GetTransform();
            // Vector Pt = probot->GetLinks()[traylinkindex]->GetGlobalCOM();
            Vector nz_tray = Htray.rotate(worldz); //normal vector of the tray described in world's frame

            // calculate Jacobians of the tray
            // probot->CalculateJacobian(traylinkindex, Pt, jacobiant_p);
            // probot->CalculateAngularVelocityJacobian(traylinkindex, jacobiant_o);

            // if (!qdiszero) {
            //  //******** Set new DOF values and compute derivatives of Jacobians **********
            //     probot->SetDOFValues(qplusdeltaqd, CLA_Nothing);
            //  //*** tray ***
            //     probot->CalculateJacobian(traylinkindex, linksvector[traylinkindex]->GetGlobalCOM(), jacobiant_pdelta);
            //     probot->CalculateAngularVelocityJacobian(traylinkindex, jacobiant_odelta);
            //     MatrixAdd(jacobiant_pdelta, jacobiant_p, jacobiant_pdiff, norm_qd/delta, -norm_qd/delta);
            //     MatrixAdd(jacobiant_odelta, jacobiant_o, jacobiant_odiff, norm_qd/delta, -norm_qd/delta);
            // }

            // Vector temp1 = MatrixMultVector(jacobiant_o, qd);
            // Vector temp2 = MatrixMultVector(jacobiant_o, qdd) + MatrixMultVector(jacobiant_odiff, qd);

            // Vector Pts = MatrixMultVector(jacobiant_p, qd);
            // Vector Ptss = MatrixMultVector(jacobiant_p, qdd) + MatrixMultVector(jacobiant_pdiff, qd);

            dReal r2 = sqrt(2);

            //iterate over n bottles
            for (int i = 0; i < nbottle; i++) {

                //******** Set DOF Values back ********
                probot->SetDOFValues(q, CLA_Nothing);

                a.resize(0);
                b.resize(0);
                c.resize(0);

                int bottlelinkindex = nlink - (i + 1);
                dReal mb = mbvect[i];
                dReal dx = dxvect[i];
                dReal dy = dyvect[i];
                dReal bottleh = bottlehvect[i];

                RaveTransform<dReal> Hbottle = probot->GetLinks()[bottlelinkindex]->GetTransform();
                Vector Pb = probot->GetLinks()[bottlelinkindex]->GetGlobalCOM();
		
                Vector nx_bottle = Hbottle.rotate(worldx);
                Vector ny_bottle = Hbottle.rotate(worldy);
		Vector nz_bottle = Hbottle.rotate(worldz);

                //Calculate Jacobians of the bottle
                probot->CalculateJacobian(bottlelinkindex, Pb, jacobianb_p);
                probot->CalculateAngularVelocityJacobian(bottlelinkindex, jacobianb_o);

                if (!qdiszero) {
                    //******** Set new DOF values and compute derivatives of Jacobians **********
                    probot->SetDOFValues(qplusdeltaqd, CLA_Nothing);

                    //*** bottle ***
                    probot->CalculateJacobian(bottlelinkindex, linksvector[bottlelinkindex]->GetGlobalCOM(), jacobianb_pdelta);
                    probot->CalculateAngularVelocityJacobian(bottlelinkindex, jacobianb_odelta);
                    MatrixAdd(jacobianb_pdelta, jacobianb_p, jacobianb_pdiff, norm_qd/delta, -norm_qd/delta);
                    MatrixAdd(jacobianb_odelta, jacobianb_o, jacobianb_odiff, norm_qd/delta, -norm_qd/delta);
		    // MatrixAdd(jacobianb_pdelta, jacobianb_p, jacobianb_pdiff, 1.0/delta, -1.0/delta);
		    // MatrixAdd(jacobianb_odelta, jacobianb_o, jacobianb_odiff, 1.0/delta, -1.0/delta);
                }

                Vector temp1 = MatrixMultVector(jacobianb_o, qd);
                Vector temp2 = MatrixMultVector(jacobianb_o, qdd) + MatrixMultVector(jacobianb_odiff, qd);

                Vector Pbs = MatrixMultVector(jacobianb_p, qd);
                Vector Pbss = MatrixMultVector(jacobianb_p, qdd) + MatrixMultVector(jacobianb_pdiff, qd);

                dReal Ns = mb*nz_tray.dot3(Pbs);
                dReal Nss = mb*nz_tray.dot3(Pbss);
                dReal N0 = -mb*nz_tray.dot3(g);

                //**************** CONSTRAINT I : N >= 0 *****************
                // actually this constraint will never be saturated.
                a.push_back(-Ns);
                b.push_back(-Nss);
                c.push_back(eps - N0);
                //********************************************************

                Vector fs = mb*Pbs - Ns*nz_tray;
                Vector fss = mb*Pbss - Nss*nz_tray;
                Vector f0 = -mb*g - N0*nz_tray;

                //*************** relaxed conditions of friction cone ****************
                //************************* CONSTAINT II/I ***************************
                a.push_back(nx_bottle.dot3(fs) - (mu/r2)*Ns);
                b.push_back(nx_bottle.dot3(fss) - (mu/r2)*Nss);
                c.push_back(nx_bottle.dot3(f0) - (mu/r2)*N0);

                a.push_back(-nx_bottle.dot3(fs) - (mu/r2)*Ns);
                b.push_back(-nx_bottle.dot3(fss) - (mu/r2)*Nss);
                c.push_back(-nx_bottle.dot3(f0) - (mu/r2)*N0);

                //************************ CONSTRAINT II/II **************************
                a.push_back(ny_bottle.dot3(fs) - (mu/r2)*Ns);
                b.push_back(ny_bottle.dot3(fss) - (mu/r2)*Nss);
                c.push_back(ny_bottle.dot3(f0) - (mu/r2)*N0);

                a.push_back(-ny_bottle.dot3(fs) - (mu/r2)*Ns);
                b.push_back(-ny_bottle.dot3(fss) - (mu/r2)*Nss);
                c.push_back(-ny_bottle.dot3(f0) - (mu/r2)*N0);
                //********************************************************************

                //*************** general conditions of friction force : fz = 0 ****************
                // a.push_back(nz.dot3(fs));
                // b.push_back(nz.dot3(fss));
                // c.push_back(nz.dot3(f0) - eps);

                // a.push_back(-nz.dot3(fs));
                // b.push_back(-nz.dot3(fss));
                // c.push_back(-nz.dot3(f0) - eps);

                //******************************************************************************

                //************************ Calculate the rate of change of the momentum ******************************

                RaveTransformMatrix<dReal> matIb;
                boost::multi_array<dReal, 2> localIb(boost::extents[3][3]), Ib(boost::extents[3][3]);
                boost::multi_array<dReal, 2> Rb(boost::extents[3][3]);

                matIb = linksvector[bottlelinkindex]->GetLocalInertia();
                localIb = ExtractI(matIb); // bottle's inertia tensor in its own frame

                Rb = ExtractR(Hbottle); // rotation matrix of the bottle
                Ib = MatricesMult3(Rb, MatricesMult3(localIb, MatrixTrans(Rb))); // bottle's inertia tensor in the world's frame

                //temp1 and temp2 are defined before (outside this loop)
                Vector Ms = MatrixMultVector(Ib, temp1);
                Vector Mss = MatrixMultVector(Ib, temp2) + temp1.cross(Ms);
                //*************************************************************************************************

                Vector Es = nz_bottle.cross(Ms) + mb*bottleh*Pbs;
                Vector Ess = nz_bottle.cross(Mss) + mb*bottleh*Pbss;
                Vector E0 = -mb*bottleh*g;

                Vector Ks = bottleh*Ns*nz_bottle - Es;
                Vector Kss = bottleh*Nss*nz_bottle - Ess;
                Vector K0 = bottleh*N0*nz_bottle - E0;

                //X
                //************************* CONSTAINT III/I ***************************
                a.push_back(nx_bottle.dot3(Ks) - dx*Ns);
                b.push_back(nx_bottle.dot3(Kss) - dx*Nss);
                c.push_back(nx_bottle.dot3(K0) - dx*N0);

                a.push_back(-nx_bottle.dot3(Ks) - dx*Ns);
                b.push_back(-nx_bottle.dot3(Kss) - dx*Nss);
                c.push_back(-nx_bottle.dot3(K0) - dx*N0);

                //Y
                //************************ CONSTRAINT III/II **************************
                a.push_back(ny_bottle.dot3(Ks) - dy*Ns);
                b.push_back(ny_bottle.dot3(Kss) - dy*Nss);
                c.push_back(ny_bottle.dot3(K0) - dy*N0);

                a.push_back(-ny_bottle.dot3(Ks) - dy*Ns);
                b.push_back(-ny_bottle.dot3(Kss) - dy*Nss);
                c.push_back(-ny_bottle.dot3(K0) - dy*N0);

                //std::cout << Ns << "\t" << Nss << "\t" << N0 << "\t" << "\n";

                //************************ ACCELERATION LIMITS **************************
                for (int j = 0; j < ndof; j++) {
                    a.push_back(qd[j]);
                    b.push_back(qdd[j]);
                    c.push_back(-amax[j]);

                    a.push_back(-qd[j]);
                    b.push_back(-qdd[j]);
                    c.push_back(-amax[j]);
                }

                avect.push_back(a);
                bvect.push_back(b);
                cvect.push_back(c);
            }
        }
    }
    nconstraints = int(avect.front().size());
}


void FrictionLimits::MatrixAdd(const boost::multi_array<dReal, 2>& A, const boost::multi_array<dReal, 2>& B, boost::multi_array<dReal, 2>& C, dReal coefA, dReal coefB) {
    for (int i = 0; i < int(A.shape()[0]); i++) {
        for (int j = 0; j < int(A.shape()[1]); j++) {
            C[i][j] = coefA*A[i][j] + coefB*B[i][j];
        }
    }
}

boost::multi_array<dReal, 2> FrictionLimits::MatrixTrans(const boost::multi_array<dReal, 2>& A) {
    boost::multi_array<dReal, 2> At(boost::extents[int(A.shape()[1])][int(A.shape()[0])]);
    for (int i = 0; i < int(A.shape()[1]); i++) {
        for (int j = 0; j < int(A.shape()[0]); j++) {
            At[i][j] = A[j][i];
        }
    }
    return At;
}

boost::multi_array<dReal, 2> FrictionLimits::MatricesMult3(const boost::multi_array<dReal, 2>& A, const boost::multi_array<dReal, 2>& B) {
    BOOST_ASSERT(int(A.shape()[0] == 3));
    BOOST_ASSERT(int(A.shape()[1] == 3));
    BOOST_ASSERT(int(B.shape()[0] == 3));
    BOOST_ASSERT(int(B.shape()[1] == 3));
    boost::multi_array<dReal, 2> C(boost::extents[3][3]);
    for (int i = 0; i < 3; i++) {
        for (int j = 0; j < 3; j++) {
            C[i][j] = A[i][0]*B[0][j] + A[i][1]*B[1][j] + A[i][2]*B[2][j];
        }
    }
    return C;
}




Vector FrictionLimits::MatrixMultVector(const boost::multi_array<dReal, 2>& M, const std::vector<dReal>& v) {
    BOOST_ASSERT(M.shape()[1] == v.size());
    Vector res;

    for (int i = 0; i < int(M.shape()[0]); i++) {
        res[i] = 0;
        for (int j = 0; j < int(v.size()); j++) {
            res[i] += M[i][j]*v[j];
        }
    }
    return res;
}

Vector FrictionLimits::MatrixMultVector(const boost::multi_array<dReal, 2>& M, const Vector& v) {
    Vector res;

    for (int i = 0; i < int(M.shape()[0]); i++) {
        res[i] = 0;
        for (int j = 0; j < 3; j++) {
            res[i] += M[i][j]*v[j];
        }
    }
    return res;
}

boost::multi_array<dReal, 2> FrictionLimits::ExtractI(const RaveTransformMatrix<dReal>& T) { //H) {
    boost::multi_array<dReal, 2> I(boost::extents[3][3]);
    // Vector x(1, 0, 0), y(0, 1, 0), z(0, 0, 1);
    // Vector resx, resy, resz;

    I[0][0] = T.m[0]; I[0][1] = T.m[1]; I[0][2] = T.m[2];
    I[1][0] = T.m[4]; I[1][1] = T.m[5]; I[1][2] = T.m[6];
    I[2][0] = T.m[8]; I[2][1] = T.m[9]; I[2][2] = T.m[10];
        
    // resx = H.rotate(x);
    // resy = H.rotate(y);
    // resz = H.rotate(z);

    // I[0][0] = resx[0]; I[0][1] = resy[0]; I[0][2] = resz[0];
    // I[1][0] = resx[1]; I[1][1] = resy[1]; I[1][2] = resz[1];
    // I[2][0] = resx[2]; I[2][1] = resy[2]; I[2][2] = resz[2];
    return I;
}

boost::multi_array<dReal, 2> FrictionLimits::ExtractR(const RaveTransform<dReal>& H) {
    boost::multi_array<dReal, 2> R(boost::extents[3][3]);
    Vector x(1, 0, 0), y(0, 1, 0), z(0, 0, 1);
    Vector resx, resy, resz;

    resx = H.rotate(x);
    resy = H.rotate(y);
    resz = H.rotate(z);

    R[0][0] = resx[0]; R[0][1] = resy[0]; R[0][2] = resz[0];
    R[1][0] = resx[1]; R[1][1] = resy[1]; R[1][2] = resz[1];
    R[2][0] = resx[2]; R[2][1] = resy[2]; R[2][2] = resz[2];
    return R;
}

Vector FrictionLimits::ExtractT(const RaveTransform<dReal>& H) {
    Vector T;
    T = H.trans;
    return T;
}

}
