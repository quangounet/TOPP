#include "SO3Constraints.h"
#include <sstream>
#include <boost/python.hpp>
#include <ctime>

namespace TOPP {
void ComputeSO3Constraints(const std::string& SO3trajstring, const std::string& constraintsstring, boost::python::list& resstringlist){
    //Convert the data
    Trajectory* SO3traj = new Trajectory(SO3trajstring);
    dReal discrtimestep;
    std::vector<dReal> taumax;
    boost::multi_array<dReal, 2> inertia(boost::extents[3][3]);
    bool noInertia = false;
    std::string buff;
    std::istringstream iss(constraintsstring);
    getline(iss, buff, '\n');
    discrtimestep = atof(buff.c_str());
    getline(iss, buff, '\n');
    VectorFromString(buff, taumax);
    for (int i = 0; i < 3; i++) {
        std::string strtmp;
        std::vector<dReal> vecttmp;
        getline(iss, strtmp, '\n');
        if (strtmp == "") {
            noInertia = true;
            std::cout << "Inertia was not input, assume Inertia to be an Identity(3) matrix\n";
            break;
        }
        VectorFromString(strtmp, vecttmp);
        inertia[i][0] = vecttmp[0];
        inertia[i][1] = vecttmp[1];
        inertia[i][2] = vecttmp[2];
    }

    // Process the data
    int ndiscrsteps = int((SO3traj->duration + 1e-10)/discrtimestep) + 1;
    int ndof = 3;
    std::string a, b, c;
    std::vector<dReal> r(ndof), rd(ndof), rdd(ndof);
    dReal nr,nr2,nr3,nr4,nr5;
    for (int i = 0; i< ndiscrsteps; i++) {
        dReal t = i*discrtimestep;
        SO3traj->Eval(t,r);
        SO3traj->Evald(t,rd);
        SO3traj->Evaldd(t,rdd);
        nr = VectorNorm(r)+1e-10;
        nr2 = nr*nr;
        nr3 = nr2*nr;
        nr4 = nr3*nr;
        nr5 = nr4*nr;
        boost::multi_array<dReal, 2> R(boost::extents[3][3]);
        R = SkewFromVect(r);
        dReal snr = std::sin(nr);
        dReal cnr = std::cos(nr);
        std::vector<dReal> rcrd(3);

        rcrd = Cross(r, rd);
        dReal rdrd = DotVect(r,rd);

        boost::multi_array<dReal, 2> Amat(boost::extents[3][3]);
        Amat = MatrixAdd(MatrixAdd(Eye(),R,1,-(1-cnr)/nr2),MatricesMult3(R,R),1,(nr-snr)/nr3);
        std::vector<dReal> C12(3);
        C12 = AddVect(Cross(rd,rcrd),rcrd,(nr-snr)/nr3,-(2*cnr+nr*snr-2)/nr4);
        std::vector<dReal> C(3);
        C = AddVect(C12,Cross(r,rcrd),1,(3*snr-nr*cnr-2*nr)*rdrd/nr5);
        std::vector<dReal> Ard = MatrixMultVector(Amat,rd);

        std::vector<dReal> at;
        std::vector<dReal> bt;
        if (noInertia) {
            //Inertia was not input, assume Inertia to be an Identity matrix. In this case, angular accelerations are the same as torques
            at = Ard;
            bt = AddVect(MatrixMultVector(Amat,rdd),C,1,1);
        }
        else
        { //Inertia was input
            at = MatrixMultVector(inertia,Ard);
            std::vector<dReal> bt_firstpart;
            std::vector<dReal> bt_secondpart;
            bt_firstpart = AddVect(MatrixMultVector(inertia,MatrixMultVector(Amat,rdd)), MatrixMultVector(inertia,C),1,1);
            bt_secondpart = Cross(Ard,at);
            bt = AddVect(bt_firstpart,bt_secondpart,1,1);
        }
        if (i != 0) {
            a.append("\n");
            b.append("\n");
            c.append("\n");
        }
        a.append(VectToString(at,false));
        b.append(VectToString(bt,false));
        c.append(VectToString(taumax,true));
    }
    resstringlist.append(a);
    resstringlist.append(b);
    resstringlist.append(c);
};

std::string VectToString(const std::vector<dReal>&vect, const bool&IsCvect){
    std::stringstream ss;
    if (IsCvect == true) { // Convert c vect to a string
        ss << -vect[0] << " ";
        ss << -vect[1] << " ";
        ss << -vect[2] << " ";
    }
    else{ //Convert a or b vect to a string
        ss << vect[0] << " ";
        ss << vect[1] << " ";
        ss << vect[2] << " ";
    }
    ss << -vect[0] << " ";
    ss << -vect[1] << " ";
    ss << -vect[2];
    return ss.str();
};


boost::multi_array<dReal, 2> SkewFromVect (const std::vector<dReal>&vect){
    boost::multi_array<dReal, 2> Skew(boost::extents[3][3]);
    Skew[0][0] = 0;
    Skew[1][0] = vect[2];
    Skew[2][0] = -vect[1];
    Skew[0][1] = -vect[2];
    Skew[1][1] = 0;
    Skew[2][1] = vect[0];
    Skew[0][2] = vect[1];
    Skew[1][2] = -vect[0];
    Skew[2][2] = 0;
    return Skew;
};

boost::multi_array<dReal, 2> Eye(){
    boost::multi_array<dReal, 2> I(boost::extents[3][3]);
    I[0][0] = 1;
    I[1][0] = 0;
    I[2][0] = 0;
    I[0][1] = 0;
    I[1][1] = 1;
    I[2][1] = 0;
    I[0][2] = 0;
    I[1][2] = 0;
    I[2][2] = 1;
    return I;
};

std::vector<dReal> AddVect(const std::vector<dReal>&a, const std::vector<dReal>&b,dReal coefa, dReal coefb  ){
    std::vector<dReal> c(3);
    c[0] = a[0]*coefa + b[0]*coefb;
    c[1] = a[1]*coefa + b[1]*coefb;
    c[2] = a[2]*coefa + b[2]*coefb;
    return c;
};


std::vector<dReal> Cross(const std::vector<dReal>&a, const std::vector<dReal>&b ){
    std::vector<dReal> c(3);
    c[0] = a[1]*b[2] - a[2]*b[1];
    c[1] = a[2]*b[0] - a[0]*b[2];
    c[2] = a[0]*b[1] - a[1]*b[0];
    return c;
};

dReal DotVect(const std::vector<dReal>&a, const std::vector<dReal>&b ){
    dReal c = a[0]*b[0] + a[1]*b[1] + a[2]*b[2];
    return c;
};

boost::multi_array<dReal, 2> MatrixAdd(const boost::multi_array<dReal, 2>&A, const boost::multi_array<dReal, 2>&B, dReal coefA, dReal coefB) {
    boost::multi_array<dReal, 2> C(boost::extents[3][3]);
    for (int i = 0; i < int(A.shape()[0]); i++) {
        for (int j = 0; j < int(A.shape()[1]); j++) {
            C[i][j] = coefA*A[i][j] + coefB*B[i][j];
        }
    }
    return C;
}

std::vector<dReal> MatrixMultVector(const boost::multi_array<dReal, 2>&M, const std::vector<dReal>&v) {
    BOOST_ASSERT(M.shape()[1] == v.size());
    std::vector<dReal> res(3);

    for (int i = 0; i < int(M.shape()[0]); i++) {
        res[i] = 0;
        for (int j = 0; j < int(v.size()); j++) {
            res[i] += M[i][j]*v[j];
        }
    }
    return res;
}

boost::multi_array<dReal, 2> MatricesMult3(const boost::multi_array<dReal, 2>&A, const boost::multi_array<dReal, 2>&B) {
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

}
