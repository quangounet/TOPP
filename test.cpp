#include "KinematicLimits.h"
#include "PiecewisePolynomialTrajectory.h"



// Testing

using namespace TOPP;

int main(){

    std::vector<dReal> coefficientsvector;
    coefficientsvector.push_back(0);
    coefficientsvector.push_back(0);
    coefficientsvector.push_back(1);
    Polynomial p(coefficientsvector);
    std::vector<Polynomial> polynomialsvector;
    polynomialsvector.push_back(p);
    Chunk chunk(1.,polynomialsvector);

    std::vector<dReal> coefficientsvector2;
    coefficientsvector2.push_back(0);
    coefficientsvector2.push_back(0);
    coefficientsvector2.push_back(2);
    Polynomial p2(coefficientsvector2);
    std::vector<Polynomial> polynomialsvector2;
    polynomialsvector2.push_back(p2);
    Chunk chunk2(2.,polynomialsvector2);

    std::vector<dReal> coefficientsvector3;
    coefficientsvector3.push_back(0);
    coefficientsvector3.push_back(0);
    coefficientsvector3.push_back(3);
    Polynomial p3(coefficientsvector3);
    std::vector<Polynomial> polynomialsvector3;
    polynomialsvector3.push_back(p3);
    Chunk chunk3(3.,polynomialsvector3);

    std::list<Chunk> chunkslist;
    chunkslist.push_back(chunk);
    chunkslist.push_back(chunk2);
    chunkslist.push_back(chunk3);

    PiecewisePolynomialTrajectory trajectory(chunkslist);

    std::vector<dReal> q(1), qd(1), qdd(1);
    dReal s=6.1;
    trajectory.Eval(s,q);
    //trajectory.Evald(s,qd);
    //trajectory.Evaldd(s,qdd);

    std::cout << q[0] << "\n";
    //std::cout << qd[0] << "\n";
    //std::cout << qdd[0] << "\n";


    // Trajectory trajectory;
    // Tunings tunings;
    // KinematicLimits kinconstraints;

    // tunings.discrtimestep = 0.01;
    // tunings.integrationtimestep = 0.001;

    // std::vector<dReal> amax, vmax;


    // kinconstraints.amax = amax;
    // kinconstraints.vmax = vmax;
    //kinconstraints.Preprocess(trajectory,tunings);
}

