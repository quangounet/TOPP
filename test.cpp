#include "KinematicLimits.h"
#include "PiecewisePolynomialTrajectory.h"



// Testing

using namespace TOPP;

int main(){

    std::vector<dReal> coefficientsvector;
    std::vector<Polynomial> polynomialsvector;

    coefficientsvector.resize(0);
    coefficientsvector.push_back(1);
    coefficientsvector.push_back(-3);
    coefficientsvector.push_back(1);
    Polynomial P0(coefficientsvector);

    coefficientsvector.resize(0);
    coefficientsvector.push_back(1);
    coefficientsvector.push_back(1);
    coefficientsvector.push_back(-2);
    Polynomial P1(coefficientsvector);

    coefficientsvector.resize(0);
    coefficientsvector.push_back(-1);
    coefficientsvector.push_back(-1);
    coefficientsvector.push_back(2);
    Polynomial P2(coefficientsvector);

    coefficientsvector.resize(0);
    coefficientsvector.push_back(-5);
    coefficientsvector.push_back(-7);
    coefficientsvector.push_back(1);
    Polynomial P3(coefficientsvector);

    polynomialsvector.resize(0);
    polynomialsvector.push_back(P0);
    polynomialsvector.push_back(P1);
    Chunk chunk0(2.,polynomialsvector);

    polynomialsvector.resize(0);
    polynomialsvector.push_back(P2);
    polynomialsvector.push_back(P3);
    Chunk chunk1(3,polynomialsvector);

    std::list<Chunk> chunkslist;
    chunkslist.push_back(chunk0);
    chunkslist.push_back(chunk1);

    PiecewisePolynomialTrajectory trajectory(chunkslist);

    std::vector<dReal> q(2);

    trajectory.Evald(0,q);
    std::cout << q[0] << "," << q[1] << "\n";
    trajectory.Evald(1,q);
    std::cout << q[0] << "," << q[1] << "\n";
    trajectory.Evald(2,q);
    std::cout << q[0] << "," << q[1] << "\n";
    trajectory.Evald(3,q);
    std::cout << q[0] << "," << q[1] << "\n";
    trajectory.Evald(4,q);
    std::cout << q[0] << "," << q[1] << "\n";
    trajectory.Evald(5,q);
    std::cout << q[0] << "," << q[1] << "\n";



    Tunings tunings;
    KinematicLimits kinconstraints;

    tunings.discrtimestep = 0.01;
    tunings.integrationtimestep = 0.001;

    std::vector<dReal> amax, vmax;


    amax.push_back(1);
    amax.push_back(2);
    vmax.push_back(0);
    vmax.push_back(0);

    kinconstraints.amax = amax;
    kinconstraints.vmax = vmax;

    kinconstraints.Preprocess(trajectory,tunings);

    kinconstraints.Preprocess(trajectory,tunings);

}

