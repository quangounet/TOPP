#include "KinematicLimits.h"
#include "PiecewisePolynomialTrajectory.h"



// Testing

using namespace TOPP;

int main(){

    std::vector<dReal> coefficientsvector;
    std::vector<Polynomial> polynomialsvector;

    coefficientsvector.resize(0);
    coefficientsvector.push_back(1);
    coefficientsvector.push_back(0);
    coefficientsvector.push_back(0);
    coefficientsvector.push_back(1);
    Polynomial P0(coefficientsvector);

    coefficientsvector.resize(0);
    coefficientsvector.push_back(0);
    coefficientsvector.push_back(0);
    coefficientsvector.push_back(0);
    coefficientsvector.push_back(-1);
    Polynomial P1(coefficientsvector);

    coefficientsvector.resize(0);
    coefficientsvector.push_back(9);
    coefficientsvector.push_back(12);
    coefficientsvector.push_back(6);
    coefficientsvector.push_back(1/3.);
    Polynomial P2(coefficientsvector);

    coefficientsvector.resize(0);
    coefficientsvector.push_back(-8);
    coefficientsvector.push_back(-12);
    coefficientsvector.push_back(-6);
    coefficientsvector.push_back(1/6.);
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

    PiecewisePolynomialTrajectory* ptrajectory;
    ptrajectory = new PiecewisePolynomialTrajectory(chunkslist);

    std::vector<dReal> q(2);

    ptrajectory->Eval(1.99,q);
    std::cout << q[0] << "," << q[1] << "\n";
    ptrajectory->Eval(2.01,q);
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

    std::cout << ptrajectory->dimension << "\n";
    kinconstraints.Preprocess(*ptrajectory,tunings);

    std::cout << "SP: " << kinconstraints.mvc.size() << "\n";

    // std::vector<dReal>::iterator it = kinconstraints.mvc.begin();
    // int toto =0;
    // while(it != kinconstraints.mvc.end()) {
    //     std::cout << toto << ": " << *it << "\n";
    //     it++;
    //     toto++;
    // }

    delete ptrajectory;

}

