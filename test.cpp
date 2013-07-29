#include "KinematicLimits.h"
#include "PiecewisePolynomialTrajectory.h"



// Testing

using namespace TOPP;

int main(){

    std::vector<dReal> coefficientsvector;
    std::vector<Polynomial> polynomialsvector;

    coefficientsvector.resize(0);
    coefficientsvector.push_back(0);
    coefficientsvector.push_back(0);
    coefficientsvector.push_back(1);
    Polynomial P0(coefficientsvector);
    polynomialsvector.resize(0);
    polynomialsvector.push_back(P0);
    Chunk chunk0(1.,polynomialsvector);

    coefficientsvector.resize(0);
    coefficientsvector.push_back(1);
    coefficientsvector.push_back(0);
    coefficientsvector.push_back(-2);
    Polynomial P1(coefficientsvector);
    polynomialsvector.resize(0);
    polynomialsvector.push_back(P1);
    Chunk chunk1(2.,polynomialsvector);

    coefficientsvector.resize(0);
    coefficientsvector.push_back(-7);
    coefficientsvector.push_back(2);
    coefficientsvector.push_back(0);
    Polynomial P2(coefficientsvector);
    polynomialsvector.resize(0);
    polynomialsvector.push_back(P2);
    Chunk chunk2(0.5,polynomialsvector);

    coefficientsvector.resize(0);
    coefficientsvector.push_back(-6);
    coefficientsvector.push_back(0);
    coefficientsvector.push_back(1);
    Polynomial P3(coefficientsvector);
    polynomialsvector.resize(0);
    polynomialsvector.push_back(P3);
    Chunk chunk3(2.,polynomialsvector);


    std::list<Chunk> chunkslist;
    chunkslist.push_back(chunk0);
    chunkslist.push_back(chunk1);
    chunkslist.push_back(chunk2);
    chunkslist.push_back(chunk3);

    PiecewisePolynomialTrajectory trajectory(chunkslist);

    std::vector<dReal> q(1), qd(1), qdd(1);


    std::list<dReal> slist, sdlist, sddlist;
    slist.push_back(0);
    slist.push_back(0.5);
    slist.push_back(4);
    slist.push_back(5.5);

    sdlist.push_back(5);
    sdlist.push_back(35);
    sdlist.push_back(15);
    sdlist.push_back(0);

    sddlist.push_back(0);
    sddlist.push_back(0);
    sddlist.push_back(0);
    sddlist.push_back(0);
    Profile profile(slist,sdlist,sddlist,0.1);


    //std::cout << profile.Eval(0.4) << "\n";

    trajectory.Reparameterize(profile);

    std::cout << trajectory.chunkslist.size() << "\n";

    std::list<Chunk>::iterator it = trajectory.chunkslist.begin();
    while(it!=trajectory.chunkslist.end()) {
        std::cout << "Chunk i: " << it->duration << "\n";
        it++;
    }


    trajectory.Eval(0.0,q);
    std::cout << "0.0: " << q[0] << "\n";

    trajectory.Eval(0.1,q);
    std::cout << "0.1: " << q[0] << "\n";

    trajectory.Eval(0.2,q);
    std::cout << "0.2: " << q[0] << "\n";

    trajectory.Eval(0.3,q);
    std::cout << "0.3: " << q[0] << "\n";





    //dReal sol;
    //assert(SolveQuadraticEquation(3,-4,0,0,1,sol));
    //std::cout << sol << "\n";

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

