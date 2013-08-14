#include "KinematicLimits.h"
#include "PiecewisePolynomialTrajectory.h"



// Testing

using namespace TOPP;



void PrintVector1d(const std::vector<dReal>& v){
    std::vector<dReal>::const_iterator it = v.begin();
    int count =0;
    while(it != v.end()) {
        std::cout << count << ": " << *it << "\n";
        it++;
        count++;
    }
}

void PrintPair(const std::pair<dReal,dReal>& p){
    std::cout << "(" << p.first << "," << p.second << ")"  << "\n";
}


int main(){

    std::vector<dReal> coefficientsvector;
    std::vector<Polynomial> polynomialsvector;

    coefficientsvector.resize(0);
    coefficientsvector.push_back(1);
    coefficientsvector.push_back(1);
    coefficientsvector.push_back(0);
    coefficientsvector.push_back(1);
    Polynomial P0(coefficientsvector);

    coefficientsvector.resize(0);
    coefficientsvector.push_back(0);
    coefficientsvector.push_back(2);
    coefficientsvector.push_back(0);
    coefficientsvector.push_back(-1);
    Polynomial P1(coefficientsvector);

    coefficientsvector.resize(0);
    coefficientsvector.push_back(11);
    coefficientsvector.push_back(13);
    coefficientsvector.push_back(6);
    coefficientsvector.push_back(1/6.);
    Polynomial P2(coefficientsvector);

    coefficientsvector.resize(0);
    coefficientsvector.push_back(-4);
    coefficientsvector.push_back(-10);
    coefficientsvector.push_back(-6);
    coefficientsvector.push_back(0.5);
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

    ptrajectory->Evaldd(1.9999,q);
    //    std::cout << q[0] << "," << q[1] << "\n";
    ptrajectory->Evaldd(2.0001,q);
    //    std::cout << q[0] << "," << q[1] << "\n";



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

    //    std::cout << ptrajectory->dimension << "\n";
    kinconstraints.Preprocess(*ptrajectory,tunings);

    std::cout << "Switch: " << kinconstraints.switchpointslist.size() << "\n";

    //    PrintVector1d(kinconstraints.mvc);

    //    std::pair<dReal,dReal> p = kinconstraints.SddLimits(0.1,10);
    //    std::cout << "(" << p.first << "," << p.second << ")"  << "\n";


    Profile profile;
    std::cout << "Integrate: " << IntegrateForward(kinconstraints,0,1,profile) << "\n";
    std::cout << "Duration: " << profile.duration << "\n";

    dReal tres;
    int startindex = 0;
    profile.currentindex = profile.nsteps;
    assert(profile.Invert(0.11,tres,true));

    std::cout << "Invert: " << tres << "\n";
    std::cout << "Invert back: " << profile.Eval(tres) << "\n";

    //PrintVector1d(profile.sdvect);

    PiecewisePolynomialTrajectory newtrajectory;
    ptrajectory->Reparameterize(profile,newtrajectory);

    std::cout << newtrajectory.duration << "\n";

    // std::cout << "\n\n\nq\n";
    // for(dReal t=0; t<newtrajectory.duration; t+=0.001) {
    //     newtrajectory.Eval(t,q);
    //     std::cout<< "--\n";
    //     PrintVector1d(q);
    // }
    // std::cout << "\n\n\nqd\n";
    // for(dReal t=0; t<newtrajectory.duration; t+=0.001) {
    //     newtrajectory.Evald(t,q);
    //     std::cout<< "--\n";
    //     PrintVector1d(q);
    // }
    // std::cout << "\n\n\nqdd\n";
    // for(dReal t=0; t<newtrajectory.duration; t+=0.001) {
    //     newtrajectory.Evaldd(t,q);
    //     std::cout<< "--\n";
    //     PrintVector1d(q);
    // }


    delete ptrajectory;

}

