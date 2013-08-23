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
    tunings.integrationtimestep = 0.01;
    tunings.threshold = 0.01;
    tunings.passswitchpointnsteps = 10;

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

    std::list<SwitchPoint>::iterator itsw = kinconstraints.switchpointslist.begin();
    while(itsw!=kinconstraints.switchpointslist.end()) {
        std::cout << "Type " << itsw->switchpointtype << ": (" << itsw->s <<"," << itsw->sd << ")\n";
        itsw++;
    }


    //    PrintVector1d(kinconstraints.mvc);

    //    std::pair<dReal,dReal> p = kinconstraints.SddLimits(0.1,10);
    //    std::cout << "(" << p.first << "," << p.second << ")"  << "\n";


    Profile profile;
    std::list<Profile> resprofileslist;
    ComputeLimitingCurves(kinconstraints,resprofileslist);
    IntegrateForward(kinconstraints,0,1e-4,profile,1e5,resprofileslist);
    resprofileslist.push_back(profile);
    IntegrateBackward(kinconstraints,ptrajectory->duration,1e-4,profile,1e5,resprofileslist);
    resprofileslist.push_back(profile);


    std::cout << "CLC: " << resprofileslist.size() << "\n";

    std::list<Profile>::iterator itprof = resprofileslist.begin();
    while(itprof != resprofileslist.end()) {
        std::cout << itprof->nsteps << "\n";
        if(itprof->nsteps > 10000) {
            //itprof->Print();
        }
        itprof++;
    }

    Profile profilefinal(resprofileslist,0.1);

    std::cout << "\n\n\n";

    dReal dt = 0.001;

    std::list<Profile>::iterator itp = resprofileslist.begin();
    itp++;
    itp++;
    itp++;
    itp++;
    Profile p = *itp;
    dReal scur,sdcur,sdd,s,sd;

    scur = p.Eval(0);
    sdcur = p.Evald(0);

    // for(dReal t=0; t<0.1; t+=dt) {
    //     s = p.Eval(t);
    //     sd = p.Evald(t);
    //     std::cout << s << "/"  << scur << " **  "<< sd << "/" << sdcur << "\n";
    //     sdd = p.Evaldd(t);
    //     scur += sdcur*dt + sdd*dt*dt*0.5;
    //     sdcur += sdd*dt;
    // }


    std::cout << "Profilefinal: " << profilefinal.duration << "\n";

    //PrintVector1d(profile.sdvect);

    PiecewisePolynomialTrajectory newtrajectory;
    ptrajectory->Reparameterize(profilefinal,newtrajectory);
    PiecewisePolynomialTrajectory newtrajectory2;
    newtrajectory.Reintegrate(0.001,newtrajectory2);



    std::cout << "\n\n\nq\n";
    for(dReal t=0; t<newtrajectory2.duration; t+=0.1) {
        newtrajectory2.Eval(t,q);
        std::cout<< "--\n";
        PrintVector1d(q);
    }
    std::cout << "\n\n\nqd\n";
    for(dReal t=0; t<newtrajectory2.duration; t+=0.1) {
        newtrajectory2.Evald(t,q);
        std::cout<< "--\n";
        PrintVector1d(q);
    }
    std::cout << "\n\n\nqdd\n";
    for(dReal t=0; t<newtrajectory2.duration; t+=0.1) {
        newtrajectory2.Evaldd(t,q);
        std::cout<< "--\n";
        PrintVector1d(q);
    }


    // std::list<Chunk> reslist;
    // std::vector<dReal> q0(2), qd0(2), qdd0(2), q1(2), qd1(2), qdd1(2);
    // ptrajectory->Eval(0,q0);
    // ptrajectory->Evald(0,qd0);
    // ptrajectory->Evaldd(0,qdd0);
    // ptrajectory->Eval(ptrajectory->duration,q1);
    // ptrajectory->Evald(ptrajectory->duration,qd1);
    // ptrajectory->Evaldd(ptrajectory->duration,qdd1);

    // Interpolate3(ptrajectory->duration,q0,qd0,qdd0,q1,qd1,qdd1,reslist);
    // PiecewisePolynomialTrajectory newtrajectory2(reslist);

    // std::vector<dReal>  qn(2), qd(2);
    // dReal dt = 0.01;
    // dReal diff;

    // for(dReal t=0; t<0.5; t+=dt) {
    //     newtrajectory.Evald(t,q);
    //     newtrajectory.Evald(t+dt,qn);
    //     newtrajectory.Evaldd(t+dt,qd);
    //     diff = (qn[0]-q[0])/dt;
    //     std::cout << diff << "----" << qd[0] << "\n";
    // }

    // std::cout << "\n\n\nqd\n";
    // for(dReal t=0; t<0.1; t+=0.01) {
    //     newtrajectory.Evald(t,q);
    //     PrintVector1d(q);
    //     std::cout << "----\n";
    //     newtrajectory2.Evald(t,q);
    //     PrintVector1d(q);
    //     std::cout << "\n\n";
    // }

    // std::cout << "\n\n\nqdd\n";
    // for(dReal t=0; t<newtrajectory.duration; t+=0.5) {
    //     newtrajectory.Evaldd(t,q);
    //     PrintVector1d(q);
    //     std::cout << "----\n";
    //     newtrajectory2.Evaldd(t,q);
    //     PrintVector1d(q);
    //     std::cout << "\n\n";
    // }


    // std::cout << newtrajectory2.duration << "\n";


    delete ptrajectory;

}

