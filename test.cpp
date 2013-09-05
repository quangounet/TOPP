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

    std::string s="2 \n 2\n 1 1 0 1\n 0 2 0 -1\n 3\n 2\n 11 13 6 0.1666666666666\n -4 -10 -6 0.5";

    PiecewisePolynomialTrajectory* ptrajectory;
    ptrajectory = new PiecewisePolynomialTrajectory(s);

    std::stringstream ss;
    ptrajectory->Write(ss);
    std::cout << "Trajectory:\n" << ss.str() << "\n";

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



    PiecewisePolynomialTrajectory newtrajectory;
    ptrajectory->Reparameterize2(resprofileslist,0.01,newtrajectory);

    std::cout << newtrajectory.duration << " " << ptrajectory->duration << "\n";

    std::stringstream ss2;
    newtrajectory.Write(ss2);
    std::cout << "Trajectory:\n" << ss2.str() << "\n";




    delete ptrajectory;

}

