#include "KinematicLimits.h"



// Testing

using namespace TOPP;

int main(){
    Trajectory trajectory;
    Tunings tunings;
    KinematicLimits kinconstraints;

    tunings.discrtimestep = 0.01;
    tunings.integrationtimestep = 0.001;

    std::vector<dReal> amax, vmax;

    kinconstraints.amax = amax;
    kinconstraints.vmax = vmax;
    kinconstraints.Preprocess(trajectory,tunings);
}

