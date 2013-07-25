#include "TOPP.h"


using namespace TOPP;

class TorqueLimits : public Constraints {
public:
    TorqueLimits(Trajectory& trajectory, Tunings& tunings) : Constraints(trajectory,tunings){
    }

    void DiscretizeDynamics();
    void ComputeMVC();
    void ComputeSwitchPoints();
};
