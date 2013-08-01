#include "TOPP.h"


namespace TOPP {

class KinematicLimits : public Constraints {
public:
    KinematicLimits() : Constraints(){

    }
    std::vector<dReal> amax, vmax;
    void Preprocess(Trajectory& trajectory, const Tunings& tunings);
    std::pair<dReal,dReal> SddLimits(dReal s, dReal sd);
    dReal SdLimitMVC(dReal s);
    void FindSingularSwitchPoints();

};
}
