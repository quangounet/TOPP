#include "TOPP.h"


namespace TOPP {

class KinematicLimits : public Constraints {
public:
    KinematicLimits() : Constraints(){
    }
    KinematicLimits(const std::string& constraintsstring);
    std::vector<dReal> amax, vmax;
    void Preprocess(Trajectory& trajectory, const Tunings& tunings);
    std::pair<dReal,dReal> SddLimits(dReal s, dReal sd);
    dReal SdLimitMVC(dReal s);
    void FindSingularSwitchPoints();

};
}
