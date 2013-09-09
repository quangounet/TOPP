#include "TOPP.h"


namespace TOPP {

class KinematicLimits : public Constraints {
public:
    KinematicLimits() : Constraints(){
    }
    KinematicLimits(const std::string& constraintsstring);
    std::vector<dReal> amax, vmax;
    std::pair<dReal,dReal> SddLimits(dReal s, dReal sd);
    void DiscretizeDynamics(){
    };
    dReal SdLimitMVC(dReal s);
    void FindSingularSwitchPoints();
};
}
