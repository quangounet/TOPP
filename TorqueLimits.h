#include "TOPP.h"


namespace TOPP {

class TorqueLimits : public Constraints {
public:
    TorqueLimits() : Constraints(){
    }
    TorqueLimits(const std::string& constraintsstring);
    std::vector<dReal> taumin, taumax;
    std::vector<std::vector<dReal> > avect, bvect, cvect;
    std::pair<dReal,dReal> SddLimits(dReal s, dReal sd);
    int dimension;
    void Interpolate(dReal s, std::vector<dReal>& a, std::vector<dReal>& b, std::vector<dReal>& c);
    void DiscretizeDynamics();
    dReal SdLimitMVC(dReal s);
    void FindSingularSwitchPoints();
};
}
