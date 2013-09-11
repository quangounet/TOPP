#include "KinematicLimits.h"


namespace TOPP {

KinematicLimits::KinematicLimits(const std::string& constraintsstring){
    int buffsize = 255;
    char buff[buffsize];
    std::istringstream iss(constraintsstring);
    iss.getline(buff,buffsize);
    VectorFromString(std::string(buff),amax);
    iss.getline(buff,buffsize);
    VectorFromString(std::string(buff),vmax);
}



std::pair<dReal,dReal> KinematicLimits::SddLimits(dReal s, dReal sd){
    dReal alpha = -INF;
    dReal beta = INF;
    dReal sdsq = sd*sd;
    dReal a_alpha_i, a_beta_i, alpha_i, beta_i;
    std::vector<dReal> qd(trajectory.dimension), qdd(trajectory.dimension);
    trajectory.Evald(s, qd);
    trajectory.Evaldd(s, qdd);
    for(int i=0; i<trajectory.dimension; i++) {
        if(qd[i]>0) {
            a_alpha_i = -amax[i];
            a_beta_i = amax[i];
        }
        else{
            a_alpha_i = amax[i];
            a_beta_i = -amax[i];
        }
        alpha_i = (a_alpha_i-sdsq*qdd[i])/qd[i];
        beta_i = (a_beta_i-sdsq*qdd[i])/qd[i];
        alpha = std::max(alpha,alpha_i);
        beta = std::min(beta,beta_i);
    }
    std::pair<dReal,dReal> result(alpha,beta);
    return result;
}



dReal KinematicLimits::SdLimitMVC(dReal s){
    std::pair<dReal,dReal> sddlimits = KinematicLimits::SddLimits(s,0);
    if(sddlimits.first > sddlimits.second) {
        return 0;
    }
    std::vector<dReal> a_alpha(trajectory.dimension), a_beta(trajectory.dimension);
    std::vector<dReal> qd(trajectory.dimension), qdd(trajectory.dimension);
    trajectory.Evald(s, qd);
    trajectory.Evaldd(s, qdd);
    for(int i=0; i<trajectory.dimension; i++) {
        if(qd[i] > 0) {
            a_alpha[i] = -amax[i];
            a_beta[i] = amax[i];
        }
        else{
            a_alpha[i] = amax[i];
            a_beta[i] = -amax[i];
        }
    }
    dReal sdmin = INF;
    for(int k=0; k<trajectory.dimension; k++) {
        for(int m=k+1; m<trajectory.dimension; m++) {
            dReal num, denum, r;
            num = qd[m]*a_alpha[k]-qd[k]*a_beta[m];
            denum = qd[m]*qdd[k]-qd[k]*qdd[m];
            if(std::abs(denum) >= TINY) {
                r = num/denum;
                if(r>=0) {
                    sdmin = std::min(sdmin,sqrt(r));
                }
            }
            num = qd[k]*a_alpha[m]-qd[m]*a_beta[k];
            denum = -denum;
            if(std::abs(denum) >= TINY) {
                r = num/denum;
                if(r>=0) {
                    sdmin = std::min(sdmin,sqrt(r));
                }
            }
        }
    }
    return sdmin;
}


void KinematicLimits::FindSingularSwitchPoints(){
    if(ndiscrsteps<3) {
        return;
    }
    int i = 0;
    std::vector<dReal> qd(trajectory.dimension),qdprev(trajectory.dimension);
    trajectory.Evald(discrsvect[i],qdprev);

    for(int i=1; i<ndiscrsteps-1; i++) {
        trajectory.Evald(discrsvect[i],qd);
        for(int j=0; j<trajectory.dimension; j++) {
            if(qd[j]*qdprev[j]<0) {
                AddSwitchPoint(i,SP_SINGULAR);
                continue;
            }
        }
        qdprev = qd;
    }
}
}
