#include "TorqueLimits.h"


namespace TOPP {

TorqueLimits::TorqueLimits(const std::string& constraintsstring){
    int buffsize = 255;
    std::vector<dReal> tmpvect;
    char buff[buffsize];
    std::istringstream iss(constraintsstring);
    iss.getline(buff,buffsize);
    VectorFromString(std::string(buff),taumin);
    iss.getline(buff,buffsize);
    VectorFromString(std::string(buff),taumax);
    while(iss.good()) {
        iss.getline(buff,buffsize);
        VectorFromString(std::string(buff),tmpvect);
        avect.push_back(tmpvect);
        iss.getline(buff,buffsize);
        VectorFromString(std::string(buff),tmpvect);
        bvect.push_back(tmpvect);
        iss.getline(buff,buffsize);
        VectorFromString(std::string(buff),tmpvect);
        bvect.push_back(tmpvect);
    }
}



void TorqueLimits::DiscretizeDynamics(){

}



void TorqueLimits::Interpolate(dReal s, std::vector<dReal>& a, std::vector<dReal>& b, std::vector<dReal>& c){
    assert(s>=0 && s<=trajectory.duration);
    int n = (int) s/tunings.discrtimestep;
    dReal coef = (s-n*tunings.discrtimestep)/tunings.discrtimestep;
    a.resize(trajectory.dimension);
    b.resize(trajectory.dimension);
    c.resize(trajectory.dimension);
    for(int i=0; i<trajectory.dimension; i++) {
        a[i] = (1-coef)*avect[n][i] + coef*avect[n+1][i];
        b[i] = (1-coef)*bvect[n][i] + coef*bvect[n+1][i];
        c[i] = (1-coef)*cvect[n][i] + coef*cvect[n+1][i];
    }
}




std::pair<dReal,dReal> TorqueLimits::SddLimits(dReal s, dReal sd){
    dReal alpha = -INF;
    dReal beta = INF;
    std::vector<dReal> a, b, c;
    Interpolate(s,a,b,c);


    dReal a_alpha_i, a_beta_i, alpha_i, beta_i;
    std::vector<dReal> qd(trajectory.dimension), qdd(trajectory.dimension);
    trajectory.Evald(s, qd);
    trajectory.Evaldd(s, qdd);
    for(int i=0; i<trajectory.dimension; i++) {
        if(qd[i]>0) {
            //            a_alpha_i = -taumax[i];
            //a_beta_i = amax[i];
        }
        else{
            //a_alpha_i = amax[i];
            //a_beta_i = -amax[i];
        }
        alpha_i = (a_alpha_i-sd*sd*qdd[i])/qd[i];
        beta_i = (a_beta_i-sd*sd*qdd[i])/qd[i];
        alpha = std::max(alpha,alpha_i);
        beta = std::min(beta,beta_i);
    }
    std::pair<dReal,dReal> result(alpha,beta);
    return result;
}



dReal TorqueLimits::SdLimitMVC(dReal s){
    std::pair<dReal,dReal> sddlimits = TorqueLimits::SddLimits(s,0);
    if(sddlimits.first > sddlimits.second) {
        return 0;
    }
    std::vector<dReal> a_alpha(trajectory.dimension), a_beta(trajectory.dimension);
    std::vector<dReal> qd(trajectory.dimension), qdd(trajectory.dimension);
    trajectory.Evald(s, qd);
    trajectory.Evaldd(s, qdd);
    for(int i=0; i<trajectory.dimension; i++) {
        if(qd[i] > 0) {
            //            a_alpha[i] = -amax[i];
            //a_beta[i] = amax[i];
        }
        else{
            //a_alpha[i] = amax[i];
            //a_beta[i] = -amax[i];
        }
    }
    dReal sdmin = INF;
    for(int k=0; k<trajectory.dimension; k++) {
        for(int m=k+1; m<trajectory.dimension; m++) {
            dReal num, denum, r;
            num = qd[m]*a_alpha[k]-qd[k]*a_beta[m];
            denum = qd[m]*qdd[k]-qd[k]*qdd[m];
            if(abs(denum) >= TINY) {
                r = num/denum;
                if(r>=0) {
                    sdmin = std::min(sdmin,sqrt(r));
                }
            }
            num = qd[k]*a_alpha[m]-qd[m]*a_beta[k];
            denum = qd[k]*qdd[m]-qd[m]*qdd[k];
            if(abs(denum) >= TINY) {
                r = num/denum;
                if(r>=0) {
                    sdmin = std::min(sdmin,sqrt(r));
                }
            }
        }
    }
    return sdmin;
}


void TorqueLimits::FindSingularSwitchPoints(){
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
