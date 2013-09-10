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
        cvect.push_back(tmpvect);
    }
}



void TorqueLimits::DiscretizeDynamics(){

}



void TorqueLimits::Interpolate(dReal s, std::vector<dReal>& a, std::vector<dReal>& b, std::vector<dReal>& c){
    a.resize(trajectory.dimension);
    b.resize(trajectory.dimension);
    c.resize(trajectory.dimension);
    assert(s>=-TINY && s<=trajectory.duration+TINY);
    if(s<0) {
        s=0;
    }
    if(s>=trajectory.duration) {
        int n = ndiscrsteps-1;
        for(int i=0; i<trajectory.dimension; i++) {
            a[i]= avect[n][i];
            b[i]= bvect[n][i];
            c[i]= cvect[n][i];
        }
        return;
    }
    int n = (int) s/tunings.discrtimestep;
    dReal coef = (s-n*tunings.discrtimestep)/tunings.discrtimestep;
    for(int i=0; i<trajectory.dimension; i++) {
        a[i] = (1-coef)*avect[n][i] + coef*avect[n+1][i];
        b[i] = (1-coef)*bvect[n][i] + coef*bvect[n+1][i];
        c[i] = (1-coef)*cvect[n][i] + coef*cvect[n+1][i];
    }
}



std::pair<dReal,dReal> TorqueLimits::SddLimits(dReal s, dReal sd){
    dReal alpha = -INF;
    dReal beta = INF;
    dReal sdsq = sd*sd;
    std::vector<dReal> a, b, c;

    dReal taumin_i, taumax_i, alpha_i, beta_i;
    Interpolate(s,a,b,c);

    for(int i=0; i<trajectory.dimension; i++) {
        if(a[i]>0) {
            taumin_i = taumin[i];
            taumax_i = taumax[i];
        }
        else{
            taumin_i = taumax[i];
            taumax_i = taumin[i];
        }
        alpha_i = (taumin_i-sdsq*b[i]-c[i])/a[i];
        beta_i = (taumax_i-sdsq*b[i]-c[i])/a[i];
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
    std::vector<dReal> tau_alpha(trajectory.dimension), tau_beta(trajectory.dimension);
    std::vector<dReal> a, b, c;
    Interpolate(s,a,b,c);

    for(int i=0; i<trajectory.dimension; i++) {
        if(a[i] > 0) {
            tau_alpha[i] = taumin[i];
            tau_beta[i] = taumax[i];
        }
        else{
            tau_alpha[i] = taumax[i];
            tau_beta[i] = taumin[i];
        }
    }
    dReal sdmin = INF;
    for(int k=0; k<trajectory.dimension; k++) {
        for(int m=k+1; m<trajectory.dimension; m++) {
            dReal num, denum, r;
            num = a[k]*(tau_alpha[m]-c[m])-a[m]*(tau_beta[k]-c[k]);
            denum = a[k]*b[m]-a[m]*b[k];
            if(abs(denum) >= TINY) {
                r = num/denum;
                if(r>=0) {
                    sdmin = std::min(sdmin,sqrt(r));
                }
            }
            num = a[m]*(tau_alpha[k]-c[k])-a[k]*(tau_beta[m]-c[m]);
            denum = -denum;
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
    std::vector<dReal> a,aprev,b,c;
    Interpolate(discrsvect[i],aprev,b,c);

    for(int i=1; i<ndiscrsteps-1; i++) {
        Interpolate(discrsvect[i],a,b,c);
        for(int j=0; j<trajectory.dimension; j++) {
            if(a[j]*aprev[j]<0) {
                AddSwitchPoint(i,SP_SINGULAR);
                continue;
            }
        }
        aprev = a;
    }
}
}
