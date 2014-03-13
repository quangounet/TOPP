// -*- coding: utf-8 -*-
// Copyright (C) 2013 Quang-Cuong Pham <cuong.pham@normalesup.org> & Rosen Diankov <rosen.diankov@gmail.com>
//
// This file is part of the Time-Optimal Path Parameterization (TOPP) library.
// TOPP is free software: you can redistribute it and/or modify
// it under the terms of the GNU Lesser General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// at your option) any later version.
//
// This program is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU Lesser General Public License for more details.
//
// You should have received a copy of the GNU Lesser General Public License
// along with this program.  If not, see <http://www.gnu.org/licenses/>.
#ifdef WITH_OPENRAVE

#include "TorqueLimitsRave.h"

using namespace OpenRAVE;

namespace TOPP {

TorqueLimitsRave::TorqueLimitsRave(RobotBasePtr probot, std::string& constraintsstring, Trajectory* ptraj){
    trajectory = *ptraj;
    int ndof = trajectory.dimension;
    std::istringstream iss(constraintsstring);
    iss >> discrtimestep;
    ReadVectorFromStream(iss, ndof, vmax);
    ReadVectorFromStream(iss, ndof, taumin);
    ReadVectorFromStream(iss, ndof, taumax);
    hasvelocitylimits = false;
    FOREACH(itv, vmax) {
        if( std::abs(*itv) > TINY ) {
            hasvelocitylimits = true;
            break;
        }
    }

    // Define the avect, bvect, cvect
    int ndiscrsteps = int((trajectory.duration+1e-10)/discrtimestep)+1;
    std::vector<dReal> q(ndof), qd(ndof), qdd(ndof), tmp0(ndof), tmp1(ndof), torquesimple;
    boost::array< std::vector< dReal >, 3 > torquecomponents;
    {
        avect.resize(ndiscrsteps);
        bvect.resize(ndiscrsteps);
        cvect.resize(ndiscrsteps);
        EnvironmentMutex::scoped_lock lock(probot->GetEnv()->GetMutex()); // lock environment
        for(int i = 0; i<ndiscrsteps; i++) {
            dReal s = i*discrtimestep;
            trajectory.Eval(s,q);
            trajectory.Evald(s,qd);
            trajectory.Evaldd(s,qdd);
            probot->SetDOFValues(q,KinBody::CLA_Nothing);
            probot->SetDOFVelocities(qd,KinBody::CLA_Nothing);
            probot->ComputeInverseDynamics(torquesimple,qd);
            probot->ComputeInverseDynamics(torquecomponents,qdd);
            VectorAdd(torquesimple,torquecomponents[1], bvect[i], 1, -1);
            VectorAdd(bvect[i], torquecomponents[2], avect[i], 1, -1);
            VectorAdd(torquecomponents[0],torquecomponents[1], bvect[i]);
            cvect[i] = torquecomponents[2];
        }
    }
}

void ConvertToTOPPTrajectory(OpenRAVE::TrajectoryBaseConstPtr pintraj, const OpenRAVE::ConfigurationSpecification& posspec, Trajectory& outtraj)
{
    int N = pintraj->GetNumWaypoints();
    if( N < 2 ) {
        throw TOPPException("openrave trajectory has less than 2 waypoints");
    }
    if( pintraj->GetDuration() <= 0 ) {
        throw TOPPException("openrave trajectory is not retimed");
    }
    OpenRAVE::ConfigurationSpecification timespec;
    timespec.AddDeltaTimeGroup();
    std::vector<OpenRAVE::ConfigurationSpecification::Group>::const_iterator itposgroup = pintraj->GetConfigurationSpecification().FindCompatibleGroup(posspec._vgroups.at(0).name, false);
    BOOST_ASSERT(itposgroup!= pintraj->GetConfigurationSpecification()._vgroups.end());

    int degree = 0;
    if( itposgroup->interpolation == "linear" ) {
        degree = 1;
    }
    else if( itposgroup->interpolation == "quadratic" ) {
        degree = 2;
    }
    else if( itposgroup->interpolation == "cubic" ) {
        degree = 3;
    }
    else if( itposgroup->interpolation == "quadric" ) {
        degree = 4;
    }
    else if( itposgroup->interpolation == "quintic" ) {
        degree = 5;
    }
    else {
        throw TOPP_EXCEPTION_FORMAT("unknown interpolation method '%s'", itposgroup->interpolation, 0);
    }

    // go up to accelerations
    std::vector<ConfigurationSpecification> derivspecs(std::min(degree,3));
    derivspecs[0] = posspec;
    for(size_t i = 1; i < derivspecs.size(); ++i) {
        derivspecs[i] = posspec.ConvertToDerivativeSpecification(i);
    }

    std::list<Chunk> listchunks;
    std::vector< std::vector<dReal> > vprevpoint(derivspecs.size()), vnewpoint(derivspecs.size());
    std::vector<dReal> vdeltatime;
    std::vector<Polynomial> polynomialsvector(posspec.GetDOF());
    std::vector<dReal> coefficientsvector(degree+1);
    for(size_t ideriv = 0; ideriv < derivspecs.size(); ++ideriv) {
        pintraj->GetWaypoint(0, vprevpoint[ideriv], derivspecs[ideriv]);
    }

    for(size_t iwaypoint = 1; iwaypoint < pintraj->GetNumWaypoints(); ++iwaypoint) {
        for(size_t ideriv = 0; ideriv < derivspecs.size(); ++ideriv) {
            pintraj->GetWaypoint(iwaypoint, vnewpoint[ideriv], derivspecs[ideriv]);
        }
        pintraj->GetWaypoint(iwaypoint, vdeltatime, timespec);
        dReal deltatime = vdeltatime.at(0);
        dReal ideltatime = 1/deltatime;
        for(int idof = 0; idof < posspec.GetDOF(); ++idof) {
            // try not to use the acceleration
            dReal px = vnewpoint.at(0).at(idof) - vprevpoint.at(0).at(idof);
            if( degree == 2 ) {
                coefficientsvector[0] = vprevpoint.at(0).at(idof);
                coefficientsvector[1] = vprevpoint.at(1).at(idof);
                coefficientsvector[2] = (px*ideltatime - vprevpoint.at(1).at(idof))*ideltatime;
            }
            else if( degree == 3 ) {
                coefficientsvector[0] = vprevpoint.at(0).at(idof);
                coefficientsvector[1] = vprevpoint.at(1).at(idof);
                coefficientsvector[2] = (3*px*ideltatime - 2*vprevpoint.at(1).at(idof) - vnewpoint.at(1).at(idof))*ideltatime;
                coefficientsvector[3] = (-2*px*ideltatime + vprevpoint.at(1).at(idof) + vnewpoint.at(1).at(idof))*ideltatime*ideltatime;
            }
            else {
                throw TOPP_EXCEPTION_FORMAT("do not degree %d yet", degree, 0);
            }
            polynomialsvector[idof].InitFromCoefficientsVector(coefficientsvector);
        }
        listchunks.push_back(Chunk(deltatime, polynomialsvector));
        vprevpoint.swap(vnewpoint);
    }
    outtraj.InitFromChunksList(listchunks);
}

void ConvertToOpenRAVETrajectory(const Trajectory& intraj, OpenRAVE::TrajectoryBasePtr pouttraj, const OpenRAVE::ConfigurationSpecification& posspec)
{
    OpenRAVE::ConfigurationSpecification newposspec = posspec;
    if( intraj.degree == 1  ) {
        FOREACH(itgroup, newposspec._vgroups) {
            itgroup->interpolation = "linear";
        }
    }
    else if( intraj.degree == 2  ) {
        FOREACH(itgroup, newposspec._vgroups) {
            itgroup->interpolation = "quadratic";
        }
    }
    else if( intraj.degree == 3  ) {
        FOREACH(itgroup, newposspec._vgroups) {
            itgroup->interpolation = "cubic";
        }
    }
    else if( intraj.degree == 4  ) {
        FOREACH(itgroup, newposspec._vgroups) {
            itgroup->interpolation = "quadric";
        }
    }
    else if( intraj.degree == 5  ) {
        FOREACH(itgroup, newposspec._vgroups) {
            itgroup->interpolation = "quintic";
        }
    }
    else if( intraj.degree == 5  ) {
        FOREACH(itgroup, newposspec._vgroups) {
            itgroup->interpolation = "quintic";
        }
    }
    else if( intraj.degree == 6  ) {
        FOREACH(itgroup, newposspec._vgroups) {
            itgroup->interpolation = "sextic";
        }
    }
    else {
        FOREACH(itgroup, newposspec._vgroups) {
            itgroup->interpolation = str(boost::format("degree%d")%intraj.degree);
        }
    }

    // go up to accelerations
    std::vector<ConfigurationSpecification> derivspecs(std::min(intraj.degree,3));
    derivspecs[0] = newposspec;
    ConfigurationSpecification totalspec = newposspec;
    for(size_t i = 1; i < derivspecs.size(); ++i) {
        derivspecs[i] = newposspec.ConvertToDerivativeSpecification(i);
        totalspec += derivspecs[i];
    }
    totalspec.AddDeltaTimeGroup();
    
    pouttraj->Init(totalspec);
    std::vector<dReal> v(totalspec.GetDOF(),0);
    std::vector<dReal> q(intraj.dimension), qd(intraj.dimension), qdd(intraj.dimension);
    intraj.Eval(0, q);
    intraj.Evald(0, qd);
    intraj.Evaldd(0, qdd);
    std::copy(q.begin(), q.end(), v.begin());
    std::copy(qd.begin(), qd.end(), v.begin()+intraj.dimension);
    std::copy(qdd.begin(), qdd.end(), v.begin()+2*intraj.dimension);
    pouttraj->Insert(0, v);
    FOREACH(itchunk, intraj.chunkslist) {
        itchunk->Eval(0, q);
        itchunk->Evald(0, qd);
        itchunk->Evaldd(0, qdd);
        std::copy(q.begin(), q.end(), v.begin());
        std::copy(qd.begin(), qd.end(), v.begin()+intraj.dimension);
        std::copy(qdd.begin(), qdd.end(), v.begin()+2*intraj.dimension);
        v.at(3*intraj.dimension) = itchunk->duration;
        pouttraj->Insert(pouttraj->GetNumWaypoints(), v);
    }
//    dReal prevtime = 0;
//    FOREACH(ittime, intraj.chunkcumulateddurationslist) {
//        dReal deltatime = *ittime - prevtime;
//        intraj.Eval(*ittime, q);
//        intraj.Evald(*ittime, qd);
//        intraj.Evaldd(*ittime, qdd);
//        std::copy(q.begin(), q.end(), v.begin());
//        std::copy(qd.begin(), qd.end(), v.begin()+intraj.dimension);
//        std::copy(qdd.begin(), qdd.end(), v.begin()+2*intraj.dimension);
//        v.at(3*intraj.dimension) = deltatime;
//        
//        pouttraj->Insert(pouttraj->GetNumWaypoints(), v);
//        prevtime = *ittime;
//    }
}

TorqueLimitsRave2::TorqueLimitsRave2(RobotBasePtr probot, OpenRAVE::TrajectoryBaseConstPtr ptraj, dReal discrtimestep)
{
    EnvironmentMutex::scoped_lock lock(probot->GetEnv()->GetMutex()); // lock environment
    _probot = probot;
    RobotBase::RobotStateSaver robotsaver(probot, KinBody::Save_LinkTransformation|KinBody::Save_LinkVelocities);
    OPENRAVE_ASSERT_OP((int)probot->GetActiveDOFIndices().size(),==,probot->GetActiveDOF()); // don't allow affine dofs
    this->discrtimestep = discrtimestep;
    ConvertToTOPPTrajectory(ptraj, probot->GetActiveConfigurationSpecification(), trajectory);
    int ndof = trajectory.dimension;
    probot->GetActiveDOFVelocityLimits(vmax);
    hasvelocitylimits = false;
    FOREACH(itv, vmax) {
        if( std::abs(*itv) > TINY ) {
            hasvelocitylimits = true;
            break;
        }
    }

    std::vector<dReal> vgearratios(probot->GetActiveDOF());
    taumin.resize(probot->GetActiveDOF());
    taumax.resize(probot->GetActiveDOF());
    for(int i = 0; i < probot->GetActiveDOF(); ++i) {
        KinBody::JointPtr pjoint = probot->GetJointFromDOFIndex(probot->GetActiveDOFIndices()[i]);
        int iaxis = probot->GetActiveDOFIndices()[i] - pjoint->GetDOFIndex();
        dReal maxtorque = pjoint->GetMaxTorque(iaxis);
        dReal maxinertia = pjoint->GetMaxInertia(iaxis);
        // ElectricMotorActuatorInfo is a new spec in openrave
        ElectricMotorActuatorInfoPtr infoElectricMotor = pjoint->GetInfo()._infoElectricMotor;
        if( !!infoElectricMotor ) {
            dReal gear_ratio = infoElectricMotor->gear_ratio;
            //maxtorque =
        }
        else {
            std::cout << "could not find electric motor definition for joint " << pjoint->GetName() << std::endl;
        }
        taumax[i] = maxtorque;
        taumin[i] = -taumax[i];
    }

    // Define the avect, bvect, cvect
    int ndiscrsteps = int((trajectory.duration+1e-10)/discrtimestep)+1;
    std::vector<dReal> q(ndof), qd(ndof), qdd(ndof), vfullvalues(probot->GetDOF()), torquesimple;
    probot->GetDOFValues(vfullvalues);
    boost::array< std::vector< dReal >, 3 > torquecomponents;
    avect.resize(ndiscrsteps);
    bvect.resize(ndiscrsteps);
    cvect.resize(ndiscrsteps);
    for(int i = 0; i<ndiscrsteps; i++) {
        dReal s = i*discrtimestep;
        trajectory.Eval(s,q);
        trajectory.Evald(s,qd);
        trajectory.Evaldd(s,qdd);
        probot->SetActiveDOFValues(q,KinBody::CLA_Nothing);
        probot->SetActiveDOFVelocities(qd,KinBody::CLA_Nothing);
        for(int idof = 0; idof < ndof; ++idof) {
            vfullvalues[probot->GetActiveDOFIndices()[idof]] = qd[idof];
        }
        probot->ComputeInverseDynamics(torquesimple,vfullvalues);
        for(int idof = 0; idof < ndof; ++idof) {
            vfullvalues[probot->GetActiveDOFIndices()[idof]] = qdd[idof];
        }
        probot->ComputeInverseDynamics(torquecomponents,vfullvalues);
        avect[i].resize(ndof);
        bvect[i].resize(ndof);
        cvect[i].resize(ndof);
        for(int idof = 0; idof < ndof; ++idof) {
            avect[i][idof] = torquesimple[idof] - torquecomponents[1][idof] - torquecomponents[2][idof];
            bvect[i][idof] = torquecomponents[0][idof] + torquecomponents[1][idof];
            cvect[i][idof] = torquecomponents[2][idof];
        }
    }
}

void TorqueLimitsRave2::InterpolateDynamics(dReal s, std::vector<dReal>& a, std::vector<dReal>& b, std::vector<dReal>& c){
    a.resize(trajectory.dimension);
    b.resize(trajectory.dimension);
    c.resize(trajectory.dimension);
    assert(s >= -TINY && s <= trajectory.duration + TINY);
    if(s < 0)
        s = 0;
    if(s >= trajectory.duration - TINY) {
        int n = ndiscrsteps - 1;
        for(int i = 0; i < trajectory.dimension; i++) {
            a[i] = avect[n][i];
            b[i] = bvect[n][i];
            c[i] = cvect[n][i];
        }
        return;
    }
    int n = int(s / discrtimestep);
    dReal coef = (s - n * discrtimestep) / discrtimestep;
    for(int i = 0; i < trajectory.dimension; i++) {
        a[i] = (1-coef)*avect[n][i] + coef*avect[n+1][i];
        b[i] = (1-coef)*bvect[n][i] + coef*bvect[n+1][i];
        c[i] = (1-coef)*cvect[n][i] + coef*cvect[n+1][i];
    }
}

void TorqueLimitsRave2::ComputeSlopeDynamicSingularity(dReal s, dReal sd, std::vector<dReal>& slopesvector) {
    dReal delta = TINY2, s2, ap, bp, cp, slope;
    std::vector<dReal> a, b, c, a2, b2, c2;
    if(s>delta) {
        delta = -delta;
    }
    s2 = s + delta;
    InterpolateDynamics(s,a,b,c);
    InterpolateDynamics(s2,a2,b2,c2);

    slopesvector.resize(0);
    for(int i=0; i<trajectory.dimension; i++) {
        ap = (a2[i]-a[i])/delta;
        bp = (b2[i]-b[i])/delta;
        cp = (c2[i]-c[i])/delta;
        slope = (-bp*sd*sd-cp)/((2*b[i]+ap)*sd);
        slopesvector.push_back(slope);
    }
}

std::pair<dReal,dReal> TorqueLimitsRave2::SddLimits(dReal s, dReal sd){
    dReal dtsq = integrationtimestep;
    dtsq = dtsq*dtsq;
    dReal safetybound = discrtimestep/dtsq;
    dReal alpha = -safetybound;
    dReal beta = safetybound;
    dReal sdsq = sd*sd;
    std::vector<dReal> a, b, c;

    dReal taumin_i, taumax_i, alpha_i, beta_i;
    InterpolateDynamics(s,a,b,c);

    for(int i=0; i<trajectory.dimension; i++) {
        if(std::abs(a[i])<TINY) {
            continue;
        }
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
    return std::make_pair(alpha,beta);
}


dReal TorqueLimitsRave2::SdLimitBobrowInit(dReal s){
    std::pair<dReal,dReal> sddlimits = TorqueLimitsRave2::SddLimits(s,0);
    if(sddlimits.first > sddlimits.second) {
        return 0;
    }
    std::vector<dReal> tau_alpha(trajectory.dimension), tau_beta(trajectory.dimension);
    std::vector<dReal> a, b, c;
    InterpolateDynamics(s,a,b,c);

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
            if(std::abs(denum) > TINY) {
                r = num/denum;
                if(r>=0) {
                    sdmin = std::min(sdmin,sqrt(r));
                }
            }
            num = a[m]*(tau_alpha[k]-c[k])-a[k]*(tau_beta[m]-c[m]);
            denum = -denum;
            if(std::abs(denum) > TINY) {
                r = num/denum;
                if(r>=0) {
                    sdmin = std::min(sdmin,sqrt(r));
                }
            }
        }
    }
    return sdmin;
}

void TorqueLimitsRave2::FindSingularSwitchPoints(){
    if(ndiscrsteps<3) {
        return;
    }
    int i = 0;
    std::vector<dReal> a,aprev,b,c;

    InterpolateDynamics(discrsvect[i],aprev,b,c);

    for(int i=1; i<ndiscrsteps-1; i++) {
        InterpolateDynamics(discrsvect[i],a,b,c);
        dReal minsd = mvcbobrow[i];
        bool found = false;
        for(int j=0; j<trajectory.dimension; j++) {
            if(a[j]*aprev[j]<0) {
                dReal r = (taumin[j]-c[j])/b[j];
                if(r<0) {
                    r = (taumax[j]-c[j])/b[j];
                }
                if(r>0) {
                    found = true;
                    minsd = std::min(minsd,sqrt(r));
                }
            }
        }
        if(found) {
            //std::cout << discrsvect[i] << "," << minsd << "\n";
            AddSwitchPoint(i,SP_SINGULAR,minsd);
        }
        aprev = a;
    }
}

} // end namespace TOPP

#endif
