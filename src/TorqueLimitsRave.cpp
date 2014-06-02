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

#include <deque>
#include <fstream>

using namespace OpenRAVE;

namespace TOPP {

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
    else if( itposgroup->interpolation == "sextic" ) {
        degree = 6;
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
        if( deltatime <= TINY ) {
            // just skip since too small and will cause more problems with other epsilons later
            if( deltatime <= 0 ) {
            }
            continue;
        }
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
        if( itchunk->duration > 0 ) {
            itchunk->Eval(0, q);
            itchunk->Evald(0, qd);
            itchunk->Evaldd(0, qdd);
            std::copy(q.begin(), q.end(), v.begin());
            std::copy(qd.begin(), qd.end(), v.begin()+intraj.dimension);
            std::copy(qdd.begin(), qdd.end(), v.begin()+2*intraj.dimension);
            v.at(3*intraj.dimension) = itchunk->duration;
            pouttraj->Insert(pouttraj->GetNumWaypoints(), v);
        }
        else if( itchunk->duration < 0 ) {
            std::cerr << "chuck duration is negative! " << itchunk->duration << std::endl;
        }
    }
}

//OpenRAVE::TrajectoryBaseConstPtr intraj, 
bool ExtractOpenRAVETrajectoryFromProfiles(const Constraints& constraints, dReal smax, OpenRAVE::TrajectoryBasePtr outtraj)
{
    if (constraints.resprofileslist.size() < 1) {
        return false;
    }

    if( smax == 0 ) {
        smax = constraints.trajectory.duration;//intraj->GetDuration();
    }
    
    //std::deque<dReal> vsampledpoints; // s, sd, sdd, deltatime
    std::vector<dReal> vsampledpoints;
    vsampledpoints.reserve(4*10000); // no idea how many points, but guessing a lot if integrationstep is 0.001 and duration is 10s
    
    ProfileSample sample = FindLowestProfileFast(0, 1e30, constraints.resprofileslist);
    if( sample.itprofile == constraints.resprofileslist.end() ) {
        RAVELOG_ERROR("failed to find first profile\n");
        return false;
    }
    vsampledpoints.push_back(sample.s);
    vsampledpoints.push_back(sample.sd);
    vsampledpoints.push_back(sample.sdd);
    vsampledpoints.push_back(0);
    
    bool bsuccess = false;
    while(!bsuccess) {
        const Profile& profile = *sample.itprofile;
        
        bool badded = false;
        dReal sprev = sample.s, sdprev = sample.sd, sddprev = sample.sdd;
        dReal tprev = sample.t;
        int sindex = sample.sindex+1;
        ProfileSample checksample;
        checksample.itprofile = constraints.resprofileslist.end();
        while(sindex < (int)profile.svect.size() && sprev < smax-TINY) {
            dReal s = profile.svect.at(sindex);
            dReal sd = profile.sdvect.at(sindex);
            dReal sdd = profile.sddvect.at(sindex);
            
            // check if there's a lower profile at s
            checksample = FindLowestProfileFast(s, sd-TINY2, constraints.resprofileslist);
            if( checksample.itprofile != constraints.resprofileslist.end() ) {
                if( sample.itprofile == checksample.itprofile ) {
                    RAVELOG_ERROR("got sample profile, unexpected!\n");
                    return false;
                }
                // found new profile, so need to use it instead.
                // have to merge with old profile by finding an intersection point.
                // s(sd) = sprev + (sd*sd - sdprev*sdprev)/(2*sddprev) => aprev * sd*sd + cprev
                if( fabs(sddprev) <= TINY ) {
                    RAVELOG_ERROR_FORMAT("sprev=%.%15e, sdprev=%.15e, sddprev=%.15e, sdd is close to 0, don't know that to do", sprev%sdprev%sddprev);
                    return false;
                }
                else if( fabs(checksample.sdd) <= TINY ) {
                    RAVELOG_ERROR("check sddprev is close to 0, don't know that to do\n");
                    return false;
                }
                else {
                    dReal aprev = 1/(2*sddprev), cprev = (sprev - sdprev*sdprev/(2*sddprev));
                    dReal anext = 1/(2*checksample.sdd), cnext = (checksample.s - checksample.sd*checksample.sd/(2*checksample.sdd));
                    dReal ad = aprev - anext;
                    if( fabs(ad) <= TINY ) {
                        RAVELOG_ERROR("two profiles have same accel, don't know that to do\n");
                        return false;
                    }
                    else {
                        dReal sdintersect2 = (cnext-cprev)/ad;
                        if( sdintersect2 < 0 ) {
                            RAVELOG_ERROR("two profiles never intersect!\n");
                            return false;
                        }
                        dReal sdintersect = sqrt(sdintersect2);
                        dReal sintersect = aprev*sdintersect2 + cprev;
                        dReal sddintersect = sddprev;
                        dReal tintersect = (sdintersect-sdprev)/sddprev;
                        BOOST_ASSERT(sprev < sintersect);
                        vsampledpoints.push_back(sintersect);
                        vsampledpoints.push_back(checksample.sd);  // needs to be checksample since interpolation goes forward
                        vsampledpoints.push_back(checksample.sdd); // needs to be checksample since interpolation goes forward
                        vsampledpoints.push_back(tintersect);
                        dReal t2 = (sdintersect-checksample.sd)/checksample.sdd;
                        if( t2 > TINY ) {
                            RAVELOG_ERROR_FORMAT("unexpted got time %f", t2);
                            return false;
                        }
                        checksample.t += t2; // go back in time
//                        if( t2 > TINY ) {
//                            BOOST_ASSERT(sintersect < checksample.s);
//                            vsampledpoints.push_back(checksample.s);
//                            vsampledpoints.push_back(checksample.sd);
//                            vsampledpoints.push_back(checksample.sdd);
//                            vsampledpoints.push_back(t2);
//                        }
                        checksample.s = sintersect;
                        sprev = checksample.s;
                        sdprev = checksample.sd;
                        sddprev = checksample.sdd;
                        badded = true;
                        break;
                    }
                }
                RAVELOG_ERROR("unexpected point\n");
                return false;
            }

            BOOST_ASSERT(sprev < s);
            vsampledpoints.push_back(s);
            vsampledpoints.push_back(sd);
            vsampledpoints.push_back(sdd);
            vsampledpoints.push_back(profile.integrationtimestep-tprev);
            badded = true;
            
            // increase the point and try again
            tprev = 0;
            sprev = s;
            sdprev = sd;
            sddprev = sdd;
            ++sindex;
        }
        
        if( sprev >= smax-TINY ) {
            bsuccess = true;
            break;
        }
        
        if( !badded ) {
            RAVELOG_ERROR("failed to add a new sample, something is wrong\n");
            return false;
        }
        
        if( checksample.itprofile == constraints.resprofileslist.end() ) {
            dReal sstart = std::min(smax, sprev-TINY);
            dReal sdthresh = 0.01;
            //sample.itprofile = constraints.resprofileslist.end();
            // got to end of ramp and another ramp was not found so search for the next ramp.
            std::list<Profile>::const_iterator itprofilemin = constraints.resprofileslist.end();
            for(std::list<Profile>::const_iterator itprofile = constraints.resprofileslist.begin(); itprofile != constraints.resprofileslist.end(); ++itprofile) {
                if( itprofile->svect.size() >= 2 && itprofile->svect.at(0) >= sstart && fabs(itprofile->sdvect.at(0)-sdprev) < sdthresh ) {
                    if( itprofilemin == constraints.resprofileslist.end() || itprofile->svect.at(0) < itprofilemin->svect.at(0) ) {
                        itprofilemin = itprofile;
                    }
                }
            }
            if( itprofilemin == constraints.resprofileslist.end() ) {
                RAVELOG_ERROR_FORMAT("failed to find next profile s=%.15e", sstart);
                return false;
            }
            
            // in sample.itprofile, take the second to last point since that should have sdd != 0
            vsampledpoints.resize(vsampledpoints.size()-4);    
            if( vsampledpoints.size() >= 4 ) {
                sprev = vsampledpoints[vsampledpoints.size()-4];
                sdprev = vsampledpoints[vsampledpoints.size()-3];
                sddprev = vsampledpoints[vsampledpoints.size()-2];
            }
            else {
                sprev = sample.itprofile->svect.at(sample.itprofile->svect.size()-2);
                sdprev = sample.itprofile->sdvect.at(sample.itprofile->svect.size()-2);
                sddprev = sample.itprofile->sddvect.at(sample.itprofile->svect.size()-2);
            }
            
            dReal snext = itprofilemin->svect.at(0);
            dReal sdnext = itprofilemin->sdvect.at(0);
            dReal sddnext = itprofilemin->sddvect.at(0);
            // s(sd) = sprev + (sd*sd - sdprev*sdprev)/(2*sddprev) => aprev * sd*sd + cprev
            if( fabs(sddprev) <= TINY ) {
                RAVELOG_ERROR("sddprev is close to 0, don't know that to do\n");
                return false;
            }
            else if( fabs(sddnext) <= TINY ) {
                RAVELOG_ERROR("check sddprev is close to 0, don't know that to do\n");
                return false;
            }
            else {
                dReal aprev = 1/(2*sddprev), cprev = (sprev - sdprev*sdprev/(2*sddprev));
                dReal anext = 1/(2*sddnext), cnext = (snext - sdnext*sdnext/(2*sddnext));
                dReal ad = aprev - anext;
                if( fabs(ad) <= TINY ) {
                    RAVELOG_ERROR("two profiles have same accel, don't know that to do\n");
                    return false;
                }
                else {
                    dReal sdintersect2 = (cnext-cprev)/ad;
                    if( sdintersect2 < 0 ) {
                        RAVELOG_ERROR("two profiles never intersect!\n");
                        return false;
                    }
                    dReal sdintersect = sqrt(sdintersect2);
                    dReal sintersect = aprev*sdintersect2 + cprev;
                    dReal sddintersect = sddprev;
                    dReal tintersect = (sdintersect-sdprev)/sddprev;
//                    if(vsampledpoints[vsampledpoints.size()-4] > sintersect-TINY) {
//                        // intersection happens between previous arc so remove last point
//                        vsampledpoints.resize(vsampledpoints.size()-4);
//                    }
//                    else {
//                        BOOST_ASSERT(tintersect >= vsampledpoints.back());
//                        tintersect -= vsampledpoints.back(); // have to subtract the last time
//                    }
                    
                    vsampledpoints.push_back(sintersect);
                    vsampledpoints.push_back(sdnext); // speed needs to be sddnext since interpolating forward sddintersect);
                    vsampledpoints.push_back(sddnext); // acceleration needs to be sddnext since interpolating forward sddintersect);
                    vsampledpoints.push_back(tintersect);

                    dReal t2 = (sdintersect-sdnext)/sddnext;
                    checksample.itprofile = itprofilemin;
                    checksample.sindex = 0;
                    checksample.s = sintersect;
                    checksample.sd = sdnext;
                    checksample.sdd = sddnext;
                    checksample.t = t2;
                }
            }
        }
        
        sample = checksample;
    }

//    RAVELOG_INFO_FORMAT("success in extracting profiles (%d)!", vsampledpoints.size());
//    std::ofstream f("points.txt");
//    f << std::setprecision(std::numeric_limits<OpenRAVE::dReal>::digits10+1);
//    for(size_t i =0; i < vsampledpoints.size(); i += 4 ) {
//        f << vsampledpoints[i] << " " << vsampledpoints[i+1] << " " << vsampledpoints[i+2] << " " << vsampledpoints[i+3] << std::endl;
//    }
    return true;
}

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
    int ndiscrsteps = int((trajectory.duration+TINY)/discrtimestep)+1;
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

TorqueLimitsRave2::TorqueLimitsRave2(RobotBasePtr probot, OpenRAVE::TrajectoryBaseConstPtr ptraj, dReal _discrtimestep)
{
    EnvironmentMutex::scoped_lock lock(probot->GetEnv()->GetMutex()); // lock environment
    _probot = probot;
    RobotBase::RobotStateSaver robotsaver(probot, KinBody::Save_LinkTransformation|KinBody::Save_LinkVelocities);
    OPENRAVE_ASSERT_OP((int)probot->GetActiveDOFIndices().size(),==,probot->GetActiveDOF()); // don't allow affine dofs
    discrtimestep = _discrtimestep;
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
    int ndiscrsteps = int((trajectory.duration+TINY)/discrtimestep)+1;
    if( ndiscrsteps > 1 ) {
        discrtimestep = trajectory.duration/(ndiscrsteps-1);
    }
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
    BOOST_ASSERT(s >= -TINY && s <= trajectory.duration + TINY);
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
    dReal coef = s - n * discrtimestep;
    if( std::abs(coef) <= TINY ) {
        a = avect[n];
        b = bvect[n];
        c = cvect[n];
    }
    else {
        BOOST_ASSERT(n+1 < (int)avect.size());
        coef /= discrtimestep;
        for(int i = 0; i < trajectory.dimension; i++) {
            a[i] = (1-coef)*avect[n][i] + coef*avect[n+1][i];
            b[i] = (1-coef)*bvect[n][i] + coef*bvect[n+1][i];
            c[i] = (1-coef)*cvect[n][i] + coef*cvect[n+1][i];
        }
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

    // at the limit, the alpha and beta curves of two dofs should meet. compute those points
    dReal sdmin = INF;
    for(int k=0; k<trajectory.dimension; k++) {
        for(int m=k+1; m<trajectory.dimension; m++) {
            dReal num, denum, r;
            denum = a[k]*b[m]-a[m]*b[k];
            if(std::abs(denum) > TINY) {
                num = a[k]*(tau_alpha[m]-c[m])-a[m]*(tau_beta[k]-c[k]);
                r = num/denum;
                if(r>=0) {
                    sdmin = std::min(sdmin,sqrt(r));
                }
                
                num = a[m]*(tau_alpha[k]-c[k])-a[k]*(tau_beta[m]-c[m]);
                denum = -denum;
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
            //std::cout << "singular switch: " << discrsvect[i] << " minsid=" << minsd << std::endl;
            AddSwitchPoint(i,SP_SINGULAR,minsd);
        }
        aprev.swap(a);
    }
}

void TorqueLimitsRave2::FindDiscontinuousSwitchPoints() {
    if(ndiscrsteps<3)
        return;
    int i = 0;
    dReal sd, sdn, sdnn;
    sd = SdLimitBobrow(discrsvect[i]);
    sdn = SdLimitBobrow(discrsvect[i+1]);
    // also look for the start of the chucks for the trajectory
    std::list<dReal>::const_iterator itchuckstart = trajectory.chunkcumulateddurationslist.begin();
    int nLastAddedSwitchIndex = -1;
    for(int i=0; i<ndiscrsteps-3; i++) {
        sdnn = SdLimitBobrow(discrsvect[i+2]);
        if(std::abs(sdnn-sdn)>100*std::abs(sdn-sd)) {
            if(sdn<sdnn) {
                AddSwitchPoint(i+1,SP_DISCONTINUOUS);
                nLastAddedSwitchIndex = i+1;
            }
            else{
                AddSwitchPoint(i+2,SP_DISCONTINUOUS);
                nLastAddedSwitchIndex = i+2;
            }
        }
        if( trajectory.degree <= 3 ) {
            // if the trajectory degree is <= 3, then the accelerations will not be differentiable at the trajectory chunk edges.
            // therefore add those discontinuity points.
            // perhaps there's a better way to compute this, but the above threshold doesn't catch it.
            if( itchuckstart != trajectory.chunkcumulateddurationslist.end() && *itchuckstart < discrsvect[i+2] ) {
                if( nLastAddedSwitchIndex < i+1 ) {
                    AddSwitchPoint(i+1,SP_DISCONTINUOUS);
                    nLastAddedSwitchIndex = i+1;
                }
                ++itchuckstart;
            }
        }
        sd = sdn;
        sdn = sdnn;
    }
}

} // end namespace TOPP

#endif
