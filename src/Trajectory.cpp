// -*- coding: utf-8 -*-
// Copyright (C) 2013 Quang-Cuong Pham <cuong.pham@normalesup.org>
//
// This file is part of the Time-Optimal Path Parameterization (TOPP) library.
// TOPP is free software: you can redistribute it and/or modify
// it under the terms of the GNU Lesser General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// at your option, any later version.
//
// This program is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU Lesser General Public License for more details.
//
// You should have received a copy of the GNU Lesser General Public License
// along with this program.  If not, see <http://www.gnu.org/licenses/>.


#include "TOPP.h"


namespace TOPP {


///////////////////////// Polynomial ///////////////////////


void Polynomial::InitFromCoefficientsVector(const std::vector<dReal>&coefficientsvector0) {
    coefficientsvector = coefficientsvector0;
    degree = coefficientsvector.size() - 1;
    // Construct first- and second-order derivative polynomials
    coefficientsvectord.resize(0);
    coefficientsvectordd.resize(0);
    for(int i = 1; i <= degree; i++) {
        coefficientsvectord.push_back(i * coefficientsvector[i]);
    }
    for(int i = 1; i <= degree - 1; i++) {
        coefficientsvectordd.push_back(i * coefficientsvectord[i]);
    }
}


Polynomial::Polynomial(const std::vector<dReal>& coefficientsvector0) {
    InitFromCoefficientsVector(coefficientsvector0);
}


Polynomial::Polynomial(std::string& s) {
    VectorFromString(s,coefficientsvector);
    InitFromCoefficientsVector(coefficientsvector);

}


// Evaluate polynomials using Horner's method
dReal Polynomial::Eval(dReal s) const {
    dReal res = 0;
    for(int i = degree; i >= 0; i--)
        res = res * s + coefficientsvector[i];
    return res;
}


dReal Polynomial::Evald(dReal s) const {
    dReal res = 0;
    for(int i = degree - 1; i >= 0; i--)
        res = res * s + coefficientsvectord[i];
    return res;
}


dReal Polynomial::Evaldd(dReal s) const {
    dReal res = 0;
    for(int i = degree - 2; i>=0; i--)
        res = res*s + coefficientsvectordd[i];
    return res;
}

dReal Polynomial::Evalddd(dReal s) const
{
    dReal res = 0;
    for(int i = degree - 2; i>=1; i--) {
        res = res*s + coefficientsvectordd[i]*i;
    }
    return res;
}

dReal Polynomial::Evaldddd(dReal s) const
{
    dReal res = 0;
    for(int i = degree - 2; i>=2; i--) {
        res = res*s + coefficientsvectordd[i]*i*(i-1);
    }
    return res;
}

void Polynomial::Write(std::stringstream& ss) {
    ss << std::setprecision(17) << coefficientsvector[0];
    for(int i = 1; i <= degree; i++)
        ss << " " << std::setprecision(17) << coefficientsvector[i];
}


///////////////////////// Chunk ///////////////////////


Chunk::Chunk(dReal duration0, const std::vector<Polynomial>& polynomialsvector0) {
    polynomialsvector = polynomialsvector0;
    dimension = polynomialsvector.size();
    duration = duration0;
    BOOST_ASSERT(dimension > 0);
    degree = polynomialsvector[0].degree;
    for(int i = 1; i < dimension; i++)
        if (polynomialsvector[i].degree > degree)
            degree = polynomialsvector[i].degree;
    // All polynomials must have the same degree
    for(int i = 1; i < dimension; i++)
        BOOST_ASSERT(degree == polynomialsvector[i].degree);
}


void Chunk::Eval(dReal s, std::vector<dReal>&q) const {
    BOOST_ASSERT(s >= -TINY2);
    BOOST_ASSERT(s <= duration+TINY2);
    for(int i = 0; i < dimension; i++)
        q[i] = polynomialsvector[i].Eval(s);
}


void Chunk::Evald(dReal s, std::vector<dReal>&qd) const {
    BOOST_ASSERT(s >= -TINY2);
    BOOST_ASSERT(s <= duration+TINY2);
    for(int i = 0; i < dimension; i++)
        qd[i] = polynomialsvector[i].Evald(s);
}


void Chunk::Evaldd(dReal s, std::vector<dReal>&qdd) const {
    BOOST_ASSERT(s >= -TINY2);
    BOOST_ASSERT(s <= duration+TINY2);
    for(int i = 0; i < dimension; i++)
        qdd[i] = polynomialsvector[i].Evaldd(s);
}

void Chunk::Evalddd(dReal s, std::vector<dReal>&qddd) const {
    BOOST_ASSERT(s >= -TINY2);
    BOOST_ASSERT(s <= duration+TINY2);
    for(int i = 0; i < dimension; i++)
        qddd[i] = polynomialsvector[i].Evalddd(s);
}

void Chunk::Evaldddd(dReal s, std::vector<dReal>&qddd) const {
    BOOST_ASSERT(s >= -TINY2);
    BOOST_ASSERT(s <= duration+TINY2);
    for(int i = 0; i < dimension; i++)
        qddd[i] = polynomialsvector[i].Evaldddd(s);
}

void Chunk::Write(std::stringstream& ss) {
    ss << std::setprecision(17) <<  duration << "\n";
    ss << dimension << "\n";
    for(int i = 0; i < dimension; i++) {
        polynomialsvector[i].Write(ss);
        ss << "\n";
    }
}


/////////////// Trajectory ///////////////////////


void Trajectory::InitFromChunksList(const std::list<Chunk>&chunkslist0) {
    chunkslist = chunkslist0;
    BOOST_ASSERT(chunkslist.size()>0);
    dimension = chunkslist.front().dimension;
    degree = chunkslist.front().degree;

    duration = 0;
    chunkdurationslist.resize(0);
    chunkcumulateddurationslist.resize(0);
    std::list<Chunk>::iterator itchunk = chunkslist.begin();
    while(itchunk != chunkslist.end()) {
        dReal chunkduration = itchunk->duration;
        if(chunkduration > TINY) {
            chunkdurationslist.push_back(chunkduration);
            chunkcumulateddurationslist.push_back(duration);  // Cumulated durations list starts with a 0
            itchunk->sbegin = duration;
            duration += chunkduration;
            itchunk->send = duration;
        }
        itchunk++;
    }
    chunkcumulateddurationslist.push_back(duration);
}


Trajectory::Trajectory(const std::list<Chunk>& chunkslist0) {
    InitFromChunksList(chunkslist0);
}


Trajectory::Trajectory(const std::string& trajectorystring) {
    std::string buff;
    std::istringstream iss(trajectorystring);
    int dimension;
    dReal duration;
    std::vector<Polynomial> polynomialsvector;
    std::list<Chunk> chunkslist0;
    while(iss.good()) {
        getline(iss, buff, '\n');
        duration = atof(buff.c_str());
        getline(iss, buff, '\n');
        dimension = atoi(buff.c_str());
        polynomialsvector.resize(0);
        for(int i = 0; i < dimension; i++) {
            getline(iss, buff, '\n');
            polynomialsvector.push_back(Polynomial(buff));
        }
        if(duration>TINY) {
            chunkslist0.push_back(Chunk(duration,polynomialsvector));
        }
    }
    InitFromChunksList(chunkslist0);
}


void Trajectory::FindChunkIndex(dReal s, int& index, dReal& remainder) const {

    std::list<dReal>::const_iterator it = chunkcumulateddurationslist.begin();
    if(s <= TINY) {
        index = 0;
        remainder = 0;
        return;
    }
    if(s >= chunkcumulateddurationslist.back()) {
        index = int(chunkslist.size())-1;
        remainder = chunkslist.back().duration;
        return;
    }
    index = 0;
    while(it != chunkcumulateddurationslist.end() && s > *it) {
        index++;
        it++;
    }
    index--;
    BOOST_ASSERT(index<=int(chunkslist.size())-1);
    it--;
    remainder = s-*it;
}


void Trajectory::Eval(dReal s, std::vector<dReal>&q) const {
    BOOST_ASSERT(s >= -TINY);
    BOOST_ASSERT(s <= duration+TINY);
    BOOST_ASSERT(dimension == int(q.size()));
    int index;
    dReal remainder;
    FindChunkIndex(s,index,remainder);
    std::list<Chunk>::const_iterator itchunk = chunkslist.begin();
    advance(itchunk,index);
    itchunk->Eval(remainder,q);
}


void Trajectory::Evald(dReal s, std::vector<dReal>&qd) const {
    BOOST_ASSERT(s >= 0-TINY);
    BOOST_ASSERT(s <= duration+TINY);
    BOOST_ASSERT(dimension == int(qd.size()));
    int index;
    dReal remainder;
    FindChunkIndex(s,index,remainder);
    std::list<Chunk>::const_iterator itchunk = chunkslist.begin();
    advance(itchunk,index);
    itchunk->Evald(remainder,qd);
}


void Trajectory::Evaldd(dReal s, std::vector<dReal>&qdd) const {
    BOOST_ASSERT(s >= -TINY);
    BOOST_ASSERT(s <= duration+TINY);
    int index;
    dReal remainder;
    FindChunkIndex(s,index,remainder);
    std::list<Chunk>::const_iterator itchunk = chunkslist.begin();
    advance(itchunk,index);
    itchunk->Evaldd(remainder,qdd);
}

void Trajectory::ComputeChunk(dReal t0, dReal tnext, dReal s, dReal sd, dReal sdd,
                              const Chunk& currentchunk, Chunk& newchunk) {

    int n = currentchunk.degree;
    int ndof = currentchunk.dimension;
    std::vector<dReal> a, rescoeffs;
    // coeffsvect[i] contains coefficients list for s^i
    std::vector<std::vector<dReal> > coeffsvects;
    std::vector<Polynomial> polynomialsvector;

    // currentchunk : sum_{i = 0}^{n} (u_i)(s^i)
    // profile : s + sd*t + 0.5*sdd*t^2
    // new chunk : sum_{j = 0}^{2n} (v_j)(t^j)

    a.resize(0);
    a.push_back(s + sd*t0 + 0.5*sdd*t0*t0);
    a.push_back(sd + sdd*t0);
    a.push_back(0.5*sdd);
    coeffsvects.push_back(a);

    std::vector<dReal> tmpvect;

    if (n >= 1) {
        for (int i = 1; i < n; i++) {
            tmpvect.resize(0);
            tmpvect.resize(2*i + 3, 0);
            for (int j = 0; j < 2*i + 1; j++) {
                tmpvect[j] += a[0]*coeffsvects[i - 1][j];
                tmpvect[j + 1] += a[1]*coeffsvects[i - 1][j];
                tmpvect[j + 2] += a[2]*coeffsvects[i - 1][j];
            }
            coeffsvects.push_back(tmpvect);
        }
    }

    for (int i = 0; i < ndof; i++) {
        rescoeffs.resize(0);
        rescoeffs.push_back(currentchunk.polynomialsvector[i].coefficientsvector[0]);
        rescoeffs.resize(2*n + 1, 0);
        for (int k = 1; k <= n; k++) {
            dReal u = currentchunk.polynomialsvector[i].coefficientsvector[k];
            int l = 2*k + 1;
            for (int j = 0; j < l; j++) {
                rescoeffs[j] += u*coeffsvects[k - 1][j];
            }
        }
        polynomialsvector.push_back(Polynomial(rescoeffs));
    }
    newchunk = Chunk(tnext - t0, polynomialsvector);
}


void Trajectory::SPieceToChunks(dReal s, dReal sd, dReal sdd, dReal T, int&
                                currentchunkindex, dReal& processedcursor, std::list<Chunk>::iterator&
                                itcurrentchunk, std::list<Chunk>& chunkslist) {
    dReal t = 0, tnext;
    dReal snext = s + T*sd + 0.5*T*T*sdd;
    int chunkindex;
    dReal remainder;
    Chunk newchunk;

    FindChunkIndex(snext,chunkindex,remainder);

    // Process all chunks that have been overpassed
    while(currentchunkindex<chunkindex) {
        if(itcurrentchunk->duration-processedcursor>=TINY) {
            bool res = SolveQuadraticEquation(s-itcurrentchunk->send,sd,0.5*sdd,tnext,t,T);
            BOOST_ASSERT(res);
            ComputeChunk(t,tnext,s-itcurrentchunk->sbegin,sd,sdd,*itcurrentchunk,newchunk);
            chunkslist.push_back(newchunk);
            t = tnext;
        }
        currentchunkindex++;
        itcurrentchunk++;
        processedcursor = 0;
    }

    // Process current chunk
    bool res = SolveQuadraticEquation((s-itcurrentchunk->sbegin)-remainder,sd,0.5*sdd,tnext,t,T);
    BOOST_ASSERT(res);
    ComputeChunk(t,tnext,s-itcurrentchunk->sbegin,sd,sdd,*itcurrentchunk,newchunk);
    chunkslist.push_back(newchunk);
    processedcursor = remainder;
}


int Trajectory::Reparameterize(Constraints& constraints, Trajectory& restrajectory, dReal smax) {

    if (constraints.resprofileslist.size() < 1) {
        return -1;
    }

    dReal scur, sdcur, snext, sdnext, sdnext2, sdd;
    dReal dt = constraints.reparamtimestep;

    // Set the reparam timestep automatically if it is initially set to 0
    if(dt == 0) {
        if(smax == 0 && constraints.resduration>TINY) {
            dt = constraints.discrtimestep*constraints.resduration/duration;
        }
        else{
            dt = constraints.discrtimestep;
        }
    }

    if (smax == 0) {
        smax = duration;
    }

    dReal dtsq = dt*dt;
    dReal dtmod;
    Profile profile;
    dReal tres;
    std::list<Chunk> newchunkslist;

    std::list<Chunk>::iterator itcurrentchunk = chunkslist.begin();
    int currentchunkindex = 0;
    dReal processedcursor = 0;

    scur = 0;
    dReal t = 0;
    FindLowestProfile(scur, profile, tres, constraints.resprofileslist);
    sdcur = profile.Evald(tres);

    while(scur<smax) {
        sdd = profile.Evaldd(tres);
        sdnext = sdcur + dt*sdd;
        snext = scur + dt*sdcur + 0.5*dtsq*sdd;
        if(snext >= scur+TINY && snext<= smax && FindLowestProfile(snext,profile,tres,constraints.resprofileslist)) {
            sdnext2 = profile.Evald(tres);
            dtmod = dt;
            // If discrepancy between integrated sd and profile's sd then
            // follow profile's sd, which requires changing dt
            if(std::abs(sdnext-sdnext2)>TINY2) {
                dtmod = 2*(snext-scur)/(sdnext2+sdcur);
                sdd = (sdnext2-sdcur)/dtmod;
            }
            SPieceToChunks(scur, sdcur, sdd, dtmod, currentchunkindex,
                           processedcursor, itcurrentchunk, newchunkslist);
        }
        else
            break;
        t+= dtmod;
        scur = snext;
        sdcur = sdnext2;
    }

    if (newchunkslist.size() < 1) {
        return -1;
    }
    restrajectory = Trajectory(newchunkslist);
    return 1;
}


int Trajectory::Reparameterize2(Constraints& constraints, Trajectory& restrajectory,
				dReal smax) {
    std::cout << "Running Reparameterize2" << "\n";
    
    std::string message;
    bool _print = true;
    
    if (constraints.resprofileslist.size() == 0) {
	return -1;
    }

    Trajectory& intraj = constraints.trajectory;
    if (intraj.chunkslist.size() == 0) {
	return -1;
    }

    if (smax == 0) {
	smax = intraj.duration;
    }

    // vsampledpoints contains s, sd, sdd, deltatime for each reparam step
    std::vector<dReal> vsampledpoints;
    vsampledpoints.reserve(80000);

    ProfileSample sample = FindLowestProfileFast(0, 1e30, constraints.resprofileslist);
    if (sample.itprofile == constraints.resprofileslist.end()) {
	message = "[Reparameterize2] Failed to find the lowest profile at s = 0.";
	std::cout << message << std::endl;
	return -1;
    }
    vsampledpoints.push_back(sample.s);
    vsampledpoints.push_back(sample.sd);
    vsampledpoints.push_back(sample.sdd);
    vsampledpoints.push_back(0);

    // Now extracting s, sd, sdd, delta, from the lowest profile and each reparam step
    bool bsuccess = false;
    while (!bsuccess) {
	const Profile& profile = *sample.itprofile;
	
	dReal sprev = sample.s, sdprev = sample.sd, sddprev = sample.sdd;
	dReal tprev = sample.t;
	int sindex = sample.sindex + 1; // next index
	ProfileSample checksample, checksample2;
	
	// Reset checksample
	checksample.itprofile = constraints.resprofileslist.end();
	
	// Step along the current profile (sample)
	bool badded = false; // set to true is the stepping is successful
	while (sindex < (int)profile.svect.size() && sprev < smax - TINY) {
	    dReal s = profile.svect.at(sindex);
	    dReal sd = profile.sdvect.at(sindex);
	    dReal sdd = profile.sddvect.at(sindex);
	    
	    // Check if there is any lower profile at s
	    checksample = FindLowestProfileFast(s, sd - 10*TINY,
						constraints.resprofileslist);
	    if (checksample.itprofile != constraints.resprofileslist.end()) {
		message = str(boost::format("[Reparameterize2]"
					    " Found a new lower profile at s = %.15e, sd = %.15e")
			      %s%sd);


		if (sample.itprofile == checksample.itprofile) {
		    message = "[Reparameterize2] Checking lower profiles got sample profile";
		    std::cout << message << std::endl;
		    return -1;
		}
		
		// We have found a new lower profile.
		// There must be an intersection somewhere.
		dReal sintersect = 0, sdintersect = 0, sddintersect = 0, tintersect = 0;
		// busechecksample tells us whether or not to use the
		// newly found profile.
		bool busechecksample = true;
		
		if (fabs(sddprev) <= TINY) { // why checking this ?
		    if (_print) {
			message = str(boost::format("[Reparameterize2]"
						    " sprev = %.15e, sdprev = %.15e,"
						    " sddprev = %.15e; sdd is close to 0.")
				      %sprev%sdprev%sddprev);
			std::cout << message << std::endl;
		    }
		    return -1;
		}
		else {
		    if (fabs(checksample.sdd) <= TINY) {
			// The new profile has sdd = 0
			tintersect = (checksample.sd - sdprev)/sddprev;
			sintersect = sprev + tintersect*(sdprev + 0.5*tintersect*sddprev);
			sdintersect = checksample.sd; // Rosen put sdprev here. I think
						      // that is not correct.
			// sddintersect always needs to be
			// checksample.sdd because the interpolation
			// goes forward.
			sddintersect = checksample.sdd;
		    }
		    else {
			/* Consider the equation 
			     v^2 = u^2 + 2ad,          --(1)
			 where u, v are initial, final velocities, a
			 is the (constant) acceleration, and d is the
			 distance.

			 For profile, we substitute v = sd, u = sdprev, 
			 a = sddprev, and d = (s - sprev) into the above equation. 
			 We have 
			 
			     s(sd) = (1/2sddprev)sd^2 + (sprev - sdprev^2/2sddprev). 
			     
			 Note that s is a function of sd.

			 For checksample, we substitute v = csd, u = sd,
			 a = csdd, and d = (cs - s) into Equation (1) to get

			     s(sd) = (1/2csddprev)sd^2 + (cs - csd^2),
			     
			 where cX means checksample.X.

			 We will use the two equations of the form
			 
			     s(sd) = A*sd^2 + C

			 to find the intersection between profile and
			 checksample.
			 */
			dReal Aprofile = 1/(2*sddprev);
			dReal Cprofile = sprev - sdprev*sdprev/(2*sddprev);
			
			dReal Acheck = 1/(2*checksample.sdd);
			dReal Ccheck = checksample.s - 
			    checksample.sd*checksample.sd/(2*checksample.sdd);

			dReal Adiff = Aprofile - Acheck;
			if (fabs(Adiff) <= TINY) {
			    // I am not sure how we should handle the
			    // following two cases.
			    // Just follow what Rosen does.
			    if (fabs(Ccheck - Cprofile) <= TINY) {
				if (_print) {
				    message = "[Reparameterize2] Adiff = Cdiff = 0";
				    std::cout << message << std::endl;
				}
				// most lokely they intersect at the same place
				sintersect = checksample.s;
				sdintersect = checksample.sd;
				sddintersect = checksample.sdd;
				tintersect = (sdintersect - sdprev)/sddprev;
			    }
			    else {
				if (_print)  {
				    message = "[Reparameterize2] Adiff = 0 but Cdiff != 0";
				    std::cout << message << std::endl;
				}
			    }
			}
			else {
			    dReal sdintersect2 = (Ccheck - Cprofile)/Adiff;
			    // I don't know why Rosen checks if
			    // sdintersect2 >= s*s here.
			    if (sdintersect2 >= s*s) {
				sdintersect = sqrt(sdintersect2);
				sintersect = Aprofile*sdintersect2 + Cprofile;
				tintersect = (sdintersect - sdprev)/sddprev;
				// sddintersect always needs to be
				// checksample.sdd because the interpolation
				// goes forward.
				sddintersect = checksample.sdd;
			    }
			}
		    }
		    
		    if (tintersect > 0 && sintersect <= s) {
			if (sprev <= sintersect + TINY) {
			    // Now we obtain a valid intersection point.
			    vsampledpoints.push_back(sintersect);
			    vsampledpoints.push_back(sdintersect);
			    vsampledpoints.push_back(sddintersect);
			    vsampledpoints.push_back(tintersect);

			    // What is this for ?
			    dReal t2 = (sdintersect - checksample.sd)/checksample.sdd;
			    checksample.t += t2; // go back in time
			    checksample.s = sintersect;
			}
			else {
			    // sintersect < sprev
			    // what does this actually mean ?
			    if (_print) {
				message = "[Reparameterize2] sintersect < sprev";
				std::cout << message << std::endl;
			    }
			    busechecksample = false;
			}
		    }
		    else {
			// Here we encounter a negative intersection 
			// (or something like that). 
			// It's either tintersect <= 0 or sintersect > s.

			// Try again.
			dReal tintersect2 = 1e30;
			dReal tdelta = 0; // time offset
			dReal sddstart = vsampledpoints.at(vsampledpoints.size() - 2);
			dReal sdstart = (vsampledpoints.at(vsampledpoints.size() - 3) + 
					 tdelta*sddstart);
			dReal sstart = (vsampledpoints.at(vsampledpoints.size() - 4) + 
					tdelta*(sdstart + 0.5*tdelta*sddstart));

			checksample2 = FindEarliestProfileIntersection
			    (sstart, sdstart, sddstart, profile.integrationtimestep,
			     constraints.resprofileslist, sample.itprofile, tintersect2);

			if (checksample2.itprofile != constraints.resprofileslist.end() &&
			    checksample2.s > sprev && tintersect2 > TINY) {
			    // Here we can resolve the profile. Found a
			    // profile with a valid intersection.
			    vsampledpoints.push_back(checksample2.s);
			    vsampledpoints.push_back(checksample2.sd);
			    vsampledpoints.push_back(checksample2.sdd);
			    vsampledpoints.push_back(tintersect2);
			    checksample = checksample2;
			}
			else {
			    // Rosen says we should go on using the same profile.
			    busechecksample = false;
			    if (_print) {
				message = "[Reparameterize2] invalid intersection cannot be fixed.";
				std::cout << message << std::endl;
			    }
			}
		    }
		}
	    
		if (busechecksample) {
		    // We just encountered a new lowest profile.
		    // Stop stepping along sample.
		    if (_print) {
			message = "[Reparameterize2] Encountered a new valid lowest profile.";
			std::cout << message << std::endl;
		    }
		    
		    sprev = checksample.s;
		    sdprev = checksample.sd;
		    sddprev = checksample.sdd;
		    badded = true;
		    break;
		}
		else {
		    // This case means FindLowestProfileFast returned
		    // something, but when we tried to find an intersect,
		    // it somehow did not work out.
		    
		    // Reset checksample
		    checksample.itprofile = constraints.resprofileslist.end();
		    if (_print) {
			message = "[Reparameterize2] Newly found lowest profile cannot be used.";
			std::cout << message << std::endl;
		    }
		}
	    }
	    
	    // If we reach here, FindLowestProfileFast did not find
	    // any lower profile, or found but it did not work.
	    // Keep stepping along the same profile (sample).
	    BOOST_ASSERT(sprev <= s);
	    BOOST_ASSERT(profile.integrationtimestep - tprev > 0);
	    
	    vsampledpoints.push_back(s);
	    vsampledpoints.push_back(sd);
	    vsampledpoints.push_back(sdd);
	    vsampledpoints.push_back(profile.integrationtimestep - tprev);

	    // Successfully stepped. Keep going forward.
	    badded = true;
	    tprev = 0;
	    sprev = s;
	    sdprev = sd;
	    sddprev = sdd;
	    ++sindex;
	}
	
	if (!badded) {
	    message = "[Reparameterize2] badded is false. Something is wrong.";
	    std::cout << message << std::endl;
	    return -1;
	}

	if (sprev >= smax - TINY) {
	    // We have successfully reached the end of this trajectory
	    bsuccess = true;
	    break;
	}

	// If we reach here, it means sindex >= (int)profile.svect.size(). 

	// This implies that we have already stepped to the end of the
	// current profile (profile) but never found any valid
	// intersection with any other profile. Note that we have not
	// reached smax either.
	if (checksample.itprofile == constraints.resprofileslist.end()) {
	    bool bfindconnection = true;
	    do {
		bfindconnection = false;
		
		dReal sstart = std::min(smax, sprev - TINY);
		dReal sdthresh = 0.01;
		dReal sddistmin = 1e30;

		std::list<Profile>::const_iterator itprofilemin = constraints.resprofileslist.end();
		size_t sconnectindexmin = 0;
		
		for (std::list<Profile>::const_iterator itprofile = constraints.resprofileslist.begin(); 
		     itprofile != constraints.resprofileslist.end();
		     ++itprofile) {
		    // Skip unrelated profiles (including *sample.itprofile)
		    if ((itprofile == sample.itprofile) || (itprofile->svect.back() <= sstart + TINY)) {
			continue;
		    }

		    // Find the first index that has an s value greater than sstart
		    std::vector<dReal>::const_iterator its = 
			std::lower_bound(itprofile->svect.begin(), itprofile->svect.end(), sstart);
		    if (its == itprofile->svect.end()) {
			continue;
		    }

		    size_t sconnectindex = its - itprofile->svect.begin();
		    
		    if (sconnectindex += 1 >= itprofile->svect.size()) {
			if (*its - sstart <= TINY) {
			    // We have got the same point as sstart. Not use it.
			    continue;
			}
		    }

		    dReal sddist = fabs(itprofile->sdvect.at(sconnectindex) - sdprev);
		    bool bisclosetosd = sddist < sdthresh;
		    bool bisminclosetosd = sddistmin < sdthresh;
		    if (itprofilemin == constraints.resprofileslist.end() ||
			*its < itprofilemin->svect.at(sconnectindexmin) ||
			(!bisminclosetosd && bisclosetosd) ||
			(!bisminclosetosd && !bisclosetosd && sddist < sddistmin && 
			 itprofile->sdvect.at(sconnectindex) < itprofilemin->sdvect.at(sconnectindexmin)) ||
			(bisminclosetosd && bisclosetosd && sddist < sddistmin)) {
			itprofilemin = itprofile;
			sddistmin = sddist;
			sconnectindexmin = sconnectindex;
		    }
		}

		if (itprofilemin == constraints.resprofileslist.end()) {
		    message = str(boost::format("[Reparameterize2]"
						" Failed to find next profile at s = %.15e")
				  %sstart);
		    std::cout << message << std::endl;
		    return -1;
		}
		
		if (sconnectindexmin > 0) {
		    // itprofilemin->svect.at(sconnectindexmin) > sstart and
		    // itprofilemin->svect.at(sconnectindexmin - 1) <= sstart
		    // This is where the ramp starts.
		    sconnectindexmin--;
		}

		dReal sintersect = 0, sdintersect = 0, sddintersect = 0, tintersect = 0;
		dReal snext, sdnext, sddnext;
		size_t sconnectindex = sconnectindexmin;
		while(true) {
		    // In sample.itprofile, take the second to last point since that should have nonzero sdd
		    vsampledpoints.resize(vsampledpoints.size() - 4);
		    if (vsampledpoints.size() >= 4) {
			sprev = vsampledpoints[vsampledpoints.size() - 4];
			sdprev = vsampledpoints[vsampledpoints.size() - 3];
			sddprev = vsampledpoints[vsampledpoints.size() - 2];
		    }
		    else {
			sprev = sample.itprofile->svect.at(sample.itprofile->svect.size() - 2);
			sdprev = sample.itprofile->sdvect.at(sample.itprofile->svect.size() - 2);
			sddprev = sample.itprofile->sddvect.at(sample.itprofile->svect.size() - 2);
		    }
		    
		    sconnectindex = sconnectindexmin;

		    while (sconnectindex < itprofilemin->svect.size()) {
			snext = itprofilemin->svect.at(sconnectindex);
			sdnext = itprofilemin->sdvect.at(sconnectindex);
			sddnext = itprofilemin->sddvect.at(sconnectindex);

			if (fabs(sddprev) <= TINY) {
			    if (fabs(sddnext) <= TINY) {
				message = str(boost::format("[Reparameterize2]"
							    " sddprev & sddnext are both close to 0 at s = %.15e;"
							    "Don't know what to do.")%sprev);
				std::cout << message << std::endl;
				return -1;
			    }
			    dReal t = (sdprev - sdnext)/sddnext;
			    sintersect = snext + t*(sdnext + t*sddnext*0.5);
			    sdintersect = sdprev;
			    sddintersect = sddprev;
			    tintersect = (sintersect - sprev)/sdprev;
			}
			else if (fabs(sddnext) <= TINY) {
			    tintersect = (sdnext - sdprev)/sddprev;
                            sintersect = sprev + tintersect*(sdprev + tintersect*sddprev*0.5);
                            sddintersect = sddprev;
                            sdintersect = sdnext;
			}
			else {
			    dReal Aprofile = 1/(2*sddprev);
			    dReal Cprofile = sprev - sdprev*sdprev/(2*sddprev);
			
			    dReal Acheck = 1/(2*sddnext);
			    dReal Ccheck = snext - sdnext*sdnext/(2*sddnext);
			    
			    dReal Adiff = Aprofile - Acheck;
			    if (fabs(Adiff) <= TINY) {
				if (fabs(Ccheck - Cprofile) <= TINY) {
				    if (_print) {
					message = "[Reparameterize2] Adiff = Cdiff = 0";
					std::cout << message << std::endl;
				    }
				    // most lokely they intersect at the same place
				    sintersect = snext;
				    sdintersect = sdnext;
				    sddintersect = sddprev;
				    tintersect = (sdintersect - sdprev)/sddprev;
				}
				else {
				    if (_print) {
					message = "[Reparameterize2] Adiff = 0 but Cdiff != 0";
					std::cout << message << std::endl;
				    }
				    dReal tintersect2 = 1e30;
				    checksample2 = FindEarliestProfileIntersection\
					(sprev, sdprev, sddprev, profile.integrationtimestep*10, 
					 constraints.resprofileslist, sample.itprofile, tintersect2);
				    
				    if (checksample2.itprofile != constraints.resprofileslist.end() &&
					checksample2.s > sprev && tintersect2 > TINY) {
					sintersect = checksample2.s;
                                        sdintersect = checksample2.sd;
                                        sddintersect = checksample2.sdd;
                                        tintersect = tintersect2;
                                        itprofilemin = checksample2.itprofile;
                                        sconnectindex = checksample2.sindex;
                                        snext = itprofilemin->svect.at(sconnectindex);
                                        sdnext = itprofilemin->sdvect.at(sconnectindex);
                                        sddnext = itprofilemin->sddvect.at(sconnectindex);
				    }
				    else {
					return -1;
				    }
				}
			    }
			    else {
				dReal sdintersect2 = (Ccheck - Cprofile)/Adiff;
				if (sdintersect2 < 0) {
				    dReal tintersect2 = 1e30;
				    checksample2 = FindEarliestProfileIntersection\
					(sprev, sdprev, sddprev, profile.integrationtimestep*10, 
					 constraints.resprofileslist, sample.itprofile, tintersect2);
				    
				    if (checksample2.itprofile != constraints.resprofileslist.end() &&
					checksample2.s > sprev && tintersect2 > TINY) {
					sintersect = checksample2.s;
                                        sdintersect = checksample2.sd;
                                        sddintersect = checksample2.sdd;
                                        tintersect = tintersect2;
                                        itprofilemin = checksample2.itprofile;
                                        sconnectindex = checksample2.sindex;
                                        snext = itprofilemin->svect.at(sconnectindex);
                                        sdnext = itprofilemin->sdvect.at(sconnectindex);
                                        sddnext = itprofilemin->sddvect.at(sconnectindex);
				    }
				    else {
					return -1;
				    }
				}
				
				sdintersect = sqrt(sdintersect2);
				sintersect = Aprofile*sdintersect2 + Cprofile;
				tintersect = (sdintersect - sdprev)/sddprev;
				sddintersect = sddprev;
			    }
			}
			
			if (tintersect > 0) {
			    break;
			}
			
			++sconnectindex;
		    }

		    if (sconnectindex < itprofilemin->svect.size()) {
			break;
		    }
		    
		    if (vsampledpoints.size() == 0) {
			message = "[Reparameterize2] Still cannot find any intersection between profiles";
			std::cout << message << std::endl;
			return -1;
		    }
		}

		vsampledpoints.push_back(sintersect);
		vsampledpoints.push_back(sdintersect);
		vsampledpoints.push_back(sddnext);
		BOOST_ASSERT(tintersect > 0);
		vsampledpoints.push_back(tintersect);

		// If sconnectindexmin is the last value in the profile, we need to add it directly
		if (sconnectindex + 1 >= itprofilemin->svect.size()) {
		    dReal t2;
                    if( fabs(sddnext) > TINY ) {
                        t2 = (sdnext - sdintersect)/sddnext;
                        BOOST_ASSERT(t2 >= 0);
                        vsampledpoints.push_back(snext);
                        vsampledpoints.push_back(sdnext);
                        vsampledpoints.push_back(sddnext);
                        vsampledpoints.push_back(t2);
                    }
                    else {
                        t2 = (snext - sintersect)/sdnext;
			std::cout << snext << " " << sintersect << " " << sdnext << "\n";
                        BOOST_ASSERT(t2 >= 0);
                        vsampledpoints.push_back(snext);
                        vsampledpoints.push_back(sdnext);
                        vsampledpoints.push_back(sddnext);
                        vsampledpoints.push_back(t2);
                    }

                    if( snext >= smax - TINY ) {
                        bsuccess = true;
                        break;
                    }
		    
		    bfindconnection = true;
                    sample.itprofile = itprofilemin;
                    sample.sindex = sconnectindex;
                    sample.s = snext;
                    sample.sd = sdnext;
                    sample.sdd = sddnext;
                    sample.t = t2;
                    sprev = sample.s;
                    sdprev = sample.sd;
                    sddprev = sample.sdd;
		}
		else {
		    dReal t2 = (sdintersect - sdnext)/sddnext;
                    checksample.itprofile = itprofilemin;
                    checksample.sindex = sconnectindex;
                    checksample.s = sintersect;
                    checksample.sd = sdnext;
                    checksample.sdd = sddnext;
                    checksample.t = t2;
                    bool bsuccess = false;
                    // there are cases where sintersect is greater
                    // than itprofilemin->svect[sconnectindex+1],
                    // which causes asserts
                    while (checksample.sindex + 1 < (int)itprofilemin->svect.size()) {
                        if (sintersect + TINY <= itprofilemin->svect[checksample.sindex + 1]) {
                            if (checksample.t <= itprofilemin->integrationtimestep) {
                                bsuccess = true;
                                break;
                            }
                        }
                        // get the next index
                        checksample.sindex++;
                        checksample.t -= itprofilemin->integrationtimestep;
                    }
                    if( !bsuccess ) {
                        message = "ramp failed, so try again (infinite loop?)\n";
			std::cout << message << std::endl;
                    }
		}

	    } while (bfindconnection);
	}

	sample = checksample;

    }
    
    // Now we successfully extracted s, sd, sdd, t from all the lowest profiles.
    // Start reparameterizing!
    
    {
	std::list<Chunk>::iterator itcurrentchunk = intraj.chunkslist.begin();
	std::list<Chunk> newchunkslist;
	int currentchunkindex = 0;
	dReal processedcursor = 0;
	size_t sindex = 0;

	while (sindex < vsampledpoints.size()) {
	    dReal s = vsampledpoints[sindex];
	    dReal sd = vsampledpoints[sindex + 1];
	    dReal sdd = vsampledpoints[sindex + 2];
	    dReal tdelta = vsampledpoints[sindex + 3]; // this will be newchunk.duration
	    
	    SPieceToChunks(s, sd, sdd, tdelta, currentchunkindex, processedcursor,
			   itcurrentchunk, newchunkslist);
	    
	    sindex += 4;
	}
	
	if (newchunkslist.size() < 1) {
	    message = "[Reparameterize2] newchunkslist is empty.";
	    std::cout << message << std::endl;
	    return -1;
	}
	
	restrajectory = Trajectory(newchunkslist);
    }
    
    // Successful!
    message = "[Reparameterize2] Successful.";
    std::cout << message << std::endl;
    return 1;
}


void Trajectory::Write(std::stringstream& ss) {
    std::list<Chunk>::iterator itchunk = chunkslist.begin();
    while(itchunk!=chunkslist.end()) {
        itchunk->Write(ss);
        itchunk++;
    }
}


} // End namespace TOPP
