#include "PiecewisePolynomialTrajectory.h"


namespace TOPP {


///////////////////////// Polynomial ///////////////////////

Polynomial::Polynomial(const std::vector<dReal>& coefficientsvector0){
    coefficientsvector = coefficientsvector0;
    degree = coefficientsvector.size()-1;
    // Construct first- and second-order derivative polynomials
    coefficientsvectord.resize(0);
    coefficientsvectordd.resize(0);
    for(int i=1; i<=degree; i++) {
        coefficientsvectord.push_back(i*coefficientsvector[i]);
    }
    for(int i=1; i<=degree-1; i++) {
        coefficientsvectordd.push_back(i*coefficientsvectord[i]);
    }
}

// Evaluate polynomials using Horner's method
dReal Polynomial::Eval(dReal s){
    dReal res = 0;
    for(int i=degree; i>=0; i--) {
        res = res*s + coefficientsvector[i];
    }
    return res;
}

dReal Polynomial::Evald(dReal s){
    dReal res = 0;
    for(int i=degree-1; i>=0; i--) {
        res = res*s + coefficientsvectord[i];
    }
    return res;
}


dReal Polynomial::Evaldd(dReal s){
    dReal res = 0;
    for(int i=degree-2; i>=0; i--) {
        res = res*s + coefficientsvectordd[i];
    }
    return res;
}



///////////////////////// Chunk ///////////////////////


Chunk::Chunk(dReal duration0, const std::vector<Polynomial>& polynomialsvector0){
    polynomialsvector = polynomialsvector0;
    dimension = polynomialsvector.size();
    duration = duration0;
    assert(dimension>0);
    degree = polynomialsvector[0].degree;
    // All polynomials must have the same degree
    for(int i=1; i<dimension; i++) {
        assert(degree == polynomialsvector[i].degree);
    }
}


void Chunk::Eval(dReal s, std::vector<dReal>&q){
    assert(s >= 0);
    assert(s <= duration);
    assert(dimension == q.size());
    for(int i=0; i<dimension; i++) {
        q[i] = polynomialsvector[i].Eval(s);
    }
}

void Chunk::Evald(dReal s, std::vector<dReal>&qd){
    assert(s >= 0);
    assert(s <= duration);
    assert(dimension == qd.size());
    for(int i=0; i<dimension; i++) {
        qd[i] = polynomialsvector[i].Evald(s);
    }
}

void Chunk::Evaldd(dReal s, std::vector<dReal>&qdd){
    assert(s >= 0);
    assert(s <= duration);
    assert(dimension == qdd.size());
    for(int i=0; i<dimension; i++) {
        qdd[i] = polynomialsvector[i].Evaldd(s);
    }
}



/////////////// PiecewisePolynomialTrajectory ///////////////////////

PiecewisePolynomialTrajectory::PiecewisePolynomialTrajectory(const std::list<Chunk>& chunkslist0){
    chunkslist = chunkslist0;
    assert(chunkslist.size()>0);
    dimension = chunkslist.front().dimension;
    degree = chunkslist.front().degree;

    duration = 0;
    chunkdurationslist.resize(0);
    chunkcumulateddurationslist.resize(0);
    std::list<Chunk>::iterator itchunk = chunkslist.begin();
    while(itchunk != chunkslist.end()) {
        assert(degree == itchunk->degree);
        dReal chunkduration = itchunk->duration;
        if(chunkduration > TINY) {
            chunkdurationslist.push_back(chunkduration);
            chunkcumulateddurationslist.push_back(duration);
            duration += chunkduration;
        }
        itchunk++;
    }
    chunkcumulateddurationslist.push_back(duration);
}

void PiecewisePolynomialTrajectory::FindChunkIndex(dReal s, int& index, dReal& remainder){
    std::list<dReal>::iterator it = chunkcumulateddurationslist.begin();
    if(s <= TINY) {
        index = 0;
        remainder = 0;
        return;
    }
    index = 0;
    while(it != chunkcumulateddurationslist.end() && s > *it) {
        index++;
        it++;
    }
    index--;
    assert(index<=chunkslist.size()-1);
    it--;
    remainder = s-*it;
}

void PiecewisePolynomialTrajectory::Eval(dReal s, std::vector<dReal>&q){
    assert(s >= 0);
    assert(s <= duration);
    assert(dimension == q.size());
    int index;
    dReal remainder;
    FindChunkIndex(s,index,remainder);
    std::list<Chunk>::iterator itchunk = chunkslist.begin();
    advance(itchunk,index);
    itchunk->Eval(remainder,q);
}

void PiecewisePolynomialTrajectory::Evald(dReal s, std::vector<dReal>&qd){
    assert(s >= 0);
    assert(s <= duration);
    assert(dimension == qd.size());
    int index;
    dReal remainder;
    FindChunkIndex(s,index,remainder);
    std::list<Chunk>::iterator itchunk = chunkslist.begin();
    advance(itchunk,index);
    itchunk->Evald(remainder,qd);
}

void PiecewisePolynomialTrajectory::Evaldd(dReal s, std::vector<dReal>&qdd){
    assert(s >= 0);
    assert(s <= duration);
    assert(dimension == qdd.size());
    int index;
    dReal remainder;
    FindChunkIndex(s,index,remainder);
    std::list<Chunk>::iterator itchunk = chunkslist.begin();
    advance(itchunk,index);
    itchunk->Evaldd(remainder,qdd);
}


void PiecewisePolynomialTrajectory::ComputeChunk(dReal t0, dReal tnext, dReal s, dReal sd, dReal sdd, const Chunk& currentchunk, Chunk& newchunk){
    assert(currentchunk.degree == 2);
    dReal a0, a1, a2, b0, b1, b2, b3, b4, u0, u1, u2;
    std::vector<dReal> coefficientsvector;
    std::vector<Polynomial> polynomialsvector;
    // current chunk: u0 + u1*s + u2*s^2
    // new chunk : c0 + c1*t + c2*t^2 + c3*t^3 + c4*t^4
    a0 = s + sd*t0 + 0.5*sdd*t0*t0;
    a1 = sd + sdd*t0;
    a2 = 0.5*sdd;
    b0 = a0*a0;
    b1 = 2*a0*a1;
    b2 = 2*a0*a2+a1*a1;
    b3 = 2*a1*a2;
    b4 = a2*a2;
    for(int i=0; i<currentchunk.dimension; i++) {
        u0 = currentchunk.polynomialsvector[i].coefficientsvector[0];
        u1 = currentchunk.polynomialsvector[i].coefficientsvector[1];
        u2 = currentchunk.polynomialsvector[i].coefficientsvector[2];
        coefficientsvector.resize(0);
        coefficientsvector.push_back(u0 + u1*a0 + u2*b0); //c0
        coefficientsvector.push_back(u1*a1 + u2*b1);      //c1
        coefficientsvector.push_back(u1*a2 + u2*b2);      //c2
        coefficientsvector.push_back(u2*b3);              //c3
        coefficientsvector.push_back(u2*b4);              //c4
        polynomialsvector.push_back(Polynomial(coefficientsvector));
    }
    newchunk = Chunk(tnext-t0, polynomialsvector);
}


void PiecewisePolynomialTrajectory::Reparameterize(const Profile& profile, PiecewisePolynomialTrajectory& newtrajectory){
    // For now, supports only reparameterization of 2nd order polynomial trajectories
    assert(degree == 2);

    if(chunkslist.size() == 0) {
        return;
    }

    dReal ps, psd, psdd, psnext, s0, t, tnext;
    dReal sbeginchunk;
    int currentchunkindex, chunkindex;
    dReal processedcursor, remainingduration, remainder;
    Chunk newchunk;
    std::list<Chunk> newchunkslist;

    std::list<Chunk>::iterator itcurrentchunk = chunkslist.begin();
    std::list<dReal>::iterator itcumduration = chunkcumulateddurationslist.begin();
    sbeginchunk = *itcumduration;
    itcumduration++;
    currentchunkindex = 0;
    processedcursor = 0;
    t = 0;

    for(int i=0; i<profile.nsteps-1; i++) {
        t = 0;
        ps = profile.svect[i];
        psd = profile.sdvect[i];
        psdd = profile.sddvect[i];
        psnext = profile.svect[i+1];
        FindChunkIndex(psnext,chunkindex,remainder);
        // Process all chunks that have been overpassed
        while(currentchunkindex<chunkindex) {
            if(itcurrentchunk->duration-processedcursor>=TINY) {
                assert(SolveQuadraticEquation(ps-*itcumduration,psd,0.5*psdd,t,profile.integrationtimestep,tnext));
                ComputeChunk(t,tnext,ps-sbeginchunk,psd,psdd,*itcurrentchunk,newchunk);
                newchunkslist.push_back(newchunk);
                t = tnext;
            }
            currentchunkindex++;
            itcurrentchunk++;
            sbeginchunk = *itcumduration;
            itcumduration++;
            processedcursor = 0;
        }

        // Process current chunk
        assert(SolveQuadraticEquation((ps-sbeginchunk)-remainder,psd,0.5*psdd,t,profile.integrationtimestep,tnext));
        ComputeChunk(t,tnext,ps-sbeginchunk,psd,psdd,*itcurrentchunk,newchunk);
        newchunkslist.push_back(newchunk);
        processedcursor = remainder;
    }

    // Update trajectory
    newtrajectory = PiecewisePolynomialTrajectory(newchunkslist);
}


} // End namespace TOPP
