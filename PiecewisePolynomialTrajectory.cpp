#include "PiecewisePolynomialTrajectory.h"


namespace TOPP {


///////////////////////// Polynomial ///////////////////////

void Polynomial::InitFromCoefficientsVector(const std::vector<dReal>&coefficientsvector0){
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



Polynomial::Polynomial(const std::vector<dReal>& coefficientsvector0){
    InitFromCoefficientsVector(coefficientsvector0);
}

Polynomial::Polynomial(const std::string& s){
    VectorFromString(s,coefficientsvector);
    InitFromCoefficientsVector(coefficientsvector);

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

void Polynomial::Write(std::stringstream& ss){
    for(int i=0; i<=degree; i++) {
        ss << std::setprecision(17) << coefficientsvector[i] << " ";
    }
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
    assert(s >= -TINY);
    assert(s <= duration+TINY);
    assert(dimension == q.size());
    for(int i=0; i<dimension; i++) {
        q[i] = polynomialsvector[i].Eval(s);
    }
}

void Chunk::Evald(dReal s, std::vector<dReal>&qd){
    assert(s >= -TINY);
    assert(s <= duration+TINY);
    assert(dimension == qd.size());
    for(int i=0; i<dimension; i++) {
        qd[i] = polynomialsvector[i].Evald(s);
    }
}

void Chunk::Evaldd(dReal s, std::vector<dReal>&qdd){
    assert(s >= -TINY);
    assert(s <= duration+TINY);
    assert(dimension == qdd.size());
    for(int i=0; i<dimension; i++) {
        qdd[i] = polynomialsvector[i].Evaldd(s);
    }
}

void Chunk::Write(std::stringstream& ss){
    ss << std::setprecision(17) <<  duration << "\n";
    ss << dimension << "\n";
    for(int i=0; i<dimension; i++) {
        polynomialsvector[i].Write(ss);
        ss << "\n";
    }
}



/////////////// PiecewisePolynomialTrajectory ///////////////////////


void PiecewisePolynomialTrajectory::InitFromChunksList(const std::list<Chunk>&chunkslist0){
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
            chunkcumulateddurationslist.push_back(duration);  // Cumulated durations list starts with a 0
            itchunk->sbegin = duration;
            duration += chunkduration;
            itchunk->send = duration;
        }
        itchunk++;
    }
    chunkcumulateddurationslist.push_back(duration);

}


PiecewisePolynomialTrajectory::PiecewisePolynomialTrajectory(const std::list<Chunk>& chunkslist0){
    InitFromChunksList(chunkslist0);
}

PiecewisePolynomialTrajectory::PiecewisePolynomialTrajectory(const std::string& trajectorystring){
    int buffsize = 255;
    char buff[buffsize];
    std::istringstream iss(trajectorystring);
    int dimension;
    dReal duration;
    std::vector<Polynomial> polynomialsvector;
    std::list<Chunk> chunkslist0;
    while(iss.good()) {
        iss.getline(buff,buffsize);
        duration = atof(buff);
        iss.getline(buff,buffsize);
        dimension = atoi(buff);
        polynomialsvector.resize(0);
        for(int i=0; i<dimension; i++) {
            iss.getline(buff,buffsize);
            polynomialsvector.push_back(Polynomial(std::string(buff)));
        }
        chunkslist0.push_back(Chunk(duration,polynomialsvector));
    }
    InitFromChunksList(chunkslist0);
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
    assert(s >= -TINY);
    assert(s <= duration+TINY);
    assert(dimension == q.size());
    int index;
    dReal remainder;
    FindChunkIndex(s,index,remainder);
    std::list<Chunk>::iterator itchunk = chunkslist.begin();
    advance(itchunk,index);
    itchunk->Eval(remainder,q);
}

void PiecewisePolynomialTrajectory::Evald(dReal s, std::vector<dReal>&qd){
    assert(s >= 0-TINY);
    assert(s <= duration+TINY);
    assert(dimension == qd.size());
    int index;
    dReal remainder;
    FindChunkIndex(s,index,remainder);
    std::list<Chunk>::iterator itchunk = chunkslist.begin();
    advance(itchunk,index);
    itchunk->Evald(remainder,qd);
}

void PiecewisePolynomialTrajectory::Evaldd(dReal s, std::vector<dReal>&qdd){
    assert(s >= -TINY);
    assert(s <= duration+TINY);
    assert(dimension == qdd.size());
    int index;
    dReal remainder;
    FindChunkIndex(s,index,remainder);
    std::list<Chunk>::iterator itchunk = chunkslist.begin();
    advance(itchunk,index);
    itchunk->Evaldd(remainder,qdd);
}


void PiecewisePolynomialTrajectory::ComputeChunk(dReal t0, dReal tnext, dReal s, dReal sd, dReal sdd, const Chunk& currentchunk, Chunk& newchunk){
    assert(currentchunk.degree <= 3);
    dReal a0, a1, a2, b0, b1, b2, b3, b4, c0, c1, c2, c3, c4, c5, c6, u0, u1, u2, u3;
    std::vector<dReal> coefficientsvector;
    std::vector<Polynomial> polynomialsvector;
    // current chunk : u0 + u1*s + u2*s^2
    // profile : s + sd*t + 0.5*sdd*t^2
    // new chunk : v0 + v1*t + v2*t^2 + v3*t^3 + v4*t^4 + v5*t^5 + v6*t^6;
    a0 = s + sd*t0 + 0.5*sdd*t0*t0;
    a1 = sd + sdd*t0;
    a2 = 0.5*sdd;
    b0 = a0*a0;
    b1 = 2*a0*a1;
    b2 = 2*a0*a2+a1*a1;
    b3 = 2*a1*a2;
    b4 = a2*a2;
    c0 = b0*a0;
    c1 = b1*a0 + b0*a1;
    c2 = b2*a0 + b1*a1 + b0*a2;
    c3 = b3*a0 + b2*a1 + b1*a2;
    c4 = b4*a0 + b3*a1 + b2*a2;
    c5 = b4*a1 + b3*a2;
    c6 = b4*a2;
    for(int i=0; i<currentchunk.dimension; i++) {
        u0 = currentchunk.polynomialsvector[i].coefficientsvector[0];
        u1 = 0;
        if(currentchunk.degree>=1) {
            u1 = currentchunk.polynomialsvector[i].coefficientsvector[1];
        }
        u2 = 0;
        if(currentchunk.degree>=2) {
            u2 = currentchunk.polynomialsvector[i].coefficientsvector[2];
        }
        u3 = 0;
        if(currentchunk.degree>=3) {
            u3 = currentchunk.polynomialsvector[i].coefficientsvector[3];
        }
        coefficientsvector.resize(0);
        coefficientsvector.push_back(u3*c0 + u2*b0 + u1*a0 + u0); //v0
        coefficientsvector.push_back(u3*c1 + u2*b1 + u1*a1);      //v1
        coefficientsvector.push_back(u3*c2 + u2*b2 + u1*a2);      //v2
        coefficientsvector.push_back(u3*c3 + u2*b3);              //v3
        coefficientsvector.push_back(u3*c4 + u2*b4);              //v4
        coefficientsvector.push_back(u3*c5);                      //v5
        coefficientsvector.push_back(u3*c6);                      //v6
        polynomialsvector.push_back(Polynomial(coefficientsvector));
    }
    newchunk = Chunk(tnext-t0, polynomialsvector);
}



void PiecewisePolynomialTrajectory::SPieceToChunks(dReal s, dReal sd, dReal sdd, dReal T, int& currentchunkindex, dReal& processedcursor, std::list<Chunk>::iterator& itcurrentchunk, std::list<Chunk>& chunkslist){

    dReal t = 0, tnext;
    dReal snext = s + T*sd + 0.5*T*T*sdd;
    int chunkindex;
    dReal remainder;
    Chunk newchunk;

    FindChunkIndex(snext,chunkindex,remainder);

    // Process all chunks that have been overpassed
    while(currentchunkindex<chunkindex) {
        if(itcurrentchunk->duration-processedcursor>=TINY) {
            assert(SolveQuadraticEquation(s-itcurrentchunk->send,sd,0.5*sdd,tnext,t,T));
            ComputeChunk(t,tnext,s-itcurrentchunk->sbegin,sd,sdd,*itcurrentchunk,newchunk);
            chunkslist.push_back(newchunk);
            t = tnext;
        }
        currentchunkindex++;
        itcurrentchunk++;
        processedcursor = 0;
    }

    // Process current chunk
    assert(SolveQuadraticEquation((s-itcurrentchunk->sbegin)-remainder,sd,0.5*sdd,tnext,t,T));
    ComputeChunk(t,tnext,s-itcurrentchunk->sbegin,sd,sdd,*itcurrentchunk,newchunk);
    chunkslist.push_back(newchunk);
    processedcursor = remainder;
}



void PiecewisePolynomialTrajectory::Reparameterize(std::list<Profile>& profileslist, dReal reparamtimestep, PiecewisePolynomialTrajectory& newtrajectory){

    dReal scur, sdcur, snext, sdnext, sdnext2, sdd;
    dReal dt = reparamtimestep;
    dReal dtsq = dt*dt;
    dReal dtmod;
    Profile profile;
    dReal tres;
    std::list<Chunk> newchunkslist;

    std::list<Chunk>::iterator itcurrentchunk = chunkslist.begin();
    int currentchunkindex = 0;
    dReal processedcursor = 0;


    // Reset currentindex
    std::list<Profile>::iterator it = profileslist.begin();
    while(it != profileslist.end()) {
        it->currentindex = 0;
        it++;
    }

    scur = 0;
    dReal t = 0;
    FindLowestProfile(scur,profile,tres,profileslist);
    sdcur = profile.Evald(tres);

    while(scur<duration) {
        sdd = profile.Evaldd(tres);
        sdnext = sdcur + dt*sdd;
        snext = scur + dt*sdcur + 0.5*dtsq*sdd;
        if(snext >= scur+TINY && FindLowestProfile(snext,profile,tres,profileslist)) {
            sdnext2 = profile.Evald(tres);
            dtmod = dt;
            if(std::abs(sdnext-sdnext2)>TINY2) {
                dtmod = 2*(snext-scur)/(sdnext2+sdcur);
                //std::cout << t << ": " << sdcur << " " << sdnext << "/" << sdnext2 << "   " << sdd << " ";
                sdd = (sdnext2-sdcur)/dtmod;
                //std::cout << sdd << "\n";
            }
            else{
                //std::cout << t << ": " << sdcur << " " << sdnext << "   " << sdd << "---------\n";
            }
            SPieceToChunks(scur,sdcur,sdd,dtmod,currentchunkindex,processedcursor,itcurrentchunk,newchunkslist);
        }
        else{
            break;
        }
        t+= dtmod;
        scur = snext;
        sdcur = sdnext2;
    }
    newtrajectory = PiecewisePolynomialTrajectory(newchunkslist);
}





void PiecewisePolynomialTrajectory::Write(std::stringstream& ss){
    std::list<Chunk>::iterator itchunk = chunkslist.begin();
    while(itchunk!=chunkslist.end()) {
        itchunk->Write(ss);
        itchunk++;
    }
}






} // End namespace TOPP
