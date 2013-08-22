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


void PiecewisePolynomialTrajectory::Reparameterize(const Profile& profile, PiecewisePolynomialTrajectory& newtrajectory){
    // For now, supports only reparameterization of 3rd order polynomial trajectories
    assert(degree <= 3);

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


void PiecewisePolynomialTrajectory::Reintegrate(dReal reintegrationtimestep, PiecewisePolynomialTrajectory& newtrajectory){

    dReal dt = reintegrationtimestep;
    dReal dtsq = dt*dt;
    dReal tcur = 0, a;

    std::vector<dReal> qcur(dimension), qdcur(dimension), qdd(dimension);
    std::vector<dReal> coefficientsvector(3);
    std::vector<Polynomial> polynomialsvector(dimension);
    std::list<Chunk> chunkslist;

    Eval(0,qcur);
    Evald(0,qdcur);

    while(tcur<=duration) {
        Evaldd(tcur,qdd);
        for(int i=0; i<dimension; i++) {
            a = qdd[i]*0.5;
            coefficientsvector[0] = qcur[i];
            coefficientsvector[1] = qdcur[i];
            coefficientsvector[2] = a;
            qcur[i] += qdcur[i]*dt+ 0.5*a*dtsq;
            qdcur[i] += a*dt;
            polynomialsvector[i] = Polynomial(coefficientsvector);
        }
        chunkslist.push_back(Chunk(dt,polynomialsvector));
        tcur += dt;
    }

    // One last time
    tcur = tcur-dt;
    dt = duration-tcur;
    dtsq = dt*dt;
    Evaldd(tcur,qdd);
    for(int i=0; i<dimension; i++) {
        a = qdd[i]*0.5;
        coefficientsvector[0] = qcur[i];
        coefficientsvector[1] = qdcur[i];
        coefficientsvector[2] = a;
        qcur[i] += qdcur[i]*dt+ 0.5*a*dtsq;
        qdcur[i] += a*dt;
        polynomialsvector[i]=Polynomial(coefficientsvector);
    }
    chunkslist.push_back(Chunk(dt,polynomialsvector));


    newtrajectory = PiecewisePolynomialTrajectory(chunkslist);


}







/////////////////////////////////////////////////////////////////////

void PiecewisePolynomialTrajectory::Convert3(dReal chunklength, PiecewisePolynomialTrajectory& newtrajectory){
    assert(duration>3*chunklength);
    int n3chunks = (int) duration/3/chunklength;
    dReal T = duration/n3chunks;
    std::list<Chunk> reslist, finalreslist;

    std::vector<dReal> q0(dimension), qd0(dimension), qdd0(dimension), q1(dimension), qd1(dimension), qdd1(dimension);

    Eval(0,q0);
    Evald(0,qd0);
    Evaldd(0,qdd0);
    dReal tnext;

    for(int i=0; i<n3chunks; i++) {
        tnext = (i+1)*T;
        Eval(tnext,q1);
        Evald(tnext,qd1);
        Evaldd(tnext,qdd1);
        Interpolate3(T,q0,qd0,qdd0,q1,qd1,qdd1,reslist);
        finalreslist.splice(finalreslist.end(),reslist);
        std::swap(q0,q1);
    }

    newtrajectory = PiecewisePolynomialTrajectory(finalreslist);
}

void Interpolate3(dReal T, const std::vector<dReal>& q0,const std::vector<dReal>& qd0, const std::vector<dReal>& qdd0, const std::vector<dReal>& q1,const std::vector<dReal>& qd1, const std::vector<dReal>& qdd1, std::list<Chunk>& reslist){

    std::vector<Polynomial> polyvector0, polyvector1, polyvector2;
    std::vector<dReal> coefficientsvector(4);

    //Chunk0 : a0*t^3+b0*t^2+c0*t+d0
    arma::mat A,B,X;
    dReal t = T/3;
    dReal t2 = t*t;
    dReal t3 = t*t2;

    //  a0   b0   c0   d0   a1   b1   c1   d1   a2   b2   c2   d2
    A << 0 << 0 << 0 << 1 << 0 << 0 << 0 << 0 << 0 << 0 << 0 << 0 << arma::endr
      << 0 << 0 << 1 << 0 << 0 << 0 << 0 << 0 << 0 << 0 << 0 << 0 << arma::endr
      << 0 << 2 << 0 << 0 << 0 << 0 << 0 << 0 << 0 << 0 << 0 << 0 << arma::endr
      << t3<< t2<< t << 1 << 0 << 0 << 0 << -1<< 0 << 0 << 0 << 0 << arma::endr
      <<3*t2<<2*t<<1 << 0 << 0 << 0 << -1 << 0 << 0 << 0 << 0 << 0 << arma::endr
      <<6*t<< 2 << 0 << 0 << 0 << -2 << 0 << 0 << 0 << 0 << 0 << 0 << arma::endr
      << 0 << 0 << 0 << 0 << t3<< t2<< t << 1 << 0 << 0 << 0 << -1<< arma::endr
      << 0 << 0 << 0 << 0 <<3*t2<<2*t<<1 << 0 << 0 << 0 << -1 << 0 << arma::endr
      << 0 << 0 << 0 << 0 <<6*t<< 2 << 0 << 0 << 0 << -2 << 0 << 0 << arma::endr
      << 0 << 0 << 0 << 0 << 0 << 0 << 0 << 0 << t3<< t2<< t << 1 << arma::endr
      << 0 << 0 << 0 << 0 << 0 << 0 << 0 << 0 << 3*t2<<2*t<<1 << 0<< arma::endr
      << 0 << 0 << 0 << 0 << 0 << 0 << 0 << 0 << 6*t<< 2 << 0 << 0<< arma::endr;

    for(int i=0; i<q0.size(); i++) {
        B.set_size(0);
        B << q0[i] << qd0[i] << qdd0[i] << 0 << 0 << 0 << 0 << 0 << 0 << q1[i] << qd1[i] << qdd1[i] << arma::endr;
        X = arma::solve(A,arma::trans(B));
        coefficientsvector[0] = X[3];
        coefficientsvector[1] = X[2];
        coefficientsvector[2] = X[1];
        coefficientsvector[3] = X[0];
        polyvector0.push_back(Polynomial(coefficientsvector));
        coefficientsvector[0] = X[7];
        coefficientsvector[1] = X[6];
        coefficientsvector[2] = X[5];
        coefficientsvector[3] = X[4];
        polyvector1.push_back(Polynomial(coefficientsvector));
        coefficientsvector[0] = X[11];
        coefficientsvector[1] = X[10];
        coefficientsvector[2] = X[9];
        coefficientsvector[3] = X[8];
        polyvector2.push_back(Polynomial(coefficientsvector));
    }

    reslist.resize(0);
    reslist.push_back(Chunk(t,polyvector0));
    reslist.push_back(Chunk(t,polyvector1));
    reslist.push_back(Chunk(t,polyvector2));

}



} // End namespace TOPP
