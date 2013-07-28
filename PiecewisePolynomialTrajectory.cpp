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
	if(chunkduration > TINY){
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
    index = 0;
    while(it != chunkcumulateddurationslist.end() && s >= *it) {
        index++;
        it++;
    }
    index--;
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



void Reparameterize(const Profile& profile){
  if(chunkslist.size() == 0){
    return;
  }
  
  dReal s, sd, sdd, si, t, tnext;
  int currentchunkindex, chunkindex;
  dReal processedcursor, remainder;

  std::list<Chunk>::iterator itcurrentchunk = chunkslist.begin();
  std::list<dReal>::iterator itcumduration = chunkcumulateddurationslist.begin();
  itcumduration++;
  currentchunkindex = 0;
  processedcursor = 0;
  t = 0;

  for(int i=0; i<profile.nsteps-1; i++){
    s = profile.svect[i];
    sd = profile.sdvect[i];
    halfsdd = 0.5*profile.sddvect[i];
    snext = profile.svect[i+1];
    FindChunkIndex(snext,chunkindex,remainder);    
    // Process all chunks that have been overpassed
    while(currentchunkindex<chunkindex){
      // Make a new chunk from currentprocessedcursor to the end of the chunk
      tnext = SolveEq(s-*itcumduration,sd,halfsdd);
      // Make chunk from t to tnext

      currentchunkindex++;
      itcurrentchunk++;
      itcumduration++;
      processedcursor = 0;
      t = tnext;
    }

    // Process current chunk
    // Make a new chunk from currentprocessedcursor to remainder
    itcumduration--;
    tnext = SolveEq(s-(*itcumduration+remainder),sd,halfsdd);
    itcumduration++;
    // Make chunk from t to tnext


    t = tnext;
    processedcursor = remainder;
  }


}
