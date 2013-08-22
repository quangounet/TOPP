#include "TOPP.h"


namespace TOPP {

class Polynomial {
public:
    Polynomial(const std::vector<dReal>& coefficientsvector);
    Polynomial(){
    }
    int degree;
    std::vector<dReal> coefficientsvector;
    std::vector<dReal> coefficientsvectord;
    std::vector<dReal> coefficientsvectordd;
    dReal Eval(dReal s);
    dReal Evald(dReal s);
    dReal Evaldd(dReal s);
};


class Chunk {
public:
    Chunk(dReal duration, const std::vector<Polynomial>& polynomialsvector);
    Chunk(){
    };
    int dimension;
    int degree;
    dReal duration;
    std::vector<Polynomial> polynomialsvector;
    void Eval(dReal s, std::vector<dReal>&q);
    void Evald(dReal s, std::vector<dReal>&qd);
    void Evaldd(dReal s, std::vector<dReal>&qdd);

};


class PiecewisePolynomialTrajectory : public Trajectory {
public:
    PiecewisePolynomialTrajectory(const std::list<Chunk>& chunkslist);
    PiecewisePolynomialTrajectory(){
    }
    int degree;
    std::list<Chunk> chunkslist;
    std::list<dReal> chunkdurationslist;
    std::list<dReal> chunkcumulateddurationslist;

    void FindChunkIndex(dReal s, int& index, dReal& remainder);
    void Eval(dReal s, std::vector<dReal>&q);
    void Evald(dReal s, std::vector<dReal>&qd);
    void Evaldd(dReal s, std::vector<dReal>&qdd);
    void ComputeChunk(dReal t, dReal t0, dReal s, dReal sd, dReal sdd, const Chunk& currentchunk, Chunk& newchunk);
    void Reparameterize(const Profile& profile, PiecewisePolynomialTrajectory& newtrajectory);
    void Reintegrate(dReal reintegrationtimestep, PiecewisePolynomialTrajectory& newtrajectory);
    void Convert3(dReal chunklength, PiecewisePolynomialTrajectory& newtrajectory);

};


void Interpolate3(dReal T, const std::vector<dReal>& q0,const std::vector<dReal>& qd0, const std::vector<dReal>& qdd0, const std::vector<dReal>& q1,const std::vector<dReal>& qd1, const std::vector<dReal>& qdd1, std::list<Chunk>& reslist);

}
