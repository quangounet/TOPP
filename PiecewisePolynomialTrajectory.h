#include "TOPP.h"


namespace TOPP {

class Polynomial {
public:
    Polynomial(const std::vector<dReal>& coefficientsvector);
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

    int dimension;
    int degree;
    dReal duration;
    std::list<Chunk> chunkslist;
    std::list<dReal> chunkdurationslist;
    std::list<dReal> chunkcumulateddurationslist;

    void FindChunkIndex(dReal s, int& index, dReal& remainder);
    void Eval(dReal s, std::vector<dReal>&q);
    void Evald(dReal s, std::vector<dReal>&qd);
    void Evaldd(dReal s, std::vector<dReal>&qdd);
    void Reparameterize(const Profile& profile); //Reparameterize in place

};
}
