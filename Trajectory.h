#ifndef TRAJECTORY_H
#define TRAJECTORY_H

namespace TOPP {


class Constraints;  // defined in TOPP.h


class Polynomial {
public:
    Polynomial(const std::vector<dReal>& coefficientsvector);
    Polynomial(const std::string& s);
    Polynomial(){
    }
    void InitFromCoefficientsVector(const std::vector<dReal>&coefficientsvector);
    int degree;
    std::vector<dReal> coefficientsvector;
    std::vector<dReal> coefficientsvectord;
    std::vector<dReal> coefficientsvectordd;
    dReal Eval(dReal s);
    dReal Evald(dReal s);
    dReal Evaldd(dReal s);
    void Write(std::stringstream& ss);
};


class Chunk {
public:
    Chunk(dReal duration, const std::vector<Polynomial>& polynomialsvector);
    Chunk(){
    };
    int dimension;
    int degree;
    dReal duration;
    dReal sbegin, send;
    std::vector<Polynomial> polynomialsvector;
    void Eval(dReal s, std::vector<dReal>&q);
    void Evald(dReal s, std::vector<dReal>&qd);
    void Evaldd(dReal s, std::vector<dReal>&qdd);
    void Write(std::stringstream& ss);
};


class Trajectory {
public:
    Trajectory(const std::list<Chunk>& chunkslist);
    Trajectory(const std::string& trajectorystring);
    Trajectory(){
    }
    void InitFromChunksList(const std::list<Chunk>&chunkslist);

    int dimension;
    dReal duration;
    int degree;

    std::list<Chunk> chunkslist;
    std::list<dReal> chunkdurationslist;
    std::list<dReal> chunkcumulateddurationslist;

    void FindChunkIndex(dReal s, int& index, dReal& remainder);
    void Eval(dReal s, std::vector<dReal>&q);
    void Evald(dReal s, std::vector<dReal>&qd);
    void Evaldd(dReal s, std::vector<dReal>&qdd);

    // Reparameterize one chunk
    void ComputeChunk(dReal t, dReal t0, dReal s, dReal sd, dReal sdd, const
                      Chunk& currentchunk, Chunk& newchunk);

    // Reparameterize a S-piece into chunkslist
    // Inputs:
    // - s-piece : s+sd*t+0.5*sdd*t^2 for t in [0,T]
    // - index of current chunk (original trajectory)
    // - processedcursor : time instant up to which the current chunk has been
    //                     processed
    // - pointer to current chunk
    // Return the list of reparameterized chunks
    void SPieceToChunks(dReal s, dReal sd, dReal sdd, dReal T, int&
                        currentchunkindex, dReal& processedcursor,
                        std::list<Chunk>::iterator& itcurrentchunk, std::list<Chunk>& chunkslist);

    // Reparameterize the trajectory
    // The degree of the polynomials of restrajectory will be 2*d where d is
    // the degree of the polynomials in the original trajectory
    int Reparameterize(Constraints& constraints, Trajectory& restrajectory);

    // Write the trajectory to the stream
    void Write(std::stringstream& ss);
};


} // End namespace TOPP


#endif // TRAJECTORY_H
