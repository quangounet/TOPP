#ifndef TRAJECTORY_H
#define TRAJECTORY_H

namespace TOPP {

typedef double dReal;

class Constraints;  // defined in TOPP.h


class Polynomial {
public:
    Polynomial(const std::vector<dReal>& coefficientsvector);
    Polynomial(std::string& s);
    Polynomial(){
    }
    void InitFromCoefficientsVector(const std::vector<dReal>&coefficientsvector);
    int degree;
    // Weak coefficients first
    std::vector<dReal> coefficientsvector;
    std::vector<dReal> coefficientsvectord;
    std::vector<dReal> coefficientsvectordd;
    dReal Eval(dReal s) const;
    dReal Evald(dReal s) const;
    dReal Evaldd(dReal s) const;
    /// \brief use derivative of coefficientsvectordd to evalute
    dReal Evalddd(dReal s) const;
    /// \brief use second derivative of coefficientsvectordd to evalute
    dReal Evaldddd(dReal s) const;
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
    void Eval(dReal s, std::vector<dReal>&q) const;
    void Evald(dReal s, std::vector<dReal>&qd) const;
    void Evaldd(dReal s, std::vector<dReal>&qdd) const;
    /// \brief use derivative of coefficientsvectordd to evalute
    void Evalddd(dReal s, std::vector<dReal>&qddd) const;
    /// \brief use second derivative of coefficientsvectordd to evalute
    void Evaldddd(dReal s, std::vector<dReal>&qddd) const;
    void Write(std::stringstream& ss);
};


class Trajectory {
public:
    Trajectory(const std::list<Chunk>& chunkslist);
    Trajectory(const std::string& trajectorystring);
    Trajectory(){
    }
    void InitFromChunksList(const std::list<Chunk>&chunkslist);
    void InitFromString(const std::string& trajectorystring);

    int dimension;
    dReal duration;
    int degree;

    std::list<Chunk> chunkslist;
    std::list<dReal> chunkdurationslist;
    std::list<dReal> chunkcumulateddurationslist;

    // Return the index of the chunk that contains s, remainder indicates the s-time elapsed between the beginning of that chunk and s
    void FindChunkIndex(dReal s, int& index, dReal& remainder) const;
    void Eval(dReal s, std::vector<dReal>&q) const;
    void Evald(dReal s, std::vector<dReal>&qd) const;
    void Evaldd(dReal s, std::vector<dReal>&qdd) const;

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
    int Reparameterize(Constraints& constraints, Trajectory& restrajectory, dReal smax = 0);

    // Write the trajectory to the stream
    void Write(std::stringstream& ss);
};


} // End namespace TOPP


#endif // TRAJECTORY_H
