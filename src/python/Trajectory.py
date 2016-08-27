from numpy import *
import bisect
import pylab, scipy
import StringIO


from pylab import arange, array, double, plot, zeros


class Polynomial(object):
    @staticmethod
    def FromString(polynomial_string):
        s = polynomial_string.strip(" \n")
        coeff_list = [double(x) for x in s.split(' ')]
        return Polynomial(coeff_list)

    def __init__(self, coeff_list):
        # NB: we adopt the weak-term-first convention for inputs
        self.coeff_list = coeff_list
        self.q = poly1d(coeff_list[::-1])
        self.qd = polyder(self.q)
        self.qdd = polyder(self.qd)
        self.degree = self.q.order

    def pad_coeff_string(self, new_degree):
        while len(self.coeff_list) <= new_degree:
            self.coeff_list.append(0.)

    def Eval(self, s):
        return self.q(s)

    def Evald(self, s):
        return self.qd(s)

    def Evaldd(self, s):
        return self.qdd(s)

    def Evaldn(self, s, n):
        return polyder(self.q,n)(s)

    def __str__(self):
        return ' '.join(map(str, self.coeff_list))


class Chunk():
    def __init__(self, duration, poly_list):
        self.polynomialsvector = poly_list
        self.dimension = len(poly_list)
        self.duration = duration

        # TODO: current limitation in polynomials
        degrees = [poly.degree for poly in poly_list]
        self.degree = max(degrees)
        for poly in poly_list:
            poly.pad_coeff_string(self.degree)

    def Eval(self, s):
        q = zeros(self.dimension)
        for i in range(self.dimension):
            q[i] = self.polynomialsvector[i].Eval(s)
        return q

    def Evald(self, s):
        qd = zeros(self.dimension)
        for i in range(self.dimension):
            qd[i] = self.polynomialsvector[i].Evald(s)
        return qd

    def Evaldd(self, s):
        qdd = zeros(self.dimension)
        for i in range(self.dimension):
            qdd[i] = self.polynomialsvector[i].Evaldd(s)
        return qdd

    def Evaldn(self, s, n):
        qdn = zeros(self.dimension)
        for i in range(self.dimension):
            qdn[i] = self.polynomialsvector[i].Evaldn(s,n)
        return qdn


    def __str__(self):
        chunks_str = '\n'.join(map(str, self.polynomialsvector))
        return '%f\n%d\n%s' % (self.duration, self.dimension, chunks_str)


class PiecewisePolynomialTrajectory():
    def __init__(self, chunkslist):
        self.chunkslist = chunkslist
        self.dimension = self.chunkslist[0].dimension
        self.degree = self.chunkslist[0].degree
        self.duration = 0
        self.chunkcumulateddurationslist = []
        for c in chunkslist:
            self.chunkcumulateddurationslist.append(self.duration)
            self.duration += c.duration

    @staticmethod
    def FromString(trajectorystring):
        buff = StringIO.StringIO(trajectorystring)
        chunkslist = []
        while buff.pos < buff.len:
            duration = double(buff.readline())
            dimension = int(buff.readline())
            poly_vector = []
            for i in range(dimension):
                poly_vector.append(Polynomial.FromString(buff.readline()))
            chunkslist.append(Chunk(duration, poly_vector))
        return PiecewisePolynomialTrajectory(chunkslist)

    def FindChunkIndex(self, s):
        if s == 0:
            s = 1e-10
        i = bisect.bisect_left(self.chunkcumulateddurationslist, s) - 1
        remainder = s - self.chunkcumulateddurationslist[i]
        return i, remainder

    def Eval(self, s):
        i, remainder = self.FindChunkIndex(s)
        return self.chunkslist[i].Eval(remainder)

    def Evald(self, s):
        i, remainder = self.FindChunkIndex(s)
        return self.chunkslist[i].Evald(remainder)

    def Evaldd(self, s):
        i, remainder = self.FindChunkIndex(s)
        return self.chunkslist[i].Evaldd(remainder)

    def Evaldn(self, s, n):
        i, remainder = self.FindChunkIndex(s)
        return self.chunkslist[i].Evaldn(remainder,n)

    def Plot(self, dt, f=''):
        tvect = arange(0, self.duration + dt, dt)
        qvect = array([self.Eval(t) for t in tvect])
        plot(tvect, qvect, f, linewidth=2)

    def Plotd(self, dt, f=''):
        tvect = arange(0, self.duration + dt, dt)
        qdvect = array([self.Evald(t) for t in tvect])
        plot(tvect, qdvect, f, linewidth=2)

    def Plotdd(self, dt, f=''):
        tvect = arange(0, self.duration + dt, dt)
        qddvect = array([self.Evaldd(t) for t in tvect])
        plot(tvect, qddvect, f, linewidth=2)

    def __str__(self):
        return '\n'.join([str(chunk) for chunk in self.chunkslist])


def CropChunk(c, s0, s1):
    # Compute the Taylor shift at s0
    polynomialsvector = []
    for i in range(c.dimension):
        coeffs = []
        Pin = c.polynomialsvector[i].q
        for n in range(c.degree+1):
            coeffs.append(1./scipy.misc.factorial(n)*Pin(s0))
            Pin = polyder(Pin)
        polynomialsvector.append(Polynomial(coeffs))
    return Chunk(s1-s0, polynomialsvector)
            

# Assumes that i0 < i1
def InsertIntoTrajectory(traj,traj2,s0,s1):
    i0,r0 = traj.FindChunkIndex(s0)
    i1,r1 = traj.FindChunkIndex(s1)
    c0 = traj.chunkslist[i0]
    c1 = traj.chunkslist[i1]
    chunk0 = CropChunk(c0, 0, r0)
    chunk1 = CropChunk(c1, r1, c1.duration) 
    tolerance = 0.05
    if linalg.linalg.norm(traj2.Eval(0)-c0.Eval(r0))>=tolerance :
        print "Position mismatch at s0 : ", linalg.linalg.norm(traj2.Eval(0)-c0.Eval(r0))
        return None
    if linalg.linalg.norm(traj2.Eval(traj2.duration)-c1.Eval(r1))>=tolerance:
        print "Position mismatch at s1 : ", linalg.linalg.norm(traj2.Eval(traj2.duration)-c1.Eval(r1))
        return None
    if linalg.linalg.norm(traj2.Evald(0)-c0.Evald(r0)) >= tolerance:
        print "Velocity mismatch at s0 : ", linalg.linalg.norm(traj2.Evald(0)-c0.Evald(r0))
        return None
    if linalg.linalg.norm(traj2.Evald(traj2.duration)-c1.Evald(r1)) >= tolerance:
        print "Velocity mismatch at s1: ", linalg.linalg.norm(traj2.Evald(traj2.duration)-c1.Evald(r1))
        return None
    newchunkslist = list(traj.chunkslist)
    for i in range(i1-i0+1):
        newchunkslist.pop(i0)
    newchunkslist.insert(i0,chunk1)
    traj2.chunkslist.reverse()
    for chunk in traj2.chunkslist:
        newchunkslist.insert(i0,chunk)
    newchunkslist.insert(i0,chunk0)
    return PiecewisePolynomialTrajectory(newchunkslist)


def Concatenate(traj1,traj2):
    c1 = list(traj1.chunkslist)
    c1.extend(traj2.chunkslist)
    return PiecewisePolynomialTrajectory(c1)


def SubTraj(traj,s0,s1=-1):
    newchunkslist = []
    if s1 == -1:
        s1 = traj.duration
    i0,r0 = traj.FindChunkIndex(s0)
    i1,r1 = traj.FindChunkIndex(s1)
    c0 = traj.chunkslist[i0]
    c1 = traj.chunkslist[i1]
    if i0 == i1 :
        newchunkslist.append(CropChunk(c0,r0,r1))
    else:
        newchunkslist.append(CropChunk(c0, r0, c0.duration)) 
        i = i0+1
        while i < i1:
            newchunkslist.append(traj.chunkslist[i])
            i = i+1
        newchunkslist.append(CropChunk(c1, 0, r1))
    return PiecewisePolynomialTrajectory(newchunkslist)
