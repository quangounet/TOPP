from numpy import *
from pylab import *
import StringIO
import bisect



def ProfileFromLines(lines):
    l = lines[0]
    [duration,dt] = [double(x) for x in l.split(' ')]
    l = lines[1]
    sarray = array([double(x) for x in l.split(' ')])
    l = lines[2]
    sdarray = array([double(x) for x in l.split(' ')])
    return [duration,dt,sarray,sdarray]


def ProfilesFromString(s):
    s = s.strip(" \n")
    profileslist = []
    lines = [l.strip(" \n") for l in s.split('\n')]
    n = len(lines) / 3
    for i in range(n):        
        profileslist.append(ProfileFromLines(lines[3*i:3*i+3]))
    return profileslist


def VectorFromString(s):
    s = s.strip(" \n")
    return array([double(x) for x in s.split(' ')])


class Polynomial():
    def __init__(self,polynomialstring):
        self.coefficientsvector = VectorFromString(polynomialstring)
        self.degree = len(self.coefficientsvector)-1
        self.coefficientsvectord = zeros(self.degree)
        self.coefficientsvectordd = zeros(self.degree-1)
        for i in range(1,self.degree+1):
            self.coefficientsvectord[i-1] = i*self.coefficientsvector[i]
        for i in range(1,self.degree):
            self.coefficientsvectordd[i-1] = i*self.coefficientsvectord[i]

    def Eval(self,s):
        res = 0
        for i in range(self.degree,-1,-1):
            res = res*s + self.coefficientsvector[i];
        return res

    def Evald(self,s):
        res = 0
        for i in range(self.degree-1,-1,-1):
            res = res*s + self.coefficientsvectord[i];
        return res

    def Evaldd(self,s):
        res = 0
        for i in range(self.degree-2,-1,-1):
            res = res*s + self.coefficientsvectordd[i];
        return res

    def Write(self):
        ss = "";
        for i in range(0,self.degree+1):
            ss += str(self.coefficientsvector[i])
        return ss


class Chunk():
    def __init__(self,duration,polynomialsvector):
        self.polynomialsvector = polynomialsvector;
        self.dimension = len(polynomialsvector)
        self.duration = duration
        self.degree = polynomialsvector[0].degree

    def Eval(self,s):
        q = zeros(self.dimension)
        for i in range(self.dimension):
            q[i] = self.polynomialsvector[i].Eval(s)
        return q

    def Evald(self,s):
        qd = zeros(self.dimension)
        for i in range(self.dimension):
            qd[i] = self.polynomialsvector[i].Evald(s)
        return qd

    def Evaldd(self,s):
        qdd = zeros(self.dimension)
        for i in range(self.dimension):
            qdd[i] = self.polynomialsvector[i].Evaldd(s)
        return qdd

    def Write(self):
        ss = str(self.duration) + "\n"
        ss += str(self.dimension) + "\n"
        for i in range(self.dimension):
            ss += self.polynomialsvector[i].Write() + "\n"
        return ss



class PieceWisePolyTrajectory():
    
    def __init__(self,trajectorystring):
        buff = StringIO.StringIO(trajectorystring)
        chunkslist = []
        while(buff.pos<buff.len):
            duration = double(buff.readline())
            dimension = int(buff.readline())
            polynomialsvector = [];
            for i in range(dimension):
                polynomialsvector.append(Polynomial(buff.readline()))
            chunkslist.append(Chunk(duration,polynomialsvector))
        self.InitFromChunkslist(chunkslist)
            

    def InitFromChunkslist(self,chunkslist):
        self.chunkslist = chunkslist
        self.dimension = self.chunkslist[0].dimension
        self.degree = self.chunkslist[0].degree
        self.duration = 0
        self.chunkcumulateddurationslist = []
        for c in chunkslist:
            self.chunkcumulateddurationslist.append(self.duration)
            self.duration += c.duration        

    def FindChunkIndex(self,s):
        if(s==0):
            s = 1e-10
        i = bisect.bisect_left(self.chunkcumulateddurationslist,s)-1
        remainder = s - self.chunkcumulateddurationslist[i]
        return i,remainder
    
    def Eval(self,s):
        i,remainder = self.FindChunkIndex(s)        
        return self.chunkslist[i].Eval(remainder)
    
    def Evald(self,s):
        i,remainder = self.FindChunkIndex(s)
        return self.chunkslist[i].Evald(remainder)
    
    def Evaldd(self,s):
        i,remainder = self.FindChunkIndex(s)
        return self.chunkslist[i].Evaldd(remainder)
        
    def Plot(self,dt,f=''):
        tvect = arange(0,self.duration+dt,dt)
        qvect = array([self.Eval(t) for t in tvect])
        plot(tvect,qvect,f)
            
    def Plotd(self,dt,f=''):
        tvect = arange(0,self.duration+dt,dt)
        qvect = array([self.Evald(t) for t in tvect])
        plot(tvect,qvect,f)
    
    def Plotdd(self,dt,f=''):
        tvect = arange(0,self.duration+dt,dt)
        qvect = array([self.Evaldd(t) for t in tvect])
        plot(tvect,qvect,f)
            
 

