# -*- coding: utf-8 -*-
# Copyright (C) 2013 Quang-Cuong Pham <cuong.pham@normalesup.org>
#
# This file is part of the Time-Optimal Path Parameterization (TOPP) library.
# TOPP is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <http://www.gnu.org/licenses/>.


from numpy import *
from pylab import *
import string
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


def Interpolate3rdDegree(q0,q1,qd0,qd1,T):
    a=((qd1-qd0)*T-2*(q1-q0-qd0*T))/T**3
    b=(3*(q1-q0-qd0*T)-(qd1-qd0)*T)/T**2
    c=qd0
    d=q0
    return a,b,c,d

def ComputeConstraints(traj,amax,discrtimestep):
    # Sample the dynamics constraints
    ndiscrsteps = int((traj.duration+1e-10)/discrtimestep)+1;
    constraintstring = ""
    for i in range(ndiscrsteps):
        t = i*discrtimestep
        q = traj.Eval(t)
        qd = traj.Evald(t)
        qdd = traj.Evaldd(t)
        constraintstring += "\n" + string.join([str(x) for x in qd]) + " " + string.join([str(x) for x in -qd])
        constraintstring += "\n" + string.join([str(x) for x in qdd]) + " " + string.join([str(x) for x in -qdd]) 
        constraintstring += "\n" + string.join([str(x) for x in -amax]) + " " + string.join([str(x) for x in -amax])
    return constraintstring



def PlotProfiles(profileslist,figstart=0):
    figure(figstart)
    clf()
    hold('on')
    mvcbobrow = profileslist.pop(0)
    plot(mvcbobrow[2],mvcbobrow[3],'m--',linewidth=4)
    mvcdirect = profileslist.pop(0)
    plot(mvcdirect[2],mvcdirect[3],'c--',linewidth=4)
    for p in profileslist:
        plot(p[2],p[3],linewidth=2)
    axis([0,mvcbobrow[0],0,2*max([max(p[3]) for p in profileslist])])
    title('MVCs and profiles')


def PlotKinematics(traj0,traj1,dt=0.01,vmax=[],amax=[],figstart=0):
    colorcycle = ['r','g','b','m','c','y']
    colorcycle = colorcycle[0:traj0.dimension]
    Tmax = max(traj0.duration,traj1.duration)
    # Joint angles
    figure(figstart)
    clf()
    hold('on')
    ax=gca()        
    ax.set_color_cycle(colorcycle)
    traj0.Plot(dt,'--')
    traj1.Plot(dt)
    title('Joint values')
    # Velocity
    figure(figstart+1)
    clf()
    hold('on')
    ax=gca()        
    ax.set_color_cycle(colorcycle)
    traj0.Plotd(dt,'--')
    traj1.Plotd(dt)
    for v in vmax:
        plot([0,Tmax],[v,v],'-.')
    for v in vmax:
        plot([0,Tmax],[-v,-v],'-.')
    if(len(vmax)>0):
        Vmax = 1.2*max(vmax)
        axis([0,Tmax,-Vmax,Vmax])
    title('Joint velocities')
    # Acceleration
    figure(figstart+2)
    ax=gca()        
    ax.set_color_cycle(colorcycle)
    clf()
    hold('on')
    traj0.Plotdd(dt,'--')
    traj1.Plotdd(dt)
    for a in amax:
        plot([0,Tmax],[a,a],'-.')
    for a in amax:
        plot([0,Tmax],[-a,-a],'-.')
    if(len(amax)>0):
        Amax = 1.2*max(amax)
        axis([0,Tmax,-Amax,Amax])
    title('Joint accelerations')




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



class PiecewisePolynomialTrajectory():
    
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
        plot(tvect,qvect,f,linewidth=2)
            
    def Plotd(self,dt,f=''):
        tvect = arange(0,self.duration+dt,dt)
        qdvect = array([self.Evald(t) for t in tvect])
        plot(tvect,qdvect,f,linewidth=2)
    
    def Plotdd(self,dt,f=''):
        tvect = arange(0,self.duration+dt,dt)
        qddvect = array([self.Evaldd(t) for t in tvect])
        plot(tvect,qddvect,f,linewidth=2)
        

