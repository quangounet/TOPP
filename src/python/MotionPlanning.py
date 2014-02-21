# -*- coding: utf-8 -*-
# Copyright (C) 2013 Quang-Cuong Pham <cuong.pham@normalesup.org>
#
# This file is part of the Time-Optimal Path Parameterization (TOPP) library.
# TOPP is free software: you can redistribute it and/or modify
# it under the terms of the GNU Lesser General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU Lesser General Public License for more details.
#
# You should have received a copy of the GNU Lesser General Public License
# along with this program.  If not, see <http://www.gnu.org/licenses/>.

import string
from numpy import *
from pylab import *
from openravepy import *
import time
import Trajectory
import TOPPbindings
import TOPPpy
import TOPPopenravepy


REACHED = 0
ADVANCED = 1
TRAPPED = 2


class Config():
    def __init__(self,q):
        self.q = q
    def dist(self,c,robot=None):
        if hasattr(robot,"weights"):
            dqsq = [x*x for x in self.q-c.q]            
            return sqrt(dot(dqsq,robot.weights))
        else:
            return norm(self.q-c.q)
    def IsAt(self,c):
        return self.dist(c)<1e-10


class Vertex():
    def __init__(self,config):
        self.config = config
        self.parent = None

class Tree():
    def __init__(self,qroot):
        vroot = Vertex(Config(qroot))
        self.verticeslist = [vroot]
    def AddVertex(self,v,cnew):
        vnew = Vertex(cnew)
        vnew.parent = v
        self.verticeslist.append(vnew)
    def NearestNeighbor(self,c,robot):
        return min(self.verticeslist,key=lambda v:v.config.dist(c,robot))
        

class RRTInstance():
    def __init__(self,robot,qstart,qend,epsilon,nloops):
        self.robot = robot
        self.treestart = Tree(qstart)
        self.treeend = Tree(qend)
        self.epsilon = epsilon
        self.currentmin = norm(qstart-qend)
        for i in range(nloops):
            if(mod(i,2)==0):
                Ta = self.treestart
                Tb = self.treeend
            else:
                Ta = self.treeend
                Tb = self.treestart
            crand = Config(robot.RandomConfig())
            if(not robot.CheckCollisionConfig(crand.q)):
                print "Found one collision free config, current min dist: ", self.currentmin
                print i, len(self.treestart.verticeslist), len(self.treeend.verticeslist)
                if(self.Extend(Ta,crand)!=TRAPPED):
                    if(self.Connect(Tb,Ta.verticeslist[-1].config) == REACHED):
                        print "Found path at step ", i
                        return
        print "Failure"            

    def TraceToRoot(self,T,v):
        if(v.parent):
            l = self.TraceToRoot(T,v.parent)
            l.append(v.config.q)
            return l
        else:
            return [v.config.q]

    def BuildPath(self):
        l1 = self.TraceToRoot(self.treestart,self.treestart.verticeslist[-1])
        l2 = self.TraceToRoot(self.treeend,self.treeend.verticeslist[-1])
        l2.reverse()
        l2.pop(0)
        l1.extend(l2)
        return l1

    def NewConfig(self,v,c):
        if(c.dist(v.config)<=self.epsilon):
            testq = c.q
        else:
            dq = c.q-v.config.q
            testq = v.config.q + self.epsilon*dq/norm(dq)
        if(not self.robot.CheckCollisionSegment(c.q,v.config.q)):
            return True,Config(testq)
        else:
            return False,None

    def Extend(self,T,c):
        vnear = T.NearestNeighbor(c,self.robot)
        success,cnew = self.NewConfig(vnear,c)
        if(success):
            T.AddVertex(vnear,cnew)
            if(cnew.IsAt(c)):
                return REACHED
            else:
                return ADVANCED
        return TRAPPED

    def Connect(self,T,c):
        while True:
            vnear = T.NearestNeighbor(c,self.robot)
            d = norm(vnear.config.q-c.q)
            if d<self.currentmin:
                self.currentmin = d
            s = self.Extend(T,c)
            if s != ADVANCED:
                return s
                
def LinearShortcut(robot,path):
    totallength = 0
    for i in range(len(path)-1):
        totallength += norm(path[i+1]-path[i])
    s1 = rand()*totallength    
    for i1 in range(len(path)-1):
        d = norm(path[i1+1]-path[i1])
        if s1 <= d:
            break
        else:
            s1 -= d
    s2 = rand()*totallength    
    for i2 in range(len(path)-1):
        d = norm(path[i2+1]-path[i2])
        if s2 <= d:
            break
        else:
            s2 -= d
    if i2 == i1:
        return
    if i2 < i1:
        i = i1
        s = s1
        i1 = i2
        s1 = s2
        i2 = i
        s2 = s
    d1 = path[i1+1]-path[i1]
    d2 = path[i2+1]-path[i2]
    q1 = path[i1] + s1*d1/norm(d1)
    q2 = path[i2] + s2*d2/norm(d2)
    if not robot.CheckCollisionSegment(q1,q2):
        for i in range(i2-i1):
            path.pop(i1+1)
        path.insert(i1+1,q1)
        path.insert(i1+2,q2)
            
def RepeatLinearShortcut(robot,path,nshortcuts):
    i = 0
    while(i<len(path)-2):
        d0 = path[i+1]-path[i]
        d1 = path[i+2]-path[i+1]
        if abs(norm(d0)*norm(d1)-abs(dot(d0,d1)))<1e-5 or (not robot.CheckCollisionSegment(path[i],path[i+2])):
            path.pop(i+1)
        else:
            i = i+1
    for i in range(nshortcuts):
        totallength = 0
        for j in range(len(path)-1):
            totallength += norm(path[j+1]-path[j])
        print i, ":", totallength
        LinearShortcut(robot,path)


def MakeParabolicTrajectory(path,robot=None,constraintstring=None):
    chunkslist = []
    for i in range(len(path)-1):
        print "MakeParabolic:", i+1, "/", len(path)-1
        q0 = path[i]
        q1 = path[i+1]
        qd0 = 0 * q0
        qd1 = 0 * q1
        traj = Trajectory.PiecewisePolynomialTrajectory([Trajectory.MakeChunk(q0,q1,qd0,qd1,norm(q1-q0))])
        if robot:
            x = TOPPbindings.TOPPInstance(robot,"ZMPTorqueLimits",constraintstring,str(traj))
            ret = x.RunComputeProfiles(0.5,0.5)
            if(ret == 0):
                print "Failed reparameterizing segment: ", i
                x.WriteProfilesList()
                x.WriteSwitchPointsList()
                profileslist = TOPPpy.ProfilesFromString(x.resprofilesliststring)
                switchpointslist = TOPPpy.SwitchPointsFromString(x.switchpointsliststring)
                TOPPpy.PlotProfiles(profileslist,switchpointslist,10)
                axis([0,1,0,20])
                raw_input()                
                return None
            x.ReparameterizeTrajectory()
            x.WriteResultTrajectory()
            trajretimed = TOPPpy.PiecewisePolynomialTrajectory.FromString(x.restrajectorystring)
            #TOPPopenravepy.PlotZMP(robot,traj,trajretimed,baselimits,0.01,4)
            #TOPPpy.PlotKinematics(traj,trajretimed)
            chunkslist.extend(trajretimed.chunkslist)
            print "Done"
            #raw_input()
        else:
            chunkslist.extend(traj.chunkslist)
    return Trajectory.PiecewisePolynomialTrajectory(chunkslist)


def ParabolicShortcut(rrtrobot,constraintstring,orig_traj):
    s0 = orig_traj.duration*rand()
    s1 = orig_traj.duration*rand()
    if s0>s1:
        s = s1
        s1 = s0
        s0 = s
    if s1 - s0 < orig_traj.duration/20.:
        print "Too short"
        return None    
    q0 = orig_traj.Eval(s0)
    qd0 = orig_traj.Evald(s0)
    q1 = orig_traj.Eval(s1)
    qd1 = orig_traj.Evald(s1)
    d = norm(q1-q0)
    d = d*2/3. + randn()*d*1/3.
    if(d<0.1):
        return None
    traj = Trajectory.PiecewisePolynomialTrajectory([Trajectory.MakeChunk(q0,q1,qd0,qd1,d)])
    x = TOPPbindings.TOPPInstance(rrtrobot.robot,"ZMPTorqueLimits",constraintstring,str(traj))
    ret = x.RunComputeProfiles(1,1)
    #x.WriteProfilesList()
    #x.WriteSwitchPointsList()
    #profileslist = TOPPpy.ProfilesFromString(x.resprofilesliststring)
    #switchpointslist = TOPPpy.SwitchPointsFromString(x.switchpointsliststring)
    if(ret == 0):
        #print "Not parameterizable"
        #TOPPpy.PlotProfiles(profileslist,switchpointslist,10)
        #axis([0,traj.duration,0,20])
        #raw_input()
        return None
    if(x.resduration>s1-s0-1e-2):
        print "Not significant: ", x.resduration, s1-s0
        return None
    x.ReparameterizeTrajectory()
    x.WriteResultTrajectory()        
    trajretimed = TOPPpy.PiecewisePolynomialTrajectory.FromString(x.restrajectorystring)
    if(trajretimed.duration>s1-s0-1e-2):
        print "Not significant (2): ", trajretimed.duration, x.resduration, s1-s0
        return None
    for s in arange(0,trajretimed.duration,0.01):
        if rrtrobot.CheckCollisionConfig(trajretimed.Eval(s)):
            print "Collision"
            return None
    #TOPPpy.PlotProfiles(profileslist,switchpointslist,10)
    #TOPPpy.PlotKinematics(traj,trajretimed)
    return Trajectory.InsertIntoTrajectory(orig_traj,trajretimed,s0,s1)



def RepeatParabolicShortcut(rrtrobot,constraintstring,orig_traj,nshortcuts):
    nsuccess = 0
    cur_traj = orig_traj
    print "Original duration: ", orig_traj.duration
    for i in range(nshortcuts):
        print i
        new_traj = ParabolicShortcut(rrtrobot,constraintstring,cur_traj)
        if new_traj:
            nsuccess += 1
            cur_traj = new_traj
            print "*********** New duration: ", cur_traj.duration, " *****************"
    print nsuccess
    return cur_traj
        
