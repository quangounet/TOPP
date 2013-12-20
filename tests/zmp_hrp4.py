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

import sys
sys.path.append('..')

import TOPPbindings
import TOPPpy
import TOPPopenravepy
import MotionPlanning
import time
import string
from pylab import *
from numpy import *
from openravepy import *

ion()

########################### Robot ################################
env = Environment()
robotfile = "../robots/hrp4r.dae"
baselinkname = "BODY"
robot = TOPPopenravepy.LoadFloat(env,robotfile,baselinkname) #Load a robot with dummy joints that mimick floating base
n=robot.GetDOF()
dof_lim=robot.GetDOFLimits()
vel_lim=robot.GetDOFVelocityLimits()
dof_lim[1][18]=0.5
dof_lim[1][28]=1.2
dof_lim[1][29]=0.8
robot.SetDOFLimits(dof_lim[0],dof_lim[1])
robot.SetDOFVelocityLimits(100*vel_lim)

M = array([[-0.81072723, -0.02540599, -0.58487254,  1.37787211],
       [ 0.58541559, -0.04056306, -0.80971799,  1.72345638],
       [-0.00315253, -0.99885393,  0.04775864,  0.70553762],
       [ 0.        ,  0.        ,  0.        ,  1.        ]])

#env.SetViewer('qtcoin')
#env.GetViewer().SetCamera(M)

chest = RaveGetEnvironment(1).GetKinBody('HRP4R').GetLink('CHEST_Y_LINK')
rsho = RaveGetEnvironment(1).GetKinBody('HRP4R').GetLink('R_SHOULDER_R_LINK')
lsho = RaveGetEnvironment(1).GetKinBody('HRP4R').GetLink('L_SHOULDER_R_LINK')
rhipy = RaveGetEnvironment(1).GetKinBody('HRP4R').GetLink('R_HIP_Y_LINK')
rhipp = RaveGetEnvironment(1).GetKinBody('HRP4R').GetLink('R_HIP_P_LINK')
lhipy = RaveGetEnvironment(1).GetKinBody('HRP4R').GetLink('L_HIP_Y_LINK')
lhipp = RaveGetEnvironment(1).GetKinBody('HRP4R').GetLink('L_HIP_P_LINK')
ignore_list=[[chest,rsho],[chest,lsho],[rhipy,rhipp],[lhipy,lhipp]]



pole = RaveCreateKinBody(env,'')
pole.SetName('Pole')
pole.InitFromBoxes(array([array([0,0,0,0.02,0.02,0.6])]),True)
g=pole.GetLinks()[0].GetGeometries()[0]
g.SetAmbientColor([1,0,0])
g.SetDiffuseColor([1,0,0])
env.Add(pole,True)
Tr = eye(4)
Tr[0:3,3]=[0.15,-0.25,0.6]
pole.SetTransform(Tr)


############################ Tunings ############################
discrtimestep = 1e-2
integrationtimestep = 1e-2
reparamtimestep = 1e-2
passswitchpointnsteps = 5
tuningsstring = "%f %f %f %d"%(discrtimestep,integrationtimestep,reparamtimestep,passswitchpointnsteps)


############################ Trajectory ############################

q0 = array([  1.37383090e-16,   1.37383090e-16,   1.37383090e-16,
         1.37383090e-16,   1.37383090e-16,   1.37383090e-16,
         1.37383090e-16,   1.37383090e-16,   0.00000000e+00,
         0.00000000e+00,   0.00000000e+00,   0.00000000e+00,
         0.00000000e+00,   0.00000000e+00,   0.00000000e+00,
         0.00000000e+00,   0.00000000e+00,  -1.32645019e-02,
        -4.00000000e-01,   2.00000006e+00,  -3.27249286e-01,
        -7.85398163e-03,   0.00000000e+00,   2.00712867e-02,
        -3.82052575e-01,   7.19250178e-01,  -3.27074692e-01,
        -1.91986222e-02,   1.39626330e-01,   0.00000000e+00,
         0.00000000e+00,   0.00000000e+00,  -5.23598696e-02,
        -1.74532927e-01,   0.00000000e+00,   0.00000000e+00,
         0.00000000e+00,   0.00000000e+00,   0.00000000e+00,
         3.20594657e-09,   0.00000000e+00,  -5.23598696e-02,
         1.50000000e+00,   0.00000000e+00,  -5.23598789e-01,
         0.00000000e+00,   0.00000000e+00,   1.37383090e-16,
        -3.20594658e-09,   1.37383090e-16,   0.00000000e+00,
         0.00000000e+00,   0.00000000e+00,   0.00000000e+00,
         0.00000000e+00,   0.00000000e+00])
q1 = array([  0.00000000e+00,   0.00000000e+00,   0.00000000e+00,
         0.00000000e+00,   0.00000000e+00,   0.00000000e+00,
         0.00000000e+00,   0.00000000e+00,  -1.86448479e-16,
         0.00000000e+00,  -1.86448479e-16,   0.00000000e+00,
        -1.86448479e-16,   0.00000000e+00,  -1.86448479e-16,
         0.00000000e+00,  -7.26978090e-01,  -1.32645018e-02,
         4.83489482e-01,   1.08640828e+00,  -3.27249427e-01,
        -7.85400698e-03,   0.00000000e+00,   2.00712867e-02,
        -3.82052575e-01,   7.19250178e-01,  -3.27074692e-01,
        -1.91986255e-02,   1.17184231e+00,   7.87388429e-01,
         0.00000000e+00,   0.00000000e+00,  -5.23596384e-02,
        -1.25348046e+00,   0.00000000e+00,  -1.43109325e+00,
         0.00000000e+00,   0.00000000e+00,   1.13053233e-16,
         7.05967061e-08,   0.00000000e+00,  -5.23596384e-02,
         1.23972591e+00,   1.50437154e+00,  -1.35968925e+00,
         0.00000000e+00,   0.00000000e+00,   0.00000000e+00,
         6.10426328e-08,   0.00000000e+00,   0.00000000e+00,
         0.00000000e+00,   0.00000000e+00,   0.00000000e+00,
         0.00000000e+00,   0.00000000e+00])

v=0.3
ndoffull = len(q0)
qd0=[v]*ndoffull
qd1=[v]*ndoffull
T = 1.5

activedofs = zeros(ndoffull)
activedofs[16] = 1 #R_HIP_Y
activedofs[17] = 1 #R_HIP_R
activedofs[18] = 1 #R_HIP_P
activedofs[19] = 1 #R_KNEE_P
activedofs[20] = 1 #R_ANKLE_P
activedofs[21] = 1 #R_ANKLE_R
activedofs[22] = 0 #L_HIP_Y
activedofs[23] = 0 #L_HIP_R
activedofs[24] = 0 #L_HIP_P
activedofs[25] = 0 #L_KNEE_P
activedofs[26] = 0 #L_ANKLE_P
activedofs[27] = 0 #L_ANKLE_R
activedofs[28] = 1 #CHEST_P
activedofs[29] = 1 #CHEST_Y
activedofs[32] = 1 #R_SHOULDER_P
activedofs[33] = 1 #R_SHOULDER_R
activedofs[34] = 1 #R_SHOULDER_Y
activedofs[35] = 1 #R_ELBOW_P
activedofs[41] = 1 #L_SHOULDER_P
activedofs[42] = 1 #L_SHOULDER_R
activedofs[43] = 1 #L_SHOULDER_Y
activedofs[44] = 1 #L_ELBOW_P

robot.activedofs = activedofs
robot.qdefault = array(q0)
ndof = int(sum(activedofs))
trajectorystring = "%f\n%d"%(T,ndof)
for i in range(ndoffull):
    if(activedofs[i])>0.1:
        if(abs(q0[i]-q1[i])>1e-5):
            a,b,c,d = TOPPpy.Interpolate3rdDegree(q0[i],q1[i],qd0[i],qd1[i],T)
        else:
            a,b,c,d = 0,0,0,q0[i]
        trajectorystring += "\n%f %f %f %f"%(d,c,b,a)


trajectorystring = """1.500000
22
0.0 0.0 1.99242076297 -1.32221433021
0.0 0.0 2.53280855643 -1.67759722295
0.0 0.0 0.0 0.0
0.0 0.0 0.0 0.0
0.0 0.0 0.967382109132 -0.642750868261
0.0 0.0 -0.519933641176 0.348891374481
0.0 0.0 2.89830573159 -1.93488398449
0.0 0.0 -2.14170672114 1.44939278452
0.0 0.0 0.0 0.0
0.0 0.0 0.0 0.0
0.0 0.0 0.0 0.0
0.0 0.0 0.0 0.0
0.0 0.0 0.0 0.0
0.0 0.0 3.13256805344 -2.08537050225
0.0 0.0 4.59613512542 -3.04279209114
0.0 0.0 4.81882377745 -3.21312086562
0.0 0.0 -3.18430860851 2.12838441106
0.0 0.0 0.0 0.0
0.0 0.0 0.0 0.0
0.0 0.0 0.0 0.0
0.0 0.0 0.0 0.0
0.0 0.0 0.0 0.0"""

traj0 = TOPPpy.PiecewisePolynomialTrajectory.FromString(trajectorystring)


activelinks = ones(len(robot.GetLinks()))
for i in range(len(activelinks)):
    if(robot.GetLinks()[i].GetMass()<0.1):
        activelinks[i] = 0
robot.activelinks = activelinks

#tvect,xzmp,yzmp = TOPPopenravepy.ComputeZMP(traj0,robot,0.01)

############################ Constraints ############################

taumin = -ones(ndof)*200
taumax = ones(ndof)*200
baseaabb = robot.GetLink("L_FOOT_LINK").ComputeAABB()
border2 = 0
xmax2 = baseaabb.pos()[0]+baseaabb.extents()[0]-border2
xmin2 = baseaabb.pos()[0]-baseaabb.extents()[0]+border2
ymax2 = baseaabb.pos()[1]+baseaabb.extents()[1]-border2
ymin2 = baseaabb.pos()[1]-baseaabb.extents()[1]+border2
zmplimits = [xmin2,xmax2,ymin2,ymax2]
vmax = ones(ndof)*3
constraintstring = string.join([str(x) for x in activedofs]) + "\n" + string.join([str(x) for x in activelinks]) + "\n" + string.join([str(x) for x in taumin]) + "\n" + string.join([str(a) for a in taumax]) + "\n" +  string.join([str(a) for a in zmplimits]) + "\n" + string.join([str(a) for a in vmax])




############################ RRT ############################

border = 0.015 #safety border of 0.02
xmax = baseaabb.pos()[0]+baseaabb.extents()[0]-border
xmin = baseaabb.pos()[0]-baseaabb.extents()[0]+border
ymax = baseaabb.pos()[1]+baseaabb.extents()[1]-border
ymin = baseaabb.pos()[1]-baseaabb.extents()[1]+border
baselimits = [xmin,xmax,ymin,ymax]

q0trimmed = array(TOPPopenravepy.Trim(robot,q0))
q1trimmed = array(TOPPopenravepy.Trim(robot,q1))

rrtrobot = TOPPopenravepy.HRP4Robot(robot,baselimits,0.05)
epsilon = 0.2
nloops = 10000000

TOPPopenravepy.CheckCollisionStaticStabilityConfig(robot,q0trimmed,baselimits)
TOPPopenravepy.CheckCollisionStaticStabilityConfig(robot,q1trimmed,baselimits)

# print "Initial distance: ", norm(q1trimmed-q0trimmed)

def callback(report,physics):
    for linkpair in ignore_list:
        if ((linkpair[0]==report.plink1) and (linkpair[1]==report.plink2)) or ((linkpair[1]==report.plink1) and (linkpair[0]==report.plink2)):
            #print "Ignored: ", report
            return 1
    #print "Not ignored: ", report
    return 0
  
collisioncallbackhandel = env.RegisterCollisionCallback(callback)


# t0 = time.time()
# rrt = MotionPlanning.RRTInstance(rrtrobot,q0trimmed,q1trimmed,epsilon,nloops)
# t1 = time.time()
# path = rrt.BuildPath()

# nshortcuts = 500
# MotionPlanning.RepeatLinearShortcut(rrtrobot,path,nshortcuts)

#traj0 = MotionPlanning.MakeParabolicTrajectory(path)

# for q in path:
#     robot.SetDOFValues(TOPPopenravepy.Fill(robot,q))
#     pause(0.1)

path = [array([ 0.        , -0.0132645 , -0.4       ,  2.00000006, -0.32724929,
       -0.00785398,  0.13962633,  0.        , -0.05235987, -0.17453293,
        0.        ,  0.        , -0.05235987,  1.5       ,  0.        ,
       -0.52359879]),
 array([-0.19518036,  0.00928079, -0.12645753,  1.72075662, -0.36127285,
       -0.07139689,  0.22614146, -0.1304928 , -0.09054683, -0.34062555,
       -0.15434703, -0.54081369, -0.33771027,  1.37411648,  0.18445678,
       -0.66584273]),
 array([-0.29781148,  0.02942216,  0.02580322,  1.54797433, -0.3682557 ,
       -0.09857651,  0.22746339, -0.18234513, -0.24909796, -0.48228371,
       -0.2613325 , -0.81468988, -0.36873148,  1.29637404,  0.26884063,
       -0.71945978]),
 array([-0.44569057,  0.11256224,  0.24797617,  1.2520952 , -0.28837653,
       -0.09814029,  0.2297418 , -0.18519987, -1.45989102, -0.82290031,
       -0.4093418 , -1.01077852, -0.09295965,  1.19336025,  0.42976125,
       -0.71647015]),
 array([-0.37620994,  0.15718295,  0.29049646,  1.23477474, -0.29216782,
       -0.09048463,  0.26796266, -0.01274109, -1.79444697, -0.64532757,
       -0.57732454, -0.95157046, -0.17017959,  1.12245097,  0.58634814,
       -0.74349619]),
 array([-0.33375808,  0.16567347,  0.33596364,  1.22220649, -0.28027406,
       -0.05907819,  0.39743296,  0.28403037, -1.94416004, -0.42390882,
       -0.80131845, -0.94360093, -0.21952529,  1.07478601,  0.73972479,
       -0.79430197]),
 array([-0.37034838,  0.14035753,  0.35181946,  1.22958999, -0.27655137,
       -0.04769716,  0.4713232 ,  0.38430408, -1.77431685, -0.44100782,
       -0.71401707, -0.98714149, -0.2631794 ,  1.07215904,  0.80708118,
       -0.84391803]),
 array([-0.72697809, -0.0132645 ,  0.48348948,  1.08640828, -0.32724943,
       -0.00785401,  1.17184231,  0.78738843, -0.05235964, -1.25348046,
        0.        , -1.43109325, -0.05235964,  1.23972591,  1.50437154,
       -1.35968925])]





############################ TOPP ############################


t0 = time.time()

traj0 =  MotionPlanning.MakeParabolicTrajectory(path,robot,constraintstring,tuningsstring)


TOPPopenravepy.PlotTorques(robot,traj0,traj0,0.01,taumin,taumax,3)
TOPPopenravepy.PlotZMP(robot,traj0,traj0,zmplimits,0.01,4,border=0.1)

raw_input()
