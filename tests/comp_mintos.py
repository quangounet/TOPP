import sys
import os
sys.path.append('..')

import TOPPbindings
import TOPPpy
import time
import string
from pylab import *
from numpy import *

ion()

random.seed(0)

def vector2string(v):
    ndof = len(v);
    s = str(ndof)
    for a in v:
        s+= ' %f'%a
    return s


def log10(x):
    return log(x)/log(10)


ncurve = 10
ndofs = [2,5,10,20,50,100]
lowcoef = 0.2
vn = [[1.2/5,1.2,1.2*5]]*3 + [[1.5/5,1.5,1.5*5]]*3
an = [[1]]*7


# for ndof in ndofs:
#     for j in range(30):
#         p0a = TOPPpy.vector2string(rand(ndof)*2*pi-pi)
#         p0b = TOPPpy.vector2string(rand(ndof)*2*pi-pi)
#         p1a = TOPPpy.vector2string(rand(ndof)*2*pi-pi)
#         p1b = TOPPpy.vector2string(rand(ndof)*2*pi-pi)
#         s = '%d'%ncurve
#         s+= '\n1.0 ' + p0a + ' ' + p0b
#         for k in range(ncurve-1):    
#             a = rand(ndof)*2*pi-pi
#             b = rand(ndof)*2*pi-pi
#             c = 2*b-a
#             pa = TOPPpy.vector2string(a)
#             pb = TOPPpy.vector2string(b)
#             pc = TOPPpy.vector2string(c)
#             s+= ' ' + pa + ' ' + pb + '\n1.0 ' + pb + ' ' + pc
#         s+= ' ' + p1a + ' ' + p1b
#         filename = 'testfiles/traj-%d-%d'%(ndof,j)
#         handle = open(filename,'w')
#         handle.write(s)
#         handle.close()





gridres = 1000
gridres2 = 1000
discrtimestep = ncurve/gridres
integrationtimestep = 0 # auto
reparamtimestep = 0 # auto
passswitchpointnsteps = 0 #auto
tuningsstring = "%f %f %f %d"%(discrtimestep,integrationtimestep,reparamtimestep,passswitchpointnsteps)

time_low = []
time_low_mintos = []
time_lows = []
time_low_mintoss = []
diff_low = []
time_normal = []
time_normal_mintos = []
time_normals = []
time_normal_mintoss = []
diff_normal = []
time_high = []
time_high_mintos = []
time_highs = []
time_high_mintoss = []
diff_high = []


for i in range(len(ndofs)):
    ndof = ndofs[i]
    for v in vn[i]: 
        for a in an[i]:
            #print v,a
            limitfile = 'testfiles/limits-%d-%f-%f'%(ndof,v,a)
            comput_time_v=[]
            traj_dur_v=[]
            comput_time_mintos_v=[]
            traj_dur_mintos_v=[]
            for j in range(30):
                #print j
                trajfile = 'testfiles/traj-%d-%d'%(ndof,j)
                h = open(trajfile,'r')
                s = h.read()
                h.close()      
                Tv,p0v,p1v,p2v,p3v = TOPPpy.string2p(s)
                trajectorystring = TOPPpy.BezierToTrajectoryString(Tv,p0v,p1v,p2v,p3v)
                trajectorystring = TOPPpy.BezierToTrajectoryString(Tv,p0v,p1v,p2v,p3v)
                vmax = v*ones(ndof)
                amax = a*ones(ndof)
                constraintstring = string.join([str(v) for v in amax]) + "\n"
                constraintstring += string.join([str(v) for v in vmax])
                t0 = time.time()
                x = TOPPbindings.TOPPInstance("KinematicLimits",constraintstring,trajectorystring,tuningsstring,False);
                ret = x.RunComputeProfiles(0,0)
                # if(ret == 1):
                #     x.ReparameterizeTrajectory()
                comput_time = time.time() - t0
                # if(ret == 1):
                #     x.WriteResultTrajectory()
                #     traj1 = TOPPpy.PiecewisePolynomialTrajectory.FromString(x.restrajectorystring)                               
                # print 'Ndof: ', ndof, ' ; Vmax: ', vmax , '; Amax: ', amax, '; Traj num: ', j
                # print trajectorystring
                # print 'Computation time: ',comput_time
                comput_time_v.append(comput_time)
                if(ret == 1):
                #     print 'Trajectory duration: ', x.resduration
                    traj_dur_v.append(x.resduration)
                # else:
                #     print '************* Error **************'
                
                command = "./timeopt %s %s %d 1 > /tmp/res-%d-%f"%(trajfile,limitfile,gridres2,ndof,v)
                os.system(command)
                res = open("/tmp/res-%d-%f"%(ndof,v),"r").read()
                lines = [l.strip(" \n") for l in res.split('\n')]
                resline = ""
                for l in lines:
                    if(len(l)>0 and l[0]=="O"):
                        resline = l
                        break
                traj_dur_mintos_v.append(float(resline.split(" ")[2].strip(",")))
                comput_time_mintos_v.append(float(resline.split(" ")[4]))
                


            comput_time_mean = mean(log10(comput_time_v))
            comput_time_std = std(log10(comput_time_v))
            traj_dur_mean = mean(traj_dur_v)
            traj_dur_std = std(traj_dur_v)
            comput_time_mintos_mean = mean(log10(comput_time_mintos_v))
            comput_time_mintos_std = std(log10(comput_time_mintos_v))
            traj_dur_mintos_mean = mean(traj_dur_mintos_v)
            traj_dur_mintos_std = std(traj_dur_mintos_v)            
            print '\nNdof: ', ndof, ' ; Vmax: ', v , '; Amax: ', a;
            print "Computation time       : ", comput_time_mean, " +- ", comput_time_std
            print "Computation time mintos: ", comput_time_mintos_mean, " +- ", comput_time_mintos_std
            print "Traj duration       : ", traj_dur_mean, " +- ", traj_dur_std

            print "Traj duration mintos: ", traj_dur_mintos_mean, " +- ", traj_dur_mintos_std
            print "Difference: ", abs(traj_dur_mean-traj_dur_mintos_mean)/traj_dur_mean*100, "%"
            if v<a/2.:
                time_low.append(comput_time_mean)
                time_low_mintos.append(comput_time_mintos_mean)
                time_lows.append(comput_time_std)
                time_low_mintoss.append(comput_time_mintos_std)
                diff_low.append(abs(traj_dur_mean-traj_dur_mintos_mean)/traj_dur_mean)
            elif v>a*2:
                time_high.append(comput_time_mean)
                time_high_mintos.append(comput_time_mintos_mean)
                time_highs.append(comput_time_std)
                time_high_mintoss.append(comput_time_mintos_std)
                diff_high.append(abs(traj_dur_mean-traj_dur_mintos_mean)/traj_dur_mean)
            else:
                time_normal.append(comput_time_mean)
                time_normal_mintos.append(comput_time_mintos_mean)
                time_normals.append(comput_time_std)
                time_normal_mintoss.append(comput_time_mintos_std)
                diff_normal.append(abs(traj_dur_mean-traj_dur_mintos_mean)/traj_dur_mean)


X = log10(ndofs)

figure(0)
clf()
errorbar(X,time_low, yerr=time_lows, linewidth=3, fmt='r')
errorbar(X,time_low_mintos, yerr=time_low_mintoss, linewidth=2, fmt='r--')
errorbar(X,time_normal, yerr=time_normals, linewidth=3, fmt='g')
errorbar(X,time_normal_mintos, yerr=time_normal_mintoss, linewidth=2, fmt='g--')
errorbar(X,time_high, yerr=time_highs, linewidth=3, fmt='b')
errorbar(X,time_high_mintos, yerr=time_high_mintoss, linewidth=2, fmt='b--')
axis([0,2.2,-3.5,0])
title('Computation time',fontsize=20)
xlabel('Dof',fontsize=18)
ylabel('Computation time in seconds (log)',fontsize=18)    
ax=gca()
ax.set_xticks(X)
ax.set_xticklabels([str(x) for x in ndofs])
grid('on')

raw_input()
