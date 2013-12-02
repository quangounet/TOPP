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


integrationtimestep = 0 # auto
reparamtimestep = 0 # auto
passswitchpointnsteps = 0 #auto


time_low = []
time_low_mintos = []
time_lows = []
time_low_mintoss = []
diff_low = []

all_dur = []
all_dur_mintos = []

v = 1.2
a = 1
ndof = 10
gridresv = [100,200,300,400,500,600,700,800,900,1000]


for i in range(len(gridresv)):
    gridres = gridresv[i]
    print gridres
    discrtimestep = 1./gridres
    tuningsstring = "%f %f %f %d"%(discrtimestep,integrationtimestep,reparamtimestep,passswitchpointnsteps)
    limitfile = 'testfiles/limits-%d-%f-%f'%(ndof,v,a)
    comput_time_v=[]
    traj_dur_v=[]
    comput_time_mintos_v=[]
    traj_dur_mintos_v=[]
    for j in range(30):
        trajfile = 'testfiles/traj-%d-%d'%(ndof,j)
        h = open(trajfile,'r')
        s = h.read()
        h.close()      
        [p0v,p1v,p2v,p3v] = TOPPpy.string2p(s)
        Tv = [1]
        trajectorystring = TOPPpy.BezierToTrajectoryString(Tv,p0v,p1v,p2v,p3v)
        vmax = v*ones(ndof)
        amax = a*ones(ndof)
        constraintstring = string.join([str(v) for v in amax]) + "\n"
        constraintstring += string.join([str(v) for v in vmax])
        t0 = time.time()
        x = TOPPbindings.TOPPInstance("KinematicLimits",constraintstring,trajectorystring,tuningsstring,False);
        ret = x.RunComputeProfiles(0,0)
        comput_time = time.time() - t0
        comput_time_v.append(comput_time)
        if(ret == 1):
            traj_dur = x.resduration
        command = "./timeopt %s %s %d 1 > /tmp/res-%d-%f"%(trajfile,limitfile,gridres,ndof,v)
        traj_dur_mintos = 1
        os.system(command)
        res = open("/tmp/res-%d-%f"%(ndof,v),"r").read()
        lines = [l.strip(" \n") for l in res.split('\n')]
        resline = ""
        for l in lines:
            if(len(l)>0 and l[0]=="O"):
                resline = l
                break        
        traj_dur_mintos = float(resline.split(" ")[2].strip(","))
        comput_time_mintos_v.append(float(resline.split(" ")[4]))
        traj_dur_v.append(traj_dur)
        traj_dur_mintos_v.append(traj_dur_mintos)

    all_dur.append(traj_dur_v)
    all_dur_mintos.append(traj_dur_mintos_v)
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
    time_low.append(comput_time_mean)
    time_low_mintos.append(comput_time_mintos_mean)
    time_lows.append(comput_time_std)
    time_low_mintoss.append(comput_time_mintos_std)
    diff_low.append(abs(traj_dur_mean-traj_dur_mintos_mean)/traj_dur_mean)


X = gridresv

figure(0)
clf()
for j in range(30):
    plot(X, [(l[j]-all_dur[-1][j])/all_dur[-1][j]*100 for l in all_dur], 'g', linewidth=1)
    plot(X, [(l[j]-all_dur_mintos[-1][j])/all_dur_mintos[-1][j]*100 for l in all_dur_mintos], 'r' , linewidth=1)

xlabel('Grid res',fontsize=18)
ylabel('Computation time in seconds (log)',fontsize=18)    
ax=gca()
ax.set_xticks(X)
ax.set_xticklabels([str(x) for x in gridresv])


for j in range(30):
    l = all_dur[6]
    print j,abs((l[j]-all_dur[-1][j])/all_dur[-1][j]*100)




raw_input()
