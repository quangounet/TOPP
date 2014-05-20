import string,time,pickle,os
from pylab import *
from numpy import *
from openravepy import *
from TOPP import TOPPbindings
from TOPP import TOPPpy
from TOPP import TOPPopenravepy
from TOPP import Trajectory
from TOPP import Utilities

print "\n**********************************\nNB: This test file requires MINTOS\n**********************************\n"


ion()

ndof = 7

# Generate random trajectories
ntraj = 1000
ncurve = 1
for j in range(ntraj):
    p0a = Utilities.vect2str_mintos(rand(ndof)*2*pi-pi)
    p0b = Utilities.vect2str_mintos(rand(ndof)*2*pi-pi)
    p1a = Utilities.vect2str_mintos(rand(ndof)*2*pi-pi)
    p1b = Utilities.vect2str_mintos(rand(ndof)*2*pi-pi)
    s = '%d'%ncurve
    s+= '\n1.0 ' + p0a + ' ' + p0b
    for k in range(ncurve-1):    
        a = rand(ndof)*2*pi-pi
        b = rand(ndof)*2*pi-pi
        c = 2*b-a
        pa = Utilities.vect2str(a)
        pb = Utilities.vect2str(b)
        pc = Utilities.vect2str(c)
        s+= ' ' + pa + ' ' + pb + '\n1.0 ' + pb + ' ' + pc
    s+= ' ' + p1a + ' ' + p1b
    filename = 'traj-%d-%d'%(ndof,j)
    handle = open(filename,'w')
    handle.write(s)
    handle.close()


# Generate velocity and acceleration bounds files (for Mintos)    
vmax0 = 4
amax0 = 20
st = ' -%f'%vmax0
s = str(ndof) + st*ndof + '\n';
st = ' %f'%vmax0            
s += str(ndof) + st*ndof + '\n';
st = ' -%f'%amax0
s += str(ndof) + st*ndof + '\n';
st = ' %f'%amax0
s += str(ndof) + st*ndof;
handle = open("limit",'w')
handle.write(s)
handle.close()


# Test starts here
gridres = 300
discrtimestep = 1./gridres
comput_time_v = []
traj_dur_v = []
comput_time_mintos_v = []
traj_dur_mintos_v = []
nfailures = 0
for j in range(ntraj):
    print j
    # Load trajectory
    trajfile = 'traj-%d-%d'%(ndof,j)
    h = open(trajfile,'r')
    s = h.read()
    h.close()      
    Tv,p0v,p1v,p2v,p3v = TOPPpy.string2p(s)
    trajectorystring = TOPPpy.BezierToTrajectoryString(Tv,p0v,p1v,p2v,p3v)
    # TOPP
    vmax = vmax0*ones(ndof)
    amax = amax0*ones(ndof)
    constraintstring = str(discrtimestep)
    constraintstring += "\n" + string.join([str(v) for v in vmax])
    constraintstring += "\n" + string.join([str(a) for a in amax])
    t0 = time.time()
    x = TOPPbindings.TOPPInstance(None,"KinematicLimits",constraintstring,trajectorystring)
    x.integrationtimestep = discrtimestep
    x.extrareps = 5
    ret = x.RunComputeProfiles(0,0)
    comput_time = time.time() - t0
    comput_time_v.append(comput_time)
    if(ret == 1):
        traj_dur = x.resduration
        traj_dur_v.append(traj_dur)
    else:
        nfailures += 1
        print ">>>>>>>>>>>>>>>>>>> TOPP could not retime: ", [ndof, j, v, a, discrtimestep]
    # MINTOS
    command = "./timeopt %s %s %d 1 > res"%(trajfile,"limit",gridres)
    os.system(command)
    res = open("res","r").read()
    lines = [l.strip(" \n") for l in res.split('\n')]
    resline = ""
    for l in lines:
        if(len(l)>0 and l[0]=="O"):
            resline = l
            break        
    traj_dur_mintos = float(resline.split(" ")[2].strip(","))
    comput_time_mintos_v.append(float(resline.split(" ")[4]))
    traj_dur_mintos_v.append(traj_dur_mintos)

comput_time_mean = mean(comput_time_v)
comput_time_std = std(comput_time_v)
traj_dur_mean = mean(traj_dur_v)
traj_dur_std = std(traj_dur_v)
comput_time_mintos_mean = mean(comput_time_mintos_v)
comput_time_mintos_std = std(comput_time_mintos_v)
traj_dur_mintos_mean = mean(traj_dur_mintos_v)
traj_dur_mintos_std = std(traj_dur_mintos_v)            
print '\nNdof: ', ndof, ' ; Gridres: ', gridres
print "Number of TOPP failures:", nfailures
print "Computation time TOPP:", comput_time_mean, " +- ", comput_time_std
print "Computation time MINTOS:", comput_time_mintos_mean, " +- ", comput_time_mintos_std
print "Ratio:", comput_time_mintos_mean/comput_time_mean
print "Traj duration TOPP:", traj_dur_mean, " +- ", traj_dur_std
print "Traj duration MINTOS:", traj_dur_mintos_mean, " +- ", traj_dur_mintos_std
print "Difference:", abs(traj_dur_mean-traj_dur_mintos_mean)/traj_dur_mean*100, "%"
