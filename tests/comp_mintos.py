import sys
import os
sys.path.append('..')

import TOPPbindings
import TOPPpy
import time
import string
from pylab import *
from numpy import *



random.seed(0)

def vector2string(v):
    ndof = len(v);
    s = str(ndof)
    for a in v:
        s+= ' %f'%a
    return s


# for ndof in [2,5,10,20,50,100]:
#     for j in range(30):
#         p0s = vector2string(rand(ndof)*2*pi-pi)
#         p1s = vector2string(rand(ndof)*2*pi-pi)
#         p2s = vector2string(rand(ndof)*2*pi-pi)
#         p3s = vector2string(rand(ndof)*2*pi-pi)
#         s = '1\n1.0 ' + p0s + ' ' + p1s + ' ' + p2s + ' ' + p3s
#         filename = 'testfiles/traj-%d-%d'%(ndof,j)
#         handle = open(filename,'w')
#         handle.write(s)
#         handle.close()



ndofs = [2,5,10,20,50,100]
lowcoef = 0.2
vn = [[lowcoef,1]]*len(ndofs)
an = [[lowcoef*a,a] for a in [1,0.6,0.55,0.5,0.45,0.4]]

# for i in range(len(ndofs)):
#     ndof = ndofs[i]
#     for v in vn[i]: 
#         for a in an[i]:
#             filename = 'testfiles/limits-%d-%f-%f'%(ndof,v,a)
#             st = ' -%f'%v
#             s = str(ndof) + st*ndof + '\n';
#             st = ' %f'%v            
#             s += str(ndof) + st*ndof + '\n';
#             st = ' -%f'%a
#             s += str(ndof) + st*ndof + '\n';
#             st = ' %f'%a
#             s += str(ndof) + st*ndof;
#             handle = open(filename,'w')
#             handle.write(s)
#             handle.close()


gridres = 500
discrtimestep = 1./gridres
integrationtimestep = 0 # auto
reparamtimestep = 0 # auto
passswitchpointnsteps = 5
tuningsstring = "%f %f %f %d"%(discrtimestep,integrationtimestep,reparamtimestep,passswitchpointnsteps)



for i in range(len(ndofs)):
    ndof = ndofs[i]
    for v in vn[i]: 
        for a in an[i]:
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
                x = TOPPbindings.TOPPInstance("KinematicLimits",constraintstring,trajectorystring,tuningsstring);
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
                
                command = "./timeopt %s %s %d 1 > /tmp/res"%(trajfile,limitfile,gridres)
                os.system(command)
                res = open("/tmp/res","r").read()
                lines = [l.strip(" \n") for l in res.split('\n')]
                resline = ""
                for l in lines:
                    if(len(l)>0 and l[0]=="O"):
                        resline = l
                        break
                traj_dur_mintos_v.append(float(resline.split(" ")[2].strip(",")))
                comput_time_mintos_v.append(float(resline.split(" ")[4]))
                


            comput_time_mean = mean(comput_time_v)
            comput_time_std = std(comput_time_v)
            traj_dur_mean = mean(traj_dur_v)
            traj_dur_std = std(traj_dur_v)
            comput_time_mintos_mean = mean(comput_time_mintos_v)
            comput_time_mintos_std = std(comput_time_mintos_v)
            traj_dur_mintos_mean = mean(traj_dur_mintos_v)
            traj_dur_mintos_std = std(traj_dur_mintos_v)
            print '\nNdof: ', ndof, ' ; Vmax: ', v , '; Amax: ', a;
            print "Computation time       : ", comput_time_mean, " +- ", comput_time_std
            print "Computation time mintos: ", comput_time_mintos_mean, " +- ", comput_time_mintos_std
            print "Traj duration       : ", traj_dur_mean, " +- ", traj_dur_std

            print "Traj duration mintos: ", traj_dur_mintos_mean, " +- ", traj_dur_mintos_std

        
