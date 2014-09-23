import time
import string
import scipy
from pylab import *
from numpy import *
import cvxopt
import cvxopt.solvers
from openravepy import *
from TOPP import TOPPpy
from TOPP import Trajectory
from TOPP import Utilities
from TOPP import TOPPbindings
from TOPP import TOPPopenravepy
from cStringIO import StringIO
import StringIO as oldIO

cvxopt.solvers.options['show_progress'] = False


#############################################################
#################### Trajectory utils #######################
#############################################################

# Returns None if collision or singularity
# Calculate the dependent that compensates for a displacement freed
def Compensate(robot,freedofs,dependentdofs,constrainedlinks,dofvalues,freedelta,toljacobian=0.01,slidedofs=None,slidelinks=None,slidedelta=None):
    if slidelinks is None:
        Jfree = zeros((6*len(constrainedlinks),len(freedofs)))
        Jdependent = zeros((6*len(constrainedlinks),len(dependentdofs)))
    else:
        Jfree = zeros((6*len(constrainedlinks)+6*len(slidelinks),len(freedofs)))
        Jdependent = zeros((6*len(constrainedlinks)+6*len(slidelinks),len(dependentdofs)+len(slidedofs)))
    with robot:
        robot.SetDOFValues(dofvalues)
        # if robot.CheckSelfCollision():
        #     return None        
        # Jacobian of the constrainedlinks wrt the free dofs
        for i in range(len(constrainedlinks)):
            p = robot.GetLinks()[constrainedlinks[i]].GetTransform()[0:3,3]
            Jdependent[i*6:i*6+3,:len(dependentdofs)] = robot.ComputeJacobianTranslation(constrainedlinks[i],p,dependentdofs)
            Jdependent[i*6+3:i*6+6,:len(dependentdofs)] = robot.ComputeJacobianAxisAngle(constrainedlinks[i],dependentdofs)
            Jfree[i*6:i*6+3,:] = robot.ComputeJacobianTranslation(constrainedlinks[i],p,freedofs)
            Jfree[i*6+3:i*6+6,:] = robot.ComputeJacobianAxisAngle(constrainedlinks[i],freedofs)
        nstart = 6*len(constrainedlinks)
        if not (slidelinks is None):
            for i in range(len(slidelinks)):
                p = robot.GetLinks()[slidelinks[i]].GetTransform()[0:3,3]
                Jdependent[nstart+i*6:nstart+i*6+3,len(dependentdofs):] = robot.ComputeJacobianTranslation(slidelinks[i],p,slidedofs)
                Jdependent[nstart+i*6+3:nstart+i*6+6,len(dependentdofs):] = robot.ComputeJacobianAxisAngle(slidelinks[i],slidedofs)
                Jfree[nstart+i*6:nstart+i*6+3,:] = robot.ComputeJacobianTranslation(slidelinks[i],p,freedofs)
                Jfree[nstart+i*6+3:nstart+i*6+6,:] = robot.ComputeJacobianAxisAngle(slidelinks[i],freedofs)

    # compute the compensation
    #dependentdelta = -dot(pinv(Jdependent),dot(Jfree,freedelta))
    # Assume that Jpendent is square    
    if abs(det(Jdependent)) < toljacobian:
        print "Singularity detected: det(Jdependent) =", det(Jdependent)
        return None
    if slidelinks is None:
        dependentdelta = linalg.solve(Jdependent,-dot(Jfree,freedelta))
    else:
        v = zeros(6*len(constrainedlinks)+6*len(slidelinks))
        v[6*len(constrainedlinks):] = slidedelta
        dependentdelta = linalg.solve(Jdependent,-dot(Jfree,freedelta)+v)
    return dependentdelta


# Compensate a trajectory
# Returns None if collision or singularity
def CompensateTrajectory(robot,qstart,freedofs,dependentdofs,constrainedlinks,freetraj,dependent0,nchunks,chunksubdiv,toljacobian=0.01,slidedofs=None,slidelinks=None,slidedelta=None):
    # Compensate trajectory
    nsubdiv = nchunks*chunksubdiv
    dt = freetraj.duration/nsubdiv
    dependentvalues = array(dependent0)
    dependentv = []
    dependentv.append(array(dependentvalues))
    tv = [0]
    dofvalues = array(qstart)
    for ichunk in range(nsubdiv):
        t = ichunk * dt
        freevalues = freetraj.Eval(t)
        dofvalues[freedofs] = freevalues
        if slidedofs is None:
            dofvalues[dependentdofs] = dependentvalues
        else:
            dofvalues[dependentdofs+slidedofs] = dependentvalues
        freedelta = freetraj.Evald(t)
        if slidelinks is None:
            dependentdelta = Compensate(robot,freedofs,dependentdofs,constrainedlinks,dofvalues,freedelta,toljacobian)
        else:
            dependentdelta = Compensate(robot,freedofs,dependentdofs,constrainedlinks,dofvalues,freedelta,toljacobian,slidedofs,slidelinks,slidedelta)
        if dependentdelta is None:
            print "t=",t
            return None
        dependentvalues += dependentdelta*dt
        dependentv.append(array(dependentvalues))
        tv.append(t+dt)
    dependentv = array(dependentv)
    tv = array(tv)

    # Interpolate dependent values with splines
    tcklist = []
    for idof in range(0,dependentv.shape[1]):
        tcklist.append(scipy.interpolate.splrep(tv[::chunksubdiv],dependentv[::chunksubdiv,idof],s=0))
    t = tcklist[0][0]
    chunkslist = []
    for i in range(len(t)-1):        
        dependentpolylist = []
        if abs(t[i+1]-t[i])>1e-5:
            for tck in tcklist:
                a = 1/6. * scipy.interpolate.splev(t[i],tck,der=3)
                b = 0.5 * scipy.interpolate.splev(t[i],tck,der=2)
                c = scipy.interpolate.splev(t[i],tck,der=1)
                d = scipy.interpolate.splev(t[i],tck,der=0)
                dependentpolylist.append(Trajectory.Polynomial([d,c,b,a]))
            chunkslist.append(Trajectory.Chunk(t[i+1]-t[i],dependentpolylist))
    return Trajectory.PiecewisePolynomialTrajectory(chunkslist)


# Interpolate freely
def InterpolateFree(qbeg, qend, qsbeg, qsend, duration):
    pathstring = ''
    ndof = len(qbeg)
    pathstring += "%f\n%d"%(duration, ndof)
    for k in range(ndof):
        a,b,c,d = Utilities.Interpolate3rdDegree(qbeg[k], qend[k], qsbeg[k], qsend[k], duration)
        pathstring += "\n%f %f %f %f"%(d, c, b, a)
    return Trajectory.PiecewisePolynomialTrajectory.FromString(pathstring)


# Merge free and dependent trajs
def MergeTrajectories(q0,freetraj,dependenttraj,freedofs,dependentdofs):
    nchunks = len(dependenttraj.chunkslist)
    ti = 0
    chunkslist = []
    polyarray = array([Trajectory.Polynomial([q]) for q in q0])
    for i in range(nchunks):
        freepolylist = []
        dependentchunk = dependenttraj.chunkslist[i]
        freevalues = freetraj.Eval(ti)
        freed = freetraj.Evald(ti)
        freedd = freetraj.Evaldd(ti)
        freeddd = freetraj.Evaldn(ti,3)
        for ifree in range(freetraj.dimension):
            a = 1/6. * freeddd[ifree]
            b = 0.5 * freedd[ifree]
            c = freed[ifree]
            d = freevalues[ifree]
            freepolylist.append(Trajectory.Polynomial([d,c,b,a]))
        polyarray[freedofs] = freepolylist
        polyarray[dependentdofs] = dependentchunk.polynomialsvector
        chunkslist.append(Trajectory.Chunk(dependentchunk.duration,list(polyarray)))
        ti += dependentchunk.duration
    return Trajectory.PiecewisePolynomialTrajectory(chunkslist)

def InsertTrajectory(q0,freetraj,freedofs):
    polyarray = array([Trajectory.Polynomial([q]) for q in q0])
    chunkslist = []
    for chunk in freetraj.chunkslist:
        polyarray[freedofs] = chunk.polynomialsvector
        chunkslist.append(Trajectory.Chunk(chunk.duration,list(polyarray)))
    return Trajectory.PiecewisePolynomialTrajectory(chunkslist)

def SplitTrajectories(trajtotal,freedofs,dependentdofs):
    nchunks = len(trajtotal.chunkslist)
    freechunkslist = []
    dependentchunkslist = []
    for chunk in trajtotal.chunkslist:
        freepolylist = array(chunk.polynomialsvector)[freedofs]
        dependentpolylist = array(chunk.polynomialsvector)[dependentdofs]
        freechunkslist.append(Trajectory.Chunk(chunk.duration,freepolylist))
        dependentchunkslist.append(Trajectory.Chunk(chunk.duration,dependentpolylist))
    return Trajectory.PiecewisePolynomialTrajectory(freechunkslist), Trajectory.PiecewisePolynomialTrajectory(dependentchunkslist)

def ComputeVelocities(robot,q,linkindex,linkactivedofs,linkvelocity):
    with robot:
        robot.SetDOFValues(q)
        p = robot.GetLinks()[linkindex].GetGlobalCOM()
        J = robot.ComputeJacobianTranslation(linkindex,p,linkactivedofs)
        return dot(pinv(J),linkvelocity)

###########################################################
################# Polygon projection ######################
###########################################################

# Optimize in one direction
def OptimizeDirection(vdir,lp):
    lp_q, lp_Gextended, lp_hextended, lp_A, lp_b = lp
    lp_q[-2] = -vdir[0]
    lp_q[-1] = -vdir[1]
    try:
        sol = cvxopt.solvers.lp(lp_q, lp_Gextended, lp_hextended, lp_A, lp_b)
        if sol['status'] == 'optimal':
            z = sol['x']
            z = array(z).reshape((lp_q.size[0],))
            return True, z[-2:]
        else:
            return False, 0
    except Exception as inst:
        print inst
        return False,0
    


# Compute the polygon in the (sddot,sdot^2) plane using Hauser technique
def ComputePolygon(lp):

    t0 = time.time()
    res, z1 = OptimizeDirection(array([1.,0.]),lp)
    if not res:
        return False,0
    res, z2 = OptimizeDirection(array([cos(2*pi/3),sin(2*pi/3)]),lp)
    if not res:
        return False,0    
    res, z3 = OptimizeDirection(array([cos(4*pi/3),sin(4*pi/3)]),lp)
    if not res:
        return False,0
    v1 = Vertex(z1)
    v2 = Vertex(z2)
    v3 = Vertex(z3)
    P0 = Polygon()
    P0.fromVertices(v1,v2,v3)
    P0.iter_expand(lp,100)

    return True, P0



##############################################################
##################### Polygon classes ########################
##############################################################

class Vertex:
    def __init__(self,z):
        self.x = z[0]
        self.y = z[1]
        self.next = None
        self.expanded = False


    def length(self):
        return linalg.norm([self.x-self.next.x,self.y-self.next.y])        
    def expand(self,lp):
        v1 = self
        v2 = self.next
        v = array([v2.y-v1.y,v1.x-v2.x]) #orthogonal direction to edge
        v = 1/linalg.norm(v) * v
        res, z = OptimizeDirection(v,lp)
        if not res:
            self.expanded = True
            return False, None
        xopt, yopt = z
        if(abs(cross([xopt-v1.x,yopt-v1.y],[v1.x-v2.x,v1.y-v2.y]))<1e-2):
            self.expanded = True
            return False, None
        else:
            vnew = Vertex([xopt,yopt])
            vnew.next = self.next
            self.next = vnew
            self.expanded = False
            return True, vnew

    def Plot(self):
        plot([self.x,self.next.x],[self.y,self.next.y])
    
    def Print(self):
        print self.x,self.y,"to",self.next.x, self.next.y

class Polygon:
    def __init__(self):
        pass

    def fromVertices(self,v1,v2,v3):
        v1.next = v2
        v2.next = v3
        v3.next = v1
        self.vertices = [v1,v2,v3]

    def fromString(self,s):
        buff = oldIO.StringIO(s)
        self.vertices = []
        while(True):
            l = buff.readline()
            l = l.strip(" \n")
            if len(l) < 2:
                break
            x,y = [double(x) for x in l.split(' ')]
            vnew = Vertex([x,y])
            self.vertices.append(vnew)

        for i in range(len(self.vertices)-1):            
            self.vertices[i].next = self.vertices[i+1]            
        self.vertices[-1].next = self.vertices[0]
        
    def all_expanded(self):
        for v in self.vertices:
            if not v.expanded:
                return False
        return True

    # Returns true if there's a edge that can be expanded, and expands that edge, otherwise returns false
    def iter_expand(self,qpconstraints,maxiter = 10):
        niter = 0
        v = self.vertices[0]
        while not self.all_expanded() and niter <maxiter:
            if not v.expanded:
                res, vnew = v.expand(qpconstraints)
                if res:
                    self.vertices.append(vnew)
                    niter += 1
            else:
                v = v.next
    
    # Assumes every vertices are on the positive halfplane
    # Export the vertices starting from the left-most and going clockwise
    def sort_vertices(self):
        minsd = 1e10
        ibottom = 0
        for i in range(len(self.vertices)):
            v = self.vertices[i]
            if (v.y + v.next.y) < minsd:
                ibottom = i
                minsd = v.y + v.next.y
        for v in self.vertices:
            v.checked = False
        vcur = self.vertices[ibottom]
        newvertices = []
        while not vcur.checked:
            vcur.checked = True
            newvertices.append(vcur)
            vcur = vcur.next
        newvertices.reverse()
        vfirst = newvertices.pop(-1)
        newvertices.insert(0,vfirst)
        self.vertices = newvertices
        

    def export_vertices(self,threshold=1e-2):
        export_list = [self.vertices[0]]
        for i in range(1,len(self.vertices)-1):
            vcur = self.vertices[i]
            vlast = export_list[-1]
            if(linalg.norm([vcur.x-vlast.x,vcur.y-vlast.y]))>threshold:
                export_list.append(vcur)
        # Always add last vertex
        export_list.append(self.vertices[-1])
        return export_list
        
    # Plot
    def Plot(self):
        hold("on")
        for v in self.vertices:
            v.Plot()

    def Print(self):
        print "Polygon contains vertices"
        for v in self.vertices:
            # print"  "
            v.Print()
