import scipy
from pylab import *
from numpy import *
import cvxopt
import cvxopt.solvers
from scipy import interpolate
from openravepy import *
from TOPP import TOPPpy
from TOPP import Trajectory
from TOPP import Utilities

cvxopt.solvers.options['show_progress'] = False
import ClosedChain


#############################################################
######################## Kinematics #########################
#############################################################

def Getxztheta(T):
    x = T[0,3]
    z = T[2,3]
    theta = pi -arctan2(T[2,0],T[0,0])
    return array([x,z,theta])

def Getxytheta(T):
    x = T[0,3]
    y = T[1,3]
    theta = arctan2(T[1,0],T[0,0])
    return array([x,y,theta])


def ObjFunc(q,robot,desiredpose,weights):
    l = 9
    with robot:
        robot.SetDOFValues(q)
        x,y,theta = Getxztheta(robot.GetLinks()[l].GetTransform())
    cost = 0
    cost += weights[0] * (x - desiredpose[0])**2
    cost += weights[1] * (y - desiredpose[1])**2
    cost += weights[2] * (theta - desiredpose[2])**2
    return cost

def IK3(robot,desiredpose,q_start=[0.5,0.5,0.5],step=1e-2,threshold=1e-3,k=0.1):
    J = zeros((3,3))
    l = 9
    obj = robot.GetLinks()[l]
    q = q_start
    #with robot:
    while True:
        robot.SetDOFValues(q)
        pose = Getxztheta(obj.GetTransform())
        #print pose, desiredpose
        dpose = desiredpose - pose
        normdpose = linalg.norm(dpose)
        #print normdpose
        if normdpose <= threshold:
            break
        p = obj.GetTransform()[0:3,3]
        J[0:2,:] = robot.ComputeJacobianTranslation(l,p)[[0,2],:]
        J[2,:] = robot.ComputeJacobianAxisAngle(l)[1,:]
        Jstar = dot(linalg.inv(dot(J.T,J) + k*eye(3)),J.T)
        dposeunit = dpose / normdpose
        dposestep = dposeunit * min(step,normdpose)
        #print J, dposestep        
        dq = dot(Jstar,dposestep)
        q += dq
        #raw_input()
    return q, linalg.norm(dpose)


def ComputeDOFVelocities(robot,q,objectvel):
    robot1 = robot.robot1
    l = robot.constrainedlink
    J = zeros((3,3))
    with robot1:
        robot1.SetTransform(robot.T1)
        robot1.SetDOFValues(q)
        J[0:2,:] = robot1.ComputeJacobianTranslation(l,robot1.GetLinks()[l].GetGlobalCOM())[[0,2],:]
        J[2,:] = robot1.ComputeJacobianAxisAngle(l)[1,:]        
    return dot(pinv(J),objectvel)
        



# Calculate the dependent delta that compensates for a displacement free delta
def Compensate(robot,T1,T2,dofvalues1,dofvalues2,freedofs1,freedofs2,dependentdofs1,dependentdofs2,constrainedlink,freedelta,tol_jacobian):
    Jfree = zeros((3,3))
    Jdependent = zeros((3,3))
    # robot in configuration 1
    with robot:
        robot.SetTransform(T1)
        robot.SetDOFValues(dofvalues1)
        p = robot.GetLinks()[constrainedlink].GetTransform()[0:3,3]
        if(len(freedofs1)>0):            
            Jfree[0:2,range(len(freedofs1))] = robot.ComputeJacobianTranslation(constrainedlink,p)[[0,2],:][:,freedofs1]
            Jfree[2,range(len(freedofs1))] = robot.ComputeJacobianAxisAngle(constrainedlink)[1,freedofs1]
        if(len(dependentdofs1)>0):
            Jdependent[0:2,range(len(dependentdofs1))] = robot.ComputeJacobianTranslation(constrainedlink,p)[[0,2],:][:,dependentdofs1]
            Jdependent[2,range(len(dependentdofs1))] = robot.ComputeJacobianAxisAngle(constrainedlink)[1,dependentdofs1]
    # robot in configuration 2
    with robot:
        robot.SetTransform(T2)
        robot.SetDOFValues(dofvalues2)
        p = robot.GetLinks()[constrainedlink].GetTransform()[0:3,3]
        if(len(freedofs2)>0):
            Jfree[0:2,len(freedofs1)+array(range(len(freedofs2)))] = robot.ComputeJacobianTranslation(constrainedlink,p)[[0,2],:][:,freedofs2]
            Jfree[2,len(freedofs1)+array(range(len(freedofs2)))] = robot.ComputeJacobianAxisAngle(constrainedlink)[1,freedofs2] 
        if(len(dependentdofs2)>0):
            Jdependent[0:2,len(dependentdofs1)+array(range(len(dependentdofs2)))] = robot.ComputeJacobianTranslation(constrainedlink,p)[[0,2],:][:,dependentdofs2] 
            Jdependent[2,len(dependentdofs1)+array(range(len(dependentdofs2)))] = robot.ComputeJacobianAxisAngle(constrainedlink)[1,dependentdofs2] 
    # compute the compensation
    if(abs(det(Jdependent))<tol_jacobian):
        raise Exception("Jdependent non invertible")
    dependentdelta = linalg.solve(Jdependent,dot(Jfree,freedelta))
    return dependentdelta

# Compensate a trajectory
def CompensateTrajectory(robot,T1,T2,freedofs1,freedofs2,dependentdofs1,dependentdofs2,constrainedlink,freetraj,dependent0,nchunks,chunksubdiv, tol_jacobian):
    # Compensate trajectory
    nsubdiv = nchunks*chunksubdiv
    dt = freetraj.duration/nsubdiv
    dependentvalues = array(dependent0)
    dependentv = []
    dependentv.append(array(dependentvalues))
    tv = [0]
    dofvalues1 = zeros(robot.GetDOF())
    dofvalues2 = zeros(robot.GetDOF())
    for ichunk in range(nsubdiv):
        t = ichunk * dt
        freevalues = freetraj.Eval(t)
        freed = freetraj.Evald(t)
        if len(freedofs1)>0:
            dofvalues1[freedofs1] = freevalues[range(len(freedofs1))]
        if len(dependentdofs1)>0:
            dofvalues1[dependentdofs1] = dependentvalues[range(len(dependentdofs1))]
        if len(freedofs2)>0:        
            dofvalues2[freedofs2] = freevalues[len(freedofs1) + array(range(len(freedofs2)))]
        if len(dependentdofs2)>0:
            dofvalues2[dependentdofs2] = dependentvalues[len(dependentdofs1) + array(range(len(dependentdofs2)))]
        try:
            dependentd = Compensate(robot,T1,T2,dofvalues1,dofvalues2,freedofs1,freedofs2,dependentdofs1,dependentdofs2,constrainedlink,freed,tol_jacobian)
        except Exception as inst:
            print "Could not interpolate : Compensate failed"
            return None
        dependentvalues += dependentd*dt
        dependentv.append(array(dependentvalues))
        tv.append(t+dt)
    dependentv = array(dependentv)
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


# Interpolate with compensation
# qstartfull : both free and dependent dofs
# qgoal : only the free dofs of the goal config
def Interpolate(robot,tunings,qstartfull,qgoal,freestartvelocities,freegoalvelocities):
    # Extract params from robot and tunings
    robot1 = robot.robot1
    T1 = robot.T1
    T2 = robot.T2
    Gdofs = robot.Gdofs
    Sdofs = robot.Sdofs
    freedofs = robot.freedofs
    dependentdofs = robot.dependentdofs
    actuated = robot.actuated
    constrainedlink = robot.constrainedlink
    tol_jacobian = tunings.tol_jacobian

    duration = tunings.duration
    nchunks = tunings.nchunks
    chunksubdiv = tunings.chunksubdiv
    freetraj = ClosedChain.InterpolateFree(qstartfull[freedofs],qgoal,freestartvelocities,freegoalvelocities,duration)
    dependent0 = qstartfull[dependentdofs]
    freedofs1 = filter(lambda x:x<3,freedofs)
    freedofs2 = filter(lambda x:x>=3,freedofs)
    freedofs2 = [x-3 for x in freedofs2]
    dependentdofs1 = filter(lambda x:x<3,dependentdofs)
    dependentdofs2 = filter(lambda x:x>=3,dependentdofs)
    dependentdofs2 = [x-3 for x in dependentdofs2]

    dependenttraj = CompensateTrajectory(robot1,T1,T2,freedofs1,freedofs2,dependentdofs1,dependentdofs2,constrainedlink,freetraj,dependent0,nchunks,chunksubdiv,tol_jacobian)
    if dependenttraj is None:
        return None
    else:
        trajtotal = ClosedChain.MergeTrajectories(qstartfull,freetraj,dependenttraj,freedofs,dependentdofs)
        return trajtotal
    


###########################################################
######################## Dynamics #########################
###########################################################

# Common function to compute the inequality constraint matrix
def ComputeInequalityConstraintMatrices(taumin,taumax,extended = False):
    ndof = len(taumin)
    if extended:
        G = zeros((2*ndof+1,ndof+2))
        h = zeros(2*ndof+1)
    else:
        G = zeros((2*ndof,ndof))
        h = zeros(2*ndof)    
    G[0:6,0:6] = eye(6)
    G[6:12,0:6] = -eye(6)
    h[0:ndof] = taumax
    h[ndof:2*ndof] = -taumin        
    if extended:
        G[2*ndof,ndof+1] = -1
    return cvxopt.matrix(G), cvxopt.matrix(h)


# Compute the matrices W and S as in Nakamura Yamane
def ComputeSensitivityMatrices(robot1,T1,T2,q,Gdofs,Sdofs,actuated,constrainedlink):
    # actuated : 6-list of booleans indicating whether a joint is actuated or not
    JA = zeros((3,3))
    JB = zeros((3,3))
    JCm = zeros((3,6))
    with robot1:
        # Left arm
        robot1.SetTransform(T1)
        robot1.SetDOFValues(q[0:3])
        p = robot1.GetLinks()[constrainedlink].GetTransform()[0:3,3]
        JA[0:2,:] = robot1.ComputeJacobianTranslation(constrainedlink,p)[[0,2],:]
        JA[2,:] = robot1.ComputeJacobianAxisAngle(constrainedlink)[1,:]
        # Right arm
        robot1.SetTransform(T2)
        robot1.SetDOFValues(q[3:6])
        p = robot1.GetLinks()[constrainedlink].GetTransform()[0:3,3]
        JB[0:2,:] = robot1.ComputeJacobianTranslation(constrainedlink,p)[[0,2],:]
        JB[2,:] = robot1.ComputeJacobianAxisAngle(constrainedlink)[1,:]
        # Compute sensitivity matrices
        JCm[0:3,0:3] = JA
        JCm[0:3,3:6] = -JB
        JG = JCm[0:3,Gdofs]
        JS = JCm[0:3,Sdofs]
        if abs(det(JS))<1e-4 :
            raise Exception("JS non invertible")
        H = -dot(inv(JS),JG)
        W = zeros((5,3))
        for i in range(5):
            indiceG = [j for j,x in enumerate(Gdofs) if x == i]
            indiceS = [j for j,x in enumerate(Sdofs) if x == i]
            if(len(indiceG)>0):
                rowz = zeros(3)
                rowz[indiceG[0]] = 1
                W[i,:] = rowz
            else:
                W[i,:] = H[indiceS[0],:]
        nactuated = len(filter(lambda x:x, actuated))
        S = zeros((nactuated,3))
        for i in range(nactuated):
            indiceG = [j for j,x in enumerate(Gdofs) if x == i]
            indiceS = [j for j,x in enumerate(Sdofs) if x == i]
            if(len(indiceG)>0):
                rowz = zeros(3)
                rowz[indiceG[0]] = 1
                S[i,:] = rowz
            else:
                S[i,:] = H[indiceS[0],:]
    return W, S

# Compute the generalized torques tauG
def ComputeGeneralizedTorques(robot1,robot2,T1,T2,q,qd,qdd,Gdofs,Sdofs,actuated,constrainedlink):
    tauO = zeros(5)
    # Inverse dynamics on the open chains
    with robot1: # left arm, has 3 joints
        robot1.SetTransform(T1)
        robot1.SetDOFValues(q[0:3])
        robot1.SetDOFVelocities(qd[0:3])
        tauO[0:3] = robot1.ComputeInverseDynamics(qdd[0:3])
    with robot2: # right arm, has 2 joints
        robot2.SetTransform(T2)
        robot2.SetDOFValues(q[3:5])
        robot2.SetDOFVelocities(qd[3:5])
        tauO[3:5] = robot2.ComputeInverseDynamics(qdd[3:5])
    # Sensitivity matrices
    W, S = ComputeSensitivityMatrices(robot1,T1,T2,q,Gdofs,Sdofs,actuated,constrainedlink)
    tauG = dot(W.T,tauO)
    return tauG, W.T, S.T


# Find the optimal tauA such that ST*tauA = tauG
def OptimizeTorques(ST,tauG,ineqmatrices):
    qp_P = cvxopt.matrix(eye(6))
    qp_q = cvxopt.matrix(zeros(6))
    #qp_q = cvxopt.matrix(ones(6))
    qp_A = cvxopt.matrix(ST)
    qp_b = cvxopt.matrix(tauG)
    qp_G, qp_h = ineqmatrices
    try:
        #sol = cvxopt.solvers.lp(qp_q, qp_G, qp_h, qp_A, qp_b)
        sol = cvxopt.solvers.qp(qp_P, qp_q, qp_G, qp_h,qp_A,qp_b)
    except Exception as inst:
        print "Warning : ", inst
        return False,0

    if sol['status'] == 'optimal':
        tauA = sol['x']
        tauA = array(tauA).reshape((6,))
        return True,tauA
    else:
        print sol['status']
    
    return False,0

# Inverse kinematics for a trajectory
def ComputeTorquesTraj(robot,traj,dt=0.01):
    robot1 = robot.robot1
    robot2 = robot.robot2
    T1 = robot.T1
    T2 = robot.T2
    Gdofs = robot.Gdofs
    Sdofs = robot.Sdofs
    actuated = robot.actuated
    constrainedlink = robot.constrainedlink
    taumin = robot.taumin
    taumax = robot.taumax
    ineqmatrices = ComputeInequalityConstraintMatrices(taumin,taumax)
    torquesvect = []
    for t in arange(0,traj.duration,dt):
        q = traj.Eval(t)
        qd = traj.Evald(t)
        qdd = traj.Evaldd(t)
        tauG,WT,ST = ComputeGeneralizedTorques(robot1,robot2,T1,T2,q,qd,qdd,Gdofs,Sdofs,actuated,constrainedlink)
        status, tauA = OptimizeTorques(ST,tauG,ineqmatrices)
        if status:
            torquesvect.append(tauA)
        else:
            print "Warning at t =",t
            if(len(torquesvect)>0):
                torquesvect.append(torquesvect[-1])
            else:
                torquesvect.append(zeros(6))
    return array(torquesvect)


def PlotTorques(robot,trajtotal,trajtotal2,dtplot=0.01):
    colorcycle = ['r', 'g', 'b', 'm', 'c', 'y', 'k']
    colorcycle = colorcycle[0:trajtotal.dimension]
    Tmax = max(trajtotal.duration, trajtotal2.duration)
    torques = ComputeTorquesTraj(robot,trajtotal)
    torques2 = ComputeTorquesTraj(robot,trajtotal2)
    tv = linspace(0,trajtotal.duration,len(torques))
    tv2 = linspace(0,trajtotal2.duration,len(torques2))
    figure()
    clf()
    hold("on")
    ax = gca()
    ax.set_color_cycle(colorcycle)
    plot(tv,torques,"--",linewidth=2)
    ax.set_color_cycle(colorcycle)
    plot(tv2,torques2,linewidth=2)
    for tau in robot.taumax:
        plot([0, Tmax], [tau, tau], '-.')
    for tau in robot.taumin:
        plot([0, Tmax], [tau, tau], '-.')
    axis([0, Tmax, 1.2*min(robot.taumin), 1.2*max(robot.taumax)])
    title('Joint torques', fontsize=20)
    xlabel('Time (s)', fontsize=18)
    ylabel('Joint torques (N.m)', fontsize=18)


###########################################################
######################## Bobrow ###########################
###########################################################

# Compute the constraint string for PolygonConstraints
def ComputeConstraints(robot, tunings, trajtotal):
    # Extract params from robot and tunings
    robot1 = robot.robot1
    robot2 = robot.robot2
    T1 = robot.T1
    T2 = robot.T2
    Gdofs = robot.Gdofs
    Sdofs = robot.Sdofs
    actuated = robot.actuated
    constrainedlink = robot.constrainedlink
    taumin = robot.taumin
    taumax = robot.taumax
    discrtimestep = tunings.discrtimestep

    ndiscrsteps = int((trajtotal.duration + 1e-10) / discrtimestep) + 1
    constraintstring = ""

    # Prepare matrices
    ndof = 6
    nvariable = ndof + 2
    lp_q = cvxopt.matrix(zeros(nvariable))
    lp_Gextended, lp_hextended = ComputeInequalityConstraintMatrices(taumin,taumax,True)

    a, b, c = zeros(5) , zeros(5), zeros(5)

    got_error = False
    for i in range(ndiscrsteps):
    #for i in [0]:
        t = i * discrtimestep
        #print t

        q = trajtotal.Eval(t)
        qd = trajtotal.Evald(t)
        qdd = trajtotal.Evaldd(t)
        q1 = q[0:3]
        qd1 = qd[0:3]
        qdd1 = qdd[0:3]
        q2 = q[3:5]
        qd2 = qd[3:5]
        qdd2 = qdd[3:5]

        # Inverse dynamics on the open chains
        with robot1:        
            robot1.SetTransform(T1)
            robot1.SetDOFValues(q1)
            robot1.SetDOFVelocities(qd1)
            tm, tc, tg = robot1.ComputeInverseDynamics(qdd1, None, returncomponents=True)
            to = robot1.ComputeInverseDynamics(qd1) - tc - tg
            a[0:3] = to
            b[0:3] = tm + tc
            c[0:3] = tg
        with robot2:
            robot2.SetTransform(T2)
            robot2.SetDOFValues(q2)
            robot2.SetDOFVelocities(qd2)
            tm, tc, tg = robot2.ComputeInverseDynamics(qdd2, None, returncomponents=True)
            to = robot2.ComputeInverseDynamics(qd2) - tc - tg
            a[3:5] = to
            b[3:5] = tm + tc
            c[3:5] = tg
        # Sensitivity matrices
        try:
            W, S = ComputeSensitivityMatrices(robot1,T1,T2,q,Gdofs,Sdofs,actuated,constrainedlink)
        except Exception as inst:
            print inst
            return None
        astar = dot(W.T,a)
        bstar = dot(W.T,b)
        cstar = dot(W.T,c)

        # A and b matrices
        A = zeros((3,8))
        A[:,:6] = S.T
        A[:,6] = -astar
        A[:,7] = -bstar
        lp_A = cvxopt.matrix(A)
        lp_b = cvxopt.matrix(cstar)
           
        # Compute the Polygon
        lp = lp_q, lp_Gextended, lp_hextended, lp_A, lp_b 
        res, P = ClosedChain.ComputePolygon(lp)
        if not res:
            got_error = True
            print "Compute Polygone: Infeasible at t =", t
        else:
            P.sort_vertices()
            # figure(6)
            # clf()
            # P.Plot()
            # raw_input()
            #return P, lp
            vertices_list = P.export_vertices()
            constraintstring += "\n"
            for ivertex in range(len(vertices_list)):
                constraintstring += str(vertices_list[ivertex].x) + " " + str(abs(vertices_list[ivertex].y)) #Note that it's sd^2 which is returned
                if ivertex < len(vertices_list)-1:
                    constraintstring += " "

    if got_error:
        return None
    return constraintstring


# Compute the constraint string for PolygonConstraints
def ComputeConstraintsTorqueOnly(robot, traj, taumin, taumax, discrtimestep):
    """Sample the dynamics constraints."""
    ndiscrsteps = int((traj.duration + 1e-10) / discrtimestep) + 1
    ndof = sum(taumax>0.01)
    constraintstring = ""    
    got_error = False
    for i in range(ndiscrsteps):
        t = i * discrtimestep
        q = traj.Eval(t)
        qd = traj.Evald(t)
        qd2 = zeros(robot.GetDOF())
        qd2[robot.GetActiveDOFIndices()] = qd
        qdd = traj.Evaldd(t)
        qdd2 = zeros(robot.GetDOF())
        qdd2[robot.GetActiveDOFIndices()] = qdd
        with robot:
            robot.SetActiveDOFValues(q)
            robot.SetActiveDOFVelocities(qd)
            tm, tc, tg = robot.ComputeInverseDynamics(qdd2, None,
                                                      returncomponents=True)
            to = robot.ComputeInverseDynamics(qd2) - tc - tg
            avect = []
            bvect = []
            cvect = []
            for i in range(len(taumax)):
                if abs(taumax[i]) > 1e-10:
                    avect.append(to[i])
                    bvect.append(tm[i] + tc[i])
                    cvect.append(tg[i])

        a = array(avect)
        b = array(bvect)
        c = array(cvect)

        # A and b matrices
        A = zeros((ndof,ndof+2))
        A[:,:ndof] = eye(ndof)
        A[:,ndof] = -a
        A[:,ndof+1] = -b
        lp_A = cvxopt.matrix(A)
        lp_b = cvxopt.matrix(c)

        # G and h matrices
        G = zeros((2*ndof+1,ndof+2))
        h = zeros(2*ndof+1)
        G[:ndof,:ndof] = eye(ndof)
        G[ndof:2*ndof,:ndof] = -eye(ndof)
        h[:ndof] = filter(lambda x: abs(x) > 1e-10 ,taumax)
        h[ndof:2*ndof] = filter(lambda x: abs(x) > 1e-10, -taumin)
        G[2*ndof,ndof+1] = -1
        lp_G = cvxopt.matrix(G)
        lp_h = cvxopt.matrix(h)

        lp_q = cvxopt.matrix(zeros(ndof+2))

        # Compute the Polygon
        lp = lp_q, lp_G, lp_h, lp_A, lp_b 
        res, P = ClosedChain.ComputePolygon(lp)
        if not res:
            got_error = True
            print "Compute Polygone: Infeasible at t =", t
        else:
            P.sort_vertices()
            # figure(6)
            # clf()
            # P.Plot()
            # raw_input()
            #return P, lp
            vertices_list = P.export_vertices()
            constraintstring += "\n"
            for ivertex in range(len(vertices_list)):
                constraintstring += str(vertices_list[ivertex].x) + " " + str(abs(vertices_list[ivertex].y)) #Note that it's sd^2 which is returned
                if ivertex < len(vertices_list)-1:
                    constraintstring += " "

    if got_error:
        return None
    return constraintstring



class Robot():
    def __init__(self):
        pass


class Tunings():
    def __init__(self):
        pass
