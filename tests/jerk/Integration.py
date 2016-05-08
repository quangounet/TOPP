import numpy as np
import sys
import logging
import TOPP

# Python logging utilities
_logger = logging.getLogger(__name__)
_loglevel = 10
_FORMAT = "[%(module)s::%(funcName)s] %(message)s"
logging.basicConfig(format=_FORMAT)
_logger.setLevel(_loglevel)

_debug = True


################################################################################
#                                 FORWARD
################################################################################

def IntegrateFWFollowingBeta(toppinstance, sstart, send, sdstart=0, dt=1e-3, 
                             checkjerk=False, topptraj=None, jerkmax=None):
    assert(sstart <= send) # check soundness
    if checkjerk:
        assert(topptraj is not None)
        assert(jerkmax is not None)
    svect = []
    sdvect = []
    sddvect = []

    scur = sstart
    sdcur = sdstart
    cont = True
    istep = 0
    while cont:
        if (sdcur < 0):
            message = "integrated profile hits sd = 0"
            break
        elif (scur > send):
            message = "exit safely"
            break

        svect.append(scur)
        sdvect.append(sdcur)
        # Now we are sure that it is safe to compute sddcur        
        sddcur = toppinstance.GetBeta(scur, sdcur)
        sddvect.append(sddcur)
        if checkjerk:
            if istep > 0:
                sddd = sdvect[-1] * (sddvect[-1] - sddvect[-2])/\
                (svect[-1] - svect[-2])
                gamma = GetGamma(topptraj, svect[-1], sdvect[-1], sddvect[-1], jerkmax)
                delta = GetDelta(topptraj, svect[-1], sdvect[-1], sddvect[-1], jerkmax)
                if _debug:
                    print 'sddd = {0}\ngamma = {1}\ndelta = {2}'.\
                    format(sddd, gamma, delta)
                if gamma > delta:
                    message = "entering gamma > delta region at ({0}, {1}, {2})".\
                    format(scur, sdcur, sddcur)
                    break
                if sddd > delta:
                    message = "sddd exceeds delta at ({0}, {1}, {2})".\
                    format(scur, sdcur, sddcur)
                    break
                elif sddd < gamma:
                    message = "sddd exceeds gamma at ({0}, {1}, {2})".\
                    format(scur, sdcur, sddcur)
                    break
                
        beta = sddcur
        alpha = toppinstance.GetAlpha(scur, sdcur)
        if (alpha > beta):
            message = "entering alpha > beta region"
            break

        # Integrate forward following beta
        istep += 1
        sdnext = sdcur + dt*beta
        snext = scur + dt*(sdcur + 0.5*dt*beta)

        scur = snext
        sdcur = sdnext
    
    _logger.info(message)
    return [svect, sdvect, sddvect]


def IntegrateFWFollowingAlpha(toppinstance, sstart, send, sdstart=0, dt=1e-3, 
                              checkjerk=False, topptraj=None, jerkmax=None):
    assert(sstart <= send) # check soundness
    if checkjerk:
        assert(topptraj is not None)
        assert(jerkmax is not None)
    svect = []
    sdvect = []
    sddvect = []

    scur = sstart
    sdcur = sdstart
    cont = True
    istep = 0
    while cont:
        if (sdcur < 0):
            message = "integrated profile hits sd = 0"
            break
        elif (scur > send):
            message = "exit safely"
            break

        svect.append(scur)
        sdvect.append(sdcur)
        # Now we are sure that it is safe to compute sddcur        
        sddcur = toppinstance.GetAlpha(scur, sdcur)
        sddvect.append(sddcur)
        if checkjerk:
            if istep > 0:
                sddd = sdvect[-1] * (sddvect[-1] - sddvect[-2])/\
                (svect[-1] - svect[-2])
                gamma = GetGamma(topptraj, svect[-1], sdvect[-1], sddvect[-1], jerkmax)
                delta = GetDelta(topptraj, svect[-1], sdvect[-1], sddvect[-1], jerkmax)
                if gamma > delta:
                    message = "entering gamma > delta region at ({0}, {1}, {2})".\
                    format(scur, sdcur, sddcur)
                    break
                if sddd > delta:
                    message = "sddd exceeds delta at ({0}, {1}, {2})".\
                    format(scur, sdcur, sddcur)
                    break
                elif sddd < gamma:
                    message = "sddd exceeds gamma at ({0}, {1}, {2})".\
                    format(scur, sdcur, sddcur)
                    break
                
        alpha = sddcur
        beta = toppinstance.GetBeta(scur, sdcur)
        if (alpha > beta):
            message = "entering alpha > beta region"
            break

        # Integrate forward following beta
        istep += 1
        sdnext = sdcur + dt*alpha
        snext = scur + dt*(sdcur + 0.5*dt*alpha)

        scur = snext
        sdcur = sdnext
    
    _logger.info(message)
    return [svect, sdvect, sddvect]
    

def IntegrateFWFollowingDelta(toppinstance, rawtopptraj, sstart, send, 
                              sdstart=0, sddstart=0, jerkmax=100, dt=1e-3):
    assert(sstart <= send) # check soundness

    svect = [sstart]
    sdvect = [sdstart]
    sddvect = [sddstart]
    ndof = rawtopptraj.dimension
    
    if type(jerkmax) == int or type(jerkmax) == float:
        J = np.ones(ndof)*jerkmax
    elif type(jerkmax) == np.ndarray:
        if len(jerkmax) == ndof:
            J = np.array(jerkmax)
        else:
            message = "jerkmax dimension mismatch"
            assert False, message
    else:
        message = "unknown jerkmax type"
        assert False, message

    scur = svect[-1]
    sdcur = sdvect[-1]
    sddcur = sddvect[-1]
    cont = True
    while cont:
        if (sdcur < 0):
            message = "integrated profile hits sd = 0"
            break
        elif (scur > send):
            message = "exit safely"
            break
        beta = toppinstance.GetBeta(scur, sdcur)
        alpha = toppinstance.GetAlpha(scur, sdcur)
        if (alpha > beta):
            message = "entering alpha > beta region"
            break
        
        # Check sddd bounds at the current step
        gamma = -np.infty
        delta = np.infty        
        
        qs = rawtopptraj.Evald(scur)
        qss = rawtopptraj.Evaldd(scur)
        qsss = rawtopptraj.Evaldn(scur, 3)

        bound1 = +J - 3*sdcur*sddcur*qss - (sdcur**3)*qsss
        bound2 = -J - 3*sdcur*sddcur*qss - (sdcur**3)*qsss

        for i in xrange(ndof):
            if qs[i] > 0:
                curmax = bound1[i] / qs[i]
                curmin = bound2[i] / qs[i]
            else:
                curmax = bound2[i] / qs[i]
                curmin = bound1[i] / qs[i]

            if curmin > gamma:
                gamma = curmin
            if curmax < delta:
                delta = curmax

        if gamma > delta:
            message = "entering gamma > delta region at ({0}, {1}, {2})".\
            format(scur, sdcur, sddcur)
            break

        # Integrate forward following delta
        sddnext = sddcur + dt*delta
        if sddnext > beta:
            message = "sddnext exceeds beta at ({0}, {1}, {2})".\
            format(scur, sdcur, sddcur)
            break
        elif sddnext < alpha:
            message = "sddnext exceeds alpha at ({0}, {1}, {2})".\
            format(scur, sdcur, sddcur)
            break

        sdnext = sdcur + dt*(sddcur + dt*delta/2.0)
        snext = scur + dt*(sdcur + dt*(sddcur + dt*delta/3.0)/2.0)

        svect.append(snext)
        sdvect.append(sdnext)
        sddvect.append(sddnext)

        scur = snext
        sdcur = sdnext
        sddcur = sddnext

        # print "current config:"
        # print "    s = {0}".format(scur)
        # print "    sd = {0}".format(sdcur)
        # print "    sdd = {0}".format(sddcur)
        # print "    delta = {0}".format(delta)

    _logger.info(message)
    return [svect, sdvect, sddvect]


def IntegrateFWFollowingGamma(toppinstance, rawtopptraj, sstart, send, 
                              sdstart=0, sddstart=0, jerkmax=100, dt=1e-3):
    assert(sstart <= send) # check soundness

    svect = [sstart]
    sdvect = [sdstart]
    sddvect = [sddstart]
    ndof = rawtopptraj.dimension
    
    if type(jerkmax) == int or type(jerkmax) == float:
        J = np.ones(ndof)*jerkmax
    elif type(jerkmax) == np.ndarray:
        if len(jerkmax) == ndof:
            J = np.array(jerkmax)
        else:
            message = "jerkmax dimension mismatch"
            assert False, message
    else:
        message = "unknown jerkmax type"
        assert False, message

    scur = svect[-1]
    sdcur = sdvect[-1]
    sddcur = sddvect[-1]
    cont = True
    while cont:
        if (sdcur < 0):
            message = "integrated profile hits sd = 0"
            break
        elif (scur > send):
            message = "exit safely"
            break
        beta = toppinstance.GetBeta(scur, sdcur)
        alpha = toppinstance.GetAlpha(scur, sdcur)
        if (alpha > beta):
            message = "entering alpha > beta region"
            break
        
        # Check sddd bounds at the current step
        gamma = -np.infty
        delta = np.infty        
        
        qs = rawtopptraj.Evald(scur)
        qss = rawtopptraj.Evaldd(scur)
        qsss = rawtopptraj.Evaldn(scur, 3)

        bound1 = +J - 3*sdcur*sddcur*qss - (sdcur**3)*qsss
        bound2 = -J - 3*sdcur*sddcur*qss - (sdcur**3)*qsss

        for i in xrange(ndof):
            if qs[i] > 0:
                curmax = bound1[i] / qs[i]
                curmin = bound2[i] / qs[i]
            else:
                curmax = bound2[i] / qs[i]
                curmin = bound1[i] / qs[i]

            if curmin > gamma:
                gamma = curmin
            if curmax < delta:
                delta = curmax

        if gamma > delta:
            message = "entering gamma > delta region at ({0}, {1}, {2})".\
            format(scur, sdcur, sddcur)
            break

        # Integrate forward following gamma
        sddnext = sddcur + dt*gamma
        if sddnext > beta:
            message = "sddnext exceeds beta at ({0}, {1}, {2})".\
            format(scur, sdcur, sddcur)
            break
        elif sddnext < alpha:
            message = "sddnext exceeds alpha at ({0}, {1}, {2})".\
            format(scur, sdcur, sddcur)
            break

        sdnext = sdcur + dt*(sddcur + dt*gamma/2.0)
        snext = scur + dt*(sdcur + dt*(sddcur + dt*gamma/3.0)/2.0)

        svect.append(snext)
        sdvect.append(sdnext)
        sddvect.append(sddnext)

        scur = snext
        sdcur = sdnext
        sddcur = sddnext

        # print "current config:"
        # print "    s = {0}".format(scur)
        # print "    sd = {0}".format(sdcur)
        # print "    sdd = {0}".format(sddcur)
        # print "    gamma = {0}".format(gamma)

    _logger.info(message)
    return [svect, sdvect, sddvect]


################################################################################
#                                 BACKWARD
################################################################################

def IntegrateBWFollowingAlpha(toppinstance, sstart, send, sdstart=0, dt=1e-3,
                              checkjerk=False, topptraj=None, jerkmax=None):
    assert(sstart >= send) # check soundness
    if checkjerk:
        assert(topptraj is not None)
        assert(jerkmax is not None)

    svect = []
    sdvect = []
    sddvect = []

    scur = sstart
    sdcur = sdstart
    cont = True
    istep = 0
    while cont:
        if (sdcur < 0):
            message = "integrated profile hits sd = 0"
            break
        elif (scur < send):
            message = "exit safely"
            break

        svect.insert(0, scur)
        sdvect.insert(0, sdcur)
        # Now we are sure that it is safe to compute sddcur        
        sddcur = toppinstance.GetAlpha(scur, sdcur)
        sddvect.insert(0, sddcur)
        if checkjerk:
            if istep > 0:
                sddd = sdvect[0] * (sddvect[1] - sddvect[0])/\
                (svect[1] - svect[0])
                gamma = GetGamma(topptraj, svect[0], sdvect[0], sddvect[0], jerkmax)
                delta = GetDelta(topptraj, svect[0], sdvect[0], sddvect[0], jerkmax)
                if _debug:
                    print 'sddd = {0}\ngamma = {1}\ndelta = {2}'.\
                    format(sddd, gamma, delta)
                if gamma > delta:
                    message = "entering gamma > delta region"
                    break
                if sddd > delta:
                    message = "sddd exceeds delta"
                    break
                elif sddd < gamma:
                    message = "sddd exceeds gamma"
                    break

        alpha = sddcur
        beta = toppinstance.GetBeta(scur, sdcur)
        if (alpha > beta):
            message = "entering alpha > beta region"
            break

        # Integrate forward following alpha
        istep += 1
        sdprev = sdcur - dt*alpha
        sprev = scur - dt*(sdcur - 0.5*dt*alpha)
        
        scur = sprev
        sdcur = sdprev
    
    _logger.info(message)
    return [svect, sdvect, sddvect]


def IntegrateBWFollowingDelta(toppinstance, rawtopptraj, sstart, send, 
                              sdstart=0, sddstart=0, jerkmax=100, dt=1e-3):
    assert(sstart >= send) # check soundness

    svect = [sstart]
    sdvect = [sdstart]
    sddvect = [sddstart]
    ndof = rawtopptraj.dimension
    
    if type(jerkmax) == int or type(jerkmax) == float:
        J = np.ones(ndof)*jerkmax
    elif type(jerkmax) == np.ndarray:
        if len(jerkmax) == ndof:
            J = np.array(jerkmax)
        else:
            message = "jerkmax dimension mismatch"
            assert False, message
    else:
        message = "unknown jerkmax type"
        assert False, message

    scur = svect[-1]
    sdcur = sdvect[-1]
    sddcur = sddvect[-1]
    cont = True
    while cont:
        if (sdcur < 0):
            message = "integrated profile hits sd = 0"
            break
        elif (scur < send):
            message = "exit safely"
            break
        beta = toppinstance.GetBeta(scur, sdcur)
        alpha = toppinstance.GetAlpha(scur, sdcur)
        if (alpha > beta):
            message = "entering alpha > beta region"
            break
        
        # Check sddd bounds at the current step
        gamma = -np.infty
        delta = np.infty        
        
        qs = rawtopptraj.Evald(scur)
        qss = rawtopptraj.Evaldd(scur)
        qsss = rawtopptraj.Evaldn(scur, 3)

        bound1 = +J - 3*sdcur*sddcur*qss - (sdcur**3)*qsss
        bound2 = -J - 3*sdcur*sddcur*qss - (sdcur**3)*qsss

        for i in xrange(ndof):
            if qs[i] > 0:
                curmax = bound1[i] / qs[i]
                curmin = bound2[i] / qs[i]
            else:
                curmax = bound2[i] / qs[i]
                curmin = bound1[i] / qs[i]

            if curmin > gamma:
                gamma = curmin
            if curmax < delta:
                delta = curmax

        if gamma > delta:
            message = "entering gamma > delta region"
            break

        # Integrate forward following delta
        sddprev = sddcur - dt*delta
        if sddprev > beta:
            message = "sddprev exceeds beta"
            break
        elif sddprev < alpha:
            message = "sddprev exceeds alpha"
            break

        sdprev = sdcur - dt*(sddcur - dt*delta/2.0)
        sprev = scur - dt*(sdcur - dt*(sddcur - dt*delta/3.0)/2.0)

        svect.insert(0, sprev)
        sdvect.insert(0, sdprev)
        sddvect.insert(0, sddprev)

        scur = sprev
        sdcur = sdprev
        sddcur = sddprev

        # print "current config:"
        # print "    s = {0}".format(scur)
        # print "    sd = {0}".format(sdcur)
        # print "    sdd = {0}".format(sddcur)
        # print "    delta = {0}".format(delta)

    _logger.info(message)
    return [svect, sdvect, sddvect]


def IntegrateBWFollowingGamma(toppinstance, rawtopptraj, sstart, send, 
                              sdstart=0, sddstart=0, jerkmax=100, dt=1e-3):
    assert(sstart >= send) # check soundness

    svect = [sstart]
    sdvect = [sdstart]
    sddvect = [sddstart]
    ndof = rawtopptraj.dimension
    
    if type(jerkmax) == int or type(jerkmax) == float:
        J = np.ones(ndof)*jerkmax
    elif type(jerkmax) == np.ndarray:
        if len(jerkmax) == ndof:
            J = np.array(jerkmax)
        else:
            message = "jerkmax dimension mismatch"
            assert False, message
    else:
        message = "unknown jerkmax type"
        assert False, message

    scur = svect[-1]
    sdcur = sdvect[-1]
    sddcur = sddvect[-1]
    cont = True
    while cont:
        if (sdcur < 0):
            message = "integrated profile hits sd = 0"
            break
        elif (scur < send):
            message = "exit safely"
            break
        beta = toppinstance.GetBeta(scur, sdcur)
        alpha = toppinstance.GetAlpha(scur, sdcur)
        if (alpha > beta):
            message = "entering alpha > beta region"
            break
        
        # Check sddd bounds at the current step
        gamma = -np.infty
        delta = np.infty        
        
        qs = rawtopptraj.Evald(scur)
        qss = rawtopptraj.Evaldd(scur)
        qsss = rawtopptraj.Evaldn(scur, 3)

        bound1 = +J - 3*sdcur*sddcur*qss - (sdcur**3)*qsss
        bound2 = -J - 3*sdcur*sddcur*qss - (sdcur**3)*qsss

        for i in xrange(ndof):
            if qs[i] > 0:
                curmax = bound1[i] / qs[i]
                curmin = bound2[i] / qs[i]
            else:
                curmax = bound2[i] / qs[i]
                curmin = bound1[i] / qs[i]

            if curmin > gamma:
                gamma = curmin
            if curmax < delta:
                delta = curmax

        if gamma > delta:
            message = "entering gamma > delta region"
            break

        # Integrate forward following gamma
        sddprev = sddcur - dt*gamma
        if sddprev > beta:
            message = "sddprev exceeds beta"
            break
        elif sddprev < alpha:
            message = "sddprev exceeds alpha"
            break

        sdprev = sdcur - dt*(sddcur - dt*gamma/2.0)
        sprev = scur - dt*(sdcur - dt*(sddcur - dt*gamma/3.0)/2.0)

        svect.insert(0, sprev)
        sdvect.insert(0, sdprev)
        sddvect.insert(0, sddprev)

        scur = sprev
        sdcur = sdprev
        sddcur = sddprev

        # print "current config:"
        # print "    s = {0}".format(scur)
        # print "    sd = {0}".format(sdcur)
        # print "    sdd = {0}".format(sddcur)
        # print "    gamma = {0}".format(gamma)

    _logger.info(message)
    return [svect, sdvect, sddvect]
    

################################################################################
#                                Utilities
################################################################################

def GetGamma(topptraj, s, sd, sdd, jerkmax):
    ndof = topptraj.dimension
    if type(jerkmax) == int or type(jerkmax) == float:
        J = np.ones(ndof)*jerkmax
    elif type(jerkmax) == np.ndarray:
        if len(jerkmax) == ndof:
            J = np.array(jerkmax)
        else:
            message = "jerkmax dimension mismatch"
            assert False, message
    else:
        message = "unknown jerkmax type"
        assert False, message


    gamma = -np.infty
        
    qs = topptraj.Evald(s)
    qss = topptraj.Evaldd(s)
    qsss = topptraj.Evaldn(s, 3)

    bound1 = +J - 3*sd*sdd*qss - (sd**3)*qsss
    bound2 = -J - 3*sd*sdd*qss - (sd**3)*qsss

    for i in xrange(ndof):
        if qs[i] > 0:
            curmax = bound1[i] / qs[i]
            curmin = bound2[i] / qs[i]
        else:
            curmax = bound2[i] / qs[i]
            curmin = bound1[i] / qs[i]

        if curmin > gamma:
            gamma = curmin

    return gamma


def GetDelta(topptraj, s, sd, sdd, jerkmax):
    ndof = topptraj.dimension
    if type(jerkmax) == int or type(jerkmax) == float:
        J = np.ones(ndof)*jerkmax
    elif type(jerkmax) == np.ndarray:
        if len(jerkmax) == ndof:
            J = np.array(jerkmax)
        else:
            message = "jerkmax dimension mismatch"
            assert False, message
    else:
        message = "unknown jerkmax type"
        assert False, message


    delta = np.infty
        
    qs = topptraj.Evald(s)
    qss = topptraj.Evaldd(s)
    qsss = topptraj.Evaldn(s, 3)

    bound1 = +J - 3*sd*sdd*qss - (sd**3)*qsss
    bound2 = -J - 3*sd*sdd*qss - (sd**3)*qsss

    for i in xrange(ndof):
        if qs[i] > 0:
            curmax = bound1[i] / qs[i]
            curmin = bound2[i] / qs[i]
        else:
            curmax = bound2[i] / qs[i]
            curmin = bound1[i] / qs[i]

        if curmax < delta:
            delta = curmax

    return delta


def PlotVectorFields(toppinstance, topptraj, ax, jmax, prec=10, length=0.3):
    s_lim = ax.get_xlim()
    sd_lim = ax.get_ylim()
    sdd_lim = ax.get_zlim()
    s_coord = np.linspace(s_lim[0], s_lim[1], prec)
    sd_coord = np.linspace(sd_lim[0], sd_lim[1], prec)
    sdd_coord = np.linspace(sdd_lim[0], sdd_lim[1], prec)

    ds0 = s_coord[1] - s_coord[0]
    dsd0 = sd_coord[1] - sd_coord[0]
    dsdd0 = sdd_coord[1] - sdd_coord[0]
    yscl = dsd0 / ds0
    zscl = dsdd0 / ds0

    ds = ds0 / 2.0
    
    # (x, y, z)-position of the vectors
    svect = []
    sdvect = []
    sddvect = []
    
    # (x, y, z)-components of the vectors
    u_min = []
    v_min = []
    w_min = []
    u_max = []
    v_max = []
    w_max = []

    from itertools import product
    for (s, sd, sdd) in product(s_coord, sd_coord, sdd_coord):
        if abs(sd) < 1e-6:
            sd = 1e-6
        alpha = toppinstance.GetAlpha(s, sd)
        beta = toppinstance.GetBeta(s, sd)
        gamma = GetGamma(topptraj, s, sd, sdd, jmax)
        delta = GetDelta(topptraj, s, sd, sdd, jmax)

        na = 1.0/np.sqrt(1 + (alpha / yscl)**2)
        nb = 1.0/np.sqrt(1 + (beta / yscl)**2)
        ng = 1.0/np.sqrt(1 + (gamma / zscl)**2)
        nd = 1.0/np.sqrt(1 + (delta / zscl)**2)

        min_vect = np.array([1, alpha/sd, gamma/sd])
        min_vect = min_vect / np.linalg.norm(min_vect)

        max_vect = np.array([1, beta/sd, delta/sd])
        max_vect = max_vect / np.linalg.norm(max_vect)

        svect.append(s)
        sdvect.append(sd)
        sddvect.append(sdd)

        u_min.append(min_vect[0])
        v_min.append(min_vect[1])
        w_min.append(min_vect[2])

        u_max.append(max_vect[0])
        v_max.append(max_vect[1])
        w_max.append(max_vect[2])    

    ax.quiver(svect, sdvect, sddvect, u_max, v_max, w_max, 
              color='g', pivot='tail', arrow_length_ratio=0, length=length)
    ax.quiver(svect, sdvect, sddvect, u_min, v_min, w_min, 
              color='r', pivot='tail', arrow_length_ratio=0, length=length)
