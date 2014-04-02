def BezierToPolynomial(T, p0, p1, p2, p3):
    a = -p0 + 3 * p1 - 3 * p2 + p3
    b = 3 * p0 - 6 * p1 + 3 * p2
    c = -3 * p0 + 3 * p1
    d = p0
    return a / (T * T * T), b / (T * T), c / T, d


def BezierToTrajectoryString(Tv, p0v, p1v, p2v, p3v):
    nchunks = len(Tv)
    dimension = len(p0v[0])
    trajectorystring = ""
    for i in range(nchunks):
        if i > 0:
            trajectorystring += "\n"
        trajectorystring += str(Tv[i]) + "\n" + str(dimension)
        for j in range(dimension):
            a, b, c, d = BezierToPolynomial(Tv[i], p0v[i][j], p1v[i][j],
                                            p2v[i][j], p3v[i][j])
            trajectorystring += "\n%f %f %f %f" % (d, c, b, a)
    return trajectorystring


def Interpolate3rdDegree(q0, q1, qd0, qd1, T):
    a = ((qd1 - qd0) * T - 2 * (q1 - q0 - qd0 * T)) / T ** 3
    b = (3 * (q1 - q0 - qd0 * T) - (qd1 - qd0) * T) / T ** 2
    c = qd0
    d = q0
    return a, b, c, d


def Interpolate5thDegree(q0, q1, qd0, qd1, qdd0, qdd1, T):
    a = (6*(q1 - q0) - 3*(qd1 + qd0)*T + 0.5*(qdd1 - qdd0)*(T**2))/(T**5)
    b = (-15*(q1 - q0) + (7*qd1 + 8*qd0)*T - 0.5*(2*qdd1 - 3*qdd0)*(T**2))/(T**4)
    c = (10*(q1 - q0) - 4*(qd1 + 1.5*qd0)*T + 0.5*(qdd1 - 3*qdd0)*(T**2))/(T**3)
    d = 0.5*qdd0
    e = qd0
    f = q0
    return a, b, c, d, e, f


def vect2str(v):
    return ' '.join(map(str, v))
