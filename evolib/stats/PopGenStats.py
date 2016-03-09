import math
from scipy import misc

def Dxy(xa, xb, ya, yb):
    
    d = (xa * yb) + (xb * ya)
    n12 = (xa + xb) * (ya + yb)

    if d == 0:
        dxy = 0.0
    else:
        dxy = d / float(n12)

    return dxy

def WattersonsTheta(n, s):
    
    if s == 0:
        tw = 0.0
    else:
        a1 = sum(1.0 / i for i in range(1, n))
        tw = s / a1
    
    return tw

def TajimasThetaIO(n, io):
    
    sumPi = 0.0
    denom = misc.comb(n, 2)
    for site in io:
        nminor = sum(site)
        nmajor = n - nminor
        sumPi += nminor * nmajor
    pi = sumPi / denom
    
    return pi

def TajimasTheta(n, io):
    
    sumPi = 0.0
    denom = misc.comb(n, 2)
    for site in io:
        anc = site.count('0')
        der = n - anc
        sumPi += anc * der
    pi = sumPi / denom
    
    return pi
    
def TajimasD(n, s, tw, pi):
    
    if s == 0:
        D = None
    else:
        rawD = pi - tw
    
        a1 = 0.
        a2 = 0.
        for i in range(1, n):
            a1 += 1. / i
            a2 += 1. / (i * i)
        b1 = (n + 1.) / (3. * (n - 1.))
        b2 = 2. * (n * n + n + 3.) / (9. * n * (n - 1.))
        c1 = b1 - 1. / a1
        c2 = b2 - (n + 2.) / (a1 * n) + a2 / (a1 * a1)
        e1 = c1 / a1
        e2 = c2 / (a1 * a1 + a2)
        V = e1 * s + e2 * s * (s - 1.)

        D = rawD / math.sqrt(V)
    
    return D

def WakeleyHey(io, pop_nsam):
    """
    s1 - 001,000
    s2 - 000,001
    ss - 010,101
    sf - 000,111
    """
    nsites = len(io)
    assert len(pop_nsam) == 2
    n1, n2 = pop_nsam[0], pop_nsam[1]
    
    s1, s2, ss, sf = 0, 0, 0, 0
    for index in xrange(nsites):
        site = io[index]
        site1 = site[: n1]
        site2 = site[n1: n1 + n2]
        
        zero1 = sum(site1) < len(site1)
        zero2 = sum(site2) < len(site2)

        one1 = sum(site1) > 0
        one2 = sum(site2) > 0
        
        both1 = zero1 and one1
        both2 = zero2 and one2
        
        assert zero1 or zero2
        assert one1 or one2
        
        if both1 and not both2:
            s1 += 1
        elif not both1 and both2:
            s2 += 1
        elif both1 and both2:
            ss += 1
        elif not both1 and not both2:
            assert zero1 is not zero2
            assert one1 is not one2
            sf += 1
        else:
            raise TypeError, 'Unknown category'
        
    return s1, s2, ss, sf
        
def neutral_replacements(Ps, Pn, Ls, Ln):

    if Ls == 0 or Ln == 0 or Ps == 0:
        f = None
    else:
        if Pn == 0:
            f = 0.0
        else:
            N = float(Pn) / Ln
            S = float(Ps) / Ls
            f = N / S

    return f
