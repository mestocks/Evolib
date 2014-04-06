import math
from scipy import misc

def WattersonsTheta(n, s):
    
    if s == 0:
        tw = 0.0
    else:
        a1 = sum(1.0 / i for i in range(1, n))
        tw = s / a1
    
    return tw
    
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
