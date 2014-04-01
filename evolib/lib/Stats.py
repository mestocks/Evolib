import math
from scipy import misc

class PopStats():
    
    def seg_sites(self):
        s = len(self.IOtable)
        self.S = s
        return s
    
    def thetaW(self):
        
        if self.nsamples() is None:
            tw = 0.0
        else:
            s = self.seg_sites()
            n = self.nsamples()
            a1 = sum(1.0 / i for i in range(1, n))
            tw = s / a1
            
        self.tw = tw

        return tw
    
    def thetaPi(self):
        sumPi = 0.0
        nsamples = self.nsamples()
        if nsamples is None:
            pi = 0.0
        else:
            denom = misc.comb(nsamples, 2)
            for site in self.IOtable:
                anc = site.count('0')
                der = nsamples - anc
                sumPi += anc * der
            pi = sumPi / denom
            
        self.pi = pi
        
        return pi
    
    def tajD(self):
        
        try:
            S = self.S
            pi = self.pi
            tw = self.tw
        except AttributeError:
            S = self.seq_sites()
            pi = self.thetaPi()
            tw = self.thetaW()
        
        if S == 0:
            D = 'NaN'
        else:
            ns = self.nsamples()
            rawD = pi - tw
            
            a1 = 0.
            a2 = 0.
            for i in range(1, ns):
                a1 += 1. / i
                a2 += 1. / (i * i)
            b1 = (ns + 1.) / (3. * (ns - 1.))
            b2 = 2. * (ns * ns + ns + 3.) / (9. * ns * (ns - 1.))
            c1 = b1 - 1. / a1
            c2 = b2 - (ns + 2.) / (a1 * ns) + a2 / (a1 * a1)
            e1 = c1 / a1
            e2 = c2 / (a1 * a1 + a2)
            V = e1 * S + e2 * S * (S - 1.)
            D = rawD / math.sqrt(V)
    
        return D
