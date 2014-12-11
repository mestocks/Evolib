import evolib.stats.PopGenStats as PopGenStats

class IOstats(object):

    def __init__(self, IOdata):
        self.IOdata = self.IOdata

    def nsamples(self):
        
        try:
            n = len(self.IOdata[0])
        except IndexError:
            n = None

        return n

    def seg_sites(self):

        return len(self.IOdata)

    def thetaW(self):

        n = self.nsamples()
        s = self.seg_sites()
        
        return PopGenStats.WattersonsTheta(n, s)

    def thetaPi(self):
        
        pi = 0.0
        n = self.nsamples()
        if n is not None:
            pi = PopGenStats.TajimasThetaIO(n, self.IOdata)

        return pi

    def tajD(self):
        
        n = self.nsamples()
        s = self.seg_sites()
        tw = self.thetaW()
        pi = self.thetaPi()
        
        return PopGenStats.TajimasD(n, s, tw, pi)
    
    def wh97(self, pop_nsam):
        return PopGenStats.WakeleyHey(self.IOdata, pop_nsam)

class old_IOstats(object):

    def __init__(self, IOdata):
        self.IOdata = self.IOdata


    def nsamples(self):

        if self.IOdata == []:
            n = None
        else:
            n = len(self.IOdata[0])

        return n

    def seg_sites(self):

        s = 0
        if self.IOdata != []:
            for site in self.IOdata:
                nIO = len(set(site))
                if nIO > 1:
                    s += 1

        return s

    def thetaW(self):

        n = self.nsamples()
        s = self.seg_sites()

        return PopGenStats.WattersonsTheta(n, s)

    def thetaPi(self):
        
        pi = 0.0
        n = self.nsamples()
        if n is not None:
            n = self.nsamples()
            pi = PopGenStats.TajimasTheta(n, self.IOdata)

        return pi

    def tajD(self):
        
        n = self.nsamples()
        s = self.seg_sites()
        tw = self.thetaW()
        pi = self.thetaPi()
        
        return PopGenStats.TajimasD(n, s, tw, pi)
    
    def wh97(self, pop_nsam):
        return PopGenStats.WakeleyHey(self.IOdata, pop_nsam)
