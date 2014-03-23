class PopStats():
    
    def seg_sites(self):
        return len(self.IOtable)
    
    def thetaw(self):
        
        if self.nsamples() is None:
            tw = 0.0
        else:
            s = self.seg_sites()
            n = self.nsamples()
            a1 = sum(1.0 / i for i in range(1, n))
            tw = s / a1

        return tw
