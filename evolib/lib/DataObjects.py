import PopGenStats

###### ######

class BinaryTable(list):
    
    def add_sample(self, item):
        s = len(item)
        
        if self.nsamples() is not None:
            assert s == self.segregating_sites(), "Cannot add a different number of segregating sites."
        
        for i in range(s):
            try:
                self[i] += item[i]
                
            except IndexError:
                self.append(item[i])
    
    def nsamples(self):
        
        if self == []:
            n = None
        else:
            n = len(self[0])
            
        return n
                
    def segregating_sites(self):
        
        if self == []:
            s = 0
        else:
            s = len(self)
            
        return s
    
    def thetaw(self):
        n = self.nsamples()
        s = self.segregating_sites()
        
        return PopGenStats.WattersonsTheta(n, s)
        
    def thetapi(self):
        pi = 0.0
        
        if self.nsamples() is not None:
            n = self.nsamples()
            pi = PopGenStats.TajimasTheta(n, self)
        
        return pi
    
    def tajimasd(self):
        n = self.nsamples()
        s = self.segregating_sites()
        tw = self.thetaw()
        pi = self.thetapi()
        
        return PopGenStats.TajimasD(n, s, tw, pi)

###### ######

class SequenceData():
    
    def __init__(self, *args):
        pass
