class FastqRead():
    def __init__(self, seq, seqid, quality):
        self.seq = seq
        self.seqid = seqid
        self.quality = quality
        
    def __str__(self):
        return '\n'.join([self.seqid, self.seq, '+', self.quality])
