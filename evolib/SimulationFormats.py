from SequenceFormats2 import msFormat

class msClass():
    
    def __init__(self, data):
        self.data = data

    def __iter__(self):
        header = True
        iteration = ''
        for line in self.data:
            if header is False:
                if line.startswith('//'):
                    msIter = msFormat(iteration)
                    yield msIter
                    iteration = ''
                else:
                    iteration += line
            elif header is True:
                if line.startswith('//'):
                    header = False
        msIter = msFormat(iteration)
        yield msIter
