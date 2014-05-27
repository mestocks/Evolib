import sys

from evolib.SequenceFormats import msFormat

def ms_iter(fileObject):
    
    header = True
    iteration = ''
    
    for line in fileObject:
        
        if header is False:
            if line.startswith('//'):
                
                msClass = msFormat(iteration.rstrip())
                yield msClass
                iteration = ''
                
            else:
                iteration += line
                
        elif header is True:
            if line.startswith('//'):
                header = False
                
    msClass = msFormat(iteration.rstrip())
    yield msClass
