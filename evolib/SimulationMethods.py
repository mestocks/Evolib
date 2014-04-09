import sys

from evolib.SequenceFormats import msFormat

def ms_iter(fileObject):
    
    header = True
    iteration = ''
    
    for line in fileObject:
        
        if header is False:
            if line.rstrip() == '//':
                
                msClass = msFormat(iteration)
                yield msClass
                iteration = ''
                
            else:
                iteration += line
                
        elif header is True:
            if line.rstrip() == '//':
                header = False
                
    msClass = msFormat(iteration)
    yield msClass
