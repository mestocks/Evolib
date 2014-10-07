from evolib.formats.IteratorObjects import FastaAlignment

class FastaFormat(FastaAlignment):
    """
    Example usage:
    
       >>> from evolib.SequenceFormats import FastaFormat
       
       Example 1:
       >>> fileObject = open('example.fsa', 'r')
       >>> F = FastaFormat(fileObject)
       
       Example 2:
       >>> import sys
       >>> stdinObject = sys.stdin
       >>> F = FastaFormat(stdinObject)
       
    """
    pass
