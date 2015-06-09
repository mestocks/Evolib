
def parse_fasta_alignment(fileObject):
    """
    Returns a matrix where M[i][j] refers to the jth site of the ith individual.
    """        
    step1 = ''.join(map(lambda line: line, fileObject))
    step2 = step1.split('\n>')
    seq_table = [part.partition('\n')[2].replace('\n','') for part in step2]
    seq_names = [part.partition('\n')[0].replace('>', '') for part in step2]

    return seq_table, seq_names

def idsets(fileObject):
    """
    File format:

        <key1>:<id1>,<id2>,...
        <key2>:<id1>,<id2>,...

    Returns:

        {<key1>: [<id1>, <id2>, ...], 
         <key2>: [<id1>, <id2>, ...]}
    
    """
    rdict = {}
    for line in fileObject:
        ksplit = line.rstrip().split(":")
        rdict.update({ksplit[0]: ksplit[1].split(',')})

    return rdict
