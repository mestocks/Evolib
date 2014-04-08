def loopByColumn(array):
    m = len(array)
    n = len(array[0])
    
    for j in range(n):
        column = ''
        for i in range(m):
            item = array[i][j]
            column += item
        
        yield column

def block_iter(itr, size = 3, start = 1):
    """
    Iterator that is aware of all items within non-overlapping 
    blocks of length 'size'. Returns the current item value in 
    the iterator along with all other items within it's local 
    bin. This is useful for scenarios where you may wish to apply 
    a method to the current item in the iteration that is dependent 
    on other items in the near vacinity. One example would be when 
    calculating statistics based on whether a mutation is amino acid 
    changing or not. The statistic calculated for each site depends 
    on the contents of the codon in which it is located. For example:
    
        dna = iter(['A', 'T', 'G', 'C', 
                    'A', 'T', 'G', 'C', 
                    'A', 'T', 'G', 'C'])
                        
        for base in backToTheFutureIter(dna):
            print base[0], ''.join(base[1])
        
    Would give:
            
            A ATG
            T ATG
            G ATG
            C CAT
            A CAT
            T CAT
            G GCA
            C GCA
            A GCA
            T TGC
            G TGC
            C TGC
    """
    assert size > 1
    assert start > 0
    assert start <= size
    
    i = start
    if i != 1:
        p = []
        for k in range(start, size + 1):
            nxt = itr.next()
            p.append(nxt)
    
        while i <= size:
            yield (p[(i - (1 + start - 1)) % size], p)
            i += 1
        
    while True:
        
        if i % size == 1:
            p = []
            for j in range(size):
                try:
                    nxt = itr.next()
                    p.append(nxt)
                except StopIteration:
                    pass
       
        try:
            yield (p[(i - 1) % size], p)
        except IndexError:
            raise StopIteration
        i += 1
