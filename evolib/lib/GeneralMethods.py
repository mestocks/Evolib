def loopByColumn(array):
    m = len(array)
    n = len(array[0])
    
    for j in range(n):
        column = ''
        for i in range(m):
            item = array[i][j]
            column += item
        
        yield column
