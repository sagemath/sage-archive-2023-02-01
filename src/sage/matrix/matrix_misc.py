def row_iterator(A):
    for i in xrange(A.nrows()):
        yield A.row(i)
