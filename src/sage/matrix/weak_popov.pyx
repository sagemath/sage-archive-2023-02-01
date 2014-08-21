"""
Mulders-Storjohann algorithm to compute the weak popov form of polynomial matrices.

AUTHORS:

- David Mödinger (2014-08-21: initial version

"""
from collections import defaultdict

cdef leading_position(v):
    r"""
    Used to compute the leading position of a vector v.

    INPUT:

     - `v` - vector

    OUTPUT:

    Outputs the leading position of a vector v, which is the position with highest degree. For multiple positions with equal degrees the highest position i, or rightmost in the vector, is chosen.
    
    .. note::
    
        This method is used in the mulders-storjohann algorithm.

    """
    p = -1 # pos of max
    m = -1 # max
    for c in range(v.degree()):
        if(v[c].degree()>=m):
            m = v[c].degree()
            p = c
    return p


cdef simple_transformation(M,rowtochange,basisrow,LP):
    r"""
    Function to compute a simple transformation on a matrix.

    INPUT:

     - `M` - Matrix to operate on.
     
     - `rowtochange` - Integer to indicate the row M[rowtochange] that 
     is to be changed.
     
     - `basisrow` - Integer to indicate the row M[basisrow] that is used 
     to change M[rowtochange].
     
     - `LP` - Position of the leading positions in the two vectors.

    OUTPUT:

        Transforms M[rowtochange] into M[rowtochange]-a*x^d*M[basisrow], wherby d is the difference of
        the degree of the leading positions and a is the division of the
        leading coefficients.
    
    .. WARNING::
    
        The degree of the leadingposition of rowtochange has to be greater or equal
        to the degree of the leadingposition of basisrow.
    
    .. note::
    
        This method is used in the mulders-storjohann algorithm.

    """
    cdef delta = M[rowtochange][LP].degree()-M[basisrow][LP].degree()
    cdef alpha = (M[rowtochange][LP].coefficients()[-1]) / (M[basisrow][LP].coefficients()[-1])
    for i in range(M.ncols()):
        M[rowtochange,i] -= alpha*M[basisrow,i].shift(delta)
    return


cpdef mulders_storjohann(M):
    r"""
    Function to transform M into weak popov form.

    INPUT:

     - `M` - Matrix.
     
    OUTPUT:

        M transformed into weak popov form.
    
    ALGORITHM::
    
        This function uses the mulders-storjohann algorithm of [MS].
        It works as follow:
        #. As long as M is not in weak popov form do:
            #. Find two rows with conflicting leading positions.
            #. Do a simple transformation:
                #. Let x and y be indicators of rows with identical leading position
                #. Let LP be the Leading Position and LC the Leading Coefficient
                #. let a = LP(M[x]).degree() - LP(M[y]).degree()
                #. let d = LC(LP(M[x])) / LC(LP(M[y]))
                #. substitute M[x] = M[x] - a * x^d * M[y]

    EXAMPLES:
        
    
        
        
    REFERENCES::

    .. [MS] T. Mulders, A. Storjohann, "On lattice reduction for polynomial
          matrices," J. Symbolic Comput. 35 (2003), no. 4, 377--401
          
    """
    lps = defaultdict(list)
    for c in range(M.nrows()):
        lp = leading_position(M[c])
        if not M[c,lp]==-1:
            lps[lp].append(c)
    
    while len(lps)<M.nrows():
        for pos in lps:
            if len(lps[pos])>1:
                if (M[lps[pos][0]][pos].degree() >= M[lps[pos][1]][pos].degree()):
                    arownr = lps[pos][0]
                    brownr = lps[pos][1]
                else:
                    arownr = lps[pos][1]
                    brownr = lps[pos][0]
                simple_transformation(M,arownr,brownr,pos)
                lps[pos].remove(arownr)
                lps[leading_position(M[arownr])].append(arownr) 
                break
    return M
    
