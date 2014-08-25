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


cdef simple_transformation(M,rowtochange,basisrow,LP,U=None):
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

    Transforms M[rowtochange] into M[rowtochange]-a*x^d*M[basisrow], wherby d is the 
    difference of the degree of the leading positions and a is the division of the
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
    if U is not None:
        for i in range(U.ncols()):
            U[rowtochange,i] -= alpha*U[basisrow,i].shift(delta)


cpdef mulders_storjohann(M,transposition=False):
    r"""
    Function to transform M into weak popov form.

    INPUT:

     - `M` - Matrix over a polynomialring.
     
     - `transposition` - Boolean (default: False) indicating if a Matrix
     U should be computed so that U*M = M.weak_popov_form()
     
    OUTPUT:

    M transformed into weak popov form. If transposition is True, a touple
    (M,U) is returned with U*Original M = M in weak popov form.
    
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
    
    The value transposition can be used to get a second matrix to check 
    unimodular equivalence. ::
    
        sage: F.<a> = GF(2^4,'a')
        sage: PF.<x> = F[]
        sage: A = matrix(PF,[[1,a*x^17+1],[0,a*x^11+a^2*x^7+1]])
        sage: Ac = copy(A)
        sage: au = A.weak_popov_form(implementation="cython",transposition=True)
        sage: au[1]*Ac == au[0]
        True
        sage: au[1].is_invertible()
        True

    The cython implementation can be used to speed up the computation of
    a weak popov form. ::
    
        sage: B = matrix(PF,[[x^2+a,x^2+a,x^2+a], [x^3+a*x+1,x+a^2,x^5+a*x^4+a^2*x^3]])
        sage: B.weak_popov_form(implementation="cython")
        [                    x^2 + a                     x^2 + a                 
              x^2 + a]
        [x^5 + (a + 1)*x^3 + a*x + 1       x^5 + a*x^3 + x + a^2       a*x^4 + 
        (a^2 + a)*x^3]
    
    Matrices containing only zeros will return the way they are. ::
    
        sage: Z = matrix(PF,5,3)
        sage: Z.weak_popov_form(implementation="cython")
        [0 0 0]
        [0 0 0]
        [0 0 0]
        [0 0 0]
        [0 0 0]
    
    Generally matrices in weak popov form will just be returned. ::
    
        sage: F.<a> = GF(17,'a')
        sage: PF.<x> = F[]
        sage: C = matrix(PF,[[1,7,x],[x^2,x,4],[2,x,11]])
        sage: C.weak_popov_form(implementation="cython")
        [  1   7   x]
        [x^2   x   4]
        [  2   x  11]
    
    And the transposition will be the identity matrix. ::
    
        sage: C.weak_popov_form(implementation="cython",transposition=True)
        (
        [  1   7   x]  [1 0 0]
        [x^2   x   4]  [0 1 0]
        [  2   x  11], [0 0 1]
        )


    It is an error to call this function with a matrix not over a polynomial
    ring. ::

        sage: M = matrix([[1,0],[1,1]])
        sage: M.weak_popov_form(implementation="cython")
        Traceback (most recent call last):
        ...
        TypeError: the entries of M must lie in a univariate polynomial ring

    It is also an error to call this function using a matrix containing
    elements of the fraction field.

        sage: R.<t> = QQ['t']
        sage: M = matrix([[1/t,1/(t^2),t],[0,0,t]])
        sage: M.weak_popov_form(implementation="cython")
        Traceback (most recent call last):
        ...
        TypeError: the entries of M must lie in a univariate polynomial ring
            
    .. SEEALSO::

        :meth:`is_weak_popov <sage.matrix.matrix0.is_weak_popov>`
        
    REFERENCES::

    .. [MS] T. Mulders, A. Storjohann, "On lattice reduction for polynomial
          matrices," J. Symbolic Comput. 35 (2003), no. 4, 377--401
    """
    from sage.rings.polynomial.polynomial_ring import is_PolynomialRing
    if not is_PolynomialRing(M.base_ring()):
        raise TypeError("the entries of M must lie in a univariate polynomial ring")

    if transposition==True:
        from sage.matrix.constructor import identity_matrix
        U = identity_matrix(M.base_ring(),M.nrows())
    else:
        U = None
    lps = defaultdict(list)
    for c in range(M.nrows()):
        lp = leading_position(M[c])
        if M[c, lp] == 0:
            lps[-c-1].append(c)    # This enables easy success check while ignoring a zero row
        else:
            lps[lp].append(c)

    while len(lps) < M.nrows():
        for pos in lps:
            if len(lps[pos]) > 1:
                if (M[lps[pos][0]][pos].degree() >= M[lps[pos][1]][pos].degree()):
                    rowtochange = lps[pos][0]
                    basisrow = lps[pos][1]
                else:
                    rowtochange = lps[pos][1]
                    basisrow = lps[pos][0]
                simple_transformation(M,rowtochange,basisrow,pos,U)
                lps[pos].remove(rowtochange)
                
                if M[rowtochange, leading_position(M[rowtochange])] == 0:
                    lps[-rowtochange-1].append(rowtochange)    # ignore a line of pure zeros
                else:
                    lps[leading_position(M[rowtochange])].append(rowtochange)
                break
    if U is not None:
        return (M,U)
    return M
    
