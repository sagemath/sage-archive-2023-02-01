r"""
Computes upper cluster algebra elements for a given matrix. 

This file implements the combinatorial formula given in [LLM] for nontrivial elements
 of the upper cluster algebra associated to an ``m \times n`` matrix, ``B``. 
If ``B`` is acyclic then elements of this form form a basis for the (upper) cluster algebra.
If ``B`` is not acyclic then this is not true, but the elements can still be useful. 


AUTHORS::

- Matt Mills (2015-06-19): initial version

EXAMPLES::

    sage: B=matrix([[0,2],[-2,0]])
    sage: get_uca_element(B,[5,3])
    x0^-5 * x1^-3 * (x1^2 + 1)^2 * (x0^2 + x1^2 + 1)^3

    sage: get_uca_element(B,[-2,3])
    x1**-3 * x0**2 * (x0**2 + 1)**3

    sage: get_uca_element(B,[-2,-3])
    x0^2 * x1^3

    sage: LLM_gen_set(B)
    [1, x1^-1 * (x0^2 + 1), x0^-1 * (x1^2 + 1), x1^-1 * x0^-1 * (x0^2 + x1^2 + 1)]


    sage: B=matrix([[0,3,2],[-3,0,2],[-2,-2,0]])
    sage: get_uca_element(B,[1,2,3])
    x2^-3 * x1^-2 * x0^-1 * (x0^2*x1^2 + 1) * (x0^5*x1^2 + x0^3 + x2^2) * (x0^5*x1^2 + x1^3*x2^4 + x0^3 + x2^2)
    sage: LLM_gen_set(B)
    [1,\
     x2^-1 * (x0^2*x1^2 + 1),
     x1^-1 * (x0^3 + x2^2),
     x2^-1 * x1^-1 * (x0^5*x1^2 + x0^3 + x2^2),
     x0^-1 * (x1^3*x2^2 + 1),
     x2^-1 * x0^-1 * (x1^3*x2^2 + x0^2*x1^2 + 1),
     x1^-1 * x0^-1 * (x1^3*x2^4 + x0^3 + x2^2),
     x2^-1 * x1^-1 * x0^-1 * (x0^5*x1^2 + x1^3*x2^4 + x0^3 + x2^2)


    sage: B=matrix([[0,1,0,0],[-1,0,1,1],[0,-1,0,0],[0,-1,0,0]])
    sage: get_uca_element(B,[1,2,3,4])
    x3^-4 * x2^-3 * x1^-2 * x0^-1 * (x1 + 1)^4 * (x0*x1 + x2*x3 + x0) * (x0*x1^2 + 2*x0*x1 + x2*x3 + x0)
    sage: LLM_gen_set(B)
    [1,
     x3^-1 * (x1 + 1),
     x2^-1 * (x1 + 1),
     x3^-1 * x2^-1 * (x1 + 1)^2,
     x1^-1 * (x2*x3 + x0),
     x3^-1 * x1^-1 * (x0*x1 + x2*x3 + x0),
     x2^-1 * x1^-1 * (x0*x1 + x2*x3 + x0),
     x3^-1 * x2^-1 * x1^-1 * (x0*x1^2 + 2*x0*x1 + x2*x3 + x0),
     x0^-1 * (x1 + 1),
     x3^-1 * x0^-1 * (x1 + 1)^2,
     x2^-1 * x0^-1 * (x1 + 1)^2,
     x3^-1 * x2^-1 * x0^-1 * (x1 + 1)^3,
     x1^-1 * x0^-1 * (x1*x2*x3 + x2*x3 + x0),
     x3^-1 * x1^-1 * x0^-1 * (x1 + 1) * (x2*x3 + x0),
     x2^-1 * x1^-1 * x0^-1 * (x1 + 1) * (x2*x3 + x0),
     x3^-1 * x2^-1 * x1^-1 * x0^-1 * (x1 + 1) * (x0*x1 + x2*x3 + x0)]
"""

#*****************************************************************************
#       Copyright (C) 2013 Matt Mills <fi8380@wayne.edu>
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 2 of the License, or
# (at your option) any later version.
#                  http://www.gnu.org/licenses/
#*****************************************************************************

def _vector_decomposition(a,length):
    r"""
    Decomposes an integer vector 

    INPUT:

    A vector ``a \in  \mathbb{Z}^n.''

    OUTPUT:

    A decomposition of ``a`` into vectors ``b_i \in \{0,1\}^n`` such that ``a= \sum c_i b_i`` for ``c_i \in \mathbb{Z}.``
    Returns an array of tuples ``\right[b_i,c_i\left].`` 

    EXAMPLES::

        sage: _vector_decomposition([2,-1,3,-2],4)
        [[(1, 0, 1, 0), 2], [(0, 0, 1, 0), 1], [(0, 0, 0, -1), 1], [(0, -1, 0, -1), 1]]
    
        sage: _vector_decomposition([3,2,3,4],4)
        [[(1, 1, 1, 1), 2], [(1, 0, 1, 1), 1], [(0, 0, 0, 1), 1]]
    """
    #Finds the difference between the largest and smallest entry in the vector to determine the how many vectors are in the decomposition
    max = 0
    min = 0
    for i in range(len(a)):
        if a[i] > max:
            max=a[i]
        if a[i] < min:
            min =a[i]
    diff = max-min

    #Creates a copy of a that will be edited when decomposing the vector.  
    ap=copy(a)
    if max == 0 and min == 0:
        ap=[]
        for i in range(len(a)):
            ap.append(0)
        return [[ap,1]]
    #Resets the counter i and puts the integer partition of the ith component of a into an array. 
    i=0
    cols=[]
    for i in range(len(a)):
        c=[]
        for j in range(diff):
            if ap[i]>0:
                c.append(1)
                ap[i]-=1
            elif ap[i]<0:
                c.append(-1)
                ap[i]+=1
            elif ap[i]==0:
                c.append(0)
        cols.append(c)
    #Converts the integer partitions into decomposition vectors.
    i=0
    for i in range(len(cols)):
        if cols[i][0]<0:
            cols[i].reverse()
    mat=matrix(cols)
    #Adds a zero to the end of every vector for each frozen vertex. 
    froz_mat=matrix(length-mat.nrows(),mat.ncols())
    mat=mat.stack(froz_mat)
    mat=mat.transpose()
    vects=mat.rows()
    #Collects identical decomposition vectors and counts their multiplicities. 
    multiList=[]
    while(len(vects)>0):
        vect=vects[0]
        count=vects.count(vect)
        multiList.append([vect,count])
        i=0
        for i in range(count):
            vects.remove(vect)
    return multiList

def _compute_compatible_vectors(B,vd):
    r"""
    Returns a list of compatible vectors of each vector in the vector decomposition ``vd``.
    Compatibility is defined as in [LLM] with respect to the matrix ``B``.

    INPUT::

    - ``B`` -- a skew-symmetric matrix. Must have the same number of columns as the length of the vectors in ``vd``.
    - ``vd`` -- a collection of tuples ``(v,z)`` with ``v \in \{0,1\}^n`` and ``z \in \mathbb{Z}``.
                ``n`` must be the number of columns in ``B``. Taken from the output of vector_decomposition.

    OUTPUT::

    Returns an a 2-dimensional array containing all the vectors compatible with each vector in ``vd.`` 

    NOTE:

    If the vector in ``vd`` is negative it will not have any compatible vectors, so it does not contribute to the list.

    EXAMPLES::

        sage: _compute_compatible_vectors(B,v)
        [[[0, 0, 0, 0],
        [0, 0, 0, 1],
        [0, 0, 1, 0],
        [0, 0, 1, 1],
        [0, 1, 1, 1],
        [1, 1, 1, 1]],
        [[0, 0, 0, 0],
        [0, 0, 0, 1],
        [0, 0, 1, 0],
        [0, 0, 1, 1],
        [1, 0, 0, 0],
        [1, 0, 0, 1],
        [1, 0, 1, 0],
        [1, 0, 1, 1]],
        [[0, 0, 0, 0], [0, 0, 0, 1]]]
    
        sage: B=matrix([[0,1,1,0],[-1,0,1,1],[-1,-1,0,0],[0,-1,0,0]])
        sage: v=_vector_decomposition([2,-1,3,-2],4)
        sage: _compute_compatible_vectors(B,v)
        [[[0, 0, 0, 0], [0, 0, 1, 0], [1, 0, 1, 0]], [[0, 0, 0, 0], [0, 0, 1, 0]]]  
    """

    #E is the set of 'edges' in the quiver. It records the tuple of indices ``(i,j)`` if ``b_{ij}>0``.
    E=[]
    #Checks the upper triangular part of the exchange graph.
    num_cols=B.ncols()
    num_rows=B.nrows()
    for j in range(num_cols):
        for i in range(j,num_rows):
            if B[i][j] > 0:
                E.append([i,j])
            elif B[i][j] < 0:
                E.append([j,i])
    #Checks for edges to frozen vertices. 
    num_frozens=num_rows-num_cols
    for k in range(num_frozens):
        j=0
        for j in range(i,num_cols):
            if B[k+num_cols][j] > 0:
                E.append([i,j])
            elif B[i][j] < 0:
                E.append([j,i])

    #For each vector a in vd. check if a vector s in {0,1}^n is compatible.
    compatibleList=[]
    psetvect = _power_set(num_rows)
    for a in vd:
        negative=false
        for m in xrange(len(a)):
    #If the vector a in vd is non-positive it is not compatible with any vector. 0 vector will pass this check but will be handled later.
            if a[m]<0:
                negative=true
                break
        if negative == true:
            continue
        clist=[]
        for s in psetvect:
            pass1=true
    #The first possible failure for compatibility is if any entry in s is larger than the corresponding entry of a.
            for k in xrange(num_rows):
                if s[k]>a[0][k]:
                    pass1=false
                    break
    #The second possible failure is if (s_i,a_j-s_j) = (1,1).
            if pass1 == true:
                for  e in E:
                    if s[e[0]]==1 and (a[0][e[1]]-s[e[1]]) == 1:
                        pass1=false
                        break
            if pass1 == true:
                clist.append(s)
        compatibleList.append(clist)
    return compatibleList

def _produce_uca_element(B,vd,cList):
    r"""
    Takes the compatible vectors and uses them to produce a Laurent polynomial in the upper cluster algebra. 

    EXAMPLES::

        sage: B=matrix([[0,1,0,0],[-1,0,1,1],[0,-1,0,0],[0,-1,0,0]])
        sage: v=_vector_decomposition([1,2,1,2],4)
        sage: c=_compute_compatible_vectors(B,v)
        sage: _produce_uca_element(B,v,c)
        x3^-2 * x1^-2 * x2^-1 * x0^-1 * (x1 + 1) * (x0*x1 + x2*x3 + x0)^2

        sage: B=matrix([[0,1,1,0],[-1,0,1,1],[-1,-1,0,0],[0,-1,0,0]])   
        sage: v=_vector_decomposition([2,-1,3,-2],4)
        sage: c=_compute_compatible_vectors(B,v)
        sage: _produce_uca_element(B,v,c)
        x2^-3 * x0^-2 * x1 * x3^2 * (x0*x1 + 1) * (x0*x1 + x1*x2 + 1)^2
    """
    #Creates a the fraction field of a polynomial ring in which to build the Laurent polynomials.
    num_cols=B.ncols()
    num_rows=B.nrows()
    R=PolynomialRing(QQ,num_rows,'x')
    R.fraction_field()
    #Computes the Laurent Polynomial for each vector in the decomposition.
    finalP=[]
    #Laurent polynomial for each vector in {0,1}^n
    for i in range(len(vd)):  
        final=1
        numerator=0
        #If the vector in vd is negative then it did not contribute any compatible vectors. It will only contribute a Laurent monomial.
        if len(cList)>i:
        #Each compatible sequence gives a term in the numerator of the Laurent polynomial.
            for s in cList[i]:  
                term=1
                #Calulates the monomial in the term. 
                for j in range(num_rows): 
                    x=R.gen(j)
                    expn=0
                    #The exponent is determined by the vectors a,s, and the matrix B.
                    for k in range(num_cols):
                        expn+=((vd[i][0][k]-s[k])*_zero_max(B[j][k])+s[k]*_zero_max(-B[j][k]))
                    term=term*(x**expn)
                numerator+=term
        #Gives a numerator for the negative vector, or else the product would be zero.      
        else:
            numerator=1
        #Uses the vectors in vd to calculates the denominator of the Laurent.     
        denominator=1
        for l in range(num_cols):
            denominator=denominator*(R.gen(l))**vd[i][0][l]
        #Each copy of a vector in vd contributes a factor of the Laurent polynomial calculated from it. 
        final=(numerator/denominator)**vd[i][1]
        finalP.append(final)
    laurentP=1
    #The UCA element for the vector a is the product of the elements produced from the vectors in its decomposition. 
    for p in finalP:
        laurentP=laurentP*p
    return factor(laurentP)

def get_uca_element(B,a):
    r"""
    Computes an element in the upper cluster algebra of ``B`` corresponding to the vector ``a \in \mathbb{Z}^n``.

    INPUT::

    - ``B`` -- a skew-symmetric matrix. Must have the same number of columns as the length of the vectors in ``vd``.
    - ``a`` -- a vector in ``\mathbb{Z}^n`` where ``n`` is the number of columns in ``B``.

    OUTPUT::

    Returns an element in the upper cluster algebra. Depending on the input it may or may not be irreducible.

    EXAMPLES::

        sage: B=matrix([[0,3,-3],[-3,0,3],[3,-3,0]])
        sage: get_uca_element(B,[1,1,0])
        x1^-1 * x0^-1 * x2^3 * (x0^3 + x1^3 + x2^3) 
    
        sage: get_uca_element(B,[1,1,1])
        (2) * x2^2 * x1^2 * x0^2
    
    
        sage: B=matrix([[0,3,0],[-3,0,3],[0,-3,0]])
        sage: get_uca_element(B,[1,1,0])
        x1^-1 * x0^-1 * (x1^3*x2^3 + x0^3 + x2^3)
        sage: get_uca_element(B,[1,1,1])
        x2^-1 * x1^-1 * x0^-1 * (x1 + 1) * (x0 + x2) * (x1^2 - x1 + 1) * (x0^2 - x0*x2 + x2^2)
    """
    #Checks if the length of the
    if len(a) != B.ncols():
        raise ValueError('The length of the input vector must be the same as the number of columns of B.')
    #Runs helper functions.
    v=_vector_decomposition(a,B.nrows())
    c=_compute_compatible_vectors(B,v)
    return _produce_uca_element(B,v,c)

def LLM_gen_set(B):
    r"""
    Produces an list of upper cluster algebra elements corresponding to all vectors in ``\{0,1\}^n``. 

    INPUT::

    - ``B`` -- a skew-symmetric matrix.

    OUTPUT::

    An array of elements in the upper cluster algebra. 

    EXAMPLES::

        sage: B=matrix([[0,1,0],[-1,0,1],[0,-1,0]])
        [1,
        (x1 + 1)/x2,
        (x0 + x2)/x1,
        (x0*x1 + x0 + x2)/(x1*x2),
        (x1 + 1)/x0,
        (x1^2 + 2*x1 + 1)/(x0*x2),
        (x1*x2 + x0 + x2)/(x0*x1),
        (x0*x1 + x1*x2 + x0 + x2)/(x0*x1*x2)]
    """
    aSet=_power_set(B.ncols())
    genSet=[]
    for a in aSet:
        genSet.append(get_uca_element(B,a))
    return (genSet)


def _zero_max(int1):
    r"""
    Returns the max of an integer and zero.

    INPUT::

    - ``int1`` -- an integer.

    OUTPUT::

    The maximum of ``int1`` and zero. 

    EXAMPLES::

        sage: _zero_max(5)
        5
    
        sage: _zero_max(-5)
        0

    """
    return max(int1,0)

def _power_set(n):
    r"""
    Returns an array of all vectors in ``\{0,1\}^n``.

    INPUT::

    - ``n`` -- an integer.

    OUTPUT:: 

    A 2-dimensional array containing all elements of ``\{0,1\}^n``.

    EXAMPLES::

        sage: _power_set(2)
        [[0, 0], [0, 1], [1, 0], [1, 1]]
    
        sage: _power_set(5)
        [[0, 0, 0, 0, 0],
        [0, 0, 0, 0, 1],
        [0, 0, 0, 1, 0],
        [0, 0, 0, 1, 1],
        [0, 0, 1, 0, 0],
        [0, 0, 1, 0, 1],
        [0, 0, 1, 1, 0],
        [0, 0, 1, 1, 1],
        [0, 1, 0, 0, 0],
        [0, 1, 0, 0, 1],
        [0, 1, 0, 1, 0],
        [0, 1, 0, 1, 1],
        [0, 1, 1, 0, 0],
        [0, 1, 1, 0, 1],
        [0, 1, 1, 1, 0],
        [0, 1, 1, 1, 1],
        [1, 0, 0, 0, 0],
        [1, 0, 0, 0, 1],
        [1, 0, 0, 1, 0],
        [1, 0, 0, 1, 1],
        [1, 0, 1, 0, 0],
        [1, 0, 1, 0, 1],
        [1, 0, 1, 1, 0],
        [1, 0, 1, 1, 1],
        [1, 1, 0, 0, 0],
        [1, 1, 0, 0, 1],
        [1, 1, 0, 1, 0],
        [1, 1, 0, 1, 1],
        [1, 1, 1, 0, 0],
        [1, 1, 1, 0, 1],
        [1, 1, 1, 1, 0],
        [1, 1, 1, 1, 1]]

    """
    p=_multi_concatenate([[]],[0,1])
    for i in range(n-1):
        p=_multi_concatenate(p,[0,1])
    return p

def _multi_concatenate(l1,l2):
    r"""
    Each element of ``l2`` gets added to the end of a copy of each array in ``l1``.
    Used to produce the power set.

    INPUT::

    -``l1`` -- a 2-dimensional array.
    -``l2`` -- a single array.

    OUTPUT::

    A 2-dimensional array.
    
    EXAMPLES::

    sage: _car_product([[0,1,2]],[3,4,5])
    [[0, 1, 2, 3], [0, 1, 2, 4], [0, 1, 2, 5]]

    sage: _car_product([[0,1,2],[3,4,5]],[6,7,8])
    [[0, 1, 2, 6],
    [0, 1, 2, 7],
    [0, 1, 2, 8],
    [3, 4, 5, 6],
    [3, 4, 5, 7],
    [3, 4, 5, 8]]   
    """
    plist=[]
    for i in l1:
        for j in l2:
            ip=copy(i)
            ip.append(j)
            plist.append(ip)
    return plist