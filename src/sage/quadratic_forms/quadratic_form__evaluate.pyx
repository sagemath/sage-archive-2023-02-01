

def QFEvaluateVector(Q, v):
    """
    Python wrapper for our fast evaluation routine.
    """
    return QFEvaluateVector_cdef(Q, v)



cdef QFEvaluateVector_cdef(Q, v):
    """
    Evaluate this quadratic form Q on a vector or matrix of elements
    coercible to the base ring of the quadratic form.  If a vector
    is given then the output will be the ring element Q(v), but if a
    matrix is given then the output will be the quadratic form Q'
    which in matrix notation is given by:

            Q' = v^t * Q * v.

    """
    ## If we are passed a matrix A, return the quadratic form Q(A(x))
    ## (In matrix notation: A^t * Q * A)
    n = Q.dim()

    tmp_val = Q.base_ring()(0)
    for i from 0 <= i < n:
        for j from i <= j < n:
            tmp_val += Q[i,j] * v[i] * v[j]

    ## Return the value (over R)
    return Q.base_ring()._coerce_(tmp_val)



def QFEvaluateMatrix(Q, v, Q2):
    """
    Python wrapper for our fast evaluation routine.
    """
    return QFEvaluateMatrix_cdef(Q, v, Q2)



cdef QFEvaluateMatrix_cdef(Q, v, Q2):
    """
    Fast evaluation on a matrix.
    """
    ## Create the new quadratic form
    n = Q.dim()
    m = Q2.dim()
    for k from 0 <= k < m:
        for l from k <= l < m:
            tmp_sum = Q2.base_ring()(0)
            for i from 0 <= i < n:
                for j from i <= j < n:
                    if (k == l):
                        tmp_sum += Q[i,j] * (v[i,k] * v[j,l])
                    else:
                        tmp_sum += Q[i,j] * (v[i,k] * v[j,l] + v[i,l] * v[j,k])
            Q2[k,l] = tmp_sum
    return Q2
