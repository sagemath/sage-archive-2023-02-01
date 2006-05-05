from sage.rings.all import (bernoulli, sigma, QQ, Integer)

def eisenstein_series_qexp(k, prec=10, K=QQ):
    """
    Return the q-expansion of the weight k Eisenstein series
    to precision prec in the field K.

    INPUT:
        k -- even positive integer
        prec -- nonnegative integer
        K -- a ring in which B_k/(2*k) is invertible

    EXAMPLES:
sage: eisenstein_series_qexp(2,5)
-1/24 + q + 3*q^2 + 4*q^3 + 7*q^4 + O(q^5)
sage: eisenstein_series_qexp(2,0)
O(q^0)

    """
    k = Integer(k)
    if k%2 or k < 2:
        raise ValueError, "k (=%s) must be an even positive integer"%k
    prec = Integer(prec)
    if prec < 0:
        raise ValueError, "prec (=%s) must an even nonnegative integer"%prec
    R = K[['q']]
    coeffs = [-bernoulli(k) / (2*k)] + [sigma(n,k-1) for n in range(1,prec)]
    return R(coeffs, prec)



