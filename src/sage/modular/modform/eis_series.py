from sage.rings.all import (bernoulli, sigma, QQ, Integer)

def eisenstein_series(k, prec):
    k = Integer(k)
    R = QQ[['q']]
    coeffs = [-bernoulli(k) / (2*k)] + [sigma(n,k-1) for n in range(1,prec)]
    return R(coeffs, prec)


