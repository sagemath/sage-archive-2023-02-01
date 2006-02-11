from sage.libs._pari.functional import *

init_primes(2000000)

def pith(x):
    if x <= 2:
        return x == 2
    lx = log(x)
    l = floor(x/lx)
    r = floor(l * (1 + 3/(2*lx)))
    while r-l > 1:
        m = (l+r) >> 1
        if prime(m)<=x:
            l=m
        else:
            r=m
    return l;

set_series_precision(30)

s = variable('s')
ellld_ll = -Euler()*s + sum([(-1)**n * zeta(n) / n * s**n for n in range(2,31)])
