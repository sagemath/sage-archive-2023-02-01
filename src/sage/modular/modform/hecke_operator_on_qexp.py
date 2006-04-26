from sage.modular.dirichlet import DirichletGroup
from sage.rings.all import divisors, infinity, gcd

def hecke_operator_on_qexp(f, n, k, eps = DirichletGroup(1).gen(0)):
    v = []
    pr = f.prec()
    if pr == infinity:
        raise ValueError, "f must have finite precision."
    prec = pr // n + 1
    l = k-1
    for m in range(prec):
        am = sum([eps(d) * d**l * f[m*n//(d*d)] for d in divisors(gcd(n, m)) if (m*n) % (d*d) == 0])
        v.append(am)
    R = f.parent()
    return R(v, prec)
