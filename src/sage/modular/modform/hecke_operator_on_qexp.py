#########################################################################
#       Copyright (C) 2004--2006 William Stein <wstein@gmail.com>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#
#                  http://www.gnu.org/licenses/
#########################################################################

from sage.modular.dirichlet import DirichletGroup, is_DirichletCharacter
from sage.rings.all import (divisors, infinity, gcd, Integer,
                            is_PowerSeries)
from sage.matrix.all import matrix, MatrixSpace

def hecke_operator_on_qexp(f, n, k, eps = None,
                           prec=None, check=True, _return_list=False):
    r"""
    Given the $q$-expansion $f$ of a modular form with character
    $\varepsilon$, this function computes the Hecke operator $T_{n,k}$
    of weight~$k$ on~$f$.
    """
    if eps is None:
        eps = DirichletGroup(1).gen(0)
    if not check:
        if not is_PowerSeries(f):
            raise TypeError, "f (=%s) must be a power series"%f
        if not is_DirichletCharacter(eps):
            raise TypeError, "eps (=%s) must be a Dirichlet character"%eps
        k = Integer(k)
        n = Integer(n)
    v = []
    if prec is None:
        pr = f.prec()
        if pr is infinity:
            raise ValueError, "f must have finite precision."
        prec = pr // n + 1
    else:
        if f.prec() < prec:
            raise ValueError, "f is not known to sufficient precision."
    l = k-1
    for m in range(prec):
        am = sum([eps(d) * d**l * f[m*n//(d*d)] for \
                  d in divisors(gcd(n, m)) if (m*n) % (d*d) == 0])
        v.append(am)
    if _return_list:
        return v
    R = f.parent()
    return R(v, prec)


def _hecke_operator_on_basis(B, V, n, k, eps):
    prec = V.degree()
    TB = [hecke_operator_on_qexp(f, n, k, eps, prec, check=False, _return_list=True)
                for f in B]
    TB = [V.coordinate_vector(w) for w in TB]
    return matrix(V.base_ring(), len(B), len(B), TB, sparse=False)

def hecke_operator_on_basis(B, n, k, eps=None,
                            already_echelonized = False):
    r"""
    Given a basis $B$ of $q$-expansions for a space of modular forms
    with character $\varepsilon$ to precision at least $\#B\cdot n+1$,
    this function computes the matrix of $T_n$ relative to $B$.

    INPUT:
        B -- list of q-expansions
        n -- an integer >= 1
        k -- an integer
        eps -- Dirichlet character
        already_echelonized -- bool (default: False); if True, use that the
                basis is already in Echelon form, which saves a lot of time.
    """
    if eps is None:
        eps = DirichletGroup(1).gen(0)
    if not isinstance(B, (list, tuple)):
        raise TypeError, "B (=%s) must be a list or tuple"%B
    if len(B) == 0:
        return MatrixSpace(eps.base_ring(),0)(0)
    f = B[0]
    if not is_PowerSeries(f):
        raise TypeError, "each element of B must be a power series"
    n = Integer(n)
    k = Integer(k)
    R = f.base_ring()
    prec = (f.prec()-1)//n
    A = R**prec
    V = A.span_of_basis([g.padded_list(prec) for g in B],
                        already_echelonized = already_echelonized)
    return _hecke_operator_on_basis(B, V, n, k, eps)


