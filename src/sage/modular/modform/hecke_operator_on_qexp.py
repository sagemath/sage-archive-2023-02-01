"""
Hecke Operators on $q$-expansions.
"""

#########################################################################
#       Copyright (C) 2004--2006 William Stein <wstein@gmail.com>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#
#                  http://www.gnu.org/licenses/
#########################################################################

from sage.modular.dirichlet import DirichletGroup, is_DirichletCharacter
from sage.rings.all import (divisors, gcd, ZZ, Integer, is_PowerSeries, Infinity)
from sage.matrix.all import matrix, MatrixSpace
from element import is_ModularFormElement

def hecke_operator_on_qexp(f, n, k, eps = None,
                           prec=None, check=True, _return_list=False):
    r"""
    Given the `q`-expansion `f` of a modular form with character
    `\varepsilon`, this function computes the image of `f` under the
    Hecke operator `T_{n,k}` of weight `k`.

    EXAMPLES::

        sage: M = ModularForms(1,12)
        sage: hecke_operator_on_qexp(M.basis()[0], 3, 12)
        252*q - 6048*q^2 + 63504*q^3 - 370944*q^4 + O(q^5)
        sage: hecke_operator_on_qexp(M.basis()[0], 1, 12, prec=7)
        q - 24*q^2 + 252*q^3 - 1472*q^4 + 4830*q^5 - 6048*q^6 + O(q^7)
        sage: hecke_operator_on_qexp(M.basis()[0], 1, 12)
        q - 24*q^2 + 252*q^3 - 1472*q^4 + 4830*q^5 - 6048*q^6 - 16744*q^7 + 84480*q^8 - 113643*q^9 - 115920*q^10 + 534612*q^11 - 370944*q^12 - 577738*q^13 + O(q^14)

        sage: M.prec(20)
        20
        sage: hecke_operator_on_qexp(M.basis()[0], 3, 12)
        252*q - 6048*q^2 + 63504*q^3 - 370944*q^4 + 1217160*q^5 - 1524096*q^6 + O(q^7)
        sage: hecke_operator_on_qexp(M.basis()[0], 1, 12)
        q - 24*q^2 + 252*q^3 - 1472*q^4 + 4830*q^5 - 6048*q^6 - 16744*q^7 + 84480*q^8 - 113643*q^9 - 115920*q^10 + 534612*q^11 - 370944*q^12 - 577738*q^13 + 401856*q^14 + 1217160*q^15 + 987136*q^16 - 6905934*q^17 + 2727432*q^18 + 10661420*q^19 - 7109760*q^20 + O(q^21)

        sage: (hecke_operator_on_qexp(M.basis()[0], 1, 12)*252).add_bigoh(7)
        252*q - 6048*q^2 + 63504*q^3 - 370944*q^4 + 1217160*q^5 - 1524096*q^6 + O(q^7)

        sage: hecke_operator_on_qexp(M.basis()[0], 6, 12)
        -6048*q + 145152*q^2 - 1524096*q^3 + O(q^4)

    An example on a formal power series::

        sage: R.<q> = QQ[[]]
        sage: f = q + q^2 + q^3 + q^7 + O(q^8)
        sage: hecke_operator_on_qexp(f, 3, 12)
        q + O(q^3)
        sage: hecke_operator_on_qexp(delta_qexp(24), 3, 12).prec()
        8
        sage: hecke_operator_on_qexp(delta_qexp(25), 3, 12).prec()
        9
    """
    if eps is None:
        # Need to have base_ring=ZZ to work over finite fields, since
        # ZZ can coerce to GF(p), but QQ can't.
        eps = DirichletGroup(1, base_ring=ZZ).gen(0)
    if check:
        if not (is_PowerSeries(f) or is_ModularFormElement(f)):
            raise TypeError, "f (=%s) must be a power series or modular form"%f
        if not is_DirichletCharacter(eps):
            raise TypeError, "eps (=%s) must be a Dirichlet character"%eps
        k = Integer(k)
        n = Integer(n)
    v = []

    if prec is None:
        if is_ModularFormElement(f):
            # always want at least three coeffs, but not too many, unless
            # requested
            pr = max(f.prec(), f.parent().prec(), (n+1)*3)
            pr = min(pr, 100*(n+1))
            prec = pr // n + 1
        else:
            prec = (f.prec() / ZZ(n)).ceil()
            if prec == Infinity: prec = f.parent().default_prec() // n + 1

    if f.prec() < prec:
        f._compute_q_expansion(prec)

    l = k-1
    for m in range(prec):
        am = sum([eps(d) * d**l * f[m*n//(d*d)] for \
                  d in divisors(gcd(n, m)) if (m*n) % (d*d) == 0])
        v.append(am)
    if _return_list:
        return v
    if is_ModularFormElement(f):
        R = f.parent()._q_expansion_ring()
    else:
        R = f.parent()
    return R(v, prec)

def _hecke_operator_on_basis(B, V, n, k, eps):
    """
    Does the work for hecke_operator_on_basis once the input
    is normalized.

    EXAMPLES::

        sage: hecke_operator_on_basis(ModularForms(1,16).q_expansion_basis(30), 3, 16) # indirect doctest
        [   -3348        0]
        [       0 14348908]

    The following used to cause a segfault due to acidentally
    transposed second and third argument (#2107)::

        sage: B = victor_miller_basis(100,30)
        sage: t2 = hecke_operator_on_basis(B, 100, 2)
        Traceback (most recent call last):
        ...
        ValueError: The given basis vectors must be linearly independent.
    """
    prec = V.degree()
    TB = [hecke_operator_on_qexp(f, n, k, eps, prec, check=False, _return_list=True)
                for f in B]
    TB = [V.coordinate_vector(w) for w in TB]
    return matrix(V.base_ring(), len(B), len(B), TB, sparse=False)

def hecke_operator_on_basis(B, n, k, eps=None,
                            already_echelonized = False):
    r"""
    Given a basis `B` of `q`-expansions for a space of modular forms
    with character `\varepsilon` to precision at least `\#B\cdot n+1`,
    this function computes the matrix of `T_n` relative to `B`.

    .. note::

       If the elements of B are not known to sufficient precision,
       this function will report that the vectors are linearly
       dependent (since they are to the specified precision).

    INPUT:

    - ``B`` - list of q-expansions

    - ``n`` - an integer >= 1

    - ``k`` - an integer

    - ``eps`` - Dirichlet character

    - ``already_echelonized`` -- bool (default: False); if True, use that the
      basis is already in Echelon form, which saves a lot of time.

    EXAMPLES::

        sage: sage.modular.modform.constructor.ModularForms_clear_cache()
        sage: ModularForms(1,12).q_expansion_basis()
        [
        q - 24*q^2 + 252*q^3 - 1472*q^4 + 4830*q^5 + O(q^6),
        1 + 65520/691*q + 134250480/691*q^2 + 11606736960/691*q^3 + 274945048560/691*q^4 + 3199218815520/691*q^5 + O(q^6)
        ]
        sage: hecke_operator_on_basis(ModularForms(1,12).q_expansion_basis(), 3, 12)
        Traceback (most recent call last):
        ...
        ValueError: The given basis vectors must be linearly independent.

        sage: hecke_operator_on_basis(ModularForms(1,12).q_expansion_basis(30), 3, 12)
        [   252      0]
        [     0 177148]
    """
    if eps is None:
        eps = DirichletGroup(1).gen(0)
    if not isinstance(B, (list, tuple)):
        raise TypeError, "B (=%s) must be a list or tuple"%B
    if len(B) == 0:
        return MatrixSpace(eps.base_ring(),0)(0)
    all_powerseries = True
    for x in B:
        if not is_PowerSeries(x):
            all_powerseries = False
    if not all_powerseries:
        raise TypeError, "each element of B must be a power series"
    f = B[0]
    n = Integer(n)
    k = Integer(k)
    R = f.base_ring()
    prec = (f.prec()-1)//n
    A = R**prec
    V = A.span_of_basis([g.padded_list(prec) for g in B],
                        already_echelonized = already_echelonized)
    return _hecke_operator_on_basis(B, V, n, k, eps)


