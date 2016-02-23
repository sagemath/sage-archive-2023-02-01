r"""
Finding `F[x]`-roots, or modular `F[x]` roots, in polynomials over `F[x][y]`, where `F` is a (finite) field.

This module contains functions for finding two types of `F[x]` roots in a
polynomial over `F[x][y]`, where `F` is a field. Note that if `F` is an infinite
field, then these functions should work, but no measures are taken to limit
coefficient growth. The functions also assume the existence of a root finding
procedure in `F[x]`.

Given a `Q(x,y) \in F[x,y]`, the first type of root are actual `F[x]` roots,
i.e. polynomials `f(x) \in F[x]` such that `Q(x, f(x)) = 0`.

The second type of root are modular roots: given `Q(x,y) \in F[x,y]` and a
precision `d \in \ZZ_+`, then we find all `f(x) \in F[x]` such that `Q(x, f(x))
\equiv 0 \mod x^d`. Since this set is infinite, we return a succinct description
of all such polynomials: this is given as pairs `(f(x), h)` such that `f(x) +
g(x)x^h` is a modular root of `Q` for any `g \in F[x]`.


AUTHORS:

- Johan S. R. Nielsen, original implementation (see [Nielsen]_ for details)
- David Lucas, ported the original implementation in Sage
"""

#*****************************************************************************
#       Copyright (C) 2015 David Lucas <david.lucas@inria.fr>
#                     2015 Johan S. R. Nielsen <jsrn@jsrn.dk>
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 2 of the License, or
# (at your option) any later version.
#                  http://www.gnu.org/licenses/
#*****************************************************************************


from sage.rings.infinity import infinity
from sage.functions.other import binomial, floor

def _convert_Q_representation(Q):
    r"""
    Converts the bivariate polynomial ``Q(x,y)`` from any accepted type to the
    internally used `F[x]` list, for a field `F`.

    INPUT:

    - ``Q`` -- a bivariate polynomial, represented either over `F[x,y]`, `F[x][y]` or `F[x]` list.

    EXAMPLES::

        sage: from sage.coding.guruswami_sudan.rootfinding import _convert_Q_representation
        sage: F = GF(17)
        sage: Px.<x> = F[]
        sage: Pxy.<y> = Px[]
        sage: Q1 = (y - (x**2 + x + 1)) * (y**2 - x + 1) * (y - (x**3 + 4*x + 16))
        sage: _convert_Q_representation(Q1)
        [16*x^6 + 13*x^4 + 2*x^3 + 4*x + 16,
         x^4 + 4*x^2 + 12*x,
         x^5 + x^4 + 5*x^3 + 3*x^2 + 2*x,
         16*x^3 + 16*x^2 + 12*x,
         1]
    """
    if isinstance(Q, list):
        Rx = Q[0].parent()
        if not hasattr(Rx,'gen'):
            raise ValueError("Q must be given as F[x][y], F[x,y] or as F[x] list.")
        return Q
    else:
        # Find out if Q is in F[x,y] or F[x][y]
        Qorig = Q
        Ryx = Q.parent()
        #TODO: Check Ryx is a polynomial ring over a field
        if len(Ryx.gens())==1:
            # Ok, Q is in F[x][y]
            pass
        elif len(Ryx.gens())==2:
            F = Ryx.base_ring()
            (xs,ys) = Ryx.variable_names()
            Rx = F[xs]
            Ryx = Rx[ys]
            x, y = Rx.gen(), Ryx.gen()
            Q = Ryx(Q)
        else:
            raise ValueError("Q must be given as F[x][y], F[x,y] or as F[x] list.")
        # Then make sure Q is a list of F[x] elements
        return Q.list()

def _sanitise_rootfinding_input(Q, maxd, precision):
    r"""
    Verifies, converts and sanitises legal input to the root-finding procedures,
    as well as returning relevant helper variables.

    INPUT:

    - ``Q`` -- a bivariate polynomial, represented either over `F[x,y]`, `F[x][y]` or `F[x]` list.

    - ``maxd``, an integer, the maximal degree of a root of ``Q`` that we're
      interested in, possibly ``None``.

    - ``precision``, an integer, the precision asked for all monomials of ``Q``, possibly ``None``.

    OUTPUT:

    - ``Q``, a modified version of ``Q``, where all monomials have been
      truncated to ``precision``. Represented as an `F[x]` list.

    - ``Qinp``,  the original ``Q`` passed in input, represented as an `F[x]` list.

    - ``F``, the base ring of the coefficients in ``Q``'s first variable.

    - ``Rx``, the polynomial ring `F[x]`.

    - ``x``, the generator of ``Rx``.

    - ``maxd``, the maximal degree of a root of ``Q`` that we're interested in,
      possibly inferred according ``precision``.

    EXAMPLES::

        sage: from sage.coding.guruswami_sudan.rootfinding import _sanitise_rootfinding_input
        sage: from sage.coding.guruswami_sudan.rootfinding import _convert_Q_representation
        sage: F = GF(17)
        sage: Px.<x> = F[]
        sage: Pxy.<y> = Px[]
        sage: Q = (y - (x**2 + x + 1)) * (y**2 - x + 1) * (y - (x**3 + 4*x + 16))
        sage: Q = _convert_Q_representation(Q)
        sage: _sanitise_rootfinding_input(Q, None, None)
        ([16*x^6 + 13*x^4 + 2*x^3 + 4*x + 16,
          x^4 + 4*x^2 + 12*x,
          x^5 + x^4 + 5*x^3 + 3*x^2 + 2*x,
          16*x^3 + 16*x^2 + 12*x,
          1],
         [16*x^6 + 13*x^4 + 2*x^3 + 4*x + 16,
          x^4 + 4*x^2 + 12*x,
          x^5 + x^4 + 5*x^3 + 3*x^2 + 2*x,
          16*x^3 + 16*x^2 + 12*x,
          1],
         Finite Field of size 17,
         Univariate Polynomial Ring in x over Finite Field of size 17,
         x,
         3)
    """
    Q = _convert_Q_representation(Q)
    Qinp = Q
    Rx = Q[0].parent()
    F = Rx.base_ring()
    x = Rx.gen()

    if not maxd:
        if precision:
            maxd = precision-1
        else:
            #The maximal degree of a root is at most
            #(di-dl)/(l-i) d being the degree of a monomial
            maxd = 0
            l = len(Q) -1
            dl = Q[l].degree()
            for i in range(l):
                qi = Q[i]
                if not qi.is_zero():
                    tmp = floor((qi.degree() - dl) / (l - i))
                    if tmp > maxd:
                        maxd = tmp
    if precision:
        for t in range(len(Q)):
            if Q[t].degree >= precision:
                Q[t] = Q[t].truncate(precision)
    return (Q, Qinp, F, Rx, x, maxd)

def _strip_x_pows(Q):
    r"""
    Returns ``(Q', s)`` where ``Q'`` is ``Q`` whose all elements
    have been divided by the largest power ``s`` of ``x`` possible
    such that all these elements remain polynomials.

    INPUT:

    - ``Q`` -- a bivariate polynomial as a list of its monomial in its first variable

    OUTPUT:

    - ``(Q', s)`` a list of two elements:

        - ``Q'``, a polynomial and
        - ``s``, an integer.

    EXAMPLES::

        sage: from sage.coding.guruswami_sudan.rootfinding import _strip_x_pows
        sage: F = GF(17)
        sage: Px.<x> = F[]
        sage: Pxy.<y> = Px[]
        sage: Q1 = (y - (x**2 + x + 1)) * (y**2 - x + 1) * (y - (x**3 + 4*x + 16))
        sage: _strip_x_pows(Q1)
        (y^4 + (16*x^3 + 16*x^2 + 12*x)*y^3 + (x^5 + x^4 + 5*x^3 + 3*x^2 + 2*x)*y^2 + (x^4 + 4*x^2 + 12*x)*y + 16*x^6 + 13*x^4 + 2*x^3 + 4*x + 16, 0)
    """
    def lead_zeroes(p):
        if p.is_zero():
            return infinity
        i = 0
        while p[i].is_zero():
            i+=1
        return i
    strip = min([lead_zeroes(p) for p in Q])
    if strip == 0:
        return (Q, 0)
    if strip == infinity:
        return ([ Q[0].parent().zero() ], infinity)
    return ([ p.shift(-strip) for p in Q ] , strip)

def _roth_ruckenstein_i(Q, F, Rx, x, maxd, precision):
    r"""
    Returns all polynomials which are a solution to the root-finding problem


    This is the core of Roth-Ruckenstein's algorithm where all conversion,
    checks and parent-extraction, is being processed.

    INPUT:

    - ``Q`` -- a bivariate polynomial given as a list of its monomials
      in its first variable

    - ``F``, the base ring of the coefficients in ``Q``'s first variable,

    - ``Rx``, the polynomial ring where live all monomial in ``Q``'s first variable,

    - ``x``, the generator of ``Rx``,

    - ``maxd``, the maximal degree of a root of ``Q`` that we're interested in,

    - ``precision``, an integer, the precision asked for all monomials of ``Q``.

    OUTPUT:

    - a list, containing all suitable polynomials

    EXAMPLES::

        sage: from sage.coding.guruswami_sudan.rootfinding import _sanitise_rootfinding_input
        sage: from sage.coding.guruswami_sudan.rootfinding import _convert_Q_representation
        sage: from sage.coding.guruswami_sudan.rootfinding import _roth_ruckenstein_i
        sage: F = GF(17)
        sage: Px.<x> = F[]
        sage: Pxy.<y> = Px[]
        sage: Q = (y - (x**2 + x + 1)) * (y**2 - x + 1) * (y - (x**3 + 4*x + 16))
        sage: Q = _convert_Q_representation(Q)
        sage: res = _sanitise_rootfinding_input(Q, None, None)
        sage: _roth_ruckenstein_i(res[0], res[2], res[3], res[4], res[5], None)
        [x^3 + 4*x + 16, x^2 + x + 1]
    """
    solutions = []
    g = [F.zero()] * (maxd+1)

    def roth_rec(Q, lam, k):
        r"""
        Recursion of the root finding:
        Q is the remaining poly, lam is the power of x whose coefficient we are
        to determine now, and k is the remaining precision to handle (if ``precision`` is given)
        """
        if precision and k <= 0:
            solutions.append((Rx(g[:lam]), lam))
            return
        (T, strip) = _strip_x_pows(Q)
        if precision:
            k = k - strip
        Ty = Rx([ p[0] for p in T ])
        if Ty.is_zero() or (precision and k <= 0):
            if precision:
                solutions.append((Rx(g[:lam]), lam))
            else:
                assert all(p.is_zero() for p in Q) , ("Q is not zero but Ty is?:\nQ = %s" % Q)
                solutions.append(Rx(g[:lam]))
            return
        roots = Ty.roots(multiplicities=False)
        for gamma in roots:
            g[lam] = gamma
            if lam<maxd:
                # Construct T(y=x*y + gamma)
                ell = len(T)-1
                yc = [[binomial(s, t) * x**t * gamma**(s-t) for t in range(s+1)] for s in range(ell+1)]
                Tg = []
                for t in range(ell+1):
                    Tg.append(sum(yc[s][t] * T[s] for s in range(t, ell+1)))
                roth_rec(Tg , lam+1, k)
            else:
                if precision:
                    solutions.append((Rx(g[:lam+1]), lam+1))
                elif sum( Q[t] * gamma**t for t in range(len(Q)) ).is_zero():
                    solutions.append(Rx(g[:lam+1]))
        return
    roth_rec(Q, 0, precision)
    return solutions

def rootfind_roth_ruckenstein(Q, maxd=None, precision=None):
    r"""
    Returns the list of roots of a bivariate polynomial ``Q``.

    Uses the Roth-Ruckenstein algorithm to find roots or roots
    modulo-up-to-some-precision of a `Q \in \mathbb{F}[x][y]` where `\mathbb{F}` is a field.

    If ``precision = None`` then actual roots will be found, i.e. all `f \in \mathbb{F}[x]`
    such that `Q(f) = 0`. This will be returned as a list of `\mathbb{F}[x]` elements.

    If ``precision = k`` for some integer ``k``, then all `f \in \mathbb{F}[x]` such that
    `Q(f) \equiv 0 \mod x^k` will be returned. This set is infinite, and so it
    will be returned as a list of pairs in `\mathbb{F}[x] \times \mathbb{Z}_+`, where
    `(f, d)` denotes that `Q(f + x^d h) \equiv 0 \mod x^k` for any `h \in
    \mathbb{F}[x]`.

    If ``maxd`` is given, then find only `f` with `deg f \leq maxd`. In case
    `precision=k` setting `maxd` means to only find the roots up to precision
    `maxd`; otherwise, the precision will be `precision-1`.

    INPUT:

    - ``Q`` -- a bivariate polynomial,

    - ``maxd`` -- (default: ``None``) an integer degree bound, as defined above, and

    - ``precision`` -- (default: ``None``) an integer, as defined above.

    EXAMPLES::

        sage: from sage.coding.guruswami_sudan.rootfinding import rootfind_roth_ruckenstein
        sage: F = GF(17)
        sage: Px.<x> = F[]
        sage: Pxy.<y> = Px[]
        sage: Q = (y - (x**2 + x + 1)) * (y**2 - x + 1) * (y - (x**3 + 4*x + 16))
        sage: rootfind_roth_ruckenstein(Q, None, None)
        [x^3 + 4*x + 16, x^2 + x + 1]
    """
    (Q, Qinp, F, Rx, x, maxd) = _sanitise_rootfinding_input(Q, maxd, precision)
    if all(p.is_zero() for p in Q):
        if precision:
            return [(Rx.zero(), 0)]
        else:
            return ValueError("The zero polynomial has infinitely many roots.")
    return _roth_ruckenstein_i(Q, F, Rx, x, maxd, precision)
