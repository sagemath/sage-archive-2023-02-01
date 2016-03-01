r"""
Finding `F[x]`-roots for polynomials over `F[x][y]`, with`F` is a (finite) field, as used in the Guruswami-Sudan decoding.

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
        ([16*x^6 + 13*x^4 + 2*x^3 + 4*x + 16,
         x^4 + 4*x^2 + 12*x,
         x^5 + x^4 + 5*x^3 + 3*x^2 + 2*x,
         16*x^3 + 16*x^2 + 12*x,
         1], Univariate Polynomial Ring in x over Finite Field of size 17)
    """
    if isinstance(Q, list):
        if Q == []:
            return ([], None)
        Rx = Q[0].parent()
        if not hasattr(Rx,'gen'):
            raise ValueError("Q must be given as F[x][y], F[x,y] or as F[x] list.")
        return (Q, Rx)
    else:
        # Find out if Q is in F[x,y] or F[x][y]
        Qorig = Q
        Ryx = Q.parent()
        #TODO: Check Ryx is a polynomial ring over a field
        if len(Ryx.gens())==1:
            # Ok, Q is in F[x][y]
            Rx = Q.base_ring()
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
        return (Q.list(), Rx)

def _sanitise_rootfinding_input(Q, maxd, precision):
    r"""
    Verifies, converts and sanitises legal input to the root-finding procedures,
    as well as returning relevant helper variables.

    INPUT:

    - ``Q`` -- a bivariate polynomial, represented either over `F[x,y]`, `F[x][y]` or `F[x]` list.

    - ``maxd``, an integer, the maximal degree of a root of ``Q`` that we're
      interested in, possibly ``None``.

    - ``precision``, an integer, the precision asked for modular roots of ``Q``, possibly ``None``.

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
        sage: F = GF(17)
        sage: Px.<x> = F[]
        sage: Pxy.<y> = Px[]
        sage: Q = (y - (x**2 + x + 1)) * (y**2 - x + 1) * (y - (x**3 + 4*x + 16))
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
    (Q, Rx) = _convert_Q_representation(Q)
    if Q == []:
        return ([],[],None,Rx,None,0) # Q == 0 so just bail
    Qinp = Q
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
    Returns ``(Q', s)`` where ``Q'`` is ``Q`` whose elements have all been
    divided by the largest power ``s`` of ``x`` possible such that all these
    elements remain polynomials.

    INPUT:

    - ``Q`` -- a bivariate polynomial in `F[x][y]` as a list of `F[x]` polynomials.

    OUTPUT:

    - ``(Q', s)`` a list of two elements:

        - ``Q'``, the reduced bivariate polynomial, as a list of univariate ones.
        - ``s``, an integer, the power of `x` that was stripped from `Q`.

    EXAMPLES::

        sage: from sage.coding.guruswami_sudan.rootfinding import _strip_x_pows
        sage: F = GF(17)
        sage: Px.<x> = F[]
        sage: Q = [ Px(3*x^2 + 2*x),  Px(5*x^7 + x^6) ]
        sage: _strip_x_pows(Q)
        ([3*x + 2, 5*x^6 + x^5], 1)
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
    Returns all polynomials which are a solution to the, possibly modular,
    root-finding problem.

    This is the core of Roth-Ruckenstein's algorithm where all conversions,
    checks and parent-extraction have been done. Most of the inputs corresponds
    to the output of ``_sanitise_rootfinding_input``.

    INPUT::

    - ``Q``, a modified version of ``Q``, where all monomials have been
      truncated to ``precision``. Represented as an `F[x]` list.

    - ``Qinp``,  the original ``Q`` passed in input, represented as an `F[x]` list.

    - ``F``, the base ring of the coefficients in ``Q``'s first variable.

    - ``Rx``, the polynomial ring `F[x]`.

    - ``x``, the generator of ``Rx``.

    - ``maxd``, the maximal degree of a root of ``Q`` that we're interested in,
      possibly inferred according ``precision``.

    - ``precision``, a non-negative integer or `None`. If given, it is the
      sought precision for modular roots of `Q`. Otherwise, we will find
      unconditional roots.

    OUTPUT::

    - a list, containing all `F[x]` roots of `Q(x,y)`, possibly modular. If
    ``precision`` is given, we return a list of pairs `(f, h)`, where `f \in
    F[x]` and `h` is a non-negative integer, and where `f + h*g \equiv 0 \mod
    x^{d}` for any `g \in F[x]`, and where `d` is ``precision``.

    EXAMPLES::

        sage: from sage.coding.guruswami_sudan.rootfinding import _sanitise_rootfinding_input
        sage: from sage.coding.guruswami_sudan.rootfinding import _roth_ruckenstein_i
        sage: F = GF(17)
        sage: Px.<x> = F[]
        sage: Pxy.<y> = Px[]
        sage: Q = (y - (x**2 + x + 1)) * (y**2 - x + 1) * (y - (x**3 + 4*x + 16))
        sage: (Q, Qinp, F, Rx, x, maxd) = _sanitise_rootfinding_input(Q, None, None)
        sage: set(_roth_ruckenstein_i(Q, F, Rx, x, maxd, None))
        {x^2 + x + 1, x^3 + 4*x + 16}
        sage: set(_roth_ruckenstein_i(Q, F, Rx, x, maxd, precision=2))
        {(x + 1, 2), (2*x + 13, 2), (4*x + 16, 2), (15*x + 4, 2)}
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

    Uses the Roth-Ruckenstein algorithm to find roots or modular roots of a `Q
    \in \mathbb{F}[x][y]` where `\mathbb{F}` is a field.

    If ``precision = None`` then actual roots will be found, i.e. all `f \in
    \mathbb{F}[x]` such that `Q(f) = 0`. This will be returned as a list of
    `\mathbb{F}[x]` elements.

    If ``precision = d`` for some integer ``d``, then all `f \in \mathbb{F}[x]`
    such that `Q(f) \equiv 0 \mod x^d` will be returned. This set is infinite,
    and so it will be returned as a list of pairs in `\mathbb{F}[x] \times
    \mathbb{Z}_+`, where `(f, h)` denotes that `Q(f + x^h g) \equiv 0 \mod x^d`
    for any `g \in \mathbb{F}[x]`.

    If ``maxd`` is given, then find only `f` with `deg f \leq maxd`. In case
    `precision=d` setting `maxd` means to only find the roots up to precision
    `maxd`, i.e. `h \leq maxd` in the above; otherwise, this will be naturally
    bounded at `precision-1`.

    INPUT:

    - ``Q`` -- a bivariate polynomial, represented either over `F[x,y]`, `F[x][y]` or `F[x]` list.

    - ``maxd`` -- (default: ``None``) an non-negative integer degree bound, as defined above.

    - ``precision`` -- (default: ``None``) an integer, as defined above.

    EXAMPLES::

        sage: from sage.coding.guruswami_sudan.rootfinding import rootfind_roth_ruckenstein
        sage: F = GF(17)
        sage: Px.<x> = F[]
        sage: Py.<y> = Px[]
        sage: Q = (y - (x**2 + x + 1)) * (y**2 - x + 1) * (y - (x**3 + 4*x + 16))
        sage: roots = rootfind_roth_ruckenstein(Q); set(roots)
        {x^2 + x + 1, x^3 + 4*x + 16}
        sage: Q(roots[0])
        0
        sage: set(rootfind_roth_ruckenstein(Q, maxd = 2))
        {x^2 + x + 1}
        sage: modroots = rootfind_roth_ruckenstein(Q, precision = 3); set(modroots)
        {(4*x + 16, 3), (x^2 + x + 1, 3), (8*x^2 + 15*x + 4, 3), (9*x^2 + 2*x + 13, 3)}
        sage: (f,h) = modroots[0]
        sage: Q(f + x^h * Px.random_element()) % x^3
        0
        sage: modroots2 = rootfind_roth_ruckenstein(Q, maxd=1, precision = 3); set(modroots2)
        {(x + 1, 2), (2*x + 13, 2), (4*x + 16, 2), (15*x + 4, 2)}

    TESTS:

    Test that if `Q = 0`, then the appropriate response is given

        sage: F = GF(17)
        sage: R.<x,y> = F[]
        sage: rootfind_roth_ruckenstein(R.zero())
        ValueError('The zero polynomial has infinitely many roots.',)
        sage: rootfind_roth_ruckenstein(R.zero(), precision=1)
        [(0, 0)]

    """
    (Q, Qinp, F, Rx, x, maxd) = _sanitise_rootfinding_input(Q, maxd, precision)
    if all(p.is_zero() for p in Q):
        if precision:
            return [(Rx.zero() if hasattr(Rx,'zero') else 0, 0)]
        else:
            return ValueError("The zero polynomial has infinitely many roots.")
    return _roth_ruckenstein_i(Q, F, Rx, x, maxd, precision)
