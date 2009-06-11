"""
Educational Versions of Groebner Basis Algorithms.

Following [BW93]_ the original Buchberger algorithm (c.f. algorithm
GROEBNER in [BW93]_) and an improved version of Buchberger's algorithm
(c.g. algorithm GROEBNERNEW2 in [BW93]_) are implemented.

No attempt was made to optimize either algorithm as the emphasis of
these implementations is a clean and easy presentation. To compute a
Groebner basis in Sage efficently use the ``groebner_basis()``
method on multivariate polynomial objects.


.. note::

   The notion of 'term' and 'monomial' in [BW93]_ is swapped from the
   notion of those words in Sage (or the other way around, however you
   prefer it). In Sage a term is a monomial multiplied by a
   coefficient, while in [BW93]_ a monomial is a term multiplied by a
   coefficient. Also, what is called LM (the leading monomial) in
   Sage is called HT (the head term) in [BW93]_.

EXAMPLES:

Consider Katsura-6 w.r.t. a ``degrevlex`` ordering.::

    sage: from sage.rings.polynomial.toy_buchberger import *
    sage: P.<a,b,c,e,f,g,h,i,j,k> = PolynomialRing(GF(32003),10)
    sage: I = sage.rings.ideal.Katsura(P,6)

    sage: g1 = buchberger(I)
    sage: g2 = buchberger_improved(I)
    sage: g3 = I.groebner_basis()

All algorithms actually compute a Groebner basis::

    sage: Ideal(g1).basis_is_groebner()
    True
    sage: Ideal(g2).basis_is_groebner()
    True
    sage: Ideal(g3).basis_is_groebner()
    True

The results are correct::

    sage: Ideal(g1) == Ideal(g2) == Ideal(g3)
    True

If ``get_verbose()`` is `>= 1` a protocol is provided::

    sage: set_verbose(1)
    sage: P.<a,b,c> = PolynomialRing(GF(127),3)
    sage: I = sage.rings.ideal.Katsura(P)
    // sage...              [0]  ideal, 3 generator(s)

    sage: I
    Ideal (a + 2*b + 2*c - 1, a^2 + 2*b^2 + 2*c^2 - a, 2*a*b + 2*b*c - b) of Multivariate Polynomial Ring in a, b, c over Finite Field of size 127

The original Buchberger algorithm performs 15 useless reductions to
zero for this example::

    sage: buchberger(I)
    (a + 2*b + 2*c - 1, a^2 + 2*b^2 + 2*c^2 - a) => -2*b^2 - 6*b*c - 6*c^2 + b + 2*c
    G: set([a + 2*b + 2*c - 1, 2*a*b + 2*b*c - b, a^2 + 2*b^2 + 2*c^2 - a, -2*b^2 - 6*b*c - 6*c^2 + b + 2*c])
    <BLANKLINE>
    (a^2 + 2*b^2 + 2*c^2 - a, a + 2*b + 2*c - 1) => 0
    G: set([a + 2*b + 2*c - 1, 2*a*b + 2*b*c - b, a^2 + 2*b^2 + 2*c^2 - a, -2*b^2 - 6*b*c - 6*c^2 + b + 2*c])
    <BLANKLINE>
    (a + 2*b + 2*c - 1, 2*a*b + 2*b*c - b) => -5*b*c - 6*c^2 - 63*b + 2*c
    G: set([a + 2*b + 2*c - 1, 2*a*b + 2*b*c - b, -5*b*c - 6*c^2 - 63*b + 2*c, a^2 + 2*b^2 + 2*c^2 - a, -2*b^2 - 6*b*c - 6*c^2 + b + 2*c])
    <BLANKLINE>
    (2*a*b + 2*b*c - b, a + 2*b + 2*c - 1) => 0
    G: set([a + 2*b + 2*c - 1, 2*a*b + 2*b*c - b, -5*b*c - 6*c^2 - 63*b + 2*c, a^2 + 2*b^2 + 2*c^2 - a, -2*b^2 - 6*b*c - 6*c^2 + b + 2*c])
    <BLANKLINE>
    (2*a*b + 2*b*c - b, -5*b*c - 6*c^2 - 63*b + 2*c) => -22*c^3 + 24*c^2 - 60*b - 62*c
    G: set([a + 2*b + 2*c - 1, -22*c^3 + 24*c^2 - 60*b - 62*c, 2*a*b + 2*b*c - b, a^2 + 2*b^2 + 2*c^2 - a, -2*b^2 - 6*b*c - 6*c^2 + b + 2*c, -5*b*c - 6*c^2 - 63*b + 2*c])
    <BLANKLINE>
    (2*a*b + 2*b*c - b, -2*b^2 - 6*b*c - 6*c^2 + b + 2*c) => 0
    G: set([a + 2*b + 2*c - 1, -22*c^3 + 24*c^2 - 60*b - 62*c, 2*a*b + 2*b*c - b, a^2 + 2*b^2 + 2*c^2 - a, -2*b^2 - 6*b*c - 6*c^2 + b + 2*c, -5*b*c - 6*c^2 - 63*b + 2*c])
    <BLANKLINE>
    (2*a*b + 2*b*c - b, a^2 + 2*b^2 + 2*c^2 - a) => 0
    G: set([a + 2*b + 2*c - 1, -22*c^3 + 24*c^2 - 60*b - 62*c, 2*a*b + 2*b*c - b, a^2 + 2*b^2 + 2*c^2 - a, -2*b^2 - 6*b*c - 6*c^2 + b + 2*c, -5*b*c - 6*c^2 - 63*b + 2*c])
    <BLANKLINE>
    (a + 2*b + 2*c - 1, -2*b^2 - 6*b*c - 6*c^2 + b + 2*c) => 0
    G: set([a + 2*b + 2*c - 1, -22*c^3 + 24*c^2 - 60*b - 62*c, 2*a*b + 2*b*c - b, a^2 + 2*b^2 + 2*c^2 - a, -2*b^2 - 6*b*c - 6*c^2 + b + 2*c, -5*b*c - 6*c^2 - 63*b + 2*c])
    <BLANKLINE>
    (a^2 + 2*b^2 + 2*c^2 - a, 2*a*b + 2*b*c - b) => 0
    G: set([a + 2*b + 2*c - 1, -22*c^3 + 24*c^2 - 60*b - 62*c, 2*a*b + 2*b*c - b, a^2 + 2*b^2 + 2*c^2 - a, -2*b^2 - 6*b*c - 6*c^2 + b + 2*c, -5*b*c - 6*c^2 - 63*b + 2*c])
    <BLANKLINE>
    (-2*b^2 - 6*b*c - 6*c^2 + b + 2*c, -5*b*c - 6*c^2 - 63*b + 2*c) => 0
    G: set([a + 2*b + 2*c - 1, -22*c^3 + 24*c^2 - 60*b - 62*c, 2*a*b + 2*b*c - b, a^2 + 2*b^2 + 2*c^2 - a, -2*b^2 - 6*b*c - 6*c^2 + b + 2*c, -5*b*c - 6*c^2 - 63*b + 2*c])
    <BLANKLINE>
    (a + 2*b + 2*c - 1, -5*b*c - 6*c^2 - 63*b + 2*c) => 0
    G: set([a + 2*b + 2*c - 1, -22*c^3 + 24*c^2 - 60*b - 62*c, 2*a*b + 2*b*c - b, a^2 + 2*b^2 + 2*c^2 - a, -2*b^2 - 6*b*c - 6*c^2 + b + 2*c, -5*b*c - 6*c^2 - 63*b + 2*c])
    <BLANKLINE>
    (a^2 + 2*b^2 + 2*c^2 - a, -5*b*c - 6*c^2 - 63*b + 2*c) => 0
    G: set([a + 2*b + 2*c - 1, -22*c^3 + 24*c^2 - 60*b - 62*c, 2*a*b + 2*b*c - b, a^2 + 2*b^2 + 2*c^2 - a, -2*b^2 - 6*b*c - 6*c^2 + b + 2*c, -5*b*c - 6*c^2 - 63*b + 2*c])
    <BLANKLINE>
    (-5*b*c - 6*c^2 - 63*b + 2*c, -22*c^3 + 24*c^2 - 60*b - 62*c) => 0
    G: set([a + 2*b + 2*c - 1, -22*c^3 + 24*c^2 - 60*b - 62*c, 2*a*b + 2*b*c - b, a^2 + 2*b^2 + 2*c^2 - a, -2*b^2 - 6*b*c - 6*c^2 + b + 2*c, -5*b*c - 6*c^2 - 63*b + 2*c])
    <BLANKLINE>
    (a + 2*b + 2*c - 1, -22*c^3 + 24*c^2 - 60*b - 62*c) => 0
    G: set([a + 2*b + 2*c - 1, -22*c^3 + 24*c^2 - 60*b - 62*c, 2*a*b + 2*b*c - b, a^2 + 2*b^2 + 2*c^2 - a, -2*b^2 - 6*b*c - 6*c^2 + b + 2*c, -5*b*c - 6*c^2 - 63*b + 2*c])
    <BLANKLINE>
    (a^2 + 2*b^2 + 2*c^2 - a, -2*b^2 - 6*b*c - 6*c^2 + b + 2*c) => 0
    G: set([a + 2*b + 2*c - 1, -22*c^3 + 24*c^2 - 60*b - 62*c, 2*a*b + 2*b*c - b, a^2 + 2*b^2 + 2*c^2 - a, -2*b^2 - 6*b*c - 6*c^2 + b + 2*c, -5*b*c - 6*c^2 - 63*b + 2*c])
    <BLANKLINE>
    (-2*b^2 - 6*b*c - 6*c^2 + b + 2*c, -22*c^3 + 24*c^2 - 60*b - 62*c) => 0
    G: set([a + 2*b + 2*c - 1, -22*c^3 + 24*c^2 - 60*b - 62*c, 2*a*b + 2*b*c - b, a^2 + 2*b^2 + 2*c^2 - a, -2*b^2 - 6*b*c - 6*c^2 + b + 2*c, -5*b*c - 6*c^2 - 63*b + 2*c])
    <BLANKLINE>
    (2*a*b + 2*b*c - b, -22*c^3 + 24*c^2 - 60*b - 62*c) => 0
    G: set([a + 2*b + 2*c - 1, -22*c^3 + 24*c^2 - 60*b - 62*c, 2*a*b + 2*b*c - b, a^2 + 2*b^2 + 2*c^2 - a, -2*b^2 - 6*b*c - 6*c^2 + b + 2*c, -5*b*c - 6*c^2 - 63*b + 2*c])
    <BLANKLINE>
    (a^2 + 2*b^2 + 2*c^2 - a, -22*c^3 + 24*c^2 - 60*b - 62*c) => 0
    G: set([a + 2*b + 2*c - 1, -22*c^3 + 24*c^2 - 60*b - 62*c, 2*a*b + 2*b*c - b, a^2 + 2*b^2 + 2*c^2 - a, -2*b^2 - 6*b*c - 6*c^2 + b + 2*c, -5*b*c - 6*c^2 - 63*b + 2*c])
    <BLANKLINE>
    15 reductions to zero.
    [a + 2*b + 2*c - 1, -22*c^3 + 24*c^2 - 60*b - 62*c, 2*a*b + 2*b*c - b, a^2 + 2*b^2 + 2*c^2 - a, -2*b^2 - 6*b*c - 6*c^2 + b + 2*c, -5*b*c - 6*c^2 - 63*b + 2*c]

The 'improved' Buchberger algorithm in constrast only performs 3 reductions to zero:

    sage: buchberger_improved(I)
    (b^2 - 26*c^2 - 51*b + 51*c, b*c + 52*c^2 + 38*b + 25*c) => 11*c^3 - 12*c^2 + 30*b + 31*c
    G: set([a + 2*b + 2*c - 1, b^2 - 26*c^2 - 51*b + 51*c, 11*c^3 - 12*c^2 + 30*b + 31*c, b*c + 52*c^2 + 38*b + 25*c])
    <BLANKLINE>
    (11*c^3 - 12*c^2 + 30*b + 31*c, b*c + 52*c^2 + 38*b + 25*c) => 0
    G: set([a + 2*b + 2*c - 1, b^2 - 26*c^2 - 51*b + 51*c, 11*c^3 - 12*c^2 + 30*b + 31*c, b*c + 52*c^2 + 38*b + 25*c])
    <BLANKLINE>
    1 reductions to zero.
    [a + 2*b + 2*c - 1, b^2 - 26*c^2 - 51*b + 51*c, c^3 + 22*c^2 - 55*b + 49*c, b*c + 52*c^2 + 38*b + 25*c]

REFERENCES:

.. [BW93] Thomas Becker and Volker Weispfenning. *Groebner Bases - A
  Computational Approach To Commutative Algebra*. Springer, New York
  1993.

AUTHOR:

- Martin Albrecht (2007-05-24): initial version
"""

from sage.misc.misc import get_verbose
from sage.rings.arith import LCM
from sage.structure.sequence import Sequence

LCM = lambda f,g: f.parent().monomial_lcm(f,g)
LM = lambda f: f.lm()
LT = lambda f: f.lt()
spol = lambda f,g: LCM(LM(f),LM(g)) // LT(f) * f - LCM(LM(f),LM(g)) // LT(g) * g

def buchberger(F):
    """
    The original version of Buchberger's algorithm as presented in
    [BW93]_, page 214.

    INPUT:

    - ``F`` - an ideal in a multivariate polynomial ring

    OUTPUT:
        a Groebner basis for F

    .. note::

       The verbosity of this function may be controlled with a
       ``set_verbose()`` call. Any value >=1 will result in this
       function printing intermediate bases.

    """
    G = set(F.gens())
    B = set(filter(lambda (x,y): x!=y, [(g1,g2) for g1 in G for g2 in G]))

    if get_verbose() >=1:
        reductions_to_zero = 0

    while B!=set():
        g1,g2 = select(B)
        B.remove( (g1,g2) )

        h = spol(g1,g2).reduce(G)
        if h != 0:
            B = B.union( [(g,h) for g in G] )
            G.add( h )

        if get_verbose() >= 1:
            print "(%s, %s) => %s"%(g1, g2, h)
            print "G: %s\n"%(G)
            if h==0:
                reductions_to_zero +=1

    if get_verbose() >= 1:
        print "%d reductions to zero."%(reductions_to_zero)

    return Sequence(G)

def buchberger_improved(F):
    """
    An improved version of Buchberger's algorithm as presented in
    [BW93]_, page 232.

    This variant uses the Gebauer-Moeller Installation to apply
    Buchberger's first and second criterion to avoid useless pairs.

    INPUT:

    - ``F`` - an ideal in a multivariate polynomial ring

    OUTPUT:
        a Groebner basis for F

    .. note::

       The verbosity of this function may be controlled with a
       ``set_verbose()`` call. Any value ``>=1`` will result in this
       function printing intermediate Groebner bases.
    """
    F = inter_reduction(F.gens())

    G = set()
    B = set()

    if get_verbose() >=1:
        reductions_to_zero = 0

    while F != set():
        f = min(F)
        F.remove(f)
        G,B = update(G,B,f)

    while B != set():

        g1,g2 = select(B)
        B.remove((g1,g2))
        h = spol(g1,g2).reduce(G)
        if h!=0: G,B = update(G,B,h)

        if get_verbose() >= 1:
            print "(%s, %s) => %s"%(g1,g2,h)
            print "G: %s\n"%(G)
            if h==0:
                reductions_to_zero +=1

    if get_verbose() >= 1:
        print "%d reductions to zero."%(reductions_to_zero)

    return Sequence(inter_reduction(G))

def update(G,B,h):
    """
    Update ``G`` using the list of critical pairs ``B`` and the
    polynomial ``h`` as presented in [BW93]_, page 230. For this,
    Buchberger's first and second criterion are tested.

    This function implements the Gebauer-Moeller Installation.

    INPUT:

    - ``G`` - an intermediate Groebner basis
    - ``B`` - a list of critical pairs
    - ``h`` - a polynomial

    OUTPUT:
        a tuple of an intermediate Groebner basis and a list of
        critical pairs
    """
    R = h.parent()

    C = set([(h,g) for g in G])
    D = set()

    while C != set():
        (h,g) = C.pop()

        lcm_divides = lambda rhs: R.monomial_divides( LCM(LM(h),LM(rhs[1])), LCM(LM(h),LM(g)))

        if R.monomial_pairwise_prime(LM(h),LM(g)) or \
           (\
               not any( lcm_divides(f) for f in C ) \
               and
               not any( lcm_divides(f) for f in D ) \
            ):
            D.add( (h,g) )

    E = set()

    while D != set():
        (h,g) = D.pop()
        if not R.monomial_pairwise_prime(LM(h),LM(g)):
            E.add( (h,g) )

    B_new = set()

    while B != set():
        g1,g2 = B.pop()
        if not R.monomial_divides( LM(h),  LCM(LM(g1),LM(g2)) ) or \
               R.monomial_lcm(LM(g1),LM( h)) == LCM(LM(g1),LM(g2)) or \
               R.monomial_lcm(LM( h),LM(g2)) == LCM(LM(g1),LM(g2)) :
            B_new.add( (g1,g2) )

    B_new = B_new.union( E )

    G_new = set()

    while G != set():
        g = G.pop()
        if not R.monomial_divides(LM(h), LM(g)):
            G_new.add(g)

    G_new.add(h)

    return G_new,B_new

def select(P):
    """
    The normal selection strategy

    INPUT:

    - ``P`` - a list of critical pairs

    OUTPUT:
        an element of P
    """
    return min(P,key = lambda (fi,fj): LCM(LM(fi),LM(fj)).total_degree())


def inter_reduction(Q):
    """
    If ``Q`` is the set `(f_1, ..., f_n)` this method
    returns `(g_1, ..., g_s)` such that:

    - `<f_1,...,f_n> = <g_1,...,g_s>`
    - `LM(g_i) != LM(g_j)` for all `i != j`
    - `LM(g_i)` does not divide `m` for all monomials `m` of
      `\{g_1,...,g_{i-1}, g_{i+1},...,g_s\}`
    - `LC(g_i) == 1` for all `i`.

    INPUT:

    - ``Q`` - a set of polynomials

    EXAMPLE::

        sage: from sage.rings.polynomial.toy_buchberger import inter_reduction
        sage: inter_reduction(set())
        set([])

        sage: (x,y) = QQ['x,y'].gens()
        sage: reduced = inter_reduction(set([x^2-5*y^2,x^3]))
        sage: reduced == set([x*y^2, x^2-5*y^2])
        True
    """
    if not Q:
        return Q # if Q is empty we cannot get a base ring
    base_ring = iter(Q).next().base_ring()

    Q = set(Q)
    while True:
        Qbar = set(Q)
        for p in Qbar:
            p = Q.pop()
            h = p.reduce(Q)
            if h!=0:
                Q.add(h)
        if Qbar == Q:
            if base_ring.is_field():
                return set([f.lc()**(-1) * f for f in Qbar])
            else: return Qbar
