"""
Educational versions of Groebner basis algorithms.

Following

    [BW93] Thomas Becker and Volker Weispfenning. Groebner Bases - A
           Computational Approach To Commutative Algebra. Springer,
           1993.

the original Buchberger algorithm (c.f. algorithm GROEBNER in [BW93])
and an improved version of Buchberger's algorithm (c.g. algorithm
GROEBNERNEW2 in [BW93]) are presented.

No attempt was made to optimize either algorithm as the emphasis of
these implementations is a clean and easy presentation. To compute a
Groebner basis in \SAGE efficently use the \code{groebner_basis}
method on multivariate polynomial objects. Also, these functions
require a $\Q$ or $\F_p$ as base fields as they rely on a method
\code{MPolynomial.reduce} which is only implemented for those two base
rings yet.

NOTE: the notion of 'term' and 'monomial' in [BW93] is swapped
from the notion of those words in \SAGE (or the other way around,
however you prefer it). In \SAGE a term is a monomial multiplied by a
coefficient, while in [BW93] a monomial is a term multiplied by a
coefficient. Also, what is called LM (the leading monomial) in \SAGE
is called HT (the head term) in [BW93].

EXAMPLES:
    Consider Katsura-6 w.r.t. a $degrevlex$ ordering.

    sage: from sage.rings.polynomial.toy_buchberger import *
    sage: P.<a,b,c,e,f,g,h,i,j,k> = PolynomialRing(GF(32003),10)
    sage: I = sage.rings.ideal.Katsura(P,6)

    sage: g1 = buchberger(I)

    sage: g2 = buchberger_improved(I)

    sage: g3 = I.groebner_basis()

    All algorithms actually compute a Groebner basis:

    sage: Ideal(g1).basis_is_groebner()
    True
    sage: Ideal(g2).basis_is_groebner()
    True
    sage: Ideal(g3).basis_is_groebner()
    True

    The results are correct:

    sage: Ideal(g1) == Ideal(g2) == Ideal(g3)
    True

    If get_verbose() is $>= 1$ a protocol is provided:

    sage: set_verbose(1)
    sage: P.<a,b,c> = PolynomialRing(GF(127),3)
    sage: I = sage.rings.ideal.Katsura(P)

    sage: I
    Ideal (a + 2*b + 2*c - 1, a^2 + 2*b^2 + 2*c^2 - a, 2*a*b + 2*b*c - b) of Polynomial Ring in a, b, c over Finite Field of size 127

    The original Buchberger algorithm performs 15 useless reductions to zero for this example:

    sage: buchberger(I)
    (a^2 + 2*b^2 + 2*c^2 - a, a + 2*b + 2*c - 1) => 6*b^2 + 8*b*c + 6*c^2 - 2*b - 2*c
    G: set([a^2 + 2*b^2 + 2*c^2 - a, 6*b^2 + 8*b*c + 6*c^2 - 2*b - 2*c, a + 2*b + 2*c - 1, 2*a*b + 2*b*c - b])
    <BLANKLINE>
    (a + 2*b + 2*c - 1, 2*a*b + 2*b*c - b) => -44*b*c - 2*c^2 - 21*b + 43*c
    G: set([-44*b*c - 2*c^2 - 21*b + 43*c, a^2 + 2*b^2 + 2*c^2 - a, 6*b^2 + 8*b*c + 6*c^2 - 2*b - 2*c, a + 2*b + 2*c - 1, 2*a*b + 2*b*c - b])
    <BLANKLINE>
    (a + 2*b + 2*c - 1, a^2 + 2*b^2 + 2*c^2 - a) => 0
    G: set([-44*b*c - 2*c^2 - 21*b + 43*c, a^2 + 2*b^2 + 2*c^2 - a, 6*b^2 + 8*b*c + 6*c^2 - 2*b - 2*c, a + 2*b + 2*c - 1, 2*a*b + 2*b*c - b])
    <BLANKLINE>
    (2*a*b + 2*b*c - b, a + 2*b + 2*c - 1) => 0
    G: set([-44*b*c - 2*c^2 - 21*b + 43*c, a^2 + 2*b^2 + 2*c^2 - a, 6*b^2 + 8*b*c + 6*c^2 - 2*b - 2*c, a + 2*b + 2*c - 1, 2*a*b + 2*b*c - b])
    <BLANKLINE>
    (a + 2*b + 2*c - 1, -44*b*c - 2*c^2 - 21*b + 43*c) => 0
    G: set([-44*b*c - 2*c^2 - 21*b + 43*c, a^2 + 2*b^2 + 2*c^2 - a, 6*b^2 + 8*b*c + 6*c^2 - 2*b - 2*c, a + 2*b + 2*c - 1, 2*a*b + 2*b*c - b])
    <BLANKLINE>
    (2*a*b + 2*b*c - b, a^2 + 2*b^2 + 2*c^2 - a) => 44*c^3 - 48*c^2 - 7*b - 3*c
    G: set([-44*b*c - 2*c^2 - 21*b + 43*c, 44*c^3 - 48*c^2 - 7*b - 3*c, a^2 + 2*b^2 + 2*c^2 - a, a + 2*b + 2*c - 1, 6*b^2 + 8*b*c + 6*c^2 - 2*b - 2*c, 2*a*b + 2*b*c - b])
    <BLANKLINE>
    (6*b^2 + 8*b*c + 6*c^2 - 2*b - 2*c, -44*b*c - 2*c^2 - 21*b + 43*c) => 0
    G: set([-44*b*c - 2*c^2 - 21*b + 43*c, 44*c^3 - 48*c^2 - 7*b - 3*c, a^2 + 2*b^2 + 2*c^2 - a, a + 2*b + 2*c - 1, 6*b^2 + 8*b*c + 6*c^2 - 2*b - 2*c, 2*a*b + 2*b*c - b])
    <BLANKLINE>
    (2*a*b + 2*b*c - b, -44*b*c - 2*c^2 - 21*b + 43*c) => 0
    G: set([-44*b*c - 2*c^2 - 21*b + 43*c, 44*c^3 - 48*c^2 - 7*b - 3*c, a^2 + 2*b^2 + 2*c^2 - a, a + 2*b + 2*c - 1, 6*b^2 + 8*b*c + 6*c^2 - 2*b - 2*c, 2*a*b + 2*b*c - b])
    <BLANKLINE>
    (2*a*b + 2*b*c - b, 6*b^2 + 8*b*c + 6*c^2 - 2*b - 2*c) => 0
    G: set([-44*b*c - 2*c^2 - 21*b + 43*c, 44*c^3 - 48*c^2 - 7*b - 3*c, a^2 + 2*b^2 + 2*c^2 - a, a + 2*b + 2*c - 1, 6*b^2 + 8*b*c + 6*c^2 - 2*b - 2*c, 2*a*b + 2*b*c - b])
    <BLANKLINE>
    (a^2 + 2*b^2 + 2*c^2 - a, 2*a*b + 2*b*c - b) => 0
    G: set([-44*b*c - 2*c^2 - 21*b + 43*c, 44*c^3 - 48*c^2 - 7*b - 3*c, a^2 + 2*b^2 + 2*c^2 - a, a + 2*b + 2*c - 1, 6*b^2 + 8*b*c + 6*c^2 - 2*b - 2*c, 2*a*b + 2*b*c - b])
    <BLANKLINE>
    (a + 2*b + 2*c - 1, 6*b^2 + 8*b*c + 6*c^2 - 2*b - 2*c) => 0
    G: set([-44*b*c - 2*c^2 - 21*b + 43*c, 44*c^3 - 48*c^2 - 7*b - 3*c, a^2 + 2*b^2 + 2*c^2 - a, a + 2*b + 2*c - 1, 6*b^2 + 8*b*c + 6*c^2 - 2*b - 2*c, 2*a*b + 2*b*c - b])
    <BLANKLINE>
    (a^2 + 2*b^2 + 2*c^2 - a, 6*b^2 + 8*b*c + 6*c^2 - 2*b - 2*c) => 0
    G: set([-44*b*c - 2*c^2 - 21*b + 43*c, 44*c^3 - 48*c^2 - 7*b - 3*c, a^2 + 2*b^2 + 2*c^2 - a, a + 2*b + 2*c - 1, 6*b^2 + 8*b*c + 6*c^2 - 2*b - 2*c, 2*a*b + 2*b*c - b])
    <BLANKLINE>
    (-44*b*c - 2*c^2 - 21*b + 43*c, 44*c^3 - 48*c^2 - 7*b - 3*c) => 0
    G: set([-44*b*c - 2*c^2 - 21*b + 43*c, 44*c^3 - 48*c^2 - 7*b - 3*c, a^2 + 2*b^2 + 2*c^2 - a, a + 2*b + 2*c - 1, 6*b^2 + 8*b*c + 6*c^2 - 2*b - 2*c, 2*a*b + 2*b*c - b])
    <BLANKLINE>
    (a + 2*b + 2*c - 1, 44*c^3 - 48*c^2 - 7*b - 3*c) => 0
    G: set([-44*b*c - 2*c^2 - 21*b + 43*c, 44*c^3 - 48*c^2 - 7*b - 3*c, a^2 + 2*b^2 + 2*c^2 - a, a + 2*b + 2*c - 1, 6*b^2 + 8*b*c + 6*c^2 - 2*b - 2*c, 2*a*b + 2*b*c - b])
    <BLANKLINE>
    (a^2 + 2*b^2 + 2*c^2 - a, -44*b*c - 2*c^2 - 21*b + 43*c) => 0
    G: set([-44*b*c - 2*c^2 - 21*b + 43*c, 44*c^3 - 48*c^2 - 7*b - 3*c, a^2 + 2*b^2 + 2*c^2 - a, a + 2*b + 2*c - 1, 6*b^2 + 8*b*c + 6*c^2 - 2*b - 2*c, 2*a*b + 2*b*c - b])
    <BLANKLINE>
    (a^2 + 2*b^2 + 2*c^2 - a, 44*c^3 - 48*c^2 - 7*b - 3*c) => 0
    G: set([-44*b*c - 2*c^2 - 21*b + 43*c, 44*c^3 - 48*c^2 - 7*b - 3*c, a^2 + 2*b^2 + 2*c^2 - a, a + 2*b + 2*c - 1, 6*b^2 + 8*b*c + 6*c^2 - 2*b - 2*c, 2*a*b + 2*b*c - b])
    <BLANKLINE>
    (6*b^2 + 8*b*c + 6*c^2 - 2*b - 2*c, 44*c^3 - 48*c^2 - 7*b - 3*c) => 0
    G: set([-44*b*c - 2*c^2 - 21*b + 43*c, 44*c^3 - 48*c^2 - 7*b - 3*c, a^2 + 2*b^2 + 2*c^2 - a, a + 2*b + 2*c - 1, 6*b^2 + 8*b*c + 6*c^2 - 2*b - 2*c, 2*a*b + 2*b*c - b])
    <BLANKLINE>
    (2*a*b + 2*b*c - b, 44*c^3 - 48*c^2 - 7*b - 3*c) => 0
    G: set([-44*b*c - 2*c^2 - 21*b + 43*c, 44*c^3 - 48*c^2 - 7*b - 3*c, a^2 + 2*b^2 + 2*c^2 - a, a + 2*b + 2*c - 1, 6*b^2 + 8*b*c + 6*c^2 - 2*b - 2*c, 2*a*b + 2*b*c - b])
    <BLANKLINE>
    15 reductions to zero.
    [-44*b*c - 2*c^2 - 21*b + 43*c, 44*c^3 - 48*c^2 - 7*b - 3*c, a^2 + 2*b^2 + 2*c^2 - a, a + 2*b + 2*c - 1, 6*b^2 + 8*b*c + 6*c^2 - 2*b - 2*c, 2*a*b + 2*b*c - b]

    The 'improved' Buchberger algorithm in constrast only performs 3 reductions to zero:

    sage: buchberger_improved(I)
    (2*a*b + 2*b*c - b, a + 2*b + 2*c - 1) => -2*b^2 - b*c - 63*b
    G: set([a^2 + 2*b^2 + 2*c^2 - a, -2*b^2 - b*c - 63*b, a + 2*b + 2*c - 1, 2*a*b + 2*b*c - b])
    <BLANKLINE>
    (a^2 + 2*b^2 + 2*c^2 - a, a + 2*b + 2*c - 1) => 5*b*c + 6*c^2 + 63*b - 2*c
    G: set([5*b*c + 6*c^2 + 63*b - 2*c, a^2 + 2*b^2 + 2*c^2 - a, -2*b^2 - b*c - 63*b, a + 2*b + 2*c - 1, 2*a*b + 2*b*c - b])
    <BLANKLINE>
    (-2*b^2 - b*c - 63*b, 2*a*b + 2*b*c - b) => -22*c^3 + 24*c^2 - 60*b - 62*c
    G: set([-22*c^3 + 24*c^2 - 60*b - 62*c, a^2 + 2*b^2 + 2*c^2 - a, 5*b*c + 6*c^2 + 63*b - 2*c, a + 2*b + 2*c - 1, -2*b^2 - b*c - 63*b, 2*a*b + 2*b*c - b])
    <BLANKLINE>
    (5*b*c + 6*c^2 + 63*b - 2*c, 2*a*b + 2*b*c - b) => 0
    G: set([-22*c^3 + 24*c^2 - 60*b - 62*c, a^2 + 2*b^2 + 2*c^2 - a, 5*b*c + 6*c^2 + 63*b - 2*c, a + 2*b + 2*c - 1, -2*b^2 - b*c - 63*b, 2*a*b + 2*b*c - b])
    <BLANKLINE>
    (5*b*c + 6*c^2 + 63*b - 2*c, -2*b^2 - b*c - 63*b) => 0
    G: set([-22*c^3 + 24*c^2 - 60*b - 62*c, a^2 + 2*b^2 + 2*c^2 - a, 5*b*c + 6*c^2 + 63*b - 2*c, a + 2*b + 2*c - 1, -2*b^2 - b*c - 63*b, 2*a*b + 2*b*c - b])
    <BLANKLINE>
    (-22*c^3 + 24*c^2 - 60*b - 62*c, 5*b*c + 6*c^2 + 63*b - 2*c) => 0
    G: set([-22*c^3 + 24*c^2 - 60*b - 62*c, a^2 + 2*b^2 + 2*c^2 - a, 5*b*c + 6*c^2 + 63*b - 2*c, a + 2*b + 2*c - 1, -2*b^2 - b*c - 63*b, 2*a*b + 2*b*c - b])
    <BLANKLINE>
    3 reductions to zero.
    [-22*c^3 + 24*c^2 - 60*b - 62*c, a^2 + 2*b^2 + 2*c^2 - a, 5*b*c + 6*c^2 + 63*b - 2*c, a + 2*b + 2*c - 1, -2*b^2 - b*c - 63*b, 2*a*b + 2*b*c - b]

AUTHOR:
    -- Martin Albrecht (2007-05-24): initial version
"""

from sage.misc.misc import get_verbose
from sage.misc.misc import exists
from sage.rings.arith import LCM
from sage.structure.sequence import Sequence

LM = lambda f: f.lm()
LT = lambda f: f.lt()
spol = lambda f,g: LCM(LM(f),LM(g)) // LT(f) * f - LCM(LM(f),LM(g)) // LT(g) * g

def buchberger(F):
    """
    The original version of Buchberger's algorithm as presented in
    [BW93], page 214.

    INPUT:
        F -- an ideal in a multivariate polynomial ring

    OUTPUT:
        a Groebner basis for F

    NOTE: The verbosity of this function may be controlled with a
    \code{set_verbose()} call. Any value >=1 will result in this
    function printing intermediate Groebner bases.

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
    [BW93], page 232.

    INPUT:
        F -- an ideal in a multivariate polynomial ring

    OUTPUT:
        a Groebner basis for F

    NOTE: The verbosity of this function may be controlled with a
    \code{set_verbose()} call. Any value >=1 will result in this
    function printing intermediate Groebner bases.
    """
    F = set(F.gens())

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

    return Sequence(G)


def update(G,B,h):
    """
    Update $G$ using the list of critical pairs $B$ and the polynomial
    $h$ as presented in [BW93], page 230. For this, Buchberger's first
    and second criterion are tested.

    INPUT:
        G -- an intermediate Groebner basis
        B -- a list of critical pairs
        h -- a polynomial

    OUTPUT:
        a tuple of an intermediate Groebner basis and a list of
        critical pairs

    """
    R = h.parent()

    C = set([(h,g) for g in G])
    D = set()

    while C != set():
        (h,g1) = C.pop()

        lcm_divides = lambda rhs: R.monomial_is_divisible_by( LCM(LM(h),LM(g1)), LCM(LM(h),LM(rhs[1])) )

        if R.monomial_pairwise_prime(LM(h),LM(g)) or \
           (\
               not exists(C, lcm_divides )[0] \
               and \
               not exists(D, lcm_divides )[0] \
            ):
            D.add( (h,g1) )

    E = set()

    while D != set():
        (h,g) = D.pop()
        if not R.monomial_pairwise_prime(LM(h),LM(g)):
            E.add( (h,g) )

    B_new = set()

    while B != set():
        g1,g2 = B.pop()
        if not R.monomial_is_divisible_by( LCM(LM(g1),LM(g2)), LM(h) ) or \
               R.monomial_lcm(LM(g1),LM( h)) == LCM(LM(g1),LM(g2)) or \
               R.monomial_lcm(LM( h),LM(g2)) == LCM(LM(g1),LM(g2)) :
            B_new.add( (g1,g2) )

    B_new = B_new.union( E )

    G_new = set()

    while G != set():
        g = G.pop()
        if not R.monomial_is_divisible_by(LM(g),LM(h)):
            G_new.add(g)

    G_new.add(h)

    return G_new,B_new

def select(P):
    """
    The normal selection strategy

    INPUT:
        P -- a list of critical pairs

    OUTPUT:
        an element of P
    """
    return min(P,key = lambda (fi,fj): LCM(LM(fi),LM(fj)).total_degree())

