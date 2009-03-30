# -*- coding: utf-8 -*-
r"""
Ideals in multivariate polynomial rings.

Sage has a powerful system to compute with multivariate polynomial
rings. Most algorithms dealing with these ideals are centered the
computation of *Groebner bases*. Sage makes use of Singular to
implement this functionality. Singular is widely regarded as the best
open-source system for Groebner basis calculation in multivariate
polynomial rings over fields.

AUTHORS:

- William Stein

- Kiran S. Kedlaya (2006-02-12): added Macaulay2 analogues of some
  Singular features

- Martin Albrecht (2008,2007): refactoring, many Singular related
  functions

EXAMPLES:

We compute a Groebner basis for some given ideal. The type returned by
the ``groebner_basis`` method is ``Sequence``, i.e. it is not an
``MPolynomialIdeal``.

::

    sage: x,y,z = QQ['x,y,z'].gens()
    sage: I = ideal(x^5 + y^4 + z^3 - 1,  x^3 + y^3 + z^2 - 1)
    sage: B = I.groebner_basis()
    sage: type(B)
    <class 'sage.structure.sequence.Sequence'>

Groebner bases can be used to solve the ideal membership problem.

::

    sage: f,g,h = B
    sage: (2*x*f + g).reduce(B)
    0

    sage: (2*x*f + g) in I
    True

    sage: (2*x*f + 2*z*h + y^3).reduce(B)
    y^3

    sage: (2*x*f + 2*z*h + y^3) in I
    False

We compute a Groebner basis for cyclic 6, which is a standard
benchmark and test ideal.

::

    sage: R.<x,y,z,t,u,v> = QQ['x,y,z,t,u,v']
    sage: I = sage.rings.ideal.Cyclic(R,6)
    sage: B = I.groebner_basis()
    sage: len(B)
    45

We compute in a quotient of a polynomial ring over `\mathbb{Z}/17\mathbb{Z}`::

    sage: R.<x,y> = ZZ[]
    sage: S.<a,b> = R.quotient((x^2 + y^2, 17))     # optional -- macaulay2
    ...
    verbose 0 ... Warning: falling back to very slow toy implementation.
    sage: S                                         # optional -- macaulay2
    Quotient of Multivariate Polynomial Ring in x, y over Integer Ring
    by the ideal (x^2 + y^2, 17)

    sage: a^2 + b^2 == 0                            # optional -- macaulay2
    True
    sage: a^3 - b^2                                 # optional -- macaulay2
    a*b^2 - b^2
    sage: (a+b)^17                                  # optional -- macaulay2
    a*b^16 + b^17
    sage: S(17) == 0                                # optional -- macaulay2
    True

Working with a polynomial ring over `\mathbb{Z}`::

    sage: R.<x,y,z,w> = ZZ['x,y,z,w']
    sage: i = ideal(x^2 + y^2 - z^2 - w^2, x-y)
    sage: j = i^2
    sage: j.groebner_basis('macaulay2')          # optional - macaulay2
    [4*y^4 - 4*y^2*z^2 + z^4 - 4*y^2*w^2 + 2*z^2*w^2 + w^4,
     2*x*y^2 - 2*y^3 - x*z^2 + y*z^2 - x*w^2 + y*w^2,
     x^2 - 2*x*y + y^2]

    sage: y^2 - 2*x*y + x^2 in j                # optional - macaulay2
    True
    sage: 0 in j                                # optional - macaulay2
    True

We do a Groebner basis computation over a number field::

    sage: K.<zeta> = CyclotomicField(3)
    sage: R.<x,y,z> = K[]; R
    Multivariate Polynomial Ring in x, y, z over Cyclotomic Field of order 3 and degree 2

    sage: i = ideal(x - zeta*y + 1, x^3 - zeta*y^3); i
    Ideal (x + (-zeta)*y + 1, x^3 + (-zeta)*y^3) of Multivariate
    Polynomial Ring in x, y, z over Cyclotomic Field of order 3 and degree 2

    sage: i.groebner_basis()
    [y^3 + (2*zeta + 1)*y^2 + (zeta - 1)*y + (-1/3*zeta - 2/3), x + (-zeta)*y + 1]

    sage: S = R.quotient(i); S
    Quotient of Multivariate Polynomial Ring in x, y, z over
    Cyclotomic Field of order 3 and degree 2 by the ideal (x +
    (-zeta)*y + 1, x^3 + (-zeta)*y^3)

    sage: S.0  - zeta*S.1
    -1
    sage: S.0^3 - zeta*S.1^3
    0

Two examples from the Mathematica documentation (done in Sage):

    We compute a Groebner basis::

        sage: R.<x,y> = PolynomialRing(QQ, order='lex')
        sage: ideal(x^2 - 2*y^2, x*y - 3).groebner_basis()
         [x - 2/3*y^3, y^4 - 9/2]

    We show that three polynomials have no common root::

        sage: R.<x,y> = QQ[]
        sage: ideal(x+y, x^2 - 1, y^2 - 2*x).groebner_basis()
        [1]

The next example shows how we can use Groebner bases over
`\mathbb{Z}` to find the primes modulo which a system of
equations has a solution, when the system has no solutions over the
rationals.

    We first form a certain ideal `I` in
    `\mathbb{Z}[x, y, z]`, and note that the Groebner basis of
    `I` over `\mathbb{Q}` contains 1, so there are no
    solutions over `\mathbb{Q}` or an algebraic closure of it
    (this is not surprising as there are 4 equations in 3 unknowns).

    ::

        sage: P.<x,y,z> = PolynomialRing(ZZ,order='lex')
        sage: I = ideal(-y^2 - 3*y + z^2 + 3, -2*y*z + z^2 + 2*z + 1, \
                        x*z + y*z + z^2, -3*x*y + 2*y*z + 6*z^2)
        sage: I.change_ring(P.change_ring(QQ)).groebner_basis()
        [1]

    However, when we compute the Groebner basis of I (defined over
    $\ZZ$), we note that there is a certain integer in the ideal which
    is not 1.

    ::

        sage: I.groebner_basis()          # optional -- macaulay2
        ...
        verbose 0 ... Warning: falling back to very slow toy implementation.
        [x + y + z, y^2 + y + 23234, y*z + y + 26532, 2*y + 158864, z^2 + 17223, 2*z + 41856, 164878]

    Now for each prime `p` dividing this integer 164878, the
    Groebner basis of I modulo `p` will be non-trivial and will thus
    give a solution of the original system modulo `p`.

    ::


        sage: factor(164878)
        2 * 7 * 11777

        sage: I.change_ring(P.change_ring( GF(2) )).groebner_basis()
        [x + y + z, y^2 + y, y*z + y, z^2 + 1]
        sage: I.change_ring(P.change_ring( GF(7) )).groebner_basis()
        [x - 1, y + 3, z - 2]
        sage: I.change_ring(P.change_ring( GF(11777 ))).groebner_basis()
        [x + 5633, y - 3007, z - 2626]

    The Groebner basis modulo any product of the prime factors is also non-trivial.

    ::

        sage: I.change_ring(P.change_ring( IntegerModRing(2*7) )).groebner_basis()   # optional -- macaulay2
        verbose 0 (...: multi_polynomial_ideal.py, groebner_basis) Warning: falling back to very slow toy implementation.
        [x + y + z, y^2 + y + 8, y*z + y + 2, 2*y + 6, z^2 + 3, 2*z + 10]

    Modulo any other prime the Groebner basis is trivial so there are
    no other solutions. For example::

        sage: I.change_ring( P.change_ring( GF(3) ) ).groebner_basis()
        [1]


TESTS::

    sage: x,y,z = QQ['x,y,z'].gens()
    sage: I = ideal(x^5 + y^4 + z^3 - 1,  x^3 + y^3 + z^2 - 1)
    sage: I == loads(dumps(I))
    True

"""

#*****************************************************************************
#
#                               Sage
#
#       Copyright (C) 2005 William Stein <wstein@gmail.com>
#       Copyright (C) 2008 Martin Albrecht <malb@informatik.uni-bremen.de>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#
#    This code is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
#    General Public License for more details.
#
#  The full text of the GPL is available at:
#
#                  http://www.gnu.org/licenses/
#*****************************************************************************
from __future__ import with_statement

from sage.interfaces.all import (singular as singular_default,
                                 macaulay2 as macaulay2_default,
                                 magma as magma_default)

from sage.rings.ideal import Ideal_generic
from sage.rings.integer import Integer
from sage.structure.sequence import Sequence

from sage.misc.cachefunc import cached_method
from sage.misc.misc import prod, verbose
from sage.misc.sage_eval import sage_eval

from sage.rings.integer_ring import ZZ
import sage.rings.polynomial.toy_buchberger as toy_buchberger
import sage.rings.polynomial.toy_variety as toy_variety
import sage.rings.polynomial.toy_d_basis as toy_d_basis

class RedSBContext:
    """
    Within this context all Singular Groebner basis calculations are
    reduced automatically.

    AUTHORS:

    - Martin Albrecht
    """
    def __init__(self, singular=singular_default):
        r"""
        Within this context all Singular Groebner basis calculations are
        reduced automatically.

        INPUT:


        -  ``singular`` - Singular instance (default: default
           instance)


        EXAMPLE::

            sage: from sage.rings.polynomial.multi_polynomial_ideal import RedSBContext
            sage: P.<a,b,c> = PolynomialRing(QQ,3, order='lex')
            sage: I = sage.rings.ideal.Katsura(P,3)
            sage: singular.option('noredTail')
            sage: singular.option('noredThrough')
            sage: Is = I._singular_()
            sage: gb = Is.groebner()
            sage: gb
            84*c^4-40*c^3+c^2+c,
            7*b+210*c^3-79*c^2+3*c,
            a+2*b+2*c-1

        ::

            sage: from __future__ import with_statement
            sage: with RedSBContext(): rgb = Is.groebner()
            sage: rgb
            84*c^4-40*c^3+c^2+c,
            7*b+210*c^3-79*c^2+3*c,
            7*a-420*c^3+158*c^2+8*c-7

        Note that both bases are Groebner bases because they have pairwise
        prime leading monomials but that the monic version of the last
        element in ``rgb`` is smaller than the last element of
        ``gb`` with respect to the lexicographical term
        ordering.

        ::

            sage: (7*a-420*c^3+158*c^2+8*c-7)/7 < (a+2*b+2*c-1)
            True

        .. note::

           This context is used automatically internally whenever a
           Groebner basis is computed so the user does not need to use
           it manually.
        """
        self.singular = singular

    def __enter__(self):
        """
        EXAMPLE::

            sage: from sage.rings.polynomial.multi_polynomial_ideal import RedSBContext
            sage: P.<a,b,c> = PolynomialRing(QQ,3, order='lex')
            sage: I = sage.rings.ideal.Katsura(P,3)
            sage: singular.option('noredTail')
            sage: singular.option('noredThrough')
            sage: Is = I._singular_()
            sage: with RedSBContext(): rgb = Is.groebner()
            sage: rgb
            84*c^4-40*c^3+c^2+c,
            7*b+210*c^3-79*c^2+3*c,
            7*a-420*c^3+158*c^2+8*c-7
        """
        self.o = self.singular.option("get")
        self.singular.option("redSB")

    def __exit__(self, typ, value, tb):
        """
        EXAMPLE::

            sage: from sage.rings.polynomial.multi_polynomial_ideal import RedSBContext
            sage: P.<a,b,c> = PolynomialRing(QQ,3, order='lex')
            sage: I = sage.rings.ideal.Katsura(P,3)
            sage: singular.option('noredTail')
            sage: singular.option('noredThrough')
            sage: Is = I._singular_()
            sage: with RedSBContext(): rgb = Is.groebner()
            sage: rgb
            84*c^4-40*c^3+c^2+c,
            7*b+210*c^3-79*c^2+3*c,
            7*a-420*c^3+158*c^2+8*c-7
        """
        self.singular.option("set",self.o)

def redSB(func):
    """
    Decorator to force a reduced Singular groebner basis.

    .. note::

       This decorator is used automatically internally so the user
       does not need to use it manually.
    """
    def wrapper(*args, **kwds):
        """
        execute fucntion in ``RedSBContext``.
        """
        with RedSBContext():
            return func(*args, **kwds)
    wrapper.__doc__=func.__doc__
    return wrapper

def is_MPolynomialIdeal(x):
    r"""
    Return ``True`` if the provided argument
    ``x`` is an ideal in the multivariate polynomial ring.

    INPUT:


    -  ``x`` - an arbitrary object


    EXAMPLES::

        sage: from sage.rings.polynomial.all import is_MPolynomialIdeal
        sage: P.<x,y,z> = PolynomialRing(QQ)
        sage: I = [x + 2*y + 2*z - 1, x^2 + 2*y^2 + 2*z^2 - x, 2*x*y + 2*y*z - y]

    Sage distinguishes between a list of generators for an ideal and
    the ideal itself. This distinction is inconsistent with Singular but
    matches Magma's behavior.

    ::

        sage: is_MPolynomialIdeal(I)
        False

    ::

        sage: I = Ideal(I)
        sage: is_MPolynomialIdeal(I)
        True
    """
    return isinstance(x, MPolynomialIdeal)

class MPolynomialIdeal_magma_repr:
    def _magma_init_(self, magma):
        r"""
        Returns a Magma ideal matching this ideal if the base ring coerces
        to Magma and Magma is available.

        INPUT:


        -  ``magma`` - Magma instance


        EXAMPLES::

            sage: R.<a,b,c,d,e,f,g,h,i,j> = PolynomialRing(GF(127),10)
            sage: I = sage.rings.ideal.Cyclic(R,4)
            sage: magma(I)                                          # optional - magma
            Ideal of Polynomial ring of rank 10 over GF(127)
            Graded Reverse Lexicographical Order
            Variables: a, b, c, d, e, f, g, h, i, j
            Basis:
            [
            a + b + c + d,
            a*b + b*c + a*d + c*d,
            a*b*c + a*b*d + a*c*d + b*c*d,
            a*b*c*d + 126
            ]
        """
        P = magma(self.ring())
        G = magma(self.gens())
        return 'ideal<%s|%s>'%(P.name(), G._ref())

    def _groebner_basis_magma(self, magma=magma_default):
        r"""
        Computes a Groebner Basis for self using Magma if available.

        INPUT:


        -  ``magma`` - Magma instance or None (default
           instance) (default: None)


        EXAMPLES::

            sage: R.<a,b,c,d,e,f,g,h,i,j> = PolynomialRing(GF(127),10)
            sage: I = sage.rings.ideal.Cyclic(R,6)
            sage: gb = I.groebner_basis('magma:GroebnerBasis') # indirect doctest; optional - magma
            sage: len(gb)                                      # optional - magma
            45
        """
        R   = self.ring()
        mgb = magma(self).GroebnerBasis()
        # TODO: rewrite this to be much more sophisticated in multi-level nested cases.
        mgb = [str(mgb[i+1]) for i in range(len(mgb))]
        if R.base_ring().degree() > 1:
            a = str(R.base_ring().gen())
            mgb = [e.replace("$.1",a) for e in mgb]
        B = Sequence([R(e) for e in mgb], R, check=False, immutable=True)
        return B

class MPolynomialIdeal_singular_repr:
    r"""
    An ideal in a multivariate polynomial ring, which has an underlying
    Singular ring associated to it.
    """
    def _singular_(self, singular=singular_default):
        r"""
        Return Singular ideal corresponding to this ideal.

        EXAMPLES::

            sage: R.<x,y> = PolynomialRing(QQ, 2)
            sage: I = R.ideal([x^3 + y, y])
            sage: S = I._singular_()
            sage: S
            x^3+y,
            y
        """
        try:
            self.ring()._singular_(singular).set_ring()
            I = self.__singular
            if not (I.parent() is singular):
                raise ValueError
            I._check_valid()
            return I
        except (AttributeError, ValueError):
            self.ring()._singular_(singular).set_ring()
            try:
                # this may fail for quotient ring elements, but is
                # faster
                gens = [str(x) for x in self.gens()]
                if len(gens) == 0:
                    gens = ['0']
                self.__singular = singular.ideal(gens)
            except TypeError:
                gens = [str(x.lift()) for x in self.gens()]
                if len(gens) == 0:
                    gens = ['0']
                self.__singular = singular.ideal(gens)
        return self.__singular

    def plot(self, singular=singular_default):
        """
        If you somehow manage to install surf, perhaps you can use this
        function to implicitly plot the real zero locus of this ideal (if
        principal).

        INPUT:


        -  ``self`` - must be a principal ideal in 2 or 3 vars
           over QQ.


        EXAMPLES:

        Implicit plotting in 2-d::

            sage: R.<x,y> = PolynomialRing(QQ,2)
            sage: I = R.ideal([y^3 - x^2])
            sage: I.plot()        # cusp         optional - surf
            sage: I = R.ideal([y^2 - x^2 - 1])
            sage: I.plot()        # hyperbola    optional - surf
            sage: I = R.ideal([y^2 + x^2*(1/4) - 1])
            sage: I.plot()        # ellipse      optional - surf
            sage: I = R.ideal([y^2-(x^2-1)*(x-2)])
            sage: I.plot()        # elliptic curve  optional - surf

        Implicit plotting in 3-d::

            sage: R.<x,y,z> = PolynomialRing(QQ,3)
            sage: I = R.ideal([y^2 + x^2*(1/4) - z])
            sage: I.plot()          # a cone         optional - surf
            sage: I = R.ideal([y^2 + z^2*(1/4) - x])
            sage: I.plot()          # same code, from a different angle  optional - surf
            sage: I = R.ideal([x^2*y^2+x^2*z^2+y^2*z^2-16*x*y*z])
            sage: I.plot()          # Steiner surface   optional - surf

        AUTHORS:

        - David Joyner (2006-02-12)
        """
        if self.ring().characteristic() != 0:
            raise TypeError, "base ring must have characteristic 0"
        if not self.is_principal():
            raise TypeError, "self must be principal"
        singular.lib('surf')
        I = singular(self)
        I.plot()

    @redSB
    def complete_primary_decomposition(self, algorithm="sy"):
        r"""
        Return a list of primary ideals and their associated primes such
        that the intersection of the primary ideal `Q_i` is
        `I` = ``self``.

        An ideal `Q` is called primary if it is a proper ideal of
        the ring `R` and if whenever `ab \in Q` and
        `a \not\in Q` then `b^n \in Q` for some
        `n \in \mathbb{Z}`.

        If `Q` is a primary ideal of the ring `R`, then the
        radical ideal `P` of `Q`, i.e.
        `P = \{a \in R, a^n \in Q\}` for some
        `n \in \mathbb{Z}`, is called the
        *associated prime* of `Q`.

        If `I` is a proper ideal of the ring `R` then there
        exists a decomposition in primary ideals `Q_i` such that


        -  their intersection is `I`

        -  none of the `Q_i` contains the intersection of the
           rest, and

        -  the associated prime ideals of `Q_i` are pairwise
           different.


        This method returns these `Q_i` and their associated
        primes.

        INPUT:


        -  ``algorithm`` - string:

          -  ``'sy'`` - (default) use the shimoyama-yokoyama
             algorithm

          -  ``'gtz'`` - use the gianni-trager-zacharias
             algorithm


        OUTPUT:


        -  ``list`` - a list of primary ideals and their
           associated primes [(primary ideal, associated prime), ...]


        EXAMPLES::

            sage: R.<x,y,z> = PolynomialRing(QQ, 3, order='lex')
            sage: p = z^2 + 1; q = z^3 + 2
            sage: I = (p*q^2, y-z^2)*R
            sage: pd = I.complete_primary_decomposition(); pd
            [(Ideal (z^6 + 4*z^3 + 4, y - z^2) of Multivariate Polynomial Ring in x, y, z over Rational Field,
              Ideal (z^3 + 2, y - z^2) of Multivariate Polynomial Ring in x, y, z over Rational Field),
             (Ideal (z^2 + 1, y + 1) of Multivariate Polynomial Ring in x, y, z over Rational Field,
              Ideal (z^2 + 1, y + 1) of Multivariate Polynomial Ring in x, y, z over Rational Field)]

        ::

            sage: I.complete_primary_decomposition(algorithm = 'gtz')
            [(Ideal (z^6 + 4*z^3 + 4, y - z^2) of Multivariate Polynomial Ring in x, y, z over Rational Field,
              Ideal (z^3 + 2, y - z^2) of Multivariate Polynomial Ring in x, y, z over Rational Field),
             (Ideal (z^2 + 1, y - z^2) of Multivariate Polynomial Ring in x, y, z over Rational Field,
              Ideal (z^2 + 1, y - z^2) of Multivariate Polynomial Ring in x, y, z over Rational Field)]

        ::

            sage: reduce(lambda Qi,Qj: Qi.intersection(Qj), [Qi for (Qi,radQi) in pd]) == I
            True

        ::

            sage: [Qi.radical() == radQi for (Qi,radQi) in pd]
            [True, True]

        ALGORITHM: Uses Singular.

        REFERENCES:

        - Thomas Becker and Volker Weispfenning. Groebner Bases - A
          Computational Approach To Commutative Algebra. Springer, New
          York 1993.
        """
        try:
            return self.__complete_primary_decomposition[algorithm]
        except AttributeError:
            self.__complete_primary_decomposition = {}
        except KeyError:
            pass
        I = self._singular_()
        I.parent().lib('primdec.lib')
        if algorithm == 'sy':
            P = I.primdecSY()
        elif algorithm == 'gtz':
            P = I.primdecGTZ()

        R = self.ring()
        V = [(R.ideal(X[1]), R.ideal(X[2])) for X in P]
        V = Sequence(V)
        self.__complete_primary_decomposition[algorithm] = V
        return self.__complete_primary_decomposition[algorithm]

    def primary_decomposition(self, algorithm='sy'):
        r"""
        Return a list of primary ideals such that their intersection is
        `I` = ``self``.

        An ideal `Q` is called primary if it is a proper ideal of
        the ring `R` and if whenever `ab \in Q` and
        `a \not\in Q` then `b^n \in Q` for some
        `n \in \mathbb{Z}`.

        If `I` is a proper ideal of the ring `R` then there
        exists a decomposition in primary ideals `Q_i` such that


        -  their intersection is `I`

        -  none of the `Q_i` contains the intersection of the
           rest, and

        -  the associated prime ideals of `Q_i` are pairwise
           different.


        This method returns these `Q_i`.

        INPUT:


        -  ``algorithm`` - string:

          -  ``'sy'`` - (default) use the shimoyama-yokoyama
             algorithm

          -  ``'gtz'`` - use the gianni-trager-zacharias
             algorithm


        EXAMPLES::

            sage: R.<x,y,z> = PolynomialRing(QQ, 3, order='lex')
            sage: p = z^2 + 1; q = z^3 + 2
            sage: I = (p*q^2, y-z^2)*R
            sage: pd = I.primary_decomposition(); pd
            [Ideal (z^6 + 4*z^3 + 4, y - z^2) of Multivariate Polynomial Ring in x, y, z over Rational Field,
             Ideal (z^2 + 1, y + 1) of Multivariate Polynomial Ring in x, y, z over Rational Field]

        ::

            sage: reduce(lambda Qi,Qj: Qi.intersection(Qj), pd) == I
            True

        ALGORITHM: Uses Singular.

        REFERENCES:

        - Thomas Becker and Volker Weispfenning. Groebner Bases - A
          Computational Approach To Commutative Algebra. Springer, New
          York 1993.
        """
        return [I for I, _ in self.complete_primary_decomposition(algorithm)]

    @redSB
    def associated_primes(self, algorithm='sy'):
        r"""
        Return a list of primary ideals (and their associated primes) such
        that their intersection is `I` = ``self``.

        An ideal `Q` is called primary if it is a proper ideal of
        the ring `R` and if whenever `ab \in Q` and
        `a \not\in Q` then `b^n \in Q` for some
        `n \in \mathbb{Z}`.

        If `Q` is a primary ideal of the ring `R`, then the
        radical ideal `P` of `Q`, i.e.
        `P = \{a \in R, a^n \in Q\}` for some
        `n \in \mathbb{Z}`, is called the
        *associated prime* of `Q`.

        If `I` is a proper ideal of the ring `R` then there
        exists a decomposition in primary ideals `Q_i` such that


        -  their intersection is `I`

        -  none of the `Q_i` contains the intersection of the
           rest, and

        -  the associated prime ideals of `Q_i` are pairwise
           different.


        This method returns the associated primes of the `Q_i`.

        INPUT:


        -  ``algorithm`` - string:

          -  ``'sy'`` - (default) use the shimoyama-yokoyama
             algorithm

          -  ``'gtz'`` - use the gianni-trager-zacharias
             algorithm


        OUTPUT:


        -  ``list`` - a list of primary ideals and their
           associated primes [(primary ideal, associated prime), ...]


        EXAMPLES::

            sage: R.<x,y,z> = PolynomialRing(QQ, 3, order='lex')
            sage: p = z^2 + 1; q = z^3 + 2
            sage: I = (p*q^2, y-z^2)*R
            sage: pd = I.associated_primes(); pd
            [Ideal (z^3 + 2, y - z^2) of Multivariate Polynomial Ring in x, y, z over Rational Field,
             Ideal (z^2 + 1, y + 1) of Multivariate Polynomial Ring in x, y, z over Rational Field]

        ALGORITHM: Uses Singular.

        REFERENCES:

        - Thomas Becker and Volker Weispfenning. Groebner Bases - A
          Computational Approach To Commutative Algebra. Springer, New
          York 1993.
        """
        return [P for _,P in self.complete_primary_decomposition(algorithm)]

    def triangular_decomposition(self, algorithm=None, singular=singular_default):
        r"""
            Decompose zero-dimensional ideal ``self`` into
            triangular sets.

            This requires that the given basis is reduced w.r.t. to the
            lexicographical monomial ordering. If the basis of self does not
            have this property, the required Groebner basis is computed
            implicitly.

            INPUT:


            -  ``algorithm`` - string or None (default: None)

               ALGORITHMS:

                 -  ``singular:triangL`` - decomposition of self into
                    triangular systems (Lazard).

                 -  ``singular:triangLfak`` - decomp. of self into tri.
                    systems plus factorization.

                 -  ``singular:triangM`` - decomposition of self into
                    triangular systems (Moeller).


            OUTPUT: a list `T` of lists `t` such that the
            variety of ``self`` is the union of the varieties of
            `t` in `L` and each `t` is in triangular
            form.

            EXAMPLE::

                sage: P.<e,d,c,b,a> = PolynomialRing(QQ,5,order='lex')
                sage: I = sage.rings.ideal.Cyclic(P)
                sage: GB = Ideal(I.groebner_basis('singular:stdfglm'))
                sage: GB.triangular_decomposition('singular:triangLfak')
                [Ideal (a - 1, b - 1, c - 1, d^2 + 3*d + 1, e + d + 3) of Multivariate Polynomial Ring in e, d, c, b, a over Rational Field,
                Ideal (a - 1, b - 1, c^2 + 3*c + 1, d + c + 3, e - 1) of Multivariate Polynomial Ring in e, d, c, b, a over Rational Field,
                Ideal (a - 1, b^2 + 3*b + 1, c + b + 3, d - 1, e - 1) of Multivariate Polynomial Ring in e, d, c, b, a over Rational Field,
                Ideal (a - 1, b^4 + b^3 + b^2 + b + 1, c - b^2, d - b^3, e + b^3 + b^2 + b + 1) of Multivariate Polynomial Ring in e, d, c, b, a over Rational Field,
                Ideal (a^2 + 3*a + 1, b - 1, c - 1, d - 1, e + a + 3) of Multivariate Polynomial Ring in e, d, c, b, a over Rational Field,
                Ideal (a^2 + 3*a + 1, b + a + 3, c - 1, d - 1, e - 1) of Multivariate Polynomial Ring in e, d, c, b, a over Rational Field,
                Ideal (a^4 - 4*a^3 + 6*a^2 + a + 1, 11*b^2 - 6*b*a^3 + 26*b*a^2 - 41*b*a + 4*b + 8*a^3 - 31*a^2 + 40*a + 24, 11*c + 3*a^3 - 13*a^2 + 26*a - 2, 11*d + 3*a^3 - 13*a^2 + 26*a - 2, 11*e + 11*b - 6*a^3 + 26*a^2 - 41*a + 4) of Multivariate Polynomial Ring in e, d, c, b, a over Rational Field,
                Ideal (a^4 + a^3 + a^2 + a + 1, b - 1, c + a^3 + a^2 + a + 1, d - a^3, e - a^2) of Multivariate Polynomial Ring in e, d, c, b, a over Rational Field,
                Ideal (a^4 + a^3 + a^2 + a + 1, b - a, c - a, d^2 + 3*d*a + a^2, e + d + 3*a) of Multivariate Polynomial Ring in e, d, c, b, a over Rational Field,
                Ideal (a^4 + a^3 + a^2 + a + 1, b - a, c^2 + 3*c*a + a^2, d + c + 3*a, e - a) of Multivariate Polynomial Ring in e, d, c, b, a over Rational Field,
                Ideal (a^4 + a^3 + a^2 + a + 1, b^2 + 3*b*a + a^2, c + b + 3*a, d - a, e - a) of Multivariate Polynomial Ring in e, d, c, b, a over Rational Field,
                Ideal (a^4 + a^3 + a^2 + a + 1, b^3 + b^2*a + b^2 + b*a^2 + b*a + b + a^3 + a^2 + a + 1, c + b^2*a^3 + b^2*a^2 + b^2*a + b^2, d - b^2*a^2 - b^2*a - b^2 - b*a^2 - b*a - a^2, e - b^2*a^3 + b*a^2 + b*a + b + a^2 + a) of Multivariate Polynomial Ring in e, d, c, b, a over Rational Field,
                Ideal (a^4 + a^3 + 6*a^2 - 4*a + 1, 11*b^2 - 6*b*a^3 - 10*b*a^2 - 39*b*a - 2*b - 16*a^3 - 23*a^2 - 104*a + 24, 11*c + 3*a^3 + 5*a^2 + 25*a + 1, 11*d + 3*a^3 + 5*a^2 + 25*a + 1, 11*e + 11*b - 6*a^3 - 10*a^2 - 39*a - 2) of Multivariate Polynomial Ring in e, d, c, b, a over Rational Field]
            """

        P = self.ring()

        is_groebner = self.basis_is_groebner()

        # make sure to work w.r.t. 'lex'
        if P.term_order() != 'lex':
            Q = P.change_ring(order='lex')
        else:
            Q = P

        if is_groebner:
            if Q == P:
                I =  MPolynomialIdeal(P, self.interreduced_basis())
            else:
                I = self
                I = MPolynomialIdeal(P, I.transformed_basis('fglm')) # -> 'lex'
                I = I.change_ring(Q) # transform to 'lex' GB
        else:
            if Q == P:
                I = MPolynomialIdeal(P, self.groebner_basis())
            else:
                I = self.change_ring(Q) # transform to 'lex' GB
                I = MPolynomialIdeal(Q, I.groebner_basis())

        if I.dimension() != 0:
            raise TypeError, "dimension must be zero"

        Ibar = I._singular_()
        Ibar.attrib('isSB',1)

        singular = Ibar.parent()
        singular.lib("triang")

        if algorithm is None:
            algorithm = "singular:triangL"

        if algorithm == "singular:triangL":
            Tbar = Ibar.triangL()
        elif algorithm == "singular:triangLfak":
            Tbar = Ibar.triangLfak()
        elif algorithm == "singular:triangM":
            Tbar = Ibar.triangM()
        else:
            raise TypeError, "algorithm '%s' unknown"%algorithm

        T = Sequence([ MPolynomialIdeal(Q,[f._sage_(Q) for f in t]) for t in Tbar ])

        f = lambda x,y: cmp(x.gens(), y.gens())
        T.sort(f)

        return T

    def dimension(self, singular=singular_default):
        """
        The dimension of the ring modulo this ideal.

        EXAMPLE::

            sage: P.<x,y,z> = PolynomialRing(GF(32003),order='degrevlex')
            sage: I = ideal(x^2-y,x^3)
            sage: I.dimension()
            1

        For polynomials over a finite field of order too large for Singular,
        this falls back on a toy implementation of Buchberger to compute
        the Groebner basis, then uses the algorithm described in Chapter 9,
        Section 1 of Cox, Little, and O'Shea's "Ideals, Varieties, and Algorithms".

        EXAMPLE::

            sage: R.<x,y> = PolynomialRing(GF(2147483659),order='lex')
            sage: I = R.ideal([x*y,x*y+1])
            sage: I.dimension()
            verbose 0 (...: multi_polynomial_ideal.py, dimension) Warning: falling back to very slow toy implementation.
            0
            sage: I=ideal([x*(x*y+1),y*(x*y+1)])
            sage: I.dimension()
            verbose 0 (...: multi_polynomial_ideal.py, dimension) Warning: falling back to very slow toy implementation.
            1
            sage: I = R.ideal([x^3*y,x*y^2])
            sage: I.dimension()
            verbose 0 (...: multi_polynomial_ideal.py, dimension) Warning: falling back to very slow toy implementation.
            1
            sage: R.<x,y> = PolynomialRing(GF(2147483659),order='lex')
            sage: I = R.ideal(0)
            sage: I.dimension()
            verbose 0 (...: multi_polynomial_ideal.py, dimension) Warning: falling back to very slow toy implementation.
            2

        ALGORITHM: Uses Singular, unless the characteristic is too large.

        .. note::

           Requires computation of a Groebner basis, which can be a
           very expensive operation.
        """
        try:
            return self.__dimension
        except AttributeError:
            try:
                v = self._groebner_basis_singular_raw()
                self.__dimension = Integer(v.dim())
            except TypeError:
                if not self.base_ring().is_field():
                    raise NotImplementedError, "dimension() is implemented only over fields."
                if self.ring().term_order().is_global():
                    verbose("Warning: falling back to very slow toy implementation.", level=0)
                    # See Chapter 9, Section 1 of Cox, Little, O'Shea's "Ideals, Varieties,
                    # and Algorithms".
                    from sage.sets.set import Set
                    gb = toy_buchberger.buchberger_improved(self)
                    ring_vars = self.ring().gens()
                    n = len(ring_vars)
                    lms = [each.lm() for each in gb]
                    # compute M_j, denoted by var_lms
                    var_lms = [Set([]) for each in lms]
                    for j in xrange(len(ring_vars)):
                        for i in xrange(len(lms)):
                            if lms[i].degree(ring_vars[j]) > 0:
                                var_lms[i] += Set([j+1])
                    # compute intersections of M_j and J
                    # we assume that the iterator starts with the empty set,
                    # then iterates through all subsets of order 1,
                    # then through all subsets of order 2, etc...
                    # the way Sage currently operates
                    all_J = Set([each + 1 for each in range(n)]).subsets()
                    min_dimension = -1
                    all_J = all_J.iterator()
                    while min_dimension == -1:
                        try:
                            J = all_J.next()
                        except StopIteration:
                            min_dimension = n
                            break
                        J_intersects_all = True
                        i = 0
                        while J_intersects_all and i < len(var_lms):
                            J_intersects_all = J.intersection(var_lms[i]) != Set([])
                            i += 1
                        if J_intersects_all:
                            min_dimension = len(J)
                    return n - min_dimension
                else:
                    raise TypeError, "Local/unknown orderings not supported by 'toy_buchberger' implementation."
        return self.__dimension

    def vector_space_dimension(self):
        """
        Return the vector space dimension of the ring modulo this ideal. If
        the ideal is not zero-dimensional, a TypeError is raised.

        ALGORITHM: Uses Singular.

        EXAMPLE::

            sage: R.<u,v> = PolynomialRing(QQ)
            sage: g = u^4 + v^4 + u^3 + v^3
            sage: I = ideal(g) + ideal(g.gradient())
            sage: I.dimension()
            0
            sage: I.vector_space_dimension()
            4
        """
        gb = self.groebner_basis()

        vdim = Integer(singular_default.ideal(gb).vdim())

        if vdim == -1:
            raise TypeError, "ideal is not zero dimensional"
        else:
            return vdim

    @redSB
    def _groebner_basis_singular(self, algorithm="groebner", *args, **kwds):
        r"""
        Return the reduced Groebner basis of this ideal. If the Groebner
        basis for this ideal has been calculated before, the cached Groebner
        basis is returned regardless of the requested algorithm.

        INPUT:


        -  ``algorithm`` - see below for available algorithms


        ALGORITHMS:

        'groebner'
            use Singular's groebner heuristic to choose an algorithm (default)

        'std'
            Buchberger's algorithm

        'stdhilb'
            computes the standard basis of the homogeneous ideal in the
            basering, via a Hilbert driven standard basis computation.

        'stdfglm'
            computes the standard basis of the ideal in the basering via fglm
            (from the degrevlex ordering to the ordering of the basering).

        'slimgb'
            the *SlimGB* algorithm


        EXAMPLES:

        We compute a Groebner basis of 'cyclic 4' relative to lexicographic
        ordering.

        ::

            sage: R.<a,b,c,d> = PolynomialRing(QQ, 4, order='lex')
            sage: I = sage.rings.ideal.Cyclic(R,4); I
            Ideal (a + b + c + d, a*b + a*d + b*c + c*d, a*b*c + a*b*d
            + a*c*d + b*c*d, a*b*c*d - 1) of Multivariate Polynomial
            Ring in a, b, c, d over Rational Field

        ::

            sage: I._groebner_basis_singular()
            [c^2*d^6 - c^2*d^2 - d^4 + 1, c^3*d^2 + c^2*d^3 - c - d,
            b*d^4 - b + d^5 - d, b*c - b*d + c^2*d^4 + c*d - 2*d^2,
            b^2 + 2*b*d + d^2, a + b + c + d]

        ALGORITHM: Uses Singular.

        .. note::

           This method is called by the :meth:`.groebner_basis` method
           and the user usually doesn't need to bother with this one.
        """
        R = self.ring()
        S = self._groebner_basis_singular_raw(algorithm=algorithm, *args, **kwds)
        S =  Sequence([R(S[i+1]) for i in range(len(S))], R, check=False, immutable=True)
        return S

    def _groebner_basis_singular_raw(self, algorithm="groebner", singular=singular_default, *args, **kwds):
        r"""
        Return a Grobner basis in Singular format.

        EXAMPLE::

            sage: R.<a,b,c,d> = PolynomialRing(QQ, 4, order='lex')
            sage: I = sage.rings.ideal.Cyclic(R,4)
            sage: I._groebner_basis_singular() # indirect doctest
            [c^2*d^6 - c^2*d^2 - d^4 + 1, c^3*d^2 + c^2*d^3 - c - d,
            b*d^4 - b + d^5 - d, b*c - b*d + c^2*d^4 + c*d - 2*d^2,
            b^2 + 2*b*d + d^2, a + b + c + d]
        """
        try:
            return self.__gb_singular
        except AttributeError:
            pass
        # singular options are preserved by @redSB so we don't
        # need to do that here too
        for o,v in kwds.iteritems():
            if v:
                singular.option(o)
            else:
                singular.option("no"+o)

        if algorithm=="groebner":
            S = self._singular_().groebner()
        elif algorithm=="std":
            S = self._singular_().std()
        elif algorithm=="slimgb":
            S = self._singular_().slimgb()
        elif algorithm=="stdhilb":
            S = self._singular_().stdhilb()
        elif algorithm=="stdfglm":
            S = self._singular_().stdfglm()
        else:
            raise TypeError, "algorithm '%s' unknown"%algorithm
        self.__gb_singular = S
        return S

    def _groebner_basis_libsingular(self, algorithm="std"):
        r"""
        Return the reduced Groebner basis of this ideal. If the Groebner
        basis for this ideal has been calculated before the cached Groebner
        basis is returned regardless of the requested algorithm.

        INPUT:


        -  ``algorithm`` - see below for available algorithms


        ALGORITHMS:

        'std'
            Buchberger's algorithm

        'slimgb'
            the *SlimGB* algorithm


        EXAMPLES:

        We compute a Groebner basis of 'cyclic 4' relative to lexicographic
        ordering.

        ::

            sage: R.<a,b,c,d> = PolynomialRing(QQ, 4, order='lex')
            sage: I = sage.rings.ideal.Cyclic(R,4); I
            Ideal (a + b + c + d, a*b + a*d + b*c + c*d, a*b*c + a*b*d
            + a*c*d + b*c*d, a*b*c*d - 1) of Multivariate Polynomial
            Ring in a, b, c, d over Rational Field

        ::

            sage: I._groebner_basis_libsingular()
            [c^2*d^6 - c^2*d^2 - d^4 + 1, c^3*d^2 + c^2*d^3 - c - d,
            b*d^4 - b + d^5 - d, b*c - b*d + c^2*d^4 + c*d - 2*d^2,
            b^2 + 2*b*d + d^2, a + b + c + d]

        ALGORITHM: Uses libSINGULAR.
        """
        from sage.rings.polynomial.multi_polynomial_ideal_libsingular import std_libsingular, slimgb_libsingular

        if algorithm=="std":
            S = std_libsingular(self)
        elif algorithm=="slimgb":
            S = slimgb_libsingular(self)
        else:
            raise TypeError, "algorithm '%s' unknown"%algorithm
        return S

    def genus(self):
        """
        Return the genus of the projective curve defined by this ideal,
        which must be 1 dimensional.

        EXAMPLE: Consider the hyperelliptic curve
        `y^2 = 4x^5 - 30x^3 + 45x - 22` over
        `\mathbb{Q}`, it has genus 2::

            sage: P, x = PolynomialRing(QQ,"x").objgen()
            sage: f = 4*x^5 - 30*x^3 + 45*x - 22
            sage: C = HyperellipticCurve(f); C
            Hyperelliptic Curve over Rational Field defined by y^2 = 4*x^5 - 30*x^3 + 45*x - 22
            sage: C.genus()
            2

        ::

            sage: P.<x,y> = PolynomialRing(QQ)
            sage: f = y^2 - 4*x^5 - 30*x^3 + 45*x - 22
            sage: I = Ideal([f])
            sage: I.genus()
            2
        """
        try:
            return self.__genus
        except AttributeError:
            I = self._singular_()
            I.parent().lib('normal.lib')
            self.__genus = Integer(I.genus())
            return self.__genus

    @redSB
    def intersection(self, other):
        """
        Return the intersection of the two ideals.

        EXAMPLES::

            sage: R.<x,y> = PolynomialRing(QQ, 2, order='lex')
            sage: I = x*R
            sage: J = y*R
            sage: I.intersection(J)
            Ideal (x*y) of Multivariate Polynomial Ring in x, y over Rational Field

        The following simple example illustrates that the product need not
        equal the intersection.

        ::

            sage: I = (x^2, y)*R
            sage: J = (y^2, x)*R
            sage: K = I.intersection(J); K
            Ideal (y^2, x*y, x^2) of Multivariate Polynomial Ring in x, y over Rational Field
            sage: IJ = I*J; IJ
            Ideal (x^2*y^2, x^3, y^3, x*y) of Multivariate Polynomial Ring in x, y over Rational Field
            sage: IJ == K
            False
        """
        R = self.ring()
        if not isinstance(other, MPolynomialIdeal_singular_repr) or other.ring() != R:
            raise ValueError, "other must be an ideal in the ring of self, but it isn't."
        I = self._singular_()
        sing = I.parent()
        J = sing(other)
        K = I.intersect(J)
        return R.ideal(K)

    @redSB
    def minimal_associated_primes(self):
        r"""
        OUTPUT:


        -  ``list`` - a list of prime ideals


        EXAMPLES::

            sage: R.<x,y,z> = PolynomialRing(QQ, 3, 'xyz')
            sage: p = z^2 + 1; q = z^3 + 2
            sage: I = (p*q^2, y-z^2)*R
            sage: I.minimal_associated_primes ()
            [Ideal (z^2 + 1, -z^2 + y) of Multivariate Polynomial Ring
            in x, y, z over Rational Field, Ideal (z^3 + 2, -z^2 + y)
            of Multivariate Polynomial Ring in x, y, z over Rational
            Field]

        ALGORITHM: Uses Singular.
        """
        I = self._singular_()
        I.parent().lib('primdec.lib')
        M = I.minAssGTZ()
        R = self.ring()
        return [R.ideal(J) for J in M]

    @redSB
    def radical(self):
        r"""
        The radical of this ideal.

        EXAMPLES:

        This is an obviously not radical ideal::

            sage: R.<x,y,z> = PolynomialRing(QQ, 3)
            sage: I = (x^2, y^3, (x*z)^4 + y^3 + 10*x^2)*R
            sage: I.radical()
            Ideal (y, x) of Multivariate Polynomial Ring in x, y, z over Rational Field

        That the radical is correct is clear from the Groebner basis.

        ::

            sage: I.groebner_basis()
            [y^3, x^2]

        This is the example from the Singular manual::

            sage: p = z^2 + 1; q = z^3 + 2
            sage: I = (p*q^2, y-z^2)*R
            sage: I.radical()
            Ideal (z^2 - y, y^2*z + y*z + 2*y + 2) of Multivariate Polynomial Ring in x, y, z over Rational Field

        .. note::

           From the Singular manual: A combination of the algorithms
           of Krick/Logar and Kemper is used. Works also in positive
           characteristic (Kempers algorithm).

        ::

            sage: R.<x,y,z> = PolynomialRing(GF(37), 3)
            sage: p = z^2 + 1; q = z^3 + 2
            sage: I = (p*q^2, y - z^2)*R
            sage: I.radical()
            Ideal (z^2 - y, y^2*z + y*z + 2*y + 2) of Multivariate Polynomial Ring in x, y, z over Finite Field of size 37
        """
        S = self.ring()
        I = self._singular_()
        I.parent().lib('primdec.lib')
        r = I.radical()
        return S.ideal(r)

    @redSB
    def integral_closure(self, p=0, r=True, singular=singular_default):
        r"""
        Let `I` = ``self``.

        Returns the integral closure of `I, ..., I^p`, where
        `sI` is an ideal in the polynomial ring
        `R=k[x(1),...x(n)]`. If `p` is not given, or
        `p=0`, compute the closure of all powers up to the maximum
        degree in t occurring in the closure of `R[It]` (so this is
        the last power whose closure is not just the sum/product of the
        smaller). If `r` is given and ``r is True``,
        ``I.integral_closure()`` starts with a check whether I
        is already a radical ideal.

        INPUT:


        -  ``p`` - powers of I (default: 0)

        -  ``r`` - check whether self is a radical ideal first
           (default: True)


        EXAMPLE::

            sage: R.<x,y> = QQ[]
            sage: I = ideal([x^2,x*y^4,y^5])
            sage: I.integral_closure()
            [x^2, y^5, -x*y^3]

        ALGORITHM: Use Singular
        """
        Is = self._singular_()
        R = self.ring()
        singular =Is .parent()
        singular.load('reesclos.lib')
        ret = Sequence([ R(f) for f in Is.normalI(p,int(r))[1] ], R,
                       check=False, immutable=True)
        return ret

    def syzygy_module(self):
        r"""
        Computes the first syzygy (i.e., the module of relations of the
        given generators) of the ideal.

        EXAMPLE::

            sage: R.<x,y> = PolynomialRing(QQ)
            sage: f = 2*x^2 + y
            sage: g = y
            sage: h = 2*f + g
            sage: I = Ideal([f,g,h])
            sage: M = I.syzygy_module(); M
            [       -2        -1         1]
            [       -y 2*x^2 + y         0]
            sage: G = vector(I.gens())
            sage: M*G
            (0, 0)

        ALGORITHM: Uses Singular's syz command
        """
        return self._singular_().syz().transpose().sage_matrix(self.ring())

    @redSB
    def reduced_basis(self):
        r"""
        .. warning::

           This function is deprecated. It will be removed in a future
           release of Sage. Please use the :meth:`interreduced_basis`
           function instead.

        If this ideal is spanned by `(f_1, ..., f_n)` this method
        returns `(g_1, ..., g_s)` such that:

        - `(f_1,...,f_n) = (g_1,...,g_s)`

        - `LT(g_i) != LT(g_j)` for all `i != j`

        - `LT(g_i)` does not divide `m` for all monomials `m`
          of `\{g_1,...,g_{i-1},g_{i+1},...,g_s\}`

        - `LC(g_i) == 1` for all `i`.

        EXAMPLES::

            sage: R.<x,y,z> = PolynomialRing(QQ)
            sage: I = Ideal([z*x+y^3,z+y^3,z+x*y])
            sage: I.reduced_basis()
            doctest:...: DeprecationWarning: This function is deprecated. It will be removed in a future release of Sage. Please use the interreduced_basis() function instead.
            [x*z - z, x*y + z, y^3 + z]

            sage: R.<x,y,z> = PolynomialRing(QQ,order='negdegrevlex')
            sage: I = Ideal([z*x+y^3,z+y^3,z+x*y])
            sage: I.reduced_basis()
            [x*z + y^3, x*y - y^3, z + y^3]

        ALGORITHM:

        Uses Singular's interred command or
        ``toy_buchberger.inter_reduction`` if conversion to
        Singular fails.
        """
        from sage.misc.misc import deprecation
        deprecation("This function is deprecated. It will be removed in a future release of Sage. Please use the interreduced_basis() function instead.")
        return self.interreduced_basis()

    @redSB
    def interreduced_basis(self):
        r"""
        If this ideal is spanned by `(f_1, ..., f_n)` this method
        returns `(g_1, ..., g_s)` such that:


        -  `(f_1,...,f_n) = (g_1,...,g_s)`

        -  `LT(g_i) != LT(g_j)` for all `i != j`

        -  `LT(g_i)` does not divide `m` for all monomials
           `m` of `\{g_1,...,g_{i-1},g_{i+1},...,g_s\}`

        -  `LC(g_i) == 1` for all `i`.


        EXAMPLE::

            sage: R.<x,y,z> = PolynomialRing(QQ)
            sage: I = Ideal([z*x+y^3,z+y^3,z+x*y])
            sage: I.interreduced_basis()
            [x*z - z, x*y + z, y^3 + z]

        ::

            sage: R.<x,y,z> = PolynomialRing(QQ,order='negdegrevlex')
            sage: I = Ideal([z*x+y^3,z+y^3,z+x*y])
            sage: I.interreduced_basis()
            [x*z + y^3, x*y - y^3, z + y^3]

        ALGORITHM: Uses Singular's interred command or
        ``toy_buchberger.inter_reduction`` if conversion to
        Singular fails.
        """
        from sage.rings.polynomial.multi_polynomial_ideal_libsingular import interred_libsingular
        from sage.rings.polynomial.multi_polynomial_libsingular import MPolynomialRing_libsingular

        R = self.ring()

        if isinstance(R,MPolynomialRing_libsingular):
            return interred_libsingular(self)
        else:
            try:
                s = self._singular_().parent()
                o = s.option("get")
                s.option("redTail")
                ret = []
                for f in self._singular_().interred():
                    f = R(f)
                    ret.append(f.lc()^(-1)*f) # lead coeffs are not reduced by interred
                s.option("set",o)
            except TypeError:
                ret = toy_buchberger.inter_reduction(self.gens())

        ret = Sequence( ret, R, check=False, immutable=True)
        return ret

    def basis_is_groebner(self, singular=singular_default):
        r"""
        Returns ``True`` if the generators of
        ``self`` (``self.gens()``) form a Groebner
        basis.

        Let `I` be the set of generators of this ideal. The check
        is performed by trying to lift `Syz(LM(I))` to
        `Syz(I)` as `I` forms a Groebner basis if and only
        if for every element `S` in `Syz(LM(I))`:

        .. math::

            S \cdot G = \sum_{i=0}^{m} h_ig_i         \rightarrow_G 0.

        .

        ALGORITHM: Uses Singular

        EXAMPLE::

            sage: R.<a,b,c,d,e,f,g,h,i,j> = PolynomialRing(GF(127),10)
            sage: I = sage.rings.ideal.Cyclic(R,4)
            sage: I.basis_is_groebner()
            False
            sage: I2 = Ideal(I.groebner_basis())
            sage: I2.basis_is_groebner()
            True

        A more complicated example::

            sage: R.<U6,U5,U4,U3,U2, u6,u5,u4,u3,u2, h> = PolynomialRing(GF(7583))
            sage: l = [u6 + u5 + u4 + u3 + u2 - 3791*h, \
                       U6 + U5 + U4 + U3 + U2 - 3791*h, \
                       U2*u2 - h^2, U3*u3 - h^2, U4*u4 - h^2, \
                       U5*u4 + U5*u3 + U4*u3 + U5*u2 + U4*u2 + U3*u2 - 3791*U5*h - 3791*U4*h - 3791*U3*h - 3791*U2*h - 2842*h^2, \
                       U4*u5 + U3*u5 + U2*u5 + U3*u4 + U2*u4 + U2*u3 - 3791*u5*h - 3791*u4*h - 3791*u3*h - 3791*u2*h - 2842*h^2, \
                       U5*u5 - h^2, U4*U2*u3 + U5*U3*u2 + U4*U3*u2 + U3^2*u2 - 3791*U5*U3*h - 3791*U4*U3*h - 3791*U3^2*h - 3791*U5*U2*h \
                        - 3791*U4*U2*h + U3*U2*h - 3791*U2^2*h - 3791*U4*u3*h - 3791*U4*u2*h - 3791*U3*u2*h - 2843*U5*h^2 + 1897*U4*h^2 - 946*U3*h^2 - 947*U2*h^2 + 2370*h^3, \
                       U3*u5*u4 + U2*u5*u4 + U3*u4^2 + U2*u4^2 + U2*u4*u3 - 3791*u5*u4*h - 3791*u4^2*h - 3791*u4*u3*h - 3791*u4*u2*h + u5*h^2 - 2842*u4*h^2, \
                       U2*u5*u4*u3 + U2*u4^2*u3 + U2*u4*u3^2 - 3791*u5*u4*u3*h - 3791*u4^2*u3*h - 3791*u4*u3^2*h - 3791*u4*u3*u2*h + u5*u4*h^2 + u4^2*h^2 + u5*u3*h^2 - 2842*u4*u3*h^2, \
                       U5^2*U4*u3 + U5*U4^2*u3 + U5^2*U4*u2 + U5*U4^2*u2 + U5^2*U3*u2 + 2*U5*U4*U3*u2 + U5*U3^2*u2 - 3791*U5^2*U4*h - 3791*U5*U4^2*h - 3791*U5^2*U3*h \
                        + U5*U4*U3*h - 3791*U5*U3^2*h - 3791*U5^2*U2*h + U5*U4*U2*h + U5*U3*U2*h - 3791*U5*U2^2*h - 3791*U5*U3*u2*h - 2842*U5^2*h^2 + 1897*U5*U4*h^2 \
                        - U4^2*h^2 - 947*U5*U3*h^2 - U4*U3*h^2 - 948*U5*U2*h^2 - U4*U2*h^2 - 1422*U5*h^3 + 3791*U4*h^3, \
                       u5*u4*u3*u2*h + u4^2*u3*u2*h + u4*u3^2*u2*h + u4*u3*u2^2*h + 2*u5*u4*u3*h^2 + 2*u4^2*u3*h^2 + 2*u4*u3^2*h^2 + 2*u5*u4*u2*h^2 + 2*u4^2*u2*h^2 \
                        + 2*u5*u3*u2*h^2 + 1899*u4*u3*u2*h^2, \
                       U5^2*U4*U3*u2 + U5*U4^2*U3*u2 + U5*U4*U3^2*u2 - 3791*U5^2*U4*U3*h - 3791*U5*U4^2*U3*h - 3791*U5*U4*U3^2*h - 3791*U5*U4*U3*U2*h \
                        + 3791*U5*U4*U3*u2*h + U5^2*U4*h^2 + U5*U4^2*h^2 + U5^2*U3*h^2 - U4^2*U3*h^2 - U5*U3^2*h^2 - U4*U3^2*h^2 - U5*U4*U2*h^2 \
                        - U5*U3*U2*h^2 - U4*U3*U2*h^2 + 3791*U5*U4*h^3 + 3791*U5*U3*h^3 + 3791*U4*U3*h^3, \
                       u4^2*u3*u2*h^2 + 1515*u5*u3^2*u2*h^2 + u4*u3^2*u2*h^2 + 1515*u5*u4*u2^2*h^2 + 1515*u5*u3*u2^2*h^2 + u4*u3*u2^2*h^2 \
                        + 1521*u5*u4*u3*h^3 - 3028*u4^2*u3*h^3 - 3028*u4*u3^2*h^3 + 1521*u5*u4*u2*h^3 - 3028*u4^2*u2*h^3 + 1521*u5*u3*u2*h^3 + 3420*u4*u3*u2*h^3, \
                       U5^2*U4*U3*U2*h + U5*U4^2*U3*U2*h + U5*U4*U3^2*U2*h + U5*U4*U3*U2^2*h + 2*U5^2*U4*U3*h^2 + 2*U5*U4^2*U3*h^2 + 2*U5*U4*U3^2*h^2 \
                        + 2*U5^2*U4*U2*h^2 + 2*U5*U4^2*U2*h^2 + 2*U5^2*U3*U2*h^2 - 2*U4^2*U3*U2*h^2 - 2*U5*U3^2*U2*h^2 - 2*U4*U3^2*U2*h^2 \
                         - 2*U5*U4*U2^2*h^2 - 2*U5*U3*U2^2*h^2 - 2*U4*U3*U2^2*h^2 - U5*U4*U3*h^3 - U5*U4*U2*h^3 - U5*U3*U2*h^3 - U4*U3*U2*h^3]

            sage: Ideal(l).basis_is_groebner()
            False
            sage: gb = Ideal(l).groebner_basis()
            sage: Ideal(gb).basis_is_groebner()
            True

        .. note::

           From the Singular Manual for the reduce function we use in
           this method: 'The result may have no meaning if the second
           argument (``self``) is not a standard basis'. I (malb)
           believe this refers to the mathematical fact that the
           results may have no meaning if self is no standard basis,
           i.e., Singular doesn't 'add' any additional 'nonsense' to
           the result. So we may actually use reduce to determine if
           self is a Groebner basis.
        """
        self.ring()._singular_().set_ring()

        F = singular( self.gens(), "module" )
        LTF = singular( [f.lt() for f in self.gens()] , "module" )

        M = (F * LTF.syz()).reduce(self._singular_())

        for i in range(M.ncols()):
            if int(singular.eval("%s[1,%s+1]!=0"%(M.name(),i))):
                return False

        self._singular_().attrib('isSB',1)
        return True

    @redSB
    def transformed_basis(self, algorithm="gwalk", other_ring=None, singular=singular_default):
        r"""
        Returns a lex or ``other_ring`` Groebner Basis for
        this ideal.

        INPUT:


        -  ``algorithm`` - see below for options.

        -  ``other_ring`` - only valid for algorithm 'fglm',
           if provided conversion will be performed to this ring. Otherwise a
           lex Groebner basis will be returned.


        ALGORITHMS:

        fglm
            FGLM algorithm. The input ideal must be given with a reduced
            Groebner Basis of a zero-dimensional ideal

        gwalk
            Groebner Walk algorithm (*default*)

        awalk1
            'first alternative' algorithm

        awalk2
            'second alternative' algorithm

        twalk
            Tran algorithm

        fwalk
            Fractal Walk algorithm


        EXAMPLES::

            sage: R.<x,y,z> = PolynomialRing(QQ,3)
            sage: I = Ideal([y^3+x^2,x^2*y+x^2, x^3-x^2, z^4-x^2-y])
            sage: I = Ideal(I.groebner_basis())
            sage: S.<z,x,y> = PolynomialRing(QQ,3,order='lex')
            sage: J = Ideal(I.transformed_basis('fglm',S))
            sage: J
            Ideal (y^4 + y^3, x*y^3 - y^3, x^2 + y^3, z^4 + y^3 - y)
             of Multivariate Polynomial Ring in z, x, y over Rational Field

        ::

            sage: R.<z,y,x>=PolynomialRing(GF(32003),3,order='lex')
            sage: I=Ideal([y^3+x*y*z+y^2*z+x*z^3,3+x*y+x^2*y+y^2*z])
            sage: I.transformed_basis('gwalk')
            [y^9 - y^7*x^2 - y^7*x - y^6*x^3 - y^6*x^2 - 3*y^6 - 3*y^5*x - y^3*x^7
             - 3*y^3*x^6 - 3*y^3*x^5 - y^3*x^4 - 9*y^2*x^5 - 18*y^2*x^4 - 9*y^2*x^3
             - 27*y*x^3 - 27*y*x^2 - 27*x,
             z*x + 8297*y^8*x^2 + 8297*y^8*x + 3556*y^7 - 8297*y^6*x^4 + 15409*y^6*x^3 - 8297*y^6*x^2
             - 8297*y^5*x^5 + 15409*y^5*x^4 - 8297*y^5*x^3 + 3556*y^5*x^2 + 3556*y^5*x + 3556*y^4*x^3
             + 3556*y^4*x^2 - 10668*y^4 - 10668*y^3*x - 8297*y^2*x^9 - 1185*y^2*x^8 + 14224*y^2*x^7
             - 1185*y^2*x^6 - 8297*y^2*x^5 - 14223*y*x^7 - 10666*y*x^6 - 10666*y*x^5 - 14223*y*x^4
             + x^5 + 2*x^4 + x^3, z*y^2 + y*x^2 + y*x + 3]

        ALGORITHM: Uses Singular
        """
        from sage.rings.polynomial.polynomial_ring_constructor import PolynomialRing

        if self.basis_is_groebner():
            Is = self._singular_()
        else:
            Is = self._groebner_basis_singular_raw()

        R = self.ring()

        if algorithm in ("gwalk","awalk1","awalk2","twalk","fwalk"):
            singular.LIB("grwalk")
            gb = singular("%s(%s)"%(algorithm,Is.name()))
            return [R(f) for f in gb]
        elif algorithm == "fglm":
            Rs = self.ring()._singular_()

            # new ring
            if other_ring is None:
                nR = PolynomialRing(R.base_ring(),R.ngens(), names=R.variable_names(), order="lex")
            else:
                nR = other_ring
            nR._singular_().set_ring()

            nIs = singular.fglm(Rs,Is)

            return [nR(f) for f in nIs]

        else:
            raise TypeError, "Cannot convert basis with given algorithm"

    @redSB
    def elimination_ideal(self, variables):
        r"""
        Returns the elimination ideal this ideal with respect to the
        variables given in ``variables``.

        INPUT:


        -  ``variables`` - a list or tuple of variables in
           ``self.ring()``


        EXAMPLE::

            sage: R.<x,y,t,s,z> = PolynomialRing(QQ,5)
            sage: I = R * [x-t,y-t^2,z-t^3,s-x+y^3]
            sage: I.elimination_ideal([t,s])
            Ideal (y^2 - x*z, x*y - z, x^2 - y) of Multivariate
            Polynomial Ring in x, y, t, s, z over Rational Field

        ALGORITHM: Uses SINGULAR

        .. note::

           Requires computation of a Groebner basis, which can be a very
           expensive operation.
        """
        if not isinstance(variables, (list,tuple)):
            variables = (variables,)

        Is = self._groebner_basis_singular_raw()
        R = self.ring()
        return MPolynomialIdeal(R, [f.sage_poly(R) for f in Is.eliminate( prod(variables) ) ] )

    @redSB
    def quotient(self, J):
        r"""
        Given ideals `I` = ``self`` and `J` in
        the same polynomial ring `P`, return the ideal quotient of
        `I` by `J` consisting of the polynomials `a` of
        `P` such that `\{aJ \subset I\}`.

        This is also referred to as the colon ideal
        (`I`:`J`).

        INPUT:


        -  ``J`` - multivariate polynomial ideal


        EXAMPLE::

            sage: R.<x,y,z> = PolynomialRing(GF(181),3)
            sage: I = Ideal([x^2+x*y*z,y^2-z^3*y,z^3+y^5*x*z])
            sage: J = Ideal([x])
            sage: Q = I.quotient(J)
            sage: y*z + x in I
            False
            sage: x in J
            True
            sage: x * (y*z + x) in I
            True
        """
        R = self.ring()

        if not isinstance(J, MPolynomialIdeal):
            raise TypeError, "J needs to be a multivariate polynomial ideal"

        if not R is J.ring() and not R == J.ring():
            raise TypeError, "base rings do not match"

        return R.ideal([f.sage_poly(R) for f in self._singular_().quotient(J._singular_())])

    def variety(self, ring=None, proof=True):
        r"""
        Return the variety of ``self``.

        Given a zero-dimensional ideal `I` (== ``self``) of a
        polynomial ring P whose order is lexicographic, return the
        variety of I as a list of dictionaries with (variable, value)
        pairs. By default, the variety of the ideal over its
        coefficient field `K` is returned; ``ring`` can be specified
        to find the variety over a different ring.

        These dictionaries have cardinality equal to the number of
        variables in P and represent assignments of values to these
        variables such that all polynomials in I vanish.

        If ``ring`` is specified, then a triangular decomposition of
        ``self`` is found over the original coefficient field `K`;
        then the triangular systems are solved using root-finding over
        ``ring``. This is particularly useful when `K` is ``QQ`` (to
        allow fast symbolic computation of the triangular
        decomposition) and ``ring`` is ``RR``, ``AA``, ``CC``, or
        ``QQbar`` (to compute the whole real or complex variety of the
        ideal).

        Note that with ``ring``=``RR`` or ``CC``, computation is done
        numerically and potentially inaccurately; in particular, the
        number of points in the real variety may be miscomputed. With
        ``ring``=``AA`` or ``QQbar``, computation is done exactly
        (which may be much slower, of course).

        INPUT:

        - ``ring`` - return roots in the ``ring`` instead of the base
          ring of this ideal (default: ``None``)
        - ``proof`` - return a provably correct result (default: ``True``)

        EXAMPLE::

            sage: K.<w> = GF(27) # this example is from the MAGMA handbook
            sage: P.<x, y> = PolynomialRing(K, 2, order='lex')
            sage: I = Ideal([ x^8 + y + 2, y^6 + x*y^5 + x^2 ])
            sage: I = Ideal(I.groebner_basis()); I
            Ideal (x - y^47 - y^45 + y^44 - y^43 + y^41 - y^39 - y^38
            - y^37 - y^36 + y^35 - y^34 - y^33 + y^32 - y^31 + y^30 +
            y^28 + y^27 + y^26 + y^25 - y^23 + y^22 + y^21 - y^19 -
            y^18 - y^16 + y^15 + y^13 + y^12 - y^10 + y^9 + y^8 + y^7
            - y^6 + y^4 + y^3 + y^2 + y - 1, y^48 + y^41 - y^40 + y^37
            - y^36 - y^33 + y^32 - y^29 + y^28 - y^25 + y^24 + y^2 + y
            + 1) of Multivariate Polynomial Ring in x, y over Finite
            Field in w of size 3^3

            sage: V = I.variety(); V
            [{y: w^2 + 2, x: 2*w}, {y: w^2 + w, x: 2*w + 1}, {y: w^2 + 2*w, x: 2*w + 2}]

            sage: [f.subs(v) for f in I.gens() for v in V] # check that all polynomials vanish
            [0, 0, 0, 0, 0, 0]

        However, we only account for solutions in the ground field and not
        in the algebraic closure::

            sage: I.vector_space_dimension()
            48

        Here we compute the points of intersection of a hyperbola and a
        circle, in several fields::

            sage: K.<x, y> = PolynomialRing(QQ, 2, order='lex')
            sage: I = Ideal([ x*y - 1, (x-2)^2 + (y-1)^2 - 1])
            sage: I = Ideal(I.groebner_basis()); I
            Ideal (x + y^3 - 2*y^2 + 4*y - 4, y^4 - 2*y^3 + 4*y^2 - 4*y + 1)
            of Multivariate Polynomial Ring in x, y over Rational Field

        These two curves have one rational intersection::

            sage: I.variety()
            [{y: 1, x: 1}]

        There are two real intersections::

            sage: I.variety(ring=RR)
            [{y: 0.361103080528647, x: 2.76929235423863},
             {y: 1.00000000000000, x: 1.00000000000000}]
            sage: I.variety(ring=AA)
            [{x: 2.769292354238632?, y: 0.3611030805286474?},
             {x: 1, y: 1}]

        and a total of four intersections::

            sage: I.variety(ring=CC)
            [{y: 0.31944845973567... - 1.6331702409152...*I,
              x: 0.11535382288068... + 0.58974280502220...*I},
             {y: 0.31944845973567... + 1.6331702409152...*I,
              x: 0.11535382288068... - 0.58974280502220...*I},
             {y: 0.36110308052864..., x: 2.7692923542386...},
             {y: 1.00000000000000, x: 1.00000000000000}]
            sage: I.variety(ring=QQbar)
            [{y: 0.3194484597356763? - 1.633170240915238?*I,
              x: 0.11535382288068429? + 0.5897428050222055?*I},
             {y: 0.3194484597356763? + 1.633170240915238?*I,
              x: 0.11535382288068429? - 0.5897428050222055?*I},
             {y: 0.3611030805286474?, x: 2.769292354238632?},
             {y: 1, x: 1}]

        If the ground field's characteristic is too large for
        Singular, we resort to a toy implementation::

            sage: R.<x,y> = PolynomialRing(GF(2147483659),order='lex')
            sage: I=ideal([x^3-2*y^2,3*x+y^4])
            sage: I.variety()
            verbose 0 (...: multi_polynomial_ideal.py, variety) Warning: falling back to very slow toy implementation.
            verbose 0 (...: multi_polynomial_ideal.py, groebner_basis) Warning: falling back to very slow toy implementation.
            [{y: 0, x: 0}]

        TESTS::

            sage: K.<w> = GF(27)
            sage: P.<x, y> = PolynomialRing(K, 2, order='lex')
            sage: I = Ideal([ x^8 + y + 2, y^6 + x*y^5 + x^2 ])

        Testing the robustness of the Singular interface

        ::

            sage: T = I.triangular_decomposition('singular:triangLfak')
            sage: I.variety()
            [{y: w^2 + 2, x: 2*w}, {y: w^2 + w, x: 2*w + 1}, {y: w^2 + 2*w, x: 2*w + 2}]

        Testing that a bug is indeed fixed.

        ::

            sage: R = PolynomialRing(GF(2), 30, ['x%d'%(i+1) for i in range(30)], order='lex')
            sage: R.inject_variables()
            Defining...
            sage: I = Ideal([x1 + 1, x2, x3 + 1, x5*x10 + x10 + x18, x5*x11 + x11, \
                             x5*x18, x6, x7 + 1, x9, x10*x11 + x10 + x18, x10*x18 + x18, \
                             x11*x18, x12, x13, x14, x15, x16 + 1, x17 + x18 + 1, x19, x20, \
                             x21 + 1, x22, x23, x24, x25 + 1, x28 + 1, x29 + 1, x30, x8, \
                             x26, x1^2 + x1, x2^2 + x2, x3^2 + x3, x4^2 + x4, x5^2 + x5, \
                             x6^2 + x6, x7^2 + x7, x8^2 + x8, x9^2 + x9, x10^2 + x10, \
                             x11^2 + x11, x12^2 + x12, x13^2 + x13, x14^2 + x14, x15^2 + x15, \
                             x16^2 + x16, x17^2 + x17, x18^2 + x18, x19^2 + x19, x20^2 + x20, \
                             x21^2 + x21, x22^2 + x22, x23^2 + x23, x24^2 + x24, x25^2 + x25, \
                             x26^2 + x26, x27^2 + x27, x28^2 + x28, x29^2 + x29, x30^2 + x30])
            sage: I.basis_is_groebner()
            True
            sage: for V in I.variety():
            ...     print V
            {x14: 0, x24: 0, x16: 1, x1: 1, x3: 1, x2: 0, x5: 0, x4: 0, x19: 0, x18: 0, x7: 1, x6: 0, x10: 0, x30: 0, x28: 1, x29: 1, x13: 0, x27: 0, x11: 0, x25: 1, x9: 0, x8: 0, x20: 0, x17: 1, x23: 0, x26: 0, x15: 0, x21: 1, x12: 0, x22: 0}
            {x14: 0, x24: 0, x16: 1, x1: 1, x3: 1, x2: 0, x5: 0, x4: 1, x19: 0, x18: 0, x7: 1, x6: 0, x10: 0, x30: 0, x28: 1, x29: 1, x13: 0, x27: 0, x11: 0, x25: 1, x9: 0, x8: 0, x20: 0, x17: 1, x23: 0, x26: 0, x15: 0, x21: 1, x12: 0, x22: 0}
            {x14: 0, x24: 0, x16: 1, x1: 1, x3: 1, x2: 0, x5: 1, x4: 0, x19: 0, x18: 0, x7: 1, x6: 0, x10: 0, x30: 0, x28: 1, x29: 1, x13: 0, x27: 0, x11: 0, x25: 1, x9: 0, x8: 0, x20: 0, x17: 1, x23: 0, x26: 0, x15: 0, x21: 1, x12: 0, x22: 0}
            {x14: 0, x24: 0, x16: 1, x1: 1, x3: 1, x2: 0, x5: 1, x4: 1, x19: 0, x18: 0, x7: 1, x6: 0, x10: 0, x30: 0, x28: 1, x29: 1, x13: 0, x27: 0, x11: 0, x25: 1, x9: 0, x8: 0, x20: 0, x17: 1, x23: 0, x26: 0, x15: 0, x21: 1, x12: 0, x22: 0}
            {x14: 0, x24: 0, x16: 1, x1: 1, x3: 1, x2: 0, x5: 1, x4: 0, x19: 0, x18: 0, x7: 1, x6: 0, x10: 0, x30: 0, x28: 1, x29: 1, x13: 0, x27: 0, x11: 1, x25: 1, x9: 0, x8: 0, x20: 0, x17: 1, x23: 0, x26: 0, x15: 0, x21: 1, x12: 0, x22: 0}
            {x14: 0, x24: 0, x16: 1, x1: 1, x3: 1, x2: 0, x5: 1, x4: 1, x19: 0, x18: 0, x7: 1, x6: 0, x10: 0, x30: 0, x28: 1, x29: 1, x13: 0, x27: 0, x11: 1, x25: 1, x9: 0, x8: 0, x20: 0, x17: 1, x23: 0, x26: 0, x15: 0, x21: 1, x12: 0, x22: 0}
            {x14: 0, x24: 0, x16: 1, x1: 1, x3: 1, x2: 0, x5: 1, x4: 0, x19: 0, x18: 0, x7: 1, x6: 0, x10: 1, x30: 0, x28: 1, x29: 1, x13: 0, x27: 0, x11: 1, x25: 1, x9: 0, x8: 0, x20: 0, x17: 1, x23: 0, x26: 0, x15: 0, x21: 1, x12: 0, x22: 0}
            {x14: 0, x24: 0, x16: 1, x1: 1, x3: 1, x2: 0, x5: 1, x4: 1, x19: 0, x18: 0, x7: 1, x6: 0, x10: 1, x30: 0, x28: 1, x29: 1, x13: 0, x27: 0, x11: 1, x25: 1, x9: 0, x8: 0, x20: 0, x17: 1, x23: 0, x26: 0, x15: 0, x21: 1, x12: 0, x22: 0}
            {x14: 0, x24: 0, x16: 1, x1: 1, x3: 1, x2: 0, x5: 0, x4: 0, x19: 0, x18: 1, x7: 1, x6: 0, x10: 1, x30: 0, x28: 1, x29: 1, x13: 0, x27: 0, x11: 0, x25: 1, x9: 0, x8: 0, x20: 0, x17: 0, x23: 0, x26: 0, x15: 0, x21: 1, x12: 0, x22: 0}
            {x14: 0, x24: 0, x16: 1, x1: 1, x3: 1, x2: 0, x5: 0, x4: 1, x19: 0, x18: 1, x7: 1, x6: 0, x10: 1, x30: 0, x28: 1, x29: 1, x13: 0, x27: 0, x11: 0, x25: 1, x9: 0, x8: 0, x20: 0, x17: 0, x23: 0, x26: 0, x15: 0, x21: 1, x12: 0, x22: 0}
            {x14: 0, x24: 0, x16: 1, x1: 1, x3: 1, x2: 0, x5: 0, x4: 0, x19: 0, x18: 0, x7: 1, x6: 0, x10: 0, x30: 0, x28: 1, x29: 1, x13: 0, x27: 1, x11: 0, x25: 1, x9: 0, x8: 0, x20: 0, x17: 1, x23: 0, x26: 0, x15: 0, x21: 1, x12: 0, x22: 0}
            {x14: 0, x24: 0, x16: 1, x1: 1, x3: 1, x2: 0, x5: 0, x4: 1, x19: 0, x18: 0, x7: 1, x6: 0, x10: 0, x30: 0, x28: 1, x29: 1, x13: 0, x27: 1, x11: 0, x25: 1, x9: 0, x8: 0, x20: 0, x17: 1, x23: 0, x26: 0, x15: 0, x21: 1, x12: 0, x22: 0}
            {x14: 0, x24: 0, x16: 1, x1: 1, x3: 1, x2: 0, x5: 1, x4: 0, x19: 0, x18: 0, x7: 1, x6: 0, x10: 0, x30: 0, x28: 1, x29: 1, x13: 0, x27: 1, x11: 0, x25: 1, x9: 0, x8: 0, x20: 0, x17: 1, x23: 0, x26: 0, x15: 0, x21: 1, x12: 0, x22: 0}
            {x14: 0, x24: 0, x16: 1, x1: 1, x3: 1, x2: 0, x5: 1, x4: 1, x19: 0, x18: 0, x7: 1, x6: 0, x10: 0, x30: 0, x28: 1, x29: 1, x13: 0, x27: 1, x11: 0, x25: 1, x9: 0, x8: 0, x20: 0, x17: 1, x23: 0, x26: 0, x15: 0, x21: 1, x12: 0, x22: 0}
            {x14: 0, x24: 0, x16: 1, x1: 1, x3: 1, x2: 0, x5: 1, x4: 0, x19: 0, x18: 0, x7: 1, x6: 0, x10: 0, x30: 0, x28: 1, x29: 1, x13: 0, x27: 1, x11: 1, x25: 1, x9: 0, x8: 0, x20: 0, x17: 1, x23: 0, x26: 0, x15: 0, x21: 1, x12: 0, x22: 0}
            {x14: 0, x24: 0, x16: 1, x1: 1, x3: 1, x2: 0, x5: 1, x4: 1, x19: 0, x18: 0, x7: 1, x6: 0, x10: 0, x30: 0, x28: 1, x29: 1, x13: 0, x27: 1, x11: 1, x25: 1, x9: 0, x8: 0, x20: 0, x17: 1, x23: 0, x26: 0, x15: 0, x21: 1, x12: 0, x22: 0}
            {x14: 0, x24: 0, x16: 1, x1: 1, x3: 1, x2: 0, x5: 1, x4: 0, x19: 0, x18: 0, x7: 1, x6: 0, x10: 1, x30: 0, x28: 1, x29: 1, x13: 0, x27: 1, x11: 1, x25: 1, x9: 0, x8: 0, x20: 0, x17: 1, x23: 0, x26: 0, x15: 0, x21: 1, x12: 0, x22: 0}
            {x14: 0, x24: 0, x16: 1, x1: 1, x3: 1, x2: 0, x5: 1, x4: 1, x19: 0, x18: 0, x7: 1, x6: 0, x10: 1, x30: 0, x28: 1, x29: 1, x13: 0, x27: 1, x11: 1, x25: 1, x9: 0, x8: 0, x20: 0, x17: 1, x23: 0, x26: 0, x15: 0, x21: 1, x12: 0, x22: 0}
            {x14: 0, x24: 0, x16: 1, x1: 1, x3: 1, x2: 0, x5: 0, x4: 0, x19: 0, x18: 1, x7: 1, x6: 0, x10: 1, x30: 0, x28: 1, x29: 1, x13: 0, x27: 1, x11: 0, x25: 1, x9: 0, x8: 0, x20: 0, x17: 0, x23: 0, x26: 0, x15: 0, x21: 1, x12: 0, x22: 0}
            {x14: 0, x24: 0, x16: 1, x1: 1, x3: 1, x2: 0, x5: 0, x4: 1, x19: 0, x18: 1, x7: 1, x6: 0, x10: 1, x30: 0, x28: 1, x29: 1, x13: 0, x27: 1, x11: 0, x25: 1, x9: 0, x8: 0, x20: 0, x17: 0, x23: 0, x26: 0, x15: 0, x21: 1, x12: 0, x22: 0}



        ALGORITHM: Uses triangular decomposition.
        """
        def _variety(T, V, v=None):
            """
            Return variety ``V`` for one triangular set of
            polynomials ``T``.
            """
            if v is None: v = {}
            found = False
            for f in T:
                if f.is_univariate() and not f.is_constant():
                    T.remove(f); found = True; break

            if found is False:
                V.append(v)
                return V

            variable = f.variable(0)
            roots = f.univariate_polynomial().roots(ring=ring)

            for root,_ in roots:
                vbar = v.copy()
                vbar[variable] = root
                Tbar = [ f.subs({variable:root}) for f in T ]
                _variety(Tbar,V,vbar)

            return V

        P = self.ring()
        if ring is not None: P = P.change_ring(ring)
        try:
          TI = self.triangular_decomposition('singular:triangLfak')
          T = [list(each.gens()) for each in TI]
        except TypeError, msg: # conversion to Singular not supported
          if self.ring().term_order().is_global():
            verbose("Warning: falling back to very slow toy implementation.", level=0)
            T = toy_variety.triangular_factorization(self.groebner_basis(),proof=proof)
          else:
            raise TypeError, "Local/unknown orderings not supported by 'toy_buchberger' implementation."

        V = []
        for t in T:
            Vbar = _variety(list(t),[])
            #Vbar = _variety(list(t.gens()),[])

            for v in Vbar:
                V.append(dict([(P(var),val) for var,val in v.iteritems()]))
        V.sort()
        return Sequence(V)

    def hilbert_polynomial(self):
        r"""
        Return the Hilbert polynomial of this ideal.

        Let `I` = ``self`` be a homogeneous ideal and
        `R` = ``self.ring()`` be a graded commutative
        algebra (`R = \oplus R_d`) over a field `K`. The
        Hilbert polynomial is the unique polynomial `HP(t)` with
        rational coefficients such that `HP(d) = dim_K R_d` for
        all but finitely many positive integers `d`.

        EXAMPLE::

            sage: P.<x,y,z> = PolynomialRing(QQ)
            sage: I = Ideal([x^3*y^2 + 3*x^2*y^2*z + y^3*z^2 + z^5])
            sage: I.hilbert_polynomial()
            5*t - 5
        """
        if not self.is_homogeneous():
            raise TypeError, "Ideal must be homogeneous."

        hp = self._singular_().hilbPoly()
        t = ZZ['t'].gen()
        fp = ZZ(len(hp)-1).factorial()
        return sum([ZZ(hp[i+1])*t**i for i in xrange(len(hp))])/fp

    def hilbert_series(self, singular=singular_default):
        r"""
        Return the Hilbert series of this ideal.

        Let `I` = ``self`` be a homogeneous ideal and
        `R` = ``self.ring()`` be a graded commutative
        algebra (`R = \oplus R_d`) over a field
        `K`. Then the Hilbert function is defined as
        `H(d) = dim_K R_d` and the Hilbert series of `I`
        is defined as the formal power series
        `HS(t) = \sum_0^{\infty} H(d) t^d`.

        This power series can be expressed as
        `HS(t) = Q(t)/(1-t)^n` where `Q(t)` is a polynomial
        over `Z` and `n` the number of variables in
        `R`. This method returns `Q(t)/(1-t)^n`.

        EXAMPLE::

            sage: P.<x,y,z> = PolynomialRing(QQ)
            sage: I = Ideal([x^3*y^2 + 3*x^2*y^2*z + y^3*z^2 + z^5])
            sage: I.hilbert_series()
            (-t^4 - t^3 - t^2 - t - 1)/(-t^2 + 2*t - 1)
        """
        if not self.is_homogeneous():
            raise TypeError, "Ideal must be homogeneous."

        gb = self.groebner_basis()
        t = ZZ['t'].gen()
        n = self.ring().ngens()
        hs = singular.ideal(gb).hilb(1)
        return sum([ZZ(hs[i+1])*t**i for i in xrange(len(hs)-1)])/(1-t)**n

class MPolynomialIdeal_macaulay2_repr:
    """
    An ideal in a multivariate polynomial ring, which has an underlying
    Macaulay2 ring associated to it.

    EXAMPLES::

        sage: R.<x,y,z,w> = PolynomialRing(ZZ, 4)
        sage: I = ideal(x*y-z^2, y^2-w^2)
        sage: I
        Ideal (x*y - z^2, y^2 - w^2) of Multivariate Polynomial Ring in x, y, z, w over Integer Ring
    """
    def _macaulay2_(self, macaulay2=None):
        """
        Return Macaulay2 ideal corresponding to this ideal.
        """
        if macaulay2 is None: macaulay2 = macaulay2_default
        try:
            I = self.__macaulay2[macaulay2]
            I._check_valid()
            return I
        except KeyError:
            pass
        except AttributeError:
            self.__macaulay2 = {}
        except ValueError:
            pass

        R = self.ring()
        R._macaulay2_set_ring(macaulay2)

        gens = [repr(x) for x in self.gens()]
        if len(gens) == 0:
            gens = ['0']
        z = macaulay2.ideal(gens)
        self.__macaulay2[macaulay2] = z
        return z

    def _groebner_basis_macaulay2(self):
        r"""
        Return the Groebner basis for this ideal, computed using
        Macaulay2.

        ALGORITHM: Computed using Macaulay2. A big advantage of Macaulay2
        is that it can compute the Groebner basis of ideals in polynomial rings
        over the integers.

        EXAMPLE::

            sage: R.<x,y,z,w> = PolynomialRing(ZZ, 4)
            sage: I = ideal(x*y-z^2, y^2-w^2)
            sage: I.groebner_basis('macaulay2')                # indirect doctest; optional - macaulay2
            [z^4 - x^2*w^2, y*z^2 - x*w^2, x*y - z^2, y^2 - w^2]

        The Groebner basis can be used to compute in
        `\mathbb{Z}/n\mathbb{Z}[x,\ldots]`.

        ::

            sage: R.<x,y,z> = ZZ[]
            sage: I = ideal([y^2*z - x^3 - 19*x*z, y^2, 19^2])
            sage: I.groebner_basis('macaulay2')               # optional - macaulay2
            [x^3 + 19*x*z, y^2, 361]
            sage: I = ideal([y^2*z - x^3 - 19^2*x*z, y^2, 19^2])
            sage: I.groebner_basis('macaulay2')               # optional - macaulay2
            [x^3, y^2, 361]
        """
        I = self._macaulay2_()
        G = str(I.gb().generators().external_string()).replace('\n','')
        i = G.rfind('{{')
        j = G.rfind('}}')
        G = G[i+2:j].split(',')
        R = self.ring()
        B = [R(f) for f in G]
        B = Sequence(B, self.ring(), check=False)
        B.set_immutable()
        return B

    def _reduce_using_macaulay2(self, f):
        """
        EXAMPLES::

            sage: R.<x,y,z,w> = PolynomialRing(ZZ, 4)
            sage: I = ideal(x*y-z^2, y^2-w^2)
            sage: I._reduce_using_macaulay2(x*y-z^2 + y^2)    # optional  - macaulay2
            w^2
        """
        I = self._macaulay2_()
        M2 = I.parent()
        k = M2('(%r) %% %s'%(f, I.name()))
        R = self.ring()
        return R(k)


class MPolynomialIdeal( MPolynomialIdeal_singular_repr, \
                        MPolynomialIdeal_macaulay2_repr, \
                        MPolynomialIdeal_magma_repr, \
                        Ideal_generic ):
    def __init__(self, ring, gens, coerce=True):
        r"""
        Create an ideal in a multivariate polynomial ring.

        INPUT:


        -  ``ring`` - the ring the ideal is defined in

        -  ``gens`` - a list of generators for the ideal

        -  ``coerce`` - coerce elements to the ring
           ``ring``?


        EXAMPLES::

            sage: R.<x,y> = PolynomialRing(IntegerRing(), 2, order='lex')
            sage: R.ideal([x, y])
            Ideal (x, y) of Multivariate Polynomial Ring in x, y over Integer Ring
            sage: R.<x0,x1> = GF(3)[]
            sage: R.ideal([x0^2, x1^3])
            Ideal (x0^2, x1^3) of Multivariate Polynomial Ring in x0, x1 over Finite Field of size 3
        """
        Ideal_generic.__init__(self, ring, gens, coerce=coerce)

    def __cmp__(self, other):
        """
        EXAMPLE::

            sage: R = PolynomialRing(QQ,'x,y,z')
            sage: I = R.ideal()
            sage: I == R.ideal()
            True

        ::

            sage: R = PolynomialRing(QQ, names=[])
            sage: R.ideal(0) == R.ideal(0)
            verbose 0 (...: multi_polynomial_ideal.py, groebner_basis) Warning: falling back to very slow toy implementation.
            verbose 0 (...: multi_polynomial_ideal.py, groebner_basis) Warning: falling back to very slow toy implementation.
            True

        ::

            sage: R, (x,y) = PolynomialRing(QQ, 2, 'xy').objgens()
            sage: I = (x^3 + y, y)*R
            sage: J = (x^3 + y, y, y*x^3 + y^2)*R
            sage: I == J
            True

        ::

            sage: R = PolynomialRing(QQ, 'x,y,z', order='degrevlex')
            sage: S = PolynomialRing(QQ, 'x,y,z', order='invlex')
            sage: I = R.ideal([R.0,R.1])
            sage: J = S.ideal([S.0,S.1])
            sage: I == J
            True
            sage: cmp(I,J)
            0
            sage: I.__cmp__(J)
            0
        """
        # first check the type
        if not isinstance(other, MPolynomialIdeal):
            return 1

        # the ideals may be defined w.r.t. to different term orders
        # but are still the same.
        R = self.ring()
        S = other.ring()
        if R is not S: # rings are unique
            if type(R) == type(S) and (R.base_ring() == S.base_ring()) and (R.ngens() == S.ngens()):
                other = other.change_ring(R)
            else:
                return cmp((type(R), R.base_ring(), R.ngens()), (type(S), S.base_ring(), S.ngens()))

        # now, check whether the GBs are cached already
        if self.groebner_basis.is_in_cache() and other.groebner_basis.is_in_cache():
            l = self.groebner_basis()
            r = other.groebner_basis()
        else: # use easy GB otherwise
            try:
                l = self.change_ring(R.change_ring(order="degrevlex")).groebner_basis()
                r = other.change_ring(R.change_ring(order="degrevlex")).groebner_basis()
            except AttributeError: # e.g. quotient rings
                l = self.groebner_basis()
                r = other.groebner_basis()
        return cmp(l,r)

    def groebner_fan(self, is_groebner_basis=False, symmetry=None, verbose=False):
        r"""
        Return the Groebner fan of this ideal.

        The base ring must be `\mathbb{Q}` or a finite field
        `\mathbb{F}_p` of with `p <= 32749`.

        EXAMPLES::

            sage: P.<x,y> = PolynomialRing(QQ)
            sage: i = ideal(x^2 - y^2 + 1)
            sage: g = i.groebner_fan()
            sage: g.reduced_groebner_bases()
            [[x^2 - y^2 + 1], [-x^2 + y^2 - 1]]

        INPUT:


        -  ``is_groebner_basis`` - bool (default False). if
           True, then I.gens() must be a Groebner basis with respect to the
           standard degree lexicographic term order.

        -  ``symmetry`` - default: None; if not None, describes
           symmetries of the ideal

        -  ``verbose`` - default: False; if True, printout
           useful info during computations
        """
        import sage.rings.polynomial.groebner_fan as groebner_fan
        return groebner_fan.GroebnerFan(self, is_groebner_basis=is_groebner_basis,
                                        symmetry=symmetry, verbose=verbose)

    @cached_method
    def groebner_basis(self, algorithm='', *args, **kwds):
        r"""
        Return the reduced Groebner basis of this ideal. A Groeber basis
        `g_1,...,g_n` for an ideal `I` is a basis such
        that


        .. math::

            <LT(g_i)> = LT(I)

        ,

        i.e. the leading term ideal of `I` is spanned by the
        leading terms of `g_1,...,g_n`. Groebner bases are the
        key concept in computational ideal theory in multivariate
        polynomial rings which allows a variety of problems to be solved.
        Additionally, a *reduced* Groebner basis `G` is a unique
        representation for the ideal `<G>` with respect to the
        chosen monomial ordering.

        INPUT:


        -  ``algorithm`` - determines the algorithm to use, see
           below for available algorithms.

        -  ``*args`` - additional parameters passed to the
           respective implementations

        -  ``**kwds`` - additional keyword parameters passed
           to the respective implementations


        ALGORITHMS:

        ''
            autoselect (default)

        'singular:groebner'
            Singular's ``groebner`` command

        'singular:std'
            Singular's ``std`` command

        'singular:stdhilb'
            Singular's ``stdhib`` command

        'singular:stdfglm'
            Singular's ``stdfglm`` command

        'singular:slimgb'
            Singular's ``slimgb`` command

        'libsingular:std'
            libSingular's ``std`` command

        'libsingular:slimgb'
            libSingular's ``slimgb`` command

        'toy:buchberger'
            Sage's toy/educational buchberger without strategy

        'toy:buchberger2'
            Sage's toy/educational buchberger with strategy

        'toy:d_basis'
            Sage's toy/educational d_basis algorithm

        'macaulay2:gb'
            Macaulay2's ``gb`` command (if available)

        'magma:GroebnerBasis'
            Magma's ``Groebnerbasis`` command (if available)


        If only a system is given - e.g. 'magma' - the default algorithm is
        chosen for that system.

        .. note::

           The Singular and libSingular versions of the respective
           algorithms are identical, but the former calls an external
           Singular process while the later calls a C function, i.e. the
           calling overhead is smaller.

        EXAMPLES:

        Consider Katsura-3 over QQ with lexicographical term
        ordering. We compute the reduced Groebner basis using every
        available implementation and check their equality.

        ::

            sage: P.<a,b,c> = PolynomialRing(QQ,3, order='lex')
            sage: I = sage.rings.ideal.Katsura(P,3) # regenerate to prevent caching
            sage: I.groebner_basis()
            [a - 60*c^3 + 158/7*c^2 + 8/7*c - 1, b + 30*c^3 - 79/7*c^2 + 3/7*c, c^4 - 10/21*c^3 + 1/84*c^2 + 1/84*c]

        ::

            sage: I = sage.rings.ideal.Katsura(P,3) # regenerate to prevent caching
            sage: I.groebner_basis('singular:groebner')
            [a - 60*c^3 + 158/7*c^2 + 8/7*c - 1, b + 30*c^3 - 79/7*c^2 + 3/7*c, c^4 - 10/21*c^3 + 1/84*c^2 + 1/84*c]

        ::

            sage: I = sage.rings.ideal.Katsura(P,3) # regenerate to prevent caching
            sage: I.groebner_basis('singular:std')
            [a - 60*c^3 + 158/7*c^2 + 8/7*c - 1, b + 30*c^3 - 79/7*c^2 + 3/7*c, c^4 - 10/21*c^3 + 1/84*c^2 + 1/84*c]

        ::

            sage: I = sage.rings.ideal.Katsura(P,3) # regenerate to prevent caching
            sage: I.groebner_basis('singular:stdhilb')
            [a - 60*c^3 + 158/7*c^2 + 8/7*c - 1, b + 30*c^3 - 79/7*c^2 + 3/7*c, c^4 - 10/21*c^3 + 1/84*c^2 + 1/84*c]

        ::

            sage: I = sage.rings.ideal.Katsura(P,3) # regenerate to prevent caching
            sage: I.groebner_basis('singular:stdfglm')
            [a - 60*c^3 + 158/7*c^2 + 8/7*c - 1, b + 30*c^3 - 79/7*c^2 + 3/7*c, c^4 - 10/21*c^3 + 1/84*c^2 + 1/84*c]

        ::

            sage: I = sage.rings.ideal.Katsura(P,3) # regenerate to prevent caching
            sage: I.groebner_basis('singular:slimgb')
            [a - 60*c^3 + 158/7*c^2 + 8/7*c - 1, b + 30*c^3 - 79/7*c^2 + 3/7*c, c^4 - 10/21*c^3 + 1/84*c^2 + 1/84*c]

        Note that toy:buchberger does not return the reduced Groebner
        basis,

        ::

            sage: I = sage.rings.ideal.Katsura(P,3) # regenerate to prevent caching
            sage: I.groebner_basis('toy:buchberger')
            [a^2 - a + 2*b^2 + 2*c^2,
             a*b + b*c - 1/2*b, a + 2*b + 2*c - 1,
             b^2 + 4/3*b*c - 1/3*b + c^2 - 1/3*c,
             b*c - 1/10*b + 6/5*c^2 - 2/5*c,
             b + 30*c^3 - 79/7*c^2 + 3/7*c,
             c^4 - 10/21*c^3 + 1/84*c^2 + 1/84*c]

        but that toy:buchberger2 does.

        ::

            sage: I = sage.rings.ideal.Katsura(P,3) # regenerate to prevent caching
            sage: I.groebner_basis('toy:buchberger2')
            [a - 60*c^3 + 158/7*c^2 + 8/7*c - 1, b + 30*c^3 - 79/7*c^2 + 3/7*c, c^4 - 10/21*c^3 + 1/84*c^2 + 1/84*c]

        ::

            sage: I = sage.rings.ideal.Katsura(P,3) # regenerate to prevent caching
            sage: I.groebner_basis('macaulay2:gb') # optional - macaulay2
            [a - 60*c^3 + 158/7*c^2 + 8/7*c - 1, b + 30*c^3 - 79/7*c^2 + 3/7*c, c^4 - 10/21*c^3 + 1/84*c^2 + 1/84*c]

        ::

            sage: I = sage.rings.ideal.Katsura(P,3) # regenerate to prevent caching
            sage: I.groebner_basis('magma:GroebnerBasis') # optional - magma
            [a - 60*c^3 + 158/7*c^2 + 8/7*c - 1, b + 30*c^3 - 79/7*c^2 + 3/7*c, c^4 - 10/21*c^3 + 1/84*c^2 + 1/84*c]

        Groebner bases over `\mathbb{Z}` can be computed. However,
        the native implementation is very slow. If available Macaulay2 is
        used, which is an optional Sage package.

        ::

            sage: P.<a,b,c> = PolynomialRing(ZZ,3)
            sage: I = P * (a + 2*b + 2*c - 1, a^2 - a + 2*b^2 + 2*c^2, 2*a*b + 2*b*c - b)
            sage: I.groebner_basis()                               # optional - macaulay2
            ...
            [b^3 - b*c^2 - 24*c^3 - b^2 - 6*b*c + 2*c^2 + 2*c,
             2*b*c^2 + 36*c^3 - 9*b*c - 24*c^2 + 2*b + 4*c,
             42*c^3 + b^2 + 2*b*c - 14*c^2 + b,
             2*b^2 - 4*b*c - 6*c^2 + 2*c,
             10*b*c + 12*c^2 - b - 4*c,
             a + 2*b + 2*c - 1]

        ::

            sage: I.groebner_basis('macaulay2') # optional - macaulay2
            [b^3 + b*c^2 + 12*c^3 + b^2 + b*c - 4*c^2,
             2*b*c^2 - 6*c^3 + b^2 + 5*b*c + 8*c^2 - b - 2*c,
             42*c^3 + b^2 + 2*b*c - 14*c^2 + b,
             2*b^2 - 4*b*c - 6*c^2 + 2*c, 10*b*c + 12*c^2 - b - 4*c,
             a + 2*b + 2*c - 1]

        Sage also supports local orderings::

            sage: P.<x,y,z> = PolynomialRing(QQ,3,order='negdegrevlex')
            sage: I = P * (  x*y*z + z^5, 2*x^2 + y^3 + z^7, 3*z^5 +y ^5 )
            sage: I.groebner_basis()
            [x^2 + 1/2*y^3, x*y*z + z^5, y^5 + 3*z^5, y^4*z - 2*x*z^5, z^6]

        ALGORITHM: Uses Singular, Magma (if available), Macaulay2 (if
        available), or toy implementation.
        """
        from sage.rings.integer_mod_ring import is_IntegerModRing

        if algorithm.lower() == "magma":
            algorithm = "magma:GroebnerBasis"
        elif algorithm.lower() == "singular":
            algorithm = "singular:groebner"
        elif algorithm.lower() == "macaulay2":
            algorithm = "macaulay2:gb"
        elif algorithm.lower() == "toy":
            algorithm = "toy:buchberger2"

        if algorithm is '':
            if self.ring().base_ring() is ZZ:
#                 try:
#                     gb = self._groebner_basis_macaulay2()
#                 except TypeError:
                verbose("Warning: falling back to very slow toy implementation.", level=0)
                gb = toy_d_basis.d_basis(self, *args, **kwds)
            elif is_IntegerModRing(self.ring().base_ring()) and not self.ring().base_ring().is_field():
                try:
                    gb = self._groebner_basis_macaulay2()
                except TypeError:
                    verbose("Warning: falling back to very slow toy implementation.", level=0)

                    ch = self.ring().base_ring().characteristic()
                    R = self.ring().change_ring(ZZ)
                    I = R.ideal([R(f) for f in self.gens()+(ch,)])

                    gb = toy_d_basis.d_basis(I, *args, **kwds)

                    R = self.ring()
                    gb = filter(lambda f: f,[R(f) for f in gb])
            else:
                try:
                    gb = self._groebner_basis_singular("groebner", *args, **kwds)
                except TypeError, msg: # conversion to Singular not supported
                    if self.ring().term_order().is_global():
                        verbose("Warning: falling back to very slow toy implementation.", level=0)
                        gb = toy_buchberger.buchberger_improved(self, *args, **kwds)
                    else:
                        raise TypeError, "Local/unknown orderings not supported by 'toy_buchberger' implementation."
        elif algorithm.startswith('singular:'):
            gb = self._groebner_basis_singular(algorithm[9:])
        elif algorithm.startswith('libsingular:'):
            gb = self._groebner_basis_libsingular(algorithm[len('libsingular:'):], *args, **kwds)
        elif algorithm == 'macaulay2:gb':
            gb = self._groebner_basis_macaulay2(*args, **kwds)
        elif algorithm == 'magma:GroebnerBasis':
            gb = self._groebner_basis_magma(*args, **kwds)
        elif algorithm == 'toy:buchberger':
            gb = toy_buchberger.buchberger(self, *args, **kwds)
        elif algorithm == 'toy:buchberger2':
            gb = toy_buchberger.buchberger_improved(self, *args, **kwds)
        elif algorithm == 'toy:d_basis':
            gb = toy_d_basis.d_basis(self, *args, **kwds)
        else:
            raise TypeError, "algorithm '%s' unknown"%algorithm

        gb = sorted(gb, reverse=True)
        if self.ring().base_ring().is_field():
            gb = Sequence( [f*f.lc()**(-1) for f in gb], immutable=True, check=False)
        return gb

    def change_ring(self, P):
        r"""
        Return the ideal ``I`` in ``P`` spanned by
        the generators `g_1, ..., g_n` of self as returned by
        ``self.gens()``.

        INPUT:


        -  ``P`` - a multivariate polynomial ring


        EXAMPLE::

            sage: P.<x,y,z> = PolynomialRing(QQ,3,order='lex')
            sage: I = sage.rings.ideal.Cyclic(P)
            sage: I
            Ideal (x + y + z, x*y + x*z + y*z, x*y*z - 1) of
            Multivariate Polynomial Ring in x, y, z over Rational Field

        ::

            sage: I.groebner_basis()
            [x + y + z, y^2 + y*z + z^2, z^3 - 1]

        ::

            sage: Q.<x,y,z> = P.change_ring(order='degrevlex'); Q
            Multivariate Polynomial Ring in x, y, z over Rational Field
            sage: Q.term_order()
            Degree reverse lexicographic term order

        ::

            sage: J = I.change_ring(Q); J
            Ideal (x + y + z, x*y + x*z + y*z, x*y*z - 1) of
            Multivariate Polynomial Ring in x, y, z over Rational Field

        ::

            sage: J.groebner_basis()
            [z^3 - 1, y^2 + y*z + z^2, x + y + z]
        """
        return P.ideal([P(f) for f in self.gens()])

    def reduce(self, f):
        """
        Reduce an element modulo the reduced Groebner basis for this ideal.
        This returns 0 if and only if the element is in this ideal. In any
        case, this reduction is unique up to monomial orders.

        EXAMPLES::

            sage: R.<x,y> = PolynomialRing(QQ, 2)
            sage: I = (x^3 + y, y)*R
            sage: I.reduce(y)
            0
            sage: I.reduce(x^3)
            0
            sage: I.reduce(x - y)
            x

        ::

            sage: I = (y^2 - (x^3 + x))*R
            sage: I.reduce(x^3)
            y^2 - x
            sage: I.reduce(x^6)
            y^4 - 2*x*y^2 + x^2
            sage: (y^2 - x)^2
            y^4 - 2*x*y^2 + x^2

        .. note::

           Requires computation of a Groebner basis, which can be a very
           expensive operation.
        """
        gb = self.groebner_basis()
        return f.reduce(gb)

    def _contains_(self, f):
        r"""
        Returns ``True`` if ``f`` is in this ideal,
        ``False`` otherwise.

        EXAMPLES::

            sage: R, (x,y) = PolynomialRing(QQ, 2, 'xy').objgens()
            sage: I = (x^3 + y, y)*R
            sage: x in I # indirect doctest
            False
            sage: y in I
            True
            sage: x^3 + 2*y in I
            True

        .. note::

          Requires computation of a Groebner basis, which can be a very
          expensive operation.
        """
        g = f.reduce(self.groebner_basis())
        return self.ring()(g).is_zero()

    def homogenize(self, var='h'):
        """
        Return homogeneous ideal spanned by the homogeneous polynomials
        generated by homogenizing the generators of this ideal.

        INPUT:


        -  ``h`` - variable name or variable in cover ring
           (default: 'h')


        EXAMPLE::

            sage: P.<x,y,z> = PolynomialRing(GF(2))
            sage: I = Ideal([x^2*y + z + 1, x + y^2 + 1]); I
            Ideal (x^2*y + z + 1, y^2 + x + 1) of Multivariate
            Polynomial Ring in x, y, z over Finite Field of size 2

        ::

            sage: I.homogenize()
            Ideal (x^2*y + z*h^2 + h^3, y^2 + x*h + h^2) of
            Multivariate Polynomial Ring in x, y, z, h over Finite
            Field of size 2

        ::

            sage: I.homogenize(y)
            Ideal (x^2*y + y^3 + y^2*z, x*y) of Multivariate
            Polynomial Ring in x, y, z over Finite Field of size 2

        ::

                   sage: I = Ideal([x^2*y + z^3 + y^2*x, x + y^2 + 1])
            sage: I.homogenize()
            Ideal (x^2*y + x*y^2 + z^3, y^2 + x*h + h^2) of
            Multivariate Polynomial Ring in x, y, z, h over Finite
            Field of size 2
        """
        I = [f.homogenize(var) for f in self.gens()]
        P = max(I, key=lambda x: x.parent().ngens()).parent()
        return P.ideal([P(f) for f in I])

    def is_homogeneous(self):
        r"""
        Return ``True`` if this ideal is spanned by homogeneous
        polynomials, i.e. if it is a homogeneous ideal.

        EXAMPLE::

            sage: P.<x,y,z> = PolynomialRing(QQ,3)
            sage: I = sage.rings.ideal.Katsura(P)
            sage: I
            Ideal (x + 2*y + 2*z - 1, x^2 + 2*y^2 + 2*z^2 - x, 2*x*y +
            2*y*z - y) of Multivariate Polynomial Ring in x, y, z over
            Rational Field

        ::

            sage: I.is_homogeneous()
            False

        ::

            sage: J = I.homogenize()
            sage: J
            Ideal (x + 2*y + 2*z - h, x^2 + 2*y^2 + 2*z^2 - x*h, 2*x*y
            + 2*y*z - y*h) of Multivariate Polynomial Ring in x, y, z,
            h over Rational Field

        ::

            sage: J.is_homogeneous()
            True
        """
        for f in self.gens():
            if not f.is_homogeneous():
                return False
        return True

    def _normal_basis_libsingular(self):
        r"""
        Returns the normal basis for a given groebner basis. It will use
        the Groebner Basis as computed by
        ``MPolynomialIdeal._groebner_basis_libsingular()``.

        EXAMPLES::

            sage: R.<x,y,z> = PolynomialRing(QQ)
            sage: I = R.ideal(x^2-2*x*z+5, x*y^2+y*z+1, 3*y^2-8*x*z)
            sage: I.normal_basis()
            [z^2, y*z, x*z, z, x*y, y, x, 1]
        """
        from sage.rings.polynomial.multi_polynomial_ideal_libsingular import kbase_libsingular
        gb = self._groebner_basis_libsingular()

        return kbase_libsingular(self.ring().ideal(gb))

    def normal_basis(self, algorithm='libsingular', singular=singular_default):
        """
        Returns a vector space basis (consisting of monomials) of the
        quotient ring by the ideal, resp. of a free module by the module,
        in case it is finite dimensional and if the input is a standard
        basis with respect to the ring ordering.

        INPUT: algorithm - defaults to use libsingular, if it is anything
        else we will use the kbase() command

        EXAMPLES::

            sage: R.<x,y,z> = PolynomialRing(QQ)
            sage: I = R.ideal(x^2+y^2+z^2-4, x^2+2*y^2-5, x*z-1)
            sage: I.normal_basis()
            [y*z^2, z^2, y*z, z, x*y, y, x, 1]
            sage: I.normal_basis(algorithm='singular')
            [y*z^2, z^2, y*z, z, x*y, y, x, 1]
        """

        if algorithm == 'libsingular':
            return self._normal_basis_libsingular()
        else:
            gb = self.groebner_basis()
            return list(singular.kbase(self.ring().ideal(gb)))

    def plot(self, *args, **kwds):
        """
        Plot the real zero locus of this ideal (if principal).

        INPUT:


        -  ``self`` - must be a principal ideal in 2 variables
           over QQ.

        -  ``algorithm`` - set this to 'surf' if you want
           'surf' to plot the ideal (default: None)

        -  ``*args`` - optional tuples ``(variable,
           minimum, maximum)`` for plotting dimensions

        -  ``**kwds`` - optional keyword arguments passed on
           to ``implicit_plot``


        EXAMPLES:

        Implicit plotting in 2-d::

            sage: R.<x,y> = PolynomialRing(QQ,2)
            sage: I = R.ideal([y^3 - x^2])
            sage: I.plot()                         # cusp

        ::

            sage: I = R.ideal([y^2 - x^2 - 1])
            sage: I.plot((x,-3, 3), (y, -2, 2))  # hyperbola

        ::

            sage: I = R.ideal([y^2 + x^2*(1/4) - 1])
            sage: I.plot()                         # ellipse

        ::

            sage: I = R.ideal([y^2-(x^2-1)*(x-2)])
            sage: I.plot()                         # elliptic curve

        ::

            sage: f = ((x+3)^3 + 2*(x+3)^2 - y^2)*(x^3 - y^2)*((x-3)^3-2*(x-3)^2-y^2)
            sage: I = R.ideal(f)
            sage: I.plot()                         # the Singular logo

        This used to be trac #5267::

            sage: I = R.ideal([-x^2*y+1])
            sage: I.plot()

        AUTHORS:

        - Martin Albrecht (2008-09)
        """
        from sage.rings.rational_field import QQ
        from sage.rings.real_mpfr import RR
        from sage.plot.all import implicit_plot

        if self.ring().base_ring() is not QQ:
            raise TypeError, "Base ring must be rational field."
        if not self.is_principal():
            raise TypeError, "Ideal must be principal."


        f = self.gens()[0]

        variables = sorted(f.variables(), reverse=True)

        if len(variables) == 2 and kwds.get('algorithm','') != 'surf':
            V = [(variables[0], None, None), (variables[1], None, None)]

            if len(args) > 2:
                raise TypeError, "Expected up to 2 optional parameters but got %d."%len(args)

            # first check whether user supplied boundaries
            for e in args:
                if not isinstance(e, (tuple, list)) or len(e) != 3:
                    raise TypeError, "Optional parameter must be list or tuple or length 3."
                v,mi,ma = e

                if v not in variables:
                    raise TypeError, "Optional parameter must contain variable of ideal generator."

                vi = variables.index(v)
                V[vi] = v,mi,ma

            # now check whether we should find boundaries
            for var_index in range(2):
                if V[var_index][1] is None:
                    v, mi, ma = variables[var_index], -10, 10
                    for i in range(mi, ma):
                        roots = f.subs({v:i}).univariate_polynomial().change_ring(RR).roots()
                        if len(roots) > 0:
                            mi = i - 1
                            break

                    for i in range(ma, mi, -1):
                        roots = f.subs({v:i}).univariate_polynomial().change_ring(RR).roots()
                        if len(roots) > 0:
                            ma = i + 1
                            break
                    V[var_index] = variables[var_index], mi, ma

            kwds.setdefault("plot_points",200)
            return implicit_plot(f, V[0], V[1], **kwds)

        elif len(variables) == 3 or kwds.get('algorithm','') == 'surf':
            MPolynomialIdeal_singular_repr.plot(self, kwds.get("singular",singular_default))
        else:
            raise TypeError, "Ideal generator may not have either 2 or 3 variables."


    def weil_restriction(self):
        """
        Compute the Weil restriction of this ideal over some extension
        field.

        A Weil restriction of scalars - denoted `Res_{L/k}` - is a
        functor which, for any finite extension of fields `L/k` and
        any algebraic variety `X` over `L`, produces another
        corresponding variety `Res_{L/k}(X)`, defined over `k`. It is
        useful for reducing questions about varieties over large
        fields to questions about more complicated varieties over
        smaller fields.

        This function does not compute this Weil restriction directly
        but computes on generating sets of polynomial ideals:

        Let `d` be the degree of the field extension `L/k`, let `a` a
        generator of `L/k` and `p` the minimal polynomial of
        `L/k`. Denote this ideal by `I`.

        Specifically, this function first maps each variable `x` to
        its representation over `k`: `\sum_{i=0}^{d-1} a^i x_i`. Then
        each generator of `I` is evaulated over these representations
        and reduced modulo the minimal polynomial `p`. The result is
        interpreted as a univariate polynomial in `a` and its
        coefficients are the new generators of the returned ideal.

        If the input and the output ideals are radical, this is
        equivalent to the statement about algebraic varieties above.

        EXAMPLE::

            sage: k.<a> = GF(2^2)
            sage: P.<x,y> = PolynomialRing(k,2)
            sage: I = Ideal([x*y + 1, a*x + 1])
            sage: I.variety()
            [{y: a, x: a + 1}]
            sage: J = I.weil_restriction()
            sage: J
            Ideal (x1*y0 + x0*y1 + x1*y1, x0*y0 + x1*y1 + 1, x0 + x1, x1 + 1) of
            Multivariate Polynomial Ring in x0, x1, y0, y1 over Finite Field of size 2
            sage: J += sage.rings.ideal.FieldIdeal(J.ring()) # ensure radical ideal
            sage: J.variety()
            [{y1: 1, x1: 1, x0: 1, y0: 0}]

            sage: J.weil_restriction() # returns J
            Ideal (x1*y0 + x0*y1 + x1*y1, x0*y0 + x1*y1 + 1, x0 + x1, x1 + 1, x0^2 + x0,
                   x1^2 + x1, y0^2 + y0, y1^2 + y1)
                   of Multivariate Polynomial Ring in x0, x1, y0, y1 over Finite Field of size 2

            sage: k.<a> = GF(3^5)
            sage: P.<x,y,z> = PolynomialRing(k)
            sage: I = sage.rings.ideal.Katsura(P)
            sage: I.dimension()
            0
            sage: I.variety()
            [{y: 0, z: 0, x: 1}]

            sage: J = I.weil_restriction(); J
            Ideal (x4 - y4 - z4, x3 - y3 - z3, x2 - y2 - z2, x1 - y1 - z1, x0 - y0 - z0 - 1,
                   x2^2 - x1*x3 - x0*x4 + x4^2 - y2^2 + y1*y3 + y0*y4 - y4^2 - z2^2 + z1*z3 + z0*z4 - z4^2 - x4,
                   -x1*x2 - x0*x3 - x3*x4 - x4^2 + y1*y2 + y0*y3 + y3*y4 + y4^2 + z1*z2 + z0*z3 + z3*z4 + z4^2 - x3,
                   x1^2 - x0*x2 + x3^2 - x2*x4 + x3*x4 - y1^2 + y0*y2 - y3^2 + y2*y4 - y3*y4 - z1^2 + z0*z2 - z3^2 + z2*z4 - z3*z4 - x2,
                   -x0*x1 - x2*x3 - x3^2 - x1*x4 + x2*x4 + y0*y1 + y2*y3 + y3^2 + y1*y4 - y2*y4 + z0*z1 + z2*z3 + z3^2 + z1*z4 - z2*z4 - x1,
                   x0^2 + x2*x3 + x1*x4 - y0^2 - y2*y3 - y1*y4 - z0^2 - z2*z3 - z1*z4 - x0,
                   -x4*y0 - x3*y1 - x2*y2 - x1*y3 - x0*y4 - x4*y4 - y4*z0 - y3*z1 - y2*z2 - y1*z3 - y0*z4 - y4*z4 - y4,
                   -x3*y0 - x2*y1 - x1*y2 - x0*y3 - x4*y3 - x3*y4 + x4*y4 - y3*z0 - y2*z1 - y1*z2 - y0*z3 - y4*z3 - y3*z4 + y4*z4 - y3,
                   -x2*y0 - x1*y1 - x0*y2 - x4*y2 - x3*y3 + x4*y3 - x2*y4 + x3*y4 - y2*z0 - y1*z1 - y0*z2 - y4*z2 - y3*z3 + y4*z3 - y2*z4 + y3*z4 - y2,
                   -x1*y0 - x0*y1 - x4*y1 - x3*y2 + x4*y2 - x2*y3 + x3*y3 - x1*y4 + x2*y4 - y1*z0 - y0*z1 - y4*z1 - y3*z2 + y4*z2 - y2*z3 + y3*z3 - y1*z4 + y2*z4 - y1,
                   -x0*y0 + x4*y1 + x3*y2 + x2*y3 + x1*y4 - y0*z0 + y4*z1 + y3*z2 + y2*z3 + y1*z4 - y0) of Multivariate Polynomial Ring in
                   x0, x1, x2, x3, x4, y0, y1, y2, y3, y4, z0, z1, z2, z3, z4 over Finite Field of size 3
            sage: J += sage.rings.ideal.FieldIdeal(J.ring()) # ensure radical ideal
            sage: J.variety()
            [{y1: 0, y4: 0, z2: 0, y2: 0, x0: 1, y0: 0, x2: 0, z4: 0, z3: 0, x4: 0, x1: 0, z1: 0, z0: 0, y3: 0, x3: 0}]


        Weil restrictions are often used to study elliptic curves over
        extension fields so we give a simple example involving those::

            sage: K.<a> = QuadraticField(1/3)
            sage: E = EllipticCurve(K,[1,2,3,4,5])

        We pick a point on ``E``::

            sage: p = E.lift_x(1); p
            (1 : 2 : 1)

            sage: I = E.defining_ideal(); I
            Ideal (-x^3 - 2*x^2*z + x*y*z + y^2*z - 4*x*z^2 + 3*y*z^2 - 5*z^3)
            of Multivariate Polynomial Ring in x, y, z over Number Field in a with defining polynomial x^2 - 1/3

        Of course, the point ``p`` is a root of all generators of ``I``::

            sage: [f.subs(x=1,y=2,z=1) for f in I.gens()]
            [0]

        ``I`` is also radical::

            sage: I.radical() == I
            True

        So we compute its Weil restriction::

            sage: J = I.weil_restriction()
            sage: J
            Ideal (-3*x0^2*x1 - 1/3*x1^3 - 4*x0*x1*z0 + x1*y0*z0 + x0*y1*z0 + 2*y0*y1*z0 - 4*x1*z0^2 + 3*y1*z0^2 - 2*x0^2*z1
                    - 2/3*x1^2*z1 + x0*y0*z1 + y0^2*z1 + 1/3*x1*y1*z1 + 1/3*y1^2*z1 - 8*x0*z0*z1 + 6*y0*z0*z1 - 15*z0^2*z1 - 4/3*x1*z1^2 + y1*z1^2 - 5/3*z1^3,
                   -x0^3 - x0*x1^2 - 2*x0^2*z0 - 2/3*x1^2*z0 + x0*y0*z0 + y0^2*z0 + 1/3*x1*y1*z0 + 1/3*y1^2*z0 - 4*x0*z0^2 + 3*y0*z0^2 - 5*z0^3 - 4/3*x0*x1*z1
                    + 1/3*x1*y0*z1 + 1/3*x0*y1*z1 + 2/3*y0*y1*z1 - 8/3*x1*z0*z1 + 2*y1*z0*z1 - 4/3*x0*z1^2 + y0*z1^2 - 5*z0*z1^2)
            of Multivariate Polynomial Ring in x0, x1, y0, y1, z0, z1 over Rational Field

        We can check that the point ``p`` is still a root of all generators of ``J``::

            sage: [f.subs(x0=1,y0=2,z0=1,x1=0,y1=0,z1=0) for f in J.gens()]
            [0, 0]

        .. note::

           Based on a Singular implementation by Michael Brickenstein
        """
        from sage.rings.polynomial.polynomial_ring_constructor import PolynomialRing

        R = self.ring()
        nvars = R.ngens()
        L = R.base_ring()
        k = L.prime_subfield()

        if L.degree() == 1:
            return self

        helper = PolynomialRing(L.prime_subfield(), nvars + 1, (L.variable_name(),) + R.variable_names(), order='lex')
        myminpoly = L.polynomial().subs(helper.gen(0))

        l = [helper(str(f))  for f in self.gens()]

        r = myminpoly.degree()

        intermediate_ring = PolynomialRing(k, nvars*r+1, 'x')

        a = intermediate_ring.gen(0)

        # map e.g. x -> a^2*x_2 + a*x_1 + x_0, where x_0,..,x_2
        # represent the components of x if viewed as a vector in k^r
        map_ideal = [a]

        variables = iter(intermediate_ring.gens()[1:])
        for _ in xrange(nvars):
           map_ideal.append(sum([a**i * variables.next() for i in range(r)]))

        myminpoly = myminpoly(*map_ideal)

        l = [f(*map_ideal).reduce([myminpoly]) for f in l]

        result = []
        # split e.g. a^2*x0+a*x1+x2 to x0,x1,x2
        for f in l:
            for i in reversed(range(r)):
               g = f.coefficient(a**i)
               f =  f - a**i * g
               result.append(g)

        # eliminate parameter
        new_var_names = [str(var)+"%d"%i for var in R.gens() for i in range(r)]

        result_ring = PolynomialRing(k, nvars*r, new_var_names)

        map_ideal = (0,) + result_ring.gens()
        result = [f(*map_ideal) for f in result]

        return result_ring.ideal(result)



