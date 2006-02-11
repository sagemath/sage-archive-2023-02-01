"""
Ideals in multivariate polynomial rings.

AUTHOR:
    -- William Stein

EXAMPLES:
    sage: x,y,z = QQ['x,y,z'].gens()
    sage: I = ideal(x^5 + y^4 + z^3 - 1,  x^3 + y^3 + z^2 - 1)
    sage: B = I.groebner_basis()
    sage: len(B)
    8
    sage: [f in I for f in I.gens()]
    [True, True]

    sage: f = I.gens()[0]
    sage: I.reduce(f)
    0

    sage: g = I.gens()[1]
    sage: I.reduce(g)
    0

    sage: I.reduce(g+x^2)
    x^2

We compute a Groebner basis for cyclic 6, which is a standard
benchmark and test ideal.

    sage: x,y,z,t,u,v = QQ['x,y,z,t,u,v'].gens()
    sage: I = ideal(x + y + z + t + u + v, x*y + y*z + z*t + t*u + u*v + v*x, x*y*z + y*z*t + z*t*u + t*u*v + u*v*x + v*x*y, x*y*z*t + y*z*t*u + z*t*u*v + t*u*v*x + u*v*x*y + v*x*y*z, x*y*z*t*u + y*z*t*u*v + z*t*u*v*x + t*u*v*x*y + u*v*x*y*z + v*x*y*z*t, x*y*z*t*u*v - 1)
    sage: B = I.groebner_basis()
    sage: len(B)
    17
"""

#*****************************************************************************
#
#   SAGE: System for Algebra and Geometry Experimentation
#
#       Copyright (C) 2005 William Stein <wstein@ucsd.edu>
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

from ideal import Ideal_generic
from sage.interfaces.all import singular as singular_default, is_SingularElement
singular = singular_default
from integer import Integer
from sage.structure.sequence import Sequence

class MPolynomialIdeal(Ideal_generic):
    """
    An ideal of a multivariate polynomial ring.
    """
    def __init__(self, ring, gens, coerce=True):
        """
        Create an ideal in a multivariate polynomial ring.

        EXAMPLES:
            sage: R = PolynomialRing(IntegerRing(), 2, ['x','y']); x,y = R.gens()
            sage: R.ideal([x, y])
            Ideal (y, x) of Polynomial Ring in x, y over Integer Ring
            sage: R = PolynomialRing(GF(3), 2); x = R.gens()
            sage: R.ideal([x[0]**2, x[1]**3])
            Ideal (x_1^3, x_0^2) of Polynomial Ring in x_0, x_1 over Finite Field of size 3
        """
        Ideal_generic.__init__(self, ring, gens, coerce=coerce)

class MPolynomialIdeal_singular_repr(MPolynomialIdeal):
    """
    An ideal in a multivariate polynomial ring, which has an underlying Singular
    ring associated to it.

    EXAMPLES:
    """
    def __init__(self, ring, gens, coerce=True):
        MPolynomialIdeal.__init__(self, ring, gens, coerce=coerce)

    def _cmp_(self, other):
        # Groebner basis determine equality since ideals are in the same ring with same term order
        #c = cmp(self.gens(), other.gens())
        #if c == 0:
        #    return c
        return cmp(self.groebner_basis(), other.groebner_basis())

    def _singular_(self, singular=None):
        """
        Return Singular ideal corresponding to this ideal.

        EXAMPLES:
            sage: R, (x,y) = PolynomialRing(Q, 2, 'xy').objgens()
            sage: I = R.ideal([x^3 + y, y])
            sage: S = I._singular_()
            sage: S
            y,
            x^3+y
        """
        if singular is None: singular = singular_default
        try:
            self.ring()._singular_().set_ring()
            I = self.__singular
            if not (I.parent() is singular):
                raise ValueError
            I._check_valid()
            return I
        except (AttributeError, ValueError):
            self.ring()._singular_().set_ring()
            gens = [str(x) for x in self.gens()]
            if len(gens) == 0:
                gens = ['0']
            self.__singular = singular.ideal(gens)
        return self.__singular

    def _contains_(self, f):
        """
        EXAMPLES:
            sage: R, (x,y) = PolynomialRing(Q, 2, 'xy').objgens()
            sage: I = (x^3 + y, y)*R
            sage: x in I
            False
            sage: y in I
            True
            sage: x^3 + 2*y in I
            True
        """
        S = singular_default
        f = S(f)
        I = self._singular_(S).std()
        g = f.reduce(I, 1)  # 1 avoids tail reduction (page 67 of singular book)
        return g.is_zero()

    def complete_primary_decomposition(self, algorithm="sy"):
        r"""
        INPUT:
            algorithm -- string:
                    'sy' -- (default) use the shimoyama-yokoyama algorithm
                    'gtz' -- use the gianni-trager-zacharias algorithm
        OUTPUT:
            list -- a list of primary ideals and their associated
                    primes
                        [(primary ideal, associated prime), ...]

        ALGORITHM: Uses Singular.

        EXAMPLES:
            sage: R, (x,y,z) = PolynomialRing(Q, 3, 'xyz').objgens()
            sage: p = z^2 + 1; q = z^3 + 2
            sage: I = (p*q^2, y-z^2)*R
            sage: pd = I.complete_primary_decomposition(); pd
            [(Ideal (1 + y, 1 + z^2) of Polynomial Ring in x, y, z over Rational Field,
              Ideal (1 + y, 1 + z^2) of Polynomial Ring in x, y, z over Rational Field),
             (Ideal (4 + 4*z^3 + z^6, -1*z^2 + y) of Polynomial Ring in x, y, z over Rational Field,
              Ideal (2 + z^3, -1*z^2 + y) of Polynomial Ring in x, y, z over Rational Field)]
            sage: I.complete_primary_decomposition(algorithm = 'gtz')
            [(Ideal (-1*z^2 + y, 1 + z^2) of Polynomial Ring in x, y, z over Rational Field,
              Ideal (-1*z^2 + y, 1 + z^2) of Polynomial Ring in x, y, z over Rational Field),
             (Ideal (4 + 4*z^3 + z^6, -1*z^2 + y) of Polynomial Ring in x, y, z over Rational Field,
              Ideal (2 + z^3, -1*z^2 + y) of Polynomial Ring in x, y, z over Rational Field)]
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
        """
        EXAMPLES:
            sage: R, (x,y,z) = PolynomialRing(Q, 3, 'xyz').objgens()
            sage: p = z^2 + 1; q = z^3 + 2
            sage: I = (p*q^2, y-z^2)*R
            sage: I.primary_decomposition()
            [Ideal (1 + y, 1 + z^2) of Polynomial Ring in x, y, z over Rational Field,
            Ideal (4 + 4*z^3 + z^6, -1*z^2 + y) of Polynomial Ring in x, y, z over Rational Field]
        """
        return [I for I, _ in self.complete_primary_decomposition(algorithm)]

    def associated_primes(self, algorithm='sy'):
        """
        EXAMPLES:
            sage: R, (x,y,z) = PolynomialRing(Q, 3, 'xyz').objgens()
            sage: p = z^2 + 1; q = z^3 + 2
            sage: I = (p*q^2, y-z^2)*R
            sage: I.associated_primes()
            [Ideal (1 + y, 1 + z^2) of Polynomial Ring in x, y, z over Rational Field,
             Ideal (2 + z^3, -1*z^2 + y) of Polynomial Ring in x, y, z over Rational Field]
        """
        return [P for _,P in self.complete_primary_decomposition(algorithm)]

    def dimension(self):
        """
        The dimension of the ring modulo this ideal.
        """
        try:
            return self.__dimension
        except AttributeError:
            self.__dimension = Integer(self._singular_().std().dim())
        return self.__dimension

    def groebner_basis(self):
        """
        Return a Groebner basis of this ideal.

        ALGORITHM: Uses Singular.

        EXAMPLES:

        We compute a Groebner basis of "cyclic 4" relative to
        lexicographic ordering.

            sage: R = PolynomialRing(RationalField(), 4, ['a','b','c','d'], 'lex')
            sage: a,b,c,d = R.gens()
            sage: I = R.ideal([a+b+c+d, a*b+a*d+b*c+c*d, a*b*c + a*b*d + a*c*d + b*c*d, a*b*c*d-1])
            sage: I
            Ideal (d + c + b + a, -1 + a*b*c*d, c*d + b*c + a*d + a*b, b*c*d + a*c*d + a*b*d + a*b*c) of Polynomial Ring in a, b, c, d over Rational Field
            sage: I.groebner_basis()
             [1 - d^4 - c^2*d^2 + c^2*d^6, -1*d - c + c^2*d^3 + c^3*d^2, -1*d + d^5 - b + b*d^4, -1*d^2 - d^6 + c*d + c^2*d^4 - b*d^5 + b*c, d^2 + 2*b*d + b^2, d + c + b + a]

        \note{Some Groebner basis calculations crash on 64-bit
        opterons with \SAGE's singular build, but work fine with an
        official binary.  If you download and install a Singular
        binary from the Singular website it will not have this problem
        (you can use it with \SAGE by putting it in local/bin/).}
        """
        try:
            return self.__groebner_basis
        except AttributeError:
            S = self._singular_().groebner()
            R = self.ring()
            self.__groebner_basis = Sequence([R(S[i+1]) for i in range(len(S))], R,
                                             check=False, immutable=True)
        return self.__groebner_basis

    def genus(self):
        """
        Return the genus of the projective curve defined by this
        ideal, which must be 1 dimensional.
        """
        try:
            return self.__genus
        except AttributeError:
            I = self._singular_()
            I.parent().lib('normal.lib')
            self.__genus = Integer(I.genus())
            return self.__genus

    def intersection(self, other):
        """
        Return the intersection of the two ideals.

        EXAMPLES:
            sage: R, (x,y) = PolynomialRing(Q, 2, 'xy').objgens()
            sage: I = x*R
            sage: J = y*R
            sage: I.intersection(J)
            Ideal (x*y) of Polynomial Ring in x, y over Rational Field

        The following simple example illustrates that the product need not equal the intersection.
            sage: I = (x^2, y)*R
            sage: J = (y^2, x)*R
            sage: K = I.intersection(J); K
            Ideal (x^2, x*y, y^2) of Polynomial Ring in x, y over Rational Field
            sage: IJ = I*J; IJ
            Ideal (x^2*y^2, x^3, x*y, y^3) of Polynomial Ring in x, y over Rational Field
            sage: IJ == K
            False
        """
        R = self.ring()
        if not isinstance(other, MPolynomialIdeal_singular_repr) or other.ring() != R:
            raise ValueError, "other (=%s) must be an ideal in the ring %s"%(other, R)
        I = self._singular_()
        sing = I.parent()
        J = sing(other)
        K = I.intersect(J)
        return R.ideal(K)


    def minimal_associated_primes(self):
        r"""
        OUTPUT:
            list -- a list of prime ideals

        EXAMPLES:
            sage: R, (x,y,z) = PolynomialRing(Q, 3, 'xyz').objgens()
            sage: p = z^2 + 1; q = z^3 + 2
            sage: I = (p*q^2, y-z^2)*R
            sage: I.minimal_associated_primes ()
            [Ideal (2 + z^3, -1*z^2 + y) of Polynomial Ring in x, y, z over Rational Field,
             Ideal (-1*z^2 + y, 1 + z^2) of Polynomial Ring in x, y, z over Rational Field]

        ALGORITHM: Uses Singular.
        """
        I = self._singular_()
        I.parent().lib('primdec.lib')
        M = I.minAssGTZ()
        R = self.ring()
        return [R.ideal(J) for J in M]

    def radical(self):
        r"""
        The radical of this ideal.

        EXAMPLES:
        This is an obviously not radical ideal:
            sage: R, (x,y,z) = PolynomialRing(QQ, 3, 'xyz').objgens()
            sage: I = (x^2, y^3, (x*z)^4 + y^3 + 10*x^2)*R
            sage: I.radical()
            Ideal (y, x) of Polynomial Ring in x, y, z over Rational Field

        That the radical is correct is clear from the Groebner basis.
            sage: I.groebner_basis()
            [y^3, x^2]

        This is the example from the singular manual:
            sage: p = z^2 + 1; q = z^3 + 2
            sage: I = (p*q^2, y-z^2)*R
            sage: I.radical()
            Ideal (-1*z^2 + y, 2 + 2*z^2 + z^3 + z^5) of Polynomial Ring in x, y, z over Rational Field

        \note{(From Singular manual) A combination of the algorithms
        of Krick/Logar and Kemper is used.  Works also in positive
        characteristic (Kempers algorithm).}

            sage: R,(x,y,z) = PolynomialRing(GF(37), 3, 'xyz').objgens()
            sage: p = z^2 + 1; q = z^3 + 2
            sage: I = (p*q^2, y - z^2)*R
            sage: I.radical()
            Ideal (36*z^2 + y, 2 + 2*z^2 + z^3 + z^5) of Polynomial Ring in x, y, z over Finite Field of size 37
        """
        S = self.ring()
        I = self._singular_()
        I.parent().lib('primdec.lib')
        r = I.radical()
        return S.ideal(r)

    def reduce(self, f):
        """
        Reduce an element modulo a standard basis for this ideal.
        This returns 0 if and only if the element is in this ideal.

        EXAMPLES:
            sage: R, (x,y) = PolynomialRing(Q, 2, 'xy').objgens()
            sage: I = (x^3 + y, y)*R
            sage: I.reduce(y)
            0
            sage: I.reduce(x^3)
            0
            sage: I.reduce(x - y)
            x

            sage: I = (y^2 - (x^3 + x))*R
            sage: I.reduce(x^3)
            y^2 - x
            sage: I.reduce(x^6)
            y^4 - 2*x*y^2 + x^2
            sage: (y^2 - x)^2
            y^4 - 2*x*y^2 + x^2
        """
        try:
            f = self.ring()(f)
            S = singular_default
            I = self._singular_(S)
            g = S(f)
            h = g.reduce(I.std())
            return self.ring()(h)

        except TypeError:

            return f
