"""
Ideals in multivariate polynomial rings.

Most functionality of multivariate polynomial ideals in SAGE is
provided through SINGULAR.

AUTHOR:
    -- William Stein
    -- Kiran S. Kedlaya (2006-02-12): added Macaulay2 analogues of
              some Singular features
    -- Martin Albrecht

EXAMPLES:
    sage: x,y,z = QQ['x,y,z'].gens()
    sage: I = ideal(x^5 + y^4 + z^3 - 1,  x^3 + y^3 + z^2 - 1)
    sage: B = I.groebner_basis()
    sage: len(B)
    3
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

    sage: R.<x,y,z,t,u,v> = QQ['x,y,z,t,u,v']
    sage: I = sage.rings.ideal.Cyclic(R,6)
    sage: B = I.groebner_basis()
    sage: len(B)
    45

We compute in a quotient of a polynomial ring over Z/17*Z:
    sage: R.<x,y> = ZZ[]
    sage: S.<a,b> = R.quotient((x^2 + y^2, 17))                 # optional -- requires Macaulay2
    sage: S                                                     # optional
    Quotient of Multivariate Polynomial Ring in x, y over Integer Ring by the ideal (x^2 + y^2, 17)
    sage: a^2 + b^2 == 0                                        # optional
    True
    sage: a^3 - b^2                                             # optional
    -1*a*b^2 - b^2
    sage: (a+b)^17                                              # optional
    a*b^16 + b^17
    sage: S(17) == 0                                            # optional
    True

Working with a polynomial ring over ZZ:
    sage: R.<x,y,z,w> = ZZ['x,y,z,w']
    sage: i = ideal(x^2 + y^2 - z^2 - w^2, x-y)
    sage: j = i^2
    sage: j.groebner_basis()                                    # optional
    [x^2 - 2*x*y + y^2, 2*x*y^2 - 2*y^3 - x*z^2 + y*z^2 - x*w^2 + y*w^2, 4*y^4 - 4*y^2*z^2 + z^4 - 4*y^2*w^2 + 2*z^2*w^2 + w^4]
    sage: y^2 - 2*x*y + x^2 in j                                # optional
    True
    sage: 0 in j                                                # optional
    True

We do a Groebner basis computation over a number field:
    sage: K.<zeta> = CyclotomicField(3)
    sage: R.<x,y,z> = K[]; R
    Multivariate Polynomial Ring in x, y, z over Cyclotomic Field of order 3 and degree 2
    sage: i = ideal(x - zeta*y + 1, x^3 - zeta*y^3); i
    Ideal (x + (-zeta)*y + 1, x^3 + (-zeta)*y^3) of Multivariate Polynomial Ring in x, y, z over Cyclotomic Field of order 3 and degree 2
    sage: i.groebner_basis()
    [x + (-zeta)*y + 1, 3*y^3 + (6*zeta + 3)*y^2 + (3*zeta - 3)*y - zeta - 2]
    sage: S = R.quotient(i); S
    Quotient of Multivariate Polynomial Ring in x, y, z over Cyclotomic Field of order 3 and degree 2 by the ideal (x + (-zeta)*y + 1, x^3 + (-zeta)*y^3)
    sage: S.0  - zeta*S.1
    -1
    sage: S.0^3 - zeta*S.1^3
    0

Two examples from the Mathematica documentation (done in SAGE):
    We compute a Groebner basis:
        sage: R.<x,y> = PolynomialRing(QQ, order='lex')
        sage: ideal(x^2 - 2*y^2, x*y - 3).groebner_basis()
        [2*y^4 - 9, 3*x - 2*y^3]

    We show that three polynomials have no common root:
        sage: R.<x,y> = QQ[]
        sage: ideal(x+y, x^2 - 1, y^2 - 2*x).groebner_basis()
        [1]

TESTS:
    sage: x,y,z = QQ['x,y,z'].gens()
    sage: I = ideal(x^5 + y^4 + z^3 - 1,  x^3 + y^3 + z^2 - 1)
    sage: I == loads(dumps(I))
    True

"""

#*****************************************************************************
#
#   SAGE: System for Algebra and Geometry Experimentation
#
#       Copyright (C) 2005 William Stein <wstein@gmail.com>
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

from sage.rings.ideal import Ideal_generic
from sage.interfaces.all import singular as singular_default, is_SingularElement
from sage.interfaces.all import macaulay2 as macaulay2_default
from sage.interfaces.all import is_SingularElement
singular = singular_default
from sage.rings.integer import Integer
from sage.structure.sequence import Sequence
from sage.misc.sage_eval import sage_eval
from sage.misc.misc import prod
import sage.rings.integer_ring
import sage.rings.polynomial.toy_buchberger as toy_buchberger

def is_MPolynomialIdeal(x):
    return isinstance(x, MPolynomialIdeal)

class MPolynomialIdeal_magma_repr:
    def _magma_(self, magma=None):
        """
        Returns a MAGMA ideal matching self if the base ring coercable to MAGMA
        and MAGMA is available.

        EXAMPLES:
            sage: R.<a,b,c,d,e,f,g,h,i,j> = PolynomialRing(GF(127),10)
            sage: I = sage.rings.ideal.Cyclic(R,4)
            sage: I._magma_() #optional MAGMA
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
        if magma == None:
            import sage.interfaces.magma
            magma = sage.interfaces.magma.magma
        return magma.ideal(self.gens())

    def _magma_groebner_basis(self):
        """
        Computes a Groebner Basis for self using MAGMA if available.

        EXAMPLES:
            sage: R.<a,b,c,d,e,f,g,h,i,j> = PolynomialRing(GF(127),10)
            sage: I = sage.rings.ideal.Cyclic(R,6)
            sage: gb = I.groebner_basis("magma:GroebnerBasis") #optional MAGMA
            sage: len(gb) #optional MAGMA
            45
        """
        try:
            return self.__magma_groebner_basis
        except AttributeError:
            pass
        R = self.ring()
        mgb = self._magma_().GroebnerBasis()
        mgb = [str(mgb[i+1]) for i in range(len(mgb))]
        if R.base_ring().degree() > 1:
            a = str(R.base_ring().gen())
            mgb = [e.replace("$.1",a) for e in mgb]
        B = Sequence([R(e) for e in mgb], R, check=False, immutable=True)
        self.__magma_groebner_basis = B
        return B


class MPolynomialIdeal_singular_repr:
    """
    An ideal in a multivariate polynomial ring, which has an
    underlying Singular ring associated to it.
    """
    def __cmp__(self, other):
        # Groebner basis determine equality since ideals are in the
        # same ring with same term order

        #c = cmp(self.gens(), other.gens())
        #if c == 0:
        #    return c
        l = MPolynomialIdeal(self.ring(), self.groebner_basis()).reduced_basis()
        r = MPolynomialIdeal(self.ring(),other.groebner_basis()).reduced_basis()
        return cmp(r,l)

    def _singular_(self, singular=None):
        """
        Return Singular ideal corresponding to this ideal.

        EXAMPLES:
            sage: R, (x,y) = PolynomialRing(QQ, 2, 'xy').objgens()
            sage: I = R.ideal([x^3 + y, y])
            sage: S = I._singular_()
            sage: S
            x^3+y,
            y
        """
        if singular is None: singular = singular_default
        try:
            self.ring()._singular_(singular).set_ring()
            I = self.__singular
            if not (I.parent() is singular):
                raise ValueError
            I._check_valid()
            return I
        except (AttributeError, ValueError):
            self.ring()._singular_(singular).set_ring()
            gens = [str(x) for x in self.gens()]
            if len(gens) == 0:
                gens = ['0']
            self.__singular = singular.ideal(gens)
        return self.__singular

    def _contains_(self, f):
        """
        Returns True if f is in the ideal self.

        EXAMPLES:
            sage: R, (x,y) = PolynomialRing(QQ, 2, 'xy').objgens()
            sage: I = (x^3 + y, y)*R
            sage: x in I
            False
            sage: y in I
            True
            sage: x^3 + 2*y in I
            True

        NOTE: Requires computation of a Groebner basis, which is a
        very expensive operation.
        """

        if self.base_ring() == sage.rings.integer_ring.ZZ:
            g = self._reduce_using_macaulay2(f)
        else:
            S = singular_default
            f = S(f)
            I = self._singular_(S).groebner()
            g = f.reduce(I, 1)  # 1 avoids tail reduction (page 67 of singular book)
        return g.is_zero()

    def plot(self):
        """
        If you somehow manage to install surf, perhaps you can use
        this function to implicitly plot the real zero locus of this
        ideal (if principal).

        INPUT:
            self -- must be a principal ideal in 2 or 3 vars over QQ.

        EXAMPLES:
        Implicit plotting in 2-d:
            sage: R.<x,y> = PolynomialRing(QQ,2)
            sage: I = R.ideal([y^3 - x^2])
            sage.: I.plot()        # cusp         (optional surf)
            sage: I = R.ideal([y^2 - x^2 - 1])
            sage.: I.plot()        # hyperbola    (optional surf)
            sage: I = R.ideal([y^2 + x^2*(1/4) - 1])
            sage.: I.plot()        # ellipse      (optional surf)
            sage: I = R.ideal([y^2-(x^2-1)*(x-2)])
            sage.: I.plot()        # elliptic curve  (optional surf)

        Implicit plotting in 3-d:
            sage: R.<x,y,z> = PolynomialRing(QQ,3)
            sage: I = R.ideal([y^2 + x^2*(1/4) - z])
            sage.: I.plot()          # a cone         (optional surf)
            sage: I = R.ideal([y^2 + z^2*(1/4) - x])
            sage.: I.plot()          # same code, from a different angle  (optional surf)
            sage: I = R.ideal([x^2*y^2+x^2*z^2+y^2*z^2-16*x*y*z])
            sage.: I.plot()          # Steiner surface   (optional surf)

        AUTHOR:
            -- David Joyner (2006-02-12)
        """
        if self.ring().characteristic() != 0:
            raise TypeError, "base ring must have characteristic 0"
        if not self.is_principal():
            raise TypeError, "self must be principal"
        singular.lib('surf')
        I = singular(self)
        I.plot()

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
            sage: R.<x,y,z> = PolynomialRing(QQ, 3, order='lex')
            sage: p = z^2 + 1; q = z^3 + 2
            sage: I = (p*q^2, y-z^2)*R
            sage: pd = I.complete_primary_decomposition(); pd
            [(Ideal (z^2 + 1, y + 1) of Multivariate Polynomial Ring in x, y, z over Rational Field, Ideal (z^2 + 1, y + 1) of Multivariate Polynomial Ring in x, y, z over Rational Field), (Ideal (z^6 + 4*z^3 + 4, y - z^2) of Multivariate Polynomial Ring in x, y, z over Rational Field, Ideal (z^3 + 2, y - z^2) of Multivariate Polynomial Ring in x, y, z over Rational Field)]

            sage: I.complete_primary_decomposition(algorithm = 'gtz')
            [(Ideal (z^2 + 1, y - z^2) of Multivariate Polynomial Ring in x, y, z over Rational Field, Ideal (z^2 + 1, y - z^2) of Multivariate Polynomial Ring in x, y, z over Rational Field), (Ideal (z^6 + 4*z^3 + 4, y - z^2) of Multivariate Polynomial Ring in x, y, z over Rational Field, Ideal (z^3 + 2, y - z^2) of Multivariate Polynomial Ring in x, y, z over Rational Field)]
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
            sage: R.<x,y,z> = PolynomialRing(QQ, 3, order='lex')
            sage: p = z^2 + 1; q = z^3 + 2
            sage: I = (p*q^2, y-z^2)*R
            sage: I.primary_decomposition()
            [Ideal (z^2 + 1, y + 1) of Multivariate Polynomial Ring in x, y, z over Rational Field, Ideal (z^6 + 4*z^3 + 4, y - z^2) of Multivariate Polynomial Ring in x, y, z over Rational Field]

        """
        return [I for I, _ in self.complete_primary_decomposition(algorithm)]

    def associated_primes(self, algorithm='sy'):
        """
        EXAMPLES:
            sage: R.<x,y,z> = PolynomialRing(QQ, 3)
            sage: p = z^2 + 1; q = z^3 + 2
            sage: I = (p*q^2, y-z^2)*R
            sage: I.associated_primes()
            [Ideal (y + 1, z^2 + 1) of Multivariate Polynomial Ring in x, y, z over Rational Field, Ideal (z^2 - y, y*z + 2, y^2 + 2*z) of Multivariate Polynomial Ring in x, y, z over Rational Field]
        """
        return [P for _,P in self.complete_primary_decomposition(algorithm)]

    def dimension(self):
        """
        The dimension of the ring modulo this ideal.

        NOTE: Requires computation of a Groebner basis, which is a
        very expensive operation.
        """
        try:
            return self.__dimension
        except AttributeError:
            v = list(self.groebner_basis())
            if len(v) == 0:
                v = [0]
            self.__dimension = Integer(singular(v,"ideal").dim())
        return self.__dimension

    def _groebner_basis_using_singular(self, algorithm="groebner"):
        """
        Return a Groebner basis of this ideal. If a groebner basis for
        this ideal has been calculated before the cached Groebner
        basis is returned regardless of the requested algorithm.

        ALGORITHM: Uses Singular.

        INPUT:
            algorithm -- 'groebner' - use Singular's groebner heuristic to choose
                                      an algorithm (default)
                         'std'      - Buchberger's algorithm
                         'stdhilb'  - computes the standard basis of the homogeneous
                                      ideal in the basering, via a Hilbert driven
                                      standard basis computation.
                         'stdfglm'  - computes the standard basis of the ideal in the basering via fglm
                                      (from the degrevlex ordering to the ordering of the basering).
                         'slimgb'   - SlimGB algorithm

        EXAMPLES:

        We compute a Groebner basis of 'cyclic 4' relative to
        lexicographic ordering.

            sage: R.<a,b,c,d> = PolynomialRing(QQ, 4, order='lex')
            sage: I = sage.rings.ideal.Cyclic(R,4); I
            Ideal (a + b + c + d, a*b + a*d + b*c + c*d, a*b*c + a*b*d + a*c*d + b*c*d, a*b*c*d - 1) of Multivariate Polynomial Ring in a, b, c, d over Rational Field
            sage: I._groebner_basis_using_singular()
            [c^2*d^6 - c^2*d^2 - d^4 + 1, c^3*d^2 + c^2*d^3 - c - d, b*d^4 - b + d^5 - d, b*c - b*d^5 + c^2*d^4 + c*d - d^6 - d^2, b^2 + 2*b*d + d^2, a + b + c + d]
        """
        try:
            return self.__groebner_basis
        except AttributeError:
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
            R = self.ring()
            self.__singular_groebner_basis = S #remember this
            self.__groebner_basis = Sequence([R(S[i+1]) for i in range(len(S))], R,
                                             check=False, immutable=True)
        return self.__groebner_basis

    def _groebner_basis_using_libsingular(self, algorithm="std"):
        """
        Return a Groebner basis of this ideal. If a groebner basis for
        this ideal has been calculated before the cached Groebner
        basis is returned regardless of the requested algorithm.

        ALGORITHM: Uses libSINGULAR.

        INPUT:
            algorithm -- 'std'      - Buchberger's algorithm
                         'slimgb'   - SlimGB algorithm

        EXAMPLES:

        We compute a Groebner basis of 'cyclic 4' relative to
        lexicographic ordering.

            sage: R.<a,b,c,d> = PolynomialRing(QQ, 4, order='lex')
            sage: I = sage.rings.ideal.Cyclic(R,4); I
            Ideal (a + b + c + d, a*b + a*d + b*c + c*d, a*b*c + a*b*d + a*c*d + b*c*d, a*b*c*d - 1) of Multivariate Polynomial Ring in a, b, c, d over Rational Field
            sage: I._groebner_basis_using_libsingular()
            [c^2*d^6 - c^2*d^2 - d^4 + 1, c^3*d^2 + c^2*d^3 - c - d, b*d^4 - b + d^5 - d, b*c - b*d^5 + c^2*d^4 + c*d - d^6 - d^2, b^2 + 2*b*d + d^2, a + b + c + d]
        """
        from sage.rings.polynomial.multi_polynomial_ideal_libsingular import std_libsingular, slimgb_libsingular

        try:
            return self.__groebner_basis
        except AttributeError:
            if algorithm=="std":
                S = std_libsingular(self)
            elif algorithm=="slimgb":
                S = slimgb_libsingular(self)
            else:
                raise TypeError, "algorithm '%s' unknown"%algorithm
            self.__groebner_basis = S
        return self.__groebner_basis

    def _singular_groebner_basis(self):
        try:
            S = self.__singular_groebner_basis
        except AttributeError:
            G = self.groebner_basis()
            try:
                return self.__singular_groebner_basis
            except AttributeError:
                pass

        try:
            S._check_valid()
            return S
        except ValueError:
            G = self.groebner_basis()
            S = singular(G)
            self.__singular_groebner_basis = S
        return S



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
            sage: R.<x,y> = PolynomialRing(QQ, 2, order='lex')
            sage: I = x*R
            sage: J = y*R
            sage: I.intersection(J)
            Ideal (x*y) of Multivariate Polynomial Ring in x, y over Rational Field

        The following simple example illustrates that the product need not equal the intersection.
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


    def minimal_associated_primes(self):
        r"""
        OUTPUT:
            list -- a list of prime ideals

        EXAMPLES:
            sage: R.<x,y,z> = PolynomialRing(QQ, 3, 'xyz')
            sage: p = z^2 + 1; q = z^3 + 2
            sage: I = (p*q^2, y-z^2)*R
            sage: I.minimal_associated_primes ()
            [Ideal (z^3 + 2, -z^2 + y) of Multivariate Polynomial Ring in x, y, z over Rational Field, Ideal (z^2 + 1, -z^2 + y) of Multivariate Polynomial Ring in x, y, z over Rational Field]

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
            sage: R.<x,y,z> = PolynomialRing(QQ, 3)
            sage: I = (x^2, y^3, (x*z)^4 + y^3 + 10*x^2)*R
            sage: I.radical()
            Ideal (y, x) of Multivariate Polynomial Ring in x, y, z over Rational Field

        That the radical is correct is clear from the Groebner basis.
            sage: I.groebner_basis()
            [x^2, y^3]

        This is the example from the singular manual:
            sage: p = z^2 + 1; q = z^3 + 2
            sage: I = (p*q^2, y-z^2)*R
            sage: I.radical()
            Ideal (z^2 - y, y^2*z + y*z + 2*y + 2) of Multivariate Polynomial Ring in x, y, z over Rational Field

        \note{(From Singular manual) A combination of the algorithms
        of Krick/Logar and Kemper is used.  Works also in positive
        characteristic (Kempers algorithm).}

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

    def integral_closure(self, p=0, r=True):
        r"""
        Let I == self.

        Returns the integral closure of $I, \ldots, I^p$, where $sI$
        is an ideal in the polynomial ring $R=k[x(1),...x(n)]$. If $p$ is
        not given, or $p=0$, compute the closure of all powers up to
        the maximum degree in t occurring in the closure of $R[It]$ (so
        this is the last power whose closure is not just the
        sum/product of the smaller). If $r$ is given and \code{r is True},
        \code{I.integral_closure()} starts with a check whether I is already a
        radical ideal.

        INPUT:
            p -- powers of I (default: 0)
            r -- check whether self is a radical ideal first (default: True)

        EXAMPLE:
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

    def reduce(self, f):
        """
        Reduce an element modulo a standard basis for this ideal.
        This returns 0 if and only if the element is in this ideal.

        EXAMPLES:
            sage: R.<x,y> = PolynomialRing(QQ, 2)
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

        NOTE: Requires computation of a Groebner basis, which is a
        very expensive operation.
        """
        if self.base_ring() == sage.rings.integer_ring.ZZ:
            return self._reduce_using_macaulay2(f)

        try:
            singular = self._singular_groebner_basis().parent()
        except (AttributeError, RuntimeError):
            singular = self._singular_groebner_basis().parent()

        f = self.ring()(f)
        g = singular(f)
        try:
            self.ring()._singular_(singular).set_ring()
            h = g.reduce(self._singular_groebner_basis())
        except TypeError:
            # This is OK, since f is in the right ring -- type error
            # just means it's a rational
            return f
        try:
            return self.ring()(h)
        except (TypeError, RuntimeError):
            return self.ring()(h[1])
            # Why the above?
            # For mysterious reasons, sometimes Singular returns a length
            # one vector with the reduced polynomial in it.
            # This occurs in the following example:
            #sage: R.<x,y> = PolynomialRing(QQ, 2)
            #sage: S.<a,b> = R.quotient(x^2 + y^2)
            #sage: phi = S.hom([b,a])
            #sage: loads(dumps(phi))


    def syzygy_module(self):
        r"""
        Computes the first syzygy (i.e., the module of relations of
        the given generators) of the ideal.

        ALGORITHM: Uses Singular's syz command

        \note{The syz module is transposed before being returned}
        """
        return self._singular_().syz().transpose().sage_matrix(self.ring())

    def reduced_basis(self):
        r"""
        returns $(g_1, \dots, g_s)$ such that:

        * $(f_1,\dots,f_n) = (g_1,\dots,g_s)$
        * $LT(g_i)\neq LT(g_j)$ for all $i\neq j$
        * $LT(g_i)$ does not divide m for all monomials m of
          $\{g_1,\dots,g_{i-1},g_{i+1},\dots,g_s\}$
        * $LC(g_i) == 1$ for all $i$.

        ALGORITHM: Uses Singular's interred command

        \note{G. Pfister recommended setting option(redSB) before
        using interred for this purpose. Though the manual doesn't
        mention it.}
        """
        from sage.rings.polynomial.multi_polynomial_ideal_libsingular import interred_libsingular
        from sage.rings.polynomial.multi_polynomial_libsingular import MPolynomialRing_libsingular

        R = self.ring()

        if isinstance(R,MPolynomialRing_libsingular):
            return interred_libsingular(self)
        else:
            s = self._singular_().parent()
            o = s.option("get")
            s.option("redSB")
            s.option("redTail")
            ret = []
            for f in self._singular_().interred():
                f = R(f)
                ret.append(f/f.lc()) # lead coeffs are not reduced by interred
            ret = Sequence( ret, R, check=False, immutable=True)
            s.option("set",o)

        return ret

    def basis_is_groebner(self):
        """
        Returns true if self.gens() form a Groebner Basis. This is done by
        trying to lift Syz(LM(self)) to Syz(self) as self is a Groebner
        Basis if and only if for every element S in Syz(LM(self)):
        $$S \cdot G = \sum_{i=0}^{m} h_ig_i \rightarrow_G 0.$$.

        ALGORITHM: Uses Singular

        EXAMPLE:
            sage: R.<a,b,c,d,e,f,g,h,i,j> = PolynomialRing(GF(127),10)
            sage: I = sage.rings.ideal.Cyclic(R,4)
            sage: I.basis_is_groebner()
            False
            sage: I2 = Ideal(I.groebner_basis())
            sage: I2.basis_is_groebner()
            True

        \note{From the Singular Manualf for the reduce function we use in
        this method: 'The result may have no meaning if the second
        argument (self, malb) is not a standard basis'. I (malb) believe
        this refers to the mathematical fact that the results may have no
        meaning if self is no standard basis, i.e., Singular doesn't 'add'
        any additional 'nonsense' to the result. So we may acutally use
        reduce to determine if self is a Groebner Basis.}
        """
        from sage.matrix.constructor import matrix
        singular = self._singular_().parent()
        R = self.ring()

        F = singular( self.gens(), "module" )
        LTF = singular( [f.lt() for f in self.gens()] , "module" )

        M = (F * LTF.syz()).reduce(self._singular_())

        for i in range(M.nrows()):
            if int(singular.eval("%s[1][%s+1]!=0"%(M.name(),i))):
                return False
        return True

    def transformed_basis(self,algorithm="gwalk", other_ring=None):
        """
        Returns a lex or other_ring Groebner Basis for a given ideal
        self which must be represented through a Groebner Basis.

        INPUT:
           algorithm -- Options are:
                        * fglm - FGLM algorithm. The input ideal must be
                                 a reduced Groebner Basis of a zero-dimensional ideal
                        * gwalk (default) - Groebner Walk algorithm
                        * awalk1 - 'first alternative' algorithm
                        * awalk2 - 'second alternative' algorithm
                        * twalk  - Tran algorithm
                        * fwalk  - Fractal Walk algorithm
           other_ring  -- only valid for algorithm 'fglm', if provided conversion will
                          be performed to this ring. Otherwise a lex Groebner basis will
                          be returned.
        EXAMPLES:
           sage: # example from the Singular manual page of fglm
           sage: R.<x,y,z> = PolynomialRing(QQ,3)
           sage: I = Ideal([y^3+x^2,x^2*y+x^2, x^3-x^2, z^4-x^2-y])
           sage: singular.option('redSB')
           sage: I = Ideal(I.groebner_basis())
           sage: singular.option('noredSB') #reset
           sage: S.<z,x,y> = PolynomialRing(QQ,3,order='lex')
           sage: J = Ideal(I.transformed_basis('fglm',S))
           sage: J
           Ideal (y^4 + y^3, x*y^3 - y^3, x^2 + y^3, z^4 + y^3 - y) of Multivariate Polynomial Ring in z, x, y over Rational Field
           sage: # example from the Singular manual page of gwalk
           sage: R.<z,y,x>=PolynomialRing(GF(32003),3,order='lex')
           sage: I=Ideal([y^3+x*y*z+y^2*z+x*z^3,3+x*y+x^2*y+y^2*z])
           sage: I.transformed_basis('gwalk')
           [y^9 - y^7*x^2 - y^7*x - y^6*x^3 - y^6*x^2 - 3*y^6 - 3*y^5*x - y^3*x^7 - 3*y^3*x^6 - 3*y^3*x^5 - y^3*x^4 - 9*y^2*x^5 - 18*y^2*x^4 - 9*y^2*x^3 - 27*y*x^3 - 27*y*x^2 - 27*x,
           z*x + 8297*y^8*x^2 + 8297*y^8*x + 3556*y^7 - 8297*y^6*x^4 + 15409*y^6*x^3 - 8297*y^6*x^2 - 8297*y^5*x^5 + 15409*y^5*x^4 - 8297*y^5*x^3 + 3556*y^5*x^2 + 3556*y^5*x + 3556*y^4*x^3 + 3556*y^4*x^2 - 10668*y^4 - 10668*y^3*x - 8297*y^2*x^9 - 1185*y^2*x^8 + 14224*y^2*x^7 - 1185*y^2*x^6 - 8297*y^2*x^5 - 14223*y*x^7 - 10666*y*x^6 - 10666*y*x^5 - 14223*y*x^4 + x^5 + 2*x^4 + x^3,
           z*y^2 + y*x^2 + y*x + 3]


        ALGORITHM: Uses Singular
        """
        from sage.rings.polynomial.multi_polynomial_ring import TermOrder,MPolynomialRing
        from sage.rings.quotient_ring import is_QuotientRing

        Is = self._singular_()
        R = self.ring()

        if algorithm in ("gwalk","awalk1","awalk2","twalk","fwalk"):
            singular.LIB("grwalk")
            gb = singular("%s(%s)"%(algorithm,Is.name()))
            return [R(f) for f in gb]
        elif algorithm == "fglm":
            Rs = self.ring()._singular_()

            # new ring
            if other_ring==None:
                nR = MPolynomialRing(R.base_ring(),R.ngens(), names=R.variable_names(), order="lex")
            else:
                nR = other_ring
            nR._singular_().set_ring()

            nIs = singular.fglm(Rs,Is)

            return [nR(f) for f in nIs]

        else:
            raise TypeError, "Cannot convert basis with given algorithm"

    def elimination_ideal(self, variables):
        """
        Returns the elimination ideal of self with respect to the
        variables given in 'variables'.

        INPUT:
            variables -- a list or tuple of variables in self.ring()

        EXAMPLE:
            sage: R.<x,y,t,s,z> = PolynomialRing(QQ,5)
            sage: I = R * [x-t,y-t^2,z-t^3,s-x+y^3]
            sage: I.elimination_ideal([t,s])
            Ideal (y^2 - x*z, x*y - z, x^2 - y) of Multivariate Polynomial Ring in x, y, t, s, z over Rational Field

        ALGORITHM: Uses SINGULAR

        NOTE: Requires computation of a Groebner basis, which is a
        very expensive operation.
        """
        if not isinstance(variables, (list,tuple)):
            variables = (variables,)

        try:
            Is = self.__singular_groebner_basis
        except AttributeError:
            Is = self._singular_()

        R = self.ring()
        return MPolynomialIdeal(R, [f.sage_poly(R) for f in Is.eliminate( prod(variables) ) ] )

    def quotient(self, J):
        r"""
        Given ideals I (==self) and J of the same polynomial ring P,
        return the ideal quotient of I by J consisting of the
        polynomials a of P such that $\{aJ \subset I\}$.

        This is also referred to as the colon ideal (I:J).

        INPUT:
            J -- multivariate polynomial ideal

        EXAMPLE:
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

class MPolynomialIdeal_macaulay2_repr:
    """
    An ideal in a multivariate polynomial ring, which has an underlying
    Macaulay2 ring associated to it.

    EXAMPLES:
        sage: R.<x,y,z,w> = PolynomialRing(ZZ, 4) # optional
        sage: I = ideal(x*y-z^2, y^2-w^2)       # optional
        sage: I                                 # optional
        Ideal (x*y - z^2, y^2 - w^2) of Multivariate Polynomial Ring in x, y, z, w over Integer Ring
    """
    #def __init__(self, ring, gens, coerce=True):
    #    MPolynomialIdeal.__init__(self, ring, gens, coerce=coerce)

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

    def _macaulay2_groebner_basis(self):
        r"""
        Return the Groebner basis for this ideal, computed using Macaulay2.

        ALGORITHM: Computed using Macaulay2.  A big advantage of
        Macaulay2 is that it can compute Groebner basis of ideals in
        polynomial rings over the integers.

        EXAMPLE:
            sage: R.<x,y,z,w> = PolynomialRing(ZZ, 4)
            sage: I = ideal(x*y-z^2, y^2-w^2)
            sage: I.groebner_basis()                                     # optional -- requires macaulay2
            [y^2 - w^2, x*y - z^2, y*z^2 - x*w^2, z^4 - x^2*w^2]

        Groebner basis can be used to compute in $\Z/n\Z[x,\ldots]$.

            sage: R.<x,y,z> = ZZ[]
            sage: I = ideal([y^2*z - x^3 - 19*x*z, y^2, 19^2])
            sage: I.groebner_basis()                                     # optional -- requires macaulay2
            [361, y^2, x^3 + 19*x*z]
            sage: I = ideal([y^2*z - x^3 - 19^2*x*z, y^2, 19^2])
            sage: I.groebner_basis()                                     # optional -- requires macaulay2
            [361, y^2, x^3]
        """
        try:
            return self.__groebner_basis
        except AttributeError:
            I = self._macaulay2_()
            G = str(I.gb().generators().str()).replace('\n','')
            i = G.rfind('{{')
            j = G.rfind('}}')
            G = G[i+2:j].split(',')
            L = self.ring().gens_dict()
            B = [sage_eval(f, L) for f in G]
            B = Sequence(B, self.ring(), check=False)
            B.sort()
            B.set_immutable()
            self.__groebner_basis = B
            return B

    def _reduce_using_macaulay2(self, f):
        """
        EXAMPLES:
            sage: R.<x,y,z,w> = PolynomialRing(ZZ, 4)
            sage: I = ideal(x*y-z^2, y^2-w^2)
            sage: I._reduce_using_macaulay2(x*y-z^2 + y^2)    # optional
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
    """
    An ideal of a multivariate polynomial ring.
    """
    def __init__(self, ring, gens, coerce=True):
        """
        Create an ideal in a multivariate polynomial ring.

        EXAMPLES:
            sage: R.<x,y> = PolynomialRing(IntegerRing(), 2, order='lex')
            sage: R.ideal([x, y])
            Ideal (x, y) of Multivariate Polynomial Ring in x, y over Integer Ring
            sage: R.<x0,x1> = GF(3)[]
            sage: R.ideal([x0^2, x1^3])
            Ideal (x0^2, x1^3) of Multivariate Polynomial Ring in x0, x1 over Finite Field of size 3
        """
        Ideal_generic.__init__(self, ring, gens, coerce=coerce)

    def groebner_fan(self, is_groebner_basis=False, symmetry=None, verbose=False):
        r"""
        Return the Groebner fan of this ideal.

        The base ring must be $\Q$ or a finite field $\F_p$ of with
        $p \leq 32749$.

        INPUT:
            is_groebner_basis -- bool (default False).  if True, then I.gens() must be
                                 a Groebner basis with respect to the standard
                                 degree lexicographic term order.
            symmetry -- default: None; if not None, describes symmetries of the ideal
            verbose -- default: False; if True, printout useful info during computations
        """
        import groebner_fan
        return groebner_fan.GroebnerFan(self, is_groebner_basis=is_groebner_basis,
                                        symmetry=symmetry, verbose=verbose)

    def groebner_basis(self, algorithm=None):
        """
        Return a Groebner basis of this ideal.

        INPUT:
            algorithm -- determines the algorithm to use, available are:
                         * None - autoselect (default)
                         * 'singular:groebner' - Singular's groebner command
                         * 'singular:std' - Singular's std command
                         * 'singular:stdhilb' - Singular's stdhib command
                         * 'singular:stdfglm' - Singular's stdfglm command
                         * 'singular:slimgb' - Singular's slimgb command
                         * 'libsingular:std' - libSINGULAR's std command
                         * 'libsingular:slimgb' - libSINGULAR's slimgb command
                         * 'toy:buchberger' - SAGE's toy/educational buchberger without strategy
                         * 'toy:buchberger2' - SAGE's toy/educational buchberger with strategy
                         * 'macaulay2:gb' (if available) - Macaulay2's gb command
                         * 'magma:GroebnerBasis' (if available) - MAGMA's Groebnerbasis command

        EXAMPLES:
            Consider Katsura-3 over QQ with term ordering 'degrevlex'

            sage: P.<a,b,c> = PolynomialRing(QQ,3, order='lex')
            sage: I = sage.rings.ideal.Katsura(P,3) # regenerate to prevent caching
            sage: I.groebner_basis()
            [84*c^4 - 40*c^3 + c^2 + c, 7*b + 210*c^3 - 79*c^2 + 3*c, 7*a - 420*c^3 + 158*c^2 + 8*c - 7]


            sage: I = sage.rings.ideal.Katsura(P,3) # regenerate to prevent caching
            sage: I.groebner_basis('singular:groebner')
            [84*c^4 - 40*c^3 + c^2 + c, 7*b + 210*c^3 - 79*c^2 + 3*c, 7*a - 420*c^3 + 158*c^2 + 8*c - 7]

            sage: I = sage.rings.ideal.Katsura(P,3) # regenerate to prevent caching
            sage: I.groebner_basis('singular:std')
            [84*c^4 - 40*c^3 + c^2 + c, 7*b + 210*c^3 - 79*c^2 + 3*c, 7*a - 420*c^3 + 158*c^2 + 8*c - 7]

            sage: I = sage.rings.ideal.Katsura(P,3) # regenerate to prevent caching
            sage: I.groebner_basis('singular:stdhilb')
            [84*c^4 - 40*c^3 + c^2 + c, 7*b + 210*c^3 - 79*c^2 + 3*c, 7*a - 420*c^3 + 158*c^2 + 8*c - 7]

            sage: I = sage.rings.ideal.Katsura(P,3) # regenerate to prevent caching
            sage: I.groebner_basis('singular:stdfglm')
            [84*c^4 - 40*c^3 + c^2 + c, 7*b + 210*c^3 - 79*c^2 + 3*c, 7*a - 420*c^3 + 158*c^2 + 8*c - 7]

            sage: I = sage.rings.ideal.Katsura(P,3) # regenerate to prevent caching
            sage: I.groebner_basis('singular:slimgb')
            [84*c^4 - 40*c^3 + c^2 + c, 7*b + 210*c^3 - 79*c^2 + 3*c, 7*a - 420*c^3 + 158*c^2 + 8*c - 7]

            sage: I = sage.rings.ideal.Katsura(P,3) # regenerate to prevent caching
            sage: I.groebner_basis('toy:buchberger')
            [-30*c^4 + 100/7*c^3 - 5/14*c^2 - 5/14*c, -2*b^2 - b*c + 1/2*b, -7/125*b - 42/25*c^3 + 79/125*c^2 - 3/125*c,
            a + 2*b + 2*c - 1, a^2 - a + 2*b^2 + 2*c^2, -5*b*c + 1/2*b - 6*c^2 + 2*c, 2*a*b + 2*b*c - b]

            sage: I = sage.rings.ideal.Katsura(P,3) # regenerate to prevent caching
            sage: I.groebner_basis('toy:buchberger2')
            [30*c^4 - 100/7*c^3 + 5/14*c^2 + 5/14*c, a + 2*b + 2*c - 1,
             7/125*b + 42/25*c^3 - 79/125*c^2 + 3/125*c, a^2 - a + 2*b^2 + 2*c^2]

            sage: I = sage.rings.ideal.Katsura(P,3) # regenerate to prevent caching
            sage: I.groebner_basis('macaulay2:gb') # optional requires Macaulay2
            [84*c^4 - 40*c^3 + c^2 + c, 7*b + 210*c^3 - 79*c^2 + 3*c, 7*a - 420*c^3 + 158*c^2 + 8*c - 7]

            sage: I = sage.rings.ideal.Katsura(P,3) # regenerate to prevent caching
            sage: I.groebner_basis('magma:GroebnerBasis') # optional requires MAGMA
            [a - 60*c^3 + 158/7*c^2 + 8/7*c - 1, b + 30*c^3 - 79/7*c^2 + 3/7*c, c^4 - 10/21*c^3 + 1/84*c^2 + 1/84*c]

            If Macaulay2 is installed, Groebner bases over ZZ can be computed.

            sage: P.<a,b,c> = PolynomialRing(ZZ,3)
            sage: I = P * (a + 2*b + 2*c - 1, a^2 - a + 2*b^2 + 2*c^2, 2*a*b + 2*b*c - b)
            sage: I.groebner_basis() #optional requires Macaulay2
            [a + 2*b + 2*c - 1, 10*b*c + 12*c^2 - b - 4*c, 2*b^2 - 4*b*c - 6*c^2 + 2*c,
             42*c^3 + b^2 + 2*b*c - 14*c^2 + b, 2*b*c^2 - 6*c^3 + b^2 + 5*b*c + 8*c^2 - b - 2*c,
             b^3 + b*c^2 + 12*c^3 + b^2 + b*c - 4*c^2]

            SAGE also supports local orderings (via SINGULAR):

            sage: P.<x,y,z> = PolynomialRing(QQ,3,order='negdegrevlex')
            sage: I = P * (  x*y*z + z^5, 2*x^2 + y^3 + z^7, 3*z^5 +y ^5 )
            sage: I.groebner_basis()
            [2*x^2 + y^3, x*y*z + z^5, y^5 + 3*z^5, y^4*z - 2*x*z^5, z^6]

        ALGORITHM: Uses Singular, MAGMA (if available), Macaulay2 (if
        available), or toy implementation.

        """
        if algorithm is None:
            if self.ring().base_ring() == sage.rings.integer_ring.ZZ:
                return self._macaulay2_groebner_basis()
            else:
                return self._groebner_basis_using_singular("groebner")
        elif algorithm.startswith('singular:'):
            return self._groebner_basis_using_singular(algorithm[9:])
        elif algorithm.startswith('libsingular:'):
            return self._groebner_basis_using_libsingular(algorithm[len('libsingular:'):])
        elif algorithm == 'macaulay2:gb':
            return self._macaulay2_groebner_basis()
        elif algorithm == 'magma:GroebnerBasis':
            return self._magma_groebner_basis()
        elif algorithm == 'toy:buchberger':
            return toy_buchberger.buchberger(self)
        elif algorithm == 'toy:buchberger2':
            return toy_buchberger.buchberger_improved(self)
        else:
            raise TypeError, "algorithm '%s' unknown"%algorithm


    def change_ring(self, P):
        r"""
        Return the ideal $I$ in $P$ spanned by the generators $g_1,
        ..., g_n$ of self as returned by self.gens().

        INPUT:
            P -- a multivariate polynomial ring

        EXAMPLE:
           sage: P.<x,y,z> = PolynomialRing(QQ,3,order='lex')
           sage: I = sage.rings.ideal.Cyclic(P)
           sage: I
           Ideal (x + y + z, x*y + x*z + y*z, x*y*z - 1) of Multivariate Polynomial Ring in x, y, z over Rational Field
           sage: I.groebner_basis()
           [z^3 - 1, y^2 + y*z + z^2, x + y + z]

           sage: Q.<x,y,z> = P.new_ring(order='degrevlex'); Q
           Multivariate Polynomial Ring in x, y, z over Rational Field
           sage: Q.term_order()
           Degree reverse lexicographic term order

           sage: J = I.change_ring(Q); J
           Ideal (x + y + z, x*y + x*z + y*z, x*y*z - 1) of Multivariate Polynomial Ring in x, y, z over Rational Field
           sage: J.groebner_basis()
           [x + y + z, y^2 + y*z + z^2, z^3 - 1]
        """
        return P.ideal([P(f) for f in self.gens()])
