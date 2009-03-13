r"""
Isomorphisms of number fields.  Provides isomorphisms between relative and
absolute presentations, to and form vector spaces, name changing maps, etc.

EXAMPLES:
    sage: L.<cuberoot2, zeta3> = CyclotomicField(3).extension(x^3 - 2)
    sage: K = L.absolute_field('a')
    sage: from_K, to_K = K.structure()
    sage: from_K
    Isomorphism map:
      From: Number Field in a with defining polynomial x^6 - 3*x^5 + 6*x^4 - 11*x^3 + 12*x^2 + 3*x + 1
      To:   Number Field in cuberoot2 with defining polynomial x^3 - 2 over its base field
    sage: to_K
    Isomorphism map:
      From: Number Field in cuberoot2 with defining polynomial x^3 - 2 over its base field
      To:   Number Field in a with defining polynomial x^6 - 3*x^5 + 6*x^4 - 11*x^3 + 12*x^2 + 3*x + 1
"""

#*****************************************************************************
#       Copyright (C) 2008 William Stein <wstein@gmail.com>
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

from sage.structure.sage_object import SageObject
from sage.categories.map import Map
from sage.categories.homset import Hom

import sage.rings.rational_field as rational_field

from sage.libs.pari.all import pari

import number_field_element

QQ = rational_field.RationalField()

class NumberFieldIsomorphism(Map):
    r"""
    A base class for various isomorphisms between number fields and
    vector spaces.
    """
    def _repr_type(self):
        return "Isomorphism"

    def is_injective(self):
        return True

    def is_surjective(self):
        return True

class MapVectorSpaceToNumberField(NumberFieldIsomorphism):
    r"""
    EXAMPLES:
        sage: K.<a> = NumberField(x^4 + 3*x + 1)
        sage: V, fr, to = K.vector_space()
        sage: V
        Vector space of dimension 4 over Rational Field
        sage: fr
        Isomorphism map:
          From: Vector space of dimension 4 over Rational Field
          To:   Number Field in a with defining polynomial x^4 + 3*x + 1
        sage: to
        Isomorphism map:
          From: Number Field in a with defining polynomial x^4 + 3*x + 1
          To:   Vector space of dimension 4 over Rational Field
        sage: type(fr), type(to)
        (<class 'sage.rings.number_field.maps.MapVectorSpaceToNumberField'>,
         <class 'sage.rings.number_field.maps.MapNumberFieldToVectorSpace'>)

        sage: fr.is_injective(), fr.is_surjective()
        (True, True)

        sage: fr.domain(), to.codomain()
        (Vector space of dimension 4 over Rational Field, Vector space of dimension 4 over Rational Field)
        sage: to.domain(), fr.codomain()
        (Number Field in a with defining polynomial x^4 + 3*x + 1, Number Field in a with defining polynomial x^4 + 3*x + 1)
        sage: fr * to
        Composite map:
          From: Number Field in a with defining polynomial x^4 + 3*x + 1
          To:   Number Field in a with defining polynomial x^4 + 3*x + 1
          Defn:   Isomorphism map:
                  From: Number Field in a with defining polynomial x^4 + 3*x + 1
                  To:   Vector space of dimension 4 over Rational Field
                then
                  Isomorphism map:
                  From: Vector space of dimension 4 over Rational Field
                  To:   Number Field in a with defining polynomial x^4 + 3*x + 1
        sage: to * fr
        Composite map:
          From: Vector space of dimension 4 over Rational Field
          To:   Vector space of dimension 4 over Rational Field
          Defn:   Isomorphism map:
                  From: Vector space of dimension 4 over Rational Field
                  To:   Number Field in a with defining polynomial x^4 + 3*x + 1
                then
                  Isomorphism map:
                  From: Number Field in a with defining polynomial x^4 + 3*x + 1
                  To:   Vector space of dimension 4 over Rational Field

        sage: to(a), to(a + 1)
        ((0, 1, 0, 0), (1, 1, 0, 0))
        sage: fr(to(a)), fr(V([0, 1, 2, 3]))
        (a, 3*a^3 + 2*a^2 + a)
    """

    def __init__(self, V, K):
        self.__V = V
        self.__K = K
        self.__R = K.polynomial_ring()
        NumberFieldIsomorphism.__init__(self, Hom(V, K))

    def _call_(self, v):
        f = self.__R(self.__V(v).list())
        return self.__K._element_class(self.__K, f)

class MapNumberFieldToVectorSpace(Map):
    def __init__(self, K, V):
        self.__V = V
        self.__K = K
        self.__zero = QQ(0)
        self.__n = K.degree()
        NumberFieldIsomorphism.__init__(self, Hom(K, V))

    def _repr_type(self):
        return "Isomorphism"

    def _call_(self, x):
        y = self.__K(x)
        v = y._coefficients()
        w = v + [self.__zero]*(self.__n - len(v))
        return self.__V(w)

class MapRelativeVectorSpaceToRelativeNumberField(NumberFieldIsomorphism):
    r"""
    EXAMPLES:
        sage: L.<b> = NumberField(x^4 + 3*x^2 + 1)
        sage: K = L.relativize(L.subfields(2)[0][1], 'a'); K
        Number Field in a0 with defining polynomial x^2 - b0*x + 1 over its base field
        sage: V, fr, to = K.relative_vector_space()
        sage: V
        Vector space of dimension 2 over Number Field in b0 with defining polynomial x^2 + 1
        sage: fr
        Isomorphism map:
          From: Vector space of dimension 2 over Number Field in b0 with defining polynomial x^2 + 1
          To:   Number Field in a0 with defining polynomial x^2 - b0*x + 1 over its base field
        sage: to
        Isomorphism map:
          From: Number Field in a0 with defining polynomial x^2 - b0*x + 1 over its base field
          To:   Vector space of dimension 2 over Number Field in b0 with defining polynomial x^2 + 1
        sage: type(fr), type(to)
        (<class 'sage.rings.number_field.maps.MapRelativeVectorSpaceToRelativeNumberField'>,
         <class 'sage.rings.number_field.maps.MapRelativeNumberFieldToRelativeVectorSpace'>)

        sage: a0 = K.gen(); b0 = K.base_field().gen()
        sage: to(a0), to(a0 + 1), to(a0 + b0)
        ((0, 1), (1, 1), (b0, 1))
        sage: fr(to(a0 + 2*b0)), fr(V([0, 1])), fr(V([b0, 2*b0]))
        (a0 + 2*b0, a0, 2*b0*a0 + b0)

        sage: (fr * to)(K.gen()) == K.gen()
        True
        sage: (to * fr)(V([1, 2])) == V([1, 2])
        True
    """
    def __init__(self, V, K):
        self.__V = V
        self.__K = K
        self.__R = K.polynomial_ring()
        self.__rnf = K.pari_rnf()
        self.__B = K.base_field().absolute_field('a')
        NumberFieldIsomorphism.__init__(self, Hom(V, K))

    def _call_(self, v):
        # Given a relative vector space element, we have to
        # compute the corresponding number field element, in terms
        # of an absolute generator.
        w = self.__V(v).list()

        # First, construct from w a PARI polynomial in x with coefficients
        # that are polynomials in y:
        B = self.__B
        _, to_B = B.structure()
        # Apply to_B, so now each coefficient is in an absolute field,
        # and is expressed in terms of a polynomial in y, then make
        # the PARI poly in x.
        w = [pari(to_B(a).polynomial('y')) for a in w]
        h = pari(w).Polrev()

        # Next we write the poly in x over a poly in y in terms
        # of an absolute polynomial for the rnf.
        g = self.__R(self.__rnf.rnfeltreltoabs(h))
        return self.__K._element_class(self.__K, g)

class MapRelativeNumberFieldToRelativeVectorSpace(NumberFieldIsomorphism):
    def __init__(self, K, V):
        self.__V = V
        self.__K = K
        self.__rnf = K.pari_rnf()
        self.__zero = QQ(0)
        self.__n = K.relative_degree()
        self.__x = pari('x')
        self.__y = pari('y')
        self.__B = K.absolute_base_field()
        NumberFieldIsomorphism.__init__(self, Hom(K, V))

    def _call_(self, alpha):
        # An element of a relative number field is represented
        # internally by an absolute polynomial over QQ.
        alpha = self.__K(alpha)
        f = alpha.polynomial('x')
        # f is the absolute polynomial that defines this number field element
        if self.__K.relative_degree() == 1:
            # Pari doesn't return a valid relative polynomial from an
            # absolute one in the case of relative degree one. (Bug
            # submitted to Pari, fixed in svn unstable as of 1/22/08
            # according to Karim Belabas. Trac #1891 is a reminder to
            # remove this workaround once we update Pari in Sage.)
            # However, in this case, the polynomial we want for g is
            # just the polynomial we have for f, but as a polynomial
            # in the generator of the base field of this relative
            # extension. That generator is just -self.__rnf[0][0], so
            # we can plug it in and it works.
            g = f(-1*self.__rnf[0][0])
        else:
            g = self.__rnf.rnfeltabstorel(pari(f))
        # Now g is a relative polynomial that defines this element.
        # This g is a polynomial in a pari variable x with
        # coefficients polynomials in a variable y.
        # These coefficients define the coordinates of the
        # vector we are constructing.

        # The list v below has the coefficients that are the
        # components of the vector we are constructing, but each is
        # converted into polynomials in a variable x, which we will
        # use to define elements of the base field.
        (x, y) = (self.__x, self.__y)
        v = [g.polcoeff(i).subst(x,y) for i in range(self.__n)]
        B,from_B, _ = self.__B
        w = [from_B(B(z)) for z in v]

        # Now w gives the coefficients.
        return self.__V(w)


class IdentityMap(NumberFieldIsomorphism):
    def __init__(self, K):
        self.__K = K
        NumberFieldIsomorphism.__init__(self, Hom(K, K))

    def _repr_type(self):
        return "Identity map"

    def _call_(self, x):
        return self.__K(x)

class NameChangeMap(NumberFieldIsomorphism):
    def __init__(self, K, L):
        self.__K = K
        self.__L = L
        NumberFieldIsomorphism.__init__(self, Hom(K, L))

    def _repr_type(self):
        return "Isomorphism given by variable name change"

    def _call_(self, x):
        y = self.__K(x)
        z = y.__copy__()
        z._set_parent(self.__L)
        return z

class MapRelativeToAbsoluteNumberField(NumberFieldIsomorphism):
    r"""
    EXAMPLES:
        sage: K.<a> = NumberField(x^6 + 4*x^2 + 200)
        sage: L = K.relativize(K.subfields(3)[0][1], 'b'); L
        Number Field in b0 with defining polynomial x^2 + a0 over its base field
        sage: fr, to = L.structure()
        sage: fr
        Relative number field morphism:
          From: Number Field in b0 with defining polynomial x^2 + a0 over its base field
          To:   Number Field in a with defining polynomial x^6 + 4*x^2 + 200
          Defn: b0 |--> a
                a0 |--> -a^2
        sage: to
        Ring morphism:
          From: Number Field in a with defining polynomial x^6 + 4*x^2 + 200
          To:   Number Field in b0 with defining polynomial x^2 + a0 over its base field
          Defn: a |--> b0
        sage: type(fr), type(to)
        (<class 'sage.rings.number_field.morphism.RelativeNumberFieldHomomorphism_from_abs'>,
         <class 'sage.rings.number_field.morphism.NumberFieldHomomorphism_im_gens'>)

        sage: M.<c> = L.absolute_field(); M
        Number Field in c with defining polynomial x^6 + 4*x^2 + 200
        sage: fr, to = M.structure()
        sage: fr
        Isomorphism map:
          From: Number Field in c with defining polynomial x^6 + 4*x^2 + 200
          To:   Number Field in b0 with defining polynomial x^2 + a0 over its base field
        sage: to
        Isomorphism map:
          From: Number Field in b0 with defining polynomial x^2 + a0 over its base field
          To:   Number Field in c with defining polynomial x^6 + 4*x^2 + 200
        sage: type(fr), type(to)
        (<class 'sage.rings.number_field.maps.MapAbsoluteToRelativeNumberField'>,
         <class 'sage.rings.number_field.maps.MapRelativeToAbsoluteNumberField'>)
        sage: fr(M.gen()), to(fr(M.gen())) == M.gen()
        (b0, True)
        sage: to(L.gen()), fr(to(L.gen())) == L.gen()
        (c, True)
        sage: (to * fr)(M.gen()) == M.gen(), (fr * to)(L.gen()) == L.gen()
        (True, True)
    """

    def __init__(self, R, A):
        self.__R = R          # relative field
        self.__A = A          # absolute field
        self.__poly_ring = self.__A.polynomial_ring()
        self.__zero = QQ(0)
        self.__n = A.degree()
        NumberFieldIsomorphism.__init__(self, Hom(R, A))

    def _call_(self, x):
        f = self.__R(x).polynomial()
        return self.__A._element_class(self.__A, f)

class MapAbsoluteToRelativeNumberField(NumberFieldIsomorphism):
    def __init__(self, A, R):
        self.__A = A          # absolute field
        self.__R = R          # relative field
        self.__poly_ring = self.__A.polynomial_ring()
        self.__zero = QQ(0)
        self.__n = A.degree()
        NumberFieldIsomorphism.__init__(self, Hom(A, R))

    def _call_(self, x):
        f = self.__A(x).polynomial()
        return self.__R._element_class(self.__R, f)

class MapVectorSpaceToRelativeNumberField(NumberFieldIsomorphism):
    def __init__(self, V, L, from_V, from_K):
        self.__V = V
        self.__L = L
        self.__from_V = from_V
        self.__from_K = from_K
        NumberFieldIsomorphism.__init__(self, Hom(V, L))

    def _call_(self, x):
        return self.__from_K(self.__from_V(x))

class MapRelativeNumberFieldToVectorSpace(NumberFieldIsomorphism):
    r"""
    EXAMPLES:
        sage: K.<a> = NumberField(x^8 + 100*x^6 + x^2 + 5)
        sage: L = K.relativize(K.subfields(4)[0][1], 'b'); L
        Number Field in b0 with defining polynomial x^2 + a0 over its base field
        sage: L_to_K, K_to_L = L.structure()

        sage: V, fr, to = L.absolute_vector_space()
        sage: V
        Vector space of dimension 8 over Rational Field
        sage: fr
        Isomorphism map:
          From: Vector space of dimension 8 over Rational Field
          To:   Number Field in b0 with defining polynomial x^2 + a0 over its base field
        sage: to
        Isomorphism map:
          From: Number Field in b0 with defining polynomial x^2 + a0 over its base field
          To:   Vector space of dimension 8 over Rational Field
        sage: type(fr), type(to)
        (<class 'sage.rings.number_field.maps.MapVectorSpaceToRelativeNumberField'>,
         <class 'sage.rings.number_field.maps.MapRelativeNumberFieldToVectorSpace'>)

        sage: v = V([1, 1, 1, 1, 0, 1, 1, 1])
        sage: fr(v), to(fr(v)) == v
        ((-a0^3 + a0^2 - a0 + 1)*b0 - a0^3 - a0 + 1, True)
        sage: to(L.gen()), fr(to(L.gen())) == L.gen()
        ((0, 1, 0, 0, 0, 0, 0, 0), True)
    """
    def __init__(self, L, V, to_K, to_V):
        self.__L = L
        self.__V = V
        self.__to_K = to_K
        self.__to_V = to_V
        NumberFieldIsomorphism.__init__(self, Hom(L, V))

    def _call_(self, x):
        return self.__to_V(self.__to_K(x))

