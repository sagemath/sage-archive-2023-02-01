"""
Sets of homomorphisms between number fields
"""

#*****************************************************************************
#       Copyright (C) 2007 William Stein <wstein@gmail.com>
#       Copyright (C) 2020 Peter Bruin <P.J.Bruin@math.leidenuniv.nl>
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 2 of the License, or
# (at your option) any later version.
#                  http://www.gnu.org/licenses/
#*****************************************************************************

from sage.misc.cachefunc import cached_method
from sage.misc.superseded import deprecation

from sage.rings.homset import RingHomset_generic
from sage.rings.number_field.morphism import (NumberFieldHomomorphism_im_gens,
                                              RelativeNumberFieldHomomorphism_from_abs,
                                              CyclotomicFieldHomomorphism_im_gens)
from sage.rings.integer import Integer
from sage.rings.finite_rings.integer_mod_ring import Zmod
from sage.structure.sequence import Sequence


class NumberFieldHomset(RingHomset_generic):
    """
    Set of homomorphisms with domain a given number field.

    TESTS::

        sage: H = Hom(QuadraticField(-1, 'a'), QuadraticField(-1, 'b'))
        sage: TestSuite(H).run()
          Failure in _test_category:
        ...
        The following tests failed: _test_elements
    """
    def __init__(self, R, S, category=None):
        """
        TESTS:

        Check that :trac:`23647` is fixed::

            sage: K.<a, b> = NumberField([x^2 - 2, x^2 - 3])
            sage: e, u, v, w = End(K)
            sage: e.abs_hom().parent().category()
            Category of homsets of number fields
            sage: (v*v).abs_hom().parent().category()
            Category of homsets of number fields
        """
        if category is None:
            from sage.categories.all import Fields, NumberFields
            if S in NumberFields:
                category = NumberFields()
            elif S in Fields:
                category = Fields()
        RingHomset_generic.__init__(self, R, S, category)

    def __call__(self, im_gens, check=True):
        """
        Create the homomorphism sending the generators to ``im_gens``.

        EXAMPLES::

            sage: H = Hom(QuadraticField(-1, 'a'), QuadraticField(-1, 'b'))
            sage: phi = H([H.domain().gen()]); phi # indirect doctest
            Ring morphism:
              From: Number Field in a with defining polynomial x^2 + 1 with a = 1*I
              To:   Number Field in b with defining polynomial x^2 + 1 with b = 1*I
              Defn: a |--> b
        """
        if isinstance(im_gens, NumberFieldHomomorphism_im_gens):
            return self._coerce_impl(im_gens)
        try:
            return NumberFieldHomomorphism_im_gens(self, im_gens, check=check)
        except (NotImplementedError, ValueError):
            try:
                return self._coerce_impl(im_gens)
            except TypeError:
                raise TypeError("images do not define a valid homomorphism")

    def _coerce_impl(self, x):
        r"""
        Canonical coercion of ``x`` into this homset. The only things that
        coerce canonically into self are elements of self and of homsets equal
        to self.

        EXAMPLES::

            sage: H1 = End(QuadraticField(-1, 'a'))
            sage: H1.coerce(loads(dumps(H1[1]))) # indirect doctest
            Ring endomorphism of Number Field in a with defining polynomial x^2 + 1 with a = 1*I
              Defn: a |--> -a

        TESTS:

        We can move morphisms between categories::

            sage: f = H1.an_element()
            sage: g = End(H1.domain(), category=Rings())(f)
            sage: f == End(H1.domain(), category=NumberFields())(g)
            True
        """
        if not isinstance(x, NumberFieldHomomorphism_im_gens):
            raise TypeError
        if x.parent() is self:
            return x
        from sage.categories.all import NumberFields, Rings
        if (x.parent() == self or
            (x.domain() == self.domain() and x.codomain() == self.codomain() and
             # This would be the better check, however it returns False currently:
             # self.homset_category().is_full_subcategory(x.category_for())
             # So we check instead that this is a morphism anywhere between
             # Rings and NumberFields where the hom spaces do not change.
             NumberFields().is_subcategory(self.homset_category()) and
             self.homset_category().is_subcategory(Rings()) and
             NumberFields().is_subcategory(x.category_for()) and
             x.category_for().is_subcategory(Rings()))):
            return NumberFieldHomomorphism_im_gens(self, x.im_gens(), check=False)
        raise TypeError

    def _an_element_(self):
        r"""
        Return an element of this set of embeddings.

        EXAMPLES::

            sage: H = Hom(QuadraticField(-1, 'a'), QuadraticField(-1, 'b'))
            sage: H.an_element() # indirect doctest
            Ring morphism:
              From: Number Field in a with defining polynomial x^2 + 1 with a = 1*I
              To:   Number Field in b with defining polynomial x^2 + 1 with b = 1*I
              Defn: a |--> b

            sage: H = Hom(QuadraticField(-1, 'a'), QuadraticField(-2, 'b'))
            sage: H.an_element()
            Traceback (most recent call last):
            ...
            EmptySetError: There is no morphism from Number Field in a with defining polynomial x^2 + 1 with a = 1*I to Number Field in b with defining polynomial x^2 + 2 with b = 1.414213562373095?*I
        """
        L = self.list()
        if len(L) != 0:
            return L[0]
        else:
            from sage.categories.sets_cat import EmptySetError
            raise EmptySetError("There is no morphism from {} to {}".format(
                                              self.domain(), self.codomain()))

    def _repr_(self):
        r"""
        String representation of this homset.

        EXAMPLES::

            sage: repr(Hom(QuadraticField(-1, 'a'), QuadraticField(-1, 'b'))) # indirect doctest
            'Set of field embeddings from Number Field in a with defining polynomial x^2 + 1 with a = 1*I to Number Field in b with defining polynomial x^2 + 1 with b = 1*I'
            sage: repr(Hom(QuadraticField(-1, 'a'), QuadraticField(-1, 'a'))) # indirect doctest
            'Automorphism group of Number Field in a with defining polynomial x^2 + 1 with a = 1*I'
        """
        D = self.domain()
        C = self.codomain()
        if C == D:
            return "Automorphism group of {}".format(D)
        else:
            return "Set of field embeddings from {} to {}".format(D, C)

    def order(self):
        """
        Return the order of this set of field homomorphism.

        EXAMPLES::

            sage: k.<a> = NumberField(x^2 + 1)
            sage: End(k)
            Automorphism group of Number Field in a with defining polynomial x^2 + 1
            sage: End(k).order()
            2
            sage: k.<a> = NumberField(x^3 + 2)
            sage: End(k).order()
            1

            sage: K.<a> = NumberField( [x^3 + 2, x^2 + x + 1] )
            sage: End(K).order()
            6
        """
        return Integer(len(self.list()))

    cardinality = order

    @cached_method
    def list(self):
        """
        Return a list of all the elements of self.

        EXAMPLES::

            sage: K.<a> = NumberField(x^3 - 3*x + 1)
            sage: End(K).list()
            [
            Ring endomorphism of Number Field in a with defining polynomial x^3 - 3*x + 1
              Defn: a |--> a,
            Ring endomorphism of Number Field in a with defining polynomial x^3 - 3*x + 1
              Defn: a |--> a^2 - 2,
            Ring endomorphism of Number Field in a with defining polynomial x^3 - 3*x + 1
              Defn: a |--> -a^2 - a + 2
            ]
            sage: Hom(K, CyclotomicField(9))[0] # indirect doctest
            Ring morphism:
              From: Number Field in a with defining polynomial x^3 - 3*x + 1
              To:   Cyclotomic Field of order 9 and degree 6
              Defn: a |--> -zeta9^4 + zeta9^2 - zeta9

        An example where the codomain is a relative extension::

            sage: K.<a> = NumberField(x^3 - 2)
            sage: L.<b> = K.extension(x^2 + 3)
            sage: Hom(K, L).list()
            [
            Ring morphism:
              From: Number Field in a with defining polynomial x^3 - 2
              To:   Number Field in b with defining polynomial x^2 + 3 over its base field
              Defn: a |--> a,
            Ring morphism:
              From: Number Field in a with defining polynomial x^3 - 2
              To:   Number Field in b with defining polynomial x^2 + 3 over its base field
              Defn: a |--> -1/2*a*b - 1/2*a,
            Ring morphism:
              From: Number Field in a with defining polynomial x^3 - 2
              To:   Number Field in b with defining polynomial x^2 + 3 over its base field
              Defn: a |--> 1/2*a*b - 1/2*a
            ]
        """
        D = self.domain()
        C = self.codomain()
        if D.degree().divides(C.absolute_degree()):
            roots = D.polynomial().roots(ring=C, multiplicities=False)
            v = [D.hom([r], codomain=C, check=False) for r in roots]
        else:
            v = []
        return Sequence(v, universe=self, check=False, immutable=True, cr=v!=[])

    def __getitem__(self, n):
        r"""
        Return the ``n``th element of ``self.list()``.

        EXAMPLES::

            sage: End(CyclotomicField(37))[3] # indirect doctest
            Ring endomorphism of Cyclotomic Field of order 37 and degree 36
              Defn: zeta37 |--> zeta37^4
        """
        return self.list()[n]


class RelativeNumberFieldHomset(NumberFieldHomset):
    """
    Set of homomorphisms with domain a given relative number field.

    EXAMPLES:

    We construct a homomorphism from a relative field by giving
    the image of a generator::

        sage: L.<cuberoot2, zeta3> = CyclotomicField(3).extension(x^3 - 2)
        sage: phi = L.hom([cuberoot2 * zeta3]); phi
        Relative number field endomorphism of Number Field in cuberoot2 with defining polynomial x^3 - 2 over its base field
          Defn: cuberoot2 |--> zeta3*cuberoot2
                zeta3 |--> zeta3
        sage: phi(cuberoot2 + zeta3)
        zeta3*cuberoot2 + zeta3

    In fact, this phi is a generator for the Kummer Galois group of this
    cyclic extension::

        sage: phi(phi(cuberoot2 + zeta3))
        (-zeta3 - 1)*cuberoot2 + zeta3
        sage: phi(phi(phi(cuberoot2 + zeta3)))
        cuberoot2 + zeta3
    """
    def __call__(self, im_gen, base_map=None, base_hom=None, check=True):
        r"""
        Create a homomorphism in this homset from the given data, which can be:

        - A homomorphism from this number field.
        - A homomorphism from the absolute number field corresponding to this
          relative number field.
        - An element (specifying the image of the generator) of a ring into
          which the base ring coerces.
        - A pair consisting of an element of a ring R and a homomorphism from
          the base ring to R.

        EXAMPLES::

            sage: K.<a> = NumberField(x^2 + 1)
            sage: L.<b> = K.extension(x^4 - 2)
            sage: E = End(L)
            sage: E(E[0]) # indirect doctest
            Relative number field endomorphism of Number Field in b with defining polynomial x^4 - 2 over its base field
              Defn: b |--> b
                    a |--> a
            sage: E(L.absolute_field('c').hom(b+a, L)) # indirect doctest
            Relative number field endomorphism of Number Field in b with defining polynomial x^4 - 2 over its base field
              Defn: b |--> b
                    a |--> -a
            sage: E(-b*a) # indirect doctest
            Relative number field endomorphism of Number Field in b with defining polynomial x^4 - 2 over its base field
              Defn: b |--> -a*b
                    a |--> a
            sage: E(-a*b, K.hom([-a])) # indirect doctest
            Relative number field endomorphism of Number Field in b with defining polynomial x^4 - 2 over its base field
              Defn: b |--> -a*b
                    a |--> -a

        You can specify a map on the base field::

            sage: R.<x> = ZZ[]
            sage: K.<i> = NumberField(x^2 + 1)
            sage: L.<b> = K.extension(x^2-17)
            sage: cc = K.hom([-i])
            sage: phi = L.hom([-b],base_map=cc); phi
            Relative number field endomorphism of Number Field in b with defining polynomial x^2 - 17 over its base field
              Defn: b |--> -b
                    i |--> -i

        Using check=False, it is possible to construct homomorphisms into fields such as CC
        where calculations are only approximate.

            sage: K.<a> = QuadraticField(-7)
            sage: f = K.hom([CC(sqrt(-7))], check=False)
            sage: x = polygen(K)
            sage: L.<b> = K.extension(x^2 - a - 5)
            sage: L.Hom(CC)(f(a + 5).sqrt(), f, check=False)
            Relative number field morphism:
              From: Number Field in b with defining polynomial x^2 - a - 5 over its base field
              To:   Complex Field with 53 bits of precision
              Defn: b |--> 2.30833860703888 + 0.573085617291335*I
              a |--> -8.88178419700125e-16 + 2.64575131106459*I
        """
        if base_hom is not None:
            deprecation(26105, "Use base_map rather than base_hom")
            base_map = base_hom
        if isinstance(im_gen, NumberFieldHomomorphism_im_gens):
            # Then it must be a homomorphism from the corresponding
            # absolute number field
            abs_hom = im_gen
            K = abs_hom.domain()
            if K != self.domain().absolute_field(K.variable_name()):
                raise TypeError("domain of morphism must be absolute field of domain.")
            from_K, to_K = K.structure()
            if abs_hom.domain() != K:
                raise ValueError("domain of absolute homomorphism must be absolute field of domain.")
            if abs_hom.codomain() != self.codomain():
                raise ValueError("codomain of absolute homomorphism must be codomain of this homset.")
            return RelativeNumberFieldHomomorphism_from_abs(self, abs_hom)
        if isinstance(im_gen, RelativeNumberFieldHomomorphism_from_abs):
            return self._coerce_impl(im_gen)
        if base_map is None:
            base_map = self.default_base_hom()
        if isinstance(im_gen, (list, tuple)) and len(im_gen) == 1:
            im_gen = im_gen[0]
        if check:
            im_gen = self.codomain()(im_gen)
        return self._from_im(im_gen, base_map=base_map, check=check)

    def _coerce_impl(self, x):
        r"""
        Canonically coerce ``x`` into this homset. This will only work if ``x``
        is already in the homset.

        EXAMPLES::

            sage: L.<a, b> = NumberField([x^3 - x + 1, x^2 + 23])
            sage: E = End(L)
            sage: E.coerce(loads(dumps(E[0])))  # indirect doctest
            Relative number field endomorphism of Number Field in a with defining polynomial x^3 - x + 1 over its base field
              Defn: a |--> a
                    b |--> b
        """
        if not isinstance(x, RelativeNumberFieldHomomorphism_from_abs):
            raise TypeError
        if x.parent() is self:
            return x
        if x.parent() == self:
            return RelativeNumberFieldHomomorphism_from_abs(self, x.abs_hom())
        raise TypeError

    def _from_im(self, im_gen, base_map, check=True):
        """
        Return the homomorphism that acts on the base as given and
        sends the generator of the domain to im_gen.

        EXAMPLES::

            sage: K.<a> = NumberField(x^2 + 23)
            sage: L.<b> = K.extension(x^3 - x + 1)
            sage: End(L)._from_im( -3/23*a*b^2 + (-9/46*a - 1/2)*b + 2/23*a, K.hom([-a], K))
            Relative number field endomorphism of Number Field in b with defining polynomial x^3 - x + 1 over its base field
              Defn: b |--> -3/23*a*b^2 + (-9/46*a - 1/2)*b + 2/23*a
                    a |--> -a
        """
        K = self.domain().absolute_field('a')
        from_K, to_K = K.structure()
        a = from_K(K.gen())
        # We just have to figure out where a goes to
        # under the morphism defined by im_gen and base_map.
        L = self.codomain()
        R = L['x']
        f = R([base_map(x) for x in a.list()])
        b = f(im_gen)
        abs_hom = K.hom([b], check=check)
        return RelativeNumberFieldHomomorphism_from_abs(self, abs_hom)

    @cached_method
    def default_base_hom(self):
        r"""
        Pick an embedding of the base field of self into the codomain of this
        homset. This is done in an essentially arbitrary way.

        EXAMPLES::

            sage: L.<a, b> = NumberField([x^3 - x + 1, x^2 + 23])
            sage: M.<c> = NumberField(x^4 + 80*x^2 + 36)
            sage: Hom(L, M).default_base_hom()
            Ring morphism:
              From: Number Field in b with defining polynomial x^2 + 23
              To:   Number Field in c with defining polynomial x^4 + 80*x^2 + 36
              Defn: b |--> 1/12*c^3 + 43/6*c
        """
        v = self.domain().base_field().embeddings(self.codomain())
        if len(v) == 0:
            raise ValueError("no way to map base field to codomain.")
        return v[0]

    @cached_method
    def list(self):
        """
        Return a list of all the elements of self (for which the domain
        is a relative number field).

        EXAMPLES::

            sage: K.<a, b> = NumberField([x^2 + x + 1, x^3 + 2])
            sage: End(K).list()
            [
            Relative number field endomorphism of Number Field in a with defining polynomial x^2 + x + 1 over its base field
              Defn: a |--> a
                    b |--> b,
            ...
            Relative number field endomorphism of Number Field in a with defining polynomial x^2 + x + 1 over its base field
              Defn: a |--> a
                    b |--> -b*a - b
            ]

        An example with an absolute codomain::

            sage: K.<a, b> = NumberField([x^2 - 3, x^2 + 2])
            sage: Hom(K, CyclotomicField(24, 'z')).list()
            [
            Relative number field morphism:
              From: Number Field in a with defining polynomial x^2 - 3 over its base field
              To:   Cyclotomic Field of order 24 and degree 8
              Defn: a |--> z^6 - 2*z^2
                    b |--> -z^5 - z^3 + z,
            ...
            Relative number field morphism:
              From: Number Field in a with defining polynomial x^2 - 3 over its base field
              To:   Cyclotomic Field of order 24 and degree 8
              Defn: a |--> -z^6 + 2*z^2
                    b |--> z^5 + z^3 - z
            ]
        """
        D = self.domain()
        C = self.codomain()
        D_abs = D.absolute_field('a')
        v = [self(f, check=False) for f in D_abs.Hom(C).list()]
        return Sequence(v, universe=self, check=False, immutable=True, cr=v!=[])


class CyclotomicFieldHomset(NumberFieldHomset):
    """
    Set of homomorphisms with domain a given cyclotomic field.

    EXAMPLES::

        sage: End(CyclotomicField(16))
        Automorphism group of Cyclotomic Field of order 16 and degree 8
    """
    def __call__(self, im_gens, check=True):
        """
        Create an element of this homset.

        EXAMPLES::

            sage: K.<z> = CyclotomicField(16)
            sage: E = End(K)
            sage: E(E[0]) # indirect doctest
            Ring endomorphism of Cyclotomic Field of order 16 and degree 8
              Defn: z |--> z
            sage: E(z^5) # indirect doctest
            Ring endomorphism of Cyclotomic Field of order 16 and degree 8
              Defn: z |--> z^5
            sage: E(z^6) # indirect doctest
            Traceback (most recent call last):
            ...
            TypeError: images do not define a valid homomorphism
        """
        if isinstance(im_gens, CyclotomicFieldHomomorphism_im_gens):
            return self._coerce_impl(im_gens)
        try:
            return CyclotomicFieldHomomorphism_im_gens(self, im_gens, check=check)
        except (NotImplementedError, ValueError):
            try:
                return self._coerce_impl(im_gens)
            except TypeError:
                raise TypeError("images do not define a valid homomorphism")

    def _coerce_impl(self, x):
        r"""
        Coerce ``x`` into self. This will only work if ``x`` is already in self.

        EXAMPLES::

            sage: E = End(CyclotomicField(16))
            sage: E.coerce(E[0]) # indirect doctest
            Ring endomorphism of Cyclotomic Field of order 16 and degree 8
              Defn: zeta16 |--> zeta16
            sage: E.coerce(17) # indirect doctest
            Traceback (most recent call last):
            ...
            TypeError: no canonical coercion from Integer Ring to Automorphism group of Cyclotomic Field of order 16 and degree 8
        """
        if not isinstance(x, CyclotomicFieldHomomorphism_im_gens):
            raise TypeError
        if x.parent() is self:
            return x
        if x.parent() == self:
            return CyclotomicFieldHomomorphism_im_gens(self, x.im_gens())
        raise TypeError

    @cached_method
    def list(self):
        """
        Return a list of all the elements of self (for which the domain
        is a cyclotomic field).

        EXAMPLES::

            sage: K.<z> = CyclotomicField(12)
            sage: G = End(K); G
            Automorphism group of Cyclotomic Field of order 12 and degree 4
            sage: [g(z) for g in G]
            [z, z^3 - z, -z, -z^3 + z]
            sage: L.<a, b> = NumberField([x^2 + x + 1, x^4 + 1])
            sage: L
            Number Field in a with defining polynomial x^2 + x + 1 over its base field
            sage: Hom(CyclotomicField(12), L)[3]
            Ring morphism:
              From: Cyclotomic Field of order 12 and degree 4
              To:   Number Field in a with defining polynomial x^2 + x + 1 over its base field
              Defn: zeta12 |--> -b^2*a
            sage: list(Hom(CyclotomicField(5), K))
            []
            sage: Hom(CyclotomicField(11), L).list()
            []
        """
        D = self.domain()
        C = self.codomain()
        z = D.gen()
        n = z.multiplicative_order()
        if not n.divides(C.zeta_order()):
            v =[]
        else:
            if D == C:
                w = z
            else:
                w = C.zeta(n)
            v = [self([w**k], check=False) for k in Zmod(n) if k.is_unit()]
        return Sequence(v, universe=self, check=False, immutable=True, cr=v!=[])
