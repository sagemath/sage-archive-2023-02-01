from sage.rings.homset import RingHomset_generic
from sage.rings.morphism import RingHomomorphism_im_gens, RingHomomorphism
from sage.rings.integer import Integer
from sage.rings.finite_rings.integer_mod_ring import Zmod
from sage.structure.sequence import Sequence

class NumberFieldHomset(RingHomset_generic):
    """
    Set of homomorphisms with domain a given number field.
    """
    def __call__(self, im_gens, check=True):
        """
        EXAMPLES:
            sage: H = Hom(ZZ, QQ)
            sage: phi = H([])
            Traceback (most recent call last):
            ...
            TypeError: images do not define a valid homomorphism

        TESTS:
            sage: H = Hom(ZZ, QQ)
            sage: H == loads(dumps(H))
            True
        """
        if isinstance(im_gens, NumberFieldHomomorphism_im_gens):
            return self._coerce_impl(im_gens)
        try:
            return NumberFieldHomomorphism_im_gens(self, im_gens, check=check)
        except (NotImplementedError, ValueError), err:
            try:
                return self._coerce_impl(im_gens)
            except TypeError:
                raise TypeError, "images do not define a valid homomorphism"

    def _coerce_impl(self, x):
        if not isinstance(x, NumberFieldHomomorphism_im_gens):
            raise TypeError
        if x.parent() is self:
            return x
        if x.parent() == self:
            return NumberFieldHomomorphism_im_gens(self, x.im_gens())
        raise TypeError

    def _repr_(self):
        D = self.domain()
        C = self.codomain()
        if C == D:
            return "Automorphism group of %s"%D
        else:
            return "Set of field embeddings from %s to %s"%(D, C)

    def is_aut(self):
        return self.domain() == self.codomain()

    def order(self):
        """
        Return the order of this set of field homomorphism.

        EXAMPLES:
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
        try:
            return self.__order
        except AttributeError:
            pass
        n = len(self.list())
        self.__order = n
        return n

    def list(self):
        """
        Return a list of all the elements of self.

        EXAMPLES:
            sage: K.<a> = NumberField( [x^2 + x + 1, x^3 + 2] )
            sage: L = K.absolute_field('b')
            sage: G = End(L); G
            Automorphism group of Number Field in b with defining polynomial x^6 + 3*x^5 + 6*x^4 + 3*x^3 + 9*x + 9
            sage: G.order()
            6
            sage: G.list()
            [
            Ring endomorphism of Number Field in b with defining polynomial x^6 + 3*x^5 + 6*x^4 + 3*x^3 + 9*x + 9
              Defn: b |--> b,
            ...
            Ring endomorphism of Number Field in b with defining polynomial x^6 + 3*x^5 + 6*x^4 + 3*x^3 + 9*x + 9
              Defn: b |--> -5/9*b^5 - b^4 - 2*b^3 + 2/3*b^2 - b - 5
            ]
            sage: Hom(L, CyclotomicField(3)).list()
            []
        """
        try:
            return self.__list
        except AttributeError:
            pass
        D = self.domain()
        C = self.codomain()
        v = []
        if Integer(D.absolute_degree()).divides(Integer(C.absolute_degree())):
            Dabs = D.absolute_field(D.variable_name())
            from_Dabs, to_Dabs = Dabs.structure()
            Cabs = C.absolute_field(C.variable_name())
            from_Cabs, to_Cabs = Cabs.structure()
            f = Dabs.polynomial()
            g = Cabs['x'](f)
            r = g.roots()
            for a, _ in r:
                im = Dabs.hom([from_Cabs(a)])(to_Dabs(D.gen()))
                v.append(self([im]))
        v = Sequence(v, immutable=True, cr=v!=[])
        self.__list = v
        return v

    def __getitem__(self, n):
        return self.list()[n]


class NumberFieldHomomorphism_im_gens(RingHomomorphism_im_gens):
    def __invert__(self):
        r"""

        Return the inverse of an isomorphism of absolute number fields

        EXAMPLES::

            sage: K.<a> = NumberField(x^2 + 5)
            sage: tau1, tau2 = K.automorphisms(); tau1, tau2
            (Ring endomorphism of Number Field in a with defining polynomial x^2 + 5
              Defn: a |--> a,
             Ring endomorphism of Number Field in a with defining polynomial x^2 + 5
              Defn: a |--> -a)
            sage: ~tau1
            Ring endomorphism of Number Field in a with defining polynomial x^2 + 5
             Defn: a |--> a
            sage: ~tau2
            Ring endomorphism of Number Field in a with defining polynomial x^2 + 5
             Defn: a |--> -a

            sage: L.<z> = CyclotomicField(5)
            sage: tau1, tau2, tau3, tau4 = L.automorphisms()
            sage: (tau1, ~tau1)
            (Ring endomorphism of Cyclotomic Field of order 5 and degree 4
              Defn: z |--> z,
             Ring endomorphism of Cyclotomic Field of order 5 and degree 4
              Defn: z |--> z)
            sage: (tau2, ~tau2)
            (Ring endomorphism of Cyclotomic Field of order 5 and degree 4
              Defn: z |--> z^2,
             Ring endomorphism of Cyclotomic Field of order 5 and degree 4
              Defn: z |--> z^3)
            sage: (tau4, ~tau4)
            (Ring endomorphism of Cyclotomic Field of order 5 and degree 4
              Defn: z |--> z^3,
             Ring endomorphism of Cyclotomic Field of order 5 and degree 4
              Defn: z |--> z^2)

             sage: M.<w> = NumberField(x^4 - 5*x + 5)
             sage: phi = M.hom([z - z^2]); phi
             Ring morphism:
               From: Number Field in w with defining polynomial x^4 - 5*x + 5
               To:   Cyclotomic Field of order 5 and degree 4
               Defn: w |--> -z^2 + z
             sage: phi^-1
             Ring morphism:
               From: Cyclotomic Field of order 5 and degree 4
               To:   Number Field in w with defining polynomial x^4 - 5*x + 5
               Defn: z |--> 3/11*w^3 + 4/11*w^2 + 9/11*w - 14/11

        """
        K = self.domain()
        L = self.codomain()
        if K.degree() != L.degree():
            raise TypeError, "Can only invert isomorphisms"
        V, V_into_K, _ = K.vector_space()
        _, _, L_into_W = L.vector_space()
        linear_inverse = ~V.hom(map(L_into_W*self*V_into_K, V.basis()))
        return L.hom(map(V_into_K*linear_inverse*L_into_W, [L.gen()]))


class RelativeNumberFieldHomset(NumberFieldHomset):
    """
    Set of homomorphisms with domain a given relative number field.

    EXAMPLES:
    We construct a homomorphism from a relative field by giving
    the image of a generator:

        sage: L.<cuberoot2, zeta3> = CyclotomicField(3).extension(x^3 - 2)
        sage: phi = L.hom([cuberoot2 * zeta3]); phi
        Relative number field endomorphism of Number Field in cuberoot2 with defining polynomial x^3 - 2 over its base field
          Defn: cuberoot2 |--> zeta3*cuberoot2
                zeta3 |--> zeta3
        sage: phi(cuberoot2 + zeta3)
        zeta3*cuberoot2 + zeta3

    In fact, this phi is a generator for the Kummer Galois group of this
    cyclic extension:
        sage: phi(phi(cuberoot2 + zeta3))
        (-zeta3 - 1)*cuberoot2 + zeta3
        sage: phi(phi(phi(cuberoot2 + zeta3)))
        cuberoot2 + zeta3
    """
    def __call__(self, im_gen, base_hom=None, check=True):
        if isinstance(im_gen, NumberFieldHomomorphism_im_gens):
            # Then it must be a homomorphism from the corresponding
            # absolute number field
            abs_hom = im_gen
            K = abs_hom.domain()
            if K != self.domain().absolute_field(K.variable_name()):
                raise TypeError, "domain of morphism must be absolute field of domain."
            from_K, to_K = K.structure()
            if abs_hom.domain() != K:
                raise ValueError, "domain of absolute homomorphism must be absolute field of domain."
            if abs_hom.codomain() != self.codomain():
                raise ValueError, "codomain of absolute homomorphism must be codomain of this homset."
            return RelativeNumberFieldHomomorphism_from_abs(self, abs_hom)
        if isinstance(im_gen, RelativeNumberFieldHomomorphism_from_abs):
            return self._coerce_impl(im_gen)
        if base_hom is None:
            base_hom = self.default_base_hom()
        if isinstance(im_gen, (list, tuple)) and len(im_gen) == 1:
            im_gen = im_gen[0]
        if check:
            im_gen = self.codomain()(im_gen)
        return self._from_im(im_gen, base_hom)

    def _coerce_impl(self, x):
        if not isinstance(x, RelativeNumberFieldHomomorphism_from_abs):
            raise TypeError
        if x.parent() is self:
            return x
        if x.parent() == self:
            return RelativeNumberFieldHomomorphism_from_abs(self, x.abs_hom())
        raise TypeError

    def _from_im(self, im_gen, base_hom):
        """
        Return the homomorphism that acts on the base as given and
        sends the generator of the domain to im_gen.
        """
        K = self.domain().absolute_field('a')
        from_K, to_K = K.structure()
        a = from_K(K.gen())
        # We just have to figure out where a goes to
        # under the morphism defined by im_gen and base_hom.
        L = self.codomain()
        R = L['x']
        f = R([base_hom(x) for x in a.list()])
        b = f(im_gen)
        abs_hom = K.hom([b])
        return RelativeNumberFieldHomomorphism_from_abs(self, abs_hom)

    def default_base_hom(self):
        try:
            return self.__default_base_hom
        except AttributeError:
            pass
        v = self.domain().base_field().embeddings(self.codomain())
        if len(v) == 0:
            raise ValueError, "no way to map base field to codomain."
        self.__default_base_hom = v[0]
        return v[0]

    def list(self):
        """
        Return a list of all the elements of self.

        EXAMPLES:
            sage: K.<a, b> = NumberField( [x^2 + x + 1, x^3 + 2] )
            sage: G = End(K); G
            Automorphism group of Number Field in a with defining polynomial x^2 + x + 1 over its base field
            sage: v = G.list(); v
            [
            Relative number field endomorphism of Number Field in a with defining polynomial x^2 + x + 1 over its base field
              Defn: a |--> a
                    b |--> b,
            ...
            Relative number field endomorphism of Number Field in a with defining polynomial x^2 + x + 1 over its base field
              Defn: a |--> a
                    b |--> -b*a - b
            ]
        """
        try:
            return self.__list
        except AttributeError:
            pass
        D = self.domain()
        C = self.codomain()
        K = D.absolute_field('a')
        v = K.Hom(C).list()
        w = [self(phi) for phi in v]
        w = Sequence(w, immutable=True, cr=w!=[], universe=self)
        self.__list = w
        return w


class RelativeNumberFieldHomomorphism_from_abs(RingHomomorphism):
    def __init__(self, parent, abs_hom):
        RingHomomorphism.__init__(self, parent)
        self.__abs_hom = abs_hom
        K = abs_hom.domain()
        from_K, to_K = K.structure()
        self.__K = K
        self.__from_K = from_K
        self.__to_K = to_K

    def abs_hom(self):
        return self.__abs_hom

    def _repr_type(self):
        return "Relative number field"

    def im_gens(self):
        try:
            return self.__im_gens
        except AttributeError:
            pass
        D = self.domain()
        v = Sequence([self(x) for x in D.gens()], immutable=True)
        self.__im_gens = v
        return v

    def _repr_defn(self):
        D = self.domain()
        ig = self.im_gens()
        return '\n'.join(['%s |--> %s'%(D.gen(i), ig[i]) for\
                       i in range(D.ngens())])

    def _call_(self, x):
        return self.__abs_hom(self.__to_K(x))


class CyclotomicFieldHomset(NumberFieldHomset):
    """
    Set of homomorphisms with domain a given cyclotomic field.
    """
    def __call__(self, im_gens, check=True):
        """
        EXAMPLES:
            sage: End(CyclotomicField(16))
            Automorphism group of Cyclotomic Field of order 16 and degree 8
        """
        if isinstance(im_gens, CyclotomicFieldHomomorphism_im_gens):
            return self._coerce_impl(im_gens)
        try:
            return CyclotomicFieldHomomorphism_im_gens(self, im_gens, check=check)
        except (NotImplementedError, ValueError), err:
            try:
                return self._coerce_impl(im_gens)
            except TypeError:
                raise TypeError, "images do not define a valid homomorphism"

    def _coerce_impl(self, x):
        if not isinstance(x, CyclotomicFieldHomomorphism_im_gens):
            raise TypeError
        if x.parent() is self:
            return x
        if x.parent() == self:
            return CyclotomicFieldHomomorphism_im_gens(self, x.im_gens())
        raise TypeError

    def list(self):
        """
        Return a list of all the elements of self.

        EXAMPLES:
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
        try:
            return self.__list
        except AttributeError:
            pass

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
            v = [self([w**k]) for k in Zmod(n) if k.is_unit()]
        v = Sequence(v, immutable=True, cr=v!=[])
        self.__list = v
        return v


class CyclotomicFieldHomomorphism_im_gens(NumberFieldHomomorphism_im_gens):
    pass



