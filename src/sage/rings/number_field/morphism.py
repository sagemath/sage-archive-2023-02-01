from sage.rings.homset import RingHomset_generic
from sage.rings.morphism import RingHomomorphism_im_gens, RingHomomorphism
from sage.rings.integer import Integer
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
        D = self.domain()
        C = self.codomain().absolute_field()[0]
        f = D.absolute_polynomial()
        g = C['x'](f)
        r = g.roots()
        n = Integer(len(r))
        self.__order = n
        return n

    def list(self):
        """
        Return a list of all the elements of self.

        EXAMPLES:
            sage: K.<a> = NumberField( [x^2 + x + 1, x^3 + 2] )
            sage: L = K.absolute_field()[0]
            sage: G = End(L); G
            Automorphism group of Number Field in a0 with defining polynomial x^6 + 3*x^5 + 6*x^4 + 3*x^3 + 9*x + 9
            sage: G.order()
            6
            sage: G.list()
            [
            Ring endomorphism of Number Field in a0 with defining polynomial x^6 + 3*x^5 + 6*x^4 + 3*x^3 + 9*x + 9
              Defn: a0 |--> a0,
            ...
            Ring endomorphism of Number Field in a0 with defining polynomial x^6 + 3*x^5 + 6*x^4 + 3*x^3 + 9*x + 9
              Defn: a0 |--> -5/9*a0^5 - a0^4 - 2*a0^3 + 2/3*a0^2 - a0 - 5
            ]
        """
        try:
            return self.__list
        except AttributeError:
            pass
        D = self.domain()
        C = self.codomain()
        Dabs, from_Dabs, to_Dabs = self.domain().absolute_field()
        Cabs, from_Cabs, to_Cabs = self.codomain().absolute_field()
        f = Dabs.polynomial()
        g = Cabs['x'](f)
        r = g.roots()
        v = []
        for a, _ in r:
            im = Dabs.hom([from_Cabs(a)])(to_Dabs(D.gen()))
            v.append(self([im]))
        v = Sequence(v, immutable=True, cr=True)
        self.__list = v
        return v

    def __getitem__(self, n):
        return self.list()[n]


class NumberFieldHomomorphism_im_gens(RingHomomorphism_im_gens):
    pass


class RelativeNumberFieldHomset(NumberFieldHomset):
    """
    Set of homomorphisms with domain a given relative number field.

    EXAMPLES:
    We construct a homomorphism from a relative field by giving
    the image of a generator:

        sage: L.<cuberoot2, zeta3> = CyclotomicField(3).extension(x^3 - 2)
        sage: phi = L.hom([cuberoot2 * zeta3]); phi
        Relative number field endomorphism of Number Field in cuberoot2 with defining polynomial x^3 + -2 over its base field
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
            K, from_K, to_K = self.domain().absolute_field()
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
        K, from_K, to_K = self.domain().absolute_field()
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
                    b |--> (-b)*a + -b
            ]
        """
        try:
            return self.__list
        except AttributeError:
            pass
        D = self.domain()
        C = self.codomain()
        K = D.absolute_field()[0]
        v = K.Hom(C).list()
        w = [self(phi) for phi in v]
        w = Sequence(w, immutable=True, cr=True, universe=self)
        self.__list = w
        return w


class RelativeNumberFieldHomomorphism_from_abs(RingHomomorphism):
    def __init__(self, parent, abs_hom):
        RingHomomorphism.__init__(self, parent)
        self.__abs_hom = abs_hom
        K, from_K, to_K = self.domain().absolute_field()
        self.__K = K
        self.__from_K = K
        self.__to_K = K

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

    def __call__(self, x):
        return self.__abs_hom(self.__to_K(x))

