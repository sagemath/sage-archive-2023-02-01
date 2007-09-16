from sage.rings.homset import RingHomset_generic
from sage.rings.morphism import RingHomomorphism_im_gens
from sage.rings.integer import Integer
from sage.structure.sequence import Sequence

class NumberFieldHomset(RingHomset_generic):
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
            sage: K.<a> = NumberField( [x^3 + 2, x^2 + x + 1] )
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

            sage: G = End(K); G
            Automorphism group of Number Field in a0 with defining polynomial x^2 + x + 1 over its base field
            sage: v = G.list()
            sage: v[1]
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
            print im
            v.append(self([im]))
        v = Sequence(v, immutable=True, cr=True)
        self.__list = v
        return v

    def __getitem__(self, n):
        return self.list()[n]


class NumberFieldHomomorphism_im_gens(RingHomomorphism_im_gens):
    pass
