from sage.rings.homset import RingHomset_generic
from sage.rings.morphism import RingHomomorphism_im_gens
from sage.rings.integer import Integer
from sage.structure.sequence import Sequence

class FiniteFieldHomset(RingHomset_generic):
    """
    Set of homomorphisms with domain a given finite field.
    """
    def __call__(self, im_gens, check=True):
        """
        EXAMPLES:
            sage: R.<t> = ZZ[]
            sage: E.<a> = GF(25, modulus = t^2 - 2)
            sage: F.<b> = GF(625)
            sage: End(E)
            Automorphism group of Finite Field in a of size 5^2
            sage: list(Hom(E, F))
            [Ring morphism:
              From: Finite Field in a of size 5^2
              To:   Finite Field in b of size 5^4
              Defn: a |--> 4*b^3 + 4*b^2 + 4*b,
             Ring morphism:
              From: Finite Field in a of size 5^2
              To:   Finite Field in b of size 5^4
              Defn: a |--> b^3 + b^2 + b]
            sage: [phi(2*a)^2 for phi in Hom(E, F)]
            [3, 3]
            sage: End(GF(7))[0]
            Ring endomorphism of Finite Field of size 7
              Defn: 1 |--> 1
            sage: H = Hom(GF(7), GF(49, 'c'))
            sage: H[0](2)
            2
            sage: Hom(GF(49, 'c'), GF(7)).list()
            []
            sage: Hom(GF(49, 'c'), GF(81, 'd')).list()
            []
            sage: H = Hom(GF(9, 'a'), GF(81, 'b'))
            sage: H == loads(dumps(H))
            True
        """
        if isinstance(im_gens, FiniteFieldHomomorphism_im_gens):
            return self._coerce_impl(im_gens)
        try:
            return FiniteFieldHomomorphism_im_gens(self, im_gens, check=check)
        except (NotImplementedError, ValueError), err:
            try:
                return self._coerce_impl(im_gens)
            except TypeError:
                raise TypeError, "images do not define a valid homomorphism"

    def _coerce_impl(self, x):
        if not isinstance(x, FiniteFieldHomomorphism_im_gens):
            raise TypeError
        if x.parent() is self:
            return x
        if x.parent() == self:
            return FiniteFieldHomomorphism_im_gens(self, x.im_gens())
        raise TypeError

    def _repr_(self):
        """
        EXAMPLES:
            sage: Hom(GF(4, 'a'), GF(16, 'b'))
            Set of field embeddings from Finite Field in a of size 2^2 to Finite Field in b of size 2^4
            sage: Hom(GF(4, 'a'), GF(4, 'c'))
            Set of field embeddings from Finite Field in a of size 2^2 to Finite Field in c of size 2^2
            sage: Hom(GF(4, 'a'), GF(4, 'a'))
            Automorphism group of Finite Field in a of size 2^2
        """
        D = self.domain()
        C = self.codomain()
        if C == D:
            return "Automorphism group of %s"%D
        else:
            return "Set of field embeddings from %s to %s"%(D, C)

    def is_aut(self):
        """
        EXAMPLES:
            sage: Hom(GF(4, 'a'), GF(16, 'b')).is_aut()
            False
            sage: Hom(GF(4, 'a'), GF(4, 'c')).is_aut()
            False
            sage: Hom(GF(4, 'a'), GF(4, 'a')).is_aut()
            True
        """
        return self.domain() == self.codomain()

    def order(self):
        """
        Return the order of this set of field homomorphisms.

        EXAMPLES:
            sage: K.<a> = GF(125)
            sage: End(K)
            Automorphism group of Finite Field in a of size 5^3
            sage: End(K).order()
            3
            sage: L.<b> = GF(25)
            sage: Hom(L, K).order() == Hom(K, L).order() == 0
            True
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
        Return a list of all the elements in this set of field homomorphisms.

        EXAMPLES:
            sage: K.<a> = GF(25)
            sage: End(K)
            Automorphism group of Finite Field in a of size 5^2
            sage: list(End(K))
            [Ring endomorphism of Finite Field in a of size 5^2
              Defn: a |--> 4*a + 1,
             Ring endomorphism of Finite Field in a of size 5^2
              Defn: a |--> a]
            sage: L.<z> = GF(7^6)
            sage: [g for g in End(L) if (g^3)(z) == z]
            [Ring endomorphism of Finite Field in z of size 7^6
              Defn: z |--> 5*z^4 + 5*z^3 + 4*z^2 + 3*z + 1,
             Ring endomorphism of Finite Field in z of size 7^6
              Defn: z |--> 3*z^5 + 5*z^4 + 5*z^2 + 2*z + 3,
             Ring endomorphism of Finite Field in z of size 7^6
              Defn: z |--> z]
        """
        try:
            return self.__list
        except AttributeError:
            pass
        D = self.domain()
        C = self.codomain()
        if D.characteristic() == C.characteristic() and Integer(D.degree()).divides(Integer(C.degree())):
            f = D.modulus()
            g = C['x'](f)
            r = g.roots()
            v = [self(D.hom(a, C)) for a, _ in r]
            v = Sequence(v, immutable=True, cr=True)
        else:
            v = Sequence([], immutable=True, cr=False)
        self.__list = v
        return v

    def __getitem__(self, n):
        """
        EXAMPLES:
            sage: H = Hom(GF(32, 'a'), GF(1024, 'b'))
            sage: H[1]
            Ring morphism:
              From: Finite Field in a of size 2^5
              To:   Finite Field in b of size 2^10
              Defn: a |--> b^7 + b^5
            sage: H[2:4]
            [
            Ring morphism:
              From: Finite Field in a of size 2^5
              To:   Finite Field in b of size 2^10
              Defn: a |--> b^8 + b^6 + b^2,
            Ring morphism:
              From: Finite Field in a of size 2^5
              To:   Finite Field in b of size 2^10
              Defn: a |--> b^9 + b^7 + b^6 + b^5 + b^4
            ]
        """
        return self.list()[n]

    def index(self, item):
        """
        EXAMPLES:
            sage: K.<z> = GF(1024)
            sage: g = End(K)[3]
            sage: End(K).index(g) == 3
            True
        """
        return self.list().index(item)

class FiniteFieldHomomorphism_im_gens(RingHomomorphism_im_gens):
    pass


