from sage.structure.unique_representation import UniqueRepresentation
from sage.misc.decorators import sage_wraps
from sage.categories.pushout import pushout
from sage.categories.algebras import Algebras
from sage.rings.integer_ring import IntegerRing
from sage.rings.ring import CommutativeRing, CommutativeAlgebra
from sage.rings.morphism import RingHomomorphism

from sage.categories.action import Action
import operator


def coerce_algebras(method):
    @sage_wraps(method)
    def new_method(*algebras, **kwargs):
        other_factors = [ ]
        if base is None:
            base = self.base_ring()
        for i in len(algebras):
            algebra = algebras[i]
            if not isinstance(algebra, CommutativeAlgebra):
                raise TypeError("%s is not a commutative algebra" % algebra)
            if not isinstance(algebra, AlgebraFromMorphism):
                algebra = RingExtension(algebra)
            base_factor = algebra.base_ring()
            if base_factor is base:
                pass
            elif base_factor.has_coerce_map_from(base):
                morphism = base_factor.coerce_map_from(base)
                algebra = self.scalar_restriction(morphism)
            else:
                raise TypeError("No coercion map from %s to %s" % (over, base))
            if i == 0:
                first_factor = algebra
            else:
                other_factors.append(algebra)
        if first_factor is algebras[0]:
            return method(first_factor, *other_factors, **kwargs)
        else:
            return getattr(first_factor, method.__name__)(*other_factors, **kwargs)
    return new_method


# move the next class to element.pyx ???
class BaseActionOnRing(Action):
    def __init__(self, algebra):
        if not isinstance(algebra, AlgebraFromMorphism):
            raise TypeError("%s is not an instance of AlgebraFromMorphism" % algebra)
        self._base = algebra.base_ring()
        self._algebra = algebra
        self._defining_morphism = algebra.defining_morphism()
        Action.__init__(self, self._base, algebra, op=operator.mul)

    def _call_(self, scalar, element):
        result = self._defining_morphism(scalar) * element.element_in_ring()
        return self._algebra(result)


# Constructor for AlgebraFromMorphism
def RingExtension(ring, base=None, defining_morphism=None):
    r"""
    Create a ring extension

    INPUT:

    - ``ring`` -- a commutative ring

    - ``base`` (optional) -- a commutative ring lying below ``ring``

    - ``defining_morphism`` (optional) -- an homomorphism of rings
    from ``base`` to ``ring``

    OUTPUT:

    The ring ``ring`` viewed as an algebra over ``base`` through
    the homomorphism ``defining_morphism``.

    If ``base`` is not set the canonical base of ``ring`` (which
    is available via ``ring.base()``) is used.

    If ``defining_morphism`` is not set and ``base`` coerces to
    ``ring`` the coercion map is used for defining the extension.

    CREATION OF EXTENSIONS:

    We create an extension of finite fields::

        sage: K = GF(5^2); z2 = K.gen()
        sage: L = GF(5^4); z4 = L.gen()
        sage: E = RingExtension(L,K); E
        Finite Field in z4 of size 5^4 viewed as an algebra over Finite Field in z2 of size 5^2

    A shortcut for creating a ring extension is given by
    the ``/`` operator::

        sage: L/K
        Finite Field in z4 of size 5^4 viewed as an algebra over Finite Field in z2 of size 5^2
        sage: E is L/K
        True

    The ring and the base of ``E`` are accessible as follows::

        sage: E.base()
        Finite Field in z2 of size 5^2
        sage: E.base() is K
        True

        sage: E.ring()
        Finite Field in z4 of size 5^4
        sage: E.ring() is L
        True

    Here is an example where the base is implicit::

        sage: R.<x> = PolynomialRing(K)
        sage: RingExtension(R)
        Univariate Polynomial Ring in x over Finite Field in z2 of size 5^2 viewed as an algebra over Finite Field in z2 of size 5^2


    ELEMENTS IN EXTENSIONS:

    We can create elements in the extension `E` using standard methods::

        sage: E.zero()
        0
        sage: E.one()
        1
        sage: E.an_element()
        0
        sage: E.random_element()  # random
        4*z4^2 + 4*z4 + 2

    Conversion also works::

        sage: a = z4 + 1
        sage: a.parent()
        Finite Field in z4 of size 5^4
        sage: aE = E(a); aE
        z4 + 1
        sage: aE.parent()
        Finite Field in z4 of size 5^4 viewed as an algebra over Finite Field in z2 of size 5^2
        sage: aE.parent() is E
        True

    Of course all ring operations are available in E::

        sage: bE = (aE + 1) * (aE + 2); bE
        z4^2 + 1

    And the result stays in the extension::

        sage: bE.parent() 
        Finite Field in z4 of size 5^4 viewed as an algebra over Finite Field in z2 of size 5^2
        sage: bE.parent() is E
        True


    COERCION:

    A coercion map going from the extension to the ring (acting as the 
    identity) is set::

        sage: L.has_coerce_map_from(E)
        True
        sage: L.coerce_map_from(E)
        Ring Coercion morphism:
          From: Finite Field in z4 of size 5^4 viewed as an algebra over Finite Field in z2 of size 5^2
          To:   Finite Field in z4 of size 5^4

    Therefore we can safely add an element of the ring with an element 
    of the extension, obtaining this way an element in the ring::

        sage: c = z4^2 + 3
        sage: s = aE + c; s
        z4^2 + z4 + 4
        sage: s.parent()
        Finite Field in z4 of size 5^4
        sage: s.parent() is L
        True

    More generally, adding an element of ``E`` with an element lying in
    another finite field of caracteristic 5 works::

        sage: F = GF(5^3); z3 = F.gen()
        sage: aE + z3
        2*z12^11 + 2*z12^9 + 3*z12^8 + 4*z12^7 + z12^6 + 3*z12^5 + 3*z12^4 + 3*z12^3 + z12^2 + 3*z12 + 4

    Moreover the extension is equipped with a left action of the base.
    Thus the result of the multiplication of an element of the base by
    an element of the extension is an element of the extension::

        sage: l = z2*aE; l
        2*z4^3 + 3*z4^2 + 1
        sage: l.parent()
        Finite Field in z4 of size 5^4 viewed as an algebra over Finite Field in z2 of size 5^2

    On the other hand, the right action is not defined. Therefore when
    multiplying an element of the extension (on the left) by an element
    of the base (on the right) we get a result in the ring::

        sage: r = aE*z2; r
        2*z4^3 + 3*z4^2 + 1
        sage: r.parent()
        Finite Field in z4 of size 5^4

    However the two results of course agree after coercion::

        sage: l == r
        True

    Coercions between extensions are implemented as follows: an extension
    L1/K1 coerces to another extension L2/K2 when L1 coerces to L2 and K2
    coerces to K1 (be careful at the direction) and the defining morphisms
    agree. For example::

        sage: K1 = GF(3^4); L1 = GF(3^8);  E1 = L1/K1
        sage: K2 = GF(3^2); L2 = GF(3^16); E2 = L2/K2
        sage: E2.has_coerce_map_from(E1)
        True

        sage: x1 = E1.random_element()
        sage: x2 = E2.random_element()
        sage: s = x1 + x2
        sage: s.parent()
        Finite Field in z16 of size 3^16 viewed as an algebra over Finite Field in z2 of size 3^2
        sage: s.parent() is E2
        True


    DEFINING MORPHISM:

    When creating the extension E = L/K, we have not specified the defining
    morphism. In that case, the coercion map from the base to the ring (if
    it exists) is used by default::

        sage: E.defining_morphism()
        Ring morphism:
          From: Finite Field in z2 of size 5^2
          To:   Finite Field in z4 of size 5^4
          Defn: z2 |--> z4^3 + z4^2 + z4 + 3
        sage: E.defining_morphism() == L.coerce_map_from(K)
        True

    When there is no coercion map from the base to the ring, an error is
    raised::

        sage: L2 = GF(7^2)
        sage: RingExtension(L2,K)
        Traceback (most recent call last):
        ...
        ValueError: No coercion map from Finite Field in z2 of size 5^2 to Finite Field in z2 of size 7^2

    Of course it is also possible to use a custom defining morphism.
    As an example, let us create the algebra L/K viewed through the 
    Frobenius endomorphism.
    We first create the Frobenius map phi : K -> L::

        sage: FrobK = K.frobenius_endomorphism()
        sage: phi = FrobK.post_compose(L.coerce_map_from(K))
        sage: phi
        Composite map:
          From: Finite Field in z2 of size 5^2
          To:   Finite Field in z4 of size 5^4
          Defn:   Frobenius endomorphism z2 |--> z2^5 on Finite Field in z2 of size 5^2
                then
                  Ring morphism:
                  From: Finite Field in z2 of size 5^2
                  To:   Finite Field in z4 of size 5^4
                  Defn: z2 |--> z4^3 + z4^2 + z4 + 3

    And then the twisted extension::

        sage: ExtFrob = RingExtension(L, K, phi); ExtFrob
        Finite Field in z4 of size 5^4 viewed as an algebra over Finite Field in z2 of size 5^2

    Explicitely composing with coercion maps (on the left and on the right)
    is not mandatory: Sage does it automatically for you if necessary.
    Consequenty the extension ``ExtFrob`` can be created in a simpler way
    as follows::

        sage: ExtFrob2 = RingExtension(L, K, FrobK)
        sage: ExtFrob is ExtFrob2
        True

    We draw the attention of the user that non canonical defining morphisms 
    should be used with extreme care! We insist in particular on the following 
    tricky behaviour::

        sage: x = ExtFrob(z4)
        sage: l = z2*x; l
        4*z4^3 + 3*z4^2 + 2*z4 + 2
        sage: r = x*z2; r
        z4^3 + 2*z4^2 + 4*z4 + 3
        sage: l == r
        False

    Here is the explanatation.
    On the one hand, when computing the product ``z2*x`` the left action of the 
    base on the extension is used and ``z2`` is mapped to the extension through 
    the defining morphism.
    On the other hand, when computing the product ``x*z2``, both factors x and
    z2 are mapped to L through coercion maps and the product is computed in L.
    The defining morphism here never interferes.
    """
    if defining_morphism is None:
        coerce = True
        if base is None:
            base = ring.base_ring()
            if base is None:
                # ZZ is initial in the category of rings
                base = IntegerRing()
        if ring.has_coerce_map_from(base):
            defining_morphism = ring.coerce_map_from(base)
        else:
            raise ValueError("No coercion map from %s to %s" % (base,ring))
    else:
        domain = defining_morphism.domain()
        codomain = defining_morphism.codomain()
        coerce = (defining_morphism == codomain.coerce_map_from(domain))
        if base is None:
            base = domain
        elif base is domain:
            pass
        elif domain.has_coerce_map_from(base):
            coercion = domain.coerce_map_from(base)
            defining_morphism = defining_morphism.pre_compose(coercion)
        else:
            raise ValueError("No coercion map from %s to %s" % (base,domain))
        if ring is codomain:
            pass
        elif ring.has_coerce_map_from(codomain):
            coercion = ring.coerce_map_from(codomain)
            defining_morphism = defining_morphism.post_compose(coercion)
        else:
            raise ValueError("No coercion map from %s to %s" % (codomain,ring))
    return AlgebraFromMorphism(defining_morphism, coerce)


class AlgebraFromMorphism(CommutativeAlgebra, UniqueRepresentation):
    def __init__(self, defining_morphism, coerce):
        base = defining_morphism.domain()
        ring = defining_morphism.codomain()

        # To avoid problems...
        if isinstance(base, AlgebraFromMorphism) or isinstance(ring, AlgebraFromMorphism):
            raise NotImplementedError

        # We do not want to have a coercion map base -> self
        # So we first set the base to ring.base() (which indeed coerces to self)
        CommutativeAlgebra.__init__(self, ring.base(), category=Algebras(base)) 
        # and then update the base
        self._base = base
        self._ring = ring
        self._defining_morphism = defining_morphism
        self._coerce = coerce

        # Set left action of base on self
        # and a coercion map self -> ring
        from sage.rings.morphism import AlgebraToRing_coercion
        self._unset_coercions_used()
        self._populate_coercion_lists_(
            action_list = [BaseActionOnRing(self)],
            embedding = AlgebraToRing_coercion(self.Hom(ring), check=False))

        from sage.structure.element import AlgebraFMElement
        self.element_class = AlgebraFMElement

    def _element_constructor_(self, x, *args, **kwargs):
        from sage.structure.element import AlgebraFMElement
        if isinstance(x, AlgebraFMElement):
            x = x.element_in_ring()
        elt = self._ring(x, *args, **kwargs)
        return self.element_class(self, elt)

    def from_base_ring(self, x):
        elt = self._ring(x)
        return self.element_class(self, elt)

    def _pushout_(self, other):
        base = None
        if isinstance(other,AlgebraFromMorphism):
            # TODO: implement pushout when:
            # - defining_morphism are not coercion maps
            #   (this implies to check the equality of two morphisms)
            # - there is no direct coercion map between the bases
            #   but there exists a third ring which coerces to eash base
            #   Question: how can one find the greatest such third ring?
            if self._coerce and other._coerce:
                sbase = self._base
                obase = other._base
                if sbase.has_coerce_map_from(obase):
                    base = obase
                elif obase.has_coerce_map_from(sbase):
                    base = sbase
            ring = pushout(self._ring, other.ring())
        else:
            ring = pushout(self._ring, other)
        if base is None:
            return ring
        else:
            return RingExtension(ring,base)

    def _repr_(self):
        return "%s viewed as an algebra over %s" % (self._ring, self._base)

    def _coerce_map_from_(self, other):
        if isinstance(other, AlgebraFromMorphism):
            if self._coerce and other._coerce:
                return other.base().has_coerce_map_from(self._base) and self._ring.has_coerce_map_from(other.ring())                

    def defining_morphism(self):
        return self._defining_morphism

    def ring(self):
        return self._ring

    def _an_element_(self):
        elt = self._ring.an_element()
        return self.element_class(self, elt)

    def random_element(self):
        elt = self._ring.random_element()
        return self.element_class(self, elt)

    def scalar_restriction(self, morphism):
        coerce = self._coerce
        if isinstance(morphism, AlgebraFromMorphism):
            if coerce:
                coerce = morphism._coerce
            morphism = morphism.defining_morphism()
        elif isinstance(morphism, CommutativeRing):
            if self._base.has_coerce_map_from(morphism):
                morphism = self._base.coerce_map_from(morphism)
            else:
                raise TypeError("No coercion map from %s to %s" % (morphism, self._base))
        else:
            domain = morphism.domain()
            codomain = morphism.codomain()
            if coerce:
                coerce = (morphism == codomain.coerce_map_from(domain))
            if self._base.has_coerce_map_from(codomain):
                morphism = morphism.post_compose(self._base.coerce_map_from(codomain))
            else:
                raise TypeError("No coercion map from %s to %s" % (codomain, self._base))
        defining_morphism = self._defining_morphism.pre_compose(morphism)
        return AlgebraFromMorphism(defining_morphism, coerce)

    @coerce_algebras
    def tensor_product(self, other, base=None):
        # Implement this method in subclasses
        raise NotImplementedError("Tensor product is not implemented")

    def scalar_extension(self, morphism):
        if isinstance(morphism, AlgebraFromMorphism):
            extension = morphism
        elif isinstance(morphism, CommutativeRing):
            extension = RingExtension(morphism, self._base)
        else:
            domain = morphism.domain()
            codomain = morphism.codomain()
            coerce = (morphism == codomain.coerce_map_from(domain))
            if domain.has_coerce_map_from(self._base):
                morphism = morphism.pre_compose(domain.coerce_map_from(self._base))
            else:
                raise TypeError("No coercion map from %s to %s" % (self._base, domain))
            extension = AlgebraFromMorphism(morphism, coerce)
        algebra, morphisms = self.tensor_product(extension)
        return AlgebraFromMorphism(morphisms[1].ring(), self._coerce)
