r"""
Extension of rings

See :func:`RingExtension` for documentation.

AUTHOR:

- Xavier Caruso (2016)
"""

#############################################################################
#    Copyright (C) 2019 Xavier Caruso <xavier.caruso@normalesup.org>
#
#    This program is free software: you can redistribute it and/or modify
#    it under the terms of the GNU General Public License as published by
#    the Free Software Foundation, either version 2 of the License, or
#    (at your option) any later version.
#                  http://www.gnu.org/licenses/
#****************************************************************************


from sage.misc.cachefunc import cached_method

from sage.structure.factory import UniqueFactory
from sage.structure.category_object import normalize_names
from sage.categories.map import Map
from sage.categories.commutative_rings import CommutativeRings
from sage.categories.commutative_algebras import CommutativeAlgebras
from sage.categories.fields import Fields
from sage.categories.pushout import pushout
from sage.rings.integer_ring import ZZ
from sage.rings.infinity import Infinity
from sage.rings.ring import CommutativeRing, CommutativeAlgebra

from sage.rings.ring_extension_element import RingExtensionElement
from sage.rings.ring_extension_element import RingExtensionWithBasisElement


# Helper function

def _common_base(K, L, degree=False):
    def tower_bases(ring):
        bases = [ ]
        base = ring
        deg = 1
        degrees = [ ]
        while True:
            bases.append(base)
            if degree:
                degrees.append(deg)
                try:
                    d = base.relative_degree()
                except AttributeError:
                    try:
                        d = base.degree()
                    except AttributeError:
                        break
                if d is Infinity: break
                deg *= d
            newbase = base.base_ring()
            if newbase is base: break
            base = newbase
        return bases, degrees
    bases_K, degrees_K = tower_bases(K)
    bases_L, degrees_L = tower_bases(L)
    base = None
    for iL in range(len(bases_L)):
        try:
            iK = bases_K.index(bases_L[iL])
            base = bases_L[iL]
            break
        except ValueError:
            pass
    if base is None:
        raise NotImplementedError("unable to find a common base")
    if degree:
        return base, degrees_K[iK], degrees_L[iL]
    else:
        return base


# Factory
#########

class RingExtensionFactory(UniqueFactory):
    r"""
    Create a ring extension

    INPUT:

    - ``ring`` -- a commutative ring

    - ``base`` (optional) -- a commutative ring lying below ``ring``

    - ``defining_morphism`` (optional) -- an homomorphism of rings
    from ``base`` to ``ring``

    OUTPUT:

    The ring ``ring`` over ``base`` through
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
        Finite Field in z4 of size 5^4 over its base

    The base of ``E`` is accessible as follows::

        sage: E.base_ring()
        Finite Field in z2 of size 5^2
        sage: E.base_ring() is K
        True

    Here is an example where the base is implicit::

        sage: R.<x> = PolynomialRing(K)
        sage: RingExtension(R)
        Univariate Polynomial Ring in x over Finite Field in z2 of size 5^2 over its base


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
        Finite Field in z4 of size 5^4 over its base
        sage: aE.parent() is E
        True

    Of course all ring operations are available in E::

        sage: bE = (aE + 1) * (aE + 2); bE
        z4^2 + 1

    And the result stays in the extension::

        sage: bE.parent()
        Finite Field in z4 of size 5^4 over its base
        sage: bE.parent() is E
        True


    COERCION:

    A coercion map going from the base to the extension is set::

        sage: E.has_coerce_map_from(K)
        True
        sage: E.coerce_map_from(K)
        Ring morphism:
          From: Finite Field in z2 of size 5^2
          To:   Finite Field in z4 of size 5^4 over its base

    Therefore we can safely add an element of the base with an element
    of the extension, obtaining this way an element in the extension::

        sage: c = z2 + 3
        sage: s = aE + c; s
        z4^3 + z4^2 + 2*z4 + 2
        sage: s.parent()
        Finite Field in z4 of size 5^4 over its base
        sage: s.parent() is E
        True

    When the defining morphism of the extension ``E = L/K`` is a coercion
    map (see below), a coercion from ``E`` to ``L`` (acting as the identity)
    is set::

        sage: L.has_coerce_map_from(E)
        True
        sage: L.coerce_map_from(E)
        Ring morphism:
          From: Finite Field in z4 of size 5^4 over its base
          To:   Finite Field in z4 of size 5^4


    Coercions between extensions are implemented as follows: an extension
    L1/K1 coerces to another extension L2/K2 when L1 coerces to L2 and K2
    coerces to K1 (be careful at the direction) and the defining morphisms
    agree. For example::

        sage: K1 = GF(3^4); L1 = GF(3^8);  E1 = RingExtension(L1,K1)
        sage: K2 = GF(3^2); L2 = GF(3^16); E2 = RingExtension(L2,K2)
        sage: E2.has_coerce_map_from(E1)
        True

        sage: x1 = E1.random_element()
        sage: x2 = E2.random_element()
        sage: s = x1 + x2
        sage: s.parent()
        Finite Field in z16 of size 3^16 over its base
        sage: s.parent() is E2
        True


    DEFINING MORPHISM:

    When creating the extension ``E = L/K``, we have not specified the defining
    morphism. In that case, the coercion map from the base to the ring (if it
    exists) is used by default::

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

    It is also possible to use a custom defining morphism.
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
        Finite Field in z4 of size 5^4 over its base

    Explicitely composing with coercion maps (on the left and on the right)
    is not mandatory: Sage does it automatically for you if necessary.
    Consequenty the extension ``ExtFrob`` can be created in a simpler way
    as follows::

        sage: ExtFrob2 = RingExtension(L, K, FrobK)
        sage: ExtFrob is ExtFrob2
        True

    We draw the attention of the user that non canonical defining morphisms
    should be used with extreme care!
    As an example, in the example below, Sage defines a coercion map from
    ``K`` to ``ExtFrob`` but no coercion map from ``ExtFrob`` to ``L`` (since
    otherwise the Frobenius morphism from ``K`` to ``L`` would be promoted as
    a coercion map as well).

        sage: ExtFrob.has_coerce_map_from(K)
        True
        sage: ExtFrob.coerce_map_from(K)
        Ring morphism:
          From: Finite Field in z2 of size 5^2
          To:   Finite Field in z4 of size 5^4 over its base

        sage: L.has_coerce_map_from(ExtFrob)
        False

        sage: aE = ExtFrob(a)
        sage: aE + c
        4*z4^3 + 4*z4^2 + 2
        sage: a + c^5
        4*z4^3 + 4*z4^2 + 2


    AUTHOR:

    - Xavier Caruso (2016)
    """
    def create_key_and_extra_args(self, *args, **kwargs):
        from sage.rings.ring_extension_morphism import RingExtensionHomomorphism

        rings = [ ]
        defining_morphism = None
        for i in range(len(args)):
            if args[i] in CommutativeRings():
                rings.append(args[i])
            elif isinstance(args[i], Map):
                if i < len(args) - 1:
                    raise TypeError("the defining morphism must be the last argument")
                defining_morphism = args[i]
            else:
                raise TypeError
        if len(rings) == 0 and defining_morphism is not None:
            base = defining_morphism.domain()
            ring = defining_morphism.codomain()
        elif len(rings) == 1:
            ring = rings[0]
            base = ring.base_ring()
        elif len(rings) == 2:
            ring, base = rings
        else:
            raise TypeError("you must specify at least a ring or a defining morphism")
        if defining_morphism is None:
            if isinstance(base, RingExtension_class):
                backend_base = base._backend
                if ring.has_coerce_map_from(backend_base):
                    defining_morphism = RingExtensionHomomorphism(base.Hom(ring), ring.coerce_map_from(backend_base))
            else:
                if ring.has_coerce_map_from(base):
                    defining_morphism = ring.coerce_map_from(base)
            if defining_morphism is None:
                raise ValueError("No coercion map from %s to %s" % (base,ring))
        else:
            if defining_morphism.domain() is not base:
                defining_morphism = defining_morphism.extend_domain(base)
            if defining_morphism.codomain() is not ring:
                defining_morphism = defining_morphism.extend_codomain(ring)
        # if not isinstance(base, RingExtension_class) and base is not base.base_ring():
        #     base = RingExtension(base)
        #     defining_morphism = RingExtensionHomomorphism(base.Hom(ring), defining_morphism)
        if isinstance(ring, RingExtension_class):
            from sage.rings.ring_extension_morphism import backend_morphism
            defining_morphism = backend_morphism(defining_morphism, forget="codomain")
            ring = ring._backend

        gen = basis = names = None
        if 'names' in kwargs:
            names = kwargs['names']
        elif 'name' in kwargs:
            names = kwargs['name']
        if 'gens' in kwargs or 'gen' in kwargs:
            if 'gens' in kwargs:
                gens = kwargs['gens']
                if not isinstance(gens, (list, tuple)):
                    raise TypeError("gens must be a list or a tuple")
                if len(gens) > 1:
                    raise NotImplementedError("only ring extensions with a single generator are implemented")
                gen = ring(gens[0])
            if 'gen' in kwargs:
                gen = ring(kwargs['gen'])
            if names is None:
                raise TypeError("you must specify the name of the generator")
            names = normalize_names(1, names)
        elif 'basis' in kwargs:
            basis = kwargs['basis']
            if not isinstance(basis, (list, tuple)):
                raise TypeError("basis must be a list or a tuple")
            if names is not None and len(basis) != len(names):
                raise ValueError("the number of names does not match the length of the basis")
            basis = [ ring(x) for x in kwargs['basis'] ]
            names = tuple(names)
        if gen is None and basis is None:
            if ring.ngens() == 1:
                try:
                    if isinstance(base, RingExtension_class):
                        _ = _common_base(base._backend, ring, degree=True)
                    else:
                        _ = _common_base(base, ring, degree=True)
                    gen = ring.gen()
                    if names is None:
                        names = ring.variable_names()
                    elif isinstance(names, (list, tuple)):
                        if len(names) != 1:
                            raise ValueError("the number of names does not match the number of generators")
                        names = tuple(names)
                    else:
                        names = (names,)
                except NotImplementedError:
                    pass
        extra_args = { }
        if 'repr' in kwargs:
            extra_args['repr'] = str(kwargs['repr'])
        return (defining_morphism, gen, basis, names), extra_args

    def create_object(self, version, key, **extra_args):
        (defining_morphism, gen, basis, names) = key
        if gen is not None:
            return RingExtensionWithGen(defining_morphism, gen, names[0])
        elif basis is not None:
            return RingExtensionWithBasis(defining_morphism, basis, names)
        else:
            return RingExtension_class(defining_morphism, **extra_args)


RingExtension = RingExtensionFactory("RingExtension")



# General extensions
####################

class RingExtension_class(CommutativeAlgebra):
    r"""
    Create a ring extension

    Do not call this function directly
    Use :func:`RingExtension` instead

    INPUT:

    - ``defining_morphism`` -- a ring homomorphism

    - ``coerce`` -- boolean (specify whether defining_morphism
    is a coercion map or not)

    OUTPUT:

    The extension defined by ``defining_morphism``

    EXAMPLES::

        sage: K = GF(5^2)
        sage: L = GF(5^4)
        sage: E = RingExtension(L,K); E
        Finite Field in z4 of size 5^4 over its base

        sage: from sage.rings.ring_extension import RingExtension_class
        sage: isinstance(E, RingExtension_class)
        True

    See :func:`RingExtension` for more documentation

    AUTHOR:

    - Xavier Caruso (2016)
    """
    Element = RingExtensionElement

    def __init__(self, defining_morphism, repr=None):
        r"""
        TESTS::

        The attribute _coerce indicates if the defining morphism is
        a coercion map::

            sage: K = GF(5^2)
            sage: E = RingExtension(K, K, K.frobenius_endomorphism())
            sage: E._coerce
            False

            sage: L = GF(5^4)
            sage: E1 = RingExtension(L,K)
            sage: E1._coerce
            True

        """
        from sage.rings.ring_extension_morphism import RingExtensionHomomorphism
        from sage.rings.ring_extension_morphism import RingExtensionBackendIsomorphism
        from sage.rings.ring_extension_morphism import RingExtensionBackendReverseIsomorphism
        from sage.rings.ring_extension_morphism import backend_morphism, are_equal_morphisms

        base = defining_morphism.domain()
        ring = defining_morphism.codomain()
        CommutativeAlgebra.__init__(self, ZZ, category=CommutativeAlgebras(base))
        self._base = base
        self._backend = ring
        self._repr = repr
        self._backend_defining_morphism = defining_morphism
        self._defining_morphism = RingExtensionHomomorphism(self._base.Hom(self), defining_morphism)

        # Some checkings
        if (base not in CommutativeRings()
         or ring not in CommutativeRings()
         or not defining_morphism.category_for().is_subcategory(CommutativeRings())):
            raise TypeError("only commutative rings are accepted")
        f = ring.Hom(ring).identity()
        b = self
        while isinstance(b, RingExtension_class):
            f *= backend_morphism(b._backend_defining_morphism)
            b = b.base_ring()
            if isinstance(b, RingExtension_class):
                backend = b._backend
            else:
                backend = b
            if ring.has_coerce_map_from(backend) and not are_equal_morphisms(f, None):
                # TODO: find a better message
                raise ValueError("exotic defining morphism between two rings in the tower; consider using another variable name")

        # We register coercion/conversion maps
        self._unset_coercions_used()
        self.register_coercion(self._defining_morphism.__copy__())
        self.register_coercion(RingExtensionBackendIsomorphism(ring.Hom(self)))
        ring._unset_coercions_used()
        ring.register_conversion(RingExtensionBackendReverseIsomorphism(self.Hom(ring)))


    def __hash__(self):
        return id(self)

    def from_base_ring(self, r):
        r"""
        Return the canonical embedding of ``r`` into ``self``.

        INPUT:

        - ``r`` -- an element of the base of the ring of this extension

        EXAMPLES::

            sage: k = GF(5)
            sage: K = GF(5^2); z2 = K.gen()
            sage: L = GF(5^4); z4 = L.gen()
            sage: E = RingExtension(L,K); E
            Finite Field in z4 of size 5^4 over its base

            sage: E.base()
            Finite Field in z2 of size 5^2
            sage: E._backend.base()
            Finite Field of size 5

            sage: x = E.from_base_ring(k(2)); x
            2
            sage: x.parent()
            Finite Field in z4 of size 5^4 over its base

            sage: x = E.from_base_ring(z2); x
            z4^3 + z4^2 + z4 + 3

        TESTS::

            sage: R.<t> = K[]
            sage: F = RingExtension(R, K, K.frobenius_endomorphism())
            sage: F.from_base_ring(z2)
            4*z2 + 1
            sage: F.defining_morphism()(z2)
            4*z2 + 1
        """
        if r not in self._base:
            raise TypeError("%s is not an element of the base of %s (= %s)" % (r, self._backend, self._base))
        return self.element_class(self, r)

    def _repr_(self):
        r"""
        Return a string representation of this extension.

        EXAMPLES::

            sage: K = GF(5^2)
            sage: L = GF(5^4)
            sage: E = RingExtension(L,K)
            sage: E
            Finite Field in z4 of size 5^4 over its base

            sage: E._repr_()
            'Finite Field in z4 of size 5^4 over its base'
        """
        if self._repr is not None:
            return self._repr
        else:
            return "%s over its base" % self._backend

    def _latex_(self):
        r"""
        Return a latex representation of this extension.

        EXAMPLES::

            sage: K = GF(5^2)
            sage: L = GF(5^4)
            sage: E = RingExtension(L,K)
            sage: latex(E)
            \Bold{F}_{5^{4}}/\Bold{F}_{5^{2}}
        """
        from sage.misc.latex import latex
        return "%s/%s" % (latex(self._backend), latex(self._base))

    def _coerce_map_from_(self, other):
        from sage.rings.ring_extension_morphism import RingExtensionHomomorphism
        from sage.rings.ring_extension_morphism import backend_morphism, are_equal_morphisms
        if isinstance(other, RingExtension_class):
           if self._backend.has_coerce_map_from(other._backend) and self._base.has_coerce_map_from(other._base):
                backend = self._backend.coerce_map_from(other._backend)
                f = backend * backend_morphism(other._defining_morphism)
                g = backend_morphism(self._defining_morphism * self._base.coerce_map_from(other._base))
                if are_equal_morphisms(f, g):
                    return RingExtensionHomomorphism(other.Hom(self), backend)

    def base(self):
        r"""
        Return the base of this extension

        EXAMPLES::

            sage: K = GF(5^2); z2 = K.gen()
            sage: L = GF(5^4); z4 = L.gen()
            sage: E = RingExtension(L,K); E
            Finite Field in z4 of size 5^4 over its base
            sage: E.base()
            Finite Field in z2 of size 5^2
            sage: E.base() is K
            True

        Note that the base of the extension is generally different
        from the base of the ring of the extension::

            sage: E._backend.base()
            Finite Field of size 5
        """
        return self._base

    def bases(self):
        L = [ self ]
        base = self
        while isinstance(base, RingExtension_class):
            base = base.base_ring()
            L.append(base)
        return L

    def is_defined_over(self, base):
        b = self
        while isinstance(b, RingExtension_class):
            if b is base: return True
            b = b.base_ring()
        return b == base

    def _check_base(self, base):
        if base is None:
            return self._base
        if not self.is_defined_over(base):
            raise ValueError("not (explicitely) defined over %s" % base)
        return base

    def defining_morphism(self, base=None):
        r"""
        Return the defining morphism of this extension

        EXAMPLES::

            sage: K = GF(5^2); z2 = K.gen()
            sage: L = GF(5^4); z4 = L.gen()
            sage: E = RingExtension(L,K); E
            Finite Field in z4 of size 5^4 over its base
            sage: E.defining_morphism()
            Ring morphism:
              From: Finite Field in z2 of size 5^2
              To:   Finite Field in z4 of size 5^4
              Defn: z2 |--> z4^3 + z4^2 + z4 + 3

            sage: E2 = RingExtension(L, K, K.frobenius_endomorphism()); E2
            Finite Field in z4 of size 5^4 over its base
            sage: E2.defining_morphism()
            Composite map:
              From: Finite Field in z2 of size 5^2
              To:   Finite Field in z4 of size 5^4
              Defn:   Frobenius endomorphism z2 |--> z2^5 on Finite Field in z2 of size 5^2
                    then
                      Ring morphism:
                      From: Finite Field in z2 of size 5^2
                      To:   Finite Field in z4 of size 5^4
                      Defn: z2 |--> z4^3 + z4^2 + z4 + 3
        """
        base = self._check_base(base)
        return self.coerce_map_from(base)

    def _an_element_(self):
        r"""
        Return an element of this extension

        This element is constructed by converting an element of
        the ring of this extension to the extension itself

        TESTS::

            sage: K = GF(5^2)
            sage: L = GF(5^4)
            sage: E = RingExtension(L,K); E
            Finite Field in z4 of size 5^4 over its base

            sage: x = E.an_element()  # indirect doctest
            sage: x
            0
            sage: x.parent()
            Finite Field in z4 of size 5^4 over its base

            sage: a = E._backend.an_element()
            sage: a == x
            True
        """
        elt = self._backend.an_element()
        return self.element_class(self, elt)

    def gens(self, base=None):
        self._check_base(base)
        return tuple([ self(x) for x in self._backend.gens() ])

    def ngens(self, base=None):
        return len(self.gens(base))

    def gen(self):
        r"""
        Return a generator of this extension

        This element is constructed by converting a generator of
        the ring of this extension to the extension itself

        TESTS::

            sage: K = GF(5^2)
            sage: L = GF(5^4)
            sage: E = RingExtension(L,K); E
            Finite Field in z4 of size 5^4 over its base

            sage: x = E.gen()
            sage: x
            z4
            sage: x.parent()
            Finite Field in z4 of size 5^4 over its base

            sage: a = E._backend.gen()
            sage: a == x
            True
        """
        return self.gens()[0]

    def random_element(self):
        r"""
        Return a random element of this extension

        This element is constructed by converting a random element of
        the ring of this extension to the extension itself

        TESTS::

            sage: K = GF(5^2)
            sage: L = GF(5^4)
            sage: E = RingExtension(L,K); E
            Finite Field in z4 of size 5^4 over its base

            sage: x = E.random_element(); x  # random
            0
            sage: x.parent()
            Finite Field in z4 of size 5^4 over its base
        """
        elt = self._backend.random_element()
        return self.element_class(self, elt)

    def degree_over(self, base):
        base = self._check_base(base)
        return self._degree_over(base)

    def _degree_over(self, base):
        if base is self:
            return ZZ(1)
        raise NotImplementedError("degree is not implemented (and maybe not defined) for this extension")

    def degree(self, base):
        return self.degree_over(base)

    def relative_degree(self):
        return self.degree_over(self._base)

    def absolute_degree(self):
        return self.relative_degree() * self._base.absolute_degree()

    def is_finite_over(self, base=None):
        base = self._check_base(base)
        return self._is_finite_over(base)

    def _is_finite_over(self, base):
        raise NotImplementedError

    def is_free_over(self, base=None):
        base = self._check_base(base)
        return self._is_free_over(base)

    def _is_free_over(self, base):
        raise NotImplementedError

    def is_field(self, proof=False):
        return self._backend.is_field(proof=proof)

    @cached_method
    def fraction_field(self, extend_base=False, **kwargs):
        if not 'repr' in kwargs:
            kwargs['repr'] = "Fraction field of %s" % self
        if extend_base:
            from sage.rings.ring_extension_morphism import RingExtensionHomomorphism
            from sage.rings.ring_extension_morphism import backend_morphism
            defining_morphism = backend_morphism(self._backend_defining_morphism)
            defining_morphism = defining_morphism.extend_to_fraction_field()
            if isinstance(self._base, RingExtension_class):
                base = self._base.fraction_field(extend_base)
                ring = defining_morphism.codomain()
                defining_morphism = RingExtensionHomomorphism(base.Hom(ring), defining_morphism)
            return RingExtension(defining_morphism, **kwargs)
        else:
            if self.is_field():
                return self
            else:
                return RingExtension(self._backend.fraction_field(), self, **kwargs)

    def _Hom_(self, other, category):
        from sage.rings.ring_extension_homset import RingExtensionHomset
        return RingExtensionHomset(self, other, category)

    def hom(self, im_gens, codomain=None, base_map=None, check=True, category=None):
        from sage.rings.ring_extension_homset import RingExtensionHomset
        from sage.rings.ring_extension_morphism import RingExtensionHomomorphism
        if codomain is None:
            from sage.structure.sequence import Sequence
            codomain = Sequence(im_gens).universe()
        parent = self.Hom(codomain, category=category)
        return RingExtensionHomomorphism(parent, im_gens, base_map, check)

    def scalar_restriction(self, newbase):
        r"""
        Return the scalar restriction of this extension to
        the new base ``newbase``.

        INPUT:

        - ``newbase`` -- it can be either
          * a ring which coerces to the base of ``self`` or
          * an extension whose ring coerces to the base of ``self`` or
          * a morphism whose codomain coerces to the base of ``self``

        OUTPUT:

        The scalar restriction of ``self`` to the ``newbase``.

        It is defined as follows.
        In each above three cases, we get a ring homorphism ``f``
        from a domain ``domain`` to the base of ``self``.
        The scalar restriction of ``self`` to ``newbase`` is the
        extension ``self._backend/domain`` with defining morphism the
        compositum of ``f`` with the defining morphism of ``self``

        EXAMPLES::

            sage: k = GF(5^2)
            sage: K = GF(5^4)
            sage: L = GF(5^8)
            sage: E = RingExtension(L,K); E
            Finite Field in z8 of size 5^8 over its base

            sage: E1 = E.scalar_restriction(k); E1
            Finite Field in z8 of size 5^8 over its base
            sage: E1.defining_morphism()
            Ring morphism:
              From: Finite Field in z2 of size 5^2
              To:   Finite Field in z8 of size 5^8
              Defn: z2 |--> z8^7 + z8^6 + z8^5 + z8^4 + 3*z8^3 + 3*z8^2 + z8 + 3

            sage: Frob = k.frobenius_endomorphism()
            sage: E2 = E.scalar_restriction(Frob); E2
            Finite Field in z8 of size 5^8 over its base
            sage: E2.defining_morphism()
            Composite map:
              From: Finite Field in z2 of size 5^2
              To:   Finite Field in z8 of size 5^8
              Defn:   Frobenius endomorphism z2 |--> z2^5 on Finite Field in z2 of size 5^2
                    then
                      Ring morphism:
                      From: Finite Field in z2 of size 5^2
                      To:   Finite Field in z4 of size 5^4
                      Defn: z2 |--> z4^3 + z4^2 + z4 + 3
                    then
                      Ring morphism:
                      From: Finite Field in z4 of size 5^4
                      To:   Finite Field in z8 of size 5^8
                      Defn: z4 |--> 4*z8^7 + z8^6 + 3*z8^4 + z8^3 + z8^2 + 4

            sage: F = RingExtension(K, k, Frob)
            sage: E3 = E.scalar_restriction(F); E3
            Finite Field in z8 of size 5^8 over its base
            sage: E3.defining_morphism()
            Composite map:
              From: Finite Field in z2 of size 5^2
              To:   Finite Field in z8 of size 5^8
              Defn:   Frobenius endomorphism z2 |--> z2^5 on Finite Field in z2 of size 5^2
                    then
                      Ring morphism:
                      From: Finite Field in z2 of size 5^2
                      To:   Finite Field in z4 of size 5^4
                      Defn: z2 |--> z4^3 + z4^2 + z4 + 3
                    then
                      Ring morphism:
                      From: Finite Field in z4 of size 5^4
                      To:   Finite Field in z8 of size 5^8
                      Defn: z4 |--> 4*z8^7 + z8^6 + 3*z8^4 + z8^3 + z8^2 + 4

            sage: E2 is E3
            True
        """
        if isinstance(newbase, RingExtension_class):
            newbase = newbase.defining_morphism()
        elif isinstance(newbase, CommutativeRing):
            if self._base.has_coerce_map_from(newbase):
                newbase = self._base.coerce_map_from(newbase)
            else:
                raise TypeError("No coercion map from %s to %s" % (newbase, self._base))
        else:
            domain = newbase.domain()
            codomain = newbase.codomain()
            if self._base.has_coerce_map_from(codomain):
                newbase = newbase.post_compose(self._base.coerce_map_from(codomain))
            else:
                raise TypeError("No coercion map from %s to %s" % (codomain, self._base))
        defining_morphism = self._backend_defining_morphism.pre_compose(newbase)
        return RingExtension_class(defining_morphism)


# Finite free extensions
########################

class RingExtensionWithBasis(RingExtension_class):
    Element = RingExtensionWithBasisElement

    def __init__(self, defining_morphism, basis, names=None, check=True, repr=None):
        RingExtension_class.__init__(self, defining_morphism, repr=repr)
        self._basis = [ self(b) for b in basis ]
        if names is None:
            names = [ ]
            for b in self._basis:
                b = b._backend
                if b == 1:
                    names.append("")
                sb = str(b)
                if b._is_atomic() or (sb[0] == "(" and sb[-1] == ")"):
                    names.append(sb)
                else:
                    names.append("(" + sb + ")")
        else:
            if len(names) != len(self._basis):
                raise ValueError("the number of names does not match the cardinality of the basis")
        self._basis_names = names
        self._names = tuple(names)
        if check:
            try:
                _ = self.free_module(map=True)
            except (ZeroDivisionError, ArithmeticError):
                raise ValueError("the given family is not a basis")

    def _degree_over(self, base):
        if base is self:
            return ZZ(1)
        elif base is self._base:
            return len(self._basis)
        else:
            return len(self._basis) * self._base._degree_over(base)

    def _is_finite_over(self, base):
        if base is self or base is self._base:
            return True
        return self._base._is_finite_over(base)

    def _is_free_over(self, base):
        if base is self or base is self._base:
            return True
        return self._base._is_free_over(base)

    def basis_over(self, base=None):
        base = self._check_base(base)
        return self._basis_over(base)

    def _basis_over(self, base):
        if base is self:
            return [ self.one() ]
        elif base is self._base:
            return self._basis[:]
        else:
            b = self._base._basis_over(base)
            return [ x*y for y in b for x in self._basis ]

    def free_module(self, base=None, map=True):
        base = self._check_base(base)
        return self._free_module(base, map)

    # @cached_method
    def _free_module(self, base, map):
        d = self._degree_over(base)
        if map:
            from sage.rings.ring_extension_morphism import MapVectorSpaceToRelativeField, MapRelativeFieldToVectorSpace
            return base**d, MapVectorSpaceToRelativeField(self, base), MapRelativeFieldToVectorSpace(self, base)
        else:
            return base**d

    @cached_method
    def fraction_field(self, extend_base=False, **kwargs):
        if extend_base and len(kwargs) == 0:
            kwargs['basis'] = self._basis
            kwargs['names'] = self._basis_names
        return RingExtension_class.fraction_field(self, extend_base, **kwargs)


class RingExtensionWithGen(RingExtensionWithBasis):
    def __init__(self, defining_morphism, gen, name, check=True, repr=None):
        self._name = str(name)
        backend_base = defining_morphism.domain()
        if isinstance(backend_base, RingExtension_class):
            backend_base = backend_base._backend
        _, deg_domain, deg_codomain = _common_base(backend_base, defining_morphism.codomain(), degree=True)
        degree = deg_codomain // deg_domain
        names = [ "" ]
        if degree == 1:
            self._name = None
        else:
            names += [ self._name ] + [ "%s^%s" % (self._name, i) for i in range(2,degree) ]
        basis = [ gen ** i for i in range(degree) ]
        RingExtensionWithBasis.__init__(self, defining_morphism, basis, names, check, repr=repr)
        self._gen = self._backend(gen)
        self._names = (self._name,)
        self._type = "Ring"
        if self._backend in Fields():
            self._type = "Field"

    def _repr_(self, base=None):
        if self._name is None:
            return "%s over its base" % self._backend
        else:
            return "%s in %s with defining polynomial %s over its base" % (self._type, self._name, self.modulus())

    def modulus(self, var='x'):
        from sage.rings.polynomial.polynomial_ring_constructor import PolynomialRing
        _, _, j = self.free_module(map=True)
        d = self.relative_degree()
        coeffs = [ -c for c in j(self._gen**d) ] + [ 1 ]
        S = PolynomialRing(self._base, name=var)
        return S(coeffs)

    def gens(self, base=None):
        if base is None:
            return (self(self._gen),)
        base = self._check_base(base)
        gens = tuple([])
        b = self
        while b is not base:
            gens += b.gens()
            b = b.base()
        return gens

    @cached_method
    def fraction_field(self, extend_base=False, **kwargs):
        if extend_base and len(kwargs) == 0:
            kwargs['gen'] = self._gen
            kwargs['name'] = self._name
        return RingExtension_class.fraction_field(self, extend_base, **kwargs)
