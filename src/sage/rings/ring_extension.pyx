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
from sage.cpython.getattr cimport AttributeErrorMessage

from sage.structure.factory import UniqueFactory
from sage.structure.parent cimport Parent
from sage.structure.element cimport Element
from sage.structure.category_object import normalize_names
from sage.categories.map cimport Map
from sage.categories.commutative_rings import CommutativeRings
from sage.categories.commutative_algebras import CommutativeAlgebras
from sage.categories.fields import Fields
from sage.rings.ring cimport CommutativeRing, CommutativeAlgebra
from sage.rings.integer_ring import ZZ
from sage.rings.infinity import Infinity

from sage.rings.ring_extension_element cimport RingExtensionElement
from sage.rings.ring_extension_element cimport RingExtensionFractionFieldElement
from sage.rings.ring_extension_element cimport RingExtensionWithBasisElement
from sage.rings.ring_extension_element cimport backend_element, from_backend_element
from sage.rings.ring_extension_morphism cimport RingExtensionHomomorphism
from sage.rings.ring_extension_morphism cimport RingExtensionBackendIsomorphism
from sage.rings.ring_extension_morphism cimport RingExtensionBackendReverseIsomorphism
from sage.rings.ring_extension_morphism cimport backend_morphism, are_equal_morphisms
from sage.rings.ring_extension_morphism cimport MapVectorSpaceToRelativeField, MapRelativeFieldToVectorSpace



# Helper functions
##################

cdef _common_base(K, L, degree):
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


cpdef backend_parent(R):
    if isinstance(R, RingExtension_class):
        return (<RingExtension_class>R)._backend
    else:
        return R

cpdef from_backend_parent(R, RingExtension_class E):
    cdef CommutativeRing base = E
    while isinstance(base, RingExtension_class):
        if (<RingExtension_class>base)._backend is R:
            return base
        base = base._base


cdef to_backend(arg):
    if isinstance(arg, list):
        return [ to_backend(x) for x in arg ]
    elif isinstance(arg, tuple):
        return ( to_backend(x) for x in arg )
    elif isinstance(arg, dict):
        return { to_backend(key): to_backend(value) for (key, value) in arg.items() }
    elif isinstance(arg, RingExtension_class):
        return (<RingExtension_class>arg)._backend
    elif isinstance(arg, RingExtensionElement):
        return (<RingExtensionElement>arg)._backend
    return arg

cdef from_backend(arg, RingExtension_class E):
    cdef ans = None
    if isinstance(arg, list):
        ans = [ from_backend(x,E) for x in arg ]
    elif isinstance(arg, tuple):
        ans = tuple(from_backend(x,E) for x in arg)
    elif isinstance(arg, dict):
        ans = { from_backend(key,E): from_backend(value,E) for (key, value) in arg.items() }
    elif isinstance(arg, Parent):
        ans = from_backend_parent(arg,E)
    elif isinstance(arg, Element):
        ans = from_backend_element(arg,E)
    if ans is None:
        return arg
    else:
        return ans


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
    def create_key_and_extra_args(self, constructors, defining_morphism, gen=None, names=None):
        key = (defining_morphism, gen, names)
        return key, { 'constructors': constructors }

    def create_object(self, version, key, **extra_args):
        (defining_morphism, gen, names) = key
        constructors = extra_args['constructors']
        for (constructor, kwargs) in constructors[:-1]:
            try:
                return constructor(defining_morphism, **kwargs)
            except (NotImplementedError, ValueError, TypeError):
                pass
        (constructor, kwargs) = constructors[-1]
        return constructor(defining_morphism, **kwargs)

_RingExtension = RingExtensionFactory("RingExtension")


def RingExtension(*args, gen=None, gens=None, name=None, names=None):
    from sage.categories.map import Map
    from sage.rings.ring_extension_morphism import RingExtensionHomomorphism

    # We parse the arguments
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
        base = None
        ring = rings[0]
    elif len(rings) == 2:
        ring, base = rings
    else:
        raise TypeError("you must specify at least a ring or a defining morphism")

    # We compute the defining morphism
    print_parent_as = None
    print_elements_as = None
    if defining_morphism is None:
        if base is None:
            base = ring.base_ring()
        if isinstance(base, RingExtension_class):
            backend_base = (<RingExtension_class>base)._backend
            if ring.has_coerce_map_from(backend_base):
                defining_morphism = RingExtensionHomomorphism(base.Hom(ring), ring.coerce_map_from(backend_base))
        else:
            if ring.has_coerce_map_from(base):
                defining_morphism = ring.coerce_map_from(base)
        if defining_morphism is None:
            raise ValueError("No coercion map from %s to %s" % (base,ring))
    else:
        if base is None:
            base = defining_morphism.domain()
        elif defining_morphism.domain() is not base:
            defining_morphism = defining_morphism.extend_domain(base)
        if defining_morphism.codomain() is not ring:
            defining_morphism = defining_morphism.extend_codomain(ring)
    if isinstance(ring, RingExtension_class):
        defining_morphism = backend_morphism(defining_morphism, forget="codomain")
        # print_parent_as = ...
        print_elements_as = ring
        ring = (<RingExtension_class>ring)._backend

    if names is None and name is not None:
        names = name
    if gens is not None:
        if gen is not None:
            raise ValueError("you cannot specify gen and gens simultaneously")
        if not isinstance(gens, (list, tuple)):
            raise TypeError("gens must be a list or a tuple")
        if len(gens) > 1:
            raise NotImplementedError("only ring extensions with a single generator are implemented")
        gen = ring(gens[0])
    if gen is not None:
        if names is None:
            raise TypeError("you must specify the name of the generator")
        names = normalize_names(1, names)
        constructors = [ (RingExtensionWithGen, {'gen': gen, 'names': names}) ]
    else:
        kwargs = { 'print_parent_as': print_parent_as, 'print_elements_as': print_elements_as }
        constructors = [ (RingExtension_class, kwargs) ]
        # we check if there is a canonical generator
        # if so, we first try the constructor RingExtensionWithGen
        if ring.ngens() == 1:
            gen = ring.gen()
            if names is None:
                names = ring.variable_names()
            elif isinstance(names, (list, tuple)):
                if len(names) != 1:
                    raise ValueError("the number of names does not match the number of generators")
                names = tuple(names)
            else:
                names = (names,)
            constructors = [ (RingExtensionWithGen, {'gen': gen, 'names': names}) ] + constructors

    return _RingExtension(constructors, defining_morphism, gen, names)


# General extensions
####################

cdef class RingExtension_class(CommutativeAlgebra):
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
    # Apparently the method __getattr__ clashes with the coercion system
    # It's why we set _element_constructor_ here
    # (Is there a better solution to fix this issue?)
    Element = _element_constructor_ = RingExtensionElement

    def __init__(self, defining_morphism, print_parent_as=None, print_elements_as=None):
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
        cdef CommutativeRing base, ring
        cdef CommutativeRing b, backend
        cdef Map f

        base = defining_morphism.domain()
        ring = defining_morphism.codomain()
        CommutativeAlgebra.__init__(self, ZZ, category=CommutativeAlgebras(base))
        self._base = base
        self._backend = ring
        self._backend_defining_morphism = defining_morphism
        self._defining_morphism = RingExtensionHomomorphism(self._base.Hom(self), defining_morphism)
        self._print_parent_as = print_parent_as
        self._print_elements_as = print_elements_as
        self._type = "Ring"
        if self._backend in Fields():
            self._type = "Field"

        # Some checkings
        if (base not in CommutativeRings()
         or ring not in CommutativeRings()
         or not defining_morphism.category_for().is_subcategory(CommutativeRings())):
            raise TypeError("only commutative rings are accepted")
        f = ring.Hom(ring).identity()
        b = self
        while isinstance(b, RingExtension_class):
            f *= backend_morphism((<RingExtension_class>b)._backend_defining_morphism)
            b = b._base
            if isinstance(b, RingExtension_class):
                backend = (<RingExtension_class>b)._backend
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

    def __getattr__(self, name):
        method = None
        if hasattr(self._backend, name):
            method = getattr(self._backend, name)
        if not callable(method):
            raise AttributeError(AttributeErrorMessage(self, name))
        def wrapper(*args, **kwargs):
            output = method(*to_backend(args), **to_backend(kwargs))
            return from_backend(output, self)
        wrapper.__doc__ = method.__doc__
        return wrapper

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

    def __repr__(self):
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
        print_as = self._print_parent_as
        if print_as is not None:
            return str(print_as)
        return self._repr_()

    def _repr_(self):
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
        cdef RingExtension_class right
        if isinstance(other, RingExtension_class):
            right = <RingExtension_class>other
            if self._backend.has_coerce_map_from(right._backend) and self._base.has_coerce_map_from(right._base):
                backend = self._backend.coerce_map_from(right._backend)
                f = backend * backend_morphism(right._defining_morphism)
                g = backend_morphism(self._defining_morphism * self._base.coerce_map_from(right._base))
                if are_equal_morphisms(f, g):
                    return RingExtensionHomomorphism(right.Hom(self), backend)

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

    cpdef is_defined_over(self, base):
        cdef CommutativeRing b
        b = self
        while isinstance(b, RingExtension_class):
            if b is base: return True
            b = (<RingExtension_class>b)._base
        return b == base

    cpdef CommutativeRing _check_base(self, CommutativeRing base):
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

    cpdef _degree_over(self, CommutativeRing base):
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

    cpdef _is_finite_over(self, CommutativeRing base):
        raise NotImplementedError

    def is_free_over(self, base=None):
        base = self._check_base(base)
        return self._is_free_over(base)

    cpdef _is_free_over(self, CommutativeRing base):
        raise NotImplementedError

    def is_field(self, proof=False):
        return self._backend.is_field(proof=proof)

    @cached_method
    def fraction_field(self, extend_base=False):
        defining_morphism = self._defining_morphism_fraction_field(extend_base)
        if defining_morphism is None:
            return self
        if extend_base:
            return _RingExtension([(RingExtensionFractionField, {'ring': self})], defining_morphism)
        else:
            return _RingExtension([(RingExtensionFractionField, {})], defining_morphism)

    cdef Map _defining_morphism_fraction_field(self, bint extend_base):
        if extend_base:
            defining_morphism = backend_morphism(self._backend_defining_morphism)
            defining_morphism = defining_morphism.extend_to_fraction_field()
            if isinstance(self._base, RingExtension_class):
                base = self._base.fraction_field(extend_base)
                ring = defining_morphism.codomain()
                defining_morphism = RingExtensionHomomorphism(base.Hom(ring), defining_morphism)
        else:
            if self.is_field():
                defining_morphism = None
            else:
                ring = self._backend.fraction_field()
                defining_morphism = RingExtensionHomomorphism(self.Hom(ring), ring.coerce_map_from(self._backend))
        return defining_morphism

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


# Fraction fields
#################

cdef class RingExtensionFractionField(RingExtension_class):
    Element = _element_constructor_ = RingExtensionFractionFieldElement

    def __init__(self, defining_morphism, ring=None, **kwargs):
        RingExtension_class.__init__(self, defining_morphism, **kwargs)
        if ring is None:
            self._ring = self._base
        else:
            self._ring = ring

    def ring(self):
        return self._ring

    def _repr_(self):
        if self._ring in Fields():
            s = str(self._ring)
        else:
            s = "Fraction field of %s" % self._ring
        if not isinstance(self._ring, RingExtension_class):
            s += " over its base"
        return s


# Finite free extensions
########################

cdef class RingExtensionWithBasis(RingExtension_class):
    Element = _element_constructor_ = RingExtensionWithBasisElement

    def __init__(self, defining_morphism, basis, names=None, check=True, **kwargs):
        RingExtension_class.__init__(self, defining_morphism, **kwargs)
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

    cpdef _degree_over(self, CommutativeRing base):
        if base is self:
            return ZZ(1)
        elif base is self._base:
            return len(self._basis)
        else:
            return len(self._basis) * self._base._degree_over(base)

    cpdef _is_finite_over(self, CommutativeRing base):
        if base is self or base is self._base:
            return True
        return self._base._is_finite_over(base)

    cpdef _is_free_over(self, CommutativeRing base):
        if base is self or base is self._base:
            return True
        return self._base._is_free_over(base)

    def basis_over(self, base=None):
        base = self._check_base(base)
        return self._basis_over(base)

    cpdef _basis_over(self, CommutativeRing base):
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

    @cached_method
    def _free_module(self, base, map):
        d = self._degree_over(base)
        if map:
            return base**d, MapVectorSpaceToRelativeField(self, base), MapRelativeFieldToVectorSpace(self, base)
        else:
            return base**d

    @cached_method
    def fraction_field(self, extend_base=False):
        defining_morphism = self._defining_morphism_fraction_field(extend_base)
        if defining_morphism is None:
            return self
        if extend_base:
            basis = self._basis
            names = self._basis_names
            constructor = RingExtensionWithBasis
            kwargs = { 'basis': basis, 'names': names, 'check': False }
        else:
            gen = names = None
            constructor = RingExtensionFractionField
            kwargs = { 'print_elements_as': self.fraction_field(extend_base=True) }
        return _RingExtension([(constructor, kwargs)], defining_morphism, gen, names)


cdef class RingExtensionWithGen(RingExtensionWithBasis):
    def __init__(self, defining_morphism, gen, names, check=True, **kwargs):
        self._name = names[0]
        backend_base = backend_parent(defining_morphism.domain())
        _, deg_domain, deg_codomain = _common_base(backend_base, defining_morphism.codomain(), True)
        degree = deg_codomain // deg_domain
        basis_names = [ "" ]
        if degree == 1:
            self._name = None
        else:
            basis_names += [ self._name ] + [ "%s^%s" % (self._name, i) for i in range(2,degree) ]
        basis = [ gen ** i for i in range(degree) ]
        RingExtensionWithBasis.__init__(self, defining_morphism, basis, basis_names, check, **kwargs)
        self._gen = self._backend(gen)
        self._names = (self._name,)

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
    def fraction_field(self, extend_base=False):
        defining_morphism = self._defining_morphism_fraction_field(extend_base)
        if defining_morphism is None:
            return self
        if extend_base:
            gen = self._gen
            names = self._names
            constructor = RingExtensionWithGen
            kwargs = { 'gen': gen, 'names': names, 'check': False }
        else:
            gen = names = None
            constructor = RingExtensionFractionField
            kwargs = { 'print_elements_as': self.fraction_field(extend_base=True) }
        return _RingExtension([(constructor, kwargs)], defining_morphism, gen, names)
