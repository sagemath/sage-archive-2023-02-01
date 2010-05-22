#*****************************************************************************
#       Copyright (C) 2008 Robert Bradshaw <robertwb@math.washington.edu>
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

include "../ext/stdsage.pxi"

import homset

from sage.structure.element import generic_power
from sage.structure.parent import Set_PythonType

def unpickle_map(_class, parent, _dict, _slots):
    """
    Auxiliary function for unpickling a map.

    TEST::

        sage: R.<x,y> = QQ[]
        sage: f = R.hom([x+y,x-y],R)
        sage: f == loads(dumps(f))  # indirect doctest
        True

    """
    # should we use slots?
    # from element.pyx
    cdef Map mor = _class.__new__(_class)
    mor._set_parent(parent)
    mor._update_slots(_slots)
    if HAS_DICTIONARY(mor):
        mor.__dict__ = _dict
    return mor

def is_Map(x):
    """
    Auxiliary function: Is the argument a map?

    EXAMPLE::

        sage: R.<x,y> = QQ[]
        sage: f = R.hom([x+y,x-y],R)
        sage: from sage.categories.map import is_Map
        sage: is_Map(f)
        True
    """
    return isinstance(x, Map)

cdef class Map(Element):
    """
    Basic class for all maps.

    NOTE:

    The call method is of course not implemented in this base class. This must
    be done in the sub classes, by overloading :meth:`_call_` and possibly also
    :meth:`_call_with_args`.

    EXAMPLES:

    Usually, instances of this class will not be constructed directly, but
    for example like this::

        sage: from sage.categories.morphism import SetMorphism
        sage: X.<x> = ZZ[]
        sage: Y = ZZ
        sage: phi = SetMorphism(Hom(X, Y, Rings()), lambda p: p[0])
        sage: phi(x^2+2*x-1)
        -1
        sage: R.<x,y> = QQ[]
        sage: f = R.hom([x+y,x-y],R)
        sage: f(x^2+2*x-1)
        x^2 + 2*x*y + y^2 + 2*x + 2*y - 1
    """
    def __init__(self, parent, codomain=None):
        """
        INPUT:

        There can be one or two arguments of this init method. If it is one argument,
        it must be a hom space. If it is two arguments, it must be two parent structures
        that will be domain and codomain of the map-to-be-created.

        TESTS::

            sage: from sage.categories.map import Map

        Using a hom space::

            sage: Map(Hom(QQ, ZZ, Rings()))
            Generic map:
              From: Rational Field
              To:   Integer Ring

        Using domain and codomain::

            sage: Map(QQ['x'], SymmetricGroup(6))
            Generic map:
              From: Univariate Polynomial Ring in x over Rational Field
              To:   SymmetricGroup(6)
        """
        if codomain is not None:
            if PY_TYPE_CHECK(parent, type):
                parent = Set_PythonType(parent)
            parent = homset.Hom(parent, codomain)
        elif not isinstance(parent, homset.Homset):
            raise TypeError, "parent (=%s) must be a Homspace"%parent
        Element.__init__(self, parent)
        self._domain = parent.domain()
        self._codomain = parent.codomain()
        if self._domain.is_exact() and self._codomain.is_exact():
            self._coerce_cost = 10 # default value.
        else:
            self._coerce_cost = 10000 # inexact morphisms are bad.

    cdef _update_slots(self, _slots):
        """
        Auxiliary method, used in pickling.

        INPUT:

        - A dictionary of slots to be updated.
          The dictionary must have the keys ``'_domain'`` and ``'_codomain'``,
          and may have the key ``'_repr_type_str'``.

        TESTS:

        Since it is a ``cdef``d method, it is tested using a dummy python method.
        ::

            sage: from sage.categories.map import Map
            sage: f = Map(Hom(QQ, ZZ, Rings()))
            sage: f._update_slots_test({"_domain": RR, "_codomain": QQ}) # indirect doctest
            sage: f.domain()
            Real Field with 53 bits of precision
            sage: f.codomain()
            Rational Field
            sage: f._repr_type_str
            sage: f._update_slots_test({"_repr_type_str": "bla", "_domain": RR, "_codomain": QQ})
            sage: f._repr_type_str
            'bla'
        """
        self._domain = _slots['_domain']
        self._codomain = _slots['_codomain']
        # Several pickles exist without a _repr_type_str, so
        # if there is none saved, we just set it to None.
        if _slots.has_key('_repr_type_str'):
            self._repr_type_str = _slots['_repr_type_str']
        else:
            self._repr_type_str = None

    def _update_slots_test(self, _slots):
        """
        A Python method to test the cdef _update_slots method.

        TESTS::

            sage: from sage.categories.map import Map
            sage: f = Map(Hom(QQ, ZZ, Rings()))
            sage: f._update_slots_test({"_domain": RR, "_codomain": QQ})
            sage: f.domain()
            Real Field with 53 bits of precision
            sage: f.codomain()
            Rational Field
            sage: f._repr_type_str
            sage: f._update_slots_test({"_repr_type_str": "bla", "_domain": RR, "_codomain": QQ})
            sage: f._repr_type_str
            'bla'
        """
        self._update_slots(_slots)

    cdef _extra_slots(self, _slots):
        """
        A Python method to test the cdef _extra_slots method.

        TESTS::

            sage: from sage.categories.map import Map
            sage: f = Map(Hom(QQ, ZZ, Rings()))
            sage: f._extra_slots_test({"bla": 1})
            {'_codomain': Integer Ring, '_domain': Rational Field, 'bla': 1, '_repr_type_str': None}
        """
        _slots['_domain'] = self._domain
        _slots['_codomain'] = self._codomain
        _slots['_repr_type_str'] = self._repr_type_str
        return _slots

    def _extra_slots_test(self, _slots):
        """
        Auxiliary method, used in pickling.

        INPUT:

        A dictionary.

        OUTPUT:

        The given dictionary, that is updated by the slots '_domain', '_codomain' and '_repr_type_str'.


        TESTS::

            sage: from sage.categories.map import Map
            sage: f = Map(Hom(QQ, ZZ, Rings()))
            sage: f._extra_slots_test({"bla": 1}) # indirect doctest
            {'_codomain': Integer Ring, '_domain': Rational Field, 'bla': 1, '_repr_type_str': None}
        """
        return self._extra_slots(_slots)

    def __reduce__(self):
        """
        TESTS::

            sage: from sage.categories.map import Map
            sage: f = Map(Hom(QQ, ZZ, Rings())); f
            Generic map:
              From: Rational Field
              To:   Integer Ring
            sage: loads(dumps(f))  # indirect doctest
            Generic map:
              From: Rational Field
              To:   Integer Ring

        """
        if HAS_DICTIONARY(self):
            _dict = self.__dict__
        else:
            _dict = {}
        return unpickle_map, (self.__class__, self._parent, _dict, self._extra_slots({}))

    def _repr_type(self):
        """
        Return a string describing the specific type of this map, to be used when printing ``self``.

        NOTE:

        By default, the string ``"Generic"`` is returned. Subclasses may overload this method.

        EXAMPLE::

            sage: from sage.categories.map import Map
            sage: f = Map(Hom(QQ, ZZ, Rings()))
            sage: print f._repr_type()
            Generic
            sage: R.<x,y> = QQ[]
            sage: phi = R.hom([x+y,x-y],R)
            sage: print phi._repr_type()
            Ring

        """
        if self._repr_type_str is None:
            return "Generic"
        else:
            return self._repr_type_str

    def _repr_defn(self):
        """
        Return a string describing the definition of ``self``, to be used when printing ``self``.

        NOTE:

        By default, the empty string is returned. Subclasses may overload this method.

        EXAMPLE::

            sage: from sage.categories.map import Map
            sage: f = Map(Hom(QQ, ZZ, Rings()))
            sage: f._repr_defn() == ''
            True
            sage: R.<x,y> = QQ[]
            sage: f = R.hom([x+y,x-y],R)
            sage: print f._repr_defn()
            x |--> x + y
            y |--> x - y

        """
        return ""

    def _repr_(self):
        """
        NOTE:

        The string representation is based on the strings returned by
        :meth:`_repr_defn` and :meth:`_repr_type`, as well as the domain
        and the codomain.

        TEST::

            sage: from sage.categories.map import Map
            sage: Map(Hom(QQ, ZZ, Rings()))    # indirect doctest
            Generic map:
              From: Rational Field
              To:   Integer Ring
            sage: R.<x,y> = QQ[]
            sage: R.hom([x+y,x-y],R)
            Ring endomorphism of Multivariate Polynomial Ring in x, y over Rational Field
              Defn: x |--> x + y
                    y |--> x - y
        """
        s = "%s map:"%self._repr_type()
        s += "\n  From: %s"%self.domain()
        s += "\n  To:   %s"%self.codomain()
        d = self._repr_defn()
        if d != '':
            s += "\n  Defn: %s"%('\n        '.join(self._repr_defn().split('\n')))
        return s

    cpdef domain(self):
        """
        Return the domain of ``self``.

        EXAMPLE::

            sage: from sage.categories.map import Map
            sage: f = Map(Hom(QQ, ZZ, Rings()))
            sage: f.domain()
            Rational Field
            sage: R.<x,y> = QQ[]
            sage: phi = R.hom([x+y,x-y],R)
            sage: phi.domain()
            Multivariate Polynomial Ring in x, y over Rational Field
        """
        return self._domain

    cpdef codomain(self):
        """
        Return the codomain of ``self``.

        EXAMPLE::

            sage: from sage.categories.map import Map
            sage: f = Map(Hom(QQ, ZZ, Rings()))
            sage: f.codomain()
            Integer Ring
            sage: R.<x,y> = QQ[]
            sage: phi = R.hom([x+y,x-y],R)
            sage: phi.codomain()
            Multivariate Polynomial Ring in x, y over Rational Field
        """
        return self._codomain

    def category_for(self):
        """
        Returns the category self is a morphism for.

        NOTE:

        This is different from the category of maps to which this
        map belongs *as an object*.

        EXAMPLES::

            sage: from sage.categories.morphism import SetMorphism
            sage: X.<x> = ZZ[]
            sage: Y = ZZ
            sage: phi = SetMorphism(Hom(X, Y, Rings()), lambda p: p[0])
            sage: phi.category_for()
            Category of rings
            sage: phi.category()
            Category of hom sets in Category of rings
            sage: R.<x,y> = QQ[]
            sage: f = R.hom([x+y,x-y],R)
            sage: f.category_for()
            Category of commutative rings
            sage: f.category()
            Category of hom sets in Category of rings

        FIXME: find a better name for this method
        """
        return self.parent().homset_category()

    def __call__(self, x, *args, **kwds):
        """
        Apply this map to ``x``.

        IMPLEMENTATION:

        - To implement the call method in a subclass of Map, implement
          :meth:`_call_` and possibly also :meth:`_call_with_args` and
          :meth:`pushout`.
        - If ``x`` can be coerced into the domain of ``self``, then the
          method ``_call_`` (single underscore) is called after coercion.
        - Otherwise, it is tried to call the method ``pushout`` to ``x``
          (which may be provided in subclasses); in that way, ``self``
          could be applied to objects like ideals.
        - If ``self`` is called with additional arguments (or keyword
          arguments) then :meth:`_call_with_args` is called.

        INPUT:

            ``x`` -- an element coercible to the domain of ``self``; also
                     objects like ideals are supported in some cases

        OUTPUT:

            an element (or ideal, etc.)

        EXAMPLES::

            sage: R.<x,y> = QQ[]; phi=R.hom([y,x])
            sage: phi(y)          # indirect doctest
            x

        We take the image of an ideal::

            sage: I = ideal(x,y); I
            Ideal (x, y) of Multivariate Polynomial Ring in x, y over Rational Field
            sage: phi(I)
            Ideal (y, x) of Multivariate Polynomial Ring in x, y over Rational Field

        TEST::

            sage: from sage.categories.map import Map
            sage: f = Map(Hom(ZZ, QQ, Rings()))
            sage: f(1/2)
            Traceback (most recent call last):
            ...
            AttributeError: 'sage.categories.map.Map' object has no attribute 'pushforward'

        """
        if len(args) == 0 and len(kwds) == 0:
            if not PY_TYPE_CHECK(x, Element):
                return self._call_(x)
            elif (<Element>x)._parent is not self._domain:
                try:
                    x = self._domain(x)
                except TypeError:
                    try:
                        return self.pushforward(x)
                    except (TypeError, NotImplementedError):
                        raise TypeError, "%s must be coercible into %s"%(x, self._domain)
            return self._call_(x)
        else:
            if PY_TYPE_CHECK(x, Element):
                if (<Element>x)._parent is not self._domain:
                    x = self._domain(x)
            return self._call_with_args(x, args, kwds)

    cpdef Element _call_(self, x):
        """
        Call method with a single argument, not implemented in the base class.

        TEST::

            sage: from sage.categories.map import Map
            sage: f = Map(Hom(QQ, ZZ, Rings()))
            sage: f(1/2)             # indirect doctest
            Traceback (most recent call last):
            ...
            NotImplementedError: <type 'sage.categories.map.Map'>
        """
        raise NotImplementedError, type(self)

    cpdef Element _call_with_args(self, x, args=(), kwds={}):
        """
        Call method with multiple arguments, not implemented in the base class.

        TEST::

            sage: from sage.categories.map import Map
            sage: f = Map(Hom(QQ, ZZ, Rings()))
            sage: f(1/2,2,foo='bar')        # indirect doctest
            Traceback (most recent call last):
            ...
            NotImplementedError: _call_with_args not overridden to accept arguments for <type 'sage.categories.map.Map'>
        """
        if len(args) == 0 and len(kwds) == 0:
            return self(x)
        else:
            raise NotImplementedError, "_call_with_args not overridden to accept arguments for %s" % type(self)

    def __mul__(self, right):
        r"""
        The multiplication * operator is operator composition

        IMPLEMENTATION:

        If you want to change the behaviour of composition for
        derived classes, please overload :meth:`_composition_`
        (but not :meth:`_composition`!) of the left factor.

        INPUT:

        - ``self`` -- Map
        - ``right`` -- Map

        OUTPUT:

        The map $x \mapsto self(right(x))$.

        EXAMPLES::

            sage: from sage.categories.morphism import SetMorphism
            sage: X.<x> = ZZ[]
            sage: Y = ZZ
            sage: Z = QQ
            sage: phi_xy = SetMorphism(Hom(X, Y, Rings()), lambda p: p[0])
            sage: phi_yz = SetMorphism(Hom(Y, Z, CommutativeAdditiveMonoids()), lambda y: QQ(y)/2)
            sage: phi_yz * phi_xy
            Composite map:
              From: Univariate Polynomial Ring in x over Integer Ring
              To:   Rational Field
              Defn:   Generic morphism:
                      From: Univariate Polynomial Ring in x over Integer Ring
                      To:   Integer Ring
                    then
                      Generic morphism:
                      From: Integer Ring
                      To:   Rational Field

        If ``right`` is a ring homomorphism given by the images of
        generators, then it is attempted to form the composition
        accordingly. Only if this fails, or if the result does not
        belong to the given homset, a formal composite map is
        returned (as above).
        ::

            sage: R.<x,y> = QQ[]
            sage: S.<a,b> = QQ[]
            sage: f = R.hom([x+y,x-y],R)
            sage: f = R.hom([a+b,a-b])
            sage: g = S.hom([x+y,x-y])
            sage: f*g
            Ring endomorphism of Multivariate Polynomial Ring in a, b over Rational Field
              Defn: a |--> 2*a
                    b |--> 2*b
            sage: h = SetMorphism(Hom(S,QQ,Rings()), lambda p: p.lc())
            sage: h*f
            Composite map:
              From: Multivariate Polynomial Ring in x, y over Rational Field
              To:   Rational Field
              Defn:   Ring morphism:
                      From: Multivariate Polynomial Ring in x, y over Rational Field
                      To:   Multivariate Polynomial Ring in a, b over Rational Field
                      Defn: x |--> a + b
                            y |--> a - b
                    then
                      Generic morphism:
                      From: Multivariate Polynomial Ring in a, b over Rational Field
                      To:   Rational Field
        """
        if not isinstance(right, Map):
            raise TypeError, "right (=%s) must be a map to multiply it by %s"%(right, self)
        if right.codomain() != self.domain():
            raise TypeError, "self (=%s) domain must equal right (=%s) codomain"%(self, right)
        return self._composition(right)

    def _composition(self, right):
        """
        Composition of maps, which generically returns a :class:`CompositeMap`.

        INPUT:

        - ``self``  -- a Map in some ``Hom(Y, Z, category_left)``
        - ``right`` -- a Map in some ``Hom(X, Y, category_right)``

        OUTPUT:

        Returns the composition of ``self`` and ``right`` as a
        morphism in ``Hom(X, Z, category)`` where ``category`` is the
        meet of ``category_left`` and ``category_right``.

        EXAMPLES::

            sage: from sage.categories.morphism import SetMorphism
            sage: X.<x> = ZZ[]
            sage: Y = ZZ
            sage: Z = QQ
            sage: phi_xy = SetMorphism(Hom(X, Y, Rings()), lambda p: p[0])
            sage: phi_yz = SetMorphism(Hom(Y, Z, CommutativeAdditiveMonoids()), lambda y: QQ(y)/2)
            sage: phi_yz._composition(phi_xy)
            Composite map:
              From: Univariate Polynomial Ring in x over Integer Ring
              To:   Rational Field
              Defn:   Generic morphism:
                      From: Univariate Polynomial Ring in x over Integer Ring
                      To:   Integer Ring
                    then
                      Generic morphism:
                      From: Integer Ring
                      To:   Rational Field
            sage: phi_yz.category_for()
            Category of commutative additive monoids
        """
        category = self.category_for()._meet_(right.category_for())
        H = homset.Hom(right.domain(), self.codomain(), category)
        return self._composition_(right, H)

    def _composition_(self, right, homset):
        """
        INPUT:

        - ``self``, ``right`` -- maps
        - homset -- a homset

        ASSUMPTION:

        The codomain of ``right`` is contained in the domain of ``self``.
        This assumption is not verified.

        OUTPUT:

        Returns a formal composite map, the composition of ``right``
        followed by ``self``, as a morphism in ``homset``.

        Classes deriving from :class:`Map` are encouraged to override
        this whenever meaningful. This is the case, e.g., for ring
        homomorphisms.

        EXAMPLES::

            sage: Rx.<x> = ZZ['x']
            sage: Ry.<y> = ZZ['y']
            sage: Rz.<z> = ZZ['z']
            sage: phi_xy = Rx.hom([y+1]); phi_xy
            Ring morphism:
              From: Univariate Polynomial Ring in x over Integer Ring
              To:   Univariate Polynomial Ring in y over Integer Ring
              Defn: x |--> y + 1
            sage: phi_yz = Ry.hom([z+1]); phi_yz
            Ring morphism:
              From: Univariate Polynomial Ring in y over Integer Ring
              To:   Univariate Polynomial Ring in z over Integer Ring
              Defn: y |--> z + 1
            sage: phi_xz = phi_yz._composition_(phi_xy, Hom(Rx, Rz, Monoids()))
            sage: phi_xz
            Composite map:
              From: Univariate Polynomial Ring in x over Integer Ring
              To:   Univariate Polynomial Ring in z over Integer Ring
              Defn:   Ring morphism:
                      From: Univariate Polynomial Ring in x over Integer Ring
                      To:   Univariate Polynomial Ring in y over Integer Ring
                      Defn: x |--> y + 1
                    then
                      Ring morphism:
                      From: Univariate Polynomial Ring in y over Integer Ring
                      To:   Univariate Polynomial Ring in z over Integer Ring
                      Defn: y |--> z + 1
            sage: phi_xz.category_for()
            Category of monoids

        TEST:

        This illustrates that it is not tested whether the maps can actually
        be composed, i.e., whether codomain and domain match.
        ::

            sage: R.<x,y> = QQ[]
            sage: S.<a,b> = QQ[]
            sage: f_R = R.hom([x+y,x-y],R)
            sage: f_S = S.hom([a+b,a-b],S)
            sage: foo_bar = f_R._composition_(f_S, Hom(S, R, Monoids()))
            sage: foo_bar(a)
            2*x

        However, it is tested when attempting to compose the maps in
        the usual multiplicative notation::

            sage: f_R*f_S
            Traceback (most recent call last):
            ...
            TypeError: self (=Ring endomorphism of Multivariate Polynomial Ring in x, y over Rational Field
              Defn: x |--> x + y
                    y |--> x - y) domain must equal right (=Ring endomorphism of Multivariate Polynomial Ring in a, b over Rational Field
              Defn: a |--> a + b
                    b |--> a - b) codomain
        """
        return FormalCompositeMap(homset, right, self)

    def pre_compose(self, right):
        """
        INPUT:
         - ``self`` -- a Map in some ``Hom(Y, Z, category_left)``
         - ``left`` -- a Map in some ``Hom(X, Y, category_right)``

        Returns the composition of ``right`` followed by ``self`` as a
        morphism in ``Hom(X, Z, category)`` where ``category`` is the
        meet of ``category_left`` and ``category_right``.

        EXAMPLES::

            sage: from sage.categories.morphism import SetMorphism
            sage: X.<x> = ZZ[]
            sage: Y = ZZ
            sage: Z = QQ
            sage: phi_xy = SetMorphism(Hom(X, Y, Rings()), lambda p: p[0])
            sage: phi_yz = SetMorphism(Hom(Y, Z, Monoids()), lambda y: QQ(y**2))
            sage: phi_xz = phi_yz.pre_compose(phi_xy); phi_xz
            Composite map:
              From: Univariate Polynomial Ring in x over Integer Ring
              To:   Rational Field
              Defn:   Generic morphism:
                      From: Univariate Polynomial Ring in x over Integer Ring
                      To:   Integer Ring
                    then
                      Generic morphism:
                      From: Integer Ring
                      To:   Rational Field
            sage: phi_xz.category_for()
            Category of monoids
        """
        if self.domain() is not right.codomain():
            right = right.extend_codomain(self.domain())
        return self._composition(right)

    def post_compose(self, left):
        """
        INPUT:
         - ``self`` -- a Map in some ``Hom(X, Y, category_right)``
         - ``left`` -- a Map in some ``Hom(Y, Z, category_left)``

        Returns the composition of ``self`` followed by ``right`` as a
        morphism in ``Hom(X, Z, category)`` where ``category`` is the
        meet of ``category_left`` and ``category_right``.

        Caveat: see the current restrictions on :method:`Category.meet`

        EXAMPLES::

            sage: from sage.categories.morphism import SetMorphism
            sage: X.<x> = ZZ[]
            sage: Y = ZZ
            sage: Z = QQ
            sage: phi_xy = SetMorphism(Hom(X, Y, Rings()), lambda p: p[0])
            sage: phi_yz = SetMorphism(Hom(Y, Z, Monoids()), lambda y: QQ(y**2))
            sage: phi_xz = phi_xy.post_compose(phi_yz); phi_xz
            Composite map:
              From: Univariate Polynomial Ring in x over Integer Ring
              To:   Rational Field
              Defn:   Generic morphism:
                      From: Univariate Polynomial Ring in x over Integer Ring
                      To:   Integer Ring
                    then
                      Generic morphism:
                      From: Integer Ring
                      To:   Rational Field
            sage: phi_xz.category_for()
            Category of monoids
        """
        return left._composition(self)

    def extend_domain(self, new_domain):
        r"""
        INPUT:
            self          -- a member of Hom(Y, Z)
            new_codomain  -- an object X such that there is a canonical
                             coercion $\phi$ in Hom(X, Y)

        OUTPUT:
            An element of Hom(X, Z) obtained by composing self with the $\phi$.
            If no canonical $\phi$ exists, a TypeError is raised.

        EXAMPLES:
            sage: mor = CDF.coerce_map_from(RDF)
            sage: mor.extend_domain(QQ)
            Composite map:
              From: Rational Field
              To:   Complex Double Field
              Defn:   Native morphism:
                      From: Rational Field
                      To:   Real Double Field
                    then
                      Native morphism:
                      From: Real Double Field
                      To:   Complex Double Field
            sage: mor.extend_domain(ZZ['x'])
            Traceback (most recent call last):
            ...
            TypeError: No coercion from Univariate Polynomial Ring in x over Integer Ring to Real Double Field
        """
        cdef Map connecting = self.domain().coerce_map_from(new_domain)
        if connecting is None:
            raise TypeError, "No coercion from %s to %s" % (new_domain, self.domain())
        elif connecting.codomain() is not self.domain():
            raise RuntimeError, "BUG: coerce_map_from should always return a map to self (%s)" % self.domain()
        else:
            return self.pre_compose(connecting)

    def extend_codomain(self, new_codomain):
        r"""
        INPUT:
            self          -- a member of Hom(X, Y)
            new_codomain  -- an object Z such that there is a canonical
                             coercion $\phi$ in Hom(Y, Z)

        OUTPUT:
            An element of Hom(X, Z) obtained by composing self with the $\phi$.
            If no canonical $\phi$ exists, a TypeError is raised.

        EXAMPLES:
            sage: mor = QQ.coerce_map_from(ZZ)
            sage: mor.extend_codomain(RDF)
            Composite map:
              From: Integer Ring
              To:   Real Double Field
              Defn:   Natural morphism:
                      From: Integer Ring
                      To:   Rational Field
                    then
                      Native morphism:
                      From: Rational Field
                      To:   Real Double Field
            sage: mor.extend_codomain(GF(7))
            Traceback (most recent call last):
            ...
            TypeError: No coercion from Rational Field to Finite Field of size 7
        """
        cdef Map connecting = new_codomain.coerce_map_from(self.codomain())
        if connecting is None:
            raise TypeError, "No coercion from %s to %s" % (self.codomain(), new_codomain)
        elif connecting.domain() is not self.codomain():
            raise RuntimeError, "BUG: coerce_map_from should always return a map from its input (%s)" % new_codomain
        else:
            return self.post_compose(connecting)

    def is_injective(self):
        """
        Tells whether the map is injective (not implemented in the base class).

        TEST::

            sage: from sage.categories.map import Map
            sage: f = Map(Hom(QQ, ZZ, Rings()))
            sage: f.is_injective()
            Traceback (most recent call last):
            ...
            NotImplementedError: <type 'sage.categories.map.Map'>
        """
        raise NotImplementedError, type(self)

    def is_surjective(self):
        """
        Tells whether the map is surjective (not implemented in the base class).

        TEST::

            sage: from sage.categories.map import Map
            sage: f = Map(Hom(QQ, ZZ, Rings()))
            sage: f.is_surjective()
            Traceback (most recent call last):
            ...
            NotImplementedError: <type 'sage.categories.map.Map'>
        """
        raise NotImplementedError, type(self)

    def __pow__(Map self, n, dummy):
        """
        TESTS:
            sage: R.<x> = ZZ['x']
            sage: phi = R.hom([x+1]); phi
            Ring endomorphism of Univariate Polynomial Ring in x over Integer Ring
              Defn: x |--> x + 1

            sage: phi^0
            Identity endomorphism of Univariate Polynomial Ring in x over Integer Ring

            sage: phi^2 == phi*phi
            True

            sage: S.<y> = QQ[]
            sage: psi = R.hom([y^2])
            sage: psi^1
            Ring morphism:
              From: Univariate Polynomial Ring in x over Integer Ring
              To:   Univariate Polynomial Ring in y over Rational Field
              Defn: x |--> y^2
            sage: psi^2
            Traceback (most recent call last):
            ...
            TypeError: self must be an endomorphism.

            sage: K.<a> = NumberField(x^4 - 5*x + 5)
            sage: C5.<z> = CyclotomicField(5)
            sage: tau = K.hom([z - z^2]); tau
            Ring morphism:
              From: Number Field in a with defining polynomial x^4 - 5*x + 5
              To:   Cyclotomic Field of order 5 and degree 4
              Defn: a |--> -z^2 + z
            sage: tau^-1
            Ring morphism:
              From: Cyclotomic Field of order 5 and degree 4
              To:   Number Field in a with defining polynomial x^4 - 5*x + 5
              Defn: z |--> 3/11*a^3 + 4/11*a^2 + 9/11*a - 14/11

        """
        if self._domain is not self._codomain and n != 1 and n != -1:
            raise TypeError, "self must be an endomorphism."
        if n == 0:
            from sage.categories.morphism import IdentityMorphism
            return IdentityMorphism(self._parent)
        return generic_power(self, n)

    def section(self):
        """
        Return a section of self.

        NOTE:

        By default, it returns ``None``. You may override it in subclasses,
        so that it returns an instance of :class:`Section`.

        TEST::

            sage: R.<x,y> = QQ[]
            sage: f = R.hom([x+y,x-y],R)
            sage: print f.section()
            None
        """
        return None

cdef class Section(Map):
    """
    A formal section of a map.

    NOTE:

    Call methods are not implemented for the base class ``Section``.


    EXAMPLE::

        sage: from sage.categories.map import Section
        sage: R.<x,y> = ZZ[]
        sage: S.<a,b> = QQ[]
        sage: f = R.hom([a+b,a-b])
        sage: sf = Section(f); sf
        Section map:
          From: Multivariate Polynomial Ring in a, b over Rational Field
          To:   Multivariate Polynomial Ring in x, y over Integer Ring
        sage: sf(a)
        Traceback (most recent call last):
        ...
        NotImplementedError: <type 'sage.categories.map.Section'>

    """

    def __init__(self, map):
        """
        INPUT:

        A map.

        TEST::

            sage: from sage.categories.map import Section
            sage: R.<x,y> = QQ[]
            sage: f = R.hom([x+y,x-y],R)
            sage: sf = Section(f); sf
            Section map:
              From: Multivariate Polynomial Ring in x, y over Rational Field
              To:   Multivariate Polynomial Ring in x, y over Rational Field
        """
        from sage.categories.homset import Hom
        from sage.categories.all import SetsWithPartialMaps
        Map.__init__(self, Hom(map.codomain(), map.domain(), SetsWithPartialMaps()))
        self._inverse = map    # TODO: Use this attribute somewhere!

    def _repr_type(self):
        """
        Return a string describing the type of this map (which is "Section").

        TEST::

            sage: from sage.categories.map import Section
            sage: R.<x,y> = QQ[]
            sage: f = R.hom([x+y,x-y],R)
            sage: sf = Section(f)
            sage: sf         # indirect doctest
            Section map:
              From: Multivariate Polynomial Ring in x, y over Rational Field
              To:   Multivariate Polynomial Ring in x, y over Rational Field
        """
        return "Section"

cdef class FormalCompositeMap(Map):
    """
    Formal composite maps.

    A formal composite map is formed by two maps, so that the codomain of the
    first map is contained in the domain of the second map.

    NOTE:

    When calling a composite with additional arguments, these arguments are
    *only* passed to the second underlying map.

    EXAMPLE::

        sage: R.<x> = QQ[]
        sage: S.<a> = QQ[]
        sage: from sage.categories.morphism import SetMorphism
        sage: f = SetMorphism(Hom(R,S,Rings()), lambda p: p[0]*a^p.degree())
        sage: g = S.hom([2*x])
        sage: f*g
        Composite map:
          From: Univariate Polynomial Ring in a over Rational Field
          To:   Univariate Polynomial Ring in a over Rational Field
          Defn:   Ring morphism:
                  From: Univariate Polynomial Ring in a over Rational Field
                  To:   Univariate Polynomial Ring in x over Rational Field
                  Defn: a |--> 2*x
                then
                  Generic morphism:
                  From: Univariate Polynomial Ring in x over Rational Field
                  To:   Univariate Polynomial Ring in a over Rational Field
        sage: g*f
        Composite map:
          From: Univariate Polynomial Ring in x over Rational Field
          To:   Univariate Polynomial Ring in x over Rational Field
          Defn:   Generic morphism:
                  From: Univariate Polynomial Ring in x over Rational Field
                  To:   Univariate Polynomial Ring in a over Rational Field
                then
                  Ring morphism:
                  From: Univariate Polynomial Ring in a over Rational Field
                  To:   Univariate Polynomial Ring in x over Rational Field
                  Defn: a |--> 2*x
        sage: (f*g)(2*a^2+5)
        5*a^2
        sage: (g*f)(2*x^2+5)
        20*x^2
    """

    def __init__(self, parent, first, second):
        """
        INPUT:

        - ``parent``: a homset
        - ``first``, ``second``: two maps

        NOTE:

        The intended use is of course that the codomain of the
        first map is contained in the domain of the second map,
        so that the two maps can be composed, and that the
        composition belongs to ``parent``. However, none of
        these conditions is verified in the init method.

        The user is advised to compose two maps ``f`` and ``g``
        in multiplicative notation, ``g*f``, since this will in
        some cases return a more efficient map object than a
        formal composite map.

        TEST::

            sage: R.<x,y> = QQ[]
            sage: S.<a,b> = QQ[]
            sage: f = R.hom([a+b,a-b])
            sage: g = S.hom([x+y,x-y])
            sage: H = Hom(R,R,Rings())
            sage: from sage.categories.map import FormalCompositeMap
            sage: m = FormalCompositeMap(H,f,g); m
            Composite map:
              From: Multivariate Polynomial Ring in x, y over Rational Field
              To:   Multivariate Polynomial Ring in x, y over Rational Field
              Defn:   Ring morphism:
                      From: Multivariate Polynomial Ring in x, y over Rational Field
                      To:   Multivariate Polynomial Ring in a, b over Rational Field
                      Defn: x |--> a + b
                            y |--> a - b
                    then
                      Ring morphism:
                      From: Multivariate Polynomial Ring in a, b over Rational Field
                      To:   Multivariate Polynomial Ring in x, y over Rational Field
                      Defn: a |--> x + y
                            b |--> x - y
            sage: m(x), m(y)
            (2*x, 2*y)

        """
        Map.__init__(self, parent)
        self.__first = first
        self.__second = second
        self._coerce_cost = (<Map>first)._coerce_cost + (<Map>second)._coerce_cost

    cdef _update_slots(self, _slots):
        """
        Used in pickling.

        TEST::

            sage: R.<x,y> = QQ[]
            sage: S.<a,b> = QQ[]
            sage: f = R.hom([a+b,a-b])
            sage: g = S.hom([x+y,x-y])
            sage: from sage.categories.map import FormalCompositeMap
            sage: H = Hom(R,R,Rings())
            sage: m = FormalCompositeMap(H,f,g)
            sage: m == loads(dumps(m))    # indirect doctest
            True

        """
        self.__first = _slots['__first']
        self.__second = _slots['__second']
        Map._update_slots(self, _slots)

    cdef _extra_slots(self, _slots):
        """
        Used in pickling.

        TEST::

            sage: R.<x,y> = QQ[]
            sage: S.<a,b> = QQ[]
            sage: f = R.hom([a+b,a-b])
            sage: g = S.hom([x+y,x-y])
            sage: from sage.categories.map import FormalCompositeMap
            sage: H = Hom(R,R,Rings())
            sage: m = FormalCompositeMap(H,f,g)
            sage: m == loads(dumps(m))    # indirect doctest
            True

        """
        _slots['__first'] = self.__first
        _slots['__second'] = self.__second
        return Map._extra_slots(self, _slots)

    def __cmp__(self, other):
        """
        TEST::

            sage: R.<x,y> = QQ[]
            sage: S.<a,b> = QQ[]
            sage: f = R.hom([a+b,a-b])
            sage: g = S.hom([x+y,x-y])
            sage: from sage.categories.map import FormalCompositeMap
            sage: H = Hom(R,R,Rings())
            sage: m = FormalCompositeMap(H,f,g)
            sage: m == loads(dumps(m))
            True

        """
        c = cmp(type(self),type(other))
        if c == 0:
            c = cmp([self.__first,self.__second],[(<FormalCompositeMap>other).__first,(<FormalCompositeMap>other).__second])
        return c

    cpdef Element _call_(self, x):
        """
        Call with a single argument

        TEST::

            sage: R.<x> = QQ[]
            sage: S.<a> = QQ[]
            sage: from sage.categories.morphism import SetMorphism
            sage: f = SetMorphism(Hom(R,S,Rings()), lambda p: p[0]*a^p.degree())
            sage: g = S.hom([2*x])
            sage: (g*f)((x+1)^2), (f*g)((a+1)^2)     # indirect doctest
            (4*x^2, a^2)
        """
        return self.__second._call_(self.__first._call_(x))

    cpdef Element _call_with_args(self, x, args=(), kwds={}):
        """
        Additional rguments are only passed to the second map

        TEST::

            sage: from sage.categories.morphism import SetMorphism
            sage: R.<x> = QQ[]
            sage: def foo(x,*args,**kwds):
            ...    print 'foo called with',args,kwds
            ...    return x
            sage: def bar(x,*args,**kwds):
            ...    print 'bar called with',args,kwds
            ...    return x
            sage: f = SetMorphism(Hom(R,R,Rings()), foo)
            sage: b = SetMorphism(Hom(R,R,Rings()), bar)
            sage: c = b*f
            sage: c(2,'hello world',test=1)       # indirect doctest
            foo called with () {}
            bar called with ('hello world',) {'test': 1}
            2
            sage: c = f*b
            sage: c(2,'hello world',test=1)
            bar called with () {}
            foo called with ('hello world',) {'test': 1}
            2

        """
        return self.__second._call_with_args(self.__first._call_(x), args, kwds)

    def _repr_type(self):
        """
        Return a string describing the type of ``self``, namely "Composite"

        TEST::

            sage: R.<x> = QQ[]
            sage: S.<a> = QQ[]
            sage: from sage.categories.morphism import SetMorphism
            sage: f = SetMorphism(Hom(R,S,Rings()), lambda p: p[0]*a^p.degree())
            sage: g = S.hom([2*x])
            sage: f*g            # indirect doctest
            Composite map:
              From: Univariate Polynomial Ring in a over Rational Field
              To:   Univariate Polynomial Ring in a over Rational Field
              Defn:   Ring morphism:
                      From: Univariate Polynomial Ring in a over Rational Field
                      To:   Univariate Polynomial Ring in x over Rational Field
                      Defn: a |--> 2*x
                    then
                      Generic morphism:
                      From: Univariate Polynomial Ring in x over Rational Field
                      To:   Univariate Polynomial Ring in a over Rational Field
        """
        return "Composite"

    def _repr_defn(self):
        """
        Return a string describing the definition of ``self``

        The return value is obtained from the string representations
        of the two constituents.

        TEST::

            sage: R.<x> = QQ[]
            sage: S.<a> = QQ[]
            sage: from sage.categories.morphism import SetMorphism
            sage: f = SetMorphism(Hom(R,S,Rings()), lambda p: p[0]*a^p.degree())
            sage: g = S.hom([2*x])
            sage: f*g            # indirect doctest
            Composite map:
              From: Univariate Polynomial Ring in a over Rational Field
              To:   Univariate Polynomial Ring in a over Rational Field
              Defn:   Ring morphism:
                      From: Univariate Polynomial Ring in a over Rational Field
                      To:   Univariate Polynomial Ring in x over Rational Field
                      Defn: a |--> 2*x
                    then
                      Generic morphism:
                      From: Univariate Polynomial Ring in x over Rational Field
                      To:   Univariate Polynomial Ring in a over Rational Field
        """
        return "  %s\nthen\n  %s"%(self.__first, self.__second)

    def first(self):
        """
        The first map in the formal composition, where the
        composition is ``x|--> second(first(x))``.

        EXAMPLE::

            sage: R.<x> = QQ[]
            sage: S.<a> = QQ[]
            sage: from sage.categories.morphism import SetMorphism
            sage: f = SetMorphism(Hom(R,S,Rings()), lambda p: p[0]*a^p.degree())
            sage: g = S.hom([2*x])
            sage: (f*g).first() is g
            True
        """
        return self.__first

    def second(self):
        """
        The second map in the formal composition, where the
        composition is x|--> second(first(x)).

        EXAMPLE::

            sage: R.<x> = QQ[]
            sage: S.<a> = QQ[]
            sage: from sage.categories.morphism import SetMorphism
            sage: f = SetMorphism(Hom(R,S,Rings()), lambda p: p[0]*a^p.degree())
            sage: g = S.hom([2*x])
            sage: (f*g).second() is f
            True
        """
        return self.__second

    def is_injective(self):
        """
        Tell whether ``self`` is injective.

        It raises ``NotImplementedError`` if it can't be determined.

        EXAMPLE::

            sage: V1 = QQ^2
            sage: V2 = QQ^3
            sage: phi1 = (QQ^1).hom(Matrix([[1],[1]]),V1)
            sage: phi2 = V1.hom(Matrix([[1,2],[3,4],[5,6]]),V2)

        If both constituents are injective, the composition is injective::

            sage: from sage.categories.map import FormalCompositeMap
            sage: c1 = FormalCompositeMap(Hom(QQ^1,V2,phi1.category_for()),phi1,phi2)
            sage: c1.is_injective()
            True

        If it can not be determined whether the composition is injective,
        an error is raised::

            sage: psi1 = V2.hom(Matrix([[1,2,3],[4,5,6]]),V1)
            sage: c2 = FormalCompositeMap(Hom(V1,V1,phi2.category_for()),phi2,psi1)
            sage: c2.is_injective()
            Traceback (most recent call last):
            ...
            NotImplementedError: Not enough information to deduce injectivity.

        If the first map is surjective and the second map is not injective,
        then the composition is not injective::

            sage: psi2 = V1.hom([1,1],QQ^1)
            sage: c3 = FormalCompositeMap(Hom(V2,QQ^1,phi2.category_for()),psi2,psi1)
            sage: c3.is_injective()
            False

        """
        if self.__first.is_injective():
            if self.__second.is_injective():
                return True
            elif self.__first.is_surjective():
                return False
            else:
                raise NotImplementedError, "Not enough information to deduce injectivity."
        else:
            return False

    def is_surjective(self):
        """
        Tell whether ``self`` is surjective.

        It raises ``NotImplementedError`` if it can't be determined.

        EXAMPLE::

            sage: from sage.categories.map import FormalCompositeMap
            sage: V3 = QQ^3
            sage: V2 = QQ^2
            sage: V1 = QQ^1

        If both maps are surjective, the composition is surjective::

            sage: phi32 = V3.hom(Matrix([[1,2,3],[4,5,6]]),V2)
            sage: phi21 = V2.hom(Matrix([[1,1]]),V1)
            sage: c_phi = FormalCompositeMap(Hom(V3,V1,phi32.category_for()),phi32,phi21)
            sage: c_phi.is_surjective()
            True

        If the second map is not surjective, the composition is not
        surjective::

            sage: FormalCompositeMap(Hom(V3,V1,phi32.category_for()),phi32,V2.hom(Matrix([[0,0]]),V1)).is_surjective()
            False

        If the second map is an isomorphism and the first map is not
        surjective, then the composition is not surjective::

            sage: FormalCompositeMap(Hom(V2,V1,phi32.category_for()),V2.hom(Matrix([[0,0]]),V1),V1.hom(Matrix([[1]]),V1)).is_surjective()
            False

        Otherwise, surjectivity of the composition can not be determined::

            sage: FormalCompositeMap(Hom(V2,V1,phi32.category_for()),V2.hom(Matrix([[1,1],[1,1]]),V2),V2.hom(Matrix([[1,1]]),V1)).is_surjective()
            Traceback (most recent call last):
            ...
            NotImplementedError: Not enough information to deduce surjectivity.

        """
        if self.__second.is_surjective():
            if self.__first.is_surjective():
                return True
            elif self.__second.is_injective():
                return False
            else:
                raise NotImplementedError, "Not enough information to deduce surjectivity."
        else:
            return False

