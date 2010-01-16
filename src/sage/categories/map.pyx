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
    # should we use slots?
    # from element.pyx
    cdef Map mor = _class.__new__(_class)
    mor._set_parent(parent)
    mor._update_slots(_slots)
    if HAS_DICTIONARY(mor):
        mor.__dict__ = _dict
    return mor

def is_Map(x):
    return isinstance(x, Map)

cdef class Map(Element):
    def __init__(self, parent, codomain=None):
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
        A Python method to test the cdef _update_slots method

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
        A Python method to test the cdef _extra_slots method

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
        return self._extra_slots(_slots)

    def __reduce__(self):
        if HAS_DICTIONARY(self):
            _dict = self.__dict__
        else:
            _dict = {}
        return unpickle_map, (self.__class__, self._parent, _dict, self._extra_slots({}))

    def _repr_type(self):
        if self._repr_type_str is None:
            return "Generic"
        else:
            return self._repr_type_str

    def _repr_defn(self):
        return ""

    def _repr_(self):
        s = "%s map:"%self._repr_type()
        s += "\n  From: %s"%self.domain()
        s += "\n  To:   %s"%self.codomain()
        d = self._repr_defn()
        if d != '':
            s += "\n  Defn: %s"%('\n        '.join(self._repr_defn().split('\n')))
        return s

    cpdef domain(self):
        return self._domain

    cpdef codomain(self):
        return self._codomain

    def category_for(self):
        """
        Returns the category self is a morphism for

        EXAMPLES::

            sage: from sage.categories.morphism import SetMorphism
            sage: X.<x> = ZZ[]
            sage: Y = ZZ
            sage: phi = SetMorphism(Hom(X, Y, Rings()), lambda p: p[0])
            sage: phi.category_for()
            Category of rings

        FIXME: find a better name for this method
        """
        return self.parent().homset_category()

    def __call__(self, x, *args, **kwds):
        """
        Apply this map to x.

        INPUT:
            x -- an element coercible to self; also objects like
                 ideals are supported in some cases

        OUTPUT:
            an element (or ideal, etc.)

        EXAMPLES:
            sage: R.<x,y> = QQ[]; phi=R.hom([y,x])
            sage: phi(y)
            x

        We take the image of an ideal:
            sage: I = ideal(x,y); I
            Ideal (x, y) of Multivariate Polynomial Ring in x, y over Rational Field
            sage: phi(I)
            Ideal (y, x) of Multivariate Polynomial Ring in x, y over Rational Field
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
        raise NotImplementedError, type(self)

    cpdef Element _call_with_args(self, x, args=(), kwds={}):
        if len(args) == 0 and len(kwds) == 0:
            return self(x)
        else:
            raise NotImplementedError, "_call_with_args not overridden to accept arguments for %s" % type(self)

    def __mul__(self, right):
        r"""
        The multiplication * operator is operator composition

        INPUT:
            self -- Map
            right -- Map

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
        """
        if not isinstance(right, Map):
            raise TypeError, "right (=%s) must be a map to multiply it by %s"%(right, self)
        if right.codomain() != self.domain():
            raise TypeError, "self (=%s) domain must equal right (=%s) codomain"%(self, right)
        return self._composition(right)

    def _composition(self, right):
        """
        INPUT:
         - ``self``  -- a Map in some ``Hom(Y, Z, category_left)``
         - ``right`` -- a Map in some ``Hom(X, Y, category_right)``

        Returns the composition of ``self`` and ``right`` as a
        morphism in ``Hom(X, Z, category)`` where ``category`` is the
        meet of ``category_left`` and ``category_right``.

        Caveat: see the current restrictions on :method:`Category.meet`

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
        category = self.category_for().meet(right.category_for())
        H = homset.Hom(right.domain(), self.codomain(), category)
        return self._composition_(right, H)

    def _composition_(self, right, homset):
        """
        INPUT:
         - ``self``, ``right`` -- maps
         - homset -- a homset

        Returns the composition of ``self`` and ``right``, as a
        morphism in ``homset``.

        By default, this returns a formal composite map. Class
        deriving from :class:`Map` are encouraged to override this
        whenever meaningful.

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
        """
        return FormalCompositeMap(homset, right, self)

    def pre_compose(self, right):
        if self.domain() is not right.codomain():
            right = right.extend_codomain(self.domain())
        return self._composition(right)

    def post_compose(self, left):
        """
        INPUT:
         - ``self`` -- a Map in some ``Hom(X, Y, category_right)``
         - ``left`` -- a Map in some ``Hom(Y, Z, category_left)``

        Returns the composition of ``self`` and ``right`` as a
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
            sage: phi_xy.post_compose(phi_yz)
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
        raise NotImplementedError, type(self)

    def is_surjective(self):
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

            sage: phi^2
            Composite map:
              From: Univariate Polynomial Ring in x over Integer Ring
              To:   Univariate Polynomial Ring in x over Integer Ring
              Defn:   Ring endomorphism of Univariate Polynomial Ring in x over Integer Ring
                      Defn: x |--> x + 1
                    then
                      Ring endomorphism of Univariate Polynomial Ring in x over Integer Ring
                      Defn: x |--> x + 1

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
        return None

cdef class Section(Map):
    def __init__(self, map):
        from sage.categories.homset import Hom
        from sage.categories.all import SetsWithPartialMaps
        Map.__init__(self, Hom(map.codomain(), map.domain(), SetsWithPartialMaps()))
        self._inverse = map

    def _repr_type(self):
        return "Section"

cdef class FormalCompositeMap(Map):
    def __init__(self, parent, first, second):
        Map.__init__(self, parent)
        self.__first = first
        self.__second = second
        self._coerce_cost = (<Map>first)._coerce_cost + (<Map>second)._coerce_cost

    cdef _update_slots(self, _slots):
        self.__first = _slots['__first']
        self.__second = _slots['__second']
        Map._update_slots(self, _slots)

    cdef _extra_slots(self, _slots):
        _slots['__first'] = self.__first
        _slots['__second'] = self.__second
        return Map._extra_slots(self, _slots)

    cpdef Element _call_(self, x):
        return self.__second._call_(self.__first._call_(x))

    cpdef Element _call_with_args(self, x, args=(), kwds={}):
        return self.__second._call_with_args(self.__first._call_(x), args, kwds)

    def _repr_type(self):
        return "Composite"

    def _repr_defn(self):
        return "  %s\nthen\n  %s"%(self.__first, self.__second)

    def first(self):
        """
        The first map in the formal composition, where the
        composition is x|--> second(first(x)).

        """
        return self.__first

    def second(self):
        """
        The second map in the formal composition, where the
        composition is x|--> second(first(x)).
        """
        return self.__second

    def is_injective(self):
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
        if self.__second.is_surjective():
            if self.__first.is_surjective():
                return True
            elif self.__second.is_injective():
                return False
            else:
                raise NotImplementedError, "Not enough information to deduce surjectivity."
        else:
            return False

