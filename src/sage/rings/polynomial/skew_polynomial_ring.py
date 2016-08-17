r"""
Skew Univariate Polynomial Rings

This module provides the :class:`~sage.rings.polynomial.skew_polynomial_ring.SkewPolynomialRing_general`
which constructs a general dense skew univariate polynomials over commutative base rings with
automorphisms over the base rings. This is the set of formal polynomials where the coefficients
are written on the left of the variable of the skew polynomial ring. The modified multiplication
operation over elements of the base ring is extended to all elements of the skew poynomial ring
by associativity and distributivity.

This module also provides :class:`~sage.rings.polynomial.skew_polynomial_ring.SkewPolynomialRing_finite_field`
which is a specialized class for skew polynomial rings over finite fields. It inherits from
:class:`~sage.rings.polynomial.skew_polynomial_ring.SkewPolynomialRing_general` but allows for
the more efficient computations in the case of finite fields.

AUTHOR:

- Xavier Caruso (2012-06-29): initial version

- Arpit Merchant (2016-08-04): improved docstrings, fixed doctests and refactored classes and methods

- Johan Rosenkilde (2016-08-03): changes for bug fixes, docstring and doctest errors

"""

#############################################################################
#    Copyright (C) 2012 Xavier Caruso <xavier.caruso@normalesup.org>
#
#    This program is free software: you can redistribute it and/or modify
#    it under the terms of the GNU General Public License as published by
#    the Free Software Foundation, either version 2 of the License, or
#    (at your option) any later version.
#                  http://www.gnu.org/licenses/
#****************************************************************************

from sage.structure.unique_representation import UniqueRepresentation
from sage.structure.element import Element
import sage.algebras.algebra
import sage.categories.basic as categories
from sage.rings.integer import Integer
from sage.structure.category_object import normalize_names
from sage.misc.prandom import randint
from sage.categories.morphism import Morphism
from sage.categories.morphism import IdentityMorphism
import sage.misc.latex as latex
from sage.rings.polynomial.skew_polynomial_element import SkewPolynomial

#########################################################################################

class SectionSkewPolynomialCenterInjection(Section):
    def _call_ (self, x):
        order = self.inverse()._order
        section = self.inverse()._embed.section()
        lx = x.list()
        l = [ ]
        mod = 0
        for c in lx:
            if mod == 0:
                l.append(section(c))
            else:
                if not c.is_zero():
                    raise ValueError("%s is not in the center" % x)
            mod += 1
            if mod == order:
                mod = 0
        return self.codomain()(l)


class SkewPolynomialCenterInjection(RingHomomorphism):
    def __init__(self,domain,codomain,embed,order):
        RingHomomorphism.__init__(self,Hom(domain,codomain))
        self._embed = embed
        self._order = order
        self._codomain = codomain
        self._section = SectionSkewPolynomialCenterInjection(self)

    def _repr_(self):
        return "Embedding of the center of %s into this ring" % self._codomain

    def _call_(self,x):
        k = self._codomain.base_ring ()
        l = [ ]
        lz = [ k(0) ] * (self._order-1)
        for c in x.list():
            l += [ self._embed(c) ] + lz
        return self._codomain (l)

    def section(self):
        return self._section


class CenterSkewPolynomialRing(PolynomialRing_general):
    """
    A specific class for the center of a skew polynomial ring.
    """

    def __init__ (self, skew_ring, names=None, sparse=False, element_class=None):
        if not isinstance (skew_ring, SkewPolynomialRing_general):
            raise TypeError("%s is not a Skew Polynomial Ring" % skew_ring)
        self._skew_ring = skew_ring
        base_ring = skew_ring.base_ring()
        kfixed, embed = skew_ring._map.fixed_points()
        self._embed_basering = embed
        order = skew_ring._map.order()
        if order == Infinity:
            raise NotImplementedError
        self.__is_sparse = sparse
        self._PolynomialRing_general__is_sparse = sparse
        if element_class:
            self._polynomial_class = element_class
        else:
            if sparse:
                raise NotImplementedError("sparse skew polynomials are not implemented")
            else:
                self._polynomial_class = sage.rings.polynomial.skew_polynomial_element.CenterSkewPolynomial_generic_dense

        # Algebra.__init__ also calls __init_extra__ of Algebras(...).parent_class, which
        # tries to provide a conversion from the base ring, if it does not exist.
        # This is for algebras that only do the generic stuff in their initialisation.
        # But here, we want to use PolynomialBaseringInjection. Hence, we need to
        # wipe the memory and construct the conversion from scratch.
        sage.algebras.algebra.Algebra.__init__(self, kfixed, names=names, normalize=True, category=None)

        if names is None:
            if order == 1:
                self._variable_name = skew_ring.variable_name()
                self._latex_variable_name = skew_ring.latex_variable_names()[0]
                self._parenthesis = False
            else:
                self._variable_name = skew_ring.variable_name () + "^" + str(order)
                self._latex_variable_name = skew_ring.latex_variable_names()[0] + "^{" + str (order) + "}"
                self._parenthesis = True
        else:
            self._variable_name = sage.algebras.algebra.Algebra.variable_name(self)
            self._latex_variable_name = sage.algebras.algebra.Algebra.latex_variable_names(self)[0]
            self._parenthesis = False
        self._names = [ self._variable_name ]
        self.__generator = self._polynomial_class (self, [0,1], is_gen=True)
        base_inject = PolynomialBaseringInjection(kfixed,self)
        center_inject = SkewPolynomialCenterInjection (self, skew_ring, embed, order)
        self._unset_coercions_used()
        self._populate_coercion_lists_(
            coerce_list = [base_inject],
            convert_list = [list, base_inject],
            embedding = center_inject)

    def _repr_ (self):
        """
        Return a string representation of this ring.
        """
        s = "Center of %s:\n" % self._skew_ring
        s += PolynomialRing_general._repr_(self)
        return s

    def _latex_ (self):
        """
        Return a latex representation of this ring.
        """
        return "%s[%s]"%(latex.latex(self.base_ring()), self._latex_variable_name)

    def gen (self,n=0):
        """
        Return the generator of this ring.

        EXAMPLES::

            sage: k.<t> = GF(2^10)
            sage: Frob = k.frobenius_endomorphism()
            sage: S.<x> = k['x',Frob]
            sage: Z = S.center()
            sage: Z.gen()
            (x^10)
        """
        if n == 0:
            return self.__generator
        raise IndexError

    def variable_name(self, parenthesis=True):
        """
        INPUT:

        -  ``parenthesis`` -- a boolean (default: True)

        OUTPUT:

        A string representation of the variable name of this ring.
        If ``parenthesis`` is true and the variable is not atomic,
        parenthesis are added around the variable name.

        EXAMPLES::

            sage: k.<t> = GF(3^5)
            sage: Frob = k.frobenius_endomorphism()
            sage: S.<x> = k['x',Frob]

            sage: Z = S.center(); Z
            Center of Skew Polynomial Ring in x over Finite Field in t of size 3^5 twisted by t |--> t^3:
            Univariate Polynomial Ring in (x^5) over Finite Field of size 3
            sage: Z.variable_name()
            '(x^5)'
            sage: Z.variable_name(parenthesis=False)
            'x^5'

            sage: Z = S.center(name='y'); Z
            Center of Skew Polynomial Ring in x over Finite Field in t of size 3^5 twisted by t |--> t^3:
            Univariate Polynomial Ring in y over Finite Field of size 3
            sage: Z.variable_name()
            'y'
        """
        if parenthesis and self._parenthesis:
            return "(" + self._variable_name + ")"
        else:
            return self._variable_name

    def latex_variable_names(self, parenthesis=True):
        """
        INPUT:

        -  ``parenthesis`` -- a boolean (default: True)

        OUTPUT:

        A list composed with just one element which is a latex
        representation of the variable name of this ring.
        If ``parenthesis`` is true and the variable is not atomic,
        parenthesis are added around the variable name.

        EXAMPLES::

            sage: k.<t> = GF(3^4)
            sage: Frob = k.frobenius_endomorphism()
            sage: S.<x> = k['x',Frob]

            sage: Z = S.center(); Z
            Center of Skew Polynomial Ring in x over Finite Field in t of size 3^4 twisted by t |--> t^3:
            Univariate Polynomial Ring in (x^4) over Finite Field of size 3
            sage: Z.latex_variable_names()
            ['(x^{4})']
            sage: Z.latex_variable_names(parenthesis=False)
            ['x^{4}']

            sage: Z = S.center(name='y'); Z
            Center of Skew Polynomial Ring in x over Finite Field in t of size 3^4 twisted by t |--> t^3:
            Univariate Polynomial Ring in y over Finite Field of size 3
            sage: Z.latex_variable_names()   # atomic variable
            ['y']
        """
        if parenthesis and self._parenthesis:
            return [ "(" + self._latex_variable_name + ")" ]
        else:
            return [ self._latex_variable_name ]

class SkewPolynomialRing_general(sage.algebras.algebra.Algebra,UniqueRepresentation):
    """
    General Skew Univariate polynomial ring over a ring.

    DEFINITION:

    Let `R` be a commutative ring and let `\sigma` be an automorphism over `R`. An
    automorphism (also called twist map) is a structure preserving map from a
    mathematical object (in this case, R) onto itself that also admits an inverse.
    The ring of skew polynomials over an indeterminate variable `X` is defined then,
    as the ring structure on the set `R` as:
    `R[X, \sigma] = { a_{n-1}X^{n-1} + ... + a_{1}X + a_{0} | a_{i} \in R and n \in N }`
    where the addition operation on `R[X, \sigma]` is given by the usual abelian
    group polynomial addition rule and the multiplication operation is defined
    by the modified rule `X*a = \sigma(a)X`.

    This ring is non-commutative and its elements are all such skew polynomials whose
    coefficients come from `R`.

    Reference: "Theory of Non-Commutative Polynomials" - Oystein Ore

    .. TODO::

        Add derivations.

    EXAMPLES::

        sage: R.<t> = ZZ[]
        sage: sigma = R.hom([t+1])
        sage: S.<x> = SkewPolynomialRing(R,sigma); S
        Skew Polynomial Ring in x over Univariate Polynomial Ring in t over Integer Ring twisted by t |--> t + 1

    One can also use a shorter syntax::

        sage: S.<x> = R['x',sigma]; S
        Skew Polynomial Ring in x over Univariate Polynomial Ring in t over Integer Ring twisted by t |--> t + 1

    Be careful, with the latter syntax, one cannot omit the name of the
    variable neither in LHS nor in RHS. If we omit it in LHS, the variable
    is not created::

        sage: Sy = R['y',sigma]; Sy
        Skew Polynomial Ring in y over Univariate Polynomial Ring in t over Integer Ring twisted by t |--> t + 1
        sage: y.parent()
        Traceback (most recent call last):
        ...
        NameError: name 'y' is not defined

    If we omit it in RHS, sage tries to create a polynomial ring and fails::

        sage: Sz.<z> = R[sigma]
        Traceback (most recent call last):
        ...
        ValueError: variable name 'Ring endomorphism of Univariate Polynomial Ring in t over Integer Ring\n
        Defn: t |--> t + 1' is not alphanumeric

    As for polynomials, skew polynomial rings with different variable names
    are not equal::

        sage: R['x',sigma] == R['y',sigma]
        False

    Of course, skew polynomial rings with different twist maps are not
    equal as well::

        sage: R['x',sigma] == R['x',sigma^2]
        False

    Saving and loading of polynomial rings works::

        sage: loads(dumps(R['x',sigma])) == R['x',sigma]
        True

    There is a coercion map from the base ring of the skew polynomial rings::

        sage: S.has_coerce_map_from(R)
        True
        sage: x.parent()
        Skew Polynomial Ring in x over Univariate Polynomial Ring in t over Integer Ring twisted by t |--> t + 1
        sage: t.parent()
        Univariate Polynomial Ring in t over Integer Ring
        sage: y = x+t; y
        x + t
        sage: y.parent() is S
        True

    .. SEE ALSO::

        :meth:`sage.rings.polynomial.skew_polynomial_ring_constructor.SkewPolynomialRing`
        :mod:`sage.rings.polynomial.skew_polynomial_element`
    """
    @staticmethod
    def __classcall__(cls, base_ring, map, name=None, sparse=False, element_class=None):
        """
        Input name mangling for `SkewPolynomialRing_general` class into
        the `SkewPolynomialRing_general_with_category` class so that it
        inherit all the methods from the super class. Sets the default
        values for `name`, `sparse` and `element_class`.

        EXAMPLES:

            sage: R.<t> = ZZ[]
            sage: sigma = R.hom([t+1])
            sage: S.<x> = SkewPolynomialRing(R,sigma)
            sage: S.__class__(R, sigma, x)
            Skew Polynomial Ring in x over Univariate Polynomial Ring in t over Integer Ring twisted by t |--> t + 1
            sage: type(S)
            <class 'sage.rings.polynomial.skew_polynomial_ring.SkewPolynomialRing_general_with_category'>
        """
        if not element_class:
            if sparse:
                raise NotImplementedError("sparse skew polynomials are not implemented")
            else:
                from sage.rings.polynomial import skew_polynomial_element
                element_class = skew_polynomial_element.SkewPolynomial_generic_dense
        return super(SkewPolynomialRing_general,cls).__classcall__(cls,base_ring,map,name,sparse,element_class)

    def __init__(self, base_ring, map, name, sparse, element_class):
        """
        This method is a constructor for a general, dense univariate skew polynomial ring.

        INPUT::

        - ``base_ring`` -- a commutative ring

        - ``map`` -- an automorphism of the base ring

        - ``name`` -- string or list of strings representing the name of the variables of ring

        - ``sparse`` -- boolean (default: ``False``)

        - ``element_class`` -- class representing the type of element to be used in ring

        EXAMPLES::

            sage: R.<t> = ZZ[]
            sage: sigma = R.hom([t+1])
            sage: S.<x> = SkewPolynomialRing(R,sigma); S
            Skew Polynomial Ring in x over Univariate Polynomial Ring in t over Integer Ring twisted by t |--> t + 1
            sage: S([1]) + S([-1])
            0
            sage: k.<t> = GF(5^3)
            sage: Frob = k.frobenius_endomorphism()
            sage: T.<x> = k['x', Frob]; T
            Skew Polynomial Ring in x over Finite Field in t of size 5^3 twisted by t |--> t^5
        """
        category = categories.Rings()
        self.__is_sparse = sparse
        self._polynomial_class = element_class
        if map is None:
            self._map = IdentityMorphism(base_ring)
        else:
            if isinstance (map, Morphism):
                if map.domain () == base_ring and map.codomain () == base_ring:
                    self._map = map
                else:
                    raise TypeError("given map is not an automorphism of %s" % base_ring)
            else:
                raise TypeError("given map is not a ring homomorphism")
        self._maps = { 0:IdentityMorphism(base_ring), 1:self._map }
        self._center = { }
        self._center_variable = None
        self._no_generic_basering_coercion = True
        sage.algebras.algebra.Algebra.__init__(self, base_ring, names=name, normalize=True, category=category)
        self.__generator = self._polynomial_class(self, [0,1], is_gen=True)
        base_inject = sage.rings.polynomial.skew_polynomial_element.SkewPolynomialBaseringInjection(base_ring,self)
        self._populate_coercion_lists_(
                coerce_list = [base_inject],
                convert_list = [list, base_inject])

    def __reduce__(self):
        """
        Return the globally unique skew polynomial ring based on
        given arguments.

        TESTS:

            sage: R.<t> = ZZ[]
            sage: sigma = R.hom([t+1])
            sage: S.<x> = SkewPolynomialRing(R,sigma); S
            Skew Polynomial Ring in x over Univariate Polynomial Ring in t over Integer Ring twisted by t |--> t + 1
            sage: T.<x> = SkewPolynomialRing(R,sigma); T
            Skew Polynomial Ring in x over Univariate Polynomial Ring in t over Integer Ring twisted by t |--> t + 1
            sage: S is T
            True
        """
        import sage.rings.polynomial.skew_polynomial_ring_constructor
        return (sage.rings.polynomial.skew_polynomial_ring_constructor.SkewPolynomialRing,
                (self.base_ring(), self.twist_map(), self.variable_name(), self.is_sparse()))

    def _element_constructor_(self, x=None, check=True, is_gen = False, construct=False, **kwds):
        """
        Convert ``x`` into an element of this univariate skew polynomial ring,
        possibly non-canonically.

        INPUT:

        - ``x`` -- an element of the base ring of ``self`` or a ring that
          has a coerce map from ``self`` (default: ``None``).

        - ``check`` -- boolean (default: ``True``)

        - ``is_gen`` -- boolean (default: ``False``)

        - ``construct`` -- boolean (default: ``False``)

        OUTPUT:

        An element of ``self``.

        EXAMPLES::

            sage: R.<t> = ZZ[]
            sage: sigma = R.hom([t+1])
            sage: S.<x> = SkewPolynomialRing(R,sigma)
            sage: S(1 + x + x^2 + x^3)
            x^3 + x^2 + x + 1
            sage: S(1 + t)
            t + 1
            sage: S(1 + t).degree()
            0
            sage: S(0).list()
            []
        """
        C = self._polynomial_class
        if isinstance(x, list):
            return C(self, x, check=check, is_gen=False,construct=construct)
        if isinstance(x, Element):
            P = x.parent()
            def build(check):
                if x.is_zero():
                    return P.zero()
                else:
                    return C(self, [x], check=check, is_gen=False, construct=construct)
            if P is self:
                return x
            elif P is self.base_ring():
                build(False)
            elif P == self.base_ring() or self.base_ring().has_coerce_map_from(P):
                build(True)
        try:
            return x._polynomial_(self)
        except AttributeError:
            pass
        if isinstance(x,str):
            try:
                from sage.misc.parser import Parser, LookupNameMaker
                R = self.base_ring()
                p = Parser(Integer, R, LookupNameMaker({self.variable_name(): self.gen()}, R))
                return self(p.parse(x))
            except NameError:
                raise TypeError("unable to coerce string")
        return C(self, x, check, is_gen, construct=construct, **kwds)

    def _coerce_map_from_(self, P):
        """
        Check whether ``self`` has a coerce map from ``P``.

        The rings that canonically coerce into this ring are:

        - this ring itself

        - any ring that canonically coerces to the base ring of this ring

        - skew polynomial rings in the same variable and automorphism over
          any base ring that canonically coerces to the base ring of this ring

        INPUT:

        - ``P`` -- a ring.

        OUTPUT:

        Return ``True`` or ``False``.

        .. NOTE::

            Sparse skew polynomials are not implemented.

        EXAMPLES::

            sage: R.<t> = ZZ[]
            sage: sigma = R.hom([t+1])
            sage: S.<x> = SkewPolynomialRing(R,sigma)
            sage: S.has_coerce_map_from(S)
            True
            sage: S.has_coerce_map_from(R)
            True
            sage: S.has_coerce_map_from(ZZ)
            True
            sage: S.has_coerce_map_from(GF(5^3))
            False

            sage: S.coerce_map_from(ZZ)
            Composite map:
              From: Integer Ring
              To:   Skew Polynomial Ring in x over Univariate Polynomial Ring in t over Integer Ring twisted by t |--> t + 1
              Defn:   Polynomial base injection morphism:
                      From: Integer Ring
                      To:   Univariate Polynomial Ring in t over Integer Ring
                    then
                      Skew Polynomial base injection morphism:
                      From: Univariate Polynomial Ring in t over Integer Ring
                      To:   Skew Polynomial Ring in x over Univariate Polynomial Ring in t over Integer Ring twisted by t |--> t + 1
            sage: S.coerce_map_from(S)
            Identity endomorphism of Skew Polynomial Ring in x over Univariate Polynomial Ring in t over Integer Ring twisted by t |--> t + 1
        """
        try:
            connecting = self.base_ring().coerce_map_from(P)
            if connecting is not None:
                return self.coerce_map_from(self.base_ring()) * connecting
        except TypeError:
            pass
        try:
            if isinstance(P, SkewPolynomialRing_general):
                if self.__is_sparse and not P.is_sparse():
                    return False
                if P.variable_name() == self.variable_name():
                    if P.base_ring() is self.base_ring() and \
                            self.base_ring() is ZZ_sage:
                       if self._implementation_names == ('NTL',):
                            return False
                    return self.base_ring().has_coerce_map_from(P.base_ring())
        except AttributeError:
            pass

    def _repr_(self):
        """
        Return string representation of ``self``.

        EXAMPLES::

            sage: R.<t> = ZZ[]
            sage: sigma = R.hom([t+1])
            sage: S.<x> = SkewPolynomialRing(R,sigma)
            sage: S
            Skew Polynomial Ring in x over Univariate Polynomial Ring in t over Integer Ring twisted by t |--> t + 1
        """
        s = "Skew Polynomial Ring in %s over %s twisted by %s"%(self.variable_name(), self.base_ring(), self._map._repr_short())
        if self.is_sparse():
            s = "Sparse " + s
        return s

    def _latex_(self):
        """
        Return latex representation of ``self``.

        EXAMPLES::

            sage: R.<t> = ZZ[]
            sage: sigma = R.hom([t+1])
            sage: S.<x> = SkewPolynomialRing(R,sigma)
            sage: latex(S)
            \Bold{Z}[t][x,\begin{array}{l}
            \text{\texttt{Ring{ }endomorphism{ }of{ }Univariate{ }Polynomial{ }Ring{ }in{ }t{ }over{ }Integer{ }Ring}}\\
            \text{\texttt{{ }{ }Defn:{ }t{ }|{-}{-}>{ }t{ }+{ }1}}
            \end{array}]
        """
        return "%s[%s,%s]"%(latex.latex(self.base_ring()), self.latex_variable_names()[0], latex.latex(self._map))

    def change_var(self, var):
        r"""
        Return the skew polynomial ring in variable ``var`` over the same base
        ring.

        INPUT:

        - ``var`` -- a string representing the name of the new variable of ``self``

        OUTPUT:

        ``self`` with variable name name changed to ``var``.

        EXAMPLES::

            sage: k.<t> = GF(5^3)
            sage: Frob = k.frobenius_endomorphism()
            sage: R.<x> = SkewPolynomialRing(k,Frob); R
            Skew Polynomial Ring in x over Finite Field in t of size 5^3 twisted by t |--> t^5
            sage: R.change_var('y')
            Skew Polynomial Ring in y over Finite Field in t of size 5^3 twisted by t |--> t^5
        """
        from sage.rings.polynomial.skew_polynomial_ring_constructor import SkewPolynomialRing
        return SkewPolynomialRing(self.base_ring(), self.twist_map(), names = var, sparse=self.is_sparse())

    def characteristic(self):
        """
        Return the characteristic of the base ring of this skew polynomial ring.

        EXAMPLES::

            sage: R.<t> = QQ[]
            sage: sigma = R.hom([t+1])
            sage: R['x',sigma].characteristic()
            0

            sage: k.<u> = GF(5^3)
            sage: Frob = k.frobenius_endomorphism()
            sage: k['y',Frob].characteristic()
            5
        """
        return self.base_ring().characteristic()

    def twist_map(self, n=1):
        """
        Return the twist map, otherwise known as the automorphism over the base ring of
        ``self``, iterated `n` times.

        INPUT:

        -  ``n`` - a relative integer (default: 1)

        OUTPUT:

        -  The `n`-th iterative of the twist map of this skew polynomial ring.

        EXAMPLES::

            sage: R.<t> = QQ[]
            sage: sigma = R.hom([t+1])
            sage: S.<x> = R['x',sigma]
            sage: S.twist_map()
            Ring endomorphism of Univariate Polynomial Ring in t over Rational Field
              Defn: t |--> t + 1
            sage: S.twist_map() == sigma
            True
            sage: S.twist_map(10)
            Ring endomorphism of Univariate Polynomial Ring in t over Rational Field
              Defn: t |--> t + 10

        If `n` in negative, Sage tries to compute the inverse of the twist map.
        Sometimes it fails (even if the twist map is actually invertible)::

            sage: S.twist_map(-1)
            Traceback (most recent call last):
            ...
            NotImplementedError: inversion of the twist map Ring endomorphism of Univariate Polynomial Ring in t over Rational Field
                  Defn: t |--> t + 1

        Sometimes it succeeds::

            sage: k.<t> = GF(5^3)
            sage: Frob = k.frobenius_endomorphism()
            sage: T.<y> = k['y',Frob]
            sage: T.twist_map(-1)
            Frobenius endomorphism t |--> t^(5^2) on Finite Field in t of size 5^3
        """
        try:
            return self._maps[n]
        except KeyError:
            if n >= 0:
                map = self._map**n
                self._maps[n] = map
                return map
            else:
                try:
                    map = self._map**n
                except TypeError:
                    raise NotImplementedError("inversion of the twist map %s" % self._map)
                self._maps[n] = map
                return map

    def gen(self, n=0):
        """
        Return the indeterminate generator of this skew polynomial ring.

        INPUT:

        - ``n`` -- an integer (default: 0)

        EXAMPLES::

            sage: R.<t> = QQ[]
            sage: sigma = R.hom([t+1])
            sage: S.<x> = R['x',sigma]; S
            Skew Polynomial Ring in x over Univariate Polynomial Ring in t over Rational Field twisted by t |--> t + 1
            sage: y = S.gen(); y
            x
            sage: y == x
            True
            sage: y is x
            True
            sage: y.is_gen()
            True
        """
        if n != 0:
            raise IndexError("generator n not defined")
        return self.__generator

    def gens_dict(self):
        """
        Return a dictionary whose entries are ``{name:variable,...}``,
        where ``name`` stands for the variable names of this
        object (as strings) and ``variable`` stands for the corresponding
        generators (as elements of this object).

        EXAMPLES::

            sage: R.<t> = ZZ[]
            sage: sigma = R.hom([t+1])
            sage: S.<x> = SkewPolynomialRing(R,sigma)
            sage: S.gens_dict()
            {'x': x}
        """
        return dict(zip(self.variable_names(), self.gens()))

    def parameter(self):
        """
        Return the generator of this skew polynomial ring.

        This is the same as ``self.gen()``.

        EXAMPLES::

            sage: R.<t> = QQ[]
            sage: sigma = R.hom([t+1])
            sage: S.<x> = R['x',sigma]; S
            Skew Polynomial Ring in x over Univariate Polynomial Ring in t over Rational Field twisted by t |--> t + 1
            sage: y = S.parameter(); y
            x
        """
        return self.gen()

    def is_finite(self):
        """
        Return ``False`` since skew polynomial rings are not finite (unless the
        base ring is 0.)

        EXAMPLES::

            sage: k.<t> = GF(5^3)
            sage: k.is_finite()
            True
            sage: Frob = k.frobenius_endomorphism()
            sage: S.<x> = k['x',Frob]
            sage: S.is_finite()
            False
        """
        R = self.base_ring()
        if R.is_finite() and R.order() == 1:
            return True
        return False

    def is_exact(self):
        """
        Return ``True`` if elements of this skew polynomial ring are exact.
        It happens if and only if elements of the base ring are exact.

        EXAMPLES::

            sage: k.<t> = GF(5^3)
            sage: Frob = k.frobenius_endomorphism()
            sage: S.<x> = k['x',Frob]
            sage: S.is_exact()
            True
            sage: S.base_ring().is_exact()
            True

            sage: R.<u> = k[[]]
            sage: sigma = R.hom([u+u^2])
            sage: T.<y> = R['y',sigma]
            sage: T.is_exact()
            False
            sage: T.base_ring().is_exact()
            False
        """
        return self.base_ring().is_exact()

    def is_sparse(self):
        """
        Return ``True`` if elements of this polynomial ring have a sparse
        representation.

        Since sparse skew polynomials are not yet implemented, this
        function always returns ``False``.

        EXAMPLES:

            sage: R.<t> = RR[]
            sage: sigma = R.hom([t+1])
            sage: S.<x> = R['x',sigma]
            sage: S.is_sparse()
            False
        """
        return self.__is_sparse

    def ngens(self):
        """
        Return the number of generators of this skew polynomial ring, which is 1.

        EXAMPLES::

            sage: R.<t> = RR[]
            sage: sigma = R.hom([t+1])
            sage: S.<x> = R['x',sigma]
            sage: S.ngens()
            1
        """
        return 1

    def random_element(self, degree=2, monic=False, *args, **kwds):
        r"""
        Return a random skew polynomial.

        INPUT:

        -  ``degree`` - Integer with degree (default: 2)
           or a tuple of integers with minimum and maximum degrees

        -  ``monic`` - if True, returns a monic skew polynomial
           (default: ``False``)

        -  ``*args, **kwds`` - Passed on to the ``random_element`` method for
           the base ring

        OUTPUT:

        -  Skew polynomial such that the coefficients of `x^i`, for `i` up
           to ``degree``, are random elements from the base ring, randomized
           subject to the arguments ``*args`` and ``**kwds``

        EXAMPLES::

            sage: k.<t> = GF(5^3)
            sage: Frob = k.frobenius_endomorphism()
            sage: S.<x> = k['x',Frob]
            sage: S.random_element()  # random
            (2*t^2 + 3)*x^2 + (4*t^2 + t + 4)*x + 2*t^2 + 2
            sage: S.random_element(monic=True)  # random
            x^2 + (2*t^2 + t + 1)*x + 3*t^2 + 3*t + 2

        If a tuple of two integers is given for the degree argument, a random
        integer will be chosen between the first and second element of the
        tuple as the degree::

            sage: S.random_element(degree=(2,7))  # random
            (3*t^2 + 1)*x^4 + (4*t + 2)*x^3 + (4*t + 1)*x^2 + (t^2 + 3*t + 3)*x + 3*t^2 + 2*t + 2

        If the minimal degree is greater than the maximal degree, sage raises
        a ValueError::

            sage: S.random_element(degree=(5,4))
            Traceback (most recent call last):
            ...
            ValueError: minimum degree must be less or equal than maximum degree

        When ``monic`` is false, the returned skew polynomial may have a degree
        less than ``degree`` (it happens when the random ``leading coefficient``
        is zero)::

            sage: S.random_element(degree=4) #random
            (3*t^2 + t)*x^3 + (2*t + 2)*x^2 + (3*t^2 + 2*t + 2)*x + t

        However, if ``monic`` is true, this can't happen::

            sage: S.random_element(degree=4,monic=True)  # random
            x^4 + (t^2 + 3*t + 3)*x^3 + (t^2 + 4*t + 3)*x^2 + (2*t^2 + 4*t + 4)*x + 2*t^2 + 2*t
        """
        R = self.base_ring()
        if isinstance(degree, (list, tuple)):
            if len(degree) != 2:
                raise ValueError("degree argument must be an integer or a tuple of 2 integers (min_degree, max_degree)")
            if degree[0] > degree[1]:
                raise ValueError("minimum degree must be less or equal than maximum degree")
            degree = randint(*degree)
        if monic:
            return self ([ R.random_element (*args, **kwds) for _ in range (degree) ] + [ R.one() ])
        else:
            return self ([ R.random_element (*args, **kwds) for _ in range (degree+1) ])

    def is_commutative(self):
        """
        Return ``True`` if this skew polynomial ring is commutative
        (i.e. if the twist map is the identity).

        EXAMPLES::

            sage: k.<t> = GF(5^3)
            sage: Frob = k.frobenius_endomorphism()
            sage: S.<x> = k['x',Frob]
            sage: S.is_commutative()
            False

            sage: T.<y> = k['y',Frob^3]
            sage: T.is_commutative()
            True
        """
        return self.twist_map().is_identity()

    def center(self,names=None,name=None):
        r"""
        Return the center of this skew polynomial ring.

        .. NOTE::

            If F denotes the subring of R fixed by `\sigma`, the center of
            `R[X,\sigma]` is `F` if `\sigma` has infinite order and `F[X^r]`
            if `\sigma` has finite order `r`.

        .. WARNING::

            This function assumes that `\sigma` has a method order() (which
            returns its order) and a method fixed_points() (which returns
            the subring `F` together with the embedding of `F` into `R`).
            The case where `\sigma` has infinite order is not implemented
            yet.

        INPUT:

        - ``name`` -- a string (or None)

        OUTPUT:

        The center of this skew polynomial ring.

        If ``name`` is given, the name of the variable of the center
        (which is a polynomial ring) is assigned to �t.
        Otherwise, the notation `(x^r)` (where `x` is the name of the
        variable of this skew polynomial ring and `r` is the order of
        `\sigma`) is used.

        EXAMPLES::

            sage: k.<t> = GF(5^3)
            sage: Frob = k.frobenius_endomorphism()
            sage: S.<x> = k['x',Frob]; S
            Skew Polynomial Ring in x over Finite Field in t of size 5^3 twisted by t |--> t^5
            sage: Z = S.center(); Z
            Center of Skew Polynomial Ring in x over Finite Field in t of size 5^3 twisted by t |--> t^5:
            Univariate Polynomial Ring in (x^3) over Finite Field of size 5
            sage: Z.gen()
            (x^3)

        We can also specify another variable name::

            sage: Zy.<y> = S.center(); Zy
            Center of Skew Polynomial Ring in x over Finite Field in t of size 5^3 twisted by t |--> t^5:
            Univariate Polynomial Ring in y over Finite Field of size 5
            sage: y.parent() == Zy
            True

        Coercion from the center into the skew polynomial ring works::

            sage: a = S.random_element(); a
            (2*t^2 + 3)*x^2 + (4*t^2 + t + 4)*x + 2*t^2 + 2
            sage: b = Z.random_element(); b
            3*(x^3) + 2
            sage: c = a + b; c
            3*x^3 + (2*t^2 + 3)*x^2 + (4*t^2 + t + 4)*x + 2*t^2 + 4
            sage: c.parent()
            Skew Polynomial Ring in x over Finite Field in t of size 5^3 twisted by t |--> t^5

        We also have a section map in the other direction::

            sage: z = x^6 + 2*x^3
            sage: z.parent()
            Skew Polynomial Ring in x over Finite Field in t of size 5^3 twisted by t |--> t^5
            sage: zz = Z(z); zz
            (x^3)^2 + 2*(x^3)
            sage: zz.parent()
            Center of Skew Polynomial Ring in x over Finite Field in t of size 5^3 twisted by t |--> t^5:
            Univariate Polynomial Ring in (x^3) over Finite Field of size 5

            sage: v = x^4
            sage: Z(v)
            Traceback (most recent call last):
            ...
            ValueError: x^4 is not in the center
        """
        if name is None:
            name = names
        if name is None:
            name = self._center_variable
        else:
            name = normalize_names(1, name)[0]
            self._center_variable = name
        try:
            return self._center[name]
        except KeyError:
            self._center[name] = CenterSkewPolynomialRing(self, sparse=self.__is_sparse, names=name)
            return self._center[name]

    def centre(self,names=None,name=None):
        r"""
        Return the centre of this skew polynomial ring.

        .. NOTE::

            If F denotes the subring of R fixed by `\sigma`, the centre of
            `R[X,\sigma]` is `F` if `\sigma` has infinite order and `F[X^r]`
            if `\sigma` has finite order `r`.

        .. WARNING::

            This function assumes that `\sigma` has a method order() (which
            returns its order) and a method fixed_points() (which returns
            the subring `F` together with the embedding of `F` into `R`).
            The case where `\sigma` has infinite order is not implemented
            yet.

        OUTPUT:

        The centre of this skew polynomial ring.

        If ``name`` is given, the name of the variable of the centre
        (which is a polynomial ring) is assigned to �t.
        Otherwise, the notation `(x^r)` (where `x` is the name of the
        variable of this skew polynomial ring and `r` is the order of
        `\sigma`) is used.

        EXAMPLES::

            sage: k.<t> = GF(7^3)
            sage: Frob = k.frobenius_endomorphism()
            sage: S.<x> = k['x',Frob]; S
            Skew Polynomial Ring in x over Finite Field in t of size 7^3 twisted by t |--> t^7
            sage: Z = S.centre(); Z
            Center of Skew Polynomial Ring in x over Finite Field in t of size 7^3 twisted by t |--> t^7:
            Univariate Polynomial Ring in (x^3) over Finite Field of size 7
            sage: Z.gen()
            (x^3)

        We can also specify another variable name::

            sage: Zy.<y> = S.centre(); Zy
            Center of Skew Polynomial Ring in x over Finite Field in t of size 7^3 twisted by t |--> t^7:
            Univariate Polynomial Ring in y over Finite Field of size 7
            sage: y.parent() == Zy
            True

        Coercion from the centre into the skew polynomial ring works::

            sage: a = (3*t^2 + 4*t + 6)*x^2 + (2*t + 1)*x + 3*t^2 + 6*t + 4
            sage: b = 4 * Z.gen() + 3; b
            4*(x^3) + 3
            sage: c = a + b; c
            4*x^3 + (3*t^2 + 4*t + 6)*x^2 + (2*t + 1)*x + 3*t^2 + 6*t
            sage: c.parent()
            Skew Polynomial Ring in x over Finite Field in t of size 7^3 twisted by t |--> t^7

        We also have a section map in the other direction::

            sage: z = x^6 + 2*x^3
            sage: z.parent()
            Skew Polynomial Ring in x over Finite Field in t of size 7^3 twisted by t |--> t^7
            sage: zz = Z(z); zz
            (x^3)^2 + 2*(x^3)
            sage: zz.parent()
            Center of Skew Polynomial Ring in x over Finite Field in t of size 7^3 twisted by t |--> t^7:
            Univariate Polynomial Ring in (x^3) over Finite Field of size 7

            sage: v = x^4
            sage: Z(v)
            Traceback (most recent call last):
            ...
            ValueError: x^4 is not in the center
        """
        return self.center(name=name,names=names)    

class SkewPolynomialRing_finite_field(SkewPolynomialRing_general):
    """
    A specialized class for skew polynomial rings over finite fields.

    .. SEEALSO::

        :meth:`sage.rings.polynomial.skew_polynomial_ring_constructor.SkewPolynomialRing`
        :class:`sage.rings.polynomial.skew_polynomial_ring.SkewPolynomialRing_general`
        :mod:`sage.rings.polynomial.skew_polynomial_finite_field`

    .. TODO::

        Add methods related to center of skew polynomial ring, irreducibility, karatsuba
        multiplication and factorization.
    """
    @staticmethod
    def __classcall__(cls, base_ring, map, name=None, sparse=False, element_class=None):
        if not element_class:
            if sparse:
                raise NotImplementedError("sparse skew polynomials are not implemented")
            else:
                from sage.rings.polynomial import skew_polynomial_finite_field
                element_class = skew_polynomial_finite_field.SkewPolynomial_finite_field_dense
                return super(SkewPolynomialRing_general,cls).__classcall__(cls,base_ring,map,name,sparse,element_class)

    def __init__(self, base_ring, map, name, sparse, element_class):
        """
        This method is a constructor for a general, dense univariate skew polynomial ring
        over a finite field.

        INPUT::

        - ``base_ring`` -- a commutative ring

        - ``map`` -- an automorphism of the base ring

        - ``name`` -- string or list of strings representing the name of the variables of ring

        - ``sparse`` -- boolean (default: ``False``)

        - ``element_class`` -- class representing the type of element to be used in ring

        ..NOTE::

            Multivariate and Sparse rings are not implemented.

        EXAMPLES::

            sage: k.<t> = GF(5^3)
            sage: Frob = k.frobenius_endomorphism()
            sage: T.<x> = k['x', Frob]; T
            Skew Polynomial Ring in x over Finite Field in t of size 5^3 twisted by t |--> t^5
        """
        self._order = -1
        try:
            self._order = map.order()
        except (AttributeError,NotImplementedError):
            pass
        if self._order < 0:
            try:
                if map.is_identity():
                    self._order = 1
            except (AttributeError,NotImplementedError):
                pass
        if self._order < 0:
            raise NotImplementedError("unable to determine the order of %s" % map)
        SkewPolynomialRing_general.__init__ (self, base_ring, map, name, sparse, element_class)
        self._maps = [ map**i for i in range(self._order) ]

    def twist_map(self, n=1):
        """
        Return the twist map, otherwise known as the automorphism over the base ring of
        ``self``, iterated `n` times.

        INPUT:

        -  ``n`` - a relative integer (default: 1)

        OUTPUT:

        -  The `n`-th iterative of the twist map of this skew polynomial ring.

        EXAMPLES::

            sage: k.<t> = GF(5^3)
            sage: Frob = k.frobenius_endomorphism()
            sage: S.<x> = k['x',Frob]
            sage: S.twist_map()
            Frobenius endomorphism t |--> t^5 on Finite Field in t of size 5^3
            sage: S.twist_map(11)
            Frobenius endomorphism t |--> t^(5^2) on Finite Field in t of size 5^3
            sage: S.twist_map(3)
            Identity endomorphism of Finite Field in t of size 5^3

        It also works if `n` is negative::

            sage: S.twist_map(-1)
            Frobenius endomorphism t |--> t^(5^2) on Finite Field in t of size 5^3
        """
        return self._maps[n%self._order]
