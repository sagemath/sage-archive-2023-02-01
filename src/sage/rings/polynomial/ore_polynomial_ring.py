r"""
Univariate Ore Polynomial Rings

AUTHOR:

- Xavier Caruso (2020-04)

"""
# ***************************************************************************
#    Copyright (C) 2020 Xavier Caruso <xavier.caruso@normalesup.org>
#
#    This program is free software: you can redistribute it and/or modify
#    it under the terms of the GNU General Public License as published by
#    the Free Software Foundation, either version 2 of the License, or
#    (at your option) any later version.
#                  https://www.gnu.org/licenses/
# ***************************************************************************

import sage

from sage.structure.richcmp import op_EQ
from sage.misc.prandom import randint
from sage.misc.cachefunc import cached_method
from sage.rings.infinity import Infinity
from sage.structure.category_object import normalize_names

from sage.structure.unique_representation import UniqueRepresentation
from sage.rings.ring import Algebra, Field
from sage.rings.integer import Integer

from sage.categories.commutative_rings import CommutativeRings
from sage.categories.algebras import Algebras
from sage.categories.fields import Fields

from sage.categories.morphism import Morphism, IdentityMorphism
from sage.rings.morphism import RingHomomorphism
from sage.rings.derivation import RingDerivation
from sage.categories.homset import Hom
from sage.categories.map import Section

from sage.rings.polynomial.polynomial_ring_constructor import PolynomialRing
from sage.rings.polynomial.ore_polynomial_element import OrePolynomialBaseringInjection

WORKING_CENTER_MAX_TRIES = 1000


# Generic implementation of skew polynomial rings
#################################################

class OrePolynomialRing(Algebra, UniqueRepresentation):
    Element = None

    @staticmethod
    def __classcall_private__(cls, base_ring, twist=None, names=None, sparse=False):
        if base_ring not in CommutativeRings():
            raise TypeError('base_ring must be a commutative ring')
        if isinstance(twist, Morphism):
            if (twist.domain() is not base_ring
             or twist.codomain() is not base_ring):
                raise TypeError("the twisting morphism must be an endomorphism of base_ring (=%s)" % base_ring)
            if twist.is_identity():
                morphism = None
            else:
                morphism = twist
            derivation = None
        elif isinstance(twist, RingDerivation):
            if (twist.domain() is not base_ring
             or twist.codomain() is not base_ring):
                raise TypeError("the twisting derivation must be an endomorphism of base_ring (=%s)" % base_ring)
            morphism = twist.parent().twisting_morphism()
            if twist:
                derivation = twist
            else:
                derivation = None
        else:
            raise TypeError("the twist map must be a ring morphism or a derivation")
        if sparse:
            raise NotImplementedError("sparse skew polynomial rings are not implemented")
        if names is None:
            raise TypeError("you must specify the name of the variable")
        try:
            names = normalize_names(1, names)[0]
        except IndexError:
            raise NotImplementedError("multivariate skew polynomials rings not supported")

        # We find the best constructor
        if derivation is None and morphism is None:
            return PolynomialRing(base_ring, names, sparse)

        from sage.rings.polynomial import skew_polynomial_ring
        constructors = [ ]
        if derivation is None:
            if base_ring in Fields():
                try:
                    order = morphism.order()
                    if order is not Infinity:
                        if base_ring.is_finite():
                            constructors.append(skew_polynomial_ring.SkewPolynomialRing_finite_field)
                        else:
                            constructors.append(skew_polynomial_ring.SkewPolynomialRing_finite_order)
                except (AttributeError, NotImplementedError):
                    pass
            constructors.append(skew_polynomial_ring.SkewPolynomialRing)
        
        for constructor in constructors:
            try:
                return constructor(base_ring, morphism, derivation, names, sparse)
            except (AttributeError, NotImplementedError):
                pass

        # We fallback to generic implementation
        return cls.__classcall__(cls, base_ring, morphism, derivation, names, sparse)

    def __init__(self, base_ring, morphism, derivation, name, sparse, category=None):
        r"""
        Initialize ``self``.

        INPUT:

        - ``base_ring`` -- a commutative ring

        - ``twisting_morphism`` -- an automorphism of the base ring

        - ``name`` -- string or list of strings representing the name of
          the variables of ring

        - ``sparse`` -- boolean (default: ``False``)

        - ``category`` -- a category

        EXAMPLES::

            sage: R.<t> = ZZ[]
            sage: sigma = R.hom([t+1])
            sage: S.<x> = SkewPolynomialRing(R,sigma)
            sage: S.category()
            Category of algebras over Univariate Polynomial Ring in t over Integer Ring
            sage: S([1]) + S([-1])
            0
            sage: TestSuite(S).run()
        """
        if self.Element is None:
            self.Element = sage.rings.polynomial.ore_polynomial_element.OrePolynomial_generic_dense
        self.__is_sparse = sparse
        self._morphism = morphism
        self._derivation = derivation
        category = Algebras(base_ring).or_subcategory(category)
        Algebra.__init__(self, base_ring, names=name, normalize=True, category=category)

    def _element_constructor_(self, a=None, check=True, construct=False, **kwds):
        r"""
        Convert a base ring element ``a`` into a constant of this univariate
        skew polynomial ring, possibly non-canonically.

        INPUT:

        - ``a`` -- (default: ``None``) an element of the base ring
          of ``self`` or a ring that has a coerce map from ``self``

        - ``check`` -- boolean (default: ``True``)

        - ``construct`` -- boolean (default: ``False``)

        OUTPUT:

        An zero-degree skew polynomial in ``self``, equal to ``a``.

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

        TESTS::

            sage: S(x, check=True)
            x
        """
        C = self.Element
        if isinstance(a, list):
            return C(self, a, check=check, construct=construct)
        if isinstance(a, sage.structure.element.Element):
            P = a.parent()
            def build(check):
                if a.is_zero():
                    return P.zero()
                else:
                    return C(self, [a], check=check, construct=construct)
            if P is self:
                return a
            elif P is self.base_ring():
                build(False)
            elif P == self.base_ring() or self.base_ring().has_coerce_map_from(P):
                build(True)
        try:
            return a._polynomial_(self)
        except AttributeError:
            pass
        if isinstance(a, str):
            try:
                from sage.misc.parser import Parser, LookupNameMaker
                R = self.base_ring()
                p = Parser(Integer, R, LookupNameMaker({self.variable_name(): self.gen()}, R))
                return self(p.parse(a))
            except NameError:
                raise TypeError("unable to coerce string")
        return C(self, a, check, construct=construct, **kwds)

    def _coerce_map_from_base_ring(self):
        """
        Return a coercion map from the base ring of ``self``.

        EXAMPLES::

            sage: R.<t> = ZZ[]
            sage: S.<x> = SkewPolynomialRing(R, R.hom([t + 1]))
            sage: S.coerce_map_from(R)
            Skew Polynomial base injection morphism:
              From: Univariate Polynomial Ring in t over Integer Ring
              To:   Skew Polynomial Ring in x over Univariate Polynomial Ring in t over Integer Ring twisted by t |--> t + 1
            sage: x.parent()
            Skew Polynomial Ring in x over Univariate Polynomial Ring in t over Integer Ring twisted by t |--> t + 1
            sage: t.parent()
            Univariate Polynomial Ring in t over Integer Ring
            sage: y = x + t  # indirect doctest
            sage: y
            x + t
            sage: y.parent() is S
            True
        """
        return OrePolynomialBaseringInjection(self.base_ring(), self)

    def _coerce_map_from_(self, P):
        r"""
        Check whether ``self`` has a coerce map from ``P``.

        The rings that canonically coerce into this ring are:

        - this ring itself

        - any ring that canonically coerces to the base ring of this ring

        - skew polynomial rings in the same variable and automorphism over
          any base ring that canonically coerces to the base ring of this ring

        INPUT:

        - ``P`` -- a ring

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
        base_ring = self.base_ring()
        try:
            connecting = base_ring.coerce_map_from(P)
            if connecting is not None:
                return self.coerce_map_from(base_ring) * connecting
        except TypeError:
            pass
        if isinstance(P, OrePolynomialRing):
            if self.__is_sparse and not P.is_sparse():
                return False
            if P.variable_name() == self.variable_name():
                return base_ring.has_coerce_map_from(P.base_ring())

    def _repr_(self):
        r"""
        Return a string representation of ``self``.

        EXAMPLES::

            sage: R.<t> = ZZ[]
            sage: sigma = R.hom([t+1])
            sage: S.<x> = SkewPolynomialRing(R,sigma)
            sage: S
            Skew Polynomial Ring in x over Univariate Polynomial Ring in t over Integer Ring twisted by t |--> t + 1
        """
        if self._derivation is None:
            s = "Skew Polynomial Ring in %s over %s twisted by %s" % (self.variable_name(),
                                                                      self.base_ring(),
                                                                      self._morphism._repr_short())
        else:
            s = "Ore Polynomial Ring in %s over %s twisted by " % (self.variable_name(), self.base_ring())
            if self._morphism is not None:
                s += "%s and " % self._morphism._repr_short()
            s += self._derivation._repr_()
        if self.is_sparse():
            s = "Sparse " + s
        return s

    def _latex_(self):
        r"""
        Return a latex representation of ``self``.

        EXAMPLES::

            sage: R.<t> = ZZ[]
            sage: sigma = R.hom([t+1])
            sage: S.<x> = SkewPolynomialRing(R,sigma)
            sage: latex(S)
            \Bold{Z}[t][x,\begin{array}{l}
            \text{\texttt{Ring{ }endomorphism...}}
            \end{array}]
        """
        from sage.misc.latex import latex
        s = "%s\\left[%s" % latex(self.base_ring(), self._latex_variable_names()[0])
        sep = ";"
        if self._morphism is not None:
            s += sep + latex(self._morphism)
            sep = ","
        if self._derivation is not None:
            s += sep + latex(self._derivation)
        return s + "\\right]"

    def change_var(self, var):
        r"""
        Return the skew polynomial ring in variable ``var`` with the same base
        ring and twist map as ``self``.

        INPUT:

        - ``var`` -- a string representing the name of the new variable.

        OUTPUT:

        ``self`` with variable name changed to ``var``.

        EXAMPLES::

            sage: k.<t> = GF(5^3)
            sage: Frob = k.frobenius_endomorphism()
            sage: R.<x> = SkewPolynomialRing(k,Frob); R
            Skew Polynomial Ring in x over Finite Field in t of size 5^3 twisted by t |--> t^5
            sage: Ry = R.change_var('y'); Ry
            Skew Polynomial Ring in y over Finite Field in t of size 5^3 twisted by t |--> t^5
            sage: Ry is R.change_var('y')
            True
        """
        if self._derivation is None:
            twist = self._morphism
        else:
            twist = self._derivation
        return OrePolynomialRing(self.base_ring(), twist, names=var,
                                 sparse=self.__is_sparse)

    def characteristic(self):
        r"""
        Return the characteristic of the base ring of ``self``.

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


    @cached_method
    def twisting_morphism(self, n=1):
        r"""
        Return the twisting endomorphism of the base ring iterated ``n`` times.

        INPUT:

        -  ``n`` - an integer (default: 1)

        EXAMPLES::

            sage: R.<t> = QQ[]
            sage: sigma = R.hom([t+1])
            sage: S.<x> = R['x',sigma]
            sage: S.twisting_morphism()
            Ring endomorphism of Univariate Polynomial Ring in t over Rational Field
              Defn: t |--> t + 1
            sage: S.twisting_morphism() == sigma
            True
            sage: S.twisting_morphism(10)
            Ring endomorphism of Univariate Polynomial Ring in t over Rational Field
              Defn: t |--> t + 10

        If ``n`` in negative, Sage tries to compute the inverse of the
        twist map::

            sage: k.<t> = GF(5^3)
            sage: Frob = k.frobenius_endomorphism()
            sage: T.<y> = k['y',Frob]
            sage: T.twisting_morphism(-1)
            Frobenius endomorphism t |--> t^(5^2) on Finite Field in t of size 5^3

        Sometimes it fails, even if the twist map is actually invertible::

            sage: S.twisting_morphism(-1)
            Traceback (most recent call last):
            ...
            NotImplementedError: inversion of the twisting morphism Ring endomorphism of Univariate Polynomial Ring in t over Rational Field
                  Defn: t |--> t + 1
        """
        if self._morphism is not None:
            try:
                return self._morphism ** n
            except TypeError as e:
                if n < 0:
                    raise NotImplementedError("inversion of the twisting morphism %s" % self._morphism)
                else:
                    raise ValueError("Unexpected error in iterating the twisting morphism: %s", e)

    def twist_map(self, n):
        # Deprecation
        return self.twisting_morphism(n)

    def twisting_derivation(self):
        return self._derivation


    @cached_method
    def gen(self, n=0):
        r"""
        Return the indeterminate generator of this skew polynomial ring.

        INPUT:

        - ``n`` -- index of generator to return (default: 0). Exists for
          compatibility with other polynomial rings.

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
            sage: S.gen(0)
            x

        This is also known as the parameter::

            sage: S.parameter() is S.gen()
            True
        """
        if n != 0:
            raise IndexError("generator %s not defined" % n)
        return self.Element(self, [0,1])

    parameter = gen

    def gens_dict(self):
        r"""
        Return a {name: variable} dictionary of the generators of ``self``.

        EXAMPLES::

            sage: R.<t> = ZZ[]
            sage: sigma = R.hom([t+1])
            sage: S.<x> = SkewPolynomialRing(R,sigma)
            sage: S.gens_dict()
            {'x': x}
        """
        return dict(zip(self.variable_names(), self.gens()))

    def is_finite(self):
        r"""
        Return ``False`` since skew polynomial rings are not finite
        (unless the base ring is `0`.)

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
        return R.is_finite() and R.order() == 1

    def is_exact(self):
        r"""
        Return ``True`` if elements of this skew polynomial ring are exact.
        This happens if and only if elements of the base ring are exact.

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
        r"""
        Return ``True`` if the elements of this polynomial ring are sparsely
        represented.

        .. WARNING::

            Since sparse skew polynomials are not yet implemented, this
            function always returns ``False``.

        EXAMPLES::

            sage: R.<t> = RR[]
            sage: sigma = R.hom([t+1])
            sage: S.<x> = R['x',sigma]
            sage: S.is_sparse()
            False
        """
        return self.__is_sparse

    def ngens(self):
        r"""
        Return the number of generators of this skew polynomial ring,
        which is 1.

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
        Return a random skew polynomial in ``self``.

        INPUT:

        - ``degree`` -- (default: 2) integer with degree
          or a tuple of integers with minimum and maximum degrees

        - ``monic`` -- (default: ``False``) if ``True``, return a monic
          skew polynomial

        - ``*args, **kwds`` -- passed on to the ``random_element`` method
          for the base ring

        OUTPUT:

        Skew polynomial such that the coefficients of `x^i`, for `i` up
        to ``degree``, are random elements from the base ring, randomized
        subject to the arguments ``*args`` and ``**kwds``.

        EXAMPLES::

            sage: k.<t> = GF(5^3)
            sage: Frob = k.frobenius_endomorphism()
            sage: S.<x> = k['x', Frob]
            sage: S.random_element()  # random
            (2*t^2 + 3)*x^2 + (4*t^2 + t + 4)*x + 2*t^2 + 2
            sage: S.random_element(monic=True)  # random
            x^2 + (2*t^2 + t + 1)*x + 3*t^2 + 3*t + 2

        Use ``degree`` to obtain polynomials of higher degree

            sage: p = S.random_element(degree=5)   # random
            (t^2 + 3*t)*x^4 + (4*t + 4)*x^3 + (4*t^2 + 4*t)*x^2 + (2*t^2 + 1)*x + 3

        When ``monic`` is ``False``, the returned skew polynomial may have
        a degree less than ``degree`` (it happens when the random leading
        coefficient is zero). However, if ``monic`` is ``True``, this can't
        happen::

            sage: p = S.random_element(degree=4, monic=True)
            sage: p.leading_coefficient() == S.base_ring().one()
            True
            sage: p.degree() == 4
            True

        If a tuple of two integers is given for the degree argument, a random
        integer will be chosen between the first and second element of the
        tuple as the degree, both inclusive::

            sage: S.random_element(degree=(2,7))  # random
            (3*t^2 + 1)*x^4 + (4*t + 2)*x^3 + (4*t + 1)*x^2
             + (t^2 + 3*t + 3)*x + 3*t^2 + 2*t + 2

        If the first tuple element is greater than the second, a a
        ``ValueError`` is raised::

            sage: S.random_element(degree=(5,4))
            Traceback (most recent call last):
            ...
            ValueError: first degree argument must be less or equal to the second
        """
        R = self.base_ring()
        if isinstance(degree, (list, tuple)):
            if len(degree) != 2:
                raise ValueError("degree argument must be an integer or a tuple of 2 integers (min_degree, max_degree)")
            if degree[0] > degree[1]:
                raise ValueError("first degree argument must be less or equal to the second")
            degree = randint(*degree)
        if monic:
            return self([R.random_element(*args, **kwds) for _ in range(degree)] + [R.one()])
        else:
            return self([R.random_element(*args, **kwds) for _ in range(degree+1)])

    def random_irreducible(self, degree=2, monic=True, *args, **kwds):
        r"""
        Return a random irreducible skew polynomial.

        .. WARNING::

            Elements of this skew polynomial ring need to have a method
            is_irreducible(). Currently, this method is implemented only
            when the base ring is a finite field.

        INPUT:

        -  ``degree`` - Integer with degree (default: 2)
           or a tuple of integers with minimum and maximum degrees

        -  ``monic`` - if True, returns a monic skew polynomial
           (default: True)

        -  ``*args, **kwds`` - Passed on to the ``random_element`` method for
           the base ring

        OUTPUT:

        -  A random skew polynomial

        EXAMPLES::

            sage: k.<t> = GF(5^3)
            sage: Frob = k.frobenius_endomorphism()
            sage: S.<x> = k['x',Frob]
            sage: A = S.random_irreducible(); A
            x^2 + (4*t^2 + 3*t + 4)*x + 4*t^2 + t
            sage: A.is_irreducible()
            True
            sage: B = S.random_irreducible(degree=3,monic=False); B  # random
            (4*t + 1)*x^3 + (t^2 + 3*t + 3)*x^2 + (3*t^2 + 2*t + 2)*x + 3*t^2 + 3*t + 1
            sage: B.is_irreducible()
            True
        """
        if isinstance(degree, (list, tuple)):
            if len(degree) != 2:
                raise ValueError("degree argument must be an integer or a tuple of 2 integers (min_degree, max_degree)")
            if degree[0] > degree[1]:
                raise ValueError("minimum degree must be less or equal than maximum degree")
            degree = randint(*degree)
        while True:
            irred = self.random_element((degree,degree), monic=monic)
            if irred.is_irreducible():
                return irred

    def is_commutative(self):
        r"""
        Return ``True`` if this skew polynomial ring is commutative, i.e. if the
        twist map is the identity.

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
        return self.twisting_morphism().is_identity()
