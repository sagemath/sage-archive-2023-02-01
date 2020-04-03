r"""
Skew Univariate Polynomial Rings

This module provides the :class:`~sage.rings.polynomial.skew_polynomial_ring.SkewPolynomialRing_general`
which constructs a general dense skew univariate polynomials over commutative base rings with
automorphisms over the base rings. This is the set of formal polynomials where the coefficients
are written on the left of the variable of the skew polynomial ring. The modified multiplication
operation over elements of the base ring is extended to all elements of the skew poynomial ring
by associativity and distributivity.

This module also provides :class:`~sage.rings.polynomial.skew_polynomial_ring.SkewPolynomialRing_finite_order`
which is a specialized class for skew polynomial rings over fields equipped with an automorphism of
finite order. It inherits from
:class:`~sage.rings.polynomial.skew_polynomial_ring.SkewPolynomialRing_general` but contains more
methods and provides better algorithms.

AUTHOR:

- Xavier Caruso (2012-06-29): initial version

- Arpit Merchant (2016-08-04): improved docstrings, fixed doctests
  and refactored classes and methods

- Johan Rosenkilde (2016-08-03): changes for bug fixes, docstring and
  doctest errors

"""
# ***************************************************************************
#    Copyright (C) 2012 Xavier Caruso <xavier.caruso@normalesup.org>
#
#    This program is free software: you can redistribute it and/or modify
#    it under the terms of the GNU General Public License as published by
#    the Free Software Foundation, either version 2 of the License, or
#    (at your option) any later version.
#                  https://www.gnu.org/licenses/
# ***************************************************************************

from sage.structure.richcmp import op_EQ
from sage.misc.prandom import randint
from sage.misc.cachefunc import cached_method
from sage.rings.infinity import Infinity
from sage.structure.category_object import normalize_names
import sage.misc.latex as latex

from sage.structure.unique_representation import UniqueRepresentation
from sage.structure.element import Element
from sage.rings.ring import Algebra
from sage.rings.ring import Field
from sage.categories.algebras import Algebras
from sage.rings.integer import Integer

from sage.categories.morphism import Morphism, IdentityMorphism
from sage.rings.morphism import RingHomomorphism
from sage.categories.homset import Hom
from sage.categories.map import Section

from sage.rings.polynomial.polynomial_element import PolynomialBaseringInjection
from sage.rings.polynomial.polynomial_ring import PolynomialRing_general
from sage.rings.polynomial.skew_polynomial_element import SkewPolynomial
from sage.rings.polynomial.skew_polynomial_element import SkewPolynomialBaseringInjection


# Helper functions

def _base_ring_to_fraction_field(S):
    """
    Return the unique skew polynomial ring over the fraction field of
    ``S.base_ring()`` which has ``S`` a sub-ring (internal method).

    INPUT:

    - ``S`` -- a skew polynomial ring.

    OUTPUT:

    - ``Q`` -- the skew polynomial ring over the fraction field of
      ``S.base_ring``.

    EXAMPLES::

        sage: from sage.rings.polynomial.skew_polynomial_ring import _base_ring_to_fraction_field
        sage: R.<t> = ZZ[]
        sage: sigma = R.hom([t+1])
        sage: S.<x> = R['x', sigma]
        sage: _base_ring_to_fraction_field(S)
        Skew Polynomial Ring in x over Fraction Field of Univariate Polynomial Ring in t over Integer Ring twisted by t |-->  t + 1
    """
    R = S.base_ring()
    if isinstance(R, Field):
        return S
    else:
        Q = R.fraction_field()
        gens = R.gens()
        sigmaS = S.twist_map()
        # try:
        sigmaQ = Q.hom([ Q(sigmaS(g)) for g in gens ])
        return Q[S.variable_name(), sigmaQ]
        # except Exception, e:
        #     raise ValueError("unable to lift the twist map to a twist map over %s (error was: %s)" % (Q, e))

def _minimal_vanishing_polynomial(R, eval_pts):
    """
    Return the minimal vanishing polynomial (internal function).

    See the documentation for
    :meth:`SkewPolynomialRing.minimal_vanishing_polynomial` for a description.

    INPUT:

    - ``R`` -- A skew polynomial ring over a field.

    - ``eval_pts`` -- a list of evaluation points

    OUTPUT:

    The minimal vanishing polynomial.

    EXAMPLES::

        sage: from sage.rings.polynomial.skew_polynomial_ring import _minimal_vanishing_polynomial
        sage: k.<t> = GF(5^3)
        sage: Frob = k.frobenius_endomorphism()
        sage: S.<x> = k['x',Frob]
        sage: eval_pts = [1, t, t^2]
        sage: b = _minimal_vanishing_polynomial(S, eval_pts); b
        doctest:...: FutureWarning: This class/method/function is marked as experimental. It, its functionality or its interface might change without a formal deprecation.
        See http://trac.sagemath.org/13215 for details.
        x^3 + 4
    """
    l = len(eval_pts)
    if l == 0:
        return R.one()
    elif l == 1:
        e = eval_pts[0]
        if e.is_zero():
            return R.one()
        else:
            return R.gen() - R.twist_map()(e)/e
    else:
        t = l//2
        A = eval_pts[:t]
        B = eval_pts[t:]
        M_A = _minimal_vanishing_polynomial(R, A)
        B_moved = M_A.multi_point_evaluation(B)
        M_at_B_moved = _minimal_vanishing_polynomial(R, B_moved)
        return M_at_B_moved * M_A


def _lagrange_polynomial(R, eval_pts, values):
    """
    Return the Lagrange polynomial of the given points if it exists.

    Otherwise return an unspecified polynomial (internal method).

    See the documentation for
    :meth:`SkewPolynomialRing.lagrange_polynomial` for a description
    of Lagrange polynomial.

    INPUT:

    - ``R`` -- a skew polynomial ring over a field

    - ``eval_pts`` -- list of evaluation points

    - ``values`` -- list of values that the Lagrange polynomial takes
        at the respective ``eval_pts``

    OUTPUT:

    - the Lagrange polynomial.

    EXAMPLES::

        sage: from sage.rings.polynomial.skew_polynomial_ring import _lagrange_polynomial
        sage: k.<t> = GF(5^3)
        sage: Frob = k.frobenius_endomorphism()
        sage: S.<x> = k['x',Frob]
        sage: eval_pts = [ t , t^2 ]
        sage: values = [ 3*t^2 + 4*t + 4 , 4*t ]
        sage: d = _lagrange_polynomial(S, eval_pts, values); d
        x + t
        sage: d.multi_point_evaluation(eval_pts) == values
        True

    The following restrictions are impossible to satisfy because the evaluation
    points are linearly dependent over the fixed field of the twist map, and the
    corresponding values do not match::

        sage: eval_pts = [ t, 2*t ]
        sage: values = [ 1, 3 ]
        sage: _lagrange_polynomial(S, eval_pts, values)
        Traceback (most recent call last):
        ...
        ValueError: the given evaluation points are linearly dependent over the fixed field of the twist map,
        so a Lagrange polynomial could not be determined (and might not exist).
    """
    l = len(eval_pts)
    if l == 1:
        if eval_pts[0].is_zero():
            # This is due to linear dependence among the eval_pts.
            raise ValueError("the given evaluation points are linearly dependent over the fixed field of the twist map, so a Lagrange polynomial could not be determined (and might not exist).")
        return (values[0]/eval_pts[0])*R.one()
    else:
        t = l//2
        A = eval_pts[:t]
        B = eval_pts[t:]
        M_A = _minimal_vanishing_polynomial(R, A)
        M_B = _minimal_vanishing_polynomial(R, B)
        A_ = M_B.multi_point_evaluation(A)
        B_ = M_A.multi_point_evaluation(B)
        I_1 = _lagrange_polynomial(R, A_, values[:t])
        I_2 = _lagrange_polynomial(R, B_, values[t:])
        return I_1 * M_B + I_2 * M_A



# Generic implementation of skew polynomial rings
#################################################

class SkewPolynomialRing_general(Algebra, UniqueRepresentation):
    r"""
    A general implementation of univariate skew polynomialring over a commutative ring.

    Let `R` be a commutative ring, and let `\sigma` be an automorphism of
    `R`.  The ring of skew polynomials `R[X, \sigma]` is the polynomial
    ring `R[X]`, where the addition is the usual polynomial addition, but
    the multiplication operation is defined by the modified rule

    .. MATH::

        X*a = \sigma(a) X.

    This means that `R[X, \sigma]` is a non-commutative ring. Skew polynomials
    were first introduced by Ore [Ore1933]_.

    EXAMPLES::

        sage: R.<t> = ZZ[]
        sage: sigma = R.hom([t+1])
        sage: S.<x> = SkewPolynomialRing(R,sigma); S
        Skew Polynomial Ring in x over Univariate Polynomial Ring in t over Integer Ring
         twisted by t |--> t + 1

    One can also use a shorter syntax::

        sage: S.<x> = R['x',sigma]; S
        Skew Polynomial Ring in x over Univariate Polynomial Ring in t over Integer Ring
         twisted by t |--> t + 1

    If we omit the diamond notation, the variable holding the indeterminate is
    not assigned::

        sage: Sy = R['y',sigma]
        sage: y
        Traceback (most recent call last):
        ...
        NameError: name 'y' is not defined
        sage: Sy.gen()
        y

    Note however that contrary to usual polynomial rings, we cannot omit the
    variable name on the RHS, since this collides with the notation for creating polynomial rings::

        sage: Sz.<z> = R[sigma]
        Traceback (most recent call last):
        ...
        ValueError: variable name 'Ring endomorphism of Univariate Polynomial Ring in t over Integer Ring\n
            Defn: t |--> t + 1' is not alphanumeric

    Of course, skew polynomial rings with different twist maps are not
    equal either::

        sage: R['x',sigma] == R['x',sigma^2]
        False

    Saving and loading of polynomial rings works::

        sage: loads(dumps(R['x',sigma])) == R['x',sigma]
        True

    .. SEEALSO::

        :meth:`sage.rings.polynomial.skew_polynomial_ring_constructor.SkewPolynomialRing`
        :mod:`sage.rings.polynomial.skew_polynomial_element`
    """
    @staticmethod
    def __classcall__(cls, base_ring, twist_map=None, name=None, sparse=False,
                      element_class=None):
        r"""
        Set the default values for ``name``, ``sparse`` and ``element_class``.

        EXAMPLES::

            sage: R.<t> = ZZ[]
            sage: sigma = R.hom([t+1])
            sage: S.<x> = SkewPolynomialRing(R, sigma)
            sage: S.__class__(R, sigma, x)
            Skew Polynomial Ring in x over Univariate Polynomial Ring in t over Integer Ring
             twisted by t |--> t + 1
        """
        if not element_class:
            if sparse:
                raise NotImplementedError("sparse skew polynomials are not implemented")
            else:
                from sage.rings.polynomial import skew_polynomial_element
                element_class = skew_polynomial_element.SkewPolynomial_generic_dense
        if twist_map is None:
            twist_map = IdentityMorphism(base_ring)
        else:
            if not isinstance(twist_map, Morphism):
                raise TypeError("given map is not a ring homomorphism")
            if twist_map.domain() != base_ring or twist_map.codomain() != base_ring:
                raise TypeError("given map is not an automorphism of %s" % base_ring)
        return super(SkewPolynomialRing_general,cls).__classcall__(cls,
                         base_ring, twist_map, name, sparse, element_class)

    def __init__(self, base_ring, twist_map, name, sparse, element_class):
        r"""
        Initialize ``self``.

        INPUT:

        - ``base_ring`` -- a commutative ring

        - ``twist_map`` -- an automorphism of the base ring

        - ``name`` -- string or list of strings representing the name of
          the variables of ring

        - ``sparse`` -- boolean (default: ``False``)

        - ``element_class`` -- class representing the type of element to
          be used in ring

        EXAMPLES::

            sage: R.<t> = ZZ[]
            sage: sigma = R.hom([t+1])
            sage: S.<x> = SkewPolynomialRing(R,sigma)
            sage: S.category()
            Category of algebras over (unique factorization domains and commutative algebras over (euclidean domains and infinite enumerated sets and metric spaces) and infinite sets)
            sage: S([1]) + S([-1])
            0
            sage: TestSuite(S).run()

            sage: k.<t> = GF(5^3)
            sage: Frob = k.frobenius_endomorphism()
            sage: T.<x> = k['x', Frob]; T
            Skew Polynomial Ring in x over Finite Field in t of size 5^3 twisted by t |--> t^5
            sage: T.category()
            Category of algebras over finite enumerated fields

        We skip the pickling tests currently because ``Frob`` does not
        pickle correctly (see note on :trac:`13215`)::

            sage: TestSuite(T).run(skip=["_test_pickling", "_test_elements"])
        """
        self.__is_sparse = sparse
        self._polynomial_class = element_class
        self._map = twist_map
        self._maps = {0: IdentityMorphism(base_ring), 1: self._map}
        Algebra.__init__(self, base_ring, names=name, normalize=True,
                         category=Algebras(base_ring.category()))

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
        C = self._polynomial_class
        if isinstance(a, list):
            return C(self, a, check=check, construct=construct)
        if isinstance(a, Element):
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
            sage: y = x + t; y
            x + t
            sage: y.parent() is S
            True
        """
        return SkewPolynomialBaseringInjection(self.base_ring(), self)

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
        if isinstance(P, SkewPolynomialRing_general):
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
        s = "Skew Polynomial Ring in %s over %s twisted by %s" % (self.variable_name(),
                                                                  self.base_ring(),
                                                                  self._map._repr_short())
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
        return "%s[%s,%s]" % (latex(self.base_ring()), self.latex_variable_names()[0],
                              latex(self._map))

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
        from sage.rings.polynomial.skew_polynomial_ring_constructor import SkewPolynomialRing
        return SkewPolynomialRing(self.base_ring(), self._map, names=var,
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
    def twist_map(self, n=1):
        r"""
        Return the twist map, the automorphism of the base ring of
        ``self``, iterated ``n`` times.

        INPUT:

        -  ``n`` - an integer (default: 1)

        OUTPUT:

        ``n``-th iterative of the twist map of this skew polynomial ring.

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

        If ``n`` in negative, Sage tries to compute the inverse of the
        twist map::

            sage: k.<t> = GF(5^3)
            sage: Frob = k.frobenius_endomorphism()
            sage: T.<y> = k['y',Frob]
            sage: T.twist_map(-1)
            Frobenius endomorphism t |--> t^(5^2) on Finite Field in t of size 5^3

        Sometimes it fails, even if the twist map is actually invertible::

            sage: S.twist_map(-1)
            Traceback (most recent call last):
            ...
            NotImplementedError: inversion of the twist map Ring endomorphism of Univariate Polynomial Ring in t over Rational Field
                  Defn: t |--> t + 1
        """
        try:
            return self._map ** n
        except TypeError as e:
            if n < 0:
                raise NotImplementedError("inversion of the twist map %s" % self._map)
            else:
                raise ValueError("Unexpected error in iterating the twist map: %s", e)


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
        return self._polynomial_class(self, [0,1])

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
        return self.twist_map().is_identity()

    def minimal_vanishing_polynomial(self, eval_pts):
        """
        Return the minimal-degree, monic skew polynomial which vanishes at all
        the given evaluation points.

        The degree of the vanishing polynomial is at most the length of
        ``eval_pts``. Equality holds if and only if the elements of ``eval_pts``
        are linearly independent over the fixed field of ``self.twist_map()``.

        - ``eval_pts`` -- list of evaluation points which are linearly
          independent over the fixed field of the twist map of the associated
          skew polynomial ring

        OUTPUT:

        The minimal vanishing polynomial.

        EXAMPLES::

            sage: k.<t> = GF(5^3)
            sage: Frob = k.frobenius_endomorphism()
            sage: S.<x> = k['x',Frob]
            sage: eval_pts = [1, t, t^2]
            sage: b = S.minimal_vanishing_polynomial(eval_pts); b
            x^3 + 4

        The minimal vanishing polynomial evaluates to 0 at each of the evaluation points::

            sage: eval = b.multi_point_evaluation(eval_pts); eval
            [0, 0, 0]

        If the evaluation points are linearly dependent over the fixed field of
        the twist map, then the returned polynomial has lower degree than the
        number of evaluation points::

            sage: S.minimal_vanishing_polynomial([t])
            x + 3*t^2 + 3*t
            sage: S.minimal_vanishing_polynomial([t, 3*t])
            x + 3*t^2 + 3*t
        """
        return _minimal_vanishing_polynomial(_base_ring_to_fraction_field(self), eval_pts)

    def lagrange_polynomial(self, points):
        r"""
        Return the minimal-degree polynomial which interpolates the given
        points.

        More precisely, given `n` pairs `(x_1, y_1), ..., (x_n, y_n) \in R^2`,
        where `R` is ``self.base_ring()``, compute a skew polynomial `p(x)` such
        that `p(x_i) = y_i` for each `i`, under the condition that the `x_i` are
        linearly independent over the fixed field of ``self.twist_map()``.

        If the `x_i` are linearly independent over the fixed field of
        ``self.twist_map()`` then such a polynomial is guaranteed to exist.
        Otherwise, it might exist depending on the `y_i`, but the algorithm used
        in this implementation does not support that, and so an error is always
        raised.

        INPUT:

        - ``points`` -- a list of pairs ``(x_1, y_1),..., (x_n, y_n)`` of
          elements of the base ring of ``self``. The `x_i` should be linearly
          independent over the fixed field of ``self.twist_map()``.

        OUTPUT:

        The Lagrange polynomial.

        EXAMPLES::

            sage: k.<t> = GF(5^3)
            sage: Frob = k.frobenius_endomorphism()
            sage: S.<x> = k['x',Frob]
            sage: points = [(t, 3*t^2 + 4*t + 4), (t^2, 4*t)]
            sage: d = S.lagrange_polynomial(points); d
            x + t

            sage: R.<t> = ZZ[]
            sage: sigma = R.hom([t+1])
            sage: T.<x> = R['x', sigma]
            sage: points = [ (1, t^2 + 3*t + 4), (t, 2*t^2 + 3*t + 1), (t^2, t^2 + 3*t + 4) ]
            sage: p = T.lagrange_polynomial(points); p
            ((-t^4 - 2*t - 3)/-2)*x^2 + (-t^4 - t^3 - t^2 - 3*t - 2)*x + (-t^4 - 2*t^3 - 4*t^2 - 10*t - 9)/-2
            sage: p.multi_point_evaluation([1, t, t^2]) == [ t^2 + 3*t + 4, 2*t^2 + 3*t + 1, t^2 + 3*t + 4 ]
            True

        If the `x_i` are linearly dependent over the fixed field of
        ``self.twist_map()``, then an error is raised::

            sage: T.lagrange_polynomial([ (t, 1), (2*t, 3) ])
            Traceback (most recent call last):
            ...
            ValueError: the given evaluation points are linearly dependent over the fixed field of the twist map,
            so a Lagrange polynomial could not be determined (and might not exist).
        """
        l = len(points)
        if not all( len(pair) == 2 for pair in points ):
            raise TypeError("supplied points must be pairs of elements of base ring")
        eval_pts = [ x for (x,_) in points ]
        values   = [ y for (_,y) in points ]

        if l > len(set(eval_pts)):
            raise TypeError("the evaluation points must be distinct")
        zero_i = [ i for i in range(l) if eval_pts[i].is_zero() ]
        if zero_i and not values[zero_i[0]].is_zero():
            raise TypeError("a skew polynomial always evaluates to 0 at 0, but a non-zero value was requested.")

        return _lagrange_polynomial(_base_ring_to_fraction_field(self), eval_pts, values)


# Special classes for twisting morphisms with finite order
##########################################################

class SectionSkewPolynomialCenterInjection(Section):
    r"""
    Section of the canonical injection of the centre of a skew
    polynomial ring into this ring

    TESTS::

        sage: k.<a> = GF(5^3)
        sage: S.<x> = SkewPolynomialRing(k, k.frobenius_endomorphism())
        sage: Z = S.centre()
        sage: iota = S.coerce_map_from(Z)
        sage: sigma = iota.section()
        sage: TestSuite(sigma).run()
    """
    def _call_ (self, x):
        r"""
        Return `x` viewed as an element of the centre

        EXAMPLES::

            sage: k.<a> = GF(5^3)
            sage: S.<x> = SkewPolynomialRing(k, k.frobenius_endomorphism())
            sage: Z = S.centre()
            sage: iota = S.coerce_map_from(Z)
            sage: sigma = iota.section()
            sage: sigma(x^3)
            (x^3)
            sage: sigma(x^2)
            Traceback (most recent call last):
            ...
            ValueError: x^2 is not in the center
        """
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

    def _richcmp_(self, right, op):
        r"""
        Compare this morphism with ``right``

        TESTS::

            sage: k.<a> = GF(5^3)
            sage: S.<x> = SkewPolynomialRing(k, k.frobenius_endomorphism())
            sage: Z = S.centre()
            sage: iota = S.coerce_map_from(Z)
            sage: sigma = iota.section()
 
            sage: s = loads(dumps(sigma))
            sage: s == sigma
            True
            sage: s is sigma
            False
        """
        if op == op_EQ:
            return (self.domain() is right.domain()) and (self.codomain() is right.codomain()) 
        return NotImplemented


class SkewPolynomialCenterInjection(RingHomomorphism):
    r"""
    Canonical injection of the centre of a skew polynomial ring
    into this ring

    TESTS::

        sage: k.<a> = GF(5^3)
        sage: S.<x> = SkewPolynomialRing(k, k.frobenius_endomorphism())
        sage: Z = S.centre()
        sage: iota = S.coerce_map_from(Z)
        sage: TestSuite(iota).run()
    """
    def __init__(self,domain,codomain,embed,order):
        r"""
        Initialize this morphism

        EXAMPLES::

            sage: k.<a> = GF(5^3)
            sage: S.<x> = SkewPolynomialRing(k, k.frobenius_endomorphism())
            sage: Z = S.centre()
            sage: S.coerce_map_from(Z)   # indirect doctest
            Embedding of the center of Skew Polynomial Ring in x over Finite Field in a of size 5^3 twisted by a |--> a^5 into this ring
        """
        RingHomomorphism.__init__(self,Hom(domain,codomain))
        self._embed = embed
        self._order = order
        self._codomain = codomain
        self._section = SectionSkewPolynomialCenterInjection(self)

    def _repr_(self):
        r"""
        Return a string representation of this morphism

        EXAMPLES::

            sage: k.<a> = GF(5^3)
            sage: S.<x> = SkewPolynomialRing(k, k.frobenius_endomorphism())
            sage: Z = S.centre()
            sage: iota = S.coerce_map_from(Z)
            sage: iota
            Embedding of the center of Skew Polynomial Ring in x over Finite Field in a of size 5^3 twisted by a |--> a^5 into this ring
            sage: iota._repr_()
            'Embedding of the center of Skew Polynomial Ring in x over Finite Field in a of size 5^3 twisted by a |--> a^5 into this ring'
        """
        return "Embedding of the center of %s into this ring" % self._codomain

    def _call_(self,x):
        r"""
        Return the image of `x` by this morphism

        TESTS::

            sage: k.<a> = GF(5^3)
            sage: S.<x> = SkewPolynomialRing(k, k.frobenius_endomorphism())
            sage: Z = S.centre()
            sage: iota = S.coerce_map_from(Z)

            sage: x3 = Z.gen(); x3
            (x^3)
            sage: iota(x3)
            x^3
        """
        k = self._codomain.base_ring ()
        l = [ ]
        lz = [ k(0) ] * (self._order-1)
        for c in x.list():
            l += [ self._embed(c) ] + lz
        return self._codomain (l)

    def _richcmp_(self, right, op):
        r"""
        Compare this morphism with ``right``

        TESTS::

            sage: k.<a> = GF(5^3)
            sage: S.<x> = SkewPolynomialRing(k, k.frobenius_endomorphism())
            sage: Z = S.centre()
            sage: iota = S.coerce_map_from(Z)
 
            sage: i = loads(dumps(iota))
            sage: i == iota
            True
            sage: i is iota
            False
        """
        if op == op_EQ:
            return (self.domain() is right.domain()) and (self.codomain() is right.codomain()) 
        return NotImplemented

    def section(self):
        r"""
        Return a section of this morphism

        EXAMPLES::

            sage: k.<a> = GF(5^3)
            sage: S.<x> = SkewPolynomialRing(k, k.frobenius_endomorphism())
            sage: Z = S.centre()
            sage: iota = S.coerce_map_from(Z)
            sage: sigma = iota.section()
            sage: sigma(x^3)
            (x^3)
        """
        return self._section


class CenterSkewPolynomialRing(PolynomialRing_general, UniqueRepresentation):
    """
    A class for the centre of a skew polynomial ring

    TESTS::

        sage: k.<a> = GF(5^3)
        sage: S.<x> = SkewPolynomialRing(k, k.frobenius_endomorphism())
        sage: Z = S.centre()
        sage: TestSuite(Z).run()  
    """
    def __init__ (self, skew_ring, names=None, sparse=False):
        r"""
        Initialize this ring

        EXAMPLES::

            sage: k.<a> = GF(5^3)
            sage: S.<x> = SkewPolynomialRing(k, k.frobenius_endomorphism())
            sage: S.centre()  # indirect doctest
            Center of Skew Polynomial Ring in x over Finite Field in a of size 5^3 twisted by a |--> a^5:
            Univariate Polynomial Ring in (x^3) over Finite Field of size 5
        """
        if not isinstance (skew_ring, SkewPolynomialRing_general):
            raise TypeError("%s is not a Skew Polynomial Ring" % skew_ring)
        self._skew_ring = skew_ring
        base_ring = skew_ring.base_ring()
        kfixed, embed = skew_ring._map.fixed_field()
        self._embed_basering = embed
        order = skew_ring._map.order()
        if order is Infinity:
            raise NotImplementedError
        self.__is_sparse = sparse
        self._PolynomialRing_general__is_sparse = sparse
        if sparse:
            raise NotImplementedError("sparse skew polynomials are not implemented")
        else:
            from sage.rings.polynomial.skew_polynomial_element import CenterSkewPolynomial_generic_dense
            self.Element = self._polynomial_class = CenterSkewPolynomial_generic_dense

        # Algebra.__init__ also calls __init_extra__ of Algebras(...).parent_class, which
        # tries to provide a conversion from the base ring, if it does not exist.
        # This is for algebras that only do the generic stuff in their initialisation.
        # But here, we want to use PolynomialBaseringInjection. Hence, we need to
        # wipe the memory and construct the conversion from scratch.
        Algebra.__init__(self, kfixed, names=names, normalize=True, category=None)

        if names is None:
            self._pickling_variable_name = None
            if order == 1:
                self._variable_name = skew_ring.variable_name()
                self._latex_variable_name = skew_ring.latex_variable_names()[0]
                self._parenthesis = False
            else:
                self._variable_name = skew_ring.variable_name() + "^" + str(order)
                self._latex_variable_name = skew_ring.latex_variable_names()[0] + "^{" + str (order) + "}"
                self._parenthesis = True
        else:
            self._pickling_variable_name = self._variable_name = Algebra.variable_name(self)
            self._latex_variable_name = Algebra.latex_variable_names(self)[0]
            self._parenthesis = False
        self._names = [ self._variable_name ]
        self.__generator = self._polynomial_class(self, [0,1], is_gen=True)
        base_inject = PolynomialBaseringInjection(kfixed, self)
        center_inject = SkewPolynomialCenterInjection (self, skew_ring, embed, order)
        self._unset_coercions_used()
        self._populate_coercion_lists_(
            coerce_list = [base_inject],
            convert_list = [list, base_inject],
            embedding = center_inject)

    def __reduce__(self):
        r"""
        Return a function and a list of arguments that reconstruct
        this object

        TESTS::

            sage: k.<a> = GF(5^3)
            sage: Frob = k.frobenius_endomorphism()
            sage: S.<x> = SkewPolynomialRing(k, Frob)
            sage: Z = S.centre()
            sage: Z.__reduce__()
            (<class 'sage.rings.polynomial.skew_polynomial_ring.CenterSkewPolynomialRing'>,
             (Skew Polynomial Ring in x over Finite Field in a of size 5^3 twisted by a |--> a^5,
              None,
              False))

            sage: loads(dumps(Z)) is Z
            True
        """
        return CenterSkewPolynomialRing, (self._skew_ring, self._pickling_variable_name, self.is_sparse())

    def _repr_ (self):
        r"""
        Return a string representation of this ring

        EXAMPLES::

            sage: k.<a> = GF(5^3)
            sage: S.<x> = SkewPolynomialRing(k, k.frobenius_endomorphism())
            sage: Z = S.centre()
            sage: Z
            Center of Skew Polynomial Ring in x over Finite Field in a of size 5^3 twisted by a |--> a^5:
            Univariate Polynomial Ring in (x^3) over Finite Field of size 5

            sage: Z._repr_()
            'Center of Skew Polynomial Ring in x over Finite Field in a of size 5^3 twisted by a |--> a^5:\nUnivariate Polynomial Ring in (x^3) over Finite Field of size 5'
        """
        s = "Center of %s:\n" % self._skew_ring
        s += PolynomialRing_general._repr_(self)
        return s

    def _latex_ (self):
        r"""
        Return a latex representation of this ring

        EXAMPLES::

            sage: k.<a> = GF(5^3)
            sage: S.<x> = SkewPolynomialRing(k, k.frobenius_endomorphism())
            sage: Z = S.centre()
            sage: Z
            Center of Skew Polynomial Ring in x over Finite Field in a of size 5^3 twisted by a |--> a^5:
            Univariate Polynomial Ring in (x^3) over Finite Field of size 5

            sage: Z._latex_()
            '\\Bold{F}_{5}[x^{3}]'
        """
        return "%s[%s]"%(latex.latex(self.base_ring()), self._latex_variable_name)

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


class SkewPolynomialRing_finite_order(SkewPolynomialRing_general):
    """
    A specialized class for skew polynomial rings over finite fields.

    .. SEEALSO::

        :meth:`sage.rings.polynomial.skew_polynomial_ring_constructor.SkewPolynomialRing`
        :class:`sage.rings.polynomial.skew_polynomial_ring.SkewPolynomialRing_general`
        :mod:`sage.rings.polynomial.skew_polynomial_finite_order`
    """
    @staticmethod
    def __classcall__(cls, base_ring, map, name=None, sparse=False, element_class=None):
        if not element_class:
            if sparse:
                raise NotImplementedError("sparse skew polynomials are not implemented")
            else:
                from sage.rings.polynomial import skew_polynomial_finite_order
                element_class = skew_polynomial_finite_order.SkewPolynomial_finite_order_dense
        return super(SkewPolynomialRing_general,cls).__classcall__(cls, base_ring, map, name, sparse, element_class)

    def __init__(self, base_ring, twist_map, name, sparse, element_class):
        SkewPolynomialRing_general.__init__(self, base_ring, twist_map, name, sparse, element_class)
        self._order = twist_map.order()

    def center(self, names=None, name=None):
        r"""
        Return the center of this skew polynomial ring.

        .. NOTE::

            If F denotes the subring of R fixed by `\sigma`, the center of
            `R[X,\sigma]` is `F` if `\sigma` has infinite order and `F[X^r]`
            if `\sigma` has finite order `r`.

        .. WARNING::

            This function assumes that `\sigma` has a method order() (which
            returns its order) and a method fixed_field() (which returns
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

            sage: a = t*x^2 + t*x + 2
            sage: b = a.reduced_norm(); b
            2*(x^3)^2 + (x^3) + 3
            sage: c = a + b; c
            2*x^6 + x^3 + t*x^2 + t*x
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
        if name is not None:
            name = normalize_names(1,name)[0]
        return CenterSkewPolynomialRing(self, name, self.is_sparse())

    def centre(self, names=None, name=None):
        r"""
        Return the centre of this skew polynomial ring.

        .. NOTE::

            If F denotes the subring of R fixed by `\sigma`, the centre of
            `R[X,\sigma]` is `F` if `\sigma` has infinite order and `F[X^r]`
            if `\sigma` has finite order `r`.

        .. WARNING::

            This function assumes that `\sigma` has a method order() (which
            returns its order) and a method fixed_field() (which returns
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
