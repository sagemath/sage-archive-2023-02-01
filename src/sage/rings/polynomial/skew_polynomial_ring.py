r"""
Skew Univariate Polynomial Rings

Sage implements dense skew univariate polynomials over commutative rings.

DEFINITION:

Given a ring `R` and a ring endomorphism `\sigma` of `R`, the ring of
skew polynomials `R[x,\sigma]` is the usual abelian group polynomial
`R[x]` equipped with the modification multiplication deduced from the
rule `X a = \sigma(a) X`.

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
    ValueError: variable name 'Ring endomorphism of Univariate Polynomial Ring in t over Integer Ring\n  Defn: t |--> t + 1' is not alphanumeric

As for polynomials, skew polynomial rings with different variable names
are not equal::

    sage: R['x',sigma] == R['y',sigma]
    False

Of course, skew polynomial rings with different twist maps are not
equal as well

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

AUTHOR:

- Xavier Caruso (2012-06-29)
"""

#############################################################################
#    Copyright (C) 2012 Xavier Caruso <xavier.caruso@normalesup.org>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#
#                  http://www.gnu.org/licenses/
#****************************************************************************

from sage.structure.unique_representation import UniqueRepresentation

from sage.rings.infinity import Infinity

from sage.structure.element import Element
import sage.algebras.algebra
import sage.categories.basic as categories
import sage.rings.ring as ring
from sage.structure.element import RingElement
from sage.rings.ring import is_Ring, PrincipalIdealDomain
import sage.rings.polynomial.polynomial_element_generic as polynomial_element_generic
import sage.rings.rational_field as rational_field
from sage.rings.integer_ring import is_IntegerRing, IntegerRing
from sage.rings.integer import Integer
from sage.libs.pari.all import pari_gen

from sage.structure.category_object import normalize_names
from sage.misc.prandom import randint

from sage.categories.morphism import Morphism
from sage.categories.morphism import IdentityMorphism
from sage.categories.homset import Hom
from sage.categories.map import Section
from sage.rings.morphism import RingHomomorphism

from sage.matrix.matrix_space import MatrixSpace

import sage.misc.latex as latex

from sage.rings.polynomial.polynomial_ring import PolynomialRing_general
from sage.rings.polynomial.polynomial_element import Polynomial_generic_dense
from sage.rings.polynomial.polynomial_element import PolynomialBaseringInjection
from sage.rings.polynomial.skew_polynomial_element import SkewPolynomial


def is_SkewPolynomialRing(S):
    """
    Return True if S is a skew polynomial ring.

    EXAMPLES::

        sage: from sage.rings.polynomial.skew_polynomial_ring import is_SkewPolynomialRing
        sage: is_SkewPolynomialRing(QQ['x'])
        False

        sage: k.<t> = GF(5^3)
        sage: Frob = k.frobenius_endomorphism()
        sage: is_SkewPolynomialRing(k['x',Frob])
        True
    """
    return isinstance(S, SkewPolynomialRing_general)


#########################################################################################

class SkewPolynomialRing_general(sage.algebras.algebra.Algebra,UniqueRepresentation):
    """
    Skew Univariate polynomial ring over a ring.
    """
    Element = SkewPolynomial
    @staticmethod
    def __classcall__(cls, base_ring, map, name=None, sparse=False, element_class=None):
        if not element_class:
            if sparse:
                raise NotImplementedError("sparse skew polynomials are not implemented")
            else:
                from sage.rings.polynomial import skew_polynomial_element
                element_class = skew_polynomial_element.SkewPolynomial_generic_dense
        return super(SkewPolynomialRing_general,cls).__classcall__(cls,base_ring,map,name,sparse,element_class)

    def __init__(self, base_ring, map, name, sparse, element_class):
        """
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
                    raise TypeError("The given map is not an endomorphism of %s" % base_ring)
            else:
                raise TypeError("The given map is not a ring homomorphism")
        self._maps = { 0:IdentityMorphism(base_ring), 1:self._map }
        self._center = { }
        self._center_variable = None
        # Algebra.__init__ also calls __init_extra__ of Algebras(...).parent_class, which
        # tries to provide a conversion from the base ring, if it does not exist.
        # This is for algebras that only do the generic stuff in their initialisation.
        # But here, we want to use PolynomialBaseringInjection. Hence, we need to
        # wipe the memory and construct the conversion from scratch.
        sage.algebras.algebra.Algebra.__init__(self, base_ring, names=name, normalize=True, category=category)
        self.__generator = self._polynomial_class(self, [0,1], is_gen=True)
        base_inject = sage.rings.polynomial.skew_polynomial_element.SkewPolynomialBaseringInjection(base_ring,self)
        self._unset_coercions_used()
        self._populate_coercion_lists_(
                coerce_list = [base_inject],
                convert_list = [list, base_inject])

    def __reduce__(self):
        import sage.rings.polynomial.skew_polynomial_ring_constructor
        return (sage.rings.polynomial.skew_polynomial_ring_constructor.SkewPolynomialRing,
                (self.base_ring(), self.twist_map(), self.variable_name(), None, self.is_sparse()))


    def _element_constructor_(self, x=None, check=True, is_gen = False, construct=False, **kwds):
        """
        Convert ``x`` into this univariate skew polynomial ring,
        possibly non-canonically.

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
        """
        C = self._polynomial_class
        if isinstance(x, list):
            return C(self, x, check=check, is_gen=False,construct=construct)
        if isinstance(x, Element):
            P = x.parent()
            if P is self:
                return x
            elif P is self.base_ring():
                # It *is* the base ring, hence, we should not need to check.
                # Moreover, if x is equal to zero then we usually need to
                # provide [] to the polynomial class, not [x], if we don't want
                # to check (normally, polynomials like to strip trailing zeroes).
                # However, in the padic case, we WANT that trailing
                # zeroes are not stripped, because O(5)==0, but still it must
                # not be forgotten. It should be the job of the __init__ method
                # to decide whether to strip or not to strip.
                return C(self, [x], check=False, is_gen=False,
                         construct=construct)
            elif P == self.base_ring():
                return C(self, [x], check=True, is_gen=False,
                         construct=construct)

            elif self.base_ring().has_coerce_map_from(P):
                return C(self, [x], check=True, is_gen=False,
                        construct=construct)
        try: #if hasattr(x, '_polynomial_'):
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
                raise TypeError("Unable to coerce string")
        return C(self, x, check, is_gen, construct=construct, **kwds)

    #def construction(self):
    #    from sage.categories.pushout import PolynomialFunctor
    #    return PolynomialFunctor(self.variable_name(), sparse=self.__is_sparse), self.base_ring()

    def _coerce_map_from_(self, P):
        """
        The rings that canonically coerce into this ring are:

        - this ring itself

        - any ring that canonically coerces to the base ring of this ring

        - skew polynomial rings in the same variable and automorphism over
          any base ring that canonically coerces to the base ring of this ring

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
                        Ring morphism:
                        From: Univariate Polynomial Ring in t over Integer Ring
                        To:   Skew Polynomial Ring in x over Univariate Polynomial Ring in t over Integer Ring twisted by t |--> t + 1
            sage: S.coerce_map_from(S)
            Identity endomorphism of Skew Polynomial Ring in x over Univariate Polynomial Ring in t over Integer Ring twisted by t |--> t + 1
        """
        # handle constants that canonically coerce into self.base_ring()
        # first, if possible
        try:
            connecting = self.base_ring().coerce_map_from(P)
            if connecting is not None:
                return self.coerce_map_from(self.base_ring()) * connecting
        except TypeError:
            pass

        # skew polynomial rings in the same variable over a base that canonically
        # coerces into self.base_ring()
        try:
            if is_SkewPolynomialRing(P):
                if self.__is_sparse and not P.is_sparse():
                    return False
                if P.variable_name() == self.variable_name():
                    if P.base_ring() is self.base_ring() and \
                            self.base_ring() is ZZ_sage:
                        # We're trying to coerce from FLINT->NTL
                        # or vice versa.  Only allow coercions from
                        # NTL->FLINT, not vice versa.
                        # Unfortunately this doesn't work, because
                        # the parents for ZZ[x]-with-NTL and
                        # ZZ[x]-with-FLINT are equal, and the coercion model
                        # believes this means that both coercions are valid;
                        # but we'll probably change that in the
                        # coercion model, at which point this code will
                        # become useful.
                        if self._implementation_names == ('NTL',):
                            return False
                    return self.base_ring().has_coerce_map_from(P.base_ring())
        except AttributeError:
            pass

    #def _sage_input_(self, sib, coerced):
    #    r"""
    #    Produce an expression which will reproduce this value when
    #    evaluated.
    #
    #    EXAMPLES::
    #
    #        sage: sage_input(GF(5)['x']['y'], verify=True)
    #        # Verified
    #        GF(5)['x']['y']
    #        sage: from sage.misc.sage_input import SageInputBuilder
    #        sage: ZZ['z']._sage_input_(SageInputBuilder(), False)
    #        {constr_parent: {subscr: {atomic:ZZ}[{atomic:'z'}]} with gens: ('z',)}
    #    """
    #    base = sib(self.base_ring())
    #    sie = base[self.variable_name()]
    #    gens_syntax = sib.empty_subscript(base)
    #    return sib.parent_with_gens(self, sie, self.variable_names(), 'R',
    #                                gens_syntax=gens_syntax)

    # Skew Polynomial rings should be unique parents. Hence,
    # no need for __cmp__. Or actually, having a __cmp__
    # method that identifies a dense with a sparse ring
    # is a bad bad idea!
    #def __cmp__(left, right):
    #    c = cmp(type(left),type(right))
    #    if c: return c
    #    return cmp((left.base_ring(), left.variable_name(), left.twist_map(), left.is_sparse()),
    #               (right.base_ring(), right.variable_name(), right.twist_map(), right.is_sparse()))

    def _repr_(self):
        """
        Return representation of self.
        
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
        Return latex representation of self.

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
        Return the skew polynomial ring in variable var over the same base
        ring.

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
        Return the characteristic of the base field of this skew polynomial ring.

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

    def twist_map(self,n=1):
        """
        Return the twist map (eventually iterated several times) used to define
        this skew polynomial ring.

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
            NotImplementedError: Inversion of the twist map Ring endomorphism of Univariate Polynomial Ring in t over Rational Field
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
                    raise NotImplementedError("Inversion of the twist map %s" % self._map)
                self._maps[n] = map
                return map
                    
    def gen(self,n=0):
        """
        Return the indeterminate generator of this skew polynomial ring.

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
            raise IndexError, "generator n not defined"
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
        """
        return self.gen()

    def is_finite(self):
        """
        Return False since skew polynomial rings are not finite (unless the
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
        Return True if elements of this skew polynomial ring are exact.
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
        Return true if elements of this polynomial ring have a sparse
        representation.

        Since sparse skew polynomials are not yet implemented, this
        function always returns False.
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
           (default: False)

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
            return self ([ R.random_element (*args, **kwds) for _ in range (degree) ] + [ R(1) ])
        else:
            return self ([ R.random_element (*args, **kwds) for _ in range (degree+1) ])

    def is_commutative(self):
        """
        Return true if this skew polynomial ring is commutative
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
