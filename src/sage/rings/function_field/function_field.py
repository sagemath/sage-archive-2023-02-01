r"""
Function Fields

A function field is an extension field of transcendental degree one over a
rational function field over arbitrary constant field. A function field in Sage
is an extension either of a rational function field or of another function
field.

EXAMPLES:

We create an extension of a rational function field, and do some
simple arithmetic in it. ::

    sage: K.<x> = FunctionField(GF(5^2,'a')); K
    Rational function field in x over Finite Field in a of size 5^2
    sage: R.<y> = K[]
    sage: L.<y> = K.extension(y^3 - (x^3 + 2*x*y + 1/x)); L
    Function field in y defined by y^3 + 3*x*y + (4*x^4 + 4)/x
    sage: y^2
    y^2
    sage: y^3
    2*x*y + (x^4 + 1)/x
    sage: a = 1/y; a
    (4*x/(4*x^4 + 4))*y^2 + 2*x^2/(4*x^4 + 4)
    sage: a * y
    1

We next make an extension of the above function field, illustrating
that arithmetic with a tower of 3 fields is fully supported. ::

    sage: S.<t> = L[]
    sage: M.<t> = L.extension(t^2 - x*y)
    sage: M
    Function field in t defined by t^2 + 4*x*y
    sage: t^2
    x*y
    sage: 1/t
    ((1/(x^4 + 1))*y^2 + 2*x/(4*x^4 + 4))*t
    sage: M.base_field()
    Function field in y defined by y^3 + 3*x*y + (4*x^4 + 4)/x
    sage: M.base_field().base_field()
    Rational function field in x over Finite Field in a of size 5^2

It is also possible to construct function fields over an imperfect base field::

    sage: N.<u> = FunctionField(K)

and function fields as inseparable extensions::

    sage: R.<v> = K[]
    sage: O.<v> = K.extension(v^5 - x)

However, most of advanced computations are available only for global function
fields as yet. A global function field is an extension field of a rational
function field over a finite field by a separable and irreducible polynomial
over the rational function field.

EXAMPLES:

We compute all the Weierstrass points of the Klein quartic over `\GF{2}` and
gap numbers for ordinary places. ::

    sage: K.<x> = FunctionField(GF(2)); _.<Y> = K[]
    sage: L.<y> = K.extension(Y^3 + x^3*Y + x)
    sage: L.genus()
    3
    sage: L.weierstrass_places()
    [Place (1/x, 1/x^3*y^2 + 1/x),
     Place (1/x, 1/x^3*y^2 + 1/x^2*y + 1),
     Place (x, y),
     Place (x + 1, (x^3 + 1)*y + x + 1),
     Place (x^3 + x + 1, y + 1),
     Place (x^3 + x + 1, y + x^2),
     Place (x^3 + x + 1, y + x^2 + 1),
     Place (x^3 + x^2 + 1, y + x),
     Place (x^3 + x^2 + 1, y + x^2 + 1),
     Place (x^3 + x^2 + 1, y + x^2 + x + 1)]
    sage: L.gaps()
    [1, 2, 3]

TESTS::

    sage: TestSuite(K).run()
    sage: TestSuite(L).run()  # long time (8s on sage.math, 2012)
    sage: TestSuite(M).run()  # long time (52s on sage.math, 2012)
    sage: TestSuite(N).run()  # long time
    sage: TestSuite(O).run()  # long time

The following two test suites do not pass ``_test_elements`` yet since
``R.an_element()`` has a ``_test_category`` method which it should not have.
It is not the fault of the function field code so this will
be fixed in another ticket::

    sage: TestSuite(R).run(skip = '_test_elements')
    sage: TestSuite(S).run(skip = '_test_elements')

AUTHORS:

- William Stein (2010): initial version

- Robert Bradshaw (2010-05-30): added is_finite()

- Julian Rueth (2011-06-08, 2011-09-14, 2014-06-23, 2014-06-24, 2016-11-13):
  fixed hom(), extension(); use @cached_method; added derivation(); added
  support for relative vector spaces; fixed conversion to base fields

- Maarten Derickx (2011-09-11): added doctests

- Syed Ahmad Lavasani (2011-12-16): added genus(), is_RationalFunctionField()

- Simon King (2014-10-29): use the same generator names for a function field
  extension and the underlying polynomial ring

- Kwankyu Lee (2016): added global function fields

"""
from __future__ import absolute_import
#*****************************************************************************
#       Copyright (C) 2010 William Stein <wstein@gmail.com>
#       Copyright (C) 2010 Robert Bradshaw <robertwb@math.washington.edu>
#       Copyright (C) 2011-2016 Julian Rueth <julian.rueth@gmail.com>
#       Copyright (C) 2011 Maarten Derickx <m.derickx.student@gmail.com>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#  as published by the Free Software Foundation; either version 2 of
#  the License, or (at your option) any later version.
#                  http://www.gnu.org/licenses/
#*****************************************************************************
from sage.misc.cachefunc import cached_method
from sage.misc.lazy_import import lazy_import

from sage.interfaces.all import singular

from sage.arith.all import lcm

from sage.rings.ring import Field
from sage.rings.polynomial.polynomial_ring_constructor import PolynomialRing
from sage.modules.free_module_element import vector

from sage.categories.homset import Hom
from sage.categories.morphism import SetMorphism
from sage.categories.function_fields import FunctionFields
CAT = FunctionFields()

from .element import (
    FunctionFieldElement,
    FunctionFieldElement_rational,
    FunctionFieldElement_polymod,
    FunctionFieldElement_global)

from .order import (
    FunctionFieldOrder_basis,
    FunctionFieldOrderInfinite_basis,
    FunctionFieldMaximalOrder_rational,
    FunctionFieldMaximalOrder_global,
    FunctionFieldMaximalOrderInfinite_rational,
    FunctionFieldMaximalOrderInfinite_global)

lazy_import('sage.rings.function_field.differential', 'DifferentialsSpace')
lazy_import('sage.matrix.constructor', 'matrix')
lazy_import('sage.modules.free_module', 'FreeModule')

def is_FunctionField(x):
    """
    Return True if x is a function field.

    EXAMPLES::

        sage: from sage.rings.function_field.function_field import is_FunctionField
        sage: is_FunctionField(QQ)
        False
        sage: is_FunctionField(FunctionField(QQ, 't'))
        True
    """
    if isinstance(x, FunctionField): return True
    return x in FunctionFields()

def is_RationalFunctionField(x):
    """
    Return True if x is a rational function field.

    EXAMPLES::

        sage: from sage.rings.function_field.function_field import is_RationalFunctionField
        sage: is_RationalFunctionField(QQ)
        False
        sage: is_RationalFunctionField(FunctionField(QQ, 't'))
        True
    """
    if isinstance(x, RationalFunctionField):
        return True
    else:
        return False

class FunctionField(Field):
    """
    The abstract base class for all function fields.

    - ``base_field`` -- a function fied, of which this function field is an
      extension

    - ``names`` -- a string that gives the name of the generator

    EXAMPLES::

        sage: K.<x> = FunctionField(QQ)
        sage: K
        Rational function field in x over Rational Field
    """
    def __init__(self, base_field, names, category=CAT):
        """
        Initialize.

        TESTS::

            sage: K.<x> = FunctionField(QQ)
            sage: from sage.rings.function_field.function_field import FunctionField
            sage: isinstance(K, FunctionField)
            True
        """
        Field.__init__(self, base_field, names=names, category=category)

        # allow conversion into the constant base field
        from sage.categories.homset import Hom
        from .maps import FunctionFieldConversionToConstantBaseField
        to_constant_base_field = FunctionFieldConversionToConstantBaseField(Hom(self, self.constant_base_field()))
        # the conversion map must not keep the field alive if that is the only reference to it
        to_constant_base_field._make_weak_references()
        self.constant_base_field().register_conversion(to_constant_base_field)

    def is_perfect(self):
        r"""
        Return whether the field is perfect, i.e., its characteristic is `p=0`
        or every element has a `p`-th root.

        EXAMPLES::

            sage: FunctionField(QQ, 'x').is_perfect()
            True
            sage: FunctionField(GF(2), 'x').is_perfect()
            False
        """
        return self.characteristic() == 0

    def some_elements(self):
         """
         Return a list of elements in the function field.

         EXAMPLES::

            sage: K.<x> = FunctionField(QQ)
            sage: elements = K.some_elements()
            sage: elements # random output
            [(x - 3/2)/(x^2 - 12/5*x + 1/18)]
            sage: False in [e in K for e in elements]
            False
         """
         return [self.random_element(), self.random_element(), self.random_element()]

    def characteristic(self):
        """
        Return the characteristic of the function field.

        EXAMPLES::

            sage: K.<x> = FunctionField(QQ)
            sage: K.characteristic()
            0
            sage: K.<x> = FunctionField(GF(7))
            sage: K.characteristic()
            7
            sage: R.<y> = K[]
            sage: L.<y> = K.extension(y^2 - x)
            sage: L.characteristic()
            7
        """
        return self.constant_base_field().characteristic()

    def is_finite(self):
        """
        Return whether the function field is finite, which it is not.

        EXAMPLES::

            sage: R.<t> = FunctionField(QQ)
            sage: R.is_finite()
            False
            sage: R.<t> = FunctionField(GF(7))
            sage: R.is_finite()
            False
        """
        return False

    def is_global(self):
        """
        Return whether the function field is global, that is, whether
        the constant base_field is finite.

        EXAMPLES::

            sage: R.<t> = FunctionField(QQ)
            sage: R.is_global()
            False
            sage: R.<t> = FunctionField(GF(7))
            sage: R.is_global()
            True
        """
        return self.constant_base_field().is_finite()

    def extension(self, f, names=None):
        """
        Create an extension `K(y)` of this function field `K` extended with
        a root `y` of the univariate polynomial `f` over `K`.

        INPUT:

        - ``f`` -- a univariate polynomial over `K`

        - ``names`` -- a string or tuple of length 1 that names the variable `y`

        OUTPUT:

        - a function field

        EXAMPLES::

            sage: K.<x> = FunctionField(QQ); R.<y> = K[]
            sage: K.extension(y^5 - x^3 - 3*x + x*y)
            Function field in y defined by y^5 + x*y - x^3 - 3*x

        A nonintegral defining polynomial::

            sage: K.<t> = FunctionField(QQ); R.<y> = K[]
            sage: K.extension(y^3 + (1/t)*y + t^3/(t+1), 'z')
            Function field in z defined by z^3 + 1/t*z + t^3/(t + 1)

        The defining polynomial need not be monic or integral::

            sage: K.extension(t*y^3 + (1/t)*y + t^3/(t+1))
            Function field in y defined by t*y^3 + 1/t*y + t^3/(t + 1)
        """
        from . import constructor
        return constructor.FunctionFieldExtension(f, names)

    def order_with_basis(self, basis, check=True):
        """
        Return the order with given basis over the maximal order of
        the base field.

        INPUT:

        - ``basis`` -- a list of elements of self

        - ``check`` -- bool (default: True); if True, check that the basis is
          really linearly independent and that the module it spans is closed
          under multiplication, and contains the identity element.

        OUTPUT:

        - an order in the function field

        EXAMPLES::

            sage: K.<x> = FunctionField(QQ); R.<y> = K[]; L.<y> = K.extension(y^3 + x^3 + 4*x + 1)
            sage: O = L.order_with_basis([1, y, y^2]); O
            Order in Function field in y defined by y^3 + x^3 + 4*x + 1
            sage: O.basis()
            (1, y, y^2)

        Note that 1 does not need to be an element of the basis, as long it is in the module spanned by it::

            sage: O = L.order_with_basis([1+y, y, y^2]); O
            Order in Function field in y defined by y^3 + x^3 + 4*x + 1
            sage: O.basis()
            (y + 1, y, y^2)

        The following error is raised when the module spanned by the basis is not closed under multiplication::

            sage: O = L.order_with_basis([1, x^2 + x*y, (2/3)*y^2]); O
            Traceback (most recent call last):
            ...
            ValueError: The module generated by basis [1, x*y + x^2, 2/3*y^2] must be closed under multiplication

        and this happens when the identity is not in the module spanned by the basis::

            sage: O = L.order_with_basis([x, x^2 + x*y, (2/3)*y^2])
            Traceback (most recent call last):
            ...
            ValueError: The identity element must be in the module spanned by basis [x, x*y + x^2, 2/3*y^2]
        """
        from .order import FunctionFieldOrder_basis
        return FunctionFieldOrder_basis([self(a) for a in basis], check=check)

    def order(self, x, check=True):
        """
        Return the order in the function field generated over the maximal order
        by x or the elements of x if x is a list.

        INPUT:

        - ``x`` -- element of self, or a list of elements of self

        - ``check`` -- bool (default: True); if True, check that x really
          generates an order

        EXAMPLES::

            sage: K.<x> = FunctionField(QQ); R.<y> = K[]
            sage: L.<y> = K.extension(y^3 + x^3 + 4*x + 1)
            sage: O = L.order(y); O
            Order in Function field in y defined by y^3 + x^3 + 4*x + 1
            sage: O.basis()
            (1, y, y^2)

            sage: Z = K.order(x); Z
            Order in Rational function field in x over Rational Field
            sage: Z.basis()
            (1,)

        Orders with multiple generators, not yet supported::

            sage: Z = K.order([x,x^2]); Z
            Traceback (most recent call last):
            ...
            NotImplementedError
        """
        if not isinstance(x, (list, tuple)):
            x = [x]
        if len(x) == 1:
            g = x[0]
            basis = [self(1)]
            for i in range(self.degree()-1):
                basis.append(basis[-1]*g)
        else:
            raise NotImplementedError
        return self.order_with_basis(basis, check=check)

    def order_infinite_with_basis(self, basis, check=True):
        """
        Return the order with given basis over the maximal infinite order of
        the base field.

        INPUT:

        - ``basis`` -- a list of elements of self

        - ``check`` -- bool (default: True); if True, check that the basis
          is really linearly independent and that the module it spans is
          closed under multiplication, and contains the identity element.

        OUTPUT:

        - an order in the function field

        EXAMPLES::

            sage: K.<x> = FunctionField(QQ); R.<y> = K[]
            sage: L.<y> = K.extension(y^3 + x^3 + 4*x + 1)
            sage: O = L.order_with_basis([1, y, y^2]); O
            Order in Function field in y defined by y^3 + x^3 + 4*x + 1
            sage: O.basis()
            (1, y, y^2)

        Note that 1 does not need to be an element of the basis, as long it is
        in the module spanned by it::

            sage: O = L.order_with_basis([1+y, y, y^2]); O
            Order in Function field in y defined by y^3 + x^3 + 4*x + 1
            sage: O.basis()
            (y + 1, y, y^2)

        The following error is raised when the module spanned by the basis is
        not closed under multiplication::

            sage: O = L.order_with_basis([1, x^2 + x*y, (2/3)*y^2]); O
            Traceback (most recent call last):
            ...
            ValueError: The module generated by basis [1, x*y + x^2, 2/3*y^2]
            must be closed under multiplication

        and this happens when the identity is not in the module spanned by the
        basis::

            sage: O = L.order_with_basis([x, x^2 + x*y, (2/3)*y^2])
            Traceback (most recent call last):
            ...
            ValueError: The identity element must be in the module spanned by
            basis [x, x*y + x^2, 2/3*y^2]
        """
        if len(basis) == 1:
            basis = basis[0]
            if not isinstance(basis, (list,tuple)):
                basis = (basis,)
        basis = [self(g) for g in basis]
        from .order import FunctionFieldOrderInfinite_basis
        return FunctionFieldOrderInfinite_basis(basis, check=check)

    def order_infinite(self, x, check=True):
        """
        Return the order in the function field generated over the
        maximal infinite order by x or the elements of x if x is a list.

        INPUT:

        - ``x`` -- an element or a list of elements of the function field

        - ``check`` -- bool (default: True); if True, check that
          x really generates an order

        EXAMPLES::

            sage: K.<x> = FunctionField(QQ); R.<y> = K[]; L.<y> = K.extension(y^3 + x^3 + 4*x + 1)
            sage: O = L.order(y); O
            Order in Function field in y defined by y^3 + x^3 + 4*x + 1
            sage: O.basis()
            (1, y, y^2)

            sage: Z = K.order(x); Z
            Order in Rational function field in x over Rational Field
            sage: Z.basis()
            (1,)

        Orders with multiple generators, not yet supported::

            sage: Z = K.order([x,x^2]); Z
            Traceback (most recent call last):
            ...
            NotImplementedError
        """
        if not isinstance(x, (list, tuple)):
            x = [x]
        if len(x) == 1:
            g = x[0]
            basis = [self(1)]
            for i in range(self.degree()-1):
                basis.append(basis[-1]*g)
        else:
            raise NotImplementedError
        return self.order_infinite_with_basis(basis, check=check)

    def _coerce_map_from_(self, R):
        """
        Return True if there is a coerce map from R to the function field.

        INPUT:

        - ``R`` -- a ring

        EXAMPLES::

            sage: K.<x> = FunctionField(QQ); R.<y> = K[]; L.<y> = K.extension(y^3 + x^3 + 4*x + 1)
            sage: L.equation_order()
            Order in Function field in y defined by y^3 + x^3 + 4*x + 1
            sage: L._coerce_map_from_(L.equation_order())
            True
            sage: L._coerce_map_from_(GF(7))
            False
        """
        from .order import FunctionFieldOrder
        if isinstance(R, FunctionFieldOrder) and R.fraction_field() == self:
            return True
        return False

    def _test_derivation(self, **options):
        """
        Test the correctness of the derivations of the function field.

        EXAMPLES::

            sage: K.<x> = FunctionField(QQ)
            sage: TestSuite(K).run() # indirect doctest
        """
        tester = self._tester(**options)
        S = tester.some_elements()
        K = self.constant_base_field().some_elements()
        try:
            d = self.derivation()
        except NotImplementedError:
            return # some function fields no not implement derivation() yet
        from itertools import product
        # Leibniz's law
        for x,y in tester.some_elements(product(S, S)):
            tester.assert_(d(x*y) == x*d(y) + d(x)*y)
        # Linearity
        for x,y in tester.some_elements(product(S, S)):
            tester.assert_(d(x+y) == d(x) + d(y))
        for c,x in tester.some_elements(product(K, S)):
            tester.assert_(d(c*x) == c*d(x))
        # Constants map to zero
        for c in tester.some_elements(K):
            tester.assert_(d(c) == 0)

    def _convert_map_from_(self, R):
        """
        Return a conversion from R to the function field or None if exists
        none.

        EXAMPLES::

            sage: K.<x> = FunctionField(QQ)
            sage: R.<y> = K[]
            sage: L.<y> = K.extension(y^3 + x^3 + 4*x + 1)
            sage: K(L(x)) # indirect doctest
            x
        """
        if isinstance(R, FunctionField_polymod):
            base_conversion = self.convert_map_from(R.base_field())
            if base_conversion is not None:
                from sage.categories.morphism import SetMorphism
                return base_conversion * SetMorphism(R.Hom(R.base_field()), R._to_base_field)

    def _intermediate_fields(self, base):
        r"""
        Return the fields which lie in between base and this field in the
        tower of function fields.

        INPUT:

        - ``base`` -- a function field, either this field or a field from which
          this field has been created as an extension

        OUTPUT:

        - a list of fields; The first entry is ``base``, the last entry is this field.

        EXAMPLES::

            sage: K.<x> = FunctionField(QQ)
            sage: K._intermediate_fields(K)
            [Rational function field in x over Rational Field]

            sage: R.<y> = K[]
            sage: L.<y> = K.extension(y^2-x)
            sage: L._intermediate_fields(K)
            [Function field in y defined by y^2 - x, Rational function field in x over Rational Field]

            sage: R.<z> = L[]
            sage: M.<z> = L.extension(z^2-y)
            sage: M._intermediate_fields(L)
            [Function field in z defined by z^2 - y, Function field in y defined by y^2 - x]
            sage: M._intermediate_fields(K)
            [Function field in z defined by z^2 - y, Function field in y defined by y^2 - x,
            Rational function field in x over Rational Field]

        TESTS::

            sage: K._intermediate_fields(M)
            Traceback (most recent call last):
            ...
            ValueError: field has not been constructed as a finite extension of base
            sage: K._intermediate_fields(QQ)
            Traceback (most recent call last):
            ...
            TypeError: base must be a function field
        """
        if not is_FunctionField(base):
            raise TypeError("base must be a function field")

        ret = [self]
        while ret[-1] is not base:
            ret.append(ret[-1].base_field())
            if ret[-1] is ret[-2]:
                raise ValueError("field has not been constructed as a finite extension of base")
        return ret

class FunctionField_polymod(FunctionField):
    """
    Function fields defined by a univariate polynomial, as an extension of the
    base field.

    EXAMPLES:

    We make a function field defined by a degree 5 polynomial over the
    rational function field over the rational numbers::

        sage: K.<x> = FunctionField(QQ)
        sage: R.<y> = K[]
        sage: L.<y> = K.extension(y^5 - (x^3 + 2*x*y + 1/x)); L
        Function field in y defined by y^5 - 2*x*y + (-x^4 - 1)/x

    We next make a function field over the above nontrivial function
    field L::

        sage: S.<z> = L[]
        sage: M.<z> = L.extension(z^2 + y*z + y); M
        Function field in z defined by z^2 + y*z + y
        sage: 1/z
        ((x/(-x^4 - 1))*y^4 - 2*x^2/(-x^4 - 1))*z - 1
        sage: z * (1/z)
        1

    We drill down the tower of function fields::

        sage: M.base_field()
        Function field in y defined by y^5 - 2*x*y + (-x^4 - 1)/x
        sage: M.base_field().base_field()
        Rational function field in x over Rational Field
        sage: M.base_field().base_field().constant_field()
        Rational Field
        sage: M.constant_base_field()
        Rational Field

    .. WARNING::

        It is not checked if the polynomial used to define the function field is irreducible
        Hence it is not guaranteed that this object really is a field!
        This is illustrated below.

    ::

        sage: K.<x>=FunctionField(QQ)
        sage: R.<y> = K[]
        sage: L.<y>=K.extension(x^2-y^2)
        sage: (y-x)*(y+x)
        0
        sage: 1/(y-x)
        1
        sage: y-x==0; y+x==0
        False
        False
    """
    def __init__(self, polynomial, names, element_class=FunctionFieldElement_polymod, category=CAT):
        """
        Create a function field defined as an extension of another
        function field by adjoining a root of a univariate polynomial.

        INPUT:

        - ``polynomial`` -- a univariate polynomial over a function field

        - ``names`` -- a tuple of length 1 or a string; variable names

        - ``category`` -- a category (defaults to category of function fields)

        EXAMPLES:

        We create an extension of a function field::

            sage: K.<x> = FunctionField(QQ); R.<y> = K[]
            sage: L = K.extension(y^5 - x^3 - 3*x + x*y); L
            Function field in y defined by y^5 + x*y - x^3 - 3*x

        Note the type::

            sage: type(L)
            <class 'sage.rings.function_field.function_field.FunctionField_polymod_with_category'>

        We can set the variable name, which doesn't have to be y::

            sage: L.<w> = K.extension(y^5 - x^3 - 3*x + x*y); L
            Function field in w defined by w^5 + x*w - x^3 - 3*x

        TESTS:

        Test that :trac:`17033` is fixed::

            sage: K.<t> = FunctionField(QQ)
            sage: R.<x> = QQ[]
            sage: M.<z> = K.extension(x^7-x-t)
            sage: M(x)
            z
            sage: M('z')
            z
            sage: M('x')
            Traceback (most recent call last):
            ...
            TypeError: unable to evaluate 'x' in Fraction Field of Univariate
            Polynomial Ring in t over Rational Field
        """
        from sage.rings.polynomial.polynomial_element import is_Polynomial
        if polynomial.parent().ngens()>1 or not is_Polynomial(polynomial):
            raise TypeError("polynomial must be univariate a polynomial")
        if names is None:
            names = (polynomial.variable_name(), )
        elif names != polynomial.variable_name():
            polynomial = polynomial.change_variable_name(names)
        if polynomial.degree() <= 0:
            raise ValueError("polynomial must have positive degree")
        base_field = polynomial.base_ring()
        if not isinstance(base_field, FunctionField):
            raise TypeError("polynomial must be over a FunctionField")
        self._element_class = element_class
        self._element_init_pass_parent = False
        self._base_field = base_field
        self._polynomial = polynomial

        FunctionField.__init__(self, base_field, names=names, category = category)

        self._hash = hash(polynomial)
        self._ring = self._polynomial.parent()
        self._populate_coercion_lists_(coerce_list=[base_field, self._ring])
        self._gen = self(self._ring.gen())

    def __reduce__(self):
        """
        Return the arguments which were used to create this instance.

        The rationale for this is explained in the documentation of
        ``UniqueRepresentation``.

        EXAMPLES::

            sage: K.<x> = FunctionField(QQ); R.<y> = K[]
            sage: L = K.extension(y^2 - x)
            sage: clazz,args = L.__reduce__()
            sage: clazz(*args)
            Function field in y defined by y^2 - x
        """
        from .constructor import FunctionFieldExtension
        return  FunctionFieldExtension, (self._polynomial, self._names)

    def __hash__(self):
        """
        Return hash of the function field.

        The hash value is equal to the hash of the defining polynomial.

        EXAMPLES::

            sage: K.<x> = FunctionField(QQ); R.<y> = K[]
            sage: L = K.extension(y^5 - x^3 - 3*x + x*y)
            sage: hash(L) == hash(L.polynomial())
            True
        """
        return self._hash

    def _element_constructor_(self, x):
        r"""
        Make x into an element of the function field, possibly not canonically.

        INPUT:

        - ``x`` -- the element

        OUTPUT:

        - ``x``, as an element of the function field

        TESTS::

            sage: K.<x> = FunctionField(QQ); R.<y> = K[]
            sage: L.<y> = K.extension(y^5 - (x^3 + 2*x*y + 1/x))
            sage: L._element_constructor_(L.polynomial_ring().gen())
            y
        """
        if isinstance(x, FunctionFieldElement):
            return self._element_class(self, self._ring(x.element()))
        return self._element_class(self, self._ring(x))

    def gen(self, n=0):
        """
        Return the n-th generator of the function field. By default, n is 0; any other
        value of n leads to an error. The generator is the class of `y`, if we view
        self as being presented as `K[y]/(f(y))`.

        EXAMPLES::

            sage: K.<x> = FunctionField(QQ); R.<y> = K[]
            sage: L.<y> = K.extension(y^5 - (x^3 + 2*x*y + 1/x))
            sage: L.gen()
            y
            sage: L.gen(1)
            Traceback (most recent call last):
            ...
            IndexError: Only one generator.
        """
        if n != 0: raise IndexError("Only one generator.")
        return self._gen

    def ngens(self):
        """
        Return the number of generators of the function field over its base
        field. This is by definition 1.

        EXAMPLES::

            sage: K.<x> = FunctionField(QQ); R.<y> = K[]
            sage: L.<y> = K.extension(y^5 - (x^3 + 2*x*y + 1/x))
            sage: L.ngens()
            1
        """
        return 1

    def _to_base_field(self, f):
        r"""
        Return f as an element of the :meth:`base_field`.

        INPUT:

        - ``f`` -- an element of the function field which lies in the base
          field.

        EXAMPLES::

            sage: K.<x> = FunctionField(QQ)
            sage: R.<y> = K[]
            sage: L.<y> = K.extension(y^2 - x)
            sage: L._to_base_field(L(x))
            x
            sage: L._to_base_field(y)
            Traceback (most recent call last):
            ...
            ValueError: y is not an element of the base field

        TESTS:

        Verify that :trac:`21872` has been resolved::

            sage: R.<z> = L[]
            sage: M.<z> = L.extension(z^2 - y)

            sage: M(1) in QQ
            True
            sage: M(y) in L
            True
            sage: M(x) in K
            True
            sage: z in K
            False
        """
        K = self.base_field()
        if f.element().is_constant():
            return K(f.element())
        raise ValueError("%r is not an element of the base field"%(f,))

    def _to_constant_base_field(self, f):
        """
        Return f as an element of the :meth:`constant_base_field`.

        INPUT:

        - ``f`` -- an element of the rational function field which is a
          constant

        EXAMPLES::

            sage: K.<x> = FunctionField(QQ)
            sage: R.<y> = K[]
            sage: L.<y> = K.extension(y^2 - x)
            sage: L._to_constant_base_field(L(1))
            1
            sage: L._to_constant_base_field(y)
            Traceback (most recent call last):
            ...
            ValueError: y is not an element of the base field

        TESTS:

        Verify that :trac:`21872` has been resolved::

            sage: L(1) in QQ
            True
            sage: y in QQ
            False
        """
        return self.base_field()._to_constant_base_field(self._to_base_field(f))

    def monic_integral_model(self, names):
        """
        Return a function field isomorphic to the function field, but with
        defining polynomial that is monic and integral over the base field.

        INPUT:

        - ``names`` -- name of the generator of the new function field

        EXAMPLES::

            sage: K.<x> = FunctionField(QQ); R.<y> = K[]
            sage: L.<y> = K.extension(x^2*y^5 - 1/x); L
            Function field in y defined by x^2*y^5 - 1/x
            sage: A, from_A, to_A = L.monic_integral_model('z')
            sage: A
            Function field in z defined by z^5 - x^12
            sage: from_A
            Function Field morphism:
              From: Function field in z defined by z^5 - x^12
              To:   Function field in y defined by x^2*y^5 - 1/x
              Defn: z |--> x^3*y
            sage: to_A
            Function Field morphism:
              From: Function field in y defined by x^2*y^5 - 1/x
              To:   Function field in z defined by z^5 - x^12
              Defn: y |--> 1/x^3*z
            sage: to_A(y)
            1/x^3*z
            sage: from_A(to_A(y))
            y
            sage: from_A(to_A(1/y))
            x^3*y^4
            sage: from_A(to_A(1/y)) == 1/y
            True
        """
        g, d = self._make_monic_integral(self.polynomial())
        R = self.base_field()
        K = R.extension(g, names=names)
        to_K = self.hom(K.gen() / d)
        from_K = K.hom(self.gen() * d)
        return K, from_K, to_K

    def _make_monic_integral(self, f):
        r"""
        Let y be a root of f. This function returns a monic integral polynomial
        g and an element d of the base field such that g(y*d)=0.

        EXAMPLES::

            sage: K.<x> = FunctionField(QQ); R.<y> = K[];
            sage: L.<y> = K.extension(x^2*y^5 - 1/x)
            sage: g, d = L._make_monic_integral(L.polynomial()); g,d
            (y^5 - x^12, x^3)
            sage: (y*d).is_integral()
            True
            sage: g.is_monic()
            True
            sage: g(y*d)
            0
        """
        R = f.base_ring()
        if not isinstance(R, RationalFunctionField):
            raise NotImplementedError

        # make f monic
        n = f.degree()
        c = f.leading_coefficient()
        if c != 1:
            f = f / c

        # find lcm of denominators
        from sage.arith.all import lcm
        # would be good to replace this by minimal...
        d = lcm([b.denominator() for b in f.list() if b])
        if d != 1:
            x = f.parent().gen()
            g = (d**n) * f(x/d)
        else:
            g = f
        return g, d

    def constant_field(self):
        """
        Return the algebraic closure of the constant field of the base
        field in the function field.

        EXAMPLES::

            sage: K.<x> = FunctionField(QQ); R.<y> = K[]
            sage: L.<y> = K.extension(y^5 - (x^3 + 2*x*y + 1/x))
            sage: L.constant_field()
            Traceback (most recent call last):
            ...
            NotImplementedError
        """
        raise NotImplementedError

    def constant_base_field(self):
        """
        Return the constant field of the base rational function field.

        EXAMPLES::

            sage: K.<x> = FunctionField(QQ); R.<y> = K[]
            sage: L.<y> = K.extension(y^5 - (x^3 + 2*x*y + 1/x)); L
            Function field in y defined by y^5 - 2*x*y + (-x^4 - 1)/x
            sage: L.constant_base_field()
            Rational Field
            sage: S.<z> = L[]
            sage: M.<z> = L.extension(z^2 - y)
            sage: M.constant_base_field()
            Rational Field
        """
        return self.base_field().constant_base_field()

    @cached_method(key=lambda self, base: self.base_field() if base is None else base)
    def degree(self, base=None):
        """
        Return the degree of the function field over the function field base.

        INPUT:

        - ``base`` -- a function field (default: ``None``), a
          function field from which this field has been constructed as a finite
          extension.

        EXAMPLES::

            sage: K.<x> = FunctionField(QQ)
            sage: R.<y> = K[]
            sage: L.<y> = K.extension(y^5 - (x^3 + 2*x*y + 1/x)); L
            Function field in y defined by y^5 - 2*x*y + (-x^4 - 1)/x
            sage: L.degree()
            5
            sage: L.degree(L)
            1

            sage: R.<z> = L[]
            sage: M.<z> = L.extension(z^2 - y)
            sage: M.degree(L)
            2
            sage: M.degree(K)
            10

        TESTS::

            sage: L.degree(M)
            Traceback (most recent call last):
            ...
            ValueError: base must be None or the rational function field

        """
        if base is None:
            base = self.base_field()
        if base is self:
            from sage.rings.integer_ring import ZZ
            return ZZ(1)
        return self._polynomial.degree() * self.base_field().degree(base)

    def _repr_(self):
        """
        Return the string representation of the function field.

        EXAMPLES::

            sage: K.<x> = FunctionField(QQ); R.<y> = K[]
            sage: L.<y> = K.extension(y^5 - (x^3 + 2*x*y + 1/x))
            sage: L._repr_()
            'Function field in y defined by y^5 - 2*x*y + (-x^4 - 1)/x'
        """
        return "Function field in %s defined by %s"%(self.variable_name(), self._polynomial)

    def base_field(self):
        """
        Return the base field of the function field.  This function
        field is presented as `L = K[y]/(f(y))`, and the base field is
        by definition the field `K`.

        EXAMPLES::

            sage: K.<x> = FunctionField(QQ); R.<y> = K[]
            sage: L.<y> = K.extension(y^5 - (x^3 + 2*x*y + 1/x))
            sage: L.base_field()
            Rational function field in x over Rational Field
        """
        return self._base_field

    def random_element(self, *args, **kwds):
        """
        Create a random element of the function field.  Parameters
        are passed onto the random_element method of the base_field.

        EXAMPLES::

            sage: K.<x> = FunctionField(QQ); R.<y> = K[]
            sage: L.<y> = K.extension(y^2 - (x^2 + x))
            sage: L.random_element() # random
            ((x^2 - x + 2/3)/(x^2 + 1/3*x - 1))*y^2 + ((-1/4*x^2 + 1/2*x - 1)/(-5/2*x + 2/3))*y
            + (-1/2*x^2 - 4)/(-12*x^2 + 1/2*x - 1/95)
        """
        return self(self._ring.random_element(degree=self.degree(), *args, **kwds))

    def polynomial(self):
        """
        Return the univariate polynomial that defines the function
        field, i.e., the polynomial f(y) so that the function field
        is of the form K[y]/(f(y)).

        EXAMPLES::

            sage: K.<x> = FunctionField(QQ); R.<y> = K[]
            sage: L.<y> = K.extension(y^5 - (x^3 + 2*x*y + 1/x))
            sage: L.polynomial()
            y^5 - 2*x*y + (-x^4 - 1)/x
        """
        return self._polynomial

    def polynomial_ring(self):
        """
        Return the polynomial ring used to represent elements of the
        function field.  If we view the function field as being presented
        as K[y]/(f(y)), then this function returns the ring K[y].

        EXAMPLES::

            sage: K.<x> = FunctionField(QQ); R.<y> = K[]
            sage: L.<y> = K.extension(y^5 - (x^3 + 2*x*y + 1/x))
            sage: L.polynomial_ring()
            Univariate Polynomial Ring in y over Rational function field in x over Rational Field
        """
        return self._ring

    @cached_method(key=lambda self, base: self.base_field() if base is None else base)
    def vector_space(self, base=None):
        """
        Return a vector space `V` and isomorphisms from the field to `V` and
        from `V` to the field.

        This function allows us to identify the elements of this field with
        elements of a vector space over the base field, which is useful for
        representation and arithmetic with orders, ideals, etc.

        INPUT:

        - ``base`` -- the base function field. The returned vector space is
          over the base function field. If not given, it defaults to the base
          field of the function field.

        OUTPUT:

        - ``V`` -- a vector space over the base function field

        - ``from_V`` -- an isomorphism from `V` to the field

        - ``to_V`` -- an isomorphism from the field to `V`

        EXAMPLES:

        We define a function field::

            sage: K.<x> = FunctionField(QQ); R.<y> = K[]
            sage: L.<y> = K.extension(y^5 - (x^3 + 2*x*y + 1/x)); L
            Function field in y defined by y^5 - 2*x*y + (-x^4 - 1)/x

        We get the vector spaces, and maps back and forth::

            sage: V, from_V, to_V = L.vector_space()
            sage: V
            Vector space of dimension 5 over Rational function field in x over Rational Field
            sage: from_V
            Isomorphism morphism:
              From: Vector space of dimension 5 over Rational function field in x over Rational Field
              To:   Function field in y defined by y^5 - 2*x*y + (-x^4 - 1)/x
            sage: to_V
            Isomorphism morphism:
              From: Function field in y defined by y^5 - 2*x*y + (-x^4 - 1)/x
              To:   Vector space of dimension 5 over Rational function field in x over Rational Field

        We convert an element of the vector space back to the function field::

            sage: from_V(V.1)
            y

        We define an interesting element of the function field::

            sage: a = 1/L.0; a
            (-x/(-x^4 - 1))*y^4 + 2*x^2/(-x^4 - 1)

        We convert it to the vector space, and get a vector over the base field::

            sage: to_V(a)
            (2*x^2/(-x^4 - 1), 0, 0, 0, -x/(-x^4 - 1))

        We convert to and back, and get the same element::

            sage: from_V(to_V(a)) == a
            True

        In the other direction::

            sage: v = x*V.0 + (1/x)*V.1
            sage: to_V(from_V(v)) == v
            True

        And we show how it works over an extension of an extension field::

            sage: R2.<z> = L[]; M.<z> = L.extension(z^2 -y)
            sage: M.vector_space()
            (Vector space of dimension 2 over Function field in y defined by y^5 - 2*x*y + (-x^4 - 1)/x, Isomorphism morphism:
              From: Vector space of dimension 2 over Function field in y defined by y^5 - 2*x*y + (-x^4 - 1)/x
              To:   Function field in z defined by z^2 - y, Isomorphism morphism:
              From: Function field in z defined by z^2 - y
              To:   Vector space of dimension 2 over Function field in y defined by y^5 - 2*x*y + (-x^4 - 1)/x)

        We can also get the vector space of ``M`` over ``K``::

            sage: M.vector_space(K)
            (Vector space of dimension 10 over Rational function field in x over Rational Field, Isomorphism morphism:
              From: Vector space of dimension 10 over Rational function field in x over Rational Field
              To:   Function field in z defined by z^2 - y, Isomorphism morphism:
              From: Function field in z defined by z^2 - y
              To:   Vector space of dimension 10 over Rational function field in x over Rational Field)

        """
        from .maps import MapVectorSpaceToFunctionField, MapFunctionFieldToVectorSpace
        if base is None:
            base = self.base_field()
        degree = self.degree(base)
        V = base**degree;
        from_V = MapVectorSpaceToFunctionField(V, self)
        to_V   = MapFunctionFieldToVectorSpace(self, V)
        return (V, from_V, to_V)

    def maximal_order(self):
        """
        Return the maximal order of the function field.

        This is to be implemented in a subclass.
        """
        raise NotImplementedError

    def maximal_order_infinite(self):
        """
        Return the maximal infinite order of the function field.

        This is to be implemented in a subclass.
        """
        raise NotImplementedError

    def equation_order(self):
        """
        Return the equation order of the function field.

        If we view self as being presented as K[y]/(f(y)), then this function
        returns the order generated by the class of y.  If f is not monic, then
        :meth:`_make_monic_integral` is called, and instead we get the order
        generated by some integral multiple of a root of f.

        EXAMPLES::

            sage: K.<x> = FunctionField(QQ); R.<y> = K[]
            sage: L.<y> = K.extension(y^5 - (x^3 + 2*x*y + 1/x))
            sage: O = L.equation_order()
            sage: O.basis()
            (1, x*y, x^2*y^2, x^3*y^3, x^4*y^4)

        We try an example, in which the defining polynomial is not
        monic and is not integral::

            sage: K.<x> = FunctionField(QQ); R.<y> = K[]
            sage: L.<y> = K.extension(x^2*y^5 - 1/x); L
            Function field in y defined by x^2*y^5 - 1/x
            sage: O = L.equation_order()
            sage: O.basis()
            (1, x^3*y, x^6*y^2, x^9*y^3, x^12*y^4)
        """
        d = self._make_monic_integral(self.polynomial())[1]
        return self.order(d*self.gen(), check=False)

    def hom(self, im_gens, base_morphism=None):
        """
        Create a homomorphism from self to another function field.

        INPUT:

        - ``im_gens`` -- a list of images of the generators of self
          and of successive base rings.

        - ``base_morphism`` -- (default: None) a homomorphism of the base ring,
          after the im_gens are used.  Thus if im_gens has length 2, then
          base_morphism should be a morphism from self.base_ring().base_ring().

        EXAMPLES:

        We create a rational function field, and a quadratic extension of it::

            sage: K.<x> = FunctionField(QQ); R.<y> = K[]
            sage: L.<y> = K.extension(y^2 - x^3 - 1)

        We make the field automorphism that sends y to -y::

            sage: f = L.hom(-y); f
            Function Field endomorphism of Function field in y defined by y^2 - x^3 - 1
              Defn: y |--> -y

        Evaluation works::

            sage: f(y*x - 1/x)
            -x*y - 1/x

        We try to define an invalid morphism::

            sage: f = L.hom(y+1)
            Traceback (most recent call last):
            ...
            ValueError: invalid morphism

        We make a morphism of the base rational function field::

            sage: phi = K.hom(x+1); phi
            Function Field endomorphism of Rational function field in x over Rational Field
              Defn: x |--> x + 1
            sage: phi(x^3 - 3)
            x^3 + 3*x^2 + 3*x - 2
            sage: (x+1)^3-3
            x^3 + 3*x^2 + 3*x - 2

        We make a morphism by specifying where the generators and the
        base generators go::

            sage: L.hom([-y, x])
            Function Field endomorphism of Function field in y defined by y^2 - x^3 - 1
              Defn: y |--> -y
                    x |--> x

        The usage of the keyword base_morphism is not implemented yet::

            sage: L.hom([-y, x-1], base_morphism=phi)
            Traceback (most recent call last):
            ...
            NotImplementedError: Function field homorphisms with optional argument base_morphism are not implemented yet. Please specify the images of the generators of the base fields manually.

        We make another extension of a rational function field::

            sage: K2.<t> = FunctionField(QQ); R2.<w> = K2[]
            sage: L2.<w> = K2.extension((4*w)^2 - (t+1)^3 - 1)

        We define a morphism, by giving the images of generators::

            sage: f = L.hom([4*w, t+1]); f
            Function Field morphism:
              From: Function field in y defined by y^2 - x^3 - 1
              To:   Function field in w defined by 16*w^2 - t^3 - 3*t^2 - 3*t - 2
              Defn: y |--> 4*w
                    x |--> t + 1

        Evaluation works, as expected::

            sage: f(y+x)
            4*w + t + 1
            sage: f(x*y + x/(x^2+1))
            (4*t + 4)*w + (t + 1)/(t^2 + 2*t + 2)

        We make another extension of a rational function field::

            sage: K3.<yy> = FunctionField(QQ); R3.<xx> = K3[]
            sage: L3.<xx> = K3.extension(yy^2 - xx^3 - 1)

        This is the function field L with the generators exchanged. We define a morphism to L::

            sage: g = L3.hom([x,y]); g
            Function Field morphism:
              From: Function field in xx defined by -xx^3 + yy^2 - 1
              To:   Function field in y defined by y^2 - x^3 - 1
              Defn: xx |--> x
                    yy |--> y

        """
        if base_morphism is not None:
            raise NotImplementedError("Function field homorphisms with optional argument base_morphism are not implemented yet. Please specify the images of the generators of the base fields manually.")

        if not isinstance(im_gens, (list,tuple)):
            im_gens = [im_gens]
        if len(im_gens) == 0:
            raise ValueError("no images specified")

        if len(im_gens) > 1:
            base_morphism = self.base_field().hom(im_gens[1:], base_morphism)

        # the codomain of this morphism is the field containing all the im_gens
        codomain = im_gens[0].parent();
        if base_morphism is not None:
            if base_morphism.codomain().has_coerce_map_from(codomain):
                codomain = base_morphism.codomain();

        from .maps import FunctionFieldMorphism_polymod
        return FunctionFieldMorphism_polymod(self.Hom(codomain), im_gens[0], base_morphism)

    @cached_method
    def genus(self):
        """
        Return the genus of the function field.

        For now, the genus is computed using singular.

        EXAMPLES::

            sage: K.<x> = FunctionField(QQ); R.<y> = K[]
            sage: L.<y> = K.extension(y^3 - (x^3 + 2*x*y + 1/x))
            sage: L.genus()
            3
        """
        # unfortunately singular can not compute the genus with the polynomial_ring()._singular_
        # object because genus method only accepts a ring of transdental degree 2 over a prime field
        # not a ring of transdental degree 1 over a rational function field of one variable

        if is_RationalFunctionField(self._base_field) and self._base_field.constant_field().is_prime_field():

            #Making the auxiliary ring which only has polynomials with integral coefficients.
            tmpAuxRing = PolynomialRing(self._base_field.constant_field(), str(self._base_field.gen())+','+str(self._ring.gen()))
            intMinPoly, d = self._make_monic_integral(self._polynomial)
            curveIdeal = tmpAuxRing.ideal(intMinPoly)

            singular.lib('normal.lib') #loading genus method in singular
            return int(curveIdeal._singular_().genus())

        else:
            raise NotImplementedError("Computation of genus over the rational function field not implemented yet")

    @cached_method
    def derivation(self):
        r"""
        Return a derivation of the function field over the constant base field.

        A derivation on `R` is a map `R\to R` satisfying
        `D(\alpha+\beta)=D(\alpha)+D(\beta)` and `D(\alpha\beta)=\beta
        D(\alpha)+\alpha D(\beta)` for all `\alpha, \beta \in R`. For a
        function field which is a finite extension of `K(x)` with `K` perfect,
        the derivations form a one-dimensional `K`-vector space generated by
        the derivation returned by this method.

        OUTPUT:

        - a derivation of the function field

        EXAMPLES::

            sage: K.<x> = FunctionField(GF(3))
            sage: R.<y> = K[]
            sage: L.<y> = K.extension(y^2 - x)
            sage: d = L.derivation(); d
            Derivation map:
                From: Function field in y defined by y^2 + 2*x
                To:   Function field in y defined by y^2 + 2*x
                Defn: y |--> 2/x*y
            sage: d(x)
            1
            sage: d(x^3)
            0
            sage: d(x*y)
            0
            sage: d(y)
            2/x*y

        Derivations are linear and satisfy Leibniz's law::

            sage: d(x+y) == d(x) + d(y)
            True
            sage: d(x*y) == x*d(y) + y*d(x)
            True

        If the field is a separable extension of the base field, the derivation
        extending a derivation of the base function field is uniquely
        determined. Proposition 11 of [GT1996]_ describes how to compute the
        extension. We apply the formula described there to the generator
        of the space of derivations on the base field.

        The general inseparable case is not implemented yet (see :trac:`16562`,
        :trac:`16564`.)`
        """
        from .maps import FunctionFieldDerivation_separable
        if self.polynomial().gcd(self.polynomial().derivative()).is_one():
            return FunctionFieldDerivation_separable(self, self.base_ring().derivation())
        else:
            raise NotImplementedError("construction of separable models not implemented")

class FunctionField_global(FunctionField_polymod):
    """
    Global function fields.
    """
    def __init__(self, polynomial, names):
        """
        Initialize the function field.

        INPUT:

        - ``polynomial`` -- a monic irreducible and separable polynomial

        - ``names`` -- name of the generator of the function field

        EXAMPLES::

            sage: K.<x>=FunctionField(GF(5)); _.<Y>=K[]
            sage: L.<y>=K.extension(Y^3-(x^3-1)/(x^3-2))
            sage: L
            Function field in y defined by y^3 + (4*x^3 + 1)/(x^3 + 3)
            sage: L.genus()
            4

        The defining equation needs not be monic::

            sage: K.<x>=FunctionField(GF(4)); _.<Y>=K[]
            sage: L.<y>=K.extension((1-x)*Y^7-x^3)
            sage: L.gaps()
            [1, 2, 3]
        """
        FunctionField_polymod.__init__(self, polynomial, names, element_class=FunctionFieldElement_global)

        self._cache = {}

    @cached_method
    def maximal_order(self):
        """
        Return the maximal order of the function field.

        EXAMPLES::

            sage: K.<x> = FunctionField(GF(2)); _.<Y> = K[]
            sage: L.<y> = K.extension(Y^2 + Y + x + 1/x)
            sage: L.maximal_order()
            Maximal order of Function field in y defined by y^2 + y + (x^2 + 1)/x
        """
        model, from_model, to_model = self.monic_integral_model('z')
        basis = [from_model(g) for g in model.maximal_order().basis()]
        return FunctionFieldMaximalOrder_global(self, basis)

    @cached_method
    def maximal_order_infinite(self):
        """
        Return the maximal infinite order of the function field.

        EXAMPLES::

            sage: K.<x> = FunctionField(GF(2)); _.<t> = K[]
            sage: F.<y> = K.extension(t^3 - x^2*(x^2+x+1)^2)
            sage: F.maximal_order_infinite()
            Maximal infinite order of Function field in y defined by y^3 + x^6 + x^4 + x^2

            sage: K.<x> = FunctionField(GF(2)); _.<Y> = K[]
            sage: L.<y> = K.extension(Y^2 + Y + x + 1/x)
            sage: L.maximal_order_infinite()
            Maximal infinite order of Function field in y defined by y^2 + y + (x^2 + 1)/x
        """
        M,from_M,to_M = self._inversion_isomorphism()
        basis = [from_M(g) for g in M.maximal_order().basis()]
        return FunctionFieldMaximalOrderInfinite_global(self, basis)

    @cached_method
    def _inversion_isomorphism(self):
        """
        Return an inverted function field and isomorphisms.

        Here we define an isomorphism of self with a function field M
        over k(1/x), mapping x :-> 1/x and t to y. We will call this
        inversion isomorphism. Also note that K = k(1/x) == k(x). We
        also require M to be canonical.
        """
        K = self.base_field()
        R = PolynomialRing(K,'T')
        x = K.gen()

        xinv = 1/x
        F_poly = R([K.hom(xinv)(c) for c in self.polynomial().list()])
        F = K.extension(F_poly)

        self2F = self.hom([F.gen(),xinv])
        F2self = F.hom([self.gen(),xinv])

        M,M2F,F2M = F.monic_integral_model('s')

        return M,F2self*M2F,F2M*self2F

    @cached_method
    def place_set(self):
        """
        Return the set of all places of the function field.

        EXAMPLES::

            sage: K.<x> = FunctionField(GF(2)); _.<Y> = K[]
            sage: L.<y> = K.extension(Y^2 + Y + x + 1/x)
            sage: L.place_set()
            Set of places of Function field in y defined by y^2 + y + (x^2 + 1)/x
        """
        from .place import PlaceSet
        return PlaceSet(self)

    @cached_method
    def space_of_differentials(self):
        """
        Return the space of differentials attached to the function field.

        EXAMPLES::

            sage: K.<x> = FunctionField(GF(5)); _.<Y> = K[]
            sage: L.<y> = K.extension(Y^3 - (x^3 - 1)/(x^3 - 2))
            sage: L.space_of_differentials()
            Space of differentials of Function field in y defined by y^3 + (4*x^3 + 1)/(x^3 + 3)
        """
        return DifferentialsSpace(self)

    @cached_method
    def divisor_group(self):
        """
        Return the group of divisors attached to the function field.

        EXAMPLES::

            sage: K.<x> = FunctionField(GF(5)); _.<Y> = K[]
            sage: L.<y> = K.extension(Y^3 - (x^3 - 1)/(x^3 - 2))
            sage: L.divisor_group()
            Divisor group of Function field in y defined by y^3 + (4*x^3 + 1)/(x^3 + 3)
        """
        from .divisor import DivisorGroup
        return DivisorGroup(self)

    @cached_method
    def valuation_ring(self, place):
        from .valuation_ring import FunctionFieldValuationRing_global
        return FunctionFieldValuationRing_global(self, place)

    @cached_method
    def residue_field(self, place):
        """
        Return the residue field associated with the place along with the maps
        from and to the residue field.

        The domain of the map to the residue field is the discrete valuation
        ring associated with the place.

        The discrete valuation ring is defined as the ring of all elements of
        the function field with nonnegative valuation at the place. The maximal
        ideal is the set of elements of positive valuation.  The residue field
        is then the quotient of the discrete valuation ring by its maximal
        ideal.

        If an element not in the valuation ring is applied to the map,
        TypeError is raised.

        EXAMPLES::

            sage: K.<x> = FunctionField(GF(2)); _.<Y> = K[]
            sage: L.<y> = K.extension(Y^2 + Y + x + 1/x)
            sage: p = L.places_finite()[0]
            sage: R, fr_R, to_R = L.residue_field(p)
            sage: R
            Finite Field of size 2
            sage: f = 1 + y
            sage: f.valuation(p)
            -1
            sage: to_R(f)
            Traceback (most recent call last):
            ...
            TypeError: ...
            sage: (1+1/f).valuation(p)
            0
            sage: to_R(1 + 1/f)
            1
            sage: [fr_R(e) for e in R]
            [0, 1]
        """
        return place.valuation_ring().residue_field()

    @cached_method
    def hasse_derivation(self):
        """
        Return the Hasse derivation for the function field.

        EXAMPLES::

            sage: K.<x>=FunctionField(GF(5)); _.<Y>=K[]
            sage: L.<y>=K.extension(Y^3 - (x^3 - 1)/(x^3 - 2))
            sage: L.hasse_derivation()
            Hasse derivation map:
              From: Function field in y defined by y^3 + (4*x^3 + 1)/(x^3 + 3)
              To:   Function field in y defined by y^3 + (4*x^3 + 1)/(x^3 + 3)
        """
        from .maps import FunctionFieldHasseDerivation_global
        return FunctionFieldHasseDerivation_global(self)

    def places(self, degree=1):
        """
        Return a list of the places of the degree of the function field.

        INPUT:

        - ``degree`` -- the degree of the places to be returned

        The infinite places come first.

        EXAMPLES::

            sage: F.<a> = GF(2)
            sage: K.<x> = FunctionField(F)
            sage: R.<t> = PolynomialRing(K)
            sage: L.<y> = K.extension(t^4 + t - x^5)
            sage: L.places(1)
            [Place (1/x, 1/x^4*y^3), Place (x, y), Place (x, y + 1)]
        """
        return self.places_infinite(degree) + self.places_finite(degree)

    def places_finite(self, degree=1):
        """
        Return a list of the finite places of the degree of the function field.

        INPUT:

        - ``degree`` -- the degree of the places to be returned

        EXAMPLES::

            sage: F.<a>=GF(2)
            sage: K.<x>=FunctionField(F)
            sage: R.<t>=PolynomialRing(K)
            sage: L.<y>=K.extension(t^4+t-x^5)
            sage: L.places_finite(1)
            [Place (x, y), Place (x, y + 1)]
        """
        return [place for place in self._places_finite(degree)]

    def _places_finite(self, degree):
        """
        Return a generator of finite places of the function field of ``degree``.

        EXAMPLES::

            sage: F.<a>=GF(2)
            sage: K.<x>=FunctionField(F)
            sage: R.<t>=PolynomialRing(K)
            sage: L.<y>=K.extension(t^4+t-x^5)
            sage: L._places_finite(1)
            <generator object ...>
        """
        O = self.maximal_order()
        K = self.base_field()

        from sage.rings.integer import Integer
        degree = Integer(degree)

        for d in degree.divisors():
            for p in K.places_finite(degree=d):
                for prime,_,_ in O.decomposition(p.prime_ideal()):
                    place = prime.place()
                    if place.degree() == degree:
                        yield place

    def places_infinite(self, degree=1):
        """
        Return a list of the infinite places of the degree of the function
        field.

        INPUT:

        - ``degree`` -- the degree of the places to be returned

        EXAMPLES::

            sage: F.<a>=GF(2)
            sage: K.<x>=FunctionField(F)
            sage: R.<t>=PolynomialRing(K)
            sage: L.<y>=K.extension(t^4+t-x^5)
            sage: L.places_infinite(1)
            [Place (1/x, 1/x^4*y^3)]
        """
        return [place for place in self._places_infinite(degree)]

    def _places_infinite(self, degree):
        """
        Return a generator of `infinite` places of the function field of ``degree``.

        EXAMPLES::

            sage: F.<a>=GF(2)
            sage: K.<x>=FunctionField(F)
            sage: R.<t>=PolynomialRing(K)
            sage: L.<y>=K.extension(t^4+t-x^5)
            sage: L._places_infinite(1)
            <generator object ...>
        """
        Oinf = self.maximal_order_infinite()
        K = self.base_field()
        for prime,_,_ in Oinf.decomposition():
            place = prime.place()
            if place.degree() == degree:
                yield place

    def different(self):
        """
        Return the different divisor of the function field.

        EXAMPLES::

            sage: K.<x> = FunctionField(GF(2)); R.<t> = PolynomialRing(K)
            sage: F.<y> = K.extension(t^3-x^2*(x^2+x+1)^2)
            sage: F.different()
            2*Place (x, (1/(x^3 + x^2 + x))*y^2)
             + 2*Place (x^2 + x + 1, (1/(x^3 + x^2 + x))*y^2)
        """
        O = self.maximal_order()
        Oinf = self.maximal_order_infinite()
        return O.different().divisor() + Oinf.different().divisor()

    def exact_constant_field(self):
        """
        Return the exact constant field and its embedding into the function field.

        EXAMPLES::

            sage: K.<x> = FunctionField(GF(3)); _.<Y> = K[]
            sage: f = Y^2 - x*Y + x^2 + 1 # irreducible but not absolutely irreducible
            sage: L.<y> = K.extension(f)
            sage: L.genus()
            0
            sage: L.exact_constant_field()
            (Finite Field in t of size 3^2, Ring morphism:
               From: Finite Field in t of size 3^2
               To:   Function field in y defined by y^2 + 2*x*y + x^2 + 1
               Defn: t |--> y + x)
            sage: (y+x).divisor()
            0
        """
        from .divisor import zero_divisor

        # A basis of the full constant field is obtained from
        # computing a Riemann-Roch basis of zero divisor.
        basis = zero_divisor(self).basis_function_space()

        dim = len(basis)

        for e in basis:
            _min_poly = e.minimal_polynomial('t')
            if _min_poly.degree() == dim:
                break
        k = self.constant_base_field()
        R = k['t']
        min_poly = R([k(c) for c in _min_poly.list()])

        k_ext = k.extension(min_poly, 't')

        if k_ext.is_prime_field():
            # The cover of the quotient ring k_ext is the integer ring
            # whose generator is 1. This is different from the generator
            # of k_ext.
            embedding = k_ext.hom([self(1)], self)
        else:
            embedding = k_ext.hom([e], self)

        return k_ext, embedding

    def genus(self):
        """
        Return the genus of the function field.

        EXAMPLES::

            sage: F.<a> = GF(16)
            sage: K.<x> = FunctionField(F); K
            Rational function field in x over Finite Field in a of size 2^4
            sage: R.<t> = PolynomialRing(K)
            sage: L.<y> = K.extension(t^4+t-x^5)
            sage: L.genus()
            6

        The genus is computed using the Hurwitz genus formula.
        """
        k,_ = self.exact_constant_field()
        different_degree = self.different().degree() # must be even
        return different_degree // 2 - self.degree() / k.degree() + 1

    def gaps(self):
        """
        Return the gaps of the function field.

        EXAMPLES::

            sage: K.<x> = FunctionField(GF(2)); _.<Y> = K[]
            sage: L.<y> = K.extension(Y^3 + x^3 * Y + x)
            sage: L.gaps()
            [1, 2, 3]
        """
        try:
            return self._cache['gaps']
        except KeyError:
            _, gaps = self._weierstrass_places()
            return gaps

    def weierstrass_places(self):
        """
        Return all Weierstrass places.

        EXAMPLES::

            sage: K.<x> = FunctionField(GF(2)); _.<Y> = K[]
            sage: L.<y> = K.extension(Y^3 + x^3 * Y + x)
            sage: L.weierstrass_places()
            [Place (1/x, 1/x^3*y^2 + 1/x),
             Place (1/x, 1/x^3*y^2 + 1/x^2*y + 1),
             Place (x, y),
             Place (x + 1, (x^3 + 1)*y + x + 1),
             Place (x^3 + x + 1, y + 1),
             Place (x^3 + x + 1, y + x^2),
             Place (x^3 + x + 1, y + x^2 + 1),
             Place (x^3 + x^2 + 1, y + x),
             Place (x^3 + x^2 + 1, y + x^2 + 1),
             Place (x^3 + x^2 + 1, y + x^2 + x + 1)]
        """
        try:
            R = self._cache['ramification_divisor']
        except KeyError:
            R,_ = self._weierstrass_places()
        return R.support()

    def _weierstrass_places(self):
        """
        Return the Weierstrass places together with the gap sequence
        for ordinary places.

        EXAMPLES::

            sage: K.<x> = FunctionField(GF(2)); _.<Y> = K[]
            sage: L.<y> = K.extension(Y^3 + x^3 * Y + x)
            sage: L.weierstrass_places()
            [Place (1/x, 1/x^3*y^2 + 1/x),
             Place (1/x, 1/x^3*y^2 + 1/x^2*y + 1),
             Place (x, y),
             Place (x + 1, (x^3 + 1)*y + x + 1),
             Place (x^3 + x + 1, y + 1),
             Place (x^3 + x + 1, y + x^2),
             Place (x^3 + x + 1, y + x^2 + 1),
             Place (x^3 + x^2 + 1, y + x),
             Place (x^3 + x^2 + 1, y + x^2 + 1),
             Place (x^3 + x^2 + 1, y + x^2 + x + 1)]

        This method implements Algorithm 30 in [Hes2002b]_.
        """
        W = self(self.base_field().gen()).differential().divisor()
        basis = W._basis()
        d = len(basis)

        if d == 0:
            return [],[]

        hasse = self.hasse_derivation()
        M = matrix([basis])
        e = 1
        gaps = [1]
        while M.nrows() < d:
            row = vector([hasse._derive(basis[i], e) for i in range(d)])
            if not row in M.row_space():
                M = matrix(M.rows() + [row])
                gaps.append(e+1)
            e += 1

        # This is faster than M.determinant(). Note that Mx
        # is a matrix over univariate polynomial ring.
        Mx = matrix(M.nrows(), [c._x for c in M.list()])
        detM = self(Mx.determinant() % self._polynomial)

        R = detM.divisor() + sum(gaps)*W # ramification divisor

        self._cache['ramification_divisor'] = R
        self._cache['gaps'] = gaps

        return R, gaps

    def valuation_ring(self, place):
        """
        Return the valuation ring at the place.

        EXAMPLES::

            sage: K.<x> = FunctionField(GF(2)); _.<Y> = K[]
            sage: L.<y> = K.extension(Y^2 + Y + x + 1/x)
            sage: p = L.places_finite()[0]
            sage: L.valuation_ring(p)
            Valuation ring at Place (x, x*y)
        """
        return place.valuation_ring()

class FunctionField_global_integral(FunctionField_global):
    """
    Global function fields defined by an irreducible and separable polynomial
    integral over the maximal order of the base rational function field over a
    finite field.
    """
    @cached_method
    def equation_order(self):
        """
        Return the equation order of the function field.

        EXAMPLES::

            sage: K.<x> = FunctionField(GF(2)); R.<t> = PolynomialRing(K)
            sage: F.<y> = K.extension(t^3-x^2*(x^2+x+1)^2)
            sage: F.equation_order()
            Order in Function field in y defined by y^3 + x^6 + x^4 + x^2
        """
        a = self.gen()
        basis = [a**i for i in range(self.degree())]
        return FunctionFieldOrder_basis(basis)

    @cached_method
    def maximal_order(self):
        """
        Return the maximal order of the function field.

        EXAMPLE::

            sage: K.<x> = FunctionField(GF(2));
            sage: R.<t> = PolynomialRing(K);
            sage: F.<y> = K.extension(t^4 + x^12*t^2 + x^18*t + x^21 + x^18);
            sage: O = F.maximal_order()
            sage: O.basis()
            (1, 1/x^4*y, 1/x^11*y^2 + 1/x^2, 1/x^15*y^3 + 1/x^6*y)

        The basis of the maximal order *always* starts with 1. This is assumed
        in some algorithms.
        """
        from sage.libs.singular.function import singular_function, lib
        from sage.env import SAGE_EXTCODE
        lib(SAGE_EXTCODE + '/singular/function_field/core.lib')
        normalize = singular_function('normalP')

        k = self.constant_base_field()
        K = self.base_field() # rational function field
        n = self.degree()

        # construct the defining polynomial of the function field
        # as a two-variate polynomial g in the ring k[y,x] where
        # k is the constant base field
        S,(y,x) = PolynomialRing(k, names='y,x', order='lex').objgens()
        v = self.polynomial().list()
        g = sum([v[i].numerator().subs(x) * y**i for i in range(len(v))])

        # Singular "normalP" algorithm assumes affine domain over
        # a prime field. So we constuct gflat lifting g as in
        # k_prime[yy,xx,zz]/(k_poly) where k = k_prime[zz]/(k_poly)
        R = PolynomialRing(k.prime_subfield(), names='yy,xx,zz')
        zz = R.gen(2)
        gflat = R.zero()
        for m in g.monomials():
            c = g.monomial_coefficient(m).polynomial(zz)
            gflat += c * R(m) # R(m) is a monomial in yy and xx

        k_poly = k.polynomial(zz)

        # invoke Singular
        pols_in_R = normalize(R.ideal([k_poly, gflat]))

        # reconstruct polynomials in S
        h = R.hom([y,x,k.gen()],S)
        pols_in_S = [h(f) for f in pols_in_R]

        # reconstruct polynomials in the function field
        x = K.gen()
        y = self.gen()
        pols = []
        for f in pols_in_S:
            p = f.polynomial(S.gen(0))
            s = 0
            for i in range(p.degree()+1):
                s += p[i].subs(x) * y**i
            pols.append(s)

        # Now if pols = [g1,g2,...gn,g0], then the g1/g0,g2/g0,...,gn/g0,
        # and g0/g0=1 are the module generators of the integral closure
        # of the equation order Sb = k[xb,yb] in its fraction field,
        # that is, the function field. The integral closure of k[x]
        # is then obtained by multiplying these generators with powers of y
        # as the equation order itself is an integral extension of k[x].
        d = ~ pols[-1]
        _basis = []
        for f in pols:
            b = d * f
            for i in range(n):
                _basis.append(b)
                b *= y

        # Finally we reduce _basis to get a basis over k[x]. This is done
        # of couse by Hermite normal form computation. Here we apply
        # a trick to get a basis that starts with 1 and is ordered in increasing
        # y-degrees. The trick is to use the reversed Hermite normal form.
        # Note that it is important that the overall denominator l lies in k[x].
        V, fr_V, to_V = self.vector_space()
        basis_V = [to_V(b) for b in _basis]
        l = lcm([v.denominator() for v in basis_V])

        # Why 'reversed' here? I don't know. But if we do not do that, it
        # dramatically increases the time it takes to get hermite_form_reversed...
        _mat = matrix([[c.numerator() for c in l*v] for v in reversed(basis_V)])

        mat = _mat.hermite_form_reversed(include_zero_rows=False)

        basis = [fr_V(v) / l for v in mat]
        return FunctionFieldMaximalOrder_global(self, basis)

    @cached_method
    def primitive_integal_element_infinite(self):
        """
        Return a primitive integral element over the base maximal infinite order.

        This element is integral over the maximal infinite order of the base rational
        function field and the function field is a simple extension by this element
        over the base order.

        EXAMPLES::

            sage: K.<x> = FunctionField(GF(2)); R.<t> = PolynomialRing(K)
            sage: F.<y> = K.extension(t^3-x^2*(x^2+x+1)^2)
            sage: b = F.primitive_integal_element_infinite(); b
            1/x^2*y
            sage: b.minimal_polynomial('t')
            t^3 + (x^4 + x^2 + 1)/x^4
        """
        f = self.polynomial()
        n = f.degree()
        y = self.gen()
        x = self.base_field().gen()

        cf = max([(f[i].numerator().degree()/(n-i)).ceil() for i in range(n)
                  if f[i] != 0])
        return y*x**(-cf)

    def equation_order_infinite(self):
        """
        Return the infinite equation order of the function field.

        This is by definition `o[b]` where `b` is the primitive integral
        element from :meth:`primitive_integal_element_infinite()` and
        `o` is the maximal infinite order of the base rational function
        field.

        EXAMPLES::

            sage: K.<x> = FunctionField(GF(2)); R.<t> = PolynomialRing(K)
            sage: F.<y> = K.extension(t^3-x^2*(x^2+x+1)^2)
            sage: F.equation_order_infinite()
            Infinite order in Function field in y defined by y^3 + x^6 + x^4 + x^2
        """
        b = self.primitive_integal_element_infinite()
        basis = [b**i for i in range(self.degree())]
        return FunctionFieldOrderInfinite_basis(basis)

class RationalFunctionField(FunctionField):
    """
    A rational function field K(t) in one variable, over an arbitrary
    base field.

    EXAMPLES::

        sage: K.<t> = FunctionField(GF(3)); K
        Rational function field in t over Finite Field of size 3
        sage: K.gen()
        t
        sage: 1/t + t^3 + 5
        (t^4 + 2*t + 1)/t

    There are various ways to get at the underlying fields and rings
    associated to a rational function field::

        sage: K.<t> = FunctionField(GF(7))
        sage: K.base_field()
        Rational function field in t over Finite Field of size 7
        sage: K.field()
        Fraction Field of Univariate Polynomial Ring in t over Finite Field of size 7
        sage: K.constant_field()
        Finite Field of size 7
        sage: K.maximal_order()
        Maximal order of Rational function field in t over Finite Field of size 7

    We define a morphism::

        sage: K.<t> = FunctionField(QQ)
        sage: L = FunctionField(QQ, 'tbar') # give variable name as second input
        sage: K.hom(L.gen())
        Function Field morphism:
          From: Rational function field in t over Rational Field
          To:   Rational function field in tbar over Rational Field
          Defn: t |--> tbar
    """
    def __init__(self, constant_field, names,
            element_class = FunctionFieldElement_rational,
            category=CAT):
        """
        Create a rational function field in one variable.

        INPUT:

            - ``constant_field`` -- an arbitrary field
            - ``names`` -- a string or tuple of length 1
            - ``category`` -- default: FunctionFields()

        EXAMPLES::

            sage: K.<t> = FunctionField(CC); K
            Rational function field in t over Complex Field with 53 bits of precision
            sage: K.category()
            Category of function fields
            sage: FunctionField(QQ[I], 'alpha')
            Rational function field in alpha over Number Field in I with defining polynomial x^2 + 1

        Must be over a field::

            sage: FunctionField(ZZ, 't')
            Traceback (most recent call last):
            ...
            TypeError: constant_field must be a field
        """
        if names is None:
            raise ValueError("variable name must be specified")
        elif not isinstance(names, tuple):
            names = (names, )
        if not constant_field.is_field():
            raise TypeError("constant_field must be a field")
        self._element_class = element_class
        self._element_init_pass_parent = False
        self._constant_field = constant_field
        FunctionField.__init__(self, self, names=names, category = category)
        R = constant_field[names[0]]
        self._hash = hash((constant_field, names))
        self._ring = R
        self._field = R.fraction_field()
        self._populate_coercion_lists_(coerce_list=[self._field])
        self._gen = self(R.gen())

    def __reduce__(self):
        """
        Returns the arguments which were used to create this instance. The rationale for this is explained in the documentation of ``UniqueRepresentation``.

        EXAMPLES::

            sage: K.<x> = FunctionField(QQ)
            sage: clazz,args = K.__reduce__()
            sage: clazz(*args)
            Rational function field in x over Rational Field
        """
        from .constructor import FunctionField
        return FunctionField, (self._constant_field, self._names)

    def __hash__(self):
        """
        Return hash of the function field.

        The hash is formed from the constant field and the variable names.

        EXAMPLES::

            sage: K.<t> = FunctionField(QQ)
            sage: hash(K) == hash((K.constant_base_field(), K.variable_names()))
            True

        """
        return self._hash

    def _repr_(self):
        """
        Return string representation of the function field.

        EXAMPLES::

            sage: K.<t> = FunctionField(QQ)
            sage: K._repr_()
            'Rational function field in t over Rational Field'
        """
        return "Rational function field in %s over %s"%(
            self.variable_name(), self._constant_field)

    def _element_constructor_(self, x):
        r"""
        Coerce ``x`` into an element of the function field, possibly not canonically.

        INPUT:

            - ``x`` -- the element

        OUTPUT:

            ``x``, as an element of the function field

        EXAMPLES::

            sage: K.<t> = FunctionField(QQ)
            sage: a = K._element_constructor_(K.maximal_order().gen()); a
            t
            sage: a.parent()
            Rational function field in t over Rational Field

        TESTS:

        Conversion of a string::

            sage: K('t')
            t
            sage: K('1/t')
            1/t

        Conversion of a constant polynomial over the function field::

            sage: K(K.polynomial_ring().one())
            1

        Some indirect test of conversion::

            sage: S.<x, y> = K[]
            sage: I = S*[x^2 - y^2, y-t]
            sage: I.groebner_basis()
            [x^2 - t^2, y - t]

        """
        if isinstance(x, FunctionFieldElement):
            return FunctionFieldElement_rational(self, self._field(x._x))
        try:
            x = self._field(x)
        except TypeError as Err:
            try:
                if x.parent() is self.polynomial_ring():
                    return x[0]
            except AttributeError:
                pass
            raise Err
        return FunctionFieldElement_rational(self, x)

    def _to_constant_base_field(self, f):
        r"""
        Return ``f`` as an element of the :meth:`constant_base_field`.

        INPUT:

        - ``f`` -- an element of the rational function field which is a
          constant of the underlying rational function field.

        EXAMPLES::

            sage: K.<x> = FunctionField(QQ)
            sage: K._to_constant_base_field(K(1))
            1
            sage: K._to_constant_base_field(K(x))
            Traceback (most recent call last):
            ...
            ValueError: only constants can be converted into the constant base field but x is not a constant

        TESTS:

        Verify that :trac:`21872` has been resolved::

            sage: K(1) in QQ
            True
            sage: x in QQ
            False

        """
        K = self.constant_base_field()
        if f.denominator() in K and f.numerator() in K:
            # When K is not exact, f.denominator() might not be an exact 1, so
            # we need to divide explicitly to get the correct precision
            return K(f.numerator()) / K(f.denominator())
        raise ValueError("only constants can be converted into the constant base field but %r is not a constant"%(f,))

    def _to_bivariate_polynomial(self, f):
        """
        Convert ``f`` from a univariate polynomial over the rational function
        field into a bivariate polynomial and a denominator.

        INPUT:

            - ``f`` -- a univariate polynomial over self.

        OUTPUT:

            - bivariate polynomial, denominator

        EXAMPLES::

            sage: R.<t> = FunctionField(GF(7))
            sage: S.<X> = R[]
            sage: f = (1/t)*(X^4 - 1/t^2)*(X^3 - t^3)
            sage: R._to_bivariate_polynomial(f)
            (X^7*t^2 - X^4*t^5 - X^3 + t^3, t^3)
        """
        v = f.list()
        from sage.arith.all import lcm
        denom = lcm([a.denominator() for a in v])
        S = denom.parent()
        x,t = S.base_ring()['%s,%s'%(f.parent().variable_name(),self.variable_name())].gens()
        phi = S.hom([t])
        return sum([phi((denom * v[i]).numerator()) * x**i for i in range(len(v))]), denom

    def _factor_univariate_polynomial(self, f, proof=True):
        """
        Factor the univariate polynomial ``f`` over self.

        EXAMPLES::

        We do a factorization over the function field over the rationals::

            sage: R.<t> = FunctionField(QQ)
            sage: S.<X> = R[]
            sage: f = (1/t)*(X^4 - 1/t^2)*(X^3 - t^3)
            sage: f.factor()             # indirect doctest
            (1/t) * (X - t) * (X^2 - 1/t) * (X^2 + 1/t) * (X^2 + t*X + t^2)
            sage: f.factor().prod() == f
            True

        We do a factorization over a finite prime field::

            sage: R.<t> = FunctionField(GF(7))
            sage: S.<X> = R[]
            sage: f = (1/t)*(X^4 - 1/t^2)*(X^3 - t^3)
            sage: f.factor()
            (1/t) * (X + 3*t) * (X + 5*t) * (X + 6*t) * (X^2 + 1/t) * (X^2 + 6/t)
            sage: f.factor().prod() == f
            True

        Factoring over a function field over a non-prime finite field::

            sage: k.<a> = GF(9)
            sage: R.<t> = FunctionField(k)
            sage: S.<X> = R[]
            sage: f = (1/t)*(X^3 - a*t^3)
            sage: f.factor()
            (1/t) * (X + (a + 2)*t)^3
            sage: f.factor().prod() == f
            True
        """
        F, d = self._to_bivariate_polynomial(f)
        fac = F.factor()
        x = f.parent().gen()
        t = f.parent().base_ring().gen()
        phi = F.parent().hom([x, t])
        v = [(phi(P),e) for P, e in fac]
        unit = phi(fac.unit())/d
        w = []
        for a, e in v:
            c = a.leading_coefficient()
            a = a/c
            unit *= (c**e)
            w.append((a,e))
        from sage.structure.factorization import Factorization
        return Factorization(w, unit=unit)

    def extension(self, f, names=None):
        """
        Create an extension L = K[y]/(f(y)) of the rational function field.

        INPUT:

        - ``f`` -- a univariate polynomial over self

        - ``names`` -- None or string or length-1 tuple

        OUTPUT:

        - a function field

        EXAMPLES::

            sage: K.<x> = FunctionField(QQ); R.<y> = K[]
            sage: K.extension(y^5 - x^3 - 3*x + x*y)
            Function field in y defined by y^5 + x*y - x^3 - 3*x

        A nonintegral defining polynomial::

            sage: K.<t> = FunctionField(QQ); R.<y> = K[]
            sage: K.extension(y^3 + (1/t)*y + t^3/(t+1))
            Function field in y defined by y^3 + 1/t*y + t^3/(t + 1)

        The defining polynomial need not be monic or integral::

            sage: K.extension(t*y^3 + (1/t)*y + t^3/(t+1))
            Function field in y defined by t*y^3 + 1/t*y + t^3/(t + 1)
        """
        from . import constructor
        return constructor.FunctionFieldExtension(f, names)

    @cached_method
    def polynomial_ring(self, var='x'):
        """
        Return a polynomial ring in one variable over the rational function field.

        INPUT:

        - ``var`` -- (default: 'x') a string; the name of the variable

        EXAMPLES::

            sage: K.<x> = FunctionField(QQ)
            sage: K.polynomial_ring()
            Univariate Polynomial Ring in x over Rational function field in x over Rational Field
            sage: K.polynomial_ring('T')
            Univariate Polynomial Ring in T over Rational function field in x over Rational Field
        """
        return self[var]

    @cached_method
    def vector_space(self):
        """
        Return a vector space V and isomorphisms self --> V and V --> self.

        OUTPUT:

            -  ``V`` -- a vector space over the rational numbers
            -  ``from_V`` -- an isomorphism from V to self
            -  ``to_V`` -- an isomorphism from self to V

        EXAMPLES::

            sage: K.<x> = FunctionField(QQ)
            sage: K.vector_space()
            (Vector space of dimension 1 over Rational function field in x over Rational Field, Isomorphism morphism:
              From: Vector space of dimension 1 over Rational function field in x over Rational Field
              To:   Rational function field in x over Rational Field, Isomorphism morphism:
              From: Rational function field in x over Rational Field
              To:   Vector space of dimension 1 over Rational function field in x over Rational Field)
        """
        V = self.base_field()**1
        from .maps import MapVectorSpaceToFunctionField, MapFunctionFieldToVectorSpace
        from_V = MapVectorSpaceToFunctionField(V, self)
        to_V   = MapFunctionFieldToVectorSpace(self, V)
        return (V, from_V, to_V)

    def random_element(self, *args, **kwds):
        """
        Create a random element of the rational function field.

        Parameters are passed to the random_element method of the
        underlying fraction field.

        EXAMPLES::

            sage: FunctionField(QQ,'alpha').random_element()   # random
            (-1/2*alpha^2 - 4)/(-12*alpha^2 + 1/2*alpha - 1/95)
        """
        return self(self._field.random_element(*args, **kwds))

    def degree(self, base=None):
        """
        Return the degree over the base field of the rational
        function field. Since the base field is the rational function
        field itself, the degree is 1.

        INPUT:

        - ``base`` -- must be this field or ``None``; this parameter is ignored
          and exists to resemble the interface of
          :meth:`FunctionField_polymod.degree`.

        EXAMPLES::

            sage: K.<t> = FunctionField(QQ)
            sage: K.degree()
            1
        """
        if base is None:
            base = self
        if base is not self:
            raise ValueError("base must be None or the rational function field")
        from sage.rings.integer_ring import ZZ
        return ZZ(1)

    def gen(self, n=0):
        """
        Return the ``n``-th generator of the function field.  If ``n`` is not
        0, then an IndexError is raised.

        EXAMPLES::

            sage: K.<t> = FunctionField(QQ); K.gen()
            t
            sage: K.gen().parent()
            Rational function field in t over Rational Field
            sage: K.gen(1)
            Traceback (most recent call last):
            ...
            IndexError: Only one generator.
        """
        if n != 0:
            raise IndexError("Only one generator.")
        return self._gen

    def ngens(self):
        """
        Return the number of generators, which is 1.

        EXAMPLES::

            sage: K.<t> = FunctionField(QQ)
            sage: K.ngens()
            1
        """
        return 1

    def base_field(self):
        """
        Return the base field of the rational function field, which is just
        the function field itself.

        EXAMPLES::

            sage: K.<t> = FunctionField(GF(7))
            sage: K.base_field()
            Rational function field in t over Finite Field of size 7
        """
        return self

    def hom(self, im_gens, base_morphism=None):
        """
        Create a homomorphism from self to another function field.

        INPUT:

        - ``im_gens`` -- exactly one element of some function field
        - ``base_morphism`` -- ignored

        OUTPUT:

        - a map between function fields

        EXAMPLES:

        We make a map from a rational function field to itself::

            sage: K.<x> = FunctionField(GF(7))
            sage: K.hom( (x^4 + 2)/x)
            Function Field endomorphism of Rational function field in x over Finite Field of size 7
              Defn: x |--> (x^4 + 2)/x

        We construct a map from a rational function field into a
        non-rational extension field::

            sage: K.<x> = FunctionField(GF(7)); R.<y> = K[]
            sage: L.<y> = K.extension(y^3 + 6*x^3 + x)
            sage: f = K.hom(y^2 + y  + 2); f
            Function Field morphism:
              From: Rational function field in x over Finite Field of size 7
              To:   Function field in y defined by y^3 + 6*x^3 + x
              Defn: x |--> y^2 + y + 2
            sage: f(x)
            y^2 + y + 2
            sage: f(x^2)
            5*y^2 + (x^3 + 6*x + 4)*y + 2*x^3 + 5*x + 4
        """
        from sage.structure.category_object import CategoryObject
        if isinstance(im_gens, CategoryObject):
            return self.Hom(im_gens).natural_map()
        if not isinstance(im_gens, (list,tuple)):
            im_gens = [im_gens]
        if len(im_gens) != 1:
            raise ValueError("there must be exactly one generator")
        x = im_gens[0]
        from .maps import FunctionFieldMorphism_rational
        return FunctionFieldMorphism_rational(self.Hom(x.parent()), x)

    def field(self):
        """
        Return the underlying field, forgetting the function field
        structure.

        EXAMPLES::

            sage: K.<t> = FunctionField(GF(7))
            sage: K.field()
            Fraction Field of Univariate Polynomial Ring in t over Finite Field of size 7
        """
        return self._field

    @cached_method
    def place_set(self):
        """
        Return the set of all places of the function field.
        """
        from .place import PlaceSet
        return PlaceSet(self)

    @cached_method
    def maximal_order(self):
        """
        Return the maximal order of the function field.

        Since this is a rational function field it is of the form `K(t)`, and the
        maximal order is by definition `K[t]`, where `K` is the constant field.

        EXAMPLES::

            sage: K.<t> = FunctionField(QQ)
            sage: K.maximal_order()
            Maximal order of Rational function field in t over Rational Field
            sage: K.equation_order()
            Maximal order of Rational function field in t over Rational Field
        """
        return FunctionFieldMaximalOrder_rational(self)

    equation_order = maximal_order

    @cached_method
    def maximal_order_infinite(self):
        """
        Return the maximal infinite order of the function field.

        By definition, this is the valuation ring of the degree valuation of
        the rational function field.

        EXAMPLES::

            sage: K.<t> = FunctionField(QQ)
            sage: K.maximal_order_infinite()
            Maximal infinite order of Rational function field in t over Rational Field
            sage: K.equation_order_infinite()
            Maximal infinite order of Rational function field in t over Rational Field
        """
        return FunctionFieldMaximalOrderInfinite_rational(self)

    equation_order_infinite = maximal_order_infinite

    def constant_base_field(self):
        """
        Return the field of which the rational function field is a
        transcendental extension.

        EXAMPLES::

            sage: K.<t> = FunctionField(QQ)
            sage: K.constant_field()
            Rational Field
        """
        return self._constant_field

    constant_field = constant_base_field

    def different(self):
        """
        Return the different ideal of the rational function field.

        In this case, the different is simply the zero divisor.

        EXAMPLES::

            sage: K.<t> = FunctionField(QQ)
            sage: K.different()
            0
        """
        return self.divisor_group()(0)

    def genus(self):
        """
        Return the genus of the function field, namely 0.

        EXAMPLES::

            sage: K.<x> = FunctionField(QQ)
            sage: K.genus()
            0
        """
        return 0

    @cached_method(key=lambda self, base: None)
    def vector_space(self, base=None):
        """
        Return a vector space `V` and isomorphisms from the field to `V` and
        from `V` to the field.

        This function allows us to identify the elements of this field with
        elements of a one-dimensional vector space over the field itself. This
        method exists so that all function fields (rational or not) have the
        same interface.

        INPUT:

        - ``base`` -- must be the field or None (default: None); this
          parameter is ignored and merely exists to have the same interface as
          :meth:`FunctionField_polymod.vector_space`.

        OUTPUT:

        - ``V`` -- a vector space over base field

        - ``from_V`` -- an isomorphism from `V` to the field

        - ``to_V`` -- an isomorphism from the field to `V`

        EXAMPLES::

            sage: K.<x> = FunctionField(QQ)
            sage: K.vector_space()
            (Vector space of dimension 1 over Rational function field in x over Rational Field, Isomorphism morphism:
              From: Vector space of dimension 1 over Rational function field in x over Rational Field
              To:   Rational function field in x over Rational Field, Isomorphism morphism:
              From: Rational function field in x over Rational Field
              To:   Vector space of dimension 1 over Rational function field in x over Rational Field)

        TESTS::

            sage: K.vector_space()
            (Vector space of dimension 1 over Rational function field in x over Rational Field, Isomorphism morphism:
              From: Vector space of dimension 1 over Rational function field in x over Rational Field
              To:   Rational function field in x over Rational Field, Isomorphism morphism:
              From: Rational function field in x over Rational Field
              To:   Vector space of dimension 1 over Rational function field in x over Rational Field)

        """
        from .maps import MapVectorSpaceToFunctionField, MapFunctionFieldToVectorSpace
        if base is None:
            base = self
        if base is not self:
            raise ValueError("base must be the rational function field or None")
        V = base**1
        from_V = MapVectorSpaceToFunctionField(V, self)
        to_V   = MapFunctionFieldToVectorSpace(self, V)
        return (V, from_V, to_V)

    @cached_method
    def derivation(self):
        r"""
        Return a derivation of the rational function field over the constant
        base field.

        OUTPUT:

        - a derivation of the rational function field

        The derivation maps the generator of the rational function field to 1.

        EXAMPLES::

            sage: K.<x> = FunctionField(GF(3))
            sage: m = K.derivation(); m
            Derivation map:
              From: Rational function field in x over Finite Field of size 3
              To:   Rational function field in x over Finite Field of size 3
            sage: m(x)
            1

        TESTS::

            sage: L.<y> = FunctionField(K)
            sage: L.derivation()
            Traceback (most recent call last):
            ...
            NotImplementedError: not implemented for non-perfect base fields
        """
        from .maps import FunctionFieldDerivation_rational
        if not self.constant_base_field().is_perfect():
            raise NotImplementedError("not implemented for non-perfect base fields")
        return FunctionFieldDerivation_rational(self, self.one())

    @cached_method
    def space_of_differentials(self):
        """
        Return the space of differentials of the rational function field.

        EXAMPLES::

            sage: K.<t> = FunctionField(QQ)
            sage: K.space_of_differentials()
            Space of differentials of Rational function field in t over Rational Field
        """
        return DifferentialsSpace(self)

    @cached_method
    def divisor_group(self):
        """
        Return the group of divisors of the rational function field.

        EXAMPLES::

            sage: K.<t> = FunctionField(QQ)
            sage: K.divisor_group()
            Divisor group of Rational function field in t over Rational Field
        """
        from .divisor import DivisorGroup
        return DivisorGroup(self)

class RationalFunctionField_global(RationalFunctionField):
    """
    Rational function field over finite fields.
    """
    def places(self, degree=1):
        """
        Return all places of the degree.

        INPUT:

        - ``degree`` -- (default: 1) the degree of the places

        EXAMPLES::

            sage: F.<x> = FunctionField(GF(5))
            sage: F.places()
            [Place (1/x),
             Place (x),
             Place (x + 1),
             Place (x + 2),
             Place (x + 3),
             Place (x + 4)]
        """
        if degree == 1:
            return [self.place_infinite()] + self.places_finite(degree)
        else:
            return self.places_finite(degree)

    def places_finite(self, degree=1):
        """
        Return the finite places of the degree.

        EXAMPLES::

            sage: F.<x> = FunctionField(GF(5))
            sage: F.places_finite()
            [Place (x), Place (x + 1), Place (x + 2), Place (x + 3), Place (x + 4)]
        """
        return [place for place in self._places_finite(degree)]

    def _places_finite(self, degree=1):
        """
        Return a generator for all monic irreducible polynomials.

        EXAMPLES::

            sage: F.<x> = FunctionField(GF(5))
            sage: F._places_finite()
            <generator object _places_finite at ...>
        """
        O = self.maximal_order()
        R = O._ring
        G = R.polynomials(of_degree=degree)
        for g in G:
            if not (g.is_monic() and g.is_irreducible()):
                continue
            yield O.ideal(g).place()

    def place_infinite(self):
        """
        Return the unique place at infinite.

        EXAMPLES::

            sage: F.<x> = FunctionField(GF(5))
            sage: F.place_infinite()
            Place (1/x)
        """
        return self.maximal_order_infinite().prime_ideal().place()

    def residue_field(self, place):
        """
        Return the residue field of the place along with the maps from
        and to it.

        INPUT:

        - ``place`` -- a place of the function field

        EXAMPLES::

            sage: F.<x> = FunctionField(GF(5))
            sage: p = F.places_finite(2)[0]
            sage: R, fr_R, to_R = F.residue_field(p)
            sage: R
            Finite Field in z2 of size 5^2
            sage: to_R(x) in R
            True
        """
        prime = place.prime_ideal()

        if place.is_infinite_place():
            K = self.constant_base_field()

            def from_K(e):
                return self(e)

            def to_K(f):
                n = f.numerator()
                d = f.denominator()

                n_deg = n.degree()
                d_deg =d.degree()

                if n_deg < d_deg:
                    return K(0)
                elif n_deg == d_deg:
                    return n.lc() / d.lc()
                else:
                    raise TypeError("not in the valuation ring")
        else:
            O = self.maximal_order()
            K, _from_K, _to_K = O.residue_field(prime)

            def from_K(e):
                return _from_K(e)

            def to_K(f):
                if f in O: # f.denominator() is 1
                    return _to_K(f.numerator())
                else:
                    d = self(f.denominator())
                    n = d * f

                    nv = prime.valuation(O.ideal(n))
                    dv = prime.valuation(O.ideal(d))

                    if nv > dv:
                        return K(0)
                    elif dv > nv:
                        raise TypeError("not in the valuation ring")

                    s = ~prime.gen()
                    rd = d * s**dv # in O but not in prime
                    rn = n * s**nv # in O but not in prime
                    return to_K(rn) / to_K(rd)

        mor_from_K = SetMorphism(Hom(K,self), from_K )
        mor_to_K = SetMorphism(Hom(self,K), to_K )
        return K, mor_from_K, mor_to_K

    @cached_method
    def hasse_derivation(self):
        """
        Return the Hasse derivation for the function field.

        EXAMPLES::

            sage: F.<x> = FunctionField(GF(5))
            sage: d = F.hasse_derivation()
            sage: [d(x^5,i) for i in range(10)]
            [x^5, 0, 0, 0, 0, 1, 0, 0, 0, 0]
            sage: [d(x^7,i) for i in range(10)]
            [x^7, 2*x^6, x^5, 0, 0, x^2, 2*x, 1, 0, 0]
        """
        from .maps import FunctionFieldHasseDerivation_rational
        return FunctionFieldHasseDerivation_rational(self)

    def valuation_ring(self, place):
        """
        Return the valuation ring at the place.

        EXAMPLES::

            sage: F.<x> = FunctionField(GF(2))
            sage: p = F.place_infinite()
            sage: F.valuation_ring(p)
            Valuation ring at Place (1/x)
        """
        return place.valuation_ring()
