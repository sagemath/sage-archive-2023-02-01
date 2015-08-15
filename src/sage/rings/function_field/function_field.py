"""
Function Fields

AUTHORS:

- William Stein (2010): initial version

- Robert Bradshaw (2010-05-30): added is_finite()

- Julian Rueth (2011-06-08, 2011-09-14, 2014-06-23, 2014-06-24, 2014-06-27):
  fixed hom(), extension(); use @cached_method; added derivation(); added
  support for relative vector spaces, added derivation(); added
  change_variable_name(), separable_model()

- Maarten Derickx (2011-09-11): added doctests

- Syed Ahmad Lavasani (2011-12-16): added genus(), is_RationalFunctionField()

- Simon King (2014-10-29): Use the same generator names for a function field
  extension and the underlying polynomial ring.

EXAMPLES:

We create an extension of a rational function fields, and do some
simple arithmetic in it::

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
that arithmetic with a tower of 3 fields is fully supported::

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

TESTS::

    sage: TestSuite(K).run()
    sage: TestSuite(L).run()  # long time (8s on sage.math, 2012)
    sage: TestSuite(M).run()  # long time (52s on sage.math, 2012)

The following two test suites do not pass ``_test_elements`` yet since
``R.an_element()`` has a ``_test_category`` method wich it should not have.
It is not the fault of the function field code so this will
be fixed in another ticket::

    sage: TestSuite(R).run(skip = '_test_elements')
    sage: TestSuite(S).run(skip = '_test_elements')
"""
#*****************************************************************************
#       Copyright (C) 2010 William Stein <wstein@gmail.com>
#       Copyright (C) 2010 Robert Bradshaw <robertwb@math.washington.edu>
#       Copyright (C) 2011-2014 Julian Rueth <julian.rueth@gmail.com>
#       Copyright (C) 2011 Maarten Derickx <m.derickx.student@gmail.com>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#  as published by the Free Software Foundation; either version 2 of
#  the License, or (at your option) any later version.
#                  http://www.gnu.org/licenses/
#*****************************************************************************

from sage.rings.ring import Field
from function_field_element import FunctionFieldElement, FunctionFieldElement_rational, FunctionFieldElement_polymod

from sage.misc.cachefunc import cached_method

#is needed for genus computation
from sage.rings.polynomial.polynomial_ring_constructor import PolynomialRing
from sage.interfaces.all import singular

from sage.categories.function_fields import FunctionFields
CAT = FunctionFields()

def is_FunctionField(x):
    """
    Return True if ``x`` is of function field type.

    EXAMPLES::

        sage: from sage.rings.function_field.function_field import is_FunctionField
        sage: is_FunctionField(QQ)
        False
        sage: is_FunctionField(FunctionField(QQ,'t'))
        True
    """
    if isinstance(x, FunctionField): return True
    return x in FunctionFields()

class FunctionField(Field):
    """
    The abstract base class for all function fields.

    EXAMPLES::

        sage: K.<x> = FunctionField(QQ)
        sage: isinstance(K, sage.rings.function_field.function_field.FunctionField)
        True
    """
    def is_perfect(self):
        r"""
        Return whether this field is perfect, i.e., its characteristic is `p=0`
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
        Return the characteristic of this function field.

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
        Return whether this function field is finite, which it is not.

        EXAMPLES::

            sage: R.<t> = FunctionField(QQ)
            sage: R.is_finite()
            False
            sage: R.<t> = FunctionField(GF(7))
            sage: R.is_finite()
            False
        """
        return False

    def extension(self, f, names=None):
        """
        Create an extension L = K[y]/(f(y)) of a function field,
        defined by a univariate polynomial in one variable over this
        function field K.

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
        from constructor import FunctionField_polymod as FunctionField_polymod_Constructor
        return FunctionField_polymod_Constructor(f, names)

    def order_with_basis(self, basis, check=True):
        """
        Return the order with given basis over the maximal order of
        the base field.

        INPUT:

           - ``basis`` -- a list of elements of self
           - ``check`` -- bool (default: True); if True, check that
             the basis is really linearly independent and that the
             module it spans is closed under multiplication, and
             contains the identity element.

        OUTPUT:

            - an order in this function field

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
        from function_field_order import FunctionFieldOrder_basis
        return FunctionFieldOrder_basis([self(a) for a in basis], check=check)

    def order(self, x, check=True):
        """
        Return the order in this function field generated over the
        maximal order by x or the elements of x if x is a list.

        INPUT:

           - ``x`` -- element of self, or a list of elements of self
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
        return self.order_with_basis(basis, check=check)

    def _coerce_map_from_(self, R):
        """
        Return True if there is a coerce map from R to self.

        EXAMPLES::

            sage: K.<x> = FunctionField(QQ); R.<y> = K[]; L.<y> = K.extension(y^3 + x^3 + 4*x + 1)
            sage: L.equation_order()
            Order in Function field in y defined by y^3 + x^3 + 4*x + 1
            sage: L._coerce_map_from_(L.equation_order())
            True
            sage: L._coerce_map_from_(GF(7))
            False
        """
        from function_field_order import FunctionFieldOrder
        if isinstance(R, FunctionFieldOrder) and R.fraction_field() == self:
            return True
        return False

    def _intermediate_fields(self, base):
        r"""
        Return the fields which lie in between ``base`` and this field in the
        tower of function fields.

        INPUT:

        - ``base`` -- a function field, either this field or a field from which
          this field has been created as an extension

        OUTPUT:

        A list of fields. The first entry is ``base``, the last entry is this field.

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
            [Function field in z defined by z^2 - y, Function field in y defined by y^2 - x, Rational function field in x over Rational Field]

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

    def rational_function_field(self):
        r"""
        Return the rational function field from which this field has been
        created as an extension.

        EXAMPLES::

            sage: K.<x> = FunctionField(QQ)
            sage: K.rational_function_field()
            Rational function field in x over Rational Field

            sage: R.<y> = K[]
            sage: L.<y> = K.extension(y^2-x)
            sage: L.rational_function_field()
            Rational function field in x over Rational Field

            sage: R.<z> = L[]
            sage: M.<z> = L.extension(z^2-y)
            sage: M.rational_function_field()
            Rational function field in x over Rational Field
        """
        return self if is_RationalFunctionField(self) else self.base_field().rational_function_field()

class FunctionField_polymod(FunctionField):
    """
    A function field defined by a univariate polynomial, as an
    extension of the base field.

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

        It is not checked if the polynomial used to define this function field is irreducible
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
    def __init__(self, polynomial, names,
            element_class = FunctionFieldElement_polymod,
            category=CAT):
        """
        Create a function field defined as an extension of another
        function field by adjoining a root of a univariate polynomial.

        INPUT:

            - ``polynomial`` -- a univariate polynomial over a function field
            - ``names`` -- variable names (as a tuple of length 1 or string)
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
            TypeError: unable to evaluate 'x' in Fraction Field of Univariate Polynomial Ring in t over Rational Field
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

        Field.__init__(self, base_field,
                                names=names, category = category)

        self._hash = hash(polynomial)
        self._ring = self._polynomial.parent()
        self._populate_coercion_lists_(coerce_list=[base_field, self._ring])
        self._gen = self(self._ring.gen())

    def __reduce__(self):
        """
        Returns the arguments which were used to create this instance. The rationale for this is explained in the documentation of ``UniqueRepresentation``.

        EXAMPLES::

            sage: K.<x> = FunctionField(QQ); R.<y> = K[]
            sage: L = K.extension(y^2 - x)
            sage: clazz,args = L.__reduce__()
            sage: clazz(*args)
            Function field in y defined by y^2 - x
        """
        from constructor import FunctionField_polymod as FunctionField_polymod_Constructor
        return  FunctionField_polymod_Constructor, (self._polynomial, self._names)

    def __hash__(self):
        """
        Return hash of this function field.

        The hash value is equal to the hash of the defining polynomial.

        EXAMPLES::

            sage: K.<x> = FunctionField(QQ); R.<y> = K[]
            sage: L = K.extension(y^5 - x^3 - 3*x + x*y)
            sage: hash(L) == hash(L.polynomial())
            True

        """
        return self._hash

    def is_separable(self, base=None):
        r"""
        Return whether this is a separable extension of ``base``.

        INPUT:

        - ``base`` -- a function field from which this field has been created
          as an extension or ``None`` (default: ``None``); if ``None``, then
          return whether this is a separable extension over its base field.

        EXAMPLES::

            sage: K.<x> = FunctionField(GF(2))
            sage: R.<y> = K[]
            sage: L.<y> = K.extension(y^2 - x)
            sage: L.is_separable()
            False
            sage: R.<z> = L[]
            sage: M.<z> = L.extension(z^3 - y)
            sage: M.is_separable()
            True
            sage: M.is_separable(K)
            False

        """
        if base is None:
            base = self.base_field()
        for k in self._intermediate_fields(base)[:-1]:
            f = k.polynomial()
            f /= f.leading_coefficient()
            if not f.gcd(f.derivative()).is_one():
                return False
        return True

    def monic_integral_model(self, names=None):
        """
        Return a function field isomorphic to this field but which is an
        extension of a rational function field with defining polynomial that is
        monic and integral over the constant base field.

        INPUT:

        - ``names`` -- a string, a tuple of up to two strings, or ``None``
          (default: ``None``), the name of the generator of the field, and
          the name of the generator of the underlying rational function
          field (if a tuple); if ``None``, then the names are chosen
          automatically.

        OUTPUT:

        A triple ``(F,f,t)`` where ``F`` is a function field, ``f`` is an
        isomorphism from ``F`` to this field, and ``t`` is the inverse of
        ``f``.

        EXAMPLES::

            sage: K.<x> = FunctionField(QQ)
            sage: R.<y> = K[]
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
                    x |--> x
            sage: to_A
            Function Field morphism:
              From: Function field in y defined by x^2*y^5 - 1/x
              To:   Function field in z defined by z^5 - x^12
              Defn: y |--> 1/x^3*z
                    x |--> x
            sage: to_A(y)
            1/x^3*z
            sage: from_A(to_A(y))
            y
            sage: from_A(to_A(1/y))
            x^3*y^4
            sage: from_A(to_A(1/y)) == 1/y
            True

        This also works for towers of function fields::

            sage: R.<z> = L[]
            sage: M.<z> = L.extension(z^2*y - 1/x)
            sage: M.monic_integral_model()
            (Function field in z_ defined by x^10 - x^18, Function Field morphism:
              From: Function field in z_ defined by x^10 - x^18
              To:   Function field in z defined by y*z^2 - 1/x
              Defn: z_ |--> x^2*z
                    x |--> x, Function Field morphism:
              From: Function field in z defined by y*z^2 - 1/x
              To:   Function field in z_ defined by x^10 - x^18
              Defn: z |--> 1/x^2*z_
                    y |--> 1/x^15*z_^8
                    x |--> x)

        TESTS:

        If the field is already a monic integral extension, then it is returned
        unchanged::

            sage: K.<x> = FunctionField(QQ)
            sage: R.<y> = K[]
            sage: L.<y> = K.extension(y^2-x)
            sage: L.monic_integral_model()
            (Function field in y defined by y^2 - x, Function Field endomorphism of Function field in y defined by y^2 - x
              Defn: y |--> y
                    x |--> x, Function Field endomorphism of Function field in y defined by y^2 - x
              Defn: y |--> y
                    x |--> x)

        Unless ``names`` does not match with the current names::

            sage: L.monic_integral_model(names=('yy','xx'))
            (Function field in yy defined by y^2 - xx, Function Field morphism:
              From: Function field in yy defined by y^2 - xx
              To:   Function field in y defined by y^2 - x
              Defn: yy |--> y
                    xx |--> x, Function Field morphism:
              From: Function field in y defined by y^2 - x
              To:   Function field in yy defined by y^2 - xx
              Defn: y |--> yy
                    x |--> xx)

        """
        if names is None:
            pass
        elif not isinstance(names, tuple):
            names = (names,)
        elif len(names) > 2:
            raise ValueErorr("names must contain at most 2 entries")

        if self.base_field() is not self.rational_function_field():
            L,from_L,to_L = self.simple_model()
            ret,ret_to_L,L_to_ret = L.monic_integral_model(names)
            from_ret = ret.hom( [from_L(ret_to_L(ret.gen())), from_L(ret_to_L(ret.base_field().gen()))] )
            to_ret = self.hom( [L_to_ret(to_L(k.gen())) for k in self._intermediate_fields(self.rational_function_field())] )
            return ret, from_ret, to_ret
        else:
            if self.polynomial().is_monic() and all([c.denominator().is_one() for c in self.polynomial()]):
                # self is already monic and integral
                if names is None or names == ():
                    names = (self.variable_name(),)
                return self.change_variable_name(names)
            else:
                if names is None or names == ():
                    names = (self.variable_name()+"_",)
                if len(names) == 1:
                    names = (names[0], self.rational_function_field().variable_name())

                g, d = self._make_monic_integral(self.polynomial())
                K,from_K,to_K = self.base_field().change_variable_name(names[1])
                g = g.map_coefficients(to_K)
                ret = K.extension(g, names=names[0])
                from_ret = ret.hom([self.gen() * d, self.base_field().gen()])
                to_ret = self.hom([ret.gen() / d, ret.base_field().gen()])
                return ret, from_ret, to_ret

    def _make_monic_integral(self, f):
        r"""
        Let y be a root of ``f``.  This function returns a monic
        integral polynomial g and an element d of the base field such
        that g(y*d)=0.

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
        from sage.rings.arith import lcm
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
        field in this function field.

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
        Return the degree of this function field over the function field
        ``base``.

        INPUT:

        - ``base`` -- a function field or ``None`` (default: ``None``), a
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
        Return string representation of this function field.

        EXAMPLES::

            sage: K.<x> = FunctionField(QQ); R.<y> = K[]
            sage: L.<y> = K.extension(y^5 - (x^3 + 2*x*y + 1/x))
            sage: L._repr_()
            'Function field in y defined by y^5 - 2*x*y + (-x^4 - 1)/x'
        """
        return "Function field in %s defined by %s"%(self.variable_name(), self._polynomial)

    def base_field(self):
        """
        Return the base field of this function field.  This function
        field is presented as L = K[y]/(f(y)), and the base field is
        by definition the field K.

        EXAMPLES::

            sage: K.<x> = FunctionField(QQ); R.<y> = K[]
            sage: L.<y> = K.extension(y^5 - (x^3 + 2*x*y + 1/x))
            sage: L.base_field()
            Rational function field in x over Rational Field
        """
        return self._base_field

    def random_element(self, *args, **kwds):
        """
        Create a random element of this function field.  Parameters
        are passed onto the random_element method of the base_field.

        EXAMPLES::

            sage: K.<x> = FunctionField(QQ); R.<y> = K[]
            sage: L.<y> = K.extension(y^2 - (x^2 + x))
            sage: L.random_element() # random
            ((x^2 - x + 2/3)/(x^2 + 1/3*x - 1))*y^2 + ((-1/4*x^2 + 1/2*x - 1)/(-5/2*x + 2/3))*y + (-1/2*x^2 - 4)/(-12*x^2 + 1/2*x - 1/95)
        """
        return self(self._ring.random_element(degree=self.degree(), *args, **kwds))

    def polynomial(self):
        """
        Return the univariate polynomial that defines this function
        field, i.e., the polynomial f(y) so that this function field
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
        Return the polynomial ring used to represent elements of this
        function field.  If we view this function field as being presented
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
        Return a vector space `V` and isomorphisms from this field to `V` and
        from `V` to this field.

        This function allows us to identify the elements of this field with
        elements of a vector space over the base field, which is useful for
        representation and arithmetic with orders, ideals, etc.

        INPUT:

        - ``base`` -- a function field or ``None`` (default: ``None``), the
          returned vector space is over ``base`` which defaults to the base
          field of this function field.

        OUTPUT:

        - ``V`` -- a vector space over base field
        - ``from_V`` -- an isomorphism from V to this field
        - ``to_V`` -- an isomorphism from this field to V

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
        from maps import MapVectorSpaceToFunctionField, MapFunctionFieldToVectorSpace
        if base is None:
            base = self.base_field()
        degree = self.degree(base)
        V = base**degree;
        from_V = MapVectorSpaceToFunctionField(V, self)
        to_V   = MapFunctionFieldToVectorSpace(self, V)
        return (V, from_V, to_V)

    def maximal_order(self):
        """
        Return the maximal_order of self.  If we view self as L =
        K[y]/(f(y)), then this is the ring of elements of L that are
        integral over K.

        EXAMPLES:

        This is not yet implemented...::

            sage: K.<x> = FunctionField(QQ); R.<y> = K[]
            sage: L.<y> = K.extension(y^5 - (x^3 + 2*x*y + 1/x))
            sage: L.maximal_order()
            Traceback (most recent call last):
            ...
            NotImplementedError
        """
        raise NotImplementedError

    def _element_constructor_(self, x):
        r"""
        Make ``x`` into an element of this function field, possibly not canonically.

        INPUT:

            - ``x`` -- the element

        OUTPUT:

            ``x``, as an element of this function field

        TESTS::

            sage: K.<x> = FunctionField(QQ); R.<y> = K[]
            sage: L.<y> = K.extension(y^5 - (x^3 + 2*x*y + 1/x))
            sage: L._element_constructor_(L.polynomial_ring().gen())
            y
        """
        if isinstance(x, FunctionFieldElement):
            return FunctionFieldElement_polymod(self, self._ring(x.element()))
        return FunctionFieldElement_polymod(self, self._ring(x))

    def gen(self, n=0):
        """
        Return the ``n``-th generator of this function field. By default ``n`` is 0; any other
        value of ``n`` leads to an error. The generator is the class of y, if we view
        self as being presented as K[y]/(f(y)).

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
        Return the number of generators of this function field over
        its base field.  This is by definition 1.

        EXAMPLES::

            sage: K.<x> = FunctionField(QQ); R.<y> = K[]
            sage: L.<y> = K.extension(y^5 - (x^3 + 2*x*y + 1/x))
            sage: L.ngens()
            1
        """
        return 1

    def equation_order(self):
        """
        If we view self as being presented as K[y]/(f(y)), then this
        function returns the order generated by the class of y.  If f
        is not monic, then :meth:`_make_monic_integral` is called, and instead we
        get the order generated by some integral multiple of a root of f.

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

           - ``base_morphism`` -- (default: None) a homomorphism of
             the base ring, after the im_gens are used.  Thus if
             im_gens has length 2, then base_morphism should be a morphism
             from self.base_ring().base_ring().

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

        from maps import FunctionFieldMorphism_polymod
        return FunctionFieldMorphism_polymod(self.Hom(codomain), im_gens[0], base_morphism)

    @cached_method
    def genus(self):
        """
        Return the genus of this function field
        For now, the genus is computed using singular

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
            raise NotImplementedError("Computation of genus over this rational function field not implemented yet")

    @cached_method
    def derivation(self):
        r"""
        Return a generator of the space of derivations over the constant base
        ring of this function field.

        A derivation on `R` is map `R\to R` with
        `D(\alpha+\beta)=D(\alpha)+D(\beta)` and `D(\alpha\beta)=\beta
        D(\alpha)+\alpha D(\beta)` for all `\alpha,\beta\in R`. For a function
        field which is a finite extension of `K(x)` with `K` perfect, the
        derivations form a one-dimensional `K`-vector space generated by the
        derivation returned by this method.

        ALGORITHM:

        If this field is a separable extension of another function field `F`,
        then Proposition 11 of [GT1996]_ describes how to compute the unique
        extension of a derivation on `F` to this field; we apply this algorithm
        to the generator of the space of derivations on `F`.
        If this field has not been generated as a separable extension, then we
        find an isomorphic field which is a separable extension of a rational
        function field, see :meth:`separable_model`.

        OUTPUT:

        An endofunction on this function field.

        REFERENCES:

        .. [GT1996] Gianni, P., & Trager, B. (1996). Square-free algorithms in
           positive characteristic. Applicable Algebra in Engineering,
           Communication and Computing, 7(1), 1-14.

        EXAMPLES::

            sage: K.<x> = FunctionField(GF(3))
            sage: R.<y> = K[]
            sage: L.<y> = K.extension(y^2 - x)
            sage: d = L.derivation(); d
            Derivation map:
                From: Function field in y defined by y^2 + 2*x
                To:   Function field in y defined by y^2 + 2*x
                Defn: y |--> 2/x*y
                      x |--> 1
            sage: d(x)
            1
            sage: d(x^3)
            0
            sage: d(x*y)
            0
            sage: d(y)
            2/x*y

        Currently the functionality for finding a separable model is not
        implemented (see :trac:`16562`, :trac:`16564`)::

            sage: R.<y> = K[]
            sage: L.<y> = K.extension(y^3 - x)
            sage: d = L.derivation(); d
            Derivation map:
                From: Function field in y defined by y^3 + 2*x
                To:   Function field in y defined by y^3 + 2*x
                Defn: y |--> 1
                      x |--> 0
            sage: d(x^2)
            0
            sage: d(y^2)
            2*y
            sage: d(x*y)
            x

        """
        from maps import FunctionFieldDerivation_separable, FunctionFieldDerivation_inseparable
        if self.is_separable():
            return FunctionFieldDerivation_separable(self, self.base_ring().derivation())
        else:
            return FunctionFieldDerivation_inseparable(self)

    def _simple_model(self, name='v'):
        r"""
        Helper method for :meth:`simple_model` which, for a tower of extensions
        `M/L/K(x)` over a rational function field `K(x)` with `K` perfect,
        finds a finite extension `N/K(x)` isomorphic to `M`.

        INPUT:

        - ``name`` -- a string, the name of the generator of `N`

        ALGORITHM:

        Since `K` is perfect, the extension `M/K(x)` is simple, i.e., generated
        by a single element [BM1940]_. Therefore, there are only finitely many
        intermediate fields (Exercise 3.6.7 in [Bosch2009]_).
        Let `a` be a generator of `M/L` and let `a` be a generator of `L/K(x)`.
        For some `i` the field `N_i=K(x)(a+x^ib)` is isomorphic to `M` and so
        it is enough to test for all terms of the form `a+x^ib` whether they
        generate a field of the right degree.
        Indeed, suppose for contradiction that for all `i` we had `N_i\neq M`.
        Then `N_i=N_j` for some `i,j`.  Thus `(a+x^ib)-(a+x^jb)=b(x^i-x^j)\in
        N_j` and so `b\in N_j`.  Similarly,
        `a+x^ib-x^{i-j}(a+x^jb)=a(1+x^{i-j})\in N_j` and so `a\in N_j`.
        Therefore, `N_j=M`.

        REFERENCES:

        .. [BM1940] Becker, M. F., and Saunders MacLane. The minimum number of
        generators for inseparable algebraic extensions. Bulletin of the
        American Mathematical Society 46, no. 2 (1940): 182-186.

        .. [Bosch2009] Bosch, S., Algebra, Springer 2009

        TESTS::

            sage: K.<x> = FunctionField(QQ)
            sage: R.<y> = K[]
            sage: L.<y> = K.extension(y^2-x)
            sage: R.<z> = L[]
            sage: M.<z> = L.extension(z^2-y)
            sage: M._simple_model()
            (Function field in v defined by x^4 - x,
             Function Field morphism:
              From: Function field in v defined by x^4 - x
              To:   Function field in z defined by z^2 - y
              Defn: v |--> z,
             Function Field morphism:
              From: Function field in z defined by z^2 - y
              To:   Function field in v defined by x^4 - x
              Defn: z |--> v
                    y |--> v^2)

        Check that this also works for inseparable extensions::

            sage: K.<x> = FunctionField(GF(2))
            sage: R.<y> = K[]
            sage: L.<y> = K.extension(y^2-x)
            sage: R.<z> = L[]
            sage: M.<z> = L.extension(z^2-y)
            sage: M._simple_model()
            (Function field in v defined by x^4 + x,
             Function Field morphism:
              From: Function field in v defined by x^4 + x
              To:   Function field in z defined by z^2 + y
              Defn: v |--> z,
             Function Field morphism:
              From: Function field in z defined by z^2 + y
              To:   Function field in v defined by x^4 + x
              Defn: z |--> v
                    y |--> v^2)

        An example where the generator of the last extension does not generate
        the extension of the rational function field::

            sage: K.<x> = FunctionField(GF(2))
            sage: R.<y> = K[]
            sage: L.<y> = K.extension(y^2-x)
            sage: R.<z> = L[]
            sage: M.<z> = L.extension(z^3-1)
            sage: M._simple_model()
            (Function field in v defined by x^6 + x*x^4 + x^2*x^2 + x^3 + 1,
             Function Field morphism:
               From: Function field in v defined by x^6 + x*x^4 + x^2*x^2 + x^3 + 1
               To:   Function field in z defined by z^3 + 1
               Defn: v |--> z + y,
             Function Field morphism:
               From: Function field in z defined by z^3 + 1
               To:   Function field in v defined by x^6 + x*x^4 + x^2*x^2 + x^3 + 1
               Defn: z |--> v^4 + x^2
                     y |--> v^4 + v + x^2)

        """
        M = self
        L = M.base_field()
        K = L.base_field()

        assert(is_RationalFunctionField(K))
        assert(K is not L)
        assert(L is not M)

        if not K.constant_field().is_perfect():
            raise NotImplementedError("simple_model() only implemented over perfect constant fields")

        x = K.gen()
        b = L.gen()
        a = M.gen()

        # using a+x^i*b tends to lead to huge powers of x in the minimal
        # polynomial of the resulting field; it is better to try terms of
        # the form a+i*b first (but in characteristic p>0 there are only
        # finitely many of these)
        # We systematically try elements of the form a+b*factor*x^exponent
        factor = self.constant_base_field().zero()
        exponent = 0
        while True:
            v = M(a+b*factor*x**exponent)
            minpoly = v.matrix(K).minpoly()
            if minpoly.degree() == M.degree()*L.degree():
                break
            factor += 1
            if factor == 0:
                factor = self.constant_base_field().one()
                exponent += 1

        N = K.extension(minpoly, names=(name,))

        # the morphism N -> M, v |-> v
        N_to_M = N.hom(v)

        # the morphism M -> N, b |-> M_b, a |-> M_a
        V, V_to_M, M_to_V = M.vector_space(K)
        V, V_to_N, N_no_V = N.vector_space(K)
        from sage.matrix.matrix_space import MatrixSpace
        MS = MatrixSpace(V.base_field(), V.dimension())
        # the power basis of v over K
        B = [M_to_V(v**i) for i in range(V.dimension())]
        B = MS(B)
        M_b = V_to_N(B.solve_left(M_to_V(b)))
        M_a = V_to_N(B.solve_left(M_to_V(a)))
        M_to_N = M.hom([M_a,M_b])

        return N, N_to_M, M_to_N

    @cached_method
    def simple_model(self, name=None):
        """
        Return a function field isomorphic to this field which is a simple
        extension of a rational function field.

        INPUT:

        - ``name`` -- a string or ``None`` (default: ``None``), the name of generator of the
          simple extension. If ``None``, then the name of the generator will be
          the same as the name of the generator of this function field.

        OUTPUT:

        A triple ``(F,f,t)`` where ``F`` is a field isomorphic to this field,
        ``f`` is an isomorphism from ``F`` to this function field and ``t`` is
        the inverse of ``f``.

        EXAMPLES:

        A tower of four function fields::

            sage: K.<x> = FunctionField(QQ); R.<z> = K[]
            sage: L.<z> = K.extension(z^2-x); R.<u> = L[]
            sage: M.<u> = L.extension(u^2-z); R.<v> = M[]
            sage: N.<v> = M.extension(v^2-u)

        The fields N and M as simple extensions of K::

            sage: N.simple_model()
            (Function field in v defined by x^8 - x,
             Function Field morphism:
              From: Function field in v defined by x^8 - x
              To:   Function field in v defined by v^2 - u
              Defn: v |--> v,
             Function Field morphism:
              From: Function field in v defined by v^2 - u
              To:   Function field in v defined by x^8 - x
              Defn: v |--> v
                    u |--> v^2
                    z |--> v^4
                    x |--> x)
            sage: M.simple_model()
            (Function field in u defined by x^4 - x,
             Function Field morphism:
              From: Function field in u defined by x^4 - x
              To:   Function field in u defined by u^2 - z
              Defn: u |--> u,
             Function Field morphism:
              From: Function field in u defined by u^2 - z
              To:   Function field in u defined by x^4 - x
              Defn: u |--> u
                    z |--> u^2
                    x |--> x)

        An optional parameter ``name`` can be used to set the name of the
        generator of the simple extension::

            sage: M.simple_model(name='t')
            (Function field in t defined by x^4 - x, Function Field morphism:
              From: Function field in t defined by x^4 - x
              To:   Function field in u defined by u^2 - z
              Defn: t |--> u, Function Field morphism:
              From: Function field in u defined by u^2 - z
              To:   Function field in t defined by x^4 - x
              Defn: u |--> t
                    z |--> t^2
                    x |--> x)

        An example with higher degrees::

            sage: K.<x> = FunctionField(GF(3)); R.<y> = K[]
            sage: L.<y> = K.extension(y^5-x); R.<z> = L[]
            sage: M.<z> = L.extension(z^3-x)
            sage: M.simple_model()
            (Function field in z defined by x^15 + x*x^12 + x^2*x^9 + 2*x^3*x^6 + 2*x^4*x^3 + 2*x^5 + 2*x^3,
             Function Field morphism:
               From: Function field in z defined by x^15 + x*x^12 + x^2*x^9 + 2*x^3*x^6 + 2*x^4*x^3 + 2*x^5 + 2*x^3
               To:   Function field in z defined by z^3 + 2*x
               Defn: z |--> z + y,
             Function Field morphism:
               From: Function field in z defined by z^3 + 2*x
               To:   Function field in z defined by x^15 + x*x^12 + x^2*x^9 + 2*x^3*x^6 + 2*x^4*x^3 + 2*x^5 + 2*x^3
               Defn: z |--> 2/x*z^6 + 2*z^3 + z + 2*x
                     y |--> 1/x*z^6 + z^3 + x
                     x |--> x)

        This also works for inseparable extensions::

            sage: K.<x> = FunctionField(GF(2)); R.<y> = K[]
            sage: L.<y> = K.extension(y^2-x); R.<z> = L[]
            sage: M.<z> = L.extension(z^2-y)
            sage: M.simple_model()
            (Function field in z defined by x^4 + x, Function Field morphism:
               From: Function field in z defined by x^4 + x
               To:   Function field in z defined by z^2 + y
               Defn: z |--> z, Function Field morphism:
               From: Function field in z defined by z^2 + y
               To:   Function field in z defined by x^4 + x
               Defn: z |--> z
                     y |--> z^2
                     x |--> x)

        """
        if name is None:
            name = self.variable_name()

        if is_RationalFunctionField(self.base_field()):
            # the extension is simple already
            if name == self.variable_name():
                from sage.categories.homset import Hom
                id = Hom(self,self).identity()
                return self, id, id
            else:
                ret = self.base_field().extension(self.polynomial(), names=(name,))
                f = ret.hom(self.gen())
                t = self.hom(ret.gen())
                return ret, f, t
        else:
            # recursively collapse the tower of fields
            base = self.base_field()
            base_, from_base_, to_base_ = base.simple_model()
            self_ = base_.extension(self.polynomial().map_coefficients(to_base_), names=(name,))
            gens_in_base_ = [to_base_(k.gen()) for k in base._intermediate_fields(base.rational_function_field())]
            to_self_ = self.hom([self_.gen()]+gens_in_base_)
            from_self_ = self_.hom([self.gen(),from_base_(base_.gen())])

            # now collapse self_/base_/K(x)
            ret, ret_to_self_, self__to_ret = self_._simple_model(name)
            ret_to_self = ret.hom(from_self_(ret_to_self_(ret.gen())))
            gens_in_ret = [self__to_ret(to_self_(k.gen())) for k in self._intermediate_fields(self.rational_function_field())]
            self_to_ret = self.hom(gens_in_ret)
            return ret, ret_to_self, self_to_ret

    @cached_method
    def primitive_element(self):
        r"""
        Return a primitive element over the underlying rational function field.

        If this is a finite extension of a rational function field `K(x)` with
        `K` perfect, then this is a simple extension of `K(x)`, i.e., there is
        a primitive element `y` which generates this field over `K(x)`. This
        method returns such an element `y`.

        EXAMPLES::

            sage: K.<x> = FunctionField(QQ)
            sage: R.<y> = K[]
            sage: L.<y> = K.extension(y^2-x)
            sage: R.<z> = L[]
            sage: M.<z> = L.extension(z^2-y)
            sage: R.<z> = L[]
            sage: N.<u> = L.extension(z^2-x-1)
            sage: N.primitive_element()
            u + y
            sage: M.primitive_element()
            z
            sage: L.primitive_element()
            y

        This also works for inseparable extensions::

            sage: K.<x> = FunctionField(GF(2))
            sage: R.<Y> = K[]
            sage: L.<y> = K.extension(Y^2-x)
            sage: R.<Z> = L[]
            sage: M.<z> = L.extension(Z^2-y)
            sage: M.primitive_element()
            z

        """
        N,f,t = self.simple_model()
        return f(N.gen())

    @cached_method
    def separable_model(self, names=None):
        r"""
        Return a function field isomorphic to this field which is a separable
        extension of a rational function field.

        INPUT:

        - ``names`` -- a tuple of two strings or ``None`` (default: ``None``);
          the second entry will be used as the variable name of the rational
          function field, the first entry will be used as the variable name of
          its separable extension. If ``None``, then the variable names will be
          chosen automatically.

        OUTPUT:

        A triple ``(F,f,t)`` where ``F`` is a function field, ``f`` is an
        isomorphism from ``F`` to this function field, and ``t`` is the inverse
        of ``f``.

        ALGORITHM:

        Suppose that the constant base field is perfect. If this is a monic
        integral inseparable extension of a rational function field, then the
        defining polynomial is separable if we swap the variables (Proposition
        4.1 in Chapter VIII of [Lang2002]_.)
        The algorithm reduces to this case with :meth:`monic_integral_model`.

        REFERENCES::

        .. [Lang2002] Serge Lang. Algebra. Springer, 2002.

        EXAMPLES::

            sage: K.<x> = FunctionField(GF(2))
            sage: R.<y> = K[]
            sage: L.<y> = K.extension(y^2 - x^3)
            sage: L.separable_model(('t','w'))
            (Function field in t defined by t^3 + w^2,
             Function Field morphism:
               From: Function field in t defined by t^3 + w^2
               To:   Function field in y defined by y^2 + x^3
               Defn: t |--> x
                     w |--> y,
             Function Field morphism:
               From: Function field in y defined by y^2 + x^3
               To:   Function field in t defined by t^3 + w^2
               Defn: y |--> w
                     x |--> t)

        This also works for non-integral polynomials::

            sage: K.<x> = FunctionField(GF(2))
            sage: R.<y> = K[]
            sage: L.<y> = K.extension(y^2/x - x^2)
            sage: L.separable_model()
            (Function field in y_ defined by y_^3 + x_^2,
             Function Field morphism:
               From: Function field in y_ defined by y_^3 + x_^2
               To:   Function field in y defined by 1/x*y^2 + x^2
               Defn: y_ |--> x
                     x_ |--> y,
             Function Field morphism:
               From: Function field in y defined by 1/x*y^2 + x^2
               To:   Function field in y_ defined by y_^3 + x_^2
               Defn: y |--> x_
                     x |--> y_)

        TESTS:

        Check that this also works in characteristic zero::

            sage: K.<x> = FunctionField(QQ)
            sage: R.<y> = K[]
            sage: L.<y> = K.extension(y^2 - x^3)
            sage: L.separable_model()
            (Function field in y defined by y^2 - x^3,
             Function Field endomorphism of Function field in y defined by y^2 - x^3
               Defn: y |--> y
                     x |--> x,
             Function Field endomorphism of Function field in y defined by y^2 - x^3
               Defn: y |--> y
                     x |--> x)

        Check that this works for towers of inseparable extensions::

            sage: K.<x> = FunctionField(GF(2))
            sage: R.<y> = K[]
            sage: L.<y> = K.extension(y^2 - x)
            sage: R.<z> = L[]
            sage: M.<z> = L.extension(z^2 - y)
            sage: M.separable_model()
            (Function field in z_ defined by z_ + x_^4,
             Function Field morphism:
               From: Function field in z_ defined by z_ + x_^4
               To:   Function field in z defined by z^2 + y
               Defn: z_ |--> x
                     x_ |--> z,
             Function Field morphism:
               From: Function field in z defined by z^2 + y
               To:   Function field in z_ defined by z_ + x_^4
               Defn: z |--> x_
                     y |--> x_^2
                     x |--> x_^4)

        Check that this also works if only the first extension is inseparable::

            sage: K.<x> = FunctionField(GF(2))
            sage: R.<y> = K[]
            sage: L.<y> = K.extension(y^2 - x)
            sage: R.<z> = L[]
            sage: M.<z> = L.extension(z^3 - y)
            sage: M.separable_model()
            (Function field in z_ defined by z_ + x_^6, Function Field morphism:
               From: Function field in z_ defined by z_ + x_^6
               To:   Function field in z defined by z^3 + y
               Defn: z_ |--> x
                     x_ |--> z, Function Field morphism:
               From: Function field in z defined by z^3 + y
               To:   Function field in z_ defined by z_ + x_^6
               Defn: z |--> x_
                     y |--> x_^3
                     x |--> x_^6)

        """
        if names is None:
            pass
        elif not isinstance(names, tuple):
            raise TypeError("names must be a tuple consisting of two strings")
        elif len(names) != 2:
            raise ValueError("must provide exactly two variable names")

        if self.base_ring() is not self.rational_function_field():
            L, from_L, to_L = self.simple_model()
            K, from_K, to_K = L.separable_model(names=names)
            f = K.hom([from_L(from_K(K.gen())), from_L(from_K(K.base_field().gen()))])
            t = self.hom([to_K(to_L(k.gen())) for k in self._intermediate_fields(self.rational_function_field())])
            return K, f, t

        if self.is_separable():
            # the model is already separable
            if names is None:
                names = self.variable_name(), self.base_field().variable_name()
            return self.change_variable_name(names)

        if not self.constant_base_field().is_perfect():
            raise NotImplementedError("constructing a separable model is only implemented for function fields over a perfect constant base field")

        if names is None:
            names = (self.variable_name()+"_", self.rational_function_field().variable_name()+"_")

        L, from_L, to_L = self.monic_integral_model()

        if L.is_separable():
            # L is separable
            ret, ret_to_L, L_to_ret = L.change_variable_name(names)
            f = ret.hom([from_L(ret_to_L(ret.gen())), from_L(ret_to_L(ret.base_field().gen()))])
            t = self.hom([L_to_ret(to_L(self.gen())), L_to_ret(to_L(self.base_field().gen()))])
            return ret, f, t
        else:
            # otherwise, the polynomial of L must be separable in the other variable
            from constructor import FunctionField
            K = FunctionField(self.constant_base_field(), names=(names[1],))
            # construct a field isomorphic to L on top of K

            # turn the minpoly of K into a bivariate polynomial
            if names[0] == names[1]:
                raise ValueError("names of generators must be distinct")
            from sage.rings.all import PolynomialRing
            R = PolynomialRing(self.constant_base_field(), names=names)
            S = R.remove_var(names[1])
            f = R( L.polynomial().change_variable_name(names[1]).map_coefficients(
                     lambda c:c.numerator().change_variable_name(names[0]), S))
            f = f.polynomial(R.gen(0)).change_ring(K)
            f /= f.leading_coefficient()
            # f must be separable in the other variable (otherwise it would factor)
            assert f.gcd(f.derivative()).is_one()

            ret = K.extension(f, names=(names[0],))
            # isomorphisms between L and ret are given by swapping generators
            ret_to_L = ret.hom( [L(L.base_field().gen()), L.gen()] )
            L_to_ret = L.hom( [ret(K.gen()), ret.gen()] )
            # compose with from_L and to_L to get the desired isomorphisms between self and ret
            f = ret.hom( [from_L(ret_to_L(ret.gen())), from_L(ret_to_L(ret.base_field().gen()))] )
            t = self.hom( [L_to_ret(to_L(self.gen())), L_to_ret(to_L(self.base_field().gen()))] )
            return ret, f, t

    def change_variable_name(self, name):
        r"""
        Return a field isomorphic to this field with variable(s) ``name``.

        INPUT:

        - ``name`` -- a string or a tuple consisting of a strings, the names of
          the new variables starting with a generator of this field and going
          down to the rational function field.

        OUTPUT:

        A triple ``F,f,t`` where ``F`` is a function field, ``f`` is an
        isomorphism from ``F`` to this field, and ``t`` is the inverse of
        ``f``.

        EXAMPLES::

            sage: K.<x> = FunctionField(QQ)
            sage: R.<y> = K[]
            sage: L.<y> = K.extension(y^2 - x)
            sage: R.<z> = L[]
            sage: M.<z> = L.extension(z^2 - y)

            sage: M.change_variable_name('zz')
            (Function field in zz defined by z^2 - y,
             Function Field morphism:
              From: Function field in zz defined by z^2 - y
              To:   Function field in z defined by z^2 - y
              Defn: zz |--> z
                    y |--> y
                    x |--> x,
             Function Field morphism:
              From: Function field in z defined by z^2 - y
              To:   Function field in zz defined by z^2 - y
              Defn: z |--> zz
                    y |--> y
                    x |--> x)
            sage: M.change_variable_name(('zz','yy'))
            (Function field in zz defined by z^2 - yy, Function Field morphism:
              From: Function field in zz defined by z^2 - yy
              To:   Function field in z defined by z^2 - y
              Defn: zz |--> z
                    yy |--> y
                    x |--> x, Function Field morphism:
              From: Function field in z defined by z^2 - y
              To:   Function field in zz defined by z^2 - yy
              Defn: z |--> zz
                    y |--> yy
                    x |--> x)
            sage: M.change_variable_name(('zz','yy','xx'))
            (Function field in zz defined by z^2 - yy,
             Function Field morphism:
              From: Function field in zz defined by z^2 - yy
              To:   Function field in z defined by z^2 - y
              Defn: zz |--> z
                    yy |--> y
                    xx |--> x,
             Function Field morphism:
              From: Function field in z defined by z^2 - y
              To:   Function field in zz defined by z^2 - yy
              Defn: z |--> zz
                    y |--> yy
                    x |--> xx)

        """
        if not isinstance(name, tuple):
            name = (name,)
        if len(name) == 0:
            raise ValueError("name must contain at least one string")
        elif len(name) == 1:
            base = self.base_field()
            from sage.categories.homset import Hom
            from_base = to_base = Hom(base,base).identity()
        else:
            base, from_base, to_base = self.base_field().change_variable_name(name[1:])

        ret = base.extension(self.polynomial().map_coefficients(to_base), names=(name[0],))
        f = ret.hom( [k.gen() for k in self._intermediate_fields(self.rational_function_field())] )
        t = self.hom( [k.gen() for k in ret._intermediate_fields(ret.rational_function_field())] )
        return ret, f, t

def is_RationalFunctionField(x):
    """
    Return ``True`` if ``x`` is of rational function field type.

    EXAMPLES::

        sage: from sage.rings.function_field.function_field import is_RationalFunctionField
        sage: is_RationalFunctionField(QQ)
        False
        sage: is_RationalFunctionField(FunctionField(QQ,'t'))
        True
    """
    return isinstance(x, RationalFunctionField)

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
        Maximal order in Rational function field in t over Finite Field of size 7

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
        Field.__init__(self, self, names=names, category = category)
        R = constant_field[names[0]]
        self._hash = hash((constant_field, names))
        self._constant_field = constant_field
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
        from constructor import FunctionField
        return FunctionField, (self._constant_field, self._names)

    def __hash__(self):
        """
        Return hash of this function field.

        The hash is formed from the constant field and the variable names.

        EXAMPLES::

            sage: K.<t> = FunctionField(QQ)
            sage: hash(K) == hash((K.constant_base_field(), K.variable_names()))
            True

        """
        return self._hash

    def _repr_(self):
        """
        Return string representation of this function field.

        EXAMPLES::

            sage: K.<t> = FunctionField(QQ)
            sage: K._repr_()
            'Rational function field in t over Rational Field'
        """
        return "Rational function field in %s over %s"%(
            self.variable_name(), self._constant_field)

    def _element_constructor_(self, x):
        r"""
        Coerce ``x`` into an element of this function field, possibly not canonically.

        INPUT:

            - ``x`` -- the element

        OUTPUT:

            ``x``, as an element of this function field

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
        from sage.rings.arith import LCM
        denom = LCM([a.denominator() for a in v])
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

    @cached_method
    def polynomial_ring(self, var='x'):
        """
        Return a polynomial ring in one variable over this rational function field.

        INPUT:

            - ``var`` -- a string (default: 'x')

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
        from maps import MapVectorSpaceToFunctionField, MapFunctionFieldToVectorSpace
        from_V = MapVectorSpaceToFunctionField(V, self)
        to_V   = MapFunctionFieldToVectorSpace(self, V)
        return (V, from_V, to_V)

    def random_element(self, *args, **kwds):
        """
        Create a random element of this rational function field.

        Parameters are passed to the random_element method of the
        underlying fraction field.

        EXAMPLES::

            sage: FunctionField(QQ,'alpha').random_element()   # random
            (-1/2*alpha^2 - 4)/(-12*alpha^2 + 1/2*alpha - 1/95)
        """
        return self(self._field.random_element(*args, **kwds))

    def degree(self, base=None):
        """
        Return the degree over the base field of this rational
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
        Return the ``n``-th generator of this function field.  If ``n`` is not
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
        Return the base field of this rational function field, which is just
        this function field itself.

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
        from maps import FunctionFieldMorphism_rational
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
    def maximal_order(self):
        """
        Return the maximal order of this function field.  Since this
        is a rational function field it is of the form K(t), and the
        maximal order is by definition K[t].

        EXAMPLES::

            sage: K.<t> = FunctionField(QQ)
            sage: K.maximal_order()
            Maximal order in Rational function field in t over Rational Field
            sage: K.equation_order()
            Maximal order in Rational function field in t over Rational Field
        """
        from function_field_order import FunctionFieldOrder_rational
        return FunctionFieldOrder_rational(self)

    equation_order = maximal_order

    def constant_base_field(self):
        """
        Return the field that this rational function field is a
        transcendental extension of.

        EXAMPLES::

            sage: K.<t> = FunctionField(QQ)
            sage: K.constant_field()
            Rational Field

        """
        return self._constant_field

    constant_field = constant_base_field

    def genus(self):
        """
        Return the genus of this function field
        This is always equal 0 for a rational function field

        EXAMPLES::

            sage: K.<x> = FunctionField(QQ);
            sage: K.genus()
            0
        """
        return 0

    @cached_method(key=lambda self, base: None)
    def vector_space(self, base=None):
        """
        Return a vector space `V` and isomorphisms from this field to `V` and
        from `V` to this field.

        This function allows us to identify the elements of this field with
        elements of a one-dimensional vector space over the field itself. This
        method exists so that all function fields (rational or not) have the
        same interface.

        INPUT:

        - ``base`` -- must be this field or ``None`` (default: ``None``); this
          parameter is ignored and merely exists to have the same interface as
          :meth:`FunctionField_polymod.vector_space`.

        OUTPUT:

        - ``V`` -- a vector space over base field
        - ``from_V`` -- an isomorphism from V to this field
        - ``to_V`` -- an isomorphism from this field to V

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
        from maps import MapVectorSpaceToFunctionField, MapFunctionFieldToVectorSpace
        if base is None:
            base = self
        if base is not self:
            raise ValueError("base must be the rational function field or None")
        V = base**1
        from_V = MapVectorSpaceToFunctionField(V, self)
        to_V   = MapFunctionFieldToVectorSpace(self, V)
        return (V, from_V, to_V)

    def change_variable_name(self, name):
        r"""
        Return a field isomorphic to this field with variable ``name``.

        INPUT:

        - ``name`` -- a string or a tuple consisting of a single string, the
          name of the new variable

        OUTPUT:

        A triple ``F,f,t`` where ``F`` is a rational function field, ``f`` is
        an isomorphism from ``F`` to this field, and ``t`` is the inverse of
        ``f``.

        EXAMPLES::

            sage: K.<x> = FunctionField(QQ)
            sage: L,f,t = K.change_variable_name('y')
            sage: L,f,t
            (Rational function field in y over Rational Field,
             Function Field morphism:
              From: Rational function field in y over Rational Field
              To:   Rational function field in x over Rational Field
              Defn: y |--> x,
             Function Field morphism:
              From: Rational function field in x over Rational Field
              To:   Rational function field in y over Rational Field
              Defn: x |--> y)
            sage: L.change_variable_name('x')[0] is K
            True

        """
        if isinstance(name, tuple):
            if len(name) != 1:
                raise ValueError("names must be a tuple with a single string")
            name = name[0]
        if name == self.variable_name():
            from sage.categories.homset import Hom
            id = Hom(self,self).identity()
            return self,id,id
        else:
            from constructor import FunctionField
            ret = FunctionField(self.constant_base_field(), name)
            return ret, ret.hom(self.gen()), self.hom(ret.gen())

    @cached_method
    def derivation(self):
        r"""
        Return a generator of the space of derivations over the constant base
        field of this function field.

        A derivation on `R` is a map `R \to R` with
        `D(\alpha + \beta) = D(\alpha) + D(\beta)` and
        `D(\alpha \beta) = \beta D(\alpha)+\alpha D(\beta)`
        for all `\alpha, \beta \in R`. For a function
        field `K(x)` with `K` perfect, the derivations form a one-dimensional
        `K`-vector space generated by the extension of the usual derivation on
        `K[x]` (cf. Proposition 10 in [GT1996]_.)

        OUTPUT:

        An endofunction on this function field.

        EXAMPLES::

            sage: K.<x> = FunctionField(GF(3))
            sage: K.derivation()
            Derivation map:
              From: Rational function field in x over Finite Field of size 3
              To:   Rational function field in x over Finite Field of size 3

        TESTS::

            sage: L.<y> = FunctionField(K)
            sage: L.derivation()
            Traceback (most recent call last):
            ...
            NotImplementedError: not implemented for non-perfect base fields

        """
        from maps import FunctionFieldDerivation_rational
        if not self.constant_base_field().is_perfect():
            raise NotImplementedError("not implemented for non-perfect base fields")
        return FunctionFieldDerivation_rational(self, self.one())
