r"""
Orders

An order of a function field is a subring that is, as a module over the base
maximal order, finitely generated and of maximal rank `n`, where `n` is the
extension degree of the function field. All orders are subrings of maximal
orders.

A rational function field has two maximal orders: maximal finite order `o` and
maximal infinite order `o_\infty`. The maximal order of a rational function
field over constant field `k` is just the polynomial ring `o=k[x]`. The
maximal infinite order is the set of rational functions whose denominator has
degree greater than or equal to that of the numerator.

EXAMPLES::

    sage: K.<x> = FunctionField(QQ)
    sage: O = K.maximal_order()
    sage: 1/x in O
    False
    sage: Oinf = K.maximal_order_infinite()
    sage: 1/x in Oinf
    True

In an extension of a rational function field, an order over the maximal finite
order is called a finite order while an order over the maximal infinite order
is called an infinite order. Thus a function field has one maximal finite order
`O` and one maximal infinite order `O_\infty`. There are other non-maximal
orders such as equation orders::

    sage: K.<x> = FunctionField(GF(3)); R.<y> = K[]
    sage: L.<y> = K.extension(y^3-y-x)
    sage: O = L.equation_order()
    sage: 1/y in O
    False
    sage: x/y in O
    True

AUTHORS:

- William Stein (2010): initial version

- Maarten Derickx (2011-09-14): fixed ideal_with_gens_over_base() for rational function fields

- Julian Rueth (2011-09-14): added check in _element_constructor_

- Kwankyu Lee (2017-04-30): added maximal orders of global function fields

"""
from __future__ import absolute_import
#*****************************************************************************
#       Copyright (C) 2010 William Stein <wstein@gmail.com>
#       Copyright (C) 2011 Maarten Derickx <m.derickx.student@gmail.com>
#       Copyright (C) 2011 Julian Rueth <julian.rueth@gmail.com>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#  as published by the Free Software Foundation; either version 2 of
#  the License, or (at your option) any later version.
#                  http://www.gnu.org/licenses/
#*****************************************************************************

from sage.misc.cachefunc import cached_method
from sage.misc.lazy_import import lazy_import

from sage.rings.ring import IntegralDomain, PrincipalIdealDomain

from sage.modules.free_module_element import vector
from sage.arith.all import lcm, gcd

from sage.structure.parent import Parent
from sage.structure.unique_representation import UniqueRepresentation

from sage.categories.rings import Rings
from sage.categories.integral_domains import IntegralDomains
from sage.categories.principal_ideal_domains import PrincipalIdealDomains
from sage.categories.euclidean_domains import EuclideanDomains

lazy_import('sage.matrix.special', 'block_matrix')
lazy_import('sage.matrix.constructor', 'matrix')

class FunctionFieldOrder(Parent):
    """
    Base class for orders in function fields.
    """
    def __init__(self, field, category=None):
        """
        Initialize.

        INPUT:

        - ``field`` -- function field

        EXAMPLES::

            sage: F = FunctionField(QQ,'y')
            sage: F.maximal_order()
            Maximal order of Rational function field in y over Rational Field
        """
        Parent.__init__(self, category=category or IntegralDomains())
        self._field = field

    def _repr_(self):
        """
        Return the string representation.

        EXAMPLES::

            sage: FunctionField(QQ,'y').maximal_order()._repr_()
            'Maximal order of Rational function field in y over Rational Field'
        """
        return "Order in {}".format(self._field)

    def is_finite(self):
        """
        Return ``False`` since orders are never finite.

        EXAMPLES::

            sage: FunctionField(QQ,'y').maximal_order().is_finite()
            False
        """
        return False

    def is_field(self):
        """
        Return ``False`` since orders are never fields.

        EXAMPLES::

            sage: FunctionField(QQ,'y').maximal_order().is_field()
            False
        """
        return False

    def is_noetherian(self):
        """
        Return ``True`` since orders in function fields are noetherian.

        EXAMPLES::

            sage: FunctionField(QQ,'y').maximal_order().is_noetherian()
            True
        """
        return True

    def function_field(self):
        """
        Return the function field to which the order belongs.

        EXAMPLES::

            sage: FunctionField(QQ,'y').maximal_order().function_field()
            Rational function field in y over Rational Field
        """
        return self._field

    fraction_field = function_field

    def is_subring(self, other):
        """
        Return ``True`` if the order is a subring of the other order.

        INPUT:

        - ``other`` -- order of the function field or the field itself

        EXAMPLES::

            sage: F = FunctionField(QQ,'y')
            sage: O = F.maximal_order()
            sage: O.is_subring(F)
            True
        """
        if other is self._field:
            return True
        else:
            raise NotImplementedError

class FunctionFieldOrder_basis(FunctionFieldOrder):
    """
    Order given by a basis over the maximal order of the base field.
    """
    def __init__(self, basis, check=True):
        """
        Initialize.

        INPUT:

        - ``basis`` -- list of elements of the function field

        - ``check`` -- (default: ``True``) if ``True``, check whether the module
          that ``basis`` generates forms an order

        EXAMPLES::

            sage: K.<x> = FunctionField(GF(7)); R.<y> = K[]
            sage: L.<y> = K.extension(y^4 + x*y + 4*x + 1)
            sage: O = L.equation_order(); O
            Order in Function field in y defined by y^4 + x*y + 4*x + 1

        The basis only defines an order if the module it generates is closed under
        multiplication and contains the identity element::

            sage: K.<x> = FunctionField(QQ)
            sage: R.<y> = K[]
            sage: L.<y> = K.extension(y^5 - (x^3 + 2*x*y + 1/x));
            sage: y.is_integral()
            False
            sage: L.order(y)
            Traceback (most recent call last):
            ...
            ValueError: The module generated by basis [1, y, y^2, y^3, y^4] must be closed under multiplication

        The basis also has to be linearly independent and of the same rank as the
        degree of the function field of its elements (only checked when ``check``
        is ``True``)::

            sage: L.order(L(x))
            Traceback (most recent call last):
            ...
            ValueError: Basis [1, x, x^2, x^3, x^4] is not linearly independent
            sage: sage.rings.function_field.order.FunctionFieldOrder_basis([y,y,y^3,y^4,y^5])
            Traceback (most recent call last):
            ...
            ValueError: Basis [y, y, y^3, y^4, 2*x*y + (x^4 + 1)/x] is not linearly independent
        """
        if len(basis) == 0:
            raise ValueError("basis must have positive length")

        field = basis[0].parent()
        if len(basis) != field.degree():
            raise ValueError("length of basis must equal degree of field")

        FunctionFieldOrder.__init__(self, field)

        self._basis = tuple(basis)
        V, fr, to = field.vector_space()
        R = field.base_field().maximal_order()
        self._from_module= fr
        self._to_module = to
        self._module = V.span([to(b) for b in basis], base_ring=R)
        self._ring = field.polynomial_ring()
        self._populate_coercion_lists_(coerce_list=[self._ring])
        if check:
            if self._module.rank() != field.degree():
                raise ValueError("Basis %s is not linearly independent"%(basis))
            if not to(field(1)) in self._module:
                raise ValueError("The identity element must be in the module spanned by basis %s"%(basis))
            if not all(to(a*b) in self._module for a in basis for b in basis):
                raise ValueError("The module generated by basis %s must be closed under multiplication"%(basis))

    def _element_constructor_(self, f, check=True):
        """
        Make ``f`` an element of the order.

        INPUT:

        - ``f`` -- the element

        - ``check`` -- check if the element is in the order

        EXAMPLES::

            sage: K.<x> = FunctionField(QQ)
            sage: K.maximal_order()._element_constructor_(x)
            x
        """
        field = self.function_field()

        if f.parent() is field:
            f = f.element()
        f = self._ring(f)
        if check:
            V, fr, to = field.vector_space()
            f_vector = to(field(f))
            if not f_vector in self._module:
                raise TypeError("%r is not an element of %r"%(f_vector,self))
        return field._element_class(self, f)

    def polynomial(self):
        """
        Return the defining polynomial of the function field of which this is an order.

        EXAMPLES::

            sage: K.<x> = FunctionField(GF(7)); R.<y> = K[]
            sage: L.<y> = K.extension(y^4 + x*y + 4*x + 1)
            sage: O = L.equation_order()
            sage: O.polynomial()
            y^4 + x*y + 4*x + 1
        """
        return self._field.polynomial()

    def basis(self):
        """
        Return a basis of the order over the maximal order of the base field.

        EXAMPLES::

            sage: K.<x> = FunctionField(GF(7)); R.<y> = K[]
            sage: L.<y> = K.extension(y^4 + x*y + 4*x + 1)
            sage: O = L.equation_order()
            sage: O.basis()
            (1, y, y^2, y^3)
        """
        return self._basis

    def free_module(self):
        """
        Return the free module formed by the basis over the maximal order
        of the base function field.

        EXAMPLES::

            sage: K.<x> = FunctionField(GF(7)); R.<y> = K[]
            sage: L.<y> = K.extension(y^4 + x*y + 4*x + 1)
            sage: O = L.equation_order()
            sage: O.free_module()
            Free module of degree 4 and rank 4 over Maximal order of Rational function field in x over Finite Field of size 7
            Echelon basis matrix:
            [1 0 0 0]
            [0 1 0 0]
            [0 0 1 0]
            [0 0 0 1]
        """
        return self._module

    def coordinate_vector(self, e):
        """
        Return the cooridinates of ``e`` with respect to the basis of the order.

        INPUT:

        - ``e`` -- element of the order or the function field

        EXAMPLES::

            sage: K.<x> = FunctionField(GF(7)); R.<y> = K[]
            sage: L.<y> = K.extension(y^4 + x*y + 4*x + 1)
            sage: O = L.equation_order()
            sage: f = (x + y)^3
            sage: O.coordinate_vector(f)
            (x^3, 3*x^2, 3*x, 1)
        """
        return self._module.coordinate_vector(self._to_module(e), check=False)

class FunctionFieldOrderInfinite(FunctionFieldOrder):
    """
    Base class for infinite orders on function fields.
    """
    def _repr_(self):
        """
        EXAMPLES::

            sage: FunctionField(QQ,'y').maximal_order()._repr_()
            'Maximal order of Rational function field in y over Rational Field'
        """
        return "Infinite order in %s"%self.function_field()

class FunctionFieldOrderInfinite_basis(FunctionFieldOrderInfinite):
    """
    Order given by a basis over the infinite maximal order of the base
    field.

    EXAMPLES::

        sage: K.<x> = FunctionField(GF(7)); R.<y> = K[]
        sage: L.<y> = K.extension(y^4 + x*y + 4*x + 1)
        sage: O = L.equation_order(); O
        Order in Function field in y defined by y^4 + x*y + 4*x + 1
        sage: type(O)
        <class 'sage.rings.function_field.order.FunctionFieldOrder_basis_with_category'>

    The basis only defines an order if the module it generates is closed under
    multiplication and contains the identity element (only checked when
    ``check`` is ``True``)::

        sage: K.<x> = FunctionField(QQ)
        sage: R.<y> = K[]
        sage: L.<y> = K.extension(y^5 - (x^3 + 2*x*y + 1/x));
        sage: y.is_integral()
        False
        sage: L.order(y)
        Traceback (most recent call last):
        ...
        ValueError: The module generated by basis [1, y, y^2, y^3, y^4] must be closed under multiplication

    The basis also has to be linearly independent and of the same rank as the
    degree of the function field of its elements (only checked when ``check``
    is ``True``)::

        sage: L.order(L(x))
        Traceback (most recent call last):
        ...
        ValueError: Basis [1, x, x^2, x^3, x^4] is not linearly independent
        sage: sage.rings.function_field.order.FunctionFieldOrder_basis([y,y,y^3,y^4,y^5])
        Traceback (most recent call last):
        ...
        ValueError: Basis [y, y, y^3, y^4, 2*x*y + (x^4 + 1)/x] is not linearly independent
    """
    def __init__(self, basis, check=True):
        """
        Initialize.

        INPUT:

        - ``basis`` -- elements of the function field

        - ``check`` -- boolean (default: ``True``); if ``True``, check the basis generates
          an order

        EXAMPLES::

            sage: K.<x> = FunctionField(GF(7)); R.<y> = K[]
            sage: L.<y> = K.extension(y^4 + x*y + 4*x + 1)
            sage: L.equation_order()
            Order in Function field in y defined by y^4 + x*y + 4*x + 1
        """
        if len(basis) == 0:
            raise ValueError("basis must have positive length")

        field = basis[0].parent()
        if len(basis) != field.degree():
            raise ValueError("length of basis must equal degree of field")

        FunctionFieldOrder.__init__(self, field)

        # The function field element f is in this order if and only if
        # W.coordinate_vector(to(f)) in M
        V, fr, to = field.vector_space()
        R = field.base_field().maximal_order_infinite()
        W = V.span_of_basis([to(v) for v in basis])
        from sage.modules.free_module import FreeModule
        M = FreeModule(R,W.dimension())
        self._basis = tuple(basis)
        self._ambient_space = W
        self._module = M

        self._ring = field.polynomial_ring()
        self._populate_coercion_lists_(coerce_list=[self._ring])

        if check:
            if self._module.rank() != field.degree():
                raise ValueError("Basis %s is not linearly independent"%(basis))
            if not W.coordinate_vector(to(field(1))) in self._module:
                raise ValueError("The identity element must be in the module spanned by basis %s"%(basis))
            if not all(W.coordinate_vector(to(a*b)) in self._module for a in basis for b in basis):
                raise ValueError("The module generated by basis %s must be closed under multiplication"%(basis))

    def _element_constructor_(self, f, check=True):
        """
        Make ``f`` an element of this order.

        INPUT:

        - ``f`` -- the element

        - ``check`` -- check if the element is in the order

        EXAMPLES::

            sage: K.<x> = FunctionField(QQ)
            sage: K.maximal_order()._element_constructor_(x)
            x
        """
        field=self.function_field()

        if not f.parent() is field:
            f = field(f)
        if check:
            V, fr, to = field.vector_space()
            W = self._ambient_space
            if not W.coordinate_vector(to(f)) in self._module:
                raise TypeError("%r is not an element of %r"%(f,self))
        return field._element_class(self, f.element())

    def function_field(self):
        """
        Return the function field in which this is an order.

        EXAMPLES::

            sage: K.<x> = FunctionField(GF(7)); R.<y> = K[]
            sage: L.<y> = K.extension(y^4 + x*y + 4*x + 1)
            sage: O = L.equation_order()
            sage: O.function_field()
            Function field in y defined by y^4 + x*y + 4*x + 1
        """
        return self._field

    fraction_field = function_field

    def polynomial(self):
        """
        Return the defining polynomial of the function field of which this is an order.

        EXAMPLES::

            sage: K.<x> = FunctionField(GF(7)); R.<y> = K[]
            sage: L.<y> = K.extension(y^4 + x*y + 4*x + 1)
            sage: O = L.equation_order()
            sage: O.polynomial()
            y^4 + x*y + 4*x + 1
        """
        return self._field.polynomial()

    def basis(self):
        """
        Return a basis of this order over the maximal order of the base field.

        EXAMPLES::

            sage: K.<x> = FunctionField(GF(7)); R.<y> = K[]
            sage: L.<y> = K.extension(y^4 + x*y + 4*x + 1)
            sage: O = L.equation_order()
            sage: O.basis()
            (1, y, y^2, y^3)
        """
        return self._basis

    def free_module(self):
        """
        Return the free module formed by the basis over the maximal order of
        the base field.

        EXAMPLES::

            sage: K.<x> = FunctionField(GF(7)); R.<y> = K[]
            sage: L.<y> = K.extension(y^4 + x*y + 4*x + 1)
            sage: O = L.equation_order()
            sage: O.free_module()
            Free module of degree 4 and rank 4 over Maximal order of Rational function field in x over Finite Field of size 7
            Echelon basis matrix:
            [1 0 0 0]
            [0 1 0 0]
            [0 0 1 0]
            [0 0 0 1]
        """
        return self._module

class FunctionFieldMaximalOrder(FunctionFieldOrder):
    """
    Base class of maximal orders of function fields.
    """
    def _repr_(self):
        """
        Return the string representation of the order.

        EXAMPLES::

            sage: FunctionField(QQ,'y').maximal_order()._repr_()
            'Maximal order of Rational function field in y over Rational Field'
        """
        return "Maximal order of %s"%(self.function_field(),)

class FunctionFieldMaximalOrder_rational(FunctionFieldMaximalOrder):
    """
    Maximal orders of rational function fields.
    """
    def __init__(self, field):
        """
        Initialize.

        INPUT:

        - ``field`` -- function field

        EXAMPLES::

            sage: K.<t> = FunctionField(GF(19)); K
            Rational function field in t over Finite Field of size 19
            sage: R = K.maximal_order(); R
            Maximal order of Rational function field in t over Finite Field of size 19
        """
        FunctionFieldMaximalOrder.__init__(self, field, category=EuclideanDomains())

        self._populate_coercion_lists_(coerce_list=[field._ring])

        self._ring = field._ring
        self._gen = self(self._ring.gen())
        self._basis = (self.one(),)

    def _element_constructor_(self, f):
        """
        Make ``f`` an element of this order.

        EXAMPLES::

            sage: K.<y> = FunctionField(QQ)
            sage: O = K.maximal_order()
            sage: O._element_constructor_(y)
            y
            sage: O._element_constructor_(1/y)
            Traceback (most recent call last):
            ...
            TypeError: 1/y is not an element of Maximal order of Rational function field in y over Rational Field
        """
        if f.parent() is self.function_field():
            if not f.denominator() in self.function_field().constant_base_field():
                raise TypeError("%r is not an element of %r"%(f,self))
            f = f.element()
        from .element import FunctionFieldElement_rational
        return FunctionFieldElement_rational(self, self._ring(f))

    def basis(self):
        """
        Return the basis (=1) of the order as a module over the polynomial ring.

        EXAMPLES::

            sage: K.<t> = FunctionField(GF(19))
            sage: O = K.maximal_order()
            sage: O.basis()
            (1,)
            sage: parent(O.basis()[0])
            Maximal order of Rational function field in t over Finite Field of size 19
        """
        return self._basis

    def gen(self, n=0):
        """
        Return the ``n``-th generator of the order. Since there is only one generator ``n`` must be 0.

        EXAMPLES::

            sage: O = FunctionField(QQ,'y').maximal_order()
            sage: O.gen()
            y
            sage: O.gen(1)
            Traceback (most recent call last):
            ...
            IndexError: Only one generator.
        """
        if n != 0: raise IndexError("Only one generator.")
        return self._gen

    def ngens(self):
        """
        Return 1 the number of generators of the order.

        EXAMPLES::

            sage: FunctionField(QQ,'y').maximal_order().ngens()
            1
        """
        return 1

class FunctionFieldMaximalOrder_global(FunctionFieldMaximalOrder):
    """
    Maximal orders of global function fields.
    """
    def __init__(self, field, basis):
        """
        Initialize.

        INPUT:

        - ``field`` -- function field to which this maximal order belongs

        - ``basis`` -- basis of this maximal order as a module over the base
          maximal order

        EXAMPLES::

            sage: K.<x> = FunctionField(GF(7)); R.<y> = K[]
            sage: L.<y> = K.extension(y^4 + x*y + 4*x + 1)
            sage: L.maximal_order()
            Maximal order of Function field in y defined by y^4 + x*y + 4*x + 1
        """
        FunctionFieldMaximalOrder.__init__(self, field)

        V, fr, to = field.vector_space()
        R = field.base_field().maximal_order()

        # This module is over R, but linear algebra over R (MaximalOrder)
        # is not well supported in Sage. So we keep it as a vector space
        # over rational function field.
        self._module = V.span_of_basis([to(b) for b in basis])
        self._module_base_ring = R
        self._basis = tuple(basis)
        self._from_module = fr
        self._to_module = to

        # multiplication table (lower triangular)
        n = len(basis)
        self._mtable = []
        for i in range(n):
            row = []
            for j in range(n):
                row.append(self._coordinate_vector(basis[i] * basis[j]))
            self._mtable.append(row)

        zero = vector(R._ring,n*[0])
        def mul_vecs(f,g):
            s = zero
            for i in range(n):
                if f[i].is_zero():
                    continue
                for j in range(n):
                    if g[j].is_zero():
                        continue
                    s += f[i] * g[j] * self._mtable[i][j]
            return s
        self._mul_vecs = mul_vecs

        # We prepare for using Kummer's theorem to decompose primes. Note
        # that Kummer's theorem applies to the most of places. Here we find
        # places for which the theorem does not apply.

        # this element is integral over k[x] and a generator of the field.
        for gen in basis[1:]:
            phi = gen.minimal_polynomial()
            if phi.degree() == n:
                break

        if phi.degree() == n:
            gen_vec = self._coordinate_vector(gen)
            g = gen_vec.parent().gen(0) # x
            gen_vec_pow = [g]
            for i in range(n):
                g = mul_vecs(g, gen_vec)
                gen_vec_pow.append(g)

            # find places where {1,gen,...,gen^(n-1)} is not integral basis
            W = V.span_of_basis([to(gen ** i) for i in range(phi.degree())])

            supp = []
            for g in basis:
                for c in W.coordinate_vector(to(g), check=False):
                    if not c.is_zero():
                        supp += [f for f,_ in c.denominator().factor()]
            supp = set(supp)

            self._kummer_gen = gen
            self._kummer_gen_vec_pow = gen_vec_pow
            self._kummer_polynomial = phi
            self._kummer_places = supp

    def _element_constructor_(self, f):
        """
        Construct an element of this order from ``f``.

        INPUT:

        - ``f`` -- element convertible to the function field

        EXAMPLES::

            sage: K.<x> = FunctionField(GF(4)); _.<Y> = K[]
            sage: L.<y> = K.extension(Y^2-x*Y+x^2+1)
            sage: O = L.maximal_order()
            sage: y in O
            True
            sage: 1/y in O
            False
            sage: x in O
            True
            sage: 1/x in O
            False
            sage: L.<y>=K.extension(Y^2+Y+x+1/x)
            sage: O = L.maximal_order()
            sage: 1 in O
            True
            sage: y in O
            False
            sage: x*y in O
            True
            sage: x^2*y in O
            True
        """
        F = self.function_field()

        if f.parent() is not F:
            f = F(f)

        # check if f is in this order
        if not all(e in self._module_base_ring for e in self.coordinate_vector(f)):
            raise TypeError( "{} is not an element of {}".format(f, self) )

        return f

    def polynomial(self):
        """
        Return the defining polynomial of the function field of which this is an order.

        EXAMPLES::

            sage: K.<x> = FunctionField(GF(7)); R.<y> = K[]
            sage: L.<y> = K.extension(y^4 + x*y + 4*x + 1)
            sage: O = L.equation_order()
            sage: O.polynomial()
            y^4 + x*y + 4*x + 1
        """
        return self._field.polynomial()

    def basis(self):
        """
        Return a basis of the order over the maximal order of the base function
        field.

        EXAMPLES::

            sage: K.<x> = FunctionField(GF(7)); R.<y> = K[]
            sage: L.<y> = K.extension(y^4 + x*y + 4*x + 1)
            sage: O = L.equation_order()
            sage: O.basis()
            (1, y, y^2, y^3)
        """
        return self._basis

    def free_module(self):
        """
        Return the free module formed by the basis over the maximal order of the base field.

        EXAMPLES::

            sage: K.<x> = FunctionField(GF(7)); R.<y> = K[]
            sage: L.<y> = K.extension(y^4 + x*y + 4*x + 1)
            sage: O = L.maximal_order()
            sage: O.free_module()
            Free module of degree 4 and rank 4 over Maximal order of Rational function field in x over Finite Field of size 7
            User basis matrix:
            [1 0 0 0]
            [0 1 0 0]
            [0 0 1 0]
            [0 0 0 1]
        """
        return self._module.change_ring(self._module_base_ring)

    def coordinate_vector(self, e):
        """
        Return the cooridinates of ``e`` with respect to the basis of this order.

        EXAMPLES::

            sage: K.<x> = FunctionField(GF(7)); R.<y> = K[]
            sage: L.<y> = K.extension(y^4 + x*y + 4*x + 1)
            sage: O = L.maximal_order()
            sage: O.coordinate_vector(y)
            (0, 1, 0, 0)
            sage: O.coordinate_vector(x*y)
            (0, x, 0, 0)
        """
        return self._module.coordinate_vector(self._to_module(e))

    def _coordinate_vector(self, e):
        """
        Return the cooridinate vector of ``e`` with respect to the basis
        of the order.

        Assuming ``e`` is in the maximal order, the coordinates are given
        as univariate polynomials in the underlying ring of the maximal
        order of the rational function field.

        EXAMPLES::

            sage: K.<x> = FunctionField(GF(7)); R.<y> = K[]
            sage: L.<y> = K.extension(y^4 + x*y + 4*x + 1)
            sage: O = L.maximal_order()
            sage: O._coordinate_vector(y)
            (0, 1, 0, 0)
            sage: O._coordinate_vector(x*y)
            (0, x, 0, 0)
        """
        v = self._module.coordinate_vector(self._to_module(e), check=False)
        return vector([c.numerator() for c in v])

class FunctionFieldMaximalOrderInfinite(FunctionFieldMaximalOrder, FunctionFieldOrderInfinite):
    """
    Base class of maximal infinite orders of function fields.
    """
    def _repr_(self):
        """
        EXAMPLES::

            sage: FunctionField(QQ,'y').maximal_order_infinite()
            Maximal infinite order of Rational function field in y over Rational Field

            sage: K.<x> = FunctionField(GF(2)); R.<t> = PolynomialRing(K)
            sage: F.<y> = K.extension(t^3-x^2*(x^2+x+1)^2)
            sage: F.maximal_order_infinite()
            Maximal infinite order of Function field in y defined by y^3 + x^6 + x^4 + x^2
        """
        return "Maximal infinite order of %s"%(self.function_field(),)

class FunctionFieldMaximalOrderInfinite_rational(FunctionFieldMaximalOrderInfinite):
    """
    Maximal infinite orders of rational function fields.
    """
    def __init__(self, field, category=None):
        """
        Initialize.

        INPUT:

        - ``field`` -- function field

        EXAMPLES::

            sage: K.<t> = FunctionField(GF(19)); K
            Rational function field in t over Finite Field of size 19
            sage: R = K.maximal_order_infinite(); R
            Maximal infinite order of Rational function field in t over Finite Field of size 19
        """
        FunctionFieldOrderInfinite.__init__(self, field,
                                            category=category or PrincipalIdealDomains())
        self._populate_coercion_lists_(coerce_list=[field.constant_base_field()])

    def basis(self):
        """
        Return the basis (=1) of the order as a module over the polynomial ring.

        EXAMPLES::

            sage: K.<t> = FunctionField(GF(19))
            sage: O = K.maximal_order()
            sage: O.basis()
            (1,)
            sage: parent(O.basis()[0])
            Maximal order of Rational function field in t over Finite Field of size 19
        """
        return 1/self.function_field().gen()

    def gen(self, n=0):
        """
        Return the ``n``-th generator of self. Since there is only one generator ``n`` must be 0.

        EXAMPLES::

            sage: O = FunctionField(QQ,'y').maximal_order()
            sage: O.gen()
            y
            sage: O.gen(1)
            Traceback (most recent call last):
            ...
            IndexError: Only one generator.
        """
        if n != 0: raise IndexError("Only one generator.")
        return self._gen

    def ngens(self):
        """
        Return 1 the number of generators of the order.

        EXAMPLES::

            sage: FunctionField(QQ,'y').maximal_order().ngens()
            1
        """
        return 1

    def _element_constructor_(self, f):
        """
        Make ``f`` an element of this order.

        EXAMPLES::

            sage: K.<y> = FunctionField(QQ)
            sage: O = K.maximal_order()
            sage: O._element_constructor_(y)
            y
            sage: O._element_constructor_(1/y)
            Traceback (most recent call last):
            ...
            TypeError: 1/y is not an element of Maximal order of Rational function field in y over Rational Field
        """
        if not f.parent() is self.function_field():
            f = self.function_field()(f)
        f = f.element()
        if f.denominator().degree() < f.numerator().degree():
            raise TypeError("%r is not an element of %r"%(f,self))
        from .element import FunctionFieldElement_rational
        return FunctionFieldElement_rational(self, f)

class FunctionFieldMaximalOrderInfinite_global(FunctionFieldMaximalOrderInfinite):
    """
    Maximal infinite orders of global function fields.
    """
    def __init__(self, field, basis):
        """
        Initialize.

        INPUT:

        - ``field`` -- function field

        - ``basis`` -- list of elements of the function field

        EXAMPLES::

            sage: K.<x> = FunctionField(GF(2)); _.<t> = PolynomialRing(K)
            sage: F.<y> = K.extension(t^3-x^2*(x^2+x+1)^2)
            sage: F.maximal_order_infinite()
            Maximal infinite order of Function field in y defined by y^3 + x^6 + x^4 + x^2

            sage: K.<x> = FunctionField(GF(2)); _.<Y> = K[]
            sage: L.<y> = K.extension(Y^2 + Y + x + 1/x)
            sage: L.maximal_order_infinite()
            Maximal infinite order of Function field in y defined by y^2 + y + (x^2 + 1)/x
        """
        FunctionFieldOrderInfinite.__init__(self, field)

        V, fr, to = field.vector_space()
        R = field.base_field().maximal_order_infinite()

        self._basis = tuple(basis)
        self._module = V.span_of_basis([to(v) for v in basis])
        self._module_base_ring = R
        self._to_module = to

    def _element_constructor_(self, f):
        """
        Make ``f`` an element of this order.

        EXAMPLES::

            sage: K.<x> = FunctionField(GF(2)); _.<Y> = K[]
            sage: L.<y> = K.extension(Y^2 + Y + x + 1/x)
            sage: Oinf = L.maximal_order_infinite()
            sage: Oinf.basis()
            (1, 1/x*y)
            sage: 1 in Oinf
            True
            sage: 1/x*y in Oinf
            True
            sage: x*y in Oinf
            False
            sage: 1/x in Oinf
            True
        """
        if not f.parent() is self.function_field():
            f = self.function_field()(f)

        oinf = self.function_field().base_field().maximal_order_infinite()
        coordinates = self.coordinate_vector(f)
        if not all(c in oinf for c in coordinates):
            raise TypeError("%r is not an element of %r"%(f,self))
        return f

    def basis(self):
        """
        Return a basis of this order as a module over the maximal order
        of the base function field.

        EXAMPLES::

            sage: K.<x> = FunctionField(GF(2)); _.<t> = K[]
            sage: L.<y> = K.extension(t^3 - x^2*(x^2 + x + 1)^2)
            sage: Oinf = L.maximal_order_infinite()
            sage: Oinf.basis()
            (1, 1/x^2*y, (1/(x^4 + x^3 + x^2))*y^2)

            sage: K.<x> = FunctionField(GF(2)); _.<Y> = K[]
            sage: L.<y> = K.extension(Y^2 + Y + x + 1/x)
            sage: Oinf = L.maximal_order_infinite()
            sage: Oinf.basis()
            (1, 1/x*y)
        """
        return self._basis

    def coordinate_vector(self, e):
        """
        Return the cooridinates of ``e`` with respect to the basis of the order.

        INPUT:

        - ``e`` -- element of the function field

        The returned coordinates are in the base maximal infinite order if and only
        if the element is in the order.

        EXAMPLES::

            sage: K.<x> = FunctionField(GF(2)); _.<Y> = K[]
            sage: L.<y> = K.extension(Y^2 + Y + x + 1/x)
            sage: Oinf = L.maximal_order_infinite()
            sage: f = 1/y^2
            sage: f in Oinf
            True
            sage: Oinf.coordinate_vector(f)
            ((x^3 + x^2 + x)/(x^4 + 1), x^3/(x^4 + 1))
        """
        return self._module.coordinate_vector(self._to_module(e))
