r"""
Coordinate Functions

In the context of a topological manifold `M` over a topological field `K`,
a *coordinate function*  is a function from a chart codomain
to `K`. In other words, a coordinate function is a `K`-valued function of
the coordinates associated to some chart.

More precisely, let `(U,\varphi)` be a chart on `M`, i.e. `U` is an open
subset of `M` and `\varphi: U \rightarrow V \subset K^n` is a homeomorphism
from `U` to an open subset `V` of `K^n`. A *coordinate function associated
to the chart* `(U,\varphi)` is a function

.. MATH::

    \begin{array}{cccc}
        f:&  V\subset K^n & \longrightarrow & K \\
          &  (x^1,\ldots, x^n) & \longmapsto & f(x^1,\ldots, x^n)
    \end{array}

Coordinate functions are implemented by derived classes of the abstract base
class :class:`CoordFunction`.

The class :class:`MultiCoordFunction` implements `K^m`-valued functions of the
coordinates of a chart, with $m$ a positive integer.

AUTHORS:

- Eric Gourgoulhon, Michal Bejger (2013-2015) : initial version

"""
#*****************************************************************************
#       Copyright (C) 2015 Eric Gourgoulhon <eric.gourgoulhon@obspm.fr>
#       Copyright (C) 2015 Michal Bejger <bejger@camk.edu.pl>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#  as published by the Free Software Foundation; either version 2 of
#  the License, or (at your option) any later version.
#                  http://www.gnu.org/licenses/
#*****************************************************************************

from sage.structure.sage_object import SageObject
from sage.misc.latex import latex

class CoordFunction(SageObject):
    r"""
    Abstract base class for coordinate functions.

    If `(U,\varphi)` is a chart on a topological manifold `M` of dimension `n`
    over a topological field `K`,  a *coordinate function* associated to
    `(U,\varphi)` is a map `f: V\subset K^n \rightarrow K`, where `V` is the
    codomain of `\varphi`. In other words, `f` is a `K`-valued function of the
    coordinates associated to the chart `(U,\varphi)`.

    The class :class:`CoordFunction` is an abstract one. Specific coordinate
    functions must be implemented by derived classes, like
    :class:`~sage.manifolds.coord_func_symb.CoordFunctionSymb` for
    symbolic coordinate functions.

    INPUT:

    - ``chart`` -- the chart `(U, \varphi)`, as an instance of class
      :class:`~sage.manifolds.chart.Chart`

    """
    def __init__(self, chart):
        r"""
        Base constructor for derived classes.

        TEST::

            sage: M = Manifold(2, 'M', structure='topological')
            sage: X.<x,y> = M.chart()
            sage: from sage.manifolds.coord_func import CoordFunction
            sage: f = CoordFunction(X)

        """
        self._chart = chart
        self._nc = len(chart[:])    # number of coordinates

    # ----------------------------------------------------------------
    # Methods that do not need to be re-implemented by derived classes
    # ----------------------------------------------------------------

    def chart(self):
        r"""
        Return the chart w.r.t. which the coordinate function is defined.

        OUTPUT:

        - an instance of :class:`~sage.manifolds.chart.Chart`

        EXAMPLE::

            sage: M = Manifold(2, 'M', structure='topological')
            sage: X.<x,y> = M.chart()
            sage: f = X.function(1+x+y^2)
            sage: f.chart()
            Chart (M, (x, y))
            sage: f.chart() is X
            True

        """
        return self._chart

    def scalar_field(self, name=None, latex_name=None):
        r"""
        Construct the scalar field that has the coordinate function as
        coordinate expression.

        The domain of the scalar field is the open subset covered by the chart
        on which the coordinate function is defined

        INPUT:

        - ``name`` -- (default: ``None``) name given to the scalar field
        - ``latex_name`` -- (default: ``None``) LaTeX symbol to denote the
          scalar field; if none is provided, the LaTeX symbol is set to ``name``

        OUTPUT:

        - instance of class
          :class:`~sage.manifolds.scalarfield.ScalarField`

        EXAMPLES:

        Construction of a scalar field on a 2-dimensional manifold::

            sage: M = Manifold(2, 'M', structure='topological')
            sage: c_xy.<x,y> = M.chart()
            sage: fc = c_xy.function(x+2*y^3)
            sage: f = fc.scalar_field() ; f
            Scalar field on the 2-dimensional topological manifold M
            sage: f.display()
            M --> R
            (x, y) |--> 2*y^3 + x
            sage: f.coord_function(c_xy) is fc
            True

        """
        alg = self._chart.domain().scalar_field_algebra()
        return alg.element_class(alg, coord_expression={self._chart: self},
                                 name=name, latex_name=latex_name)

    def __ne__(self, other):
        r"""
        Inequality operator.

        INPUT:

        - ``other`` -- another instance of :class:`CoordFunction` or a value

        OUTPUT:

        - ``True`` if ``self`` is different from ``other``,  or ``False``
          otherwise

        TESTS::

            sage: M = Manifold(2, 'M', structure='topological')
            sage: X.<x,y> = M.chart()
            sage: f = X.function(x+y)
            sage: g = X.function(x*y)
            sage: f != g
            True
            sage: h = X.function(x+y)
            sage: f != h
            False

        """
        return not (self == other)

    def __radd__(self, other):
        r"""
        Reflected addition operator: performs ``other + self``.

        INPUT:

        - ``other`` -- another instance of :class:`CoordFunction` or a value

        OUTPUT:

        - coordinate function resulting from the addition of ``self`` to
          ``other``

        TESTS::

            sage: M = Manifold(2, 'M', structure='topological')
            sage: X.<x,y> = M.chart()
            sage: f = X.function(x+y)
            sage: g = X.function(x*y)
            sage: f.__radd__(g)
            (x + 1)*y + x
            sage: f.__radd__(X.zero_function()) == f
            True
            sage: 2 + f  # indirect doctest
            x + y + 2

        """
        return self.__add__(other)  # since + is commutative

    def __iadd__(self, other):
        r"""
        In-place addition operator.

        INPUT:

        - ``other`` -- another instance of :class:`CoordFunction` or a value

        OUTPUT:

        - coordinate function resulting from the addition of ``self`` and
          ``other``

        TESTS::

            sage: M = Manifold(2, 'M', structure='topological')
            sage: X.<x,y> = M.chart()
            sage: f = X.function(x+y)
            sage: g = X.function(x*y)
            sage: f.__iadd__(g)
            (x + 1)*y + x
            sage: f += g; f
            (x + 1)*y + x
            sage: f.__iadd__(X.zero_function()) == f
            True

        """
        return self.__add__(other)

    def __rsub__(self, other):
        r"""
        Reflected subtraction operator: performs ``other - self``.

        INPUT:

        - ``other`` -- another instance of :class:`CoordFunction` or a value

        OUTPUT:

        - coordinate function resulting from the subtraction of ``self`` from
          ``other``

        TESTS::

            sage: M = Manifold(2, 'M', structure='topological')
            sage: X.<x,y> = M.chart()
            sage: f = X.function(x+y)
            sage: g = X.function(x*y)
            sage: f.__rsub__(g)
            (x - 1)*y - x
            sage: f.__rsub__(g) == -f + g
            True
            sage: 2 - f  # indirect doctest
            -x - y + 2

        """
        return self.__neg__().__add__(other)

    def __isub__(self, other):
        r"""
        In-place subtraction operator.

        INPUT:

        - ``other`` -- another instance of :class:`CoordFunction` or a value

        OUTPUT:

        - coordinate function resulting from the subtraction of ``other`` from
          ``self``

        TESTS::

            sage: M = Manifold(2, 'M', structure='topological')
            sage: X.<x,y> = M.chart()
            sage: f = X.function(x+y)
            sage: g = X.function(x*y)
            sage: f.__isub__(g)
            -(x - 1)*y + x
            sage: f -= g; f
            -(x - 1)*y + x
            sage: f.__isub__(X.zero_function()) == f
            True

        """
        return self.__sub__(other)

    def __rmul__(self, other):
        r"""
        Reflected multiplication operator: performs ``other * self``.

        INPUT:

        - ``other`` -- another instance of :class:`CoordFunction` or a value

        OUTPUT:

        - coordinate function resulting from the multiplication of ``other``
          by ``self``

        TESTS::

            sage: M = Manifold(2, 'M', structure='topological')
            sage: X.<x,y> = M.chart()
            sage: f = X.function(x+y)
            sage: g = X.function(x*y)
            sage: f.__rmul__(g)
            x^2*y + x*y^2
            sage: f.__rmul__(g) == g*f
            True
            sage: f.__rmul__(X.one_function()) == f
            True
            sage: 2*f  # indirect doctest
            2*x + 2*y

        """
        return self.__mul__(other)  # since * is commutative

    def __imul__(self, other):
        r"""
        In-place multiplication operator.

        INPUT:

        - ``other`` -- another instance of :class:`CoordFunction` or a value

        OUTPUT:

        - coordinate function resulting from the multiplication of ``self``
          by ``other``

        TESTS::

            sage: M = Manifold(2, 'M', structure='topological')
            sage: X.<x,y> = M.chart()
            sage: f = X.function(x+y)
            sage: g = X.function(x*y)
            sage: f.__imul__(g)
            x^2*y + x*y^2
            sage: f *= g; f
            x^2*y + x*y^2
            sage: f.__imul__(X.one_function()) == f
            True

        """
        return self.__mul__(other)

    def __rdiv__(self, other):
        r"""
        Reflected division operator: performs ``other / self``.

        INPUT:

        - ``other`` -- another instance of :class:`CoordFunction` or a value

        OUTPUT:

        - coordinate function resulting from the division of ``other`` by
          ``self``

        TESTS::

            sage: M = Manifold(2, 'M', structure='topological')
            sage: X.<x,y> = M.chart()
            sage: f = X.function(x+y)
            sage: g = X.function(x*y)
            sage: f.__rdiv__(g)
            x*y/(x + y)
            sage: f.__rdiv__(g) == g/f
            True
            sage: 2/f  # indirect doctest
            2/(x + y)

        """
        return self.__invert__().__mul__(other)

    def __idiv__(self, other):
        r"""
        In-place division operator.

        INPUT:

        - ``other`` -- another instance of :class:`CoordFunction` or a value

        OUTPUT:

        - coordinate function resulting from the division of ``self`` by
          ``other``

        TESTS::

            sage: M = Manifold(2, 'M', structure='topological')
            sage: X.<x,y> = M.chart()
            sage: f = X.function(x+y)
            sage: g = X.function(x*y)
            sage: f.__idiv__(g)
            (x + y)/(x*y)
            sage: f /= g; f
            (x + y)/(x*y)
            sage: f.__idiv__(X.one_function()) == f
            True

        """
        return self.__div__(other)

    # --------------------------------------------
    # Methods to be implemented by derived classes
    # --------------------------------------------

    def _repr_(self):
        r"""
        String representation of the object.

        TEST:

        This method must be implemented by derived classes; it is not
        implemented here::

            sage: M = Manifold(2, 'M', structure='topological')
            sage: X.<x,y> = M.chart()
            sage: from sage.manifolds.coord_func import CoordFunction
            sage: f = CoordFunction(X)
            sage: f._repr_()
            Traceback (most recent call last):
            ...
            NotImplementedError: CoordFunction._repr_ not implemented

        """
        raise NotImplementedError("CoordFunction._repr_ not implemented")

    def _latex_(self):
        r"""
        LaTeX representation of the object.

        TEST:

        This method must be implemented by derived classes; it is not
        implemented here::

            sage: M = Manifold(2, 'M', structure='topological')
            sage: X.<x,y> = M.chart()
            sage: from sage.manifolds.coord_func import CoordFunction
            sage: f = CoordFunction(X)
            sage: f._latex_()
            Traceback (most recent call last):
            ...
            NotImplementedError: CoordFunction._latex_ not implemented

        """
        raise NotImplementedError("CoordFunction._latex_ not implemented")

    def display(self):
        r"""
        Display the function in arrow notation.

        TEST:

        This method must be implemented by derived classes; it is not
        implemented here::

            sage: M = Manifold(2, 'M', structure='topological')
            sage: X.<x,y> = M.chart()
            sage: from sage.manifolds.coord_func import CoordFunction
            sage: f = CoordFunction(X)
            sage: f.display()
            Traceback (most recent call last):
            ...
            NotImplementedError: CoordFunction.display not implemented

        """
        raise NotImplementedError("CoordFunction.display not implemented")

    disp = display

    def expr(self):
        r"""
        Return some data that, along with the chart, is sufficient to
        reconstruct the object.

        For a symbolic coordinate function, this returns the symbol
        expression representing the function
        (see
        :meth:`sage.manifolds.coord_func_symb.CoordFunctionSymb.expr`)

        TEST:

        This method must be implemented by derived classes; it is not
        implemented here::

            sage: M = Manifold(2, 'M', structure='topological')
            sage: X.<x,y> = M.chart()
            sage: from sage.manifolds.coord_func import CoordFunction
            sage: f = CoordFunction(X)
            sage: f.expr()
            Traceback (most recent call last):
            ...
            NotImplementedError: CoordFunction.expr not implemented

        """
        raise NotImplementedError("CoordFunction.expr not implemented")

    def __call__(self, *coords, **options):
        r"""
        Compute the value of the function at specified coordinates.

        INPUT:

        - ``*coords`` -- list of coordinates `(x^1,...,x^n)` where the
          function `f` is to be evaluated
        - ``**options`` -- options to control the computation (e.g.
          simplification options)

        OUTPUT:

        - the value `f(x^1,...,x^n)`, where `f` is the current coordinate
          function.

        TEST:

        This method must be implemented by derived classes; it is not
        implemented here::

            sage: M = Manifold(2, 'M', structure='topological')
            sage: X.<x,y> = M.chart()
            sage: from sage.manifolds.coord_func import CoordFunction
            sage: f = CoordFunction(X)
            sage: f.__call__(2,-3)
            Traceback (most recent call last):
            ...
            NotImplementedError: CoordFunction.__call__ not implemented

        """
        raise NotImplementedError("CoordFunction.__call__ not implemented")

    def is_zero(self):
        r"""
        Return ``True`` if the function is zero and ``False`` otherwise.

        TEST:

        This method must be implemented by derived classes; it is not
        implemented here::

            sage: M = Manifold(2, 'M', structure='topological')
            sage: X.<x,y> = M.chart()
            sage: from sage.manifolds.coord_func import CoordFunction
            sage: f = CoordFunction(X)
            sage: f.is_zero()
            Traceback (most recent call last):
            ...
            NotImplementedError: CoordFunction.is_zero not implemented

        """
        raise NotImplementedError("CoordFunction.is_zero not implemented")

    def copy(self):
        r"""
        Return an exact copy of the object.

        OUTPUT:

        - an instance of :class:`CoordFunction`

        TEST:

        This method must be implemented by derived classes; it is not
        implemented here::

            sage: M = Manifold(2, 'M', structure='topological')
            sage: X.<x,y> = M.chart()
            sage: from sage.manifolds.coord_func import CoordFunction
            sage: f = CoordFunction(X)
            sage: f.copy()
            Traceback (most recent call last):
            ...
            NotImplementedError: CoordFunction.copy not implemented

        """
        raise NotImplementedError("CoordFunction.copy not implemented")

    def diff(self, coord):
        r"""
        Partial derivative with respect to a coordinate.

        INPUT:

        - ``coord`` -- either the coordinate `x^i` with respect
          to which the derivative of the coordinate function `f` is to be
          taken, or the index `i` labelling this coordinate (with the
          index convention defined on the chart domain via the parameter
          ``start_index``)

        OUTPUT:

        - instance of :class:`CoordFunction` representing the partial
          derivative `\frac{\partial f}{\partial x^i}`

        TEST:

        This method must be implemented by derived classes; it is not
        implemented here::

            sage: M = Manifold(2, 'M', structure='topological')
            sage: X.<x,y> = M.chart()
            sage: from sage.manifolds.coord_func import CoordFunction
            sage: f = CoordFunction(X)
            sage: f.diff(x)
            Traceback (most recent call last):
            ...
            NotImplementedError: CoordFunction.diff not implemented

        """
        raise NotImplementedError("CoordFunction.diff not implemented")

    def __eq__(self, other):
        r"""
        Comparison (equality) operator.

        INPUT:

        - ``other`` -- another instance of :class:`CoordFunction` or a value

        OUTPUT:

        - ``True`` if ``self`` is equal to ``other``,  or ``False`` otherwise

        TEST:

        This method must be implemented by derived classes; it is not
        implemented here::

            sage: M = Manifold(2, 'M', structure='topological')
            sage: X.<x,y> = M.chart()
            sage: from sage.manifolds.coord_func import CoordFunction
            sage: f = CoordFunction(X)
            sage: g = CoordFunction(X)
            sage: f == g
            Traceback (most recent call last):
            ...
            NotImplementedError: CoordFunction.__eq__ not implemented

        """
        raise NotImplementedError("CoordFunction.__eq__ not implemented")

    def __pos__(self):
        r"""
        Unary plus operator.

        OUTPUT:

        - an exact copy of ``self``

        TEST:

        This method must be implemented by derived classes; it is not
        implemented here::

            sage: M = Manifold(2, 'M', structure='topological')
            sage: X.<x,y> = M.chart()
            sage: from sage.manifolds.coord_func import CoordFunction
            sage: f = CoordFunction(X)
            sage: f.__pos__()
            Traceback (most recent call last):
            ...
            NotImplementedError: CoordFunction.__pos__ not implemented

        """
        raise NotImplementedError("CoordFunction.__pos__ not implemented")

    def __neg__(self):
        r"""
        Unary minus operator.

        OUTPUT:

        - the opposite of the coordinate function ``self``

        TEST:

        This method must be implemented by derived classes; it is not
        implemented here::

            sage: M = Manifold(2, 'M', structure='topological')
            sage: X.<x,y> = M.chart()
            sage: from sage.manifolds.coord_func import CoordFunction
            sage: f = CoordFunction(X)
            sage: f.__neg__()
            Traceback (most recent call last):
            ...
            NotImplementedError: CoordFunction.__neg__ not implemented

        """
        raise NotImplementedError("CoordFunction.__neg__ not implemented")

    def __invert__(self):
        r"""
        Inverse operator.

        If `f` denotes the current coordinate function and `K` the topological
        field over which the manifold is defined, the *inverse* of `f` is the
        coordinate function `1/f`, where `1` of the multiplicative identity
        of `K`.

        OUTPUT:

        - the inverse of the coordinate function ``self``

        TEST:

        This method must be implemented by derived classes; it is not
        implemented here::

            sage: M = Manifold(2, 'M', structure='topological')
            sage: X.<x,y> = M.chart()
            sage: from sage.manifolds.coord_func import CoordFunction
            sage: f = CoordFunction(X)
            sage: f.__invert__()
            Traceback (most recent call last):
            ...
            NotImplementedError: CoordFunction.__invert__ not implemented

        """
        raise NotImplementedError("CoordFunction.__invert__ not implemented")

    def __add__(self, other):
        r"""
        Addition operator.

        INPUT:

        - ``other`` -- another instance of :class:`CoordFunction` or a value

        OUTPUT:

        - coordinate function resulting from the addition of ``self`` and
          ``other``

        TEST:

        This method must be implemented by derived classes; it is not
        implemented here::

            sage: M = Manifold(2, 'M', structure='topological')
            sage: X.<x,y> = M.chart()
            sage: from sage.manifolds.coord_func import CoordFunction
            sage: f = CoordFunction(X)
            sage: f.__add__(2)
            Traceback (most recent call last):
            ...
            NotImplementedError: CoordFunction.__add__ not implemented

        """
        raise NotImplementedError("CoordFunction.__add__ not implemented")

    def __sub__(self, other):
        r"""
        Subtraction operator.

        INPUT:

        - ``other`` -- another instance of :class:`CoordFunction` or a value

        OUTPUT:

        - coordinate function resulting from the subtraction of ``other`` from
          ``self``

        TEST:

        This method must be implemented by derived classes; it is not
        implemented here::

            sage: M = Manifold(2, 'M', structure='topological')
            sage: X.<x,y> = M.chart()
            sage: from sage.manifolds.coord_func import CoordFunction
            sage: f = CoordFunction(X)
            sage: f.__sub__(2)
            Traceback (most recent call last):
            ...
            NotImplementedError: CoordFunction.__sub__ not implemented

        """
        raise NotImplementedError("CoordFunction.__sub__ not implemented")

    def __mul__(self, other):
        r"""
        Multiplication  operator.

        INPUT:

        - ``other`` -- another instance of :class:`CoordFunction` or a value

        OUTPUT:

        - coordinate function resulting from the multiplication of ``self`` by
          ``other``

        TEST:

        This method must be implemented by derived classes; it is not
        implemented here::

            sage: M = Manifold(2, 'M', structure='topological')
            sage: X.<x,y> = M.chart()
            sage: from sage.manifolds.coord_func import CoordFunction
            sage: f = CoordFunction(X)
            sage: f.__mul__(2)
            Traceback (most recent call last):
            ...
            NotImplementedError: CoordFunction.__mul__ not implemented

        """
        raise NotImplementedError("CoordFunction.__mul__ not implemented")

    def __div__(self, other):
        r"""
        Division  operator.

        INPUT:

        - ``other`` -- another instance of :class:`CoordFunction` or a value

        OUTPUT:

        - coordinate function resulting from the division of ``self`` by
          ``other``

        TEST:

        This method must be implemented by derived classes; it is not
        implemented here::

            sage: M = Manifold(2, 'M', structure='topological')
            sage: X.<x,y> = M.chart()
            sage: from sage.manifolds.coord_func import CoordFunction
            sage: f = CoordFunction(X)
            sage: f.__div__(2)
            Traceback (most recent call last):
            ...
            NotImplementedError: CoordFunction.__div__ not implemented

        """
        raise NotImplementedError("CoordFunction.__div__ not implemented")

    def exp(self):
        r"""
        Exponential of the coordinate function.

        OUTPUT:

        - coordinate function `\exp(f)`, where `f` is the current coordinate
          function.

        TEST:

        This method must be implemented by derived classes; it is not
        implemented here::

            sage: M = Manifold(2, 'M', structure='topological')
            sage: X.<x,y> = M.chart()
            sage: from sage.manifolds.coord_func import CoordFunction
            sage: f = CoordFunction(X)
            sage: f.exp()
            Traceback (most recent call last):
            ...
            NotImplementedError: CoordFunction.exp not implemented

        """
        raise NotImplementedError("CoordFunction.exp not implemented")


    def log(self, base=None):
        r"""
        Logarithm of the coordinate function.

        INPUT:

        - ``base`` -- (default: ``None``) base of the logarithm; if None, the
          natural logarithm (i.e. logarithm to base e) is returned

        OUTPUT:

        - coordinate function `\log_a(f)`, where `f` is the current coordinate
          function and `a` is the base

        TEST:

        This method must be implemented by derived classes; it is not
        implemented here::

            sage: M = Manifold(2, 'M', structure='topological')
            sage: X.<x,y> = M.chart()
            sage: from sage.manifolds.coord_func import CoordFunction
            sage: f = CoordFunction(X)
            sage: f.log()
            Traceback (most recent call last):
            ...
            NotImplementedError: CoordFunction.log not implemented

        """
        raise NotImplementedError("CoordFunction.log not implemented")


    def __pow__(self, exponent):
        r"""
        Power of the coordinate function.

        INPUT:

        - ``exponent`` -- the exponent

        OUTPUT:

        - coordinate function `f^a`, where `f` is the current coordinate
          function and `a` is the exponent

        TEST:

        This method must be implemented by derived classes; it is not
        implemented here::

            sage: M = Manifold(2, 'M', structure='topological')
            sage: X.<x,y> = M.chart()
            sage: from sage.manifolds.coord_func import CoordFunction
            sage: f = CoordFunction(X)
            sage: f.__pow__(2)
            Traceback (most recent call last):
            ...
            NotImplementedError: CoordFunction.__pow__ not implemented

        """
        raise NotImplementedError("CoordFunction.__pow__ not implemented")


    def sqrt(self):
        r"""
        Square root of the coordinate function.

        OUTPUT:

        - coordinate function `\sqrt{f}`, where `f` is the current coordinate
          function.

        TEST:

        This method must be implemented by derived classes; it is not
        implemented here::

            sage: M = Manifold(2, 'M', structure='topological')
            sage: X.<x,y> = M.chart()
            sage: from sage.manifolds.coord_func import CoordFunction
            sage: f = CoordFunction(X)
            sage: f.sqrt()
            Traceback (most recent call last):
            ...
            NotImplementedError: CoordFunction.sqrt not implemented

        """
        raise NotImplementedError("CoordFunction.sqrt not implemented")

    def cos(self):
        r"""
        Cosine of the coordinate function.

        OUTPUT:

        - coordinate function `\cos(f)`, where `f` is the current coordinate
          function.

        TEST:

        This method must be implemented by derived classes; it is not
        implemented here::

            sage: M = Manifold(2, 'M', structure='topological')
            sage: X.<x,y> = M.chart()
            sage: from sage.manifolds.coord_func import CoordFunction
            sage: f = CoordFunction(X)
            sage: f.cos()
            Traceback (most recent call last):
            ...
            NotImplementedError: CoordFunction.cos not implemented

        """
        raise NotImplementedError("CoordFunction.cos not implemented")

    def sin(self):
        r"""
        Sine of the coordinate function.

        OUTPUT:

        - coordinate function `\sin(f)`, where `f` is the current coordinate
          function.

        TEST:

        This method must be implemented by derived classes; it is not
        implemented here::

            sage: M = Manifold(2, 'M', structure='topological')
            sage: X.<x,y> = M.chart()
            sage: from sage.manifolds.coord_func import CoordFunction
            sage: f = CoordFunction(X)
            sage: f.sin()
            Traceback (most recent call last):
            ...
            NotImplementedError: CoordFunction.sin not implemented

        """
        raise NotImplementedError("CoordFunction.sin not implemented")

    def tan(self):
        r"""
        Tangent of the coordinate function.

        OUTPUT:

        - coordinate function `\tan(f)`, where `f` is the current coordinate
          function.

        TEST:

        This method must be implemented by derived classes; it is not
        implemented here::

            sage: M = Manifold(2, 'M', structure='topological')
            sage: X.<x,y> = M.chart()
            sage: from sage.manifolds.coord_func import CoordFunction
            sage: f = CoordFunction(X)
            sage: f.tan()
            Traceback (most recent call last):
            ...
            NotImplementedError: CoordFunction.tan not implemented

        """
        raise NotImplementedError("CoordFunction.tan not implemented")

    def arccos(self):
        r"""
        Arc cosine of the coordinate function.

        OUTPUT:

        - coordinate function `\arccos(f)`, where `f` is the current coordinate
          function.

        TEST:

        This method must be implemented by derived classes; it is not
        implemented here::

            sage: M = Manifold(2, 'M', structure='topological')
            sage: X.<x,y> = M.chart()
            sage: from sage.manifolds.coord_func import CoordFunction
            sage: f = CoordFunction(X)
            sage: f.arccos()
            Traceback (most recent call last):
            ...
            NotImplementedError: CoordFunction.arccos not implemented

        """
        raise NotImplementedError("CoordFunction.arccos not implemented")

    def arcsin(self):
        r"""
        Arc sine of the coordinate function.

        OUTPUT:

        - coordinate function `\arcsin(f)`, where `f` is the current coordinate
          function.

        TEST:

        This method must be implemented by derived classes; it is not
        implemented here::

            sage: M = Manifold(2, 'M', structure='topological')
            sage: X.<x,y> = M.chart()
            sage: from sage.manifolds.coord_func import CoordFunction
            sage: f = CoordFunction(X)
            sage: f.arcsin()
            Traceback (most recent call last):
            ...
            NotImplementedError: CoordFunction.arcsin not implemented

        """
        raise NotImplementedError("CoordFunction.arcsin not implemented")

    def arctan(self):
        r"""
        Arc tangent of the coordinate function.

        OUTPUT:

        - coordinate function `\arctan(f)`, where `f` is the current coordinate
          function.

        TEST:

        This method must be implemented by derived classes; it is not
        implemented here::

            sage: M = Manifold(2, 'M', structure='topological')
            sage: X.<x,y> = M.chart()
            sage: from sage.manifolds.coord_func import CoordFunction
            sage: f = CoordFunction(X)
            sage: f.arctan()
            Traceback (most recent call last):
            ...
            NotImplementedError: CoordFunction.arctan not implemented

        """
        raise NotImplementedError("CoordFunction.arctan not implemented")

    def cosh(self):
        r"""
        Hyperbolic cosine of the coordinate function.

        OUTPUT:

        - coordinate function `\cosh(f)`, where `f` is the current coordinate
          function.

        TEST:

        This method must be implemented by derived classes; it is not
        implemented here::

            sage: M = Manifold(2, 'M', structure='topological')
            sage: X.<x,y> = M.chart()
            sage: from sage.manifolds.coord_func import CoordFunction
            sage: f = CoordFunction(X)
            sage: f.cosh()
            Traceback (most recent call last):
            ...
            NotImplementedError: CoordFunction.cosh not implemented

        """
        raise NotImplementedError("CoordFunction.cosh not implemented")

    def sinh(self):
        r"""
        Hyperbolic sine of the coordinate function.

        OUTPUT:

        - coordinate function `\sinh(f)`, where `f` is the current coordinate
          function.

        TEST:

        This method must be implemented by derived classes; it is not
        implemented here::

            sage: M = Manifold(2, 'M', structure='topological')
            sage: X.<x,y> = M.chart()
            sage: from sage.manifolds.coord_func import CoordFunction
            sage: f = CoordFunction(X)
            sage: f.sinh()
            Traceback (most recent call last):
            ...
            NotImplementedError: CoordFunction.sinh not implemented

        """
        raise NotImplementedError("CoordFunction.sinh not implemented")

    def tanh(self):
        r"""
        Hyperbolic tangent of the coordinate function.

        OUTPUT:

        - coordinate function `\tanh(f)`, where `f` is the current coordinate
          function.

        TEST:

        This method must be implemented by derived classes; it is not
        implemented here::

            sage: M = Manifold(2, 'M', structure='topological')
            sage: X.<x,y> = M.chart()
            sage: from sage.manifolds.coord_func import CoordFunction
            sage: f = CoordFunction(X)
            sage: f.tanh()
            Traceback (most recent call last):
            ...
            NotImplementedError: CoordFunction.tanh not implemented

        """
        raise NotImplementedError("CoordFunction.tanh not implemented")

    def arccosh(self):
        r"""
        Inverse hyperbolic cosine of the coordinate function.

        OUTPUT:

        - coordinate function `\mathrm{arcosh}(f)`, where `f` is the current
          coordinate function.

        TEST:

        This method must be implemented by derived classes; it is not
        implemented here::

            sage: M = Manifold(2, 'M', structure='topological')
            sage: X.<x,y> = M.chart()
            sage: from sage.manifolds.coord_func import CoordFunction
            sage: f = CoordFunction(X)
            sage: f.arccosh()
            Traceback (most recent call last):
            ...
            NotImplementedError: CoordFunction.arccosh not implemented

        """
        raise NotImplementedError("CoordFunction.arccosh not implemented")

    def arcsinh(self):
        r"""
        Inverse hyperbolic sine of the coordinate function.

        OUTPUT:

        - coordinate function `\mathrm{arsinh}(f)`, where `f` is the current
          coordinate function.

        TEST:

        This method must be implemented by derived classes; it is not
        implemented here::

            sage: M = Manifold(2, 'M', structure='topological')
            sage: X.<x,y> = M.chart()
            sage: from sage.manifolds.coord_func import CoordFunction
            sage: f = CoordFunction(X)
            sage: f.arcsinh()
            Traceback (most recent call last):
            ...
            NotImplementedError: CoordFunction.arcsinh not implemented

        """
        raise NotImplementedError("CoordFunction.arcsinh not implemented")

    def arctanh(self):
        r"""
        Inverse hyperbolic tangent of the coordinate function.

        OUTPUT:

        - coordinate function `\mathrm{artanh}(f)`, where `f` is the current
          coordinate function.

        TEST:

        This method must be implemented by derived classes; it is not
        implemented here::

            sage: M = Manifold(2, 'M', structure='topological')
            sage: X.<x,y> = M.chart()
            sage: from sage.manifolds.coord_func import CoordFunction
            sage: f = CoordFunction(X)
            sage: f.arctanh()
            Traceback (most recent call last):
            ...
            NotImplementedError: CoordFunction.arctanh not implemented

        """
        raise NotImplementedError("CoordFunction.arctanh not implemented")


#*****************************************************************************

class MultiCoordFunction(SageObject):
    r"""
    Coordinate function to some Cartesian power of the base field.

    If `n` and `m` are two positive integers and `(U,\varphi)` is a chart on
    a topological manifold `M` of dimension `n` over a topological field `K`,
    a *multi-coordinate function* associated to `(U,\varphi)` is a map

    .. MATH::

        \begin{array}{llcl}
        f:& V \subset K^n & \longrightarrow & K^m \\
          & (x^1,\ldots,x^n) & \longmapsto & (f_1(x^1,\ldots,x^n),\ldots,
            f_m(x^1,\ldots,x^n)) ,
        \end{array}

    where `V` is the codomain of `\varphi`. In other words, `f` is a
    `K^m`-valued function of the coordinates associated to the chart
    `(U,\varphi)`. Each component `f_i` (`1\leq i \leq m`) is a coordinate
    function and is therefore stored as an instance of
    :class:`~sage.manifolds.coord_func.CoordFunction`.

    INPUT:

    - ``chart`` -- the chart `(U, \varphi)`
    - ``expressions`` -- list (or tuple) of length `m` of elements to
      construct the coordinate functions `f_i` (`1\leq i \leq m`); for
      symbolic coordinate functions, this must be symbolic expressions
      involving the chart coordinates, while for numerical coordinate functions,
      this must be data file names

    EXAMPLES:

    A function `f: V\subset \RR^2 \longrightarrow \RR^3`::

        sage: M = Manifold(2, 'M', structure='topological')
        sage: X.<x,y> = M.chart()
        sage: f = X.multifunction(x-y, x*y, cos(x)*exp(y)); f
        Coordinate functions (x - y, x*y, cos(x)*e^y) on the Chart (M, (x, y))
        sage: type(f)
        <class 'sage.manifolds.coord_func.MultiCoordFunction'>
        sage: f(x,y)
        (x - y, x*y, cos(x)*e^y)
        sage: latex(f)
        \left(x - y, x y, \cos\left(x\right) e^{y}\right)

    Each real-valued function `f_i` (`1\leq i \leq m`) composing `f` can be
    accessed via the square-bracket operator, by providing `i-1` as an
    argument::

        sage: f[0]
        x - y
        sage: f[1]
        x*y
        sage: f[2]
        cos(x)*e^y

    Each ``f[i-1]`` is an instance of
    :class:`~sage.manifolds.coord_func.CoordFunction`::

        sage: isinstance(f[0], sage.manifolds.coord_func.CoordFunction)
        True

    In the present case, ``f[i-1]`` belongs to the subclass
    :class:`~sage.manifolds.coord_func_symb.CoordFunctionSymb`::

        sage: type(f[0])
        <class 'sage.manifolds.coord_func_symb.CoordFunctionSymb'>
        sage: f[0].display()
        (x, y) |--> x - y

    An instance of class :class:`MultiCoordFunction` can represent a
    real-valued function (case `m=1`), although one should
    rather employ the class :class:`~sage.manifolds.coord_func.CoordFunction`
    for this purpose::


        sage: g = X.multifunction(x*y^2)
        sage: g(x,y)
        (x*y^2,)

    Evaluating the functions at specified coordinates::

        sage: f(1,2)
        (-1, 2, cos(1)*e^2)
        sage: var('a b')
        (a, b)
        sage: f(a,b)
        (a - b, a*b, cos(a)*e^b)
        sage: g(1,2)
        (4,)

    """
    def __init__(self, chart, expressions):
        r"""
        Construct a multi-coordinate function.

        TESTS::

            sage: M = Manifold(3, 'M', structure='topological')
            sage: X.<x,y,z> = M.chart()
            sage: f = X.multifunction(x+y+z, x*y*z); f
            Coordinate functions (x + y + z, x*y*z) on the Chart (M, (x, y, z))
            sage: type(f)
            <class 'sage.manifolds.coord_func.MultiCoordFunction'>
            sage: TestSuite(f).run()

        """
        self._chart = chart
        self._nc = len(self._chart._xx)   # number of coordinates
        self._nf = len(expressions)       # number of functions
        self._functions = tuple(chart.function(express)
                                for express in expressions)
        # Derived quantities:
        self._jacob = None
        self._jacob_det = None

    def _repr_(self):
        r"""
        String representation of the object.

        TESTS::

            sage: M = Manifold(2, 'M', structure='topological')
            sage: X.<x,y> = M.chart()
            sage: f = X.multifunction(x-y, x*y, cos(x)*exp(y))
            sage: f._repr_()
            'Coordinate functions (x - y, x*y, cos(x)*e^y) on the Chart (M, (x, y))'
            sage: repr(f)  # indirect doctest
            'Coordinate functions (x - y, x*y, cos(x)*e^y) on the Chart (M, (x, y))'

        """
        return "Coordinate functions {} on the {}".format(self._functions,
                                                          self._chart)

    def _latex_(self):
        r"""
        LaTeX representation of the object.

        TESTS::

            sage: M = Manifold(2, 'M', structure='topological')
            sage: X.<x,y> = M.chart()
            sage: f = X.multifunction(x-y, x*y, cos(x)*exp(y))
            sage: f._latex_()
            \left(x - y, x y, \cos\left(x\right) e^{y}\right)
            sage: latex(f)  # indirect doctest
            \left(x - y, x y, \cos\left(x\right) e^{y}\right)

        """
        return latex(self._functions)

    def expr(self):
        r"""
        Return a tuple of data, the item no. `i` begin sufficient to
        reconstruct the coordinate function no. `i`.

        In other words, if ``f`` is a multi-coordinate function, then
        ``f.chart().multifunction(*(f.expr()))`` results in a
        multi-coordinate function identical to ``f``.

        For a symbolic multi-coordinate function, :meth:`expr` returns the
        tuple of the symbolic expressions of the coordinate functions
        composing the object.

        EXAMPLE::

            sage: M = Manifold(2, 'M', structure='topological')
            sage: X.<x,y> = M.chart()
            sage: f = X.multifunction(x-y, x*y, cos(x)*exp(y))
            sage: f.expr()
            (x - y, x*y, cos(x)*e^y)
            sage: type(f.expr()[0])
            <type 'sage.symbolic.expression.Expression'>

        One shall always have::

            sage: f.chart().multifunction(*(f.expr())) == f
            True

        """
        return tuple(func.expr() for func in self._functions)

    def chart(self):
        r"""
        Return the chart w.r.t. which the multi-coordinate function is defined.

        OUTPUT:

        - an instance of :class:`~sage.manifolds.chart.Chart`

        EXAMPLE::

            sage: M = Manifold(2, 'M', structure='topological')
            sage: X.<x,y> = M.chart()
            sage: f = X.multifunction(x-y, x*y, cos(x)*exp(y))
            sage: f.chart()
            Chart (M, (x, y))
            sage: f.chart() is X
            True

        """
        return self._chart


    def __eq__(self, other):
        r"""
        Comparison (equality) operator.

        INPUT:

        - ``other`` -- another instance of :class:`MultiCoordFunction`

        OUTPUT:

        - ``True`` if ``self`` is equal to ``other``, ``False`` otherwise


        TESTS::

            sage: M = Manifold(2, 'M', structure='topological')
            sage: X.<x,y> = M.chart()
            sage: f = X.multifunction(x-y, x*y, cos(x*y))
            sage: f == X.multifunction(x-y, x*y)
            False
            sage: f == X.multifunction(x-y, x*y, 2)
            False
            sage: f == X.multifunction(x-y, x*y, cos(y*x))
            True
            sage: Y.<u,v> = M.chart()
            sage: f == Y.multifunction(u-v, u*v, cos(u*v))
            False

        """
        if other is self:
            return True
        if not isinstance(other, MultiCoordFunction):
            return False
        if other._chart != self._chart:
            return False
        if other._nf != self._nf:
            return False
        for i in range(self._nf):
            if other._functions[i] != self._functions[i]:
                return False
        return True

    def __ne__(self, other):
        r"""
        Inequality operator.

        INPUT:

        - ``other`` -- another instance of :class:`MultiCoordFunction`

        OUTPUT:

        - ``True`` if ``self`` is different from ``other``, ``False``
          otherwise

        TESTS::

            sage: M = Manifold(2, 'M', structure='topological')
            sage: X.<x,y> = M.chart()
            sage: f = X.multifunction(x-y, x*y, cos(x*y))
            sage: f != X.multifunction(x-y, x*y)
            True
            sage: f != X.multifunction(x, y, 2)
            True
            sage: f != X.multifunction(x-y, x*y, cos(x*y))
            False

        """
        return not (self == other)

    def __getitem__(self, index):
        r"""
        Return a specified function of the set represented by ``self``.

        INPUT:

        -- ``index`` -- index `i` of the function (`0\leq i \leq m-1`)

        OUTPUT

        -- instance of :class:`CoordFunction` representing the function

        TESTS::

            sage: M = Manifold(2, 'M', structure='topological')
            sage: X.<x,y> = M.chart()
            sage: f = X.multifunction(x-y, x*y, cos(x*y))
            sage: f.__getitem__(0)
            x - y
            sage: f.__getitem__(1)
            x*y
            sage: f.__getitem__(2)
            cos(x*y)
            sage: f[0], f[1], f[2]
            (x - y, x*y, cos(x*y))

        """
        return self._functions[index]

    def __call__(self, *coords, **options):
        r"""
        Compute the values of the functions at specified coordinates.

        INPUT:

        - ``*coords`` -- list of coordinates where the functions are to be
          evaluated
        - ``**options`` -- allows to pass some options, like
          ``simplify=False`` to disable simplification for symbolic
          coordinate functions

        OUTPUT:

        - tuple containing the values of the `m` functions.

        TESTS::

            sage: M = Manifold(2, 'M', structure='topological')
            sage: X.<x,y> = M.chart()
            sage: f = X.multifunction(x-y, x*y, cos(x*y))
            sage: f.__call__(2,3)
            (-1, 6, cos(6))
            sage: f(2,3)
            (-1, 6, cos(6))
            sage: f.__call__(x,y)
            (x - y, x*y, cos(x*y))

        """
        return tuple(func(*coords, **options) for func in self._functions)

    def jacobian(self):
        r"""
        Return the Jacobian matrix of the system of coordinate functions.

        ``jacobian()`` is a 2-dimensional array of size `m\times n`
        where `m` is the number of functions and `n` the number of coordinates,
        the generic element being `J_{ij} = \frac{\partial f_i}{\partial x^j}`
        with `1\leq i \leq m` (row index) and `1\leq j \leq n` (column index).

        Each `J_{ij}` is an instance of :class:`CoordFunction`.

        OUTPUT:

        - Jacobian matrix as a 2-dimensional array ``J`` of coordinate functions,
          ``J[i-1][j-1]`` being `J_{ij} = \frac{\partial f_i}{\partial x^j}`
          for `1\leq i \leq m` and `1\leq j \leq n`

        EXAMPLES:

        Jacobian of a set of 3 functions of 2 coordinates::

            sage: M = Manifold(2, 'M', structure='topological')
            sage: X.<x,y> = M.chart()
            sage: f = X.multifunction(x-y, x*y, y^3*cos(x))
            sage: f.jacobian()
            [[1, -1], [y, x], [-y^3*sin(x), 3*y^2*cos(x)]]

        Each element of the result is a coordinate function::

            sage: type(f.jacobian()[2][0])
            <class 'sage.manifolds.coord_func_symb.CoordFunctionSymb'>
            sage: f.jacobian()[2][0].display()
            (x, y) |--> -y^3*sin(x)

        Test of the computation::

            sage: [[f.jacobian()[i][j] == f[i].diff(j) for j in range(2)] for i in range(3)]
            [[True, True], [True, True], [True, True]]

        Test with ``start_index=1``::

            sage: M = Manifold(2, 'M', structure='topological', start_index=1)
            sage: X.<x,y> = M.chart()
            sage: f = X.multifunction(x-y, x*y, y^3*cos(x))
            sage: f.jacobian()
            [[1, -1], [y, x], [-y^3*sin(x), 3*y^2*cos(x)]]
            sage: [[f.jacobian()[i][j] == f[i].diff(j+1) for j in range(2)]  # note the j+1
            ....:                                         for i in range(3)]
            [[True, True], [True, True], [True, True]]

        """
        if self._jacob is None:
            xx = self._chart[:]  # coordinates x^j
            self._jacob = [[ self._functions[i].diff(xx[j])
                          for j in range(self._nc) ] for i in range(self._nf) ]
        return self._jacob

    def jacobian_det(self):
        r"""
        Return the Jacobian determinant of the system of functions.

        The number `m` of coordinate functions must equal the number `n` of
        coordinates.

        OUTPUT:

        - instance of :class:`CoordFunction` representing the determinant

        EXAMPLE:

        Jacobian determinant of a set of 2 functions of 2 coordinates::

            sage: M = Manifold(2, 'M', structure='topological')
            sage: X.<x,y> = M.chart()
            sage: f = X.multifunction(x-y, x*y)
            sage: f.jacobian_det()
            x + y

        The output of :meth:`jacobian_det` is an instance of
        :class:`CoordFunction` and can therefore be called on specific values
        of the coordinates, e.g. (x,y)=(1,2)::

            sage: type(f.jacobian_det())
            <class 'sage.manifolds.coord_func_symb.CoordFunctionSymb'>
            sage: f.jacobian_det().display()
            (x, y) |--> x + y
            sage: f.jacobian_det()(1,2)
            3

        The result is cached::

            sage: f.jacobian_det() is f.jacobian_det()
            True

        Test::

            sage: f.jacobian_det() == det(matrix([[f[i].diff(j).expr() for j in range(2)]
            ....:                                 for i in range(2)]))
            True

        Jacobian determinant of a set of 3 functions of 3 coordinates::

            sage: M = Manifold(3, 'M', structure='topological')
            sage: X.<x,y,z> = M.chart()
            sage: f = X.multifunction(x*y+z^2, z^2*x+y^2*z, (x*y*z)^3)
            sage: f.jacobian_det().display()
            (x, y, z) |--> 6*x^3*y^5*z^3 - 3*x^4*y^3*z^4 - 12*x^2*y^4*z^5 + 6*x^3*y^2*z^6

        Test::

            sage: f.jacobian_det() == det(matrix([[f[i].diff(j).expr() for j in range(3)]
            ....:                                 for i in range(3)]))
            True

        """
        def simple_determinant(aa):
            r"""
            Compute the determinant of a square matrix represented as an array.

            This function is based on Laplace's cofactor expansion.
            """
            n = len(aa)
            if n == 1:
                return aa[0][0]
            res = 0
            sign = True
            for i in range(n):
                b = []
                for k in range(i):
                    r = []
                    for l in range(1,n):
                       r.append(aa[k][l])
                    b.append(r)
                for k in range(i+1,n):
                    r = []
                    for l in range(1,n):
                       r.append(aa[k][l])
                    b.append(r)
                if sign:
                    res += aa[i][0] * simple_determinant(b)
                else:
                    res -= aa[i][0] * simple_determinant(b)
                sign = not sign
            return res

        if self._jacob_det is None:
            if (self._nf != self._nc):
                raise ValueError("the Jacobian matrix is not a square matrix")
            self.jacobian() # forces the computation of self._jacob
            self._jacob_det = simple_determinant(self._jacob)
        return self._jacob_det
