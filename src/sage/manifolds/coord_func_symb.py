r"""
Symbolic coordinate functions

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

This module implements symbolic coordinate functions via the class
:class:`CoordFunctionSymb`.

The class :class:`MultiCoordFunctionSymb` implements `K^m`-valued
functions of the coordinates of a chart.

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

from sage.symbolic.ring import SR
from sage.manifolds.coord_func import CoordFunction, MultiCoordFunction

class CoordFunctionSymb(CoordFunction):
    r"""
    Coordinate function with symbolic representation.

    If `(U,\varphi)` is a chart on a topological manifold `M` of dimension `n`
    over a topological field `K`,  a *coordinate function* associated to
    `(U,\varphi)` is a map

    .. MATH::

        \begin{array}{llcl}
        f:& V \subset K^n & \longrightarrow & K \\
          & (x^1,\ldots,x^n) & \longmapsto & f(x^1,\ldots,x^n),
        \end{array}

    where `V` is the codomain of `\varphi`. In other words, `f` is a
    `K`-valued function of the
    coordinates associated to the chart `(U,\varphi)`.

    INPUT:

    - ``chart`` -- the chart `(U, \varphi)`, as an instance of class
      :class:`~sage.manifolds.chart.Chart`
    - ``expression`` -- a symbolic expression representing `f(x^1,\ldots,x^n)`,
      where `(x^1,\ldots,x^n)` are the coordinates of the chart `(U, \varphi)`

    """

    nice_output = True # static flag for textbook-like output instead of the
                       # Pynac output for derivatives

    def __init__(self, chart, expression):
        r"""
        Construct a coordinate function.
        """
        CoordFunction.__init__(self, chart)
        self._express = SR(expression)

    # -------------------------------------------------------------
    # Methods to be implemented by derived classes of CoordFunction
    # -------------------------------------------------------------

    def _repr_(self):
        r"""
        String representation of the object.
        """
        if FunctionChart.nice_output:
            return str(ExpressionNice(self._express))
        else:
            return str(self._express)

    def _latex_(self):
        r"""
        LaTeX representation of the object.
        """
        raise NotImplementedError("CoordFunction._latex_ not implemented")

    def display(self):
        r"""
        Display the function in arrow notation.
        """
        raise NotImplementedError("CoordFunction.display not implemented")

    def __call__(self, *coords, **options):
        r"""
        Computes the value of the function at specified coordinates.

        INPUT:

        - ``*coords`` -- list of coordinates `(x^1,...,x^n)` where the
          function `f` is to be evaluated
        - ``**options`` -- options to control the computation (e.g.
          simplification options)

        OUTPUT:

        - the value `f(x^1,...,x^n)`, where `f` is the current coordinate
          function.

        """
        raise NotImplementedError("CoordFunction.__call__ not implemented")

    def is_zero(self):
        r"""
        Return ``True`` if the function is zero and ``False`` otherwise.

        """
        raise NotImplementedError("CoordFunction.is_zero not implemented")

    def __eq__(self, other):
        r"""
        Comparison (equality) operator.

        INPUT:

        - ``other`` -- another instance of :class:`CoordFunction` or a value

        OUTPUT:

        - ``True`` if ``self`` is equal to ``other``,  or ``False`` otherwise

        """
        raise NotImplementedError("CoordFunction.__eq__ not implemented")

    def __pos__(self):
        r"""
        Unary plus operator.

        OUTPUT:

        - an exact copy of ``self``

        """
        return NotImplementedError("CoordFunction.__pos__ not implemented")

    def __neg__(self):
        r"""
        Unary minus operator.

        OUTPUT:

        - the opposite of the coordinate function ``self``

        """
        return NotImplementedError("CoordFunction.__neg__ not implemented")

    def __invert__(self):
        r"""
        Inverse operator.

        If `f` denotes the current coordinate function and `K` the topological
        field over which the manifold is defined, the *inverse* of `f` is the
        coordinate function `1/f`, where `1` of the multiplicative identity
        of `K`.

        OUTPUT:

        - the inverse of the coordinate function ``self``

        """
        return NotImplementedError("CoordFunction.__invert__ not implemented")

    def __add__(self, other):
        r"""
        Addition operator.

        INPUT:

        - ``other`` -- another instance of :class:`CoordFunction` or a value

        OUTPUT:

        - coordinate function resulting from the addition of ``self`` and
          ``other``

        """
        return NotImplementedError("CoordFunction.__add__ not implemented")

    def __sub__(self, other):
        r"""
        Subtraction operator.

        INPUT:

        - ``other`` -- another instance of :class:`CoordFunction` or a value

        OUTPUT:

        - coordinate function resulting from the subtraction of ``other`` from
          ``self``

        """
        return NotImplementedError("CoordFunction.__sub__ not implemented")

    def __mul__(self, other):
        r"""
        Multiplication  operator.

        INPUT:

        - ``other`` -- another instance of :class:`CoordFunction` or a value

        OUTPUT:

        - coordinate function resulting from the multiplication of ``self`` by
          ``other``

        """
        return NotImplementedError("CoordFunction.__mul__ not implemented")

    def __div__(self, other):
        r"""
        Division  operator.

        INPUT:

        - ``other`` -- another instance of :class:`CoordFunction` or a value

        OUTPUT:

        - coordinate function resulting from the division of ``self`` by
          ``other``

        """
        return NotImplementedError("CoordFunction.__div__ not implemented")

    def exp(self):
        r"""
        Exponential of the coordinate function.

        OUTPUT:

        - coordinate function `\exp(f)`, where `f` is the current coordinate
          function.

        """
        return NotImplementedError("CoordFunction.exp not implemented")


    def log(self, base=None):
        r"""
        Logarithm of the coordinate function.

        INPUT:

        - ``base`` -- (default: ``None``) base of the logarithm; if None, the
          natural logarithm (i.e. logarithm to base e) is returned

        OUTPUT:

        - coordinate function `\log_a(f)`, where `f` is the current coordinate
          function and `a` is the base

        """
        return NotImplementedError("CoordFunction.log not implemented")


    def __pow__(self, exponent):
        r"""
        Power of the coordinate function.

        INPUT:

        - ``exponent`` -- the exponent

        OUTPUT:

        - coordinate function `f^a`, where `f` is the current coordinate
          function and `a` is the exponent

        """
        return NotImplementedError("CoordFunction.__pow__ not implemented")


    def sqrt(self):
        r"""
        Square root of the coordinate function.

        OUTPUT:

        - coordinate function `\sqrt{f}`, where `f` is the current coordinate
          function.

        """
        return NotImplementedError("CoordFunction.sqrt not implemented")

    def cos(self):
        r"""
        Cosine of the coordinate function.

        OUTPUT:

        - coordinate function `\cos(f)`, where `f` is the current coordinate
          function.

        """
        return NotImplementedError("CoordFunction.cos not implemented")

    def sin(self):
        r"""
        Sine of the coordinate function.

        OUTPUT:

        - coordinate function `\sin(f)`, where `f` is the current coordinate
          function.

        """
        return NotImplementedError("CoordFunction.sin not implemented")

    def tan(self):
        r"""
        Tangent of the coordinate function.

        OUTPUT:

        - coordinate function `\tan(f)`, where `f` is the current coordinate
          function.

        """
        return NotImplementedError("CoordFunction.tan not implemented")

    def arccos(self):
        r"""
        Arc cosine of the coordinate function.

        OUTPUT:

        - coordinate function `\arccos(f)`, where `f` is the current coordinate
          function.

        """
        return NotImplementedError("CoordFunction.arccos not implemented")

    def arcsin(self):
        r"""
        Arc sine of the coordinate function.

        OUTPUT:

        - coordinate function `\arcsin(f)`, where `f` is the current coordinate
          function.

        """
        return NotImplementedError("CoordFunction.arcsin not implemented")

    def arctan(self):
        r"""
        Arc tangent of the coordinate function.

        OUTPUT:

        - coordinate function `\arctan(f)`, where `f` is the current coordinate
          function.

        """
        return NotImplementedError("CoordFunction.arctan not implemented")

    def cosh(self):
        r"""
        Hyperbolic cosine of the coordinate function.

        OUTPUT:

        - coordinate function `\cosh(f)`, where `f` is the current coordinate
          function.

        """
        return NotImplementedError("CoordFunction.cosh not implemented")

    def sinh(self):
        r"""
        Hyperbolic sine of the coordinate function.

        OUTPUT:

        - coordinate function `\sinh(f)`, where `f` is the current coordinate
          function.

        """
        return NotImplementedError("CoordFunction.sinh not implemented")

    def tanh(self):
        r"""
        Hyperbolic tangent of the coordinate function.

        OUTPUT:

        - coordinate function `\tanh(f)`, where `f` is the current coordinate
          function.

        """
        return NotImplementedError("CoordFunction.tanh not implemented")

    def arccosh(self):
        r"""
        Inverse hyperbolic cosine of the coordinate function.

        OUTPUT:

        - coordinate function `\mathrm{arcosh}(f)`, where `f` is the current
          coordinate function.

        """
        return NotImplementedError("CoordFunction.arccosh not implemented")

    def arcsinh(self):
        r"""
        Inverse hyperbolic sine of the coordinate function.

        OUTPUT:

        - coordinate function `\mathrm{arsinh}(f)`, where `f` is the current
          coordinate function.

        """
        return NotImplementedError("CoordFunction.arcsinh not implemented")

    def arctanh(self):
        r"""
        Inverse hyperbolic tangent of the coordinate function.

        OUTPUT:

        - coordinate function `\mathrm{artanh}(f)`, where `f` is the current
          coordinate function.

        """
        return NotImplementedError("CoordFunction.arctanh not implemented")


#*****************************************************************************

class MultiCoordFunction(SageObject):
    r"""
    Base class for multi-coordinate functions

    If `(U,\varphi)` is a chart on a topological manifold `M` of dimension `n`
    over a topological field `K`,  a *multi-coordinate function* associated to
    `(U,\varphi)` is a map

    .. MATH::

        \begin{array}{llcl}
        f:& V \subset K^n & \longrightarrow & K^m \\
          & (x^1,\ldots,x^n) & \longmapsto & (f_1(x^1,\ldots,x^n),\ldots,
            f_m(x^1,\ldots,x^n)) ,
        \end{array}

    where `V` is the codomain of `\varphi`. In other words, `f` is a
    `K^m`-valued function of the coordinates associated to the chart
    `(U,\varphi)`. Each components `f_i` (`1\leq i \leq m`) is a coordinate
    function and is therefore stored as an instance of
    :class:`~sage.manifolds.coord_func.CoordFunction`.

    INPUT:

    - ``chart`` -- the chart `(U, \varphi)`
    - ``size`` -- the integer `m`

    """
    def __init__(self, chart, size):
        if not isinstance(chart, Chart):
            raise TypeError("The argument must be a chart.")
        self._chart = chart
        self._nc = len(self._chart._xx)  # number of coordinates
        self._nf = size        # number of functions
        self._functions = None # to be set be derived classes

    def _repr_(self):
        r"""
        String representation of the object.
        """
        return "Functions {} on the {}".format(self._functions, self._chart)

    def _latex_(self):
        r"""
        LaTeX representation of the object.
        """
        return latex(self._functions)

    def __eq__(self, other):
        r"""
        Comparison (equality) operator.

        INPUT:

        - ``other`` -- another instance of :class:`MultiCoordFunction`

        OUTPUT:

        - ``True`` if ``self`` is equal to ``other``,  or ``False`` otherwise

        """
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

        - True if ``self`` is different from ``other``,  or False otherwise

        """
        return not self.__eq__(other)

    def __getitem__(self, index):
        r"""
        Return a specified function of the set represented by ``self``.

        INPUT:

        -- ``index`` -- index `i` of the function (`0\leq i \leq m-1`)

        OUTPUT

        -- instance of :class:`CoordFunction` representing the function

        """
        return self._functions[index]

    def __call__(self, *coords, **options):
        r"""
        Compute the values of the functions at specified coordinates.

        INPUT:

        - ``*coords`` -- list of coordinates where the functions are to be
          evaluated
        - ``**options`` -- allows to pass ``simplify=False`` to disable the
          call of the simplification chain on the result

        OUTPUT:

        - the values of the `m` functions.

        """
        return tuple( self._functions[i](*coords, **options) for i in
                                                              range(self._nf) )
