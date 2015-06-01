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
from sage.structure.element import RingElement
from sage.misc.latex import latex
from sage.manifolds.coord_func import CoordFunction, MultiCoordFunction
from sage.manifolds.utilities import ExpressionNice, simplify_chain_real, \
                                     simplify_chain_generic

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

    omit_fargs  = False # static flag to govern whether or not
                        # the arguments of symbolic functions are printed

    def __init__(self, chart, expression):
        r"""
        Construct a coordinate function.
        """
        self._chart = chart
        self._nc = len(chart[:])    # number of coordinates
        self._express = SR(expression)
        # Definition of the simplification chain to be applied in
        # symbolic calculus:
        if self._chart.manifold().base_field() == 'real':
            self._simplify = simplify_chain_real
        else:
            self._simplify = simplify_chain_generic

    # -------------------------------------------------------------
    # Methods to be implemented by derived classes of CoordFunction
    # -------------------------------------------------------------

    def _repr_(self):
        r"""
        String representation of the object.
        """
        if CoordFunctionSymb.nice_output:
            return str(ExpressionNice(self._express))
        else:
            return str(self._express)

    def _latex_(self):
        r"""
        LaTeX representation of the object.
        """
        if CoordFunctionSymb.nice_output:
            return latex(ExpressionNice(self._express))
        else:
            return latex(self._express)

    def display(self):
        r"""
        Display the function in arrow notation.
        """
        from sage.tensor.modules.format_utilities import FormattedExpansion
        resu_txt = str((self._chart)[:]) + ' |--> ' + \
                   str(ExpressionNice(self._express))
        resu_latex = latex(self._chart[:]) + r' \mapsto' + \
                     latex(ExpressionNice(self._express))
        return FormattedExpansion(resu_txt, resu_latex)

    def __call__(self, *coords, **options):
        r"""
        Computes the value of the function at specified coordinates.

        INPUT:

        - ``*coords`` -- list of coordinates `(x^1,...,x^n)` where the
          function `f` is to be evaluated
        - ``**options`` -- allows to pass ``simplify=False`` to disable the
          call of the simplification chain on the result

        OUTPUT:

        - the value `f(x^1,...,x^n)`, where `f` is the current coordinate
          function.

        """
        #!# This should be the Python 2.7 form:
        # substitutions = {self._chart._xx[j]: coords[j] for j in
        #                                                      range(self._nc)}
        #
        # Here we use a form compatible with Python 2.6:
        substitutions = dict([(self._chart._xx[j], coords[j]) for j in
                                                              range(self._nc)])
        resu = self._express.subs(substitutions)
        if 'simplify' in options:
            if options['simplify']:
                return self._simplify(resu)
            else:
                return resu
        else:
            return self._simplify(resu)

    def is_zero(self):
        r"""
        Return ``True`` if the function is zero and ``False`` otherwise.

        EXAMPLES:

        Coordinate functions associated to a 2-dimensional chart::

            sage: M = TopManifold(2, 'M')
            sage: X.<x,y> = M.chart()
            sage: f = X.function_symb(x^2+3*y+1)
            sage: f.is_zero()
            False
            sage: f == 0
            False
            sage: g = X.function_symb(0)
            sage: g.is_zero()
            True
            sage: g == 0
            True
            sage: X.zero_function().is_zero()
            True
            sage: X.zero_function() == 0
            True

        """
        return self._express.is_zero()

    def __eq__(self, other):
        r"""
        Comparison (equality) operator.

        INPUT:

        - ``other`` -- another instance of :class:`CoordFunction` or a value

        OUTPUT:

        - ``True`` if ``self`` is equal to ``other``,  or ``False`` otherwise

        TESTS:

        Coordinate functions associated to a 2-dimensional chart::

            sage: M = TopManifold(2, 'M')
            sage: X.<x,y> = M.chart()
            sage: f = X.function_symb(x+y^2)
            sage: g = X.function_symb(x+y^2)
            sage: f.__eq__(g)
            True
            sage: f == g
            True
            sage: f.__eq__(1)
            False
            sage: h = X.function_symb(1)
            sage: h.__eq__(1)
            True
            sage: h.__eq__(f)
            False
            sage: h.__eq__(0)
            False
            sage: X.function_symb(0).__eq__(0)
            True
            sage: X.zero_function().__eq__(0)
            True

        """
        if isinstance(other, CoordFunctionSymb):
            if other._chart != self._chart:
                return False
            else:
                return bool(other._express == self._express)
        else:
            return bool(self._express == other)

    def __pos__(self):
        r"""
        Unary plus operator.

        OUTPUT:

        - an exact copy of ``self``

        TESTS:

        Coordinate functions associated to a 2-dimensional chart::

            sage: M = TopManifold(2, 'M')
            sage: X.<x,y> = M.chart()
            sage: f = X.function_symb(x+y^2)
            sage: g = +f; g
            y^2 + x
            sage: type(g)
            <class 'sage.manifolds.coord_func_symb.CoordFunctionSymb'>
            sage: g == f
            True

        """
        return CoordFunctionSymb(self._chart, self._express)

    def __neg__(self):
        r"""
        Unary minus operator.

        OUTPUT:

        - the opposite of the coordinate function ``self``

        TESTS:

        Coordinate functions associated to a 2-dimensional chart::

            sage: M = TopManifold(2, 'M')
            sage: X.<x,y> = M.chart()
            sage: f = X.function_symb(x+y^2)
            sage: g = -f; g
            -y^2 - x
            sage: type(g)
            <class 'sage.manifolds.coord_func_symb.CoordFunctionSymb'>
            sage: -g == f
            True

        """
        return CoordFunctionSymb(self._chart, self._simplify(-self._express))

    def __invert__(self):
        r"""
        Inverse operator.

        If `f` denotes the current coordinate function and `K` the topological
        field over which the manifold is defined, the *inverse* of `f` is the
        coordinate function `1/f`, where `1` of the multiplicative identity
        of `K`.

        OUTPUT:

        - the inverse of the coordinate function ``self``

        TESTS:

        Coordinate functions associated to a 2-dimensional chart::

            sage: M = TopManifold(2, 'M')
            sage: X.<x,y> = M.chart()
            sage: f = X.function_symb(1+x^2+y^2)
            sage: g = f.__invert__(); g
            1/(x^2 + y^2 + 1)
            <class 'sage.manifolds.coord_func_symb.CoordFunctionSymb'>
            sage: g == ~f
            True
            sage: g.__invert__() == f
            True

        """
        return CoordFunctionSymb(self._chart,
                                 self._simplify(SR(1) / self._express))

    def __add__(self, other):
        r"""
        Addition operator.

        INPUT:

        - ``other`` -- another instance of :class:`CoordFunction` or a value

        OUTPUT:

        - coordinate function resulting from the addition of ``self`` and
          ``other``

        TESTS::

            sage: M = TopManifold(2, 'M')
            sage: X.<x,y> = M.chart()
            sage: f = X.function_symb(x+y^2)
            sage: g = X.function_symb(x+1)
            sage: s = f.__add__(g); s.display()
            (x, y) |--> y^2 + 2*x + 1
            sage: type(s)
            <class 'sage.manifolds.coord_func_symb.CoordFunctionSymb'>
            sage: f.__add__(0).display()
            (x, y) |--> y^2 + x
            sage: f.__add__(1).display()
            (x, y) |--> y^2 + x + 1
            sage: f.__add__(pi).display()
            (x, y) |--> pi + y^2 + x
            sage: f.__add__(x).display()
            (x, y) |--> y^2 + 2*x
            sage: f.__add__(-f).display()
            (x, y) |--> 0

        """
        if isinstance(other, CoordFunctionSymb):
            if other._chart != self._chart:
                raise ValueError("two coordinate functions not defined on " +
                                 "the same chart cannot be added")
            res = self._simplify(self._express + other._express)
        elif isinstance(other, (int, RingElement)):
            res = self._simplify(self._express + other)
        else:
            # addition to a numerical coord. function shall fall in this case
            return other.__radd__(self)
        if res == 0:
            return self._chart._zero_function
        else:
            return CoordFunctionSymb(self._chart, res)

    def __sub__(self, other):
        r"""
        Subtraction operator.

        INPUT:

        - ``other`` -- another instance of :class:`CoordFunction` or a value

        OUTPUT:

        - coordinate function resulting from the subtraction of ``other`` from
          ``self``

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

        """
        raise NotImplementedError("CoordFunction.__div__ not implemented")

    def exp(self):
        r"""
        Exponential of the coordinate function.

        OUTPUT:

        - coordinate function `\exp(f)`, where `f` is the current coordinate
          function.

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

        """
        raise NotImplementedError("CoordFunction.__pow__ not implemented")


    def sqrt(self):
        r"""
        Square root of the coordinate function.

        OUTPUT:

        - coordinate function `\sqrt{f}`, where `f` is the current coordinate
          function.

        """
        raise NotImplementedError("CoordFunction.sqrt not implemented")

    def cos(self):
        r"""
        Cosine of the coordinate function.

        OUTPUT:

        - coordinate function `\cos(f)`, where `f` is the current coordinate
          function.

        """
        raise NotImplementedError("CoordFunction.cos not implemented")

    def sin(self):
        r"""
        Sine of the coordinate function.

        OUTPUT:

        - coordinate function `\sin(f)`, where `f` is the current coordinate
          function.

        """
        raise NotImplementedError("CoordFunction.sin not implemented")

    def tan(self):
        r"""
        Tangent of the coordinate function.

        OUTPUT:

        - coordinate function `\tan(f)`, where `f` is the current coordinate
          function.

        """
        raise NotImplementedError("CoordFunction.tan not implemented")

    def arccos(self):
        r"""
        Arc cosine of the coordinate function.

        OUTPUT:

        - coordinate function `\arccos(f)`, where `f` is the current coordinate
          function.

        """
        raise NotImplementedError("CoordFunction.arccos not implemented")

    def arcsin(self):
        r"""
        Arc sine of the coordinate function.

        OUTPUT:

        - coordinate function `\arcsin(f)`, where `f` is the current coordinate
          function.

        """
        raise NotImplementedError("CoordFunction.arcsin not implemented")

    def arctan(self):
        r"""
        Arc tangent of the coordinate function.

        OUTPUT:

        - coordinate function `\arctan(f)`, where `f` is the current coordinate
          function.

        """
        raise NotImplementedError("CoordFunction.arctan not implemented")

    def cosh(self):
        r"""
        Hyperbolic cosine of the coordinate function.

        OUTPUT:

        - coordinate function `\cosh(f)`, where `f` is the current coordinate
          function.

        """
        raise NotImplementedError("CoordFunction.cosh not implemented")

    def sinh(self):
        r"""
        Hyperbolic sine of the coordinate function.

        OUTPUT:

        - coordinate function `\sinh(f)`, where `f` is the current coordinate
          function.

        """
        raise NotImplementedError("CoordFunction.sinh not implemented")

    def tanh(self):
        r"""
        Hyperbolic tangent of the coordinate function.

        OUTPUT:

        - coordinate function `\tanh(f)`, where `f` is the current coordinate
          function.

        """
        raise NotImplementedError("CoordFunction.tanh not implemented")

    def arccosh(self):
        r"""
        Inverse hyperbolic cosine of the coordinate function.

        OUTPUT:

        - coordinate function `\mathrm{arcosh}(f)`, where `f` is the current
          coordinate function.

        """
        raise NotImplementedError("CoordFunction.arccosh not implemented")

    def arcsinh(self):
        r"""
        Inverse hyperbolic sine of the coordinate function.

        OUTPUT:

        - coordinate function `\mathrm{arsinh}(f)`, where `f` is the current
          coordinate function.

        """
        raise NotImplementedError("CoordFunction.arcsinh not implemented")

    def arctanh(self):
        r"""
        Inverse hyperbolic tangent of the coordinate function.

        OUTPUT:

        - coordinate function `\mathrm{artanh}(f)`, where `f` is the current
          coordinate function.

        """
        raise NotImplementedError("CoordFunction.arctanh not implemented")


#*****************************************************************************

class MultiCoordFunctionSymb(MultiCoordFunction):
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
