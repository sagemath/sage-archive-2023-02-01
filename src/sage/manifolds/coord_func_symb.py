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

    EXAMPLES:

    A symbolic coordinate function associated with a 2-dimensional chart::

        sage: M = TopManifold(2, 'M')
        sage: X.<x,y> = M.chart()
        sage: f = X.function(x^2+3*y+1)
        sage: type(f)
        <class 'sage.manifolds.coord_func_symb.CoordFunctionSymb'>
        sage: f.display()
        (x, y) |--> x^2 + 3*y + 1
        sage: f(x,y)
        x^2 + 3*y + 1

    The symbolic expression is returned when asking the direct display of
    the function::

        sage: f
        x^2 + 3*y + 1
        sage: latex(f)
        x^{2} + 3 \, y + 1

    A similar output is obtained by means of the method :meth:`expr`::

        sage: f.expr()
        x^2 + 3*y + 1

    The value of the function at specified coordinates is obtained by means
    of the standard parentheses notation::

        sage: f(2,-1)
        2
        sage: var('a b')
        (a, b)
        sage: f(a,b)
        a^2 + 3*b + 1

    An unspecified coordinate function::

        sage: g = X.function(function('G', x, y))
        sage: g
        G(x, y)
        sage: g.display()
        (x, y) |--> G(x, y)
        sage: g.expr()
        G(x, y)
        sage: g(2,3)
        G(2, 3)

    Coordinate functions can be compared to other values::

        sage: f = X.function(x^2+3*y+1)
        sage: f == 2
        False
        sage: f == x^2 + 3*y + 1
        True
        sage: g = X.function(x*y)
        sage: f == g
        False
        sage: h = X.function(x^2+3*y+1)
        sage: f == h
        True

    .. RUBRIC:: Differences between ``CoordFunctionSymb`` and callable symbolic
      expressions

    Callable symbolic expressions are defined directly from symbolic
    expressions of the coordinates::

        sage: f0(x,y) = x^2 + 3*y + 1
        sage: type(f0)
        <type 'sage.symbolic.expression.Expression'>
        sage: f0
        (x, y) |--> x^2 + 3*y + 1
        sage: f0(x,y)
        x^2 + 3*y + 1

    To get an output similar to that of ``f0`` for the coordinate function
    ``f``, we must use the method :meth:`display`::

        sage: f
        x^2 + 3*y + 1
        sage: f.display()
        (x, y) |--> x^2 + 3*y + 1
        sage: f(x,y)
        x^2 + 3*y + 1

    More importantly, instances of :class:`CoordFunctionSymb` differ from
    callable symbolic expression by the automatic simplifications in all
    operations. For instance, adding the two callable symbolic expressions::

        sage: f0(x,y,z) = cos(x)^2 ; g0(x,y,z) = sin(x)^2

    results in::

        sage: f0 + g0
        (x, y, z) |--> cos(x)^2 + sin(x)^2

    To get 1,  one has to call
    :meth:`~sage.symbolic.expression.Expression.simplify_trig`::

        sage: (f0 + g0).simplify_trig()
        (x, y, z) |--> 1

    On the contrary, the sum of the corresponding :class:`CoordFunctionSymb`
    instances is automatically simplified (see
    :func:`~sage.manifolds.utilities.simplify_chain_real`
    and :func:`~sage.manifolds.utilities.simplify_chain_generic` for details)::

        sage: f = X.function(cos(x)^2) ; g = X.function(sin(x)^2)
        sage: f + g
        1

    Another difference regards the display of partial derivatives: for callable
    symbolic functions, it relies on Pynac notation ``D[0]``, ``D[1]``, etc.::

        sage: g = function('g', x, y)
        sage: f0(x,y) = diff(g, x) + diff(g, y)
        sage: f0
        (x, y) |--> D[0](g)(x, y) + D[1](g)(x, y)

    while for coordinate functions, the display is more "textbook" like::

        sage: f = X.function(diff(g, x) + diff(g, y))
        sage: f
        d(g)/dx + d(g)/dy

    The difference is even more dramatic on LaTeX outputs::

        sage: latex(f0)
        \left( x, y \right) \ {\mapsto} \ D[0]\left(g\right)\left(x, y\right) + D[1]\left(g\right)\left(x, y\right)
        sage: latex(f)
        \frac{\partial\,g}{\partial x} + \frac{\partial\,g}{\partial y}

    Note that this regards only the display of coordinate functions:
    internally, the Pynac notation is still used, as we can check by asking
    for the symbolic expression stored in `f`::

        sage: f.expr()
        D[0](g)(x, y) + D[1](g)(x, y)

    One can switch to Pynac notation via the command
    :func:`~sage.manifolds.utilities.nice_derivatives`::

        sage: nice_derivatives(False)
        sage: latex(f)
        D[0]\left(g\right)\left(x, y\right) + D[1]\left(g\right)\left(x, y\right)
        sage: nice_derivatives(True)
        sage: latex(f)
        \frac{\partial\,g}{\partial x} + \frac{\partial\,g}{\partial y}

    Another difference between :class:`CoordFunctionSymb` and
    callable symbolic expression is the possibility to switch off the display
    of the arguments of unspecified functions. Consider for instance::

        sage: f = X.function(function('u', x, y) * function('v', x, y))
        sage: f
        u(x, y)*v(x, y)
        sage: f0(x,y) = function('u', x, y) * function('v', x, y)
        sage: f0
        (x, y) |--> u(x, y)*v(x, y)

    If there is a clear understanding that `u` and `v` are functions of
    `(x,y)`, the explicit mention of the latter can be cumbersome in lengthy
    tensor expressions. We can switch it off by::

        sage: omit_function_args(True)
        sage: f
        u*v

    Note that neither the callable symbolic expression ``f0`` nor the internal
    expression of ``f`` is affected by the above command::

        sage: f0
        (x, y) |--> u(x, y)*v(x, y)
        sage: f.expr()
        u(x, y)*v(x, y)

    We revert to the default behavior by::

        sage: omit_function_args(False)
        sage: f
        u(x, y)*v(x, y)

    """

    _nice_output = True # static flag for textbook-like output instead of the
                        # Pynac output for derivatives

    _omit_fargs  = False # static flag to govern whether or not
                         # the arguments of symbolic functions are printed

    def __init__(self, chart, expression):
        r"""
        Construct a coordinate function.

        TESTS:

        Coordinate function on a real manifold::

            sage: M = TopManifold(2, 'M')
            sage: X.<x,y> = M.chart()
            sage: f = X.function(1+x*y); f
            x*y + 1
            sage: type(f)
            <class 'sage.manifolds.coord_func_symb.CoordFunctionSymb'>
            sage: TestSuite(f).run()

        Coordinate function on a complex manifold::

            sage: N = TopManifold(2, 'N', field='complex')
            sage: Y.<z,w> = N.chart()
            sage: g = Y.function(i*z + 2*w); g
            2*w + I*z
            sage: TestSuite(g).run()

        """
        self._chart = chart
        self._nc = len(chart[:])  # number of coordinates
        self._express = SR(expression)  # symbolic expression enforced
        # Definition of the simplification chain to be applied in
        # symbolic calculus:
        if self._chart.manifold().base_field() == 'real':
            self._simplify = simplify_chain_real
        else:
            self._simplify = simplify_chain_generic
        # Derived quantities:
        self._der = None  # list of partial derivatives (to be set by diff()
                          # and unset by del_derived())

    # -------------------------------------------------------------
    # Methods to be implemented by derived classes of CoordFunction
    # -------------------------------------------------------------

    def _repr_(self):
        r"""
        String representation of the object.

        TESTS::

            sage: M = TopManifold(2, 'M')
            sage: X.<x,y> = M.chart()
            sage: f = X.function(1+x*y)
            sage: f._repr_()
            'x*y + 1'
            sage: repr(f)  # indirect doctest
            'x*y + 1'
            sage: f  # indirect doctest
            x*y + 1

        """
        if CoordFunctionSymb._nice_output:
            return str(ExpressionNice(self._express))
        else:
            return str(self._express)

    def _latex_(self):
        r"""
        LaTeX representation of the object.

        TESTS::

            sage: M = TopManifold(2, 'M')
            sage: X.<x,y> = M.chart()
            sage: f = X.function(cos(x*y/2))
            sage: f._latex_()
            \cos\left(\frac{1}{2} \, x y\right)
            sage: latex(f)  # indirect doctest
            \cos\left(\frac{1}{2} \, x y\right)

        """
        if CoordFunctionSymb._nice_output:
            return latex(ExpressionNice(self._express))
        else:
            return latex(self._express)

    def display(self):
        r"""
        Display the function in arrow notation.

        The output is either text-formatted (console mode) or LaTeX-formatted
        (notebook mode).

        EXAMPLES:

        Coordinate function on a 2-dimensional manifold::

            sage: M = TopManifold(2, 'M')
            sage: X.<x,y> = M.chart()
            sage: f = X.function(cos(x*y/2))
            sage: f.display()
            (x, y) |--> cos(1/2*x*y)
            sage: latex(f.display())
            \left(x, y\right) \mapsto \cos\left(\frac{1}{2} \, x y\right)

        A shortcut is ``disp()``::

            sage: f.disp()
            (x, y) |--> cos(1/2*x*y)

        Display of the zero function::

            sage: X.zero_function().display()
            (x, y) |--> 0

        """
        from sage.tensor.modules.format_utilities import FormattedExpansion
        resu_txt = str((self._chart)[:]) + ' |--> ' + \
                   str(ExpressionNice(self._express))
        resu_latex = latex(self._chart[:]) + r' \mapsto' + \
                     latex(ExpressionNice(self._express))
        return FormattedExpansion(resu_txt, resu_latex)

    disp = display

    def expr(self):
        r"""
        Return the symbolic expression representing the image of the coordinate
        function.

        OUTPUT:

        - symbolic expression, involving the chart coordinates (instance of
          :class:`sage.symbolic.expression.Expression`)

        EXAMPLES:

        Coordinate function of a 2-dimensional manifold::

            sage: M = TopManifold(2, 'M')
            sage: X.<x,y> = M.chart()
            sage: f = X.function(x^2+3*y+1)
            sage: f.expr()
            x^2 + 3*y + 1
            sage: type(f.expr())
            <type 'sage.symbolic.expression.Expression'>

        For a symbolic coordinate function, one shall always have::

            sage: bool( f.expr() == f(*(f.chart()[:])) )
            True

        The method :meth:`expr` is useful for accessing to all the
        symbolic expression functionalities in Sage; for instance::

            sage: var('a')
            a
            sage: f = X.function(a*x*y); f.display()
            (x, y) |--> a*x*y
            sage: f.expr()
            a*x*y
            sage: f.expr().subs(a=2)
            2*x*y

        Note that for substituting the value of a coordinate, the function call
        can be used as well::

            sage: f(x,3)
            3*a*x
            sage: bool( f(x,3) == f.expr().subs(y=3) )
            True

        """
        return self._express

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

        TESTS::

            sage: M = TopManifold(2, 'M')
            sage: X.<x,y> = M.chart()
            sage: f = X.function(sin(x*y))
            sage: f.__call__(-2, 3)
            -sin(6)
            sage: f(-2, 3)
            -sin(6)
            sage: var('a b')
            (a, b)
            sage: f.__call__(a, b)
            sin(a*b)
            sage: f(a,b)
            sin(a*b)
            sage: f.__call__(pi, 1)
            0
            sage: f.__call__(pi, 1/2)
            1

        """
        if len(coords) != self._nc:
            raise ValueError("bad number of coordinates")
        substitutions = dict(zip(self._chart._xx, coords))
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
            sage: f = X.function(x^2+3*y+1)
            sage: f.is_zero()
            False
            sage: f == 0
            False
            sage: g = X.function(0)
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

    def copy(self):
        r"""
        Return an exact copy of the object.

        OUTPUT:

        - an instance of :class:`CoordFunctionSymb`

        EXAMPLE::

            sage: M = TopManifold(2, 'M')
            sage: X.<x,y> = M.chart()
            sage: f = X.function(x+y^2)
            sage: g = f.copy(); g
            y^2 + x

        By construction, ``g`` is identical to ``f``::

            sage: type(g) == type(f)
            True
            sage: g == f
            True

        but it is not the same object::

            sage: g is f
            False

        """
        return CoordFunctionSymb(self._chart, self._express)


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

        - instance of :class:`CoordFunctionSymb` representing the partial
          derivative `\frac{\partial f}{\partial x^i}`

        EXAMPLES:

        Partial derivatives of a 2-dimensional coordinate function::

            sage: M = TopManifold(2, 'M')
            sage: X.<x,y> = M.chart()
            sage: f = X.function(x^2+3*y+1); f
            x^2 + 3*y + 1
            sage: f.diff(x)
            2*x
            sage: f.diff(y)
            3

        Each partial derivatives is itself a coordinate function::

            sage: type(f.diff(x))
            <class 'sage.manifolds.coord_func_symb.CoordFunctionSymb'>

        An index can be used instead of the coordinate symbol::

            sage: f.diff(0)
            2*x
            sage: f.diff(1)
            3

        The index range depends on the convention used on the chart's domain::

            sage: M = TopManifold(2, 'M', start_index=1)
            sage: X.<x,y> = M.chart()
            sage: f = X.function(x^2+3*y+1)
            sage: f.diff(0)
            Traceback (most recent call last):
            ...
            ValueError: coordinate index out of range
            sage: f.diff(1)
            2*x
            sage: f.diff(2)
            3

        """
        from sage.calculus.functional import diff
        from sage.rings.integer import Integer
        if self._der is None:
            # the list of partial derivatives has to be updated
            self._der = [CoordFunctionSymb(self._chart,
                                           self._simplify(diff(self._express,
                                                               xx)))
                         for xx in self._chart[:]]
        if isinstance(coord, (int, Integer)):
            # NB: for efficiency, we access directly to the "private" attributes
            # of other classes. A more conventional OOP writing would be
            # coordsi = coord - self._chart.domain().start_index()
            coordsi = coord - self._chart._domain._sindex
            if coordsi < 0 or coordsi >= self._nc:
                raise ValueError("coordinate index out of range")
            return self._der[coordsi]
        else:
            return self._der[self._chart[:].index(coord)]

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
            sage: f = X.function(x+y^2)
            sage: g = X.function(x+y^2)
            sage: f.__eq__(g)
            True
            sage: f == g
            True
            sage: f.__eq__(1)
            False
            sage: h = X.function(1)
            sage: h.__eq__(1)
            True
            sage: h.__eq__(f)
            False
            sage: h.__eq__(0)
            False
            sage: X.function(0).__eq__(0)
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

        This method is identical to :meth:`copy`.

        OUTPUT:

        - an exact copy of ``self``

        TESTS:

        Coordinate functions associated to a 2-dimensional chart::

            sage: M = TopManifold(2, 'M')
            sage: X.<x,y> = M.chart()
            sage: f = X.function(x+y^2)
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
            sage: f = X.function(x+y^2)
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
            sage: f = X.function(1+x^2+y^2)
            sage: g = f.__invert__(); g
            1/(x^2 + y^2 + 1)
            sage: type(g)
            <class 'sage.manifolds.coord_func_symb.CoordFunctionSymb'>
            sage: g == ~f
            True
            sage: g.__invert__() == f
            True

        """
        return CoordFunctionSymb(self._chart,
                                 self._simplify(SR(1) / self._express))
        # NB: self._express.__invert__() would return 1/self._express
        # (cf. the code of __invert__ in src/sage/symbolic/expression.pyx)
        # Here we prefer SR(1)/self._express

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
            sage: f = X.function(x+y^2)
            sage: g = X.function(x+1)
            sage: s = f.__add__(g); s.display()
            (x, y) |--> y^2 + 2*x + 1
            sage: type(s)
            <class 'sage.manifolds.coord_func_symb.CoordFunctionSymb'>
            sage: f.__add__(0).display()
            (x, y) |--> y^2 + x
            sage: f.__add__(X.zero_function()).display()
            (x, y) |--> y^2 + x
            sage: f.__add__(1).display()
            (x, y) |--> y^2 + x + 1
            sage: f.__add__(pi).display()
            (x, y) |--> pi + y^2 + x
            sage: f.__add__(x).display()
            (x, y) |--> y^2 + 2*x
            sage: f.__add__(-f).display()
            (x, y) |--> 0
            sage: f.__radd__(g) == g.__add__(f)
            True

        """
        if isinstance(other, CoordFunctionSymb):
            if other._chart != self._chart:
                raise ValueError("two coordinate functions not defined on " +
                                 "the same chart cannot be added")
            res = self._simplify(self._express + other._express)
        elif isinstance(other, (int, RingElement)):
            res = self._simplify(self._express + SR(other))
        else:
            # addition to a numerical coord. function shall fall into this case
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

        TESTS::

            sage: M = TopManifold(2, 'M')
            sage: X.<x,y> = M.chart()
            sage: f = X.function(x+y^2)
            sage: g = X.function(x+1)
            sage: s = f.__sub__(g); s.display()
            (x, y) |--> y^2 - 1
            sage: type(s)
            <class 'sage.manifolds.coord_func_symb.CoordFunctionSymb'>
            sage: f.__sub__(0).display()
            (x, y) |--> y^2 + x
            sage: f.__sub__(X.zero_function()).display()
            (x, y) |--> y^2 + x
            sage: f.__sub__(1).display()
            (x, y) |--> y^2 + x - 1
            sage: f.__sub__(x).display()
            (x, y) |--> y^2
            sage: f.__sub__(pi).display()
            (x, y) |--> -pi + y^2 + x
            sage: f.__sub__(f).display()
            (x, y) |--> 0
            sage: f.__sub__(g) == - (g.__sub__(f))
            True
            sage: f.__rsub__(g) == g.__sub__(f)
            True

        """
        if isinstance(other, CoordFunctionSymb):
            if other._chart != self._chart:
                raise ValueError("two coordinate functions not defined on " +
                                 "the same chart cannot be subtracted")
            res = self._simplify(self._express - other._express)
        elif isinstance(other, (int, RingElement)):
            res = self._simplify(self._express - SR(other))
        else:
            # subtraction of a numerical coord. function shall fall into this
            # case
            return other.__rsub__(self)
        if res == 0:
            return self._chart._zero_function
        else:
            return CoordFunctionSymb(self._chart, res)

    def __mul__(self, other):
        r"""
        Multiplication  operator.

        INPUT:

        - ``other`` -- another instance of :class:`CoordFunction` or a value

        OUTPUT:

        - coordinate function resulting from the multiplication of ``self`` by
          ``other``

        TESTS::

            sage: M = TopManifold(2, 'M')
            sage: X.<x,y> = M.chart()
            sage: f = X.function(x+y)
            sage: g = X.function(x-y)
            sage: s = f.__mul__(g); s.display()
            (x, y) |--> x^2 - y^2
            sage: type(s)
            <class 'sage.manifolds.coord_func_symb.CoordFunctionSymb'>
            sage: f.__mul__(0).display()
            (x, y) |--> 0
            sage: f.__mul__(X.zero_function()).display()
            (x, y) |--> 0
            sage: f.__mul__(2).display()
            (x, y) |--> 2*x + 2*y
            sage: f.__mul__(pi).display()
            (x, y) |--> pi*x + pi*y
            sage: f.__mul__(x).display()
            (x, y) |--> x^2 + x*y
            sage: f.__mul__(1/f).display()
            (x, y) |--> 1
            sage: f.__rmul__(g) == g.__mul__(f)
            True

        """
        if isinstance(other, CoordFunctionSymb):
            if other._chart != self._chart:
                raise ValueError("two coordinate functions not defined on " +
                                 "the same chart cannot be multiplied")
            res = self._simplify(self._express * other._express)
        elif isinstance(other, (int, RingElement)):
            res = self._simplify(self._express * SR(other))
        else:
            # multiplication by a numerical coord. function shall fall into
            # this case
            return other.__rmul__(self)
        if res == 0:
            return self._chart._zero_function
        else:
            return CoordFunctionSymb(self._chart, res)

    def __div__(self, other):
        r"""
        Division  operator.

        INPUT:

        - ``other`` -- another instance of :class:`CoordFunction` or a value

        OUTPUT:

        - coordinate function resulting from the division of ``self`` by
          ``other``

        TESTS::

            sage: M = TopManifold(2, 'M')
            sage: X.<x,y> = M.chart()
            sage: f = X.function(x+y)
            sage: g = X.function(1+x^2+y^2)
            sage: s = f.__div__(g); s.display()
            (x, y) |--> (x + y)/(x^2 + y^2 + 1)
            sage: type(s)
            <class 'sage.manifolds.coord_func_symb.CoordFunctionSymb'>
            sage: f.__div__(X.zero_function())
            Traceback (most recent call last):
            ...
            ZeroDivisionError: Division of a coordinate function by zero
            sage: f.__div__(1).display()
            (x, y) |--> x + y
            sage: f.__div__(2).display()
            (x, y) |--> 1/2*x + 1/2*y
            sage: f.__div__(pi).display()
            (x, y) |--> (x + y)/pi
            sage: f.__div__(1+x^2).display()
            (x, y) |--> (x + y)/(x^2 + 1)
            sage: f.__div__(1+x^2).display()
            (x, y) |--> (x + y)/(x^2 + 1)
            sage: f.__div__(g) == ~(g.__div__(f))
            True
            sage: f.__rdiv__(g) == g.__div__(f)
            True

        """
        if isinstance(other, CoordFunctionSymb):
            if other._chart != self._chart:
                raise ValueError("two coordinate functions not defined on " +
                                 "the same chart cannot be divided")
            if other._express.is_zero():
                raise ZeroDivisionError("Division of a coordinate function " +
                                        "by zero")
            res = self._simplify(self._express / other._express)
        elif isinstance(other, (int, RingElement)):
            res = self._simplify(self._express / SR(other))
        else:
            # division by a numerical coord. function shall fall into this
            # case
            return other.__rdiv__(self)
        if res == 0:
            return self._chart._zero_function
        else:
            return CoordFunctionSymb(self._chart, res)

    def exp(self):
        r"""
        Exponential of the coordinate function.

        OUTPUT:

        - coordinate function `\exp(f)`, where `f` is the current coordinate
          function.

        EXAMPLES::

            sage: M = TopManifold(2, 'M')
            sage: X.<x,y> = M.chart()
            sage: f = X.function(x+y)
            sage: f.exp()
            e^(x + y)
            sage: exp(f) # equivalent to f.exp()
            e^(x + y)
            sage: exp(f).display()
            (x, y) |--> e^(x + y)
            sage: exp(X.zero_function())
            1

        """
        return CoordFunctionSymb(self._chart,
                                 self._simplify(self._express.exp()))

    def log(self, base=None):
        r"""
        Logarithm of the coordinate function.

        INPUT:

        - ``base`` -- (default: ``None``) base of the logarithm; if None, the
          natural logarithm (i.e. logarithm to base e) is returned

        OUTPUT:

        - coordinate function `\log_a(f)`, where `f` is the current coordinate
          function and `a` is the base

        EXAMPLES::

            sage: M = TopManifold(2, 'M')
            sage: X.<x,y> = M.chart()
            sage: f = X.function(x+y)
            sage: f.log()
            log(x + y)
            sage: log(f) # equivalent to f.log()
            log(x + y)
            sage: log(f).display()
            (x, y) |--> log(x + y)
            sage: f.log(2)
            log(x + y)/log(2)
            sage: log(f, 2)
            log(x + y)/log(2)

        """
        return CoordFunctionSymb(self._chart,
                                 self._simplify(self._express.log(base)))

    def __pow__(self, exponent):
        r"""
        Power of the coordinate function.

        INPUT:

        - ``exponent`` -- the exponent

        OUTPUT:

        - coordinate function `f^a`, where `f` is the current coordinate
          function and `a` is the exponent

        EXAMPLES::

            sage: M = TopManifold(2, 'M')
            sage: X.<x,y> = M.chart()
            sage: f = X.function(x+y)
            sage: f.__pow__(3)
            x^3 + 3*x^2*y + 3*x*y^2 + y^3
            sage: f^3  # equivalent to f.__pow__(3)
            x^3 + 3*x^2*y + 3*x*y^2 + y^3
            sage: f.__pow__(3).display()
            (x, y) |--> x^3 + 3*x^2*y + 3*x*y^2 + y^3
            sage: pow(f,3).display()
            (x, y) |--> x^3 + 3*x^2*y + 3*x*y^2 + y^3
            sage: (f^3).display()
            (x, y) |--> x^3 + 3*x^2*y + 3*x*y^2 + y^3
            sage: pow(X.zero_function(), 3).display()
            (x, y) |--> 0

        """
        return CoordFunctionSymb(self._chart,
                                 self._simplify(pow(self._express, exponent)))

    def sqrt(self):
        r"""
        Square root of the coordinate function.

        OUTPUT:

        - coordinate function `\sqrt{f}`, where `f` is the current coordinate
          function.

        EXAMPLES::

            sage: M = TopManifold(2, 'M')
            sage: X.<x,y> = M.chart()
            sage: f = X.function(x+y)
            sage: f.sqrt()
            sqrt(x + y)
            sage: sqrt(f)  # equivalent to f.sqrt()
            sqrt(x + y)
            sage: sqrt(f).display()
            (x, y) |--> sqrt(x + y)
            sage: sqrt(X.zero_function()).display()
            (x, y) |--> 0

        """
        return CoordFunctionSymb(self._chart,
                                 self._simplify(self._express.sqrt()))

    def cos(self):
        r"""
        Cosine of the coordinate function.

        OUTPUT:

        - coordinate function `\cos(f)`, where `f` is the current coordinate
          function.

        EXAMPLES::

            sage: M = TopManifold(2, 'M')
            sage: X.<x,y> = M.chart()
            sage: f = X.function(x*y)
            sage: f.cos()
            cos(x*y)
            sage: cos(f)  # equivalent to f.cos()
            cos(x*y)
            sage: cos(f).display()
            (x, y) |--> cos(x*y)
            sage: cos(X.zero_function()).display()
            (x, y) |--> 1

        """
        return CoordFunctionSymb(self._chart,
                                 self._simplify(self._express.cos()))

    def sin(self):
        r"""
        Sine of the coordinate function.

        OUTPUT:

        - coordinate function `\sin(f)`, where `f` is the current coordinate
          function.

        EXAMPLES::

            sage: M = TopManifold(2, 'M')
            sage: X.<x,y> = M.chart()
            sage: f = X.function(x*y)
            sage: f.sin()
            sin(x*y)
            sage: sin(f)  # equivalent to f.sin()
            sin(x*y)
            sage: sin(f).display()
            (x, y) |--> sin(x*y)
            sage: sin(X.zero_function()) == X.zero_function()
            True

        """
        return CoordFunctionSymb(self._chart,
                                 self._simplify(self._express.sin()))

    def tan(self):
        r"""
        Tangent of the coordinate function.

        OUTPUT:

        - coordinate function `\tan(f)`, where `f` is the current coordinate
          function.

        EXAMPLES::

            sage: M = TopManifold(2, 'M')
            sage: X.<x,y> = M.chart()
            sage: f = X.function(x*y)
            sage: f.tan()
            sin(x*y)/cos(x*y)
            sage: tan(f)  # equivalent to f.tan()
            sin(x*y)/cos(x*y)
            sage: tan(f).display()
            (x, y) |--> sin(x*y)/cos(x*y)
            sage: tan(X.zero_function()) == X.zero_function()
            True

        """
        return CoordFunctionSymb(self._chart,
                                 self._simplify(self._express.tan()))

    def arccos(self):
        r"""
        Arc cosine of the coordinate function.

        OUTPUT:

        - coordinate function `\arccos(f)`, where `f` is the current coordinate
          function.

        EXAMPLES::

            sage: M = TopManifold(2, 'M')
            sage: X.<x,y> = M.chart()
            sage: f = X.function(x*y)
            sage: f.arccos()
            arccos(x*y)
            sage: arccos(f)  # equivalent to f.arccos()
            arccos(x*y)
            sage: acos(f)  # equivalent to f.arccos()
            arccos(x*y)
            sage: arccos(f).display()
            (x, y) |--> arccos(x*y)
            sage: arccos(X.zero_function()).display()
            (x, y) |--> 1/2*pi

        """
        return CoordFunctionSymb(self._chart,
                                 self._simplify(self._express.arccos()))

    def arcsin(self):
        r"""
        Arc sine of the coordinate function.

        OUTPUT:

        - coordinate function `\arcsin(f)`, where `f` is the current coordinate
          function.

        EXAMPLES::

            sage: M = TopManifold(2, 'M')
            sage: X.<x,y> = M.chart()
            sage: f = X.function(x*y)
            sage: f.arcsin()
            arcsin(x*y)
            sage: arcsin(f)  # equivalent to f.arcsin()
            arcsin(x*y)
            sage: asin(f)  # equivalent to f.arcsin()
            arcsin(x*y)
            sage: arcsin(f).display()
            (x, y) |--> arcsin(x*y)
            sage: arcsin(X.zero_function()) == X.zero_function()
            True

        """
        return CoordFunctionSymb(self._chart,
                                 self._simplify(self._express.arcsin()))

    def arctan(self):
        r"""
        Arc tangent of the coordinate function.

        OUTPUT:

        - coordinate function `\arctan(f)`, where `f` is the current coordinate
          function.

        EXAMPLES::

            sage: M = TopManifold(2, 'M')
            sage: X.<x,y> = M.chart()
            sage: f = X.function(x*y)
            sage: f.arctan()
            arctan(x*y)
            sage: arctan(f)  # equivalent to f.arctan()
            arctan(x*y)
            sage: atan(f)  # equivalent to f.arctan()
            arctan(x*y)
            sage: arctan(f).display()
            (x, y) |--> arctan(x*y)
            sage: arctan(X.zero_function()) == X.zero_function()
            True

        """
        return CoordFunctionSymb(self._chart,
                                 self._simplify(self._express.arctan()))

    def cosh(self):
        r"""
        Hyperbolic cosine of the coordinate function.

        OUTPUT:

        - coordinate function `\cosh(f)`, where `f` is the current coordinate
          function.

        EXAMPLES::

            sage: M = TopManifold(2, 'M')
            sage: X.<x,y> = M.chart()
            sage: f = X.function(x*y)
            sage: f.cosh()
            cosh(x*y)
            sage: cosh(f)  # equivalent to f.cosh()
            cosh(x*y)
            sage: cosh(f).display()
            (x, y) |--> cosh(x*y)
            sage: cosh(X.zero_function()).display()
            (x, y) |--> 1

        """
        return CoordFunctionSymb(self._chart,
                                 self._simplify(self._express.cosh()))

    def sinh(self):
        r"""
        Hyperbolic sine of the coordinate function.

        OUTPUT:

        - coordinate function `\sinh(f)`, where `f` is the current coordinate
          function.

        EXAMPLES::

            sage: M = TopManifold(2, 'M')
            sage: X.<x,y> = M.chart()
            sage: f = X.function(x*y)
            sage: f.sinh()
            sinh(x*y)
            sage: sinh(f)  # equivalent to f.sinh()
            sinh(x*y)
            sage: sinh(f).display()
            (x, y) |--> sinh(x*y)
            sage: sinh(X.zero_function()) == X.zero_function()
            True

        """
        return CoordFunctionSymb(self._chart,
                                 self._simplify(self._express.sinh()))

    def tanh(self):
        r"""
        Hyperbolic tangent of the coordinate function.

        OUTPUT:

        - coordinate function `\tanh(f)`, where `f` is the current coordinate
          function.

        EXAMPLES::

            sage: M = TopManifold(2, 'M')
            sage: X.<x,y> = M.chart()
            sage: f = X.function(x*y)
            sage: f.tanh()
            sinh(x*y)/cosh(x*y)
            sage: tanh(f)  # equivalent to f.tanh()
            sinh(x*y)/cosh(x*y)
            sage: tanh(f).display()
            (x, y) |--> sinh(x*y)/cosh(x*y)
            sage: tanh(X.zero_function()) == X.zero_function()
            True

        """
        return CoordFunctionSymb(self._chart,
                                 self._simplify(self._express.tanh()))

    def arccosh(self):
        r"""
        Inverse hyperbolic cosine of the coordinate function.

        OUTPUT:

        - coordinate function `\mathrm{arcosh}(f)`, where `f` is the current
          coordinate function.

        EXAMPLES::

            sage: M = TopManifold(2, 'M')
            sage: X.<x,y> = M.chart()
            sage: f = X.function(x*y)
            sage: f.arccosh()
            arccosh(x*y)
            sage: arccosh(f)  # equivalent to f.arccosh()
            arccosh(x*y)
            sage: acosh(f)  # equivalent to f.arccosh()
            arccosh(x*y)
            sage: arccosh(f).display()
            (x, y) |--> arccosh(x*y)
            sage: arccosh(X.function(1)) == X.zero_function()
            True

        """
        return CoordFunctionSymb(self._chart,
                                 self._simplify(self._express.arccosh()))

    def arcsinh(self):
        r"""
        Inverse hyperbolic sine of the coordinate function.

        OUTPUT:

        - coordinate function `\mathrm{arsinh}(f)`, where `f` is the current
          coordinate function.

        EXAMPLES::

            sage: M = TopManifold(2, 'M')
            sage: X.<x,y> = M.chart()
            sage: f = X.function(x*y)
            sage: f.arcsinh()
            arcsinh(x*y)
            sage: arcsinh(f)  # equivalent to f.arcsinh()
            arcsinh(x*y)
            sage: asinh(f)  # equivalent to f.arcsinh()
            arcsinh(x*y)
            sage: arcsinh(f).display()
            (x, y) |--> arcsinh(x*y)
            sage: arcsinh(X.zero_function()) == X.zero_function()
            True

        """
        return CoordFunctionSymb(self._chart,
                                self._simplify(self._express.arcsinh()))

    def arctanh(self):
        r"""
        Inverse hyperbolic tangent of the coordinate function.

        OUTPUT:

        - coordinate function `\mathrm{artanh}(f)`, where `f` is the current
          coordinate function.

        EXAMPLES::

            sage: M = TopManifold(2, 'M')
            sage: X.<x,y> = M.chart()
            sage: f = X.function(x*y)
            sage: f.arctanh()
            arctanh(x*y)
            sage: arctanh(f)  # equivalent to f.arctanh()
            arctanh(x*y)
            sage: atanh(f)  # equivalent to f.arctanh()
            arctanh(x*y)
            sage: arctanh(f).display()
            (x, y) |--> arctanh(x*y)
            sage: arctanh(X.zero_function()) == X.zero_function()
            True

        """
        return CoordFunctionSymb(self._chart,
                                 self._simplify(self._express.arctanh()))

    # -------------------------------------
    # Methods specific to CoordFunctionSymb
    # -------------------------------------

    def _del_derived(self):
        r"""
        Delete the derived quantities

        TESTS::

            sage: M = TopManifold(2, 'M')
            sage: X.<x,y> = M.chart()
            sage: f = X.function(cos(x*y))
            sage: f._der
            sage: f.diff(x)
            -y*sin(x*y)
            sage: f._der
            [-y*sin(x*y), -x*sin(x*y)]
            sage: f._del_derived()
            sage: f._der

        """
        self._der = None  # reset of the partial derivatives

    def simplify(self):
        r"""
        Simplify the coordinate expression of the function.

        For details about the employed chain of simplifications, see
        :func:`~sage.manifolds.utilities.simplify_chain_real` for coordinate
        functions on real manifolds and
        :func:`~sage.manifolds.utilities.simplify_chain_generic` for the
        generic case.

        OUTPUT:

        - the coordinate function, with its expression simplified.

        EXAMPLES:

        Simplification of a 2-dimension coordinate function::

            sage: M = TopManifold(2, 'M')
            sage: X.<x,y> = M.chart()
            sage: f = X.function(cos(x)^2+sin(x)^2 + sqrt(x^2))
            sage: f.display()
            (x, y) |--> cos(x)^2 + sin(x)^2 + sqrt(x^2)
            sage: f.simplify()
            abs(x) + 1

        The method ``simplify()`` has changed the expression of ``f``::

            sage: f.display()
            (x, y) |--> abs(x) + 1

        Another example::

            sage: f = X.function((x^2-1)/(x+1))
            sage: f
            (x^2 - 1)/(x + 1)
            sage: f.simplify()
            x - 1

        Examples taking into account the declared range of a coordinate::

            sage: M =  TopManifold(2, 'M_1')
            sage: X.<x,y> = M.chart('x:(0,+oo) y')
            sage: f = X.function(sqrt(x^2))
            sage: f
            sqrt(x^2)
            sage: f.simplify()
            x

        ::

            sage: forget()  # to clear the previous assumption on x
            sage: M =  TopManifold(2, 'M_2')
            sage: X.<x,y> = M.chart('x:(-oo,0) y')
            sage: f = X.function(sqrt(x^2))
            sage: f
            sqrt(x^2)
            sage: f.simplify()
            -x

        """
        self._express = self._simplify(self._express)
        self._del_derived()
        return self

    def factor(self):
        r"""
        Factorize the coordinate expression of the function.

        OUTPUT:

        - the coordinate function, with its expression factorized.

        EXAMPLES:

        Factorization of a 2-dimensional coordinate function::

            sage: M = TopManifold(2, 'M')
            sage: X.<x,y> = M.chart()
            sage: f = X.function(x^2 + 2*x*y + y^2)
            sage: f.display()
            (x, y) |--> x^2 + 2*x*y + y^2
            sage: f.factor()
            (x + y)^2

        The method ``factor()`` has changed the expression of ``f``::

            sage: f.display()
            (x, y) |--> (x + y)^2

        """
        self._express = self._express.factor()
        self._del_derived()
        return self

    def expand(self):
        r"""
        Expand the coordinate expression of the function.

        OUTPUT:

        - the coordinate function, with its expression expanded.

        EXAMPLES:

        Expanding a 2-dimensional coordinate function::

            sage: M = TopManifold(2, 'M')
            sage: X.<x,y> = M.chart()
            sage: f = X.function((x - y)^2)
            sage: f.display()
            (x, y) |--> (x - y)^2
            sage: f.expand()
            x^2 - 2*x*y + y^2

        The method ``expand()`` has changed the expression of ``f``::

            sage: f.display()
            (x, y) |--> x^2 - 2*x*y + y^2

        """
        self._express = self._express.expand()
        self._del_derived()
        return self

    def collect(self, s):
        r"""
        Collect the coefficients of `s` in the expression of the coordinate
        function into a group.

        INPUT:

        - ``s`` -- the symbol whose coefficients will be collected.

        OUTPUT:

        - the coordinate function, with the coefficients of ``s`` grouped in
          its expression.

        EXAMPLES:

        Action on a 2-dimensional coordinate function::

            sage: M = TopManifold(2, 'M')
            sage: X.<x,y> = M.chart()
            sage: f = X.function(x^2*y + x*y + (x*y)^2)
            sage: f.display()
            (x, y) |--> x^2*y^2 + x^2*y + x*y
            sage: f.collect(y)
            x^2*y^2 + (x^2 + x)*y

        The method ``collect()`` has changed the expression of ``f``::

            sage: f.display()
            (x, y) |--> x^2*y^2 + (x^2 + x)*y

        """
        self._express = self._express.collect(s)
        self._del_derived()
        return self

    def collect_common_factors(self):
        r"""
        Collect common factors in the expression of the coordinate
        function.

        This method does not perform a full factorization but only looks
        for factors which are already explicitly present.

        OUTPUT:

        - the coordinate function, with the common factors collected in
          its expression.

        EXAMPLES:

        Action on a 2-dimensional coordinate function::

            sage: M = TopManifold(2, 'M')
            sage: X.<x,y> = M.chart()
            sage: f = X.function(x/(x^2*y + x*y))
            sage: f.display()
            (x, y) |--> x/(x^2*y + x*y)
            sage: f.collect_common_factors()
            1/((x + 1)*y)


        The method ``collect_common_factors()`` has changed the expression
        of ``f``::

            sage: f.display()
            (x, y) |--> 1/((x + 1)*y)

        """
        self._express = self._express.collect_common_factors()
        self._del_derived()
        return self
