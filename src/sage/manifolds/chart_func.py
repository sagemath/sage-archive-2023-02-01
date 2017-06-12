r"""
Chart Functions


AUTHORS:

- Eric Gourgoulhon, Marco Mancini (2017) : initial version

"""
#*****************************************************************************
#       Copyright (C) 2017 Eric Gourgoulhon <eric.gourgoulhon@obspm.fr>
#       Copyright (C) 2017 Marco Mancini <marco.mancini@obspm.fr>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#  as published by the Free Software Foundation; either version 2 of
#  the License, or (at your option) any later version.
#                  http://www.gnu.org/licenses/
#*****************************************************************************
from sage.structure.element import AlgebraElement
from sage.structure.parent import Parent
from sage.structure.unique_representation import UniqueRepresentation
from sage.categories.commutative_algebras import CommutativeAlgebras
from sage.manifolds.utilities import ExpressionNice
from sage.misc.cachefunc import cached_method

from sage.manifolds.calculus_method import CalculusMethod
from sage.symbolic.ring import SR
import sympy

from sage.structure.sage_object import SageObject


class ChartFunction(AlgebraElement):
    r"""
    Class containing the dictionary of ChartFunction

    INPUT:

    - ``parent`` -- the algebra of chart functions on a given chart

    """


    def __init__(self, parent, expression=None, calc_method=None ):
        r"""
        Initialize ``self``.

        TEST::

            sage: M = Manifold(2, 'M', structure='topological')
            sage: X.<x,y> = M.chart()
            sage: from sage.manifolds.chart_func import ChartFunction
            sage: f = ChartFunction(X.function_ring(),x^3+y)
            sage: f
            x^3 + y


            sage: M = Manifold(2, 'M', structure='topological')
            sage: X.<x,y> = M.chart(calc_method='sympy')
            sage: from sage.manifolds.chart_func import ChartFunction
            sage: g = ChartFunction(X.function_ring(),x**3+y)
            sage: g
            x**3 + y

            sage: M = Manifold(2, 'M', structure='topological')
            sage: X.<x,y> = M.chart(calc_method='sympy')
            sage: f = X.function(x**3+y)
            sage: f
            x**3 + y



        """
        AlgebraElement.__init__(self, parent)

        self._nc = len(parent._chart[:])

        self._express = { }

        # set the calculation method managing
        self._calc_method = parent._chart._calc_method
        self._simplify = self._calc_method._simplify


        if calc_method is None :
            calc_method = self._calc_method._current

        if expression is not None :
            self._express[calc_method] = self._calc_method._tranf[calc_method](expression)

        # Derived quantities:
        self._der = None  # list of partial derivatives (to be set by diff()
                          # and unset by del_derived())


    def expr(self,method=None):
        r"""
        Get expression from ``self`` with a particular method.

        TEST::


            sage: M = Manifold(2, 'M', structure='topological')
            sage: X.<x,y> = M.chart()
            sage: f = X.function(x^2+y)
            sage: f.display()
            (x, y) |--> x^2 + y
            sage: a = f.expr('sympy')
            sage: type(a)
            <class 'sympy.core.add.Add'>
            sage: f._express
            {'SR': x^2 + y, 'sympy': x**2 + y}

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

        Note that for substituting the value of a coordinate, the function
        call can be used as well::

            sage: f(x,3)
            3*a*x
            sage: bool( f(x,3) == f.expr().subs(y=3) )
            True


        """
        if method is None : method = self._calc_method._current
        if method in self._express :
            return self._express[method]
        else :
            for kk,vv in self._express.items():
                self._express[method] = self._calc_method._tranf[method](vv)
                return self._express[method]


    def set_expr(self,method,expression):
        r"""
        String representation of ``self``.

        TESTS::

            sage: M = Manifold(2, 'M', structure='topological')
            sage: X.<x,y> = M.chart()
            sage: f = X.function(1+x^2)
            sage: f._repr_()
            'x^2 + 1'
            sage: f.set_expr('sympy','x**2+1')
            sage: f  # indirect doctest
            x^2 + 1

            sage: g = X.function(1+x^3)
            sage: g._repr_()
            'x^3 + 1'
            sage: g.set_expr('sympy','x**2+y')
            Traceback (most recent call last):
            ...
            ValueError: Expressions are not equal


        """
        for kk,vv in self._express.items():
            # print(expression,self._calc_method._tranf[method](vv),vv,kk)
            if not bool(self._calc_method._tranf[method](expression) == self._calc_method._tranf[method](vv)):
                raise ValueError("Expressions are not equal")

        self._express[method] = expression




    def _repr_(self):
        r"""
        String representation of ``self``.

        TESTS::

            sage: M = Manifold(2, 'M', structure='topological')
            sage: X.<x,y> = M.chart()
            sage: f = X.function(1+x*y)
            sage: f._repr_()
            'x*y + 1'
            sage: repr(f)  # indirect doctest
            'x*y + 1'
            sage: f  # indirect doctest
            x*y + 1

        """

        return str(self.expr(self._calc_method._current))


    def _latex_(self):
        r"""
        LaTeX representation of ``self``.

        TESTS::

            sage: M = Manifold(2, 'M', structure='topological')
            sage: X.<x,y> = M.chart()
            sage: f = X.function(cos(x*y/2))
            sage: f._latex_()
            \cos\left(\frac{1}{2} \, x y\right)
            sage: latex(f)  # indirect doctest
            \cos\left(\frac{1}{2} \, x y\right)

        """
        from sage.misc.latex import latex
        return latex(self._express[self._calc_method._current])


    def display(self):
        r"""
        Display ``self`` in arrow notation.

        The output is either text-formatted (console mode) or
        LaTeX-formatted (notebook mode).

        EXAMPLES:

        Coordinate function on a 2-dimensional manifold::

            sage: M = Manifold(2, 'M', structure='topological')
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
        from sage.misc.latex import latex
        if self._calc_method._current in self._express:
            method = self._calc_method._current
        elif self._calc_method._default in self._express:
            method = self._calc_method._default
        else:
            method = self._express.keys()[0]
        resu_txt = str(self.parent()._chart[:]) + ' |--> ' + \
                   str(ExpressionNice(self._express[method]))
        resu_latex = latex(self.parent()._chart[:]) + r' \mapsto' + \
                     latex(ExpressionNice(self._express[method]))
        return FormattedExpansion(resu_txt, resu_latex)

    disp = display


    def __call__(self, *coords, **options):
        r"""
        Compute the value of the function at specified coordinates.

        INPUT:

        - ``*coords`` -- list of coordinates `(x^1, \ldots, x^n)`,
          where the function `f` is to be evaluated
        - ``**options`` -- allows to pass ``simplify=False`` to disable the
          call of the simplification chain on the result

        OUTPUT:

        - the value `f(x^1, \ldots, x^n)`, where `f` is the current
          coordinate function

        TESTS::

            sage: M = Manifold(2, 'M', structure='topological')
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



        substitutions = dict(zip(self.parent()._chart._xx, coords))
        resu = self.expr('SR').subs(substitutions)
        if 'simplify' in options:
            if options['simplify']:
                return self._simplify['SR'](resu)
            else:
                return resu
        else:
            return self._simplify['SR'](resu)


    def __bool__(self):
        r"""
        Return ``True`` if ``self`` is nonzero and ``False`` otherwise.

        This method is called by :meth:`~sage.structure.element.Element.is_zero()`.

        EXAMPLES:

        Coordinate functions associated to a 2-dimensional chart::

            sage: M = Manifold(2, 'M', structure='topological')
            sage: X.<x,y> = M.chart()
            sage: f = X.function(x^2+3*y+1)
            sage: bool(f)
            True
            sage: f.is_zero()
            False
            sage: f == 0
            False
            sage: g = X.function(0)
            sage: bool(g)
            False
            sage: g.is_zero()
            True
            sage: X.set_calculus_method('sympy')
            sage: g.is_zero()
            True
            sage: g == 0
            True
            sage: X.zero_function().is_zero()
            True
            sage: X.zero_function() == 0
            True

        """
        return not self.expr(self._calc_method._current).is_zero()

    __nonzero__ = __bool__   # For Python2 compatibility

    def is_zero(self):
        r"""
        Return ``True`` if the function is zero and ``False`` otherwise.

        EXAMPLES:

        Coordinate functions associated to a 2-dimensional chart::

            sage: M = Manifold(2, 'M', structure='topological')
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
        return bool(self.expr(self._calc_method._current) == 0)


    def copy(self):
        r"""
        Return an exact copy of the object.

        OUTPUT:

        - a :class:`ChartFunctionSymb`

        EXAMPLES::

            sage: M = Manifold(2, 'M', structure='topological')
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
        resu = type(self)(self.parent())
        for kk,vv in self._express.items():
            resu._express[kk] = vv
        return resu

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

        - a :class:`ChartFunction` representing the partial
          derivative `\frac{\partial f}{\partial x^i}`

        EXAMPLES:

        Partial derivatives of a 2-dimensional coordinate function::

            sage: M = Manifold(2, 'M', structure='topological')
            sage: X.<x,y> = M.chart(calc_method='SR')
            sage: f = X.function(x^2+3*y+1); f
            x^2 + 3*y + 1
            sage: f.diff(x)
            2*x
            sage: f.diff(y)
            3

        Each partial derivatives is itself a coordinate function::

            sage: type(f.diff(x))
            <class 'sage.manifolds.chart_func.ChartFunctionRing_with_category.element_class'>

        An index can be used instead of the coordinate symbol::

            sage: f.diff(0)
            2*x
            sage: f.diff(1)
            3

        The index range depends on the convention used on the chart's domain::

            sage: M = Manifold(2, 'M', structure='topological', start_index=1)
            sage: X.<x,y> = M.chart(calc_method='sympy')
            sage: f = X.function(x**2+3*y+1)
            sage: f.diff(0)
            Traceback (most recent call last):
            ...
            ValueError: coordinate index out of range
            sage: f.diff(1)
            2*x
            sage: f.diff(2)
            3

        The same test with ``sympy``::

            sage: M = Manifold(2, 'M', structure='topological')
            sage: X.<x,y> = M.chart(calc_method='sympy')
            sage: f = X.function(x^2+3*y+1); f
            x**2 + 3*y + 1
            sage: f.diff(x)
            2*x
            sage: f.diff(y)
            3

        """
        from sage.calculus.functional import diff
        from sympy import diff as diff_sympy
        from sage.rings.integer import Integer
        if self._der is None:

            # depending on the calc_method use different derivative
            curr = self._calc_method._current
            if curr == 'SR' :
                loc_diff = diff
            elif curr == 'sympy' :
                loc_diff = diff_sympy

            # the list of partial derivatives has to be updated
            self._der = [type(self)(self.parent(),
                                    self._simplify[curr](loc_diff(self.expr(), xx)))
                         for xx in self.parent()._chart[:]]
        if isinstance(coord, (int, Integer)):
            # NB: for efficiency, we access directly to the "private" attributes
            # of other classes. A more conventional OOP writing would be
            # coordsi = coord - self.parent()._chart.domain().start_index()
            coordsi = coord - self.parent()._chart._domain._sindex
            if coordsi < 0 or coordsi >= self._nc:
                raise ValueError("coordinate index out of range")
            return self._der[coordsi]
        else:
            return self._der[self.parent()._chart[:].index(coord)]

    def __eq__(self, other):
        r"""
        Comparison (equality) operator.

        INPUT:

        - ``other`` -- a :class:`ChartFunction` or a value

        OUTPUT:

        - ``True`` if ``self`` is equal to ``other``,  or ``False`` otherwise

        TESTS:

        Coordinate functions associated to a 2-dimensional chart::

            sage: M = Manifold(2, 'M', structure='topological')
            sage: X.<x,y> = M.chart()
            sage: f = X.function(x+y^2)
            sage: g = X.function(x+y^2)
            sage: f == g
            True
            sage: f = X.function(x+y^2)
            sage: g = X.function(x+y**2,'sympy')
            sage: f._express; g._express
            {'SR': y^2 + x}
            {'sympy': x + y**2}
            sage: f == g
            True
            sage: f == 1
            False
            sage: h = X.function(1)
            sage: h == 1
            True
            sage: h == f
            False
            sage: h == 0
            False
            sage: X.function(0) == 0
            True
            sage: X.zero_function() == 0
            True

        """
        if other is self:
            return True
        if isinstance(other, ChartFunction):
            if other.parent() != self.parent():
                return False
            else:
                if self._calc_method._current in self._express.keys():
                    method = self._calc_method._current
                else :
                    method = self._express.keys()[0]

                other.expr(method)
                return bool(other._express[method] == self._express[method])
        else:
            return bool(self._express[self._calc_method._current] == other)

    def __ne__(self, other):
        r"""
        Inequality operator.

        INPUT:

        - ``other`` -- a :class:`ChartFunction`

        OUTPUT:

        - ``True`` if ``self`` is different from ``other``, ``False``
          otherwise

        TESTS::

            sage: M = Manifold(2, 'M', structure='topological')
            sage: X.<x,y> = M.chart()
            sage: f = X.function(x-y)
            sage: f != X.function(x*y)
            True
            sage: f != X.function(x)
            True
            sage: f != X.function(x-y)
            False

        """
        return not (self == other)


    def __neg__(self):
        r"""
        Unary minus operator.

        OUTPUT:

        - the opposite of ``self``

        TESTS:

        Coordinate functions associated to a 2-dimensional chart::

            sage: M = Manifold(2, 'M', structure='topological')
            sage: X.<x,y> = M.chart(calc_method='sympy')
            sage: f = X.function(x+y^2)
            sage: g = -f; g
            -x - y**2
            sage: type(g)
            <class 'sage.manifolds.chart_func.ChartFunctionRing_with_category.element_class'>
            sage: -g == f
            True

        """
        resu = type(self)(self.parent())
        for kk,vv in self._express.items():
            resu._express[kk] = -vv
        return resu

    def __invert__(self):
        r"""
        Inverse operator.

        If `f` denotes the current coordinate function and `K` the topological
        field over which the manifold is defined, the *inverse* of `f` is the
        coordinate function `1/f`, where `1` of the multiplicative identity
        of `K`.

        OUTPUT:

        - the inverse of ``self``

        TESTS:

        Coordinate functions associated to a 2-dimensional chart::

            sage: M = Manifold(2, 'M', structure='topological')
            sage: X.<x,y> = M.chart()
            sage: f = X.function(1+x^2+y^2)
            sage: g = f.__invert__(); g
            1/(x^2 + y^2 + 1)
            sage: type(g)
            <class 'sage.manifolds.chart_func.ChartFunctionRing_with_category.element_class'>
            sage: g == ~f
            True
            sage: g.__invert__() == f
            True

        The same test with ``sympy``::

            sage: M = Manifold(2, 'M', structure='topological')
            sage: X.<x,y> = M.chart(calc_method='sympy')
            sage: f = X.function(1+x^2+y^2)
            sage: g = f.__invert__(); g
            1/(x**2 + y**2 + 1)
            sage: type(g)
            <class 'sage.manifolds.chart_func.ChartFunctionRing_with_category.element_class'>
            sage: g == ~f
            True
            sage: g.__invert__() == f
            True

        """
        curr = self._calc_method._current
        return type(self)(self.parent(),
                          calc_method = curr,
                          expression = self._simplify[curr]( 1/ self.expr()))
        # NB: self._express.__invert__() would return 1/self._express
        # (cf. the code of __invert__ in src/sage/symbolic/expression.pyx)
        # Here we prefer self._calc_method(1)/self._express


    def _add_(self, other):
        r"""
        Addition operator.

        INPUT:

        - ``other`` -- a :class:`ChartFunction` or a value

        OUTPUT:

        - coordinate function resulting from the addition of ``self``
          and ``other``

        TESTS::

            sage: M = Manifold(2, 'M', structure='topological')
            sage: X.<x,y> = M.chart(calc_method='sympy')

            sage: f = X.function(x+y^2)
            sage: g = X.function(x+1)
            sage: s = f + g; s.display()
            (x, y) |--> y^2 + 2*x + 1
            sage: type(s)
            <class 'sage.manifolds.chart_func.ChartFunctionRing_with_category.element_class'>
            sage: (f + 0).display()
            (x, y) |--> y^2 + x
            sage: (f + X.zero_function()).display()
            (x, y) |--> y^2 + x
            sage: (f + 1).display()
            (x, y) |--> y^2 + x + 1
            sage: (f + pi).display()
            (x, y) |--> pi + y^2 + x
            sage: (f + x).display()
            (x, y) |--> y^2 + 2*x
            sage: (f + -f).display()
            (x, y) |--> 0

        """
        curr = self._calc_method._current
        res = self._simplify[curr](self.expr() + other.expr())
        if res == 0:
            return self.parent().zero()
        else:
            return type(self)(self.parent(), res)


    def _sub_(self, other):
        r"""
        Subtraction operator.

        INPUT:

        - ``other`` -- a :class:`ChartFunction` or a value

        OUTPUT:

        - coordinate function resulting from the subtraction of ``other``
          from ``self``

        TESTS::

            sage: M = Manifold(2, 'M', structure='topological')
            sage: X.<x,y> = M.chart()
            sage: f = X.function(x+y^2)
            sage: g = X.function(x+1)
            sage: s = f - g; s.display()
            (x, y) |--> y^2 - 1
            sage: type(s)
            <class 'sage.manifolds.chart_func.ChartFunctionRing_with_category.element_class'>
            sage: (f - 0).display()
            (x, y) |--> y^2 + x
            sage: (f - X.zero_function()).display()
            (x, y) |--> y^2 + x
            sage: (f - 1).display()
            (x, y) |--> y^2 + x - 1
            sage: (f - x).display()
            (x, y) |--> y^2
            sage: (f - pi).display()
            (x, y) |--> -pi + y^2 + x
            sage: (f - f).display()
            (x, y) |--> 0
            sage: (f - g) == -(g - f)
            True

            sage: h = X.function(2*(x+y^2),calc_method='sympy')
            sage: (h - f).display()
            (x, y) |--> y^2 + x
        """
        curr = self._calc_method._current
        res = self._simplify[curr](self.expr() - other.expr())
        if curr =='SR' and res.is_trivial_zero():
            # NB: "if res == 0" would be too
            # expensive (cf. #22859)
            return self.parent().zero()
        return type(self)(self.parent(), res)

    def _mul_(self, other):
        r"""
        Multiplication operator.

        INPUT:

        - ``other`` -- a :class:`ChartFunction` or a value

        OUTPUT:

        - coordinate function resulting from the multiplication of ``self``
          by ``other``

        TESTS::

            sage: M = Manifold(2, 'M', structure='topological')
            sage: X.<x,y> = M.chart()
            sage: f = X.function(x+y)
            sage: g = X.function(x-y)
            sage: s = f._mul_(g); s.display()
            (x, y) |--> x^2 - y^2
            sage: type(s)
            <class 'sage.manifolds.chart_func.ChartFunctionRing_with_category.element_class'>
            sage: (f * 0).display()
            (x, y) |--> 0
            sage: (f * X.zero_function()).display()
            (x, y) |--> 0
            sage: (f * (1/f)).display()
            (x, y) |--> 1

        """
        curr = self._calc_method._current
        res = self._simplify[curr](self.expr() * other.expr())

        if curr =='SR' and res.is_trivial_zero():
            # NB: "if res == 0" would be too
            # expensive (cf. #22859)
            return self.parent().zero()
        return type(self)(self.parent(), res)

    def _rmul_(self, other):
        """
        Return ``other * self``.

        EXAMPLES::

            sage: M = Manifold(2, 'M', structure='topological')
            sage: X.<x,y> = M.chart()
            sage: one = X.function_ring().one()
            sage: 2 * one
            2
            sage: f = X.function(x+y)
            sage: (f * pi).display()
            (x, y) |--> pi*x + pi*y
            sage: (x * f).display()
            (x, y) |--> x^2 + x*y
        """
        try:
            other = self._calc_method._tranf(other)
        except (TypeError, ValueError):
            return
        return type(self)(self.parent(), other * self.expr())

    def _lmul_(self, other):
        """
        Return ``self * other``.

        EXAMPLES::

            sage: M = Manifold(2, 'M', structure='topological')
            sage: X.<x,y> = M.chart()
            sage: one = X.function_ring().one()
            sage: one * 2
            2
            sage: f = X.function(x+y)
            sage: (f * 2).display()
            (x, y) |--> 2*x + 2*y
            sage: (f * pi).display()
            (x, y) |--> pi*x + pi*y
        """
        try:
            other = self._calc_method._tranf(other)
        except (TypeError, ValueError):
            return
        return type(self)(self.parent(), self.expr() * other)


    def _div_(self, other):
        r"""
        Division operator.

        INPUT:

        - ``other`` -- a :class:`CoordFunction` or a value

        OUTPUT:

        - coordinate function resulting from the division of ``self``
          by ``other``

        TESTS::

            sage: M = Manifold(2, 'M', structure='topological')
            sage: X.<x,y> = M.chart()
            sage: f = X.function(x+y)
            sage: g = X.function(1+x^2+y^2)
            sage: s = f._div_(g); s.display()
            (x, y) |--> (x + y)/(x^2 + y^2 + 1)
            sage: type(s)
            <class 'sage.manifolds.chart_func.ChartFunctionRing_with_category.element_class'>
            sage: f / X.zero_function()
            Traceback (most recent call last):
            ...
            ZeroDivisionError: division of a coordinate function by zero
            sage: (f / 1).display()
            (x, y) |--> x + y
            sage: (f / 2).display()
            (x, y) |--> 1/2*x + 1/2*y
            sage: (f / pi).display()
            (x, y) |--> (x + y)/pi
            sage: (f / (1+x^2)).display()
            (x, y) |--> (x + y)/(x^2 + 1)
            sage: (f / (1+x^2)).display()
            (x, y) |--> (x + y)/(x^2 + 1)
            sage: (f / g) == ~(g / f)
            True

        """
        if other.is_zero():
            raise ZeroDivisionError("division of a coordinate function by zero")
        curr = self._calc_method._current
        res = self._simplify[curr](self.expr() / other.expr())
        if curr =='SR' and res.is_trivial_zero():
            # NB: "if res == 0" would be too
            # expensive (cf. #22859)
            return self.parent().zero()
        return type(self)(self.parent(), res)


    def exp(self):
        r"""
        Exponential of ``self``.

        OUTPUT:

        - coordinate function `\exp(f)`, where `f` is the current
          coordinate function

        EXAMPLES::

            sage: M = Manifold(2, 'M', structure='topological')
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

        The same test with ``sympy``::

            sage: X.set_calculus_method('sympy')
            sage: f = X.function(x+y)
            sage: f.exp();
            exp(x + y)
            sage: exp(f) # equivalent to f.exp()
            exp(x + y)
            sage: exp(f).display()
            (x, y) |--> e^(x + y)
            sage: exp(X.zero_function())
            1



        """
        curr = self._calc_method._current

        if curr == 'SR':
            val = self.expr().exp()
        elif curr == 'sympy' :
            val = sympy.exp(self.expr())

        return type(self)(self.parent(), self._simplify[curr](val))

    def log(self, base=None):
        r"""
        Logarithm of ``self``.

        INPUT:

        - ``base`` -- (default: ``None``) base of the logarithm; if ``None``,
          the natural logarithm (i.e. logarithm to base `e`) is returned

        OUTPUT:

        - coordinate function `\log_a(f)`, where `f` is the current coordinate
          function and `a` is the base

        EXAMPLES::

            sage: M = Manifold(2, 'M', structure='topological')
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

        The same test with ``sympy``::

            sage: X.set_calculus_method('sympy')
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
        curr = self._calc_method._current

        if curr == 'SR':
            val = self.expr().log(base)
        elif curr == 'sympy' :
            val = sympy.log(self.expr()) if base is None else sympy.log(self.expr(),base)


        return type(self)(self.parent(), self._simplify[curr](val))

    def __pow__(self, exponent):
        r"""
        Power of ``self``.

        INPUT:

        - ``exponent`` -- the exponent

        OUTPUT:

        - coordinate function `f^a`, where `f` is the current coordinate
          function and `a` is the exponent

        EXAMPLES::

            sage: M = Manifold(2, 'M', structure='topological')
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

        The same test with ``sympy``::

            sage: X.set_calculus_method('sympy')
            sage: f = X.function(x+y)
            sage: f.__pow__(3)
            x**3 + 3*x**2*y + 3*x*y**2 + y**3
            sage: f^3  # equivalent to f.__pow__(3)
            x**3 + 3*x**2*y + 3*x*y**2 + y**3
            sage: f.__pow__(3).display()
            (x, y) |--> x^3 + 3*x^2*y + 3*x*y^2 + y^3
            sage: pow(f,3).display()
            (x, y) |--> x^3 + 3*x^2*y + 3*x*y^2 + y^3
            sage: (f^3).display()
            (x, y) |--> x^3 + 3*x^2*y + 3*x*y^2 + y^3
            sage: pow(X.zero_function(), 3).display()
            (x, y) |--> 0



        """
        curr = self._calc_method._current

        if curr == 'SR':
            val = pow(self.expr(), exponent)
        elif curr == 'sympy' :
            val = self.expr() ** exponent

        return type(self)(self.parent(), self._simplify[curr](val))


    def sqrt(self):
        r"""
        Square root of ``self``.

        OUTPUT:

        - coordinate function `\sqrt{f}`, where `f` is the current
          coordinate function

        EXAMPLES::

            sage: M = Manifold(2, 'M', structure='topological')
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

        curr = self._calc_method._current

        if curr == 'SR':
            val = self.expr().sqrt()
        elif curr == 'sympy' :
            val = sympy.sqrt(self.expr())

        return type(self)(self.parent(),self._simplify[curr](val))

    def cos(self):
        r"""
        Cosine of ``self``.

        OUTPUT:

        - coordinate function `\cos(f)`, where `f` is the current
          coordinate function

        EXAMPLES::

            sage: M = Manifold(2, 'M', structure='topological')
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

        The same tests with ``sympy`` ::

            sage: X.set_calculus_method('sympy')
            sage: f.cos()
            cos(x*y)
            sage: cos(f)  # equivalent to f.cos()
            cos(x*y)

        """
        curr = self._calc_method._current

        if curr == 'SR':
            val = self.expr().cos()
        elif curr == 'sympy' :
            val = sympy.cos(self.expr())

        return type(self)(self.parent(),self._simplify[curr](val))


    def sin(self):
        r"""
        Sine of ``self``.

        OUTPUT:

        - coordinate function `\sin(f)`, where `f` is the current
          coordinate function

        EXAMPLES::

            sage: M = Manifold(2, 'M', structure='topological')
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

        The same tests with ``sympy`` ::

            sage: X.set_calculus_method('sympy')
            sage: f.sin()
            sin(x*y)
            sage: sin(f)  # equivalent to f.sin()
            sin(x*y)

        """
        curr = self._calc_method._current

        if curr == 'SR':
            val = self.expr().sin()
        elif curr == 'sympy' :
            val = sympy.sin(self.expr())

        return type(self)(self.parent(),self._simplify[curr](val))


    def tan(self):
        r"""
        Tangent of ``self``.

        OUTPUT:

        - coordinate function `\tan(f)`, where `f` is the current
          coordinate function

        EXAMPLES::

            sage: M = Manifold(2, 'M', structure='topological')
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

        The same test with ``sympy``::

            sage: M.set_calculus_method('sympy')
            sage: g = X.function(x*y)
            sage: g.tan()
            tan(x*y)
            sage: tan(g)  # equivalent to g.tan()
            tan(x*y)
            sage: tan(g).display()
            (x, y) |--> tan(x*y)


        """
        curr = self._calc_method._current

        if curr == 'SR':
            val = self.expr().tan()
        elif curr == 'sympy' :
            val = sympy.tan(self.expr())

        return type(self)(self.parent(),self._simplify[curr](val))


    def arccos(self):
        r"""
        Arc cosine of ``self``.

        OUTPUT:

        - coordinate function `\arccos(f)`, where `f` is the current
          coordinate function

        EXAMPLES::

            sage: M = Manifold(2, 'M', structure='topological')
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

        The same test with ``sympy``::

            sage: M.set_calculus_method('sympy')
            sage: f = X.function(x*y)
            sage: f.arccos()
            acos(x*y)
            sage: arccos(f)  # equivalent to f.arccos()
            acos(x*y)
            sage: acos(f)  # equivalent to f.arccos()
            acos(x*y)
            sage: arccos(f).display()
            (x, y) |--> arccos(x*y)

        """
        curr = self._calc_method._current

        if curr == 'SR':
            val = self.expr().arccos()
        elif curr == 'sympy' :
            val = sympy.acos(self.expr())

        return type(self)(self.parent(),self._simplify[curr](val))

    def arcsin(self):
        r"""
        Arc sine of ``self``.

        OUTPUT:

        - coordinate function `\arcsin(f)`, where `f` is the current
          coordinate function

        EXAMPLES::

            sage: M = Manifold(2, 'M', structure='topological')
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

        The same tests with ``sympy`` ::

            sage: X.set_calculus_method('sympy')
            sage: f.arcsin()
            asin(x*y)
            sage: arcsin(f)  # equivalent to f.arcsin()
            asin(x*y)
            sage: asin(f)  # equivalent to f.arcsin()
            asin(x*y)

        """
        curr = self._calc_method._current

        if curr == 'SR':
            val = self.expr().arcsin()
        elif curr == 'sympy' :
            val = sympy.asin(self.expr())

        return type(self)(self.parent(),self._simplify[curr](val))

    def arctan(self):
        r"""
        Arc tangent of ``self``.

        OUTPUT:

        - coordinate function `\arctan(f)`, where `f` is the current
          coordinate function

        EXAMPLES::

            sage: M = Manifold(2, 'M', structure='topological')
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

        The same tests with ``sympy`` ::

            sage: X.set_calculus_method('sympy')
            sage: f.arctan()
            atan(x*y)
            sage: arctan(f)  # equivalent to f.arctan()
            atan(x*y)
            sage: atan(f)  # equivalent to f.arctan()
            atan(x*y)

        """
        curr = self._calc_method._current

        if curr == 'SR':
            val = self.expr().arctan()
        elif curr == 'sympy' :
            val = sympy.atan(self.expr())

        return type(self)(self.parent(),self._simplify[curr](val))


    def cosh(self):
        r"""
        Hyperbolic cosine of ``self``.

        OUTPUT:

        - coordinate function `\cosh(f)`, where `f` is the current
          coordinate function

        EXAMPLES::

            sage: M = Manifold(2, 'M', structure='topological')
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

        The same tests with ``sympy`` ::

            sage: X.set_calculus_method('sympy')
            sage: f.cosh()
            cosh(x*y)
            sage: cosh(f)  # equivalent to f.cosh()
            cosh(x*y)

        """
        curr = self._calc_method._current

        if curr == 'SR':
            val = self.expr().cosh()
        elif curr == 'sympy' :
            val = sympy.cosh(self.expr())

        return type(self)(self.parent(),self._simplify[curr](val))

    def sinh(self):
        r"""
        Hyperbolic sine of ``self``.

        OUTPUT:

        - coordinate function `\sinh(f)`, where `f` is the current
          coordinate function

        EXAMPLES::

            sage: M = Manifold(2, 'M', structure='topological')
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

        The same tests with ``sympy`` ::

            sage: X.set_calculus_method('sympy')
            sage: f.sinh()
            sinh(x*y)
            sage: sinh(f)  # equivalent to f.sinh()
            sinh(x*y)

        """
        curr = self._calc_method._current

        if curr == 'SR':
            val = self.expr().sinh()
        elif curr == 'sympy' :
            val = sympy.sinh(self.expr())

        return type(self)(self.parent(),self._simplify[curr](val))

    def tanh(self):
        r"""
        Hyperbolic tangent of ``self``.

        OUTPUT:

        - coordinate function `\tanh(f)`, where `f` is the current
          coordinate function

        EXAMPLES::

            sage: M = Manifold(2, 'M', structure='topological')
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

        The same tests with ``sympy`` ::

            sage: X.set_calculus_method('sympy')
            sage: f.tanh()
            tanh(x*y)
            sage: tanh(f)  # equivalent to f.tanh()
            tanh(x*y)

        """
        curr = self._calc_method._current

        if curr == 'SR':
            val = self.expr().tanh()
        elif curr == 'sympy' :
            val = sympy.tanh(self.expr())

        return type(self)(self.parent(),self._simplify[curr](val))

    def arccosh(self):
        r"""
        Inverse hyperbolic cosine of ``self``.

        OUTPUT:

        - coordinate function `\mathrm{arcosh}(f)`, where `f` is the current
          coordinate function

        EXAMPLES::

            sage: M = Manifold(2, 'M', structure='topological')
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

        The same tests with ``sympy`` ::

            sage: X.set_calculus_method('sympy')
            sage: f.arccosh()
            acosh(x*y)
            sage: arccosh(f)  # equivalent to f.arccosh()
            acosh(x*y)
            sage: acosh(f)  # equivalent to f.arccosh()
            acosh(x*y)


        """
        curr = self._calc_method._current

        if curr == 'SR':
            val = self.expr().arccosh()
        elif curr == 'sympy' :
            val = sympy.acosh(self.expr())

        return type(self)(self.parent(),self._simplify[curr](val))


    def arcsinh(self):
        r"""
        Inverse hyperbolic sine of ``self``.

        OUTPUT:

        - coordinate function `\mathrm{arsinh}(f)`, where `f` is the current
          coordinate function

        EXAMPLES::

            sage: M = Manifold(2, 'M', structure='topological')
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

        The same tests with ``sympy`` ::

            sage: X.set_calculus_method('sympy')
            sage: f.arcsinh()
            asinh(x*y)
            sage: arcsinh(f)  # equivalent to f.arcsinh()
            asinh(x*y)
            sage: asinh(f)  # equivalent to f.arcsinh()
            asinh(x*y)

        """
        curr = self._calc_method._current

        if curr == 'SR':
            val = self.expr().arcsinh()
        elif curr == 'sympy' :
            val = sympy.asinh(self.expr())

        return type(self)(self.parent(),self._simplify[curr](val))

    def arctanh(self):
        r"""
        Inverse hyperbolic tangent of ``self``.

        OUTPUT:

        - coordinate function `\mathrm{artanh}(f)`, where `f` is the
          current coordinate function

        EXAMPLES::

            sage: M = Manifold(2, 'M', structure='topological')
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

        The same tests with ``sympy`` ::

            sage: X.set_calculus_method('sympy')
            sage: f.arctanh()
            atanh(x*y)
            sage: arctanh(f)  # equivalent to f.arctanh()
            atanh(x*y)
            sage: atanh(f)  # equivalent to f.arctanh()
            atanh(x*y)


        """
        curr = self._calc_method._current

        if curr == 'SR':
            val = self.expr().arctanh()
        elif curr == 'sympy' :
            val = sympy.atanh(self.expr())

        return type(self)(self.parent(),self._simplify[curr](val))

    # -------------------------------------
    # Methods specific to CoordFunctionSymb
    # -------------------------------------

    def _del_derived(self):
        r"""
        Delete the derived quantities.

        TESTS::

            sage: M = Manifold(2, 'M', structure='topological')
            sage: X.<x,y> = M.chart()
            sage: f = X.function(cos(x*y))
            sage: f._der
            sage: f.diff(x)
            -y*sin(x*y)
            sage: f._der
            [-y*sin(x*y), -x*sin(x*y)]
            sage: f._del_derived()
            sage: f._der

        The same tests with ``sympy`` ::

            sage: X.set_calculus_method('sympy')
            sage: f = X.function(cos(x*y))
            sage: f._der
            sage: f.diff(x)
            -y*sin(x*y)
            sage: f._der
            [-y*sin(x*y), -x*sin(x*y)]
            sage: type(f._der[0]._express['sympy'])
            <class 'sympy.core.mul.Mul'>
            sage: f._del_derived()
            sage: f._der

        """
        self._der = None  # reset of the partial derivatives

    def simplify(self):
        r"""
        Simplify the coordinate expression of ``self``.

        For details about the employed chain of simplifications, see
        :func:`~sage.manifolds.utilities.simplify_chain_real` for coordinate
        functions on real manifolds and
        :func:`~sage.manifolds.utilities.simplify_chain_generic` for the
        generic case.

        OUTPUT:

        - ``self`` with its expression simplified

        EXAMPLES:

        Simplification of a 2-dimension coordinate function::

            sage: M = Manifold(2, 'M', structure='topological')
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


        The same tests with ``sympy`` ::

            sage: X.set_calculus_method('sympy')
            sage: f = X.function(cos(x)^2+sin(x)^2 + sqrt(x^2))
            sage: f.display()
            (x, y) |--> cos(x)^2 + sin(x)^2 + sqrt(x^2)
            sage: f.simplify()
            abs(x) + 1

            sage: f = X.function((x^2-1)/(x+1))
            sage: f
            (x**2 - 1)/(x + 1)
            sage: f.simplify()
            x - 1

        Examples taking into account the declared range of a coordinate::

            sage: M =  Manifold(2, 'M_1', structure='topological')
            sage: X.<x,y> = M.chart('x:(0,+oo) y')
            sage: f = X.function(sqrt(x^2))
            sage: f
            x
            sage: f.simplify()
            x

        ::

            sage: forget()  # to clear the previous assumption on x
            sage: M =  Manifold(2, 'M_2', structure='topological')
            sage: X.<x,y> = M.chart('x:(-oo,0) y')
            sage: f = X.function(sqrt(x^2))
            sage: f
            sqrt(x^2)
            sage: f.simplify()
            -x


        """
        curr = self._calc_method._current
        self._express[curr] = self._simplify[curr](self.expr())
        self._del_derived()
        return self

    def factor(self):
        r"""
        Factorize the coordinate expression of ``self``.

        OUTPUT:

        - ``self`` with its expression factorized

        EXAMPLES:

        Factorization of a 2-dimensional coordinate function::

            sage: M = Manifold(2, 'M', structure='topological')
            sage: X.<x,y> = M.chart()
            sage: f = X.function(x^2 + 2*x*y + y^2)
            sage: f.display()
            (x, y) |--> x^2 + 2*x*y + y^2
            sage: f.factor()
            (x + y)^2

        The method ``factor()`` has changed the expression of ``f``::

            sage: f.display()
            (x, y) |--> (x + y)^2

        The same test with ``sympy`` ::

            sage: X.set_calculus_method('sympy')
            sage: g = X.function(x^2 + 2*x*y + y^2)
            sage: g.display()
            (x, y) |--> x^2 + 2*x*y + y^2
            sage: g.factor()
            (x + y)**2


        """
        curr = self._calc_method._current
        self._express[curr] = self.expr().factor()
        self._del_derived()
        return self

    def expand(self):
        r"""
        Expand the coordinate expression of ``self``.

        OUTPUT:

        - ``self`` with its expression expanded

        EXAMPLES:

        Expanding a 2-dimensional coordinate function::

            sage: M = Manifold(2, 'M', structure='topological')
            sage: X.<x,y> = M.chart()
            sage: f = X.function((x - y)^2)
            sage: f.display()
            (x, y) |--> (x - y)^2
            sage: f.expand()
            x^2 - 2*x*y + y^2

        The method ``expand()`` has changed the expression of ``f``::

            sage: f.display()
            (x, y) |--> x^2 - 2*x*y + y^2

        The same test with ``sympy`` ::

            sage: X.set_calculus_method('sympy')
            sage: g = X.function((x - y)^2)
            sage: g.expand()
            x**2 - 2*x*y + y**2

        """
        curr = self._calc_method._current
        self._express[curr] = self.expr().expand()
        self._del_derived()
        return self

    def collect(self, s):
        r"""
        Collect the coefficients of `s` in the expression of ``self``
        into a group.

        INPUT:

        - ``s`` -- the symbol whose coefficients will be collected

        OUTPUT:

        - ``self`` with the coefficients of ``s`` grouped in
          its expression

        EXAMPLES:

        Action on a 2-dimensional coordinate function::

            sage: M = Manifold(2, 'M', structure='topological')
            sage: X.<x,y> = M.chart()
            sage: f = X.function(x^2*y + x*y + (x*y)^2)
            sage: f.display()
            (x, y) |--> x^2*y^2 + x^2*y + x*y
            sage: f.collect(y)
            x^2*y^2 + (x^2 + x)*y

        The method ``collect()`` has changed the expression of ``f``::

            sage: f.display()
            (x, y) |--> x^2*y^2 + (x^2 + x)*y

        The same test with ``sympy`` ::

            sage: X.set_calculus_method('sympy')
            sage: f = X.function(x^2*y + x*y + (x*y)^2)
            sage: f.display()
            (x, y) |--> x^2*y^2 + x^2*y + x*y
            sage: f.collect(y)
            x**2*y**2 + y*(x**2 + x)

        """
        curr = self._calc_method._current
        self._express[curr] = self.expr().collect(s)
        self._del_derived()
        return self

    def collect_common_factors(self):
        r"""
        Collect common factors in the expression of ``self``.

        This method does not perform a full factorization but only looks
        for factors which are already explicitly present.

        OUTPUT:

        - ``self`` with the common factors collected in
          its expression

        EXAMPLES:

        Action on a 2-dimensional coordinate function::

            sage: M = Manifold(2, 'M', structure='topological')
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

        The same test with ``sympy`` ::
            sage: X.set_calculus_method('sympy')
            sage: g = X.function(x/(x^2*y + x*y))
            sage: g.display()
            (x, y) |--> x/(x^2*y + x*y)
            sage: g.collect_common_factors()
            1/(y*(x + 1))


        """
        curr = self._calc_method._current
        if curr == 'sympy':
            self._express[curr] = self.expr().simplify()
        else:
            self._express[curr] = self.expr().collect_common_factors()
        self._del_derived()
        return self



class ChartFunctionRing(Parent, UniqueRepresentation):
    """
    Ring of all chart functions on a chart.

    INPUT:

    - ``chart`` -- a coordinate chart, as an instance of class
      :class:`~sage.manifolds.chart.Chart`

    EXAMPLES::

        sage: M = Manifold(2, 'M', structure='topological')
        sage: X.<x,y> = M.chart()
        sage: FR = X.function_ring(); FR
        Ring of coordinate functions on Chart (M, (x, y))
        sage: type(FR)
        <class 'sage.manifolds.chart_func.ChartFunctionRing_with_category'>
        sage: FR.category()
        Category of commutative algebras over Symbolic Ring

    """
    Element = ChartFunction

    def __init__(self,chart):
        """
        Initialize ``self``.

        EXAMPLES::

            sage: M = Manifold(2, 'M', structure='topological')
            sage: X.<x,y> = M.chart()
            sage: FR = X.function_ring()
            sage: TestSuite(FR).run()
        """

        self._chart = chart
        Parent.__init__(self, base=SR, category=CommutativeAlgebras(SR))


    def _element_constructor_(self, expression, calc_method=None):
        # if isinstance(expression, ChartFunction):
        #     return self.element_class(self, expression)
        # from sage.rings.integer_ring import ZZ
        # if expression is SR or expression is ZZ:
        #     self.element_class(self, expression, calc_method=calc_method)
        return self.element_class(self, expression, calc_method=calc_method)


    def _coerce_map_from_(self, other):
        # from sage.symbolic.ring import SR
        from sage.rings.integer_ring import ZZ
        if other is SR or other is ZZ:
            return True
        return False


    def _repr_(self):
        """
        Return a string representation of ``self``.

        EXAMPLES::

            sage: M = Manifold(2, 'M', structure='topological')
            sage: X.<x,y> = M.chart()
            sage: X.function_ring()
            Ring of coordinate functions on Chart (M, (x, y))
        """
        return "Ring of coordinate functions on {}".format(self._chart)

    def is_integral_domain(self):
        """
        Return ``False`` as ``self`` is not an integral domain.

        EXAMPLES::

            sage: M = Manifold(2, 'M', structure='topological')
            sage: X.<x,y> = M.chart()
            sage: FR = X.function_ring()
            sage: FR.is_integral_domain()
            False
            sage: FR.is_field()
            False
        """
        return False

    @cached_method
    def zero(self):
        """
        Return the constant function `0` in ``self``.

        EXAMPLES::

            sage: M = Manifold(2, 'M', structure='topological')
            sage: X.<x,y> = M.chart()
            sage: FR = X.function_ring()
            sage: FR.zero()
            0

            sage: M = Manifold(2, 'M', structure='topological', field=Qp(5))
            sage: X.<x,y> = M.chart()
            sage: X.function_ring().zero()
            0
        """
        if self._chart.manifold().base_field_type() in ['real', 'complex']:
            elt = SR.zero()
        else:
            elt = self._chart.manifold().base_field().zero()
        return self.element_class(self, elt)

    @cached_method
    def one(self):
        """
        Return the constant function `1` in ``self``.

        EXAMPLES::

            sage: M = Manifold(2, 'M', structure='topological')
            sage: X.<x,y> = M.chart()
            sage: FR = X.function_ring()
            sage: FR.one()
            1

            sage: M = Manifold(2, 'M', structure='topological', field=Qp(5))
            sage: X.<x,y> = M.chart()
            sage: X.function_ring().one()
            1 + O(5^20)
        """
        if self._chart.manifold().base_field_type() in ['real', 'complex']:
            elt = SR.one()
        else:
            elt = self._chart.manifold().base_field().one()
        return self.element_class(self, elt)

    is_field = is_integral_domain




# TODO: Make this and CoordFunction have a common ABC
class MultiCoordFunction(SageObject):
    r"""
    Coordinate function to some Cartesian power of the base field.

    If `n` and `m` are two positive integers and `(U, \varphi)` is a chart on
    a topological manifold `M` of dimension `n` over a topological field `K`,
    a *multi-coordinate function* associated to `(U, \varphi)` is a map

    .. MATH::

        \begin{array}{llcl}
        f:& V \subset K^n & \longrightarrow & K^m \\
          & (x^1, \ldots, x^n) & \longmapsto & (f_1(x^1, \ldots, x^n),
            \ldots, f_m(x^1, \ldots, x^n)),
        \end{array}

    where `V` is the codomain of `\varphi`. In other words, `f` is a
    `K^m`-valued function of the coordinates associated to the chart
    `(U, \varphi)`. Each component `f_i` (`1 \leq i \leq m`) is a coordinate
    function and is therefore stored as a
    :class:`~sage.manifolds.chart_func.CoordFunction`.

    INPUT:

    - ``chart`` -- the chart `(U, \varphi)`
    - ``expressions`` -- list (or tuple) of length `m` of elements to
      construct the coordinate functions `f_i` (`1 \leq i \leq m`); for
      symbolic coordinate functions, this must be symbolic expressions
      involving the chart coordinates, while for numerical coordinate
      functions, this must be data file names

    EXAMPLES:

    A function `f: V \subset \RR^2 \longrightarrow \RR^3`::

        sage: M = Manifold(2, 'M', structure='topological')
        sage: X.<x,y> = M.chart()
        sage: f = X.multifunction(x-y, x*y, cos(x)*exp(y)); f
        Coordinate functions (x - y, x*y, cos(x)*e^y) on the Chart (M, (x, y))
        sage: type(f)
        <class 'sage.manifolds.chart_func.MultiCoordFunction'>
        sage: f(x,y)
        (x - y, x*y, cos(x)*e^y)
        sage: latex(f)
        \left(x - y, x y, \cos\left(x\right) e^{y}\right)

    Each real-valued function `f_i` (`1 \leq i \leq m`) composing `f` can
    be accessed via the square-bracket operator, by providing `i-1` as an
    argument::

        sage: f[0]
        x - y
        sage: f[1]
        x*y
        sage: f[2]
        cos(x)*e^y

    We can give a more verbose explanation of each function::

        sage: f[0].display()
        (x, y) |--> x - y

    Each ``f[i-1]`` is an instance of
    :class:`~sage.manifolds.chart_func.ChartFunction`::

        sage: isinstance(f[0], sage.manifolds.chart_func.ChartFunction)
        True

    In the present case, ``f[i-1]`` is an instance of the subclass
    :class:`~sage.manifolds.chart_func.ChartFunction`::

        sage: type(f[0])
        <class 'sage.manifolds.chart_func.ChartFunctionRing_with_category.element_class'>

    A class :class:`MultiCoordFunction` can represent a
    real-valued function (case `m = 1`), although one should
    rather employ the class :class:`~sage.manifolds.chart_func.ChartFunction`
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
        Initialize ``self``.

        TESTS::

            sage: M = Manifold(3, 'M', structure='topological')
            sage: X.<x,y,z> = M.chart()
            sage: f = X.multifunction(x+y+z, x*y*z); f
            Coordinate functions (x + y + z, x*y*z) on the Chart (M, (x, y, z))
            sage: type(f)
            <class 'sage.manifolds.chart_func.MultiCoordFunction'>
            sage: TestSuite(f).run()

        """
        self._chart = chart
        self._nc = len(self._chart._xx)   # number of coordinates
        self._nf = len(expressions)       # number of functions
        self._functions = tuple(chart.function(express)
                                for express in expressions)

    def _repr_(self):
        r"""
        String representation of ``self``.

        TESTS::

            sage: M = Manifold(2, 'M', structure='topological')
            sage: X.<x,y> = M.chart()
            sage: f = X.multifunction(x-y, x*y, cos(x)*exp(y))
            sage: f._repr_()
            'Coordinate functions (x - y, x*y, cos(x)*e^y) on the Chart (M, (x, y))'
            sage: f
            Coordinate functions (x - y, x*y, cos(x)*e^y) on the Chart (M, (x, y))

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
            sage: latex(f)
            \left(x - y, x y, \cos\left(x\right) e^{y}\right)

        """
        from sage.misc.latex import latex
        return latex(self._functions)

    def expr(self):
        r"""
        Return a tuple of data, the item no.`i` begin sufficient to
        reconstruct the coordinate function no. `i`.

        In other words, if ``f`` is a multi-coordinate function, then
        ``f.chart().multifunction(*(f.expr()))`` results in a
        multi-coordinate function identical to ``f``.

        For a symbolic multi-coordinate function, :meth:`expr` returns the
        tuple of the symbolic expressions of the coordinate functions
        composing the object.

        EXAMPLES::

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
        Return the chart with respect to which ``self`` is defined.

        OUTPUT:

        - a :class:`~sage.manifolds.chart.Chart`

        EXAMPLES::

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

        - ``other`` -- a :class:`MultiCoordFunction`

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

        - ``other`` -- a :class:`MultiCoordFunction`

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

        - ``index`` -- index `i` of the function (`0 \leq i \leq m-1`)

        OUTPUT:

        -- a :class:`CoordFunction` representing the function

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

        - ``*coords`` -- list of coordinates where the functions are
          to be evaluated
        - ``**options`` -- allows to pass some options, e.g.,
          ``simplify=False`` to disable simplification for symbolic
          coordinate functions

        OUTPUT:

        - tuple containing the values of the `m` functions

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

    @cached_method
    def jacobian(self):
        r"""
        Return the Jacobian matrix of the system of coordinate functions.

        ``jacobian()`` is a 2-dimensional array of size `m \times n`,
        where `m` is the number of functions and `n` the number of
        coordinates, the generic element being
        `J_{ij} = \frac{\partial f_i}{\partial x^j}` with `1 \leq i \leq m`
        (row index) and `1 \leq j \leq n` (column index).

        OUTPUT:

        - Jacobian matrix as a 2-dimensional array ``J`` of
          coordinate functions with ``J[i-1][j-1]`` being
          `J_{ij} = \frac{\partial f_i}{\partial x^j}`
          for `1 \leq i \leq m` and `1 \leq j \leq n`

        EXAMPLES:

        Jacobian of a set of 3 functions of 2 coordinates::

            sage: M = Manifold(2, 'M', structure='topological')
            sage: X.<x,y> = M.chart()
            sage: f = X.multifunction(x-y, x*y, y^3*cos(x))
            sage: f.jacobian()
            [           1           -1]
            [           y            x]
            [ -y^3*sin(x) 3*y^2*cos(x)]

        Each element of the result is a
        :class:`coordinate function <CoordFunction>`::

            sage: type(f.jacobian()[2,0])
            <class 'sage.manifolds.chart_func.ChartFunctionRing_with_category.element_class'>
            sage: f.jacobian()[2,0].display()
            (x, y) |--> -y^3*sin(x)

        Test of the computation::

            sage: [[f.jacobian()[i,j] == f[i].diff(j) for j in range(2)] for i in range(3)]
            [[True, True], [True, True], [True, True]]

        Test with ``start_index = 1``::

            sage: M = Manifold(2, 'M', structure='topological', start_index=1)
            sage: X.<x,y> = M.chart()
            sage: f = X.multifunction(x-y, x*y, y^3*cos(x))
            sage: f.jacobian()
            [           1           -1]
            [           y            x]
            [ -y^3*sin(x) 3*y^2*cos(x)]
            sage: [[f.jacobian()[i,j] == f[i].diff(j+1) for j in range(2)]  # note the j+1
            ....:                                         for i in range(3)]
            [[True, True], [True, True], [True, True]]
        """
        from sage.matrix.constructor import matrix
        mat = matrix([[func.diff(coord) for coord in self._chart[:]]
                      for func in self._functions])
        mat.set_immutable()
        return mat

    @cached_method
    def jacobian_det(self):
        r"""
        Return the Jacobian determinant of the system of functions.

        The number `m` of coordinate functions must equal the number `n`
        of coordinates.

        OUTPUT:

        - a :class:`CoordFunction` representing the determinant

        EXAMPLES:

        Jacobian determinant of a set of 2 functions of 2 coordinates::

            sage: M = Manifold(2, 'M', structure='topological')
            sage: X.<x,y> = M.chart()
            sage: f = X.multifunction(x-y, x*y)
            sage: f.jacobian_det()
            x + y

        The output of :meth:`jacobian_det` is an instance of
        :class:`CoordFunction` and can therefore be called on specific
        values of the coordinates, e.g. `(x,y) = (1,2)`::

            sage: type(f.jacobian_det())
            <class 'sage.manifolds.chart_func.ChartFunctionRing_with_category.element_class'>
            sage: f.jacobian_det().display()
            (x, y) |--> x + y
            sage: f.jacobian_det()(1,2)
            3

        The result is cached::

            sage: f.jacobian_det() is f.jacobian_det()
            True

        We verify the determinant of the Jacobian::

            sage: f.jacobian_det() == det(matrix([[f[i].diff(j).expr() for j in range(2)]
            ....:                                 for i in range(2)]))
            True

        Jacobian determinant of a set of 3 functions of 3 coordinates::

            sage: M = Manifold(3, 'M', structure='topological')
            sage: X.<x,y,z> = M.chart()
            sage: f = X.multifunction(x*y+z^2, z^2*x+y^2*z, (x*y*z)^3)
            sage: f.jacobian_det().display()
            (x, y, z) |--> 6*x^3*y^5*z^3 - 3*x^4*y^3*z^4 - 12*x^2*y^4*z^5 + 6*x^3*y^2*z^6

        We verify the determinant of the Jacobian::

            sage: f.jacobian_det() == det(matrix([[f[i].diff(j).expr() for j in range(3)]
            ....:                                 for i in range(3)]))
            True

        """
        from sage.matrix.constructor import matrix
        if self._nf != self._nc:
            raise ValueError("the Jacobian matrix is not a square matrix")
        mat = self.jacobian()
        mat_expr = matrix([[mat[i,j].expr() for i in range(self._nc)]
                           for j in range(self._nc)])
        det = mat_expr.det() # the unsimplified determinant
        func = self._functions[0]
        return type(func)(func.parent(), func._simplify[self._chart._calc_method._current](det))
