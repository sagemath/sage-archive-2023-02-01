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

from sage.manifolds.symb_method import SymbMethod
from sage.symbolic.ring import SR
import sympy



class ChartFunction(AlgebraElement):
    r"""
    Class containing the dictionary of ChartFunction

    INPUT:

    - ``parent`` -- the algebra of chart functions on a given chart

    """


    def __init__(self, parent, expression=None, method=None ):
        r"""
        Initialize ``self``.

        TEST::

            sage: M = Manifold(2, 'M', structure='topological')
            sage: X.<x,y> = M.chart()
            sage: from sage.manifolds.chart_func import ChartFunction
            sage: f = ChartFunction(X.function_ring(),x+y)


        """

        AlgebraElement.__init__(self, parent)

        self._nc = len(parent._chart[:])

        self._express = { }

        if method is None :
            method = SymbMethod()._default
        else:
            SymbMethod().test(method)

        self._express[method] = expression



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
        if method is None : method = SymbMethod().current
        if method in self._express :
            return self._express[method]
        else :
            for kk,vv in self._express.items():
                self._express[method] = SymbMethod().tranf[method](vv)
                return self._express[method]


    def set_expr(self,method,expression):

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

        return str(self.expr(SymbMethod().current))


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
        return latex(self._express[SymbMethod().current])


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
        if SymbMethod().current in self._express:
            method = SymbMethod().current
        elif SymbMethod()._default in self._express:
            method = SymbMethod()._default
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
        resu = self._express[SymbMethod().current].subs(substitutions)
        # if 'simplify' in options:
        #     if options['simplify']:
        #         return self._simplify(resu)
        #     else:
        #         return resu
        # else:
        #     return self._simplify(resu)
        return resu

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
        return self._express[SymbMethod().current].is_zero()


    def copy(self):
        r"""
        Return an exact copy of the object.

        OUTPUT:

        - a :class:`CoordFunctionSymb`

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
        # resu = type(self)(self.parent(), self._express[SymbMethod().current])
        resu = type(self)(self.parent())
        for kk,vv in self._express.items():
            resu._express[kk] = vv
        return resu


    def __eq__(self, other):
        r"""
        Comparison (equality) operator.

        INPUT:

        - ``other`` -- a :class:`CoordFunction` or a value

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
            sage: g = X.function('x+y**2','sympy')
            sage: f._express; g._express
            {'SR': y^2 + x}
            {'sympy': 'x+y**2'}
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
                if SymbMethod().current in self._express.keys():
                    method = SymbMethod().current
                else :
                    method = self._express.keys()[0]

                other.expr(method)
                return bool(other._express[method] == self._express[method])

                # return bool(other._express == self._express)
        else:
                return bool(self._express[SymbMethod().current] == other)


class ChartFunctionRing(Parent, UniqueRepresentation):

    Element = ChartFunction

    def __init__(self,chart):

        self._chart = chart
        Parent.__init__(self, base=SR, category=CommutativeAlgebras(SR))


    def _element_constructor_(self, expression, method = 'SR'):
        # if isinstance(expression, ChartFunction):
        #     return self.element_class(self, expression)
        return self.element_class(self, expression, method)


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
