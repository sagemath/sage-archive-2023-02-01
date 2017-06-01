r"""
Coordinate Functions

In the context of a topological manifold `M` over a topological field `K`,
a *coordinate function*  is a function from a chart codomain
to `K`. In other words, a coordinate function is a `K`-valued function of
the coordinates associated to some chart.

More precisely, let `(U, \varphi)` be a chart on `M`, i.e. `U` is an open
subset of `M` and `\varphi: U \to V \subset K^n` is a homeomorphism
from `U` to an open subset `V` of `K^n`. A *coordinate function associated
to the chart* `(U, \varphi)` is a function

.. MATH::

    \begin{array}{cccc}
        f:&  V\subset K^n & \longrightarrow & K \\
          &  (x^1, \ldots, x^n) & \longmapsto & f(x^1, \ldots, x^n)
    \end{array}

Coordinate functions are implemented by derived classes of the abstract base
class :class:`CoordFunction`.

The class :class:`MultiCoordFunction` implements `K^m`-valued functions of
the coordinates of a chart, with `m` a positive integer.

AUTHORS:

- Eric Gourgoulhon, Michal Bejger (2013-2015) : initial version
- Travis Scrimshaw (2016) : make :class:`CoordFunction` inheritate from
  :class:`~sage.structure.element.AlgebraElement`

"""
#*****************************************************************************
#       Copyright (C) 2015 Eric Gourgoulhon <eric.gourgoulhon@obspm.fr>
#       Copyright (C) 2015 Michal Bejger <bejger@camk.edu.pl>
#       Copyright (C) 2016 Travis Scrimshaw <tscrimsh@umn.edu>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#  as published by the Free Software Foundation; either version 2 of
#  the License, or (at your option) any later version.
#                  http://www.gnu.org/licenses/
#*****************************************************************************

from sage.misc.abstract_method import abstract_method
from sage.manifolds.coord_func import CoordFunction
from sage.manifolds.utilities import ExpressionNice

class CoordFunctionGeneric(CoordFunction):
    r"""
    Abstract base class for coordinate functions.

    If `(U, \varphi)` is a chart on a topological manifold `M` of
    dimension `n` over a topological field `K`, a *coordinate function*
    associated to `(U, \varphi)` is a map `f: V \subset K^n \to K`, where
    `V` is the codomain of `\varphi`. In other words, `f` is a `K`-valued
    function of the coordinates associated to the chart `(U, \varphi)`.

    The class :class:`CoordFunction` is an abstract one. Specific
    coordinate functions must be implemented by derived classes, like
    :class:`~sage.manifolds.coord_func_symb.CoordFunctionSymb` for
    symbolic coordinate functions.

    INPUT:

    - ``parent`` -- the algebra of coordinate functions on a given chart

    """
    def __init__(self, parent):
        r"""
        Initialize ``self``.

        TEST::

            sage: M = Manifold(2, 'M', structure='topological')
            sage: X.<x,y> = M.chart()

            sage: f = CoordFunctionGeneric(X.function_ring())

        """
        CoordFunction.__init__(self, parent)

        # symbolic expression enforced :
        self._express = None

        # Derived quantities:
        self._der = None  # list of partial derivatives (to be set by diff()
                          # and unset by del_derived())


    # ----------------------------------------------------------------
    # Methods that do not need to be re-implemented by derived classes
    # but can be overloaded
    # ----------------------------------------------------------------
    def _repr_(self):
        r"""
        String representation of ``self``.

        TESTS::

            sage: M = Manifold(2, 'M', structure='topological')
            sage: X.<x,y> = M.chart()

            sage: f = CoordFunctionGeneric(X.function_ring())
            sage: f = X.function(1+x*y)
            sage: f._repr_()
            'x*y + 1'
            sage: repr(f)  # indirect doctest
            'x*y + 1'
            sage: f  # indirect doctest
            x*y + 1

        """
        return str(self._express)

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
        return latex(self._express)


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
        resu_txt = str(self.parent()._chart[:]) + ' |--> ' + \
                   str(ExpressionNice(self._express))
        resu_latex = latex(self.parent()._chart[:]) + r' \mapsto' + \
                     latex(ExpressionNice(self._express))
        return FormattedExpansion(resu_txt, resu_latex)

    disp = display


    def expr(self):
        r"""
        Return the symbolic expression representing the image of ``self``.

        OUTPUT:

        - :class:`symbolic expression <sage.symbolic.expression.Expression>`
          involving the chart coordinates

        EXAMPLES:

        Coordinate function of a 2-dimensional manifold::

            sage: M = Manifold(2, 'M', structure='topological')
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

        Note that for substituting the value of a coordinate, the function
        call can be used as well::

            sage: f(x,3)
            3*a*x
            sage: bool( f(x,3) == f.expr().subs(y=3) )
            True

        """
        return self._express

    # --------------------------------------------
    # Methods to be implemented by derived classes
    # --------------------------------------------

    @abstract_method
    def _simplify(self,expression):
        r"""
        Simplification chain.

        TESTS:

        This method must be implemented by derived classes; it is not
        implemented here::

            sage: M = Manifold(2, 'M', structure='topological')
            sage: X.<x,y> = M.chart()
            sage: from sage.manifolds.coord_func_generic import CoordFunctionGeneric
            sage: f = CoordFunctionGeneric(X.function_ring())
            sage: f._simplify()
            Traceback (most recent call last):
            ...
            NotImplementedError: <abstract method _simplify at 0x...>
        """

    @abstract_method
    def _symb_method(self,expression):
        r"""
        Method to select the symbolic representation method.

        TESTS:

        This method must be implemented by derived classes; it is not
        implemented here::

            sage: M = Manifold(2, 'M', structure='topological')
            sage: X.<x,y> = M.chart()
            sage: from sage.manifolds.coord_func_generic import CoordFunctionGeneric
            sage: f = CoordFunctionGeneric(X.function_ring())
            sage: f._symb_method()
            Traceback (most recent call last):
            ...
            NotImplementedError: <abstract method _symb_method at 0x...>
        """
