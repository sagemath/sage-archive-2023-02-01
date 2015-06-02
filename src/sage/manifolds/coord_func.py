r"""
Coordinate functions

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
coordinates of a chart.

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

            sage: M = TopManifold(2, 'M')
            sage: X.<x,y> = M.chart()
            sage: f = X.function(1+x+y^2)
            sage: f.chart()
            Chart (M, (x, y))
            sage: f.chart() is X
            True

        """
        return self._chart

    #~ def scalar_field(self, name=None, latex_name=None):
        #~ r"""
        #~ Construct the scalar field that has the coordinate function as
        #~ coordinate expression.
#~
        #~ The domain of the scalar field is the open subset covered by the chart
        #~ associated to the coordinate function.
#~
        #~ INPUT:
#~
        #~ - ``name`` -- (default: ``None``) name given to the scalar field
        #~ - ``latex_name`` -- (default: ``None``) LaTeX symbol to denote the scalar
          #~ field; if none is provided, the LaTeX symbol is set to ``name``
#~
        #~ OUTPUT:
#~
        #~ - instance of class
          #~ :class:`~sage.manifolds.scalarfield.ScalarField`
#~
        #~ EXAMPLES:
#~
        #~ Construction of a scalar field on a 2-dimensional manifold::
#~
            #~ sage: M = Manifold(2, 'M')
            #~ sage: c_xy.<x,y> = M.chart()
            #~ sage: fc = c_xy.function(x+2*y^3)
            #~ sage: f = fc.scalar_field() ; f
            #~ scalar field on the 2-dimensional manifold 'M'
            #~ sage: f.display()
            #~ M --> R
            #~ (x, y) |--> 2*y^3 + x
            #~ sage: f.function_chart(c_xy) is fc
            #~ True
#~
        #~ """
        #~ return self._chart._domain.scalar_field_algebra().element_class(
                     #~ self._chart._domain, coord_expression={self._chart: self},
                     #~ name=name, latex_name=latex_name)

    def __ne__(self, other):
        r"""
        Inequality operator.

        INPUT:

        - ``other`` -- another instance of :class:`CoordFunction` or a value

        OUTPUT:

        - ``True`` if ``self`` is different from ``other``,  or ``False``
          otherwise

        """
        return not self.__eq__(other)

    def __radd__(self, other):
        r"""
        Addition on the left with ``other``.

        INPUT:

        - ``other`` -- another instance of :class:`CoordFunction` or a value

        OUTPUT:

        - coordinate function resulting from the addition of ``self`` and
          ``other``

        """
        return self.__add__(other)

    def __iadd__(self, other):
        r"""
        In-place addition operator.

        INPUT:

        - ``other`` -- another instance of :class:`CoordFunction` or a value

        OUTPUT:

        - coordinate function resulting from the addition of ``self`` and
          ``other``

        """
        return self.__add__(other)

    def __rsub__(self, other):
        r"""
        Subtraction from ``other``.

        INPUT:

        - ``other`` -- another instance of :class:`CoordFunction` or a value

        OUTPUT:

        - coordinate function resulting from the subtraction of ``self`` from
          ``other``

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

        """
        return self.__sub__(other)

    def __rmul__(self, other):
        r"""
        Multiplication on the left with ``other``.

        INPUT:

        - ``other`` -- another instance of :class:`CoordFunction` or a value

        OUTPUT:

        - coordinate function resulting from the multiplication of ``other``
          by ``self``

        """
        return self.__mul__(other)

    def __imul__(self, other):
        r"""
        In-place multiplication operator.

        INPUT:

        - ``other`` -- another instance of :class:`CoordFunction` or a value

        OUTPUT:

        - coordinate function resulting from the multiplication of ``self``
          by ``other``

        """
        return self.__mul__(other)

    def __rdiv__(self, other):
        r"""
        Division of ``other`` by ``self``.

        INPUT:

        - ``other`` -- another instance of :class:`CoordFunction` or a value

        OUTPUT:

        - coordinate function resulting from the division of ``other`` by
          ``self``

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

        """
        return self.__div__(other)

    # --------------------------------------------
    # Methods to be implemented by derived classes
    # --------------------------------------------

    def _repr_(self):
        r"""
        String representation of the object.
        """
        raise NotImplementedError("CoordFunction._repr_ not implemented")

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

    disp = display

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

        """
        raise NotImplementedError("CoordFunction.diff not implemented")

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
        raise NotImplementedError("CoordFunction.__pos__ not implemented")

    def __neg__(self):
        r"""
        Unary minus operator.

        OUTPUT:

        - the opposite of the coordinate function ``self``

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

class MultiCoordFunction(SageObject):
    r"""
    Multi-coordinate function.

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
    - ``*expressions`` -- the list of the coordinate expressions of the `m`
      functions `f_i` (`1\leq i \leq m`)

    """
    def __init__(self, chart, size):
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
        from sage.misc.latex import latex
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
