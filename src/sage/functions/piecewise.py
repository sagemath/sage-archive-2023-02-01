r"""
Piecewise-defined Functions

This module implement piecewise functions in a single variable. See
:mod:`sage.sets.real_set` for more information about how to construct
subsets of the real line for the domains. 

EXAMPLES::

    sage: f = piecewise([((0,1), x^3), ([-1,0], -x^2)]);  f
    piecewise(x|-->x^3 on (0, 1), x|-->-x^2 on [-1, 0]; x)
    sage: 2*f
    2*piecewise(x|-->x^3 on (0, 1), x|-->-x^2 on [-1, 0]; x)
    sage: f(x=1/2)
    1/8
    sage: plot(f)    # not tested

TODO:

- Implement max/min location and values,

AUTHORS:

- David Joyner (2006-04): initial version

- David Joyner (2006-09): added __eq__, extend_by_zero_to, unextend,
  convolution, trapezoid, trapezoid_integral_approximation,
  riemann_sum, riemann_sum_integral_approximation, tangent_line fixed
  bugs in __mul__, __add__

- David Joyner (2007-03): adding Hann filter for FS, added general FS
  filter methods for computing and plotting, added options to plotting
  of FS (eg, specifying rgb values are now allowed). Fixed bug in
  documentation reported by Pablo De Napoli.

- David Joyner (2007-09): bug fixes due to behaviour of
  SymbolicArithmetic

- David Joyner (2008-04): fixed docstring bugs reported by J Morrow; added
  support for Laplace transform of functions with infinite support.

- David Joyner (2008-07): fixed a left multiplication bug reported by
  C. Boncelet (by defining __rmul__ = __mul__).

- Paul Butler (2009-01): added indefinite integration and default_variable

- Volker Braun (2013): Complete rewrite

- Ralf Stephan (2015): Rewrite of convolution() and other calculus
  functions; many doctest adaptations

TESTS::

    sage: fast_callable(f, vars=[x])(0.5)
    0.125000000000...
"""

#*****************************************************************************
#       Copyright (C) 2006 William Stein <wstein@gmail.com>
#                     2006 David Joyner <wdjoyner@gmail.com>
#                     2013 Volker Braun <vbraun.name@gmail.com>
#
#  Distributed under the terms of the GNU General Public License (GPL),
#  version 2 or any later version.  The full text of the GPL is available at:
#                  http://www.gnu.org/licenses/
#*****************************************************************************

from sage.symbolic.function import BuiltinFunction
from sage.sets.real_set import RealSet, InternalRealInterval
from sage.symbolic.ring import SR
from sage.rings.rational_field import QQ
from sage.rings.infinity import minus_infinity, infinity

class PiecewiseFunction(BuiltinFunction):

    def __init__(self):
        """
        Piecewise function

        EXAMPLES::

            sage: var('x, y')
            (x, y)
            sage: f = piecewise([((0,1), x^2*y), ([-1,0], -x*y^2)], var=x);  f
            piecewise(x|-->x^2*y on (0, 1), x|-->-x*y^2 on [-1, 0]; x)
            sage: f(1/2)
            1/4*y
            sage: f(-1/2)
            1/2*y^2
        """
        BuiltinFunction.__init__(self, "piecewise", 
                                 latex_name="piecewise",
                                 conversions=dict(), nargs=2)

    def __call__(self, function_pieces, **kwds):
        r"""
        Piecewise functions

        INPUT:
   
        - ``function_pieces`` -- a list of pairs consisting of a
          domain and a symbolic function.

        - ``var=x`` -- a symbolic variable or ``None`` (default). The
        real variable in which the function is piecewise in.

        OUTPUT:

        A piecewise-defined function. A ``ValueError`` will be raised
        if the domains of the pieces are not pairwise disjoint.
    
        EXAMPLES::
        
            sage: my_abs = piecewise([((-1, 0), -x), ([0, 1], x)], var=x);  my_abs
            piecewise(x|-->-x on (-1, 0), x|-->x on [0, 1]; x)
            sage: [ my_abs(i/5) for i in range(-4, 5)]
            [4/5, 3/5, 2/5, 1/5, 0, 1/5, 2/5, 3/5, 4/5]

        TESTS::

            sage: piecewise([([-1, 0], -x), ([0, 1], x)], var=x)
            Traceback (most recent call last):
            ...
            ValueError: domains must be pairwise disjoint

            sage: step = piecewise([((-1, 0), -1), ([0, 0], 0), ((0, 1), 1)], var=x);  step
            piecewise(x|-->-1 on (-1, 0), x|-->0 on {0}, x|-->1 on (0, 1); x)
            sage: step(-1/2), step(0), step(1/2)
            (-1, 0, 1)
        """
        #print 'pf_call', function_pieces, kwds
        from types import *
        var = kwds.pop('var', None)
        parameters = []
        domain_list = []
        for piece in function_pieces:
            domain, function = piece
            if not isinstance(domain, RealSet):
                domain = RealSet(domain)
            if domain.is_empty():
                continue
            if isinstance(function, FunctionType):
                if var is None:
                    var = SR.var('x')
                if function.func_code.co_argcount == 0:
                    function = function()
                else:
                    function = function(var)
            function = SR(function)
            if var is None and len(function.variables()) > 0:
                var = function.variables()[0]
            parameters.append((domain, function))
            domain_list.append(domain)
        if not RealSet.are_pairwise_disjoint(*domain_list):
            raise ValueError('domains must be pairwise disjoint')
        if var is None:
            var = self.default_variable()
        parameters = SR._force_pyobject(tuple(parameters), recursive=False)
        return BuiltinFunction.__call__(self, parameters, var, **kwds)

    def _print_(self, parameters, variable):
        """
        Return a string representation
        
        OUTPUT:
        
        String.

        EXAMPLES::

            sage: p = piecewise([((-2, 0), -x), ([0, 4], x)], var=x)
            sage: str(p)    # indirect doctest
            'piecewise(x|-->-x on (-2, 0), x|-->x on [0, 4]; x)'
        """
        s = 'piecewise(' 
        args = []
        for domain, func in parameters:
            args.append('{0}|-->{1} on {2}'.format(str(variable), str(func), str(domain)))
        s += ', '.join(args) + '; {0})'.format(str(variable))
        return s

    def _subs_(self, subs_map, options, parameters, x):
        """
        Callback from Pynac `subs()`

        EXAMPLES:

        If the substitution changes the piecewise variable, it must
        evaluate to a number so that we know which component we are
        on::

            sage: p = piecewise([((-2, 0), -x), ([0, 4], x)], var=x)
            sage: p.subs(x=-1)
            1
            sage: (10+p).subs(x=-1)
            11
            sage: p.subs(x=pi)
            pi

        Auxiliary variables can be substituted arbitrarily::

            sage: var('x,y')
            (x, y)
            sage: p = piecewise([((-2, 0), -x^y), ([0, 2], x-y)], var=x);  p
            piecewise(x|-->-x^y on (-2, 0), x|-->x - y on [0, 2]; x)
            sage: p.subs(y=sin(y))
            piecewise(x|-->-x^sin(y) on (-2, 0), x|-->x - sin(y) on [0, 2]; x)
        """
        point = subs_map.apply_to(x, 0)
        if point == x:
            # substitution only in auxiliary variables
            new_params = []
            for domain, func in parameters:
                new_params.append((domain, subs_map.apply_to(func, 0)))
            return piecewise(new_params, var=x)
        if ((point.is_numeric() or point.is_constant())
            and (point.is_real())):
            if hasattr(point, 'pyobject'):
                # unwrap any numeric values
                point = point.pyobject()
        else:
            raise ValueError('substituting the piecewise variable must result in real number')

        for domain, func in parameters:
            if domain.contains(point):
                return subs_map.apply_to(func, 0)
        raise ValueError('point {} is not in the domain'.format(point))

    @staticmethod
    def in_operands(ex):
        """
        Return whether a symbolic expression contains a piecewise
        function as operand

        INPUT:

        - ``ex`` -- a symbolic expression.

        OUTPUT:

        Boolean

        EXAMPLES::

            sage: f = piecewise([([0,0], sin(x)), ((0,2), cos(x))]);  f
            piecewise(x|-->sin(x) on {0}, x|-->cos(x) on (0, 2); x)
            sage: piecewise.in_operands(f)
            True
            sage: piecewise.in_operands(1+sin(f))
            True
            sage: piecewise.in_operands(1+sin(0*f))
            False
        """
        def is_piecewise(ex):
            result = ex.operator() is piecewise
            for op in ex.operands():
                result = result or is_piecewise(op)
            return result
        return is_piecewise(ex)


    @staticmethod
    def simplify(ex):
        """
        Combine piecewise operands into single piecewise function

        OUTPUT:

        A piecewise function whose operands are not piecewiese if 
        possible, that is, as long as the piecewise variable is the same.

        EXAMPLES::
            
            sage: f = piecewise([([0,0], sin(x)), ((0,2), cos(x))])
            sage: piecewise.simplify(f)
            Traceback (most recent call last):
            ...
            NotImplementedError
        """
        raise NotImplementedError


    class EvaluationMethods:

        def expression_at(cls, self, parameters, variable, point):
            """
            Return the expression defining the piecewise function at
            ``value``

            INPUT:
            
            - ``point`` -- a real number.

            OUTPUT:

            The symbolic expression defining the function value at the
            given ``point``.

            EXAMPLES::

                sage: f = piecewise([([0,0], sin(x)), ((0,2), cos(x))]);  f
                piecewise(x|-->sin(x) on {0}, x|-->cos(x) on (0, 2); x)
                sage: f.expression_at(0)
                sin(x)
                sage: f.expression_at(1)
                cos(x)
                sage: f.expression_at(2)
                Traceback (most recent call last):
                ...
                ValueError: point is not in the domain
            """
            for domain, func in parameters:
                if domain.contains(point):
                    return func
            raise ValueError('point is not in the domain')

        which_function = expression_at

        def domains(cls, self, parameters, variable):
            """
            Return the individual domains
            
            See also :meth:`~expressions`.

            OUTPUT:

            The collection of domains of the component functions as a
            tuple of :class:`~sage.sets.real_set.RealSet`.
            
            EXAMPLES::

                sage: f = piecewise([([0,0], sin(x)), ((0,2), cos(x))]);  f
                piecewise(x|-->sin(x) on {0}, x|-->cos(x) on (0, 2); x)
                sage: f.domains()
                ({0}, (0, 2))
            """
            return tuple(dom for dom, fun in parameters)

        def domain(cls, self, parameters, variable):
            """
            Return the domain
            
            OUTPUT:

            The union of the domains of the individual pieces as a
            :class:`~sage.sets.real_set.RealSet`.
            
            EXAMPLES::

                sage: f = piecewise([([0,0], sin(x)), ((0,2), cos(x))]);  f
                piecewise(x|-->sin(x) on {0}, x|-->cos(x) on (0, 2); x)
                sage: f.domain()
                [0, 2)
            """
            intervals = []
            for domain, func in parameters:
                intervals += list(domain)
            return RealSet(*intervals)

        def __len__(cls, self, parameters, variable):
            """
            Return the number of "pieces"

            OUTPUT:

            Integer.

            EXAMPLES::

                sage: f = piecewise([([0,0], sin(x)), ((0,2), cos(x))]);  f
                piecewise(x|-->sin(x) on {0}, x|-->cos(x) on (0, 2); x)
                sage: len(f)
                2
            """
            return len(parameters)

        def expressions(cls, self, parameters, variable):
            """
            Return the individual domains
            
            See also :meth:`~domains`.

            OUTPUT:

            The collection of expressions of the component functions.
            
            EXAMPLES::

                sage: f = piecewise([([0,0], sin(x)), ((0,2), cos(x))]);  f
                piecewise(x|-->sin(x) on {0}, x|-->cos(x) on (0, 2); x)
                sage: f.expressions()
                (sin(x), cos(x))
            """
            return tuple(fun for dom, fun in parameters)

        def items(cls, self, parameters, variable):
            """
            Iterate over the pieces of the piecewise function

            .. NOTE::

                You should probably use :meth:`pieces` instead, which
                offers a nicer interface.

            OUTPUT:

            This method iterates over pieces of the piecewise
            function, each represented by a pair. The first element is
            the support, and the second the function over that
            support.

            EXAMPLES::

                sage: f = piecewise([([0,0], sin(x)), ((0,2), cos(x))])
                sage: for support, function in f.items():
                ....:     print('support is {0}, function is {1}'.format(support, function))
                support is {0}, function is sin(x)
                support is (0, 2), function is cos(x)
            """
            for pair in parameters:
                yield pair

        def __call__(cls, self, parameters, variable, value=None, **kwds):
            """
            Call the piecewise function

            EXAMPLES::

                sage: f = piecewise([([0,0], sin(x)), ((0,2), cos(x))]);  f
                piecewise(x|-->sin(x) on {0}, x|-->cos(x) on (0, 2); x)
                sage: f(0)
                0
                sage: f(1)
                cos(1)
                sage: f(2)
                Traceback (most recent call last):
                ...
                ValueError: point 2 is not in the domain
            """
            self = piecewise(parameters, var=variable)
            substitution = dict()
            for k, v in kwds.items():
                substitution[SR.var(k)] = v
            if value is not None:
                substitution[variable] = value
            return self.subs(substitution)

        def _fast_float_(cls, self, *args):
            """
            Do not support the old ``fast_float``

            OUTPUT:

            This method raises ``NotImplementedError`` so that
            plotting uses the newer `fast_callable` implementation.

            EXAMPLES::
            
                sage: f = piecewise([([0,0], sin(x)), ((0,2), cos(x))])
                sage: f._fast_float_()
                Traceback (most recent call last):
                ...
                NotImplementedError
            """
            raise NotImplementedError

        def _fast_callable_(cls, self, parameters, variable, etb):
            """
            Override the ``fast_callable``

            OUTPUT:

            A :class:`~sage.ext.fast_callable.ExpressionCall`
            representing the piecewise function in the expression
            tree.
            
            EXAMPLES::

                sage: p = piecewise([((-1, 0), -x), ([0, 1], x)], var=x)
                sage: from sage.ext.fast_callable import ExpressionTreeBuilder
                sage: etb = ExpressionTreeBuilder(vars=['x'])
                sage: p._fast_callable_(etb)
                {piecewise(x|-->-x on (-1, 0), x|-->x on [0, 1]; x)}(v_0)
            """
            # print 'ev_fast_cal', parameters, variable, etb
            self = piecewise(parameters, var=variable)
            return etb.call(self, variable)

        def restriction(cls, self, parameters, variable, restricted_domain):
            """
            Restrict the domain

            INPUT:
            
            - ``restricted_domain`` -- a
              :class:`~sage.sets.real_set.RealSet` or something that
              defines one.

            OUTPUT:

            A new piecewise function obtained by restricting the domain.

            EXAMPLES::

                sage: f = piecewise([((-oo, oo), x)]);  f
                piecewise(x|-->x on (-oo, +oo); x)
                sage: f.restriction([[-1,1], [3,3]])
                piecewise(x|-->x on [-1, 1] + {3}; x)
            """
            restricted_domain = RealSet(*restricted_domain)
            new_param = []
            for domain, func in parameters:
                domain = domain.intersection(restricted_domain)
                new_param.append((domain, func))
            return piecewise(new_param, var=variable)

        def extension(cls, self, parameters, variable, extension, extension_domain=None):
            """
            Extend the function

            INPUT:

            - ``extension`` -- a symbolic expression

            - ``extension_domain`` -- a
              :class:`~sage.sets.real_set.RealSet` or ``None``
              (default). The domain of the extension. By default, the
              entire complement of the current domain.

            EXAMPLES::

                sage: f = piecewise([((-1,1), x)]);  f
                piecewise(x|-->x on (-1, 1); x)
                sage: f(3)
                Traceback (most recent call last):
                ...
                ValueError: point 3 is not in the domain

                sage: g = f.extension(0);  g
                piecewise(x|-->x on (-1, 1), x|-->0 on (-oo, -1] + [1, +oo); x)
                sage: g(3)
                0

                sage: h = f.extension(1, RealSet.unbounded_above_closed(1));  h
                piecewise(x|-->x on (-1, 1), x|-->1 on [1, +oo); x)
                sage: h(3)
                1
            """
            self = piecewise(parameters, var=variable)
            if extension_domain is None:
                extension_domain = self.domain().complement()
            ext = ((extension_domain, SR(extension)),)
            return piecewise(parameters + ext, var=variable)

        def unextend_zero(cls, self, parameters, variable):
            """
            Remove zero pieces.

            EXAMPLES::

                sage: f = piecewise([((-1,1), x)]);  f
                piecewise(x|-->x on (-1, 1); x)
                sage: g = f.extension(0);  g
                piecewise(x|-->x on (-1, 1), x|-->0 on (-oo, -1] + [1, +oo); x)
                sage: g(3)
                0
                sage: h = g.unextend_zero()
                sage: bool(h == f)
                True
            """
            result = [(domain, func) for domain,func in parameters if func<>0]
            return piecewise(result, var=variable)

        def pieces(cls, self, parameters, variable):
            """
            Return the "pieces".

            OUTPUT:

            A tuple of piecewise functions, each having only a single
            expression.

            EXAMPLES::

                sage: p = piecewise([((-1, 0), -x), ([0, 1], x)], var=x)
                sage: p.pieces()
                (piecewise(x|-->-x on (-1, 0); x), 
                 piecewise(x|-->x on [0, 1]; x))
            """
            result = []
            for domain, func in parameters:
                result.append(piecewise([(domain, func)], var=variable))
            return tuple(result)

        def end_points(cls, self, parameters, variable):
            """
            Return a list of all interval endpoints for this function.

            EXAMPLES::

                sage: f1(x) = 1
                sage: f2(x) = 1-x
                sage: f3(x) = x^2-5
                sage: f = piecewise([[(0,1),f1],[(1,2),f2],[(2,3),f3]])
                sage: f.end_points()
                [0, 1, 2, 3]
                sage: f = piecewise([([0,0], sin(x)), ((0,2), cos(x))]);  f
                piecewise(x|-->sin(x) on {0}, x|-->cos(x) on (0, 2); x)
                sage: f.end_points()
                [0, 2]
            """
            s = set()
            for domain, func in parameters:
                for interval in domain:
                    s.add(interval.lower())
                    s.add(interval.upper())
            s.discard(minus_infinity)
            s.discard(infinity)
            return sorted(s)

        def piecewise_add(cls, self, parameters, variable, other):
            """
            Return a new piecewise function with domain the union
            of the original domains and functions summed. Undefined
            intervals in the union domain get function value `0`.

            EXAMPLES::

                sage: f = piecewise([([0,1], 1), ((2,3), x)])
                sage: g = piecewise([((1/2, 2), x)])
                sage: f.piecewise_add(g).unextend_zero()
                piecewise(x|-->1 on (0, 1/2], x|-->x + 1 on (1/2, 1], x|-->x on (1, 2) + (2, 3); x)
            """
            points = ([minus_infinity] +
                      sorted(set(self.end_points() + other.end_points())) +
                      [infinity])
            domain = []
            funcs = []
            contains_lower = False
            contains_upper = False
            for i in range(len(points)-1):
                try:
                    contains_lower = (self.domain().contains(points[i]) or
                        other.domain().contains(points[i])) and not contains_upper
                    contains_upper = (self.domain().contains(points[i+1]) or
                        other.domain().contains(points[i+1]))
                    if contains_lower:
                        if contains_upper:
                            rs = RealSet.closed(points[i],points[i+1])
                        else:
                            rs = RealSet.closed_open(points[i],points[i+1])
                    else:
                        if contains_upper:
                            rs = RealSet.open_closed(points[i],points[i+1])
                        else:
                            rs = RealSet.open(points[i],points[i+1])
                    point = (points[i+1] + points[i])/2
                except ValueError:
                    if points[i] == minus_infinity and points[i+1] == infinity:
                        rs = RealSet.open(minus_infinity, infinity)
                        point = 0
                    elif points[i] == minus_infinity:
                        if contains_lower:
                            rs = RealSet.unbounded_below_closed(points[i+1])
                        else:
                            rs = RealSet.unbounded_below_open(points[i+1])
                        point = points[i+1]-1
                    elif points[i+1] == infinity:
                        if contains_upper:
                            rs = RealSet.unbounded_above_closed(points[i])
                        else:
                            rs = RealSet.unbounded_above_open(points[i])
                        point = points[i]+1
                    else:
                        raise
                try:
                    ex1 = self.expression_at(point)
                except ValueError:
                    ex1 = 0
                try:
                    ex2 = other.expression_at(point)
                except ValueError:
                    ex2 = 0
                ex = ex1 + ex2
                if i>0 and funcs[-1] == ex:
                    # extend the previous domain
                    rs += domain[-1]
                    domain[-1] = rs
                else:
                    domain += rs
                    funcs.append(ex)
            return piecewise(zip(domain, funcs))

        def integral(cls, self, parameters, variable, x=None, a=None, b=None, definite=False):
            r"""
            By default, return the indefinite integral of the function.
            If definite=True is given, returns the definite integral.

            AUTHOR:

            - Paul Butler

            EXAMPLES::

                sage: f1(x) = 1-x
                sage: f = piecewise([((0,1),1), ((1,2),f1)])
                sage: f.integral(definite=True)
                1/2

            ::

                sage: f1(x) = -1
                sage: f2(x) = 2
                sage: f = piecewise([((0,pi/2),f1), ((pi/2,pi),f2)])
                sage: f.integral(definite=True)
                1/2*pi

                sage: f1(x) = 2
                sage: f2(x) = 3 - x
                sage: f = piecewise([[(-2, 0), f1], [(0, 3), f2]])
                sage: f.integral()
                piecewise(x|-->2*x + 4 on (-2, 0), x|-->-1/2*x^2 + 3*x + 4 on (0, 3); x)

                sage: f1(y) = -1
                sage: f2(y) = y + 3
                sage: f3(y) = -y - 1
                sage: f4(y) = y^2 - 1
                sage: f5(y) = 3
                sage: f = piecewise([[[-4,-3],f1],[(-3,-2),f2],[[-2,0],f3],[(0,2),f4],[[2,3],f5]])
                sage: F = f.integral(y)
                sage: F
                piecewise(y|-->-y - 4 on [-4, -3], y|-->1/2*y^2 + 3*y + 7/2 on (-3, -2), y|-->-1/2*y^2 - y - 1/2 on [-2, 0], y|-->1/3*y^3 - y - 1/2 on (0, 2), y|-->3*y - 35/6 on [2, 3]; y)

            Ensure results are consistent with FTC::

                sage: F(-3) - F(-4)
                -1
                sage: F(-1) - F(-3)
                1
                sage: F(2) - F(0)
                2/3
                sage: f.integral(y, 0, 2)
                2/3
                sage: F(3) - F(-4)
                19/6
                sage: f.integral(y, -4, 3)
                19/6
                sage: f.integral(definite=True)
                19/6

            ::

                sage: f1(y) = (y+3)^2
                sage: f2(y) = y+3
                sage: f3(y) = 3
                sage: f = piecewise([[(-infinity, -3), f1], [(-3, 0), f2], [(0, infinity), f3]])
                sage: f.integral()
                piecewise(y|-->1/3*y^3 + 3*y^2 + 9*y + 9 on (-oo, -3), y|-->1/2*y^2 + 3*y + 9/2 on (-3, 0), y|-->3*y + 9/2 on (0, +oo); y)

            ::

                sage: f1(x) = e^(-abs(x))
                sage: f = piecewise([[(-infinity, infinity), f1]])
                sage: f.integral(definite=True)
                2
                sage: f.integral()
                piecewise(x|-->-1/2*((sgn(x) - 1)*e^(2*x) - 2*e^x*sgn(x) + sgn(x) + 1)*e^(-x) - 1 on (-oo, +oo); x)

            ::

                sage: f = piecewise([((0, 5), cos(x))])
                sage: f.integral()
                piecewise(x|-->sin(x) on (0, 5); x)


            TESTS:

            Verify that piecewise integrals of zero work (trac #10841)::

                sage: f0(x) = 0
                sage: f = piecewise([[[0,1],f0]])
                sage: f.integral(x,0,1)
                0
                sage: f = piecewise([[[0,1], 0]])
                sage: f.integral(x,0,1)
                0
                sage: f = piecewise([[[0,1], SR(0)]])
                sage: f.integral(x,0,1)
                0

            """
            if a != None and b != None:
                F = self.integral(x)
                return F(b) - F(a)

            if a != None or b != None:
                raise TypeError('only one endpoint given')

            area = 0
            new_pieces = []

            if x == None:
                x = self.default_variable()

            # The integral is computed by iterating over the pieces in order.
            # The definite integral for each piece is calculated and accumulated in `area`.
            # The indefinite integral of each piece is also calculated,
            # and the `area` before each piece is added to the piece.
            #
            # If a definite integral is requested, `area` is returned.
            # Otherwise, a piecewise function is constructed from the indefinite integrals
            # and returned.
            #
            # An exception is made if integral is called on a piecewise function
            # that starts at -infinity. In this case, we do not try to calculate the
            # definite integral of the first piece, and the value of `area` remains 0
            # after the first piece.

            from sage.symbolic.assumptions import assume, forget
            for domain, fun in parameters:
                for interval in domain:
                    start = interval.lower()
                    end = interval.upper()
                    if start == -infinity and not definite:
                        fun_integrated = fun.integral(x, end, x)
                    else:
                        try:
                            assume(start < x)
                        except ValueError: # Assumption is redundant
                            pass
                        fun_integrated = fun.integral(x, start, x) + area
                        forget(start < x)
                        if definite or end != infinity:
                            area += fun.integral(x, start, end)
                    new_pieces.append([interval, SR(fun_integrated).function(x)])

            if definite:
                return SR(area)
            else:
                return piecewise(new_pieces)

        def critical_points(cls, self, parameters, variable):
            """
            Return the critical points of this piecewise function.

            EXAMPLES::

                sage: R.<x> = QQ[]
                sage: f1 = x^0
                sage: f2 = 10*x - x^2
                sage: f3 = 3*x^4 - 156*x^3 + 3036*x^2 - 26208*x
                sage: f = piecewise([[(0,3),f1],[(3,10),f2],[(10,20),f3]])
                sage: expected = [5, 12, 13, 14]
                sage: all(abs(e-a) < 0.001 for e,a in zip(expected, f.critical_points()))
                True

            TESTS:

            Use variables other than x (:trac:`13836`)::

                sage: R.<y> = QQ[]
                sage: f1 = y^0
                sage: f2 = 10*y - y^2
                sage: f3 = 3*y^4 - 156*y^3 + 3036*y^2 - 26208*y
                sage: f = piecewise([[(0,3),f1],[(3,10),f2],[(10,20),f3]])
                sage: expected = [5, 12, 13, 14]
                sage: all(abs(e-a) < 0.001 for e,a in zip(expected, f.critical_points()))
                True
            """
            from sage.calculus.calculus import maxima
            x = QQ[self.default_variable()].gen()
            crit_pts = []
            for domain, f in parameters:
                for interval in domain:
                    a = interval.lower()
                    b = interval.upper()
                    for root in maxima.allroots(SR(f).diff(x)==0):
                        root = float(root.rhs())
                        if a < root < b:
                            crit_pts.append(root)
            return crit_pts

        def convolution(cls, self, parameters, variable, other):
            """
            Return the convolution function,
            `f*g(t)=\int_{-\infty}^\infty f(u)g(t-u)du`, for compactly
            supported `f,g`.

            EXAMPLES::

                sage: x = PolynomialRing(QQ,'x').gen()
                sage: f = piecewise([[[0,1],1]])  ## example 0
                sage: g = f.convolution(f); g
                piecewise(x|-->x on (0, 1], x|-->-x + 2 on (1, 2]; x)
                sage: h = f.convolution(g); h
                piecewise(x|-->1/2*x^2 on (0, 1], x|-->-x^2 + 3*x - 3/2 on (1, 2], x|-->1/2*x^2 - 3*x + 9/2 on (2, 3]; x)
                sage: f = piecewise([[(0,1),1],[(1,2),2],[(2,3),1]])  ## example 1
                sage: g = f.convolution(f)
                sage: h = f.convolution(g); h
                piecewise(x|-->1/2*x^2 on (0, 1], x|-->2*x^2 - 3*x + 3/2 on (1, 3], x|-->-2*x^2 + 21*x - 69/2 on (3, 4], x|-->-5*x^2 + 45*x - 165/2 on (4, 5], x|-->-2*x^2 + 15*x - 15/2 on (5, 6], x|-->2*x^2 - 33*x + 273/2 on (6, 8], x|-->1/2*x^2 - 9*x + 81/2 on (8, 9]; x)
                sage: f = piecewise([[(-1,1),1]])                             ## example 2
                sage: g = piecewise([[(0,3),x]])
                sage: f.convolution(g)
                piecewise(x|-->1/2*x^2 + x + 1/2 on (-1, 1], x|-->2*x on (1, 2], x|-->-1/2*x^2 + x + 4 on (2, 4]; x)
                sage: g = piecewise([[(0,3),1],[(3,4),2]])
                sage: f.convolution(g)
                piecewise(x|-->x + 1 on (-1, 1], x|-->2 on (1, 2], x|-->x on (2, 3], x|-->-x + 6 on (3, 4], x|-->-2*x + 10 on (4, 5]; x)

            Check that the bugs raised in :trac:`12123` are fixed::

                sage: f = piecewise([[(-2, 2), 2]])
                sage: g = piecewise([[(0, 2), 3/4]])
                sage: f.convolution(g)
                piecewise(x|-->3/2*x + 3 on (-2, 0], x|-->3 on (0, 2], x|-->-3/2*x + 6 on (2, 4]; x)
                sage: f = piecewise([[(-1, 1), 1]])
                sage: g = piecewise([[(0, 1), x], [(1, 2), -x + 2]])
                sage: f.convolution(g)
                piecewise(x|-->1/2*x^2 + x + 1/2 on (-1, 0], x|-->-1/2*x^2 + x + 1/2 on (0, 2], x|-->1/2*x^2 - 3*x + 9/2 on (2, 3]; x)
            """
            from sage.symbolic.integration.integral import definite_integral
            f = self
            g = other
            if len(f.end_points())*len(g.end_points()) == 0:
                raise ValueError('one of the piecewise functions is nowhere defined')
            M = min(min(f.end_points()),min(g.end_points()))
            N = max(max(f.end_points()),max(g.end_points()))
            tt = SR.var('tt')
            uu = SR.var('uu')
            conv = 0
            fd,f0 = parameters[0]
            gd,g0 = next(other.items())
            if len(f)==1 and len(g)==1:
                f = f.unextend_zero()
                g = g.unextend_zero()
                a1 = fd[0].lower()
                a2 = fd[0].upper()
                b1 = gd[0].lower()
                b2 = gd[0].upper()
                i1 = f0.subs({variable: uu})
                i2 = g0.subs({variable: tt-uu})
                fg1 = definite_integral(i1*i2, uu, a1, tt-b1).subs(tt = variable)
                fg2 = definite_integral(i1*i2, uu, tt-b2, tt-b1).subs(tt = variable)
                fg3 = definite_integral(i1*i2, uu, tt-b2, a2).subs(tt = variable)
                fg4 = definite_integral(i1*i2, uu, a1, a2).subs(tt = variable)
                if a1-b1<a2-b2:
                    if a2+b1!=a1+b2:
                        h = piecewise([[(a1+b1,a1+b2),fg1],[(a1+b2,a2+b1),fg2],[(a2+b1,a2+b2),fg3]])
                    else:
                        h = piecewise([[(a1+b1,a1+b2),fg1],[(a1+b2,a2+b2),fg3]])
                else:
                    if a1+b2!=a2+b1:
                        h = piecewise([[(a1+b1,a2+b1),fg1],[(a2+b1,a1+b2),fg4],[(a1+b2,a2+b2),fg3]])
                    else:
                        h = piecewise([[(a1+b1,a2+b1),fg1],[(a2+b1,a2+b2),fg3]])
                return (piecewise([[(minus_infinity,infinity),0]]).piecewise_add(h)).unextend_zero()

            if len(f)>1 or len(g)>1:
                z = piecewise([[(0,0),0]])
                for fpiece in f.pieces():
                    for gpiece in g.pieces():
                        h = gpiece.convolution(fpiece)
                        z = z.piecewise_add(h)
                return z.unextend_zero()

        def trapezoid(cls, self, parameters, variable, N):
            """
            Return the piecewise line function defined by the trapezoid rule
            for numerical integration based on a subdivision of each domain
            interval into N subintervals.

            EXAMPLES::

                sage: f = piecewise([[[0,1], x^2], [RealSet.open_closed(1,2), 5-x^2]])
                sage: f.trapezoid(2)
                piecewise(x|-->1/2*x on (0, 1/2), x|-->3/2*x - 1/2 on (1/2, 1), x|-->7/2*x - 5/2 on (1, 3/2), x|-->-7/2*x + 8 on (3/2, 2); x)
                sage: f = piecewise([[[-1,1], 1-x^2]])
                sage: f.trapezoid(4).integral(definite=True)
                5/4
                sage: f = piecewise([[[-1,1], 1/2+x-x^3]]) ## example 3
                sage: f.trapezoid(6).integral(definite=True)
                1

            TESTS:

            Use variables or rings other than x (:trac:`13836`)::

                sage: R.<y> = QQ[]
                sage: f1 = y^2
                sage: f2 = 5-y^2
                sage: f = piecewise([[[0,1],f1], [RealSet.open_closed(1,2),f2]])
                sage: f.trapezoid(2)
                piecewise(y|-->1/2*y on (0, 1/2), y|-->3/2*y - 1/2 on (1/2, 1), y|-->7/2*y - 5/2 on (1, 3/2), y|-->-7/2*y + 8 on (3/2, 2); y)
            """
            x = QQ[self.default_variable()].gen()
            def func(x0, x1):
                f0, f1 = self(x0), self(x1)
                return [[(x0,x1),f0+(f1-f0)*(x1-x0)**(-1)*(x-x0)]]
            rsum = []
            for domain, f in parameters:
                for interval in domain:
                    a = interval.lower()
                    b = interval.upper()
                    h = (b-a)/N
                    for i in range(N):
                        x0 = a+i*h
                        x1 = a+(i+1)*h
                        rsum += func(x0, x1)
            return piecewise(rsum)

        def laplace(cls, self, parameters, variable, x='x', s='t'):
            r"""
            Returns the Laplace transform of self with respect to the variable
            var.

            INPUT:

            -  ``x`` - variable of self

            -  ``s`` - variable of Laplace transform.

            We assume that a piecewise function is 0 outside of its domain and
            that the left-most endpoint of the domain is 0.

            EXAMPLES::

                sage: x, s, w = var('x, s, w')
                sage: f = piecewise([[(0,1),1],[[1,2], 1-x]])
                sage: f.laplace(x, s)
                -e^(-s)/s + (s + 1)*e^(-2*s)/s^2 + 1/s - e^(-s)/s^2
                sage: f.laplace(x, w)
                -e^(-w)/w + (w + 1)*e^(-2*w)/w^2 + 1/w - e^(-w)/w^2

            ::

                sage: y, t = var('y, t')
                sage: f = piecewise([[[1,2], 1-y]])
                sage: f.laplace(y, t)
                (t + 1)*e^(-2*t)/t^2 - e^(-t)/t^2

            ::

                sage: s = var('s')
                sage: t = var('t')
                sage: f1(t) = -t
                sage: f2(t) = 2
                sage: f = piecewise([[[0,1],f1],[(1,infinity),f2]])
                sage: f.laplace(t,s)
                (s + 1)*e^(-s)/s^2 + 2*e^(-s)/s - 1/s^2
            """
            from sage.all import assume, exp, forget
            x = SR.var(x)
            s = SR.var(s)
            assume(s>0)
            result = 0
            for domain, f in parameters:
                for interval in domain:
                    a = interval.lower()
                    b = interval.upper()
                    result += (SR(f)*exp(-s*x)).integral(x,a,b)
            forget(s>0)
            return result

        def fourier_series_cosine_coefficient(cls, self, parameters, variable, n, L):
            r"""
            Returns the n-th Fourier series coefficient of
            `\cos(n\pi x/L)`, `a_n`.

            INPUT:


            -  ``self`` - the function f(x), defined over -L x L

            -  ``n`` - an integer n=0

            -  ``L`` - (the period)/2


            OUTPUT:
            `a_n = \frac{1}{L}\int_{-L}^L f(x)\cos(n\pi x/L)dx`

            EXAMPLES::

                sage: f(x) = x^2
                sage: f = piecewise([[(-1,1),f]])
                sage: f.fourier_series_cosine_coefficient(2,1)
                pi^(-2)
                sage: f(x) = x^2
                sage: f = piecewise([[(-pi,pi),f]])
                sage: f.fourier_series_cosine_coefficient(2,pi)
                1
                sage: f1(x) = -1
                sage: f2(x) = 2
                sage: f = piecewise([[(-pi,pi/2),f1],[(pi/2,pi),f2]])
                sage: f.fourier_series_cosine_coefficient(5,pi)
                -3/5/pi
            """
            from sage.all import cos, pi
            x = SR.var('x')
            result = 0
            for domain, f in parameters:
                for interval in domain:
                    a = interval.lower()
                    b = interval.upper()
                    result += (f*cos(pi*x*n/L)/L).integrate(x, a, b)
            return SR(result).simplify_trig()

        def fourier_series_sine_coefficient(cls, self, parameters, variable, n, L):
            r"""
            Returns the n-th Fourier series coefficient of
            `\sin(n\pi x/L)`, `b_n`.

            INPUT:


            -  ``self`` - the function f(x), defined over -L x L

            -  ``n`` - an integer n0

            -  ``L`` - (the period)/2


            OUTPUT:
            `b_n = \frac{1}{L}\int_{-L}^L f(x)\sin(n\pi x/L)dx`

            EXAMPLES::

                sage: f(x) = x^2
                sage: f = piecewise([[(-1,1),f]])
                sage: f.fourier_series_sine_coefficient(2,1)  # L=1, n=2
                0
            """
            from sage.all import sin, pi
            x = SR.var('x')
            result = 0
            for domain, f in parameters:
                for interval in domain:
                    a = interval.lower()
                    b = interval.upper()
                    result += (f*sin(pi*x*n/L)/L).integrate(x, a, b)
            return SR(result).simplify_trig()

        def fourier_series_partial_sum(cls, self, parameters, variable, N, L):
            r"""
            Returns the partial sum

            .. math::

               f(x) \sim \frac{a_0}{2} + \sum_{n=1}^N [a_n\cos(\frac{n\pi x}{L}) + b_n\sin(\frac{n\pi x}{L})],

            as a string.

            EXAMPLE::

                sage: f(x) = x^2
                sage: f = piecewise([[(-1,1),f]])
                sage: f.fourier_series_partial_sum(3,1)
                cos(2*pi*x)/pi^2 - 4*cos(pi*x)/pi^2 + 1/3
                sage: f1(x) = -1
                sage: f2(x) = 2
                sage: f = piecewise([[(-pi,pi/2),f1],[(pi/2,pi),f2]])
                sage: f.fourier_series_partial_sum(3,pi)
                -3*cos(x)/pi - 3*sin(2*x)/pi + 3*sin(x)/pi - 1/4
            """
            from sage.all import pi, sin, cos, srange
            x = self.default_variable()
            a0 = self.fourier_series_cosine_coefficient(0,L)
            result = a0/2 + sum([(self.fourier_series_cosine_coefficient(n,L)*cos(n*pi*x/L) +
                                  self.fourier_series_sine_coefficient(n,L)*sin(n*pi*x/L))
                                 for n in srange(1,N)])
            return SR(result).expand()

piecewise = PiecewiseFunction()
