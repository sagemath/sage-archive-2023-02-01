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
from sage.sets.real_set import RealSet
from sage.symbolic.ring import SR




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
        var = kwds.pop('var', None)
        parameters = []
        domain_list = []
        for piece in function_pieces:
            domain, function = piece
            if not isinstance(domain, RealSet):
                domain = RealSet(domain)
            if domain.is_empty():
                continue
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

            sage: p = piecewise([((-2, 0), -x), ([0, 2], x)], var=x)
            sage: p.subs(x=-1)
            1
            sage: (10+p).subs(x=-1)
            11
    
        Auxiliary variables can be substituted arbitrarily::

            sage: var('x,y')
            (x, y)
            sage: p = piecewise([((-2, 0), -x^y), ([0, 2], x-y)], var=x);  p
            piecewise(x|-->-x^y on (-2, 0), x|-->x - y on [0, 2]; x)
            sage: p.subs(y=sin(y))
            piecewise(x|-->-x^sin(y) on (-2, 0), x|-->x - sin(y) on [0, 2]; x)
        """
        point = subs_map.apply_to(x, 0)
        # print 'point =', point
        if point == x:
            # substitution only in auxiliary variables
            new_params = []
            for domain, func in parameters:
                new_params.append((domain, subs_map.apply_to(func, 0)))
            return piecewise(new_params, var=x)
        if not (point.is_numeric() and point.is_real()):
            raise ValueError('substituting the piecewise variable must result in real number')

        for domain, func in parameters:
            if domain.contains(point):
                return subs_map.apply_to(func, 0)
        raise ValueError('point is not in the domain')

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
    def simplify():
        """
        Combine piecewise operands into single piecewise function

        OUTPUT:

        A piecewise function whose operands are not piecewiese if 
        possible, that is, as long as the piecewise variable is the same.

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

        def iteritems(cls, self, parameters, variable):
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
                ValueError: point is not in the domain
            """
            self = piecewise(parameters, var=variable)
            substitution = dict()
            for k, v in kwds.iteritems():
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
                {CommutativeRings.element_class}(v_0)
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
                ValueError: point is not in the domain

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

        def pieces(cls, self, parameters, variable):
            """
            Return the "pieces"

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




piecewise = PiecewiseFunction()

def Piecewise(*args, **kwds):
    """
    Deprecated spelling of ``piecewise``

    EXAMPLES::

        sage: Piecewise([((-1, 0), -x), ([0, 1], x)], var=x)
        doctest:...: DeprecationWarning: use lower-case piecewise instead
        See http://trac.sagemath.org/14801 for details.
        piecewise(x|-->-x on (-1, 0), x|-->x on [0, 1]; x)
    """
    from sage.misc.superseded import deprecation
    deprecation(14801, 'use lower-case piecewise instead')
    return piecewise(*args, **kwds)



