#*****************************************************************************
#  Distributed under the terms of the GNU General Public License (GPL)
#
#    This code is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
#    General Public License for more details.
#
#  The full text of the GPL is available at:
#
#                  http://www.gnu.org/licenses/
#*****************************************************************************

from functools import wraps

from sage.ext.fast_eval import fast_float, fast_float_constant, is_fast_float

from sage.structure.element import is_Vector

def setup_for_eval_on_grid(funcs, ranges, plot_points=None, return_vars=False):
    """
    Calculate the necessary parameters to construct a list of points,
    and make the functions fast_callable.

    INPUT:

    -  ``funcs`` - a function, or a list, tuple, or vector of functions

    - ``ranges`` - a list of ranges.  A range can be a 2-tuple of
      numbers specifying the minimum and maximum, or a 3-tuple giving
      the variable explicitly.

    - ``plot_points`` - a tuple of integers specifying the number of
      plot points for each range.  If a single number is specified, it
      will be the value for all ranges.  This defaults to 2.

    - ``return_vars`` - (default False) If True, return the variables,
      in order.


    OUTPUT:


    - ``fast_funcs`` - if only one function passed, then a fast
       callable function.  If funcs is a list or tuple, then a tuple
       of fast callable functions is returned.

    - ``range_specs`` - a list of range_specs: for each range, a
       tuple is returned of the form (range_min, range_max,
       range_step) such that ``srange(range_min, range_max,
       range_step, include_endpoint=True)`` gives the correct points
       for evaluation.

    EXAMPLES::

        sage: x,y,z=var('x,y,z')
        sage: f(x,y)=x+y-z
        sage: g(x,y)=x+y
        sage: h(y)=-y
        sage: sage.plot.misc.setup_for_eval_on_grid(f, [(0, 2),(1,3),(-4,1)], plot_points=5)
        (<sage.ext...>, [(0.0, 2.0, 0.5), (1.0, 3.0, 0.5), (-4.0, 1.0, 1.25)])
        sage: sage.plot.misc.setup_for_eval_on_grid([g,h], [(0, 2),(-1,1)], plot_points=5)
        ((<sage.ext...>, <sage.ext...>), [(0.0, 2.0, 0.5), (-1.0, 1.0, 0.5)])
        sage: sage.plot.misc.setup_for_eval_on_grid([sin,cos], [(-1,1)], plot_points=9)
        ((<sage.ext...>, <sage.ext...>), [(-1.0, 1.0, 0.25)])
        sage: sage.plot.misc.setup_for_eval_on_grid([lambda x: x^2,cos], [(-1,1)], plot_points=9)
        ((<function <lambda> ...>, <sage.ext...>), [(-1.0, 1.0, 0.25)])
        sage: sage.plot.misc.setup_for_eval_on_grid([x+y], [(x,-1,1),(y,-2,2)])
        ((<sage.ext...>,), [(-1.0, 1.0, 2.0), (-2.0, 2.0, 4.0)])
        sage: sage.plot.misc.setup_for_eval_on_grid(x+y, [(x,-1,1),(y,-1,1)], plot_points=[4,9])
        (<sage.ext...>, [(-1.0, 1.0, 0.66666666666666663), (-1.0, 1.0, 0.25)])
        sage: sage.plot.misc.setup_for_eval_on_grid(x+y, [(x,-1,1),(y,-1,1)], plot_points=[4,9,10])
        Traceback (most recent call last):
        ...
        ValueError: plot_points must be either an integer or a list of integers, one for each range
        sage: sage.plot.misc.setup_for_eval_on_grid(x+y, [(1,-1),(y,-1,1)], plot_points=[4,9,10])
        Traceback (most recent call last):
        ...
        ValueError: Some variable ranges specify variables while others do not
        sage: sage.plot.misc.setup_for_eval_on_grid(x+y, [(1,-1),(-1,1)], plot_points=5)
        doctest:...: DeprecationWarning: Unnamed ranges for more than one variable is deprecated and will be removed from a future release of Sage; you can used named ranges instead, like (x,0,2)
        (<sage.ext...>, [(1.0, -1.0, 0.5), (-1.0, 1.0, 0.5)])
        sage: sage.plot.misc.setup_for_eval_on_grid(x+y, [(y,1,-1),(x,-1,1)], plot_points=5)
        (<sage.ext...>, [(1.0, -1.0, 0.5), (-1.0, 1.0, 0.5)])
        sage: sage.plot.misc.setup_for_eval_on_grid(x+y, [(x,1,-1),(x,-1,1)], plot_points=5)
        Traceback (most recent call last):
        ...
        ValueError: range variables should be distinct, but there are duplicates
        sage: sage.plot.misc.setup_for_eval_on_grid(x+y, [(x,1,1),(y,-1,1)])
        Traceback (most recent call last):
        ...
        ValueError: plot start point and end point must be different
        sage: sage.plot.misc.setup_for_eval_on_grid(x+y, [(x,1,-1),(y,-1,1)], return_vars=True)
        (<sage.ext...>, [(1.0, -1.0, 2.0), (-1.0, 1.0, 2.0)], [x, y])
        sage: sage.plot.misc.setup_for_eval_on_grid(x+y, [(y,1,-1),(x,-1,1)], return_vars=True)
        (<sage.ext...>, [(1.0, -1.0, 2.0), (-1.0, 1.0, 2.0)], [y, x])
    """
    if max(map(len, ranges)) != min(map(len, ranges)):
        raise ValueError, "Some variable ranges specify variables while others do not"

    if len(ranges[0])==3:
        vars = [r[0] for r in ranges]
        ranges = [r[1:] for r in ranges]
        if len(set(vars))<len(vars):
            raise ValueError, "range variables should be distinct, but there are duplicates"
    else:
        vars, free_vars = unify_arguments(funcs)
        if len(free_vars)>1:
            from sage.misc.misc import deprecation
            deprecation("Unnamed ranges for more than one variable is deprecated and will be removed from a future release of Sage; you can used named ranges instead, like (x,0,2)")

    # pad the variables if we don't have enough
    nargs = len(ranges)
    if len(vars)<nargs:
        vars += ('_',)*(nargs-len(vars))

    ranges = [[float(z) for z in r] for r in ranges]

    if plot_points is None:
        plot_points=2

    if not isinstance(plot_points, (list, tuple)):
        plot_points = [plot_points]*len(ranges)
    elif len(plot_points)!=nargs:
        raise ValueError, "plot_points must be either an integer or a list of integers, one for each range"

    plot_points = [int(p) if p>=2 else 2 for p in plot_points]
    range_steps = [abs(range[1] - range[0])/(p-1) for range, p in zip(ranges, plot_points)]
    if min(range_steps) == float(0):
        raise ValueError, "plot start point and end point must be different"

    options={}
    if nargs==1:
        options['expect_one_var']=True

    if is_Vector(funcs):
        funcs = list(funcs)

    #TODO: raise an error if there is a function/method in funcs that takes more values than we have ranges

    if return_vars:
        return fast_float(funcs, *vars,**options), [tuple(range+[range_step]) for range,range_step in zip(ranges, range_steps)], vars
    else:
        return fast_float(funcs, *vars,**options), [tuple(range+[range_step]) for range,range_step in zip(ranges, range_steps)]


def unify_arguments(funcs):
    """
    Returns a tuple of variables of the functions, as well as the
    number of "free" variables (i.e., variables that defined in a
    callable function).

    INPUT:

    - ``funcs`` -- a list of functions; these can be symbolic
            expressions, polynomials, etc

    OUTPUT: functions, expected arguments

    - A tuple of variables in the functions

    - A tuple of variables that were "free" in the functions

    EXAMPLES:

        sage: x,y,z=var('x,y,z')
        sage: f(x,y)=x+y-z
        sage: g(x,y)=x+y
        sage: h(y)=-y
        sage: sage.plot.misc.unify_arguments((f,g,h))
        ((x, y, z), (z,))
        sage: sage.plot.misc.unify_arguments((g,h))
        ((x, y), ())
        sage: sage.plot.misc.unify_arguments((f,z))
        ((x, y, z), (z,))
        sage: sage.plot.misc.unify_arguments((h,z))
        ((y, z), (z,))
        sage: sage.plot.misc.unify_arguments((x+y,x-y))
        ((x, y), (x, y))
    """
    from sage.symbolic.callable import is_CallableSymbolicExpression

    vars=set()
    free_variables=set()
    if not isinstance(funcs, (list, tuple)):
        funcs=[funcs]

    for f in funcs:
        if is_CallableSymbolicExpression(f):
            f_args=set(f.arguments())
            vars.update(f_args)
        else:
            f_args=set()

        try:
            free_vars = set(f.variables()).difference(f_args)
            vars.update(free_vars)
            free_variables.update(free_vars)
        except AttributeError:
            # we probably have a constant
            pass
    return tuple(sorted(vars, key=lambda x: str(x))), tuple(sorted(free_variables, key=lambda x: str(x)))

#For backward compatibility -- see #9907.
from sage.misc.decorators import options, suboptions, rename_keyword

def _multiple_of_constant(n,pos,const):
    """
    Function for internal use in formatting ticks on axes with
    nice-looking multiples of various symbolic constants, such
    as `\pi` or `e`.  Should only be used via keyword argument
    `tick_formatter` in :meth:`plot.show`.  See documentation
    for the matplotlib.ticker module for more details.

    EXAMPLES:

    Here is the intended use::

        sage: plot(sin(x), (x,0,2*pi), ticks=pi/3, tick_formatter=pi)

    Here is an unintended use, which yields unexpected (and probably
    undesired) results::

        sage: plot(x^2, (x, -2, 2), tick_formatter=pi)

    We can also use more unusual constant choices::

        sage: plot(ln(x), (x,0,10), ticks=e, tick_formatter=e)
        sage: plot(x^2, (x,0,10), ticks=[sqrt(2),8], tick_formatter=sqrt(2))
    """
    from sage.misc.latex import latex
    from sage.rings.arith import convergents
    c=[i for i in convergents(n/const.n()) if i.denominator()<12]
    return '$%s$'%latex(c[-1]*const)
