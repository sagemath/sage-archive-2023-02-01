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
from copy import copy

from sage.misc.misc import verbose, deprecation

from sage.ext.fast_eval import fast_float, fast_float_constant, is_fast_float


def ensure_subs(f):
    if not hasattr(f, 'subs'):
        from sage.calculus.all import SR
        return SR(f)
    return f

class suboptions(object):
    def __init__(self, name, **options):
        """
        A decorator for functions which collects all keywords
        starting with name_ and collects them into a dictionary
        which will be passed on to the wrapped function as a
        dictionary called name_options.

        The keyword arguments passed into the constructor are taken
        to be default for the name_options dictionary.

        EXAMPLES::

            sage: from sage.plot.misc import suboptions
            sage: s = suboptions('arrow', size=2)
            sage: s.name
            'arrow_'
            sage: s.options
            {'size': 2}
        """
        self.name = name + "_"
        self.options = options

    def __call__(self, func):
        """
        Returns a wrapper around func

        EXAMPLES::

            sage: from sage.plot.misc import suboptions
            sage: def f(*args, **kwds): print list(sorted(kwds.items()))
            sage: f = suboptions('arrow', size=2)(f)
            sage: f(size=2)
            [('arrow_options', {'size': 2}), ('size', 2)]
            sage: f(arrow_size=3)
            [('arrow_options', {'size': 3})]
            sage: f(arrow_options={'size':4})
            [('arrow_options', {'size': 4})]
            sage: f(arrow_options={'size':4}, arrow_size=5)
            [('arrow_options', {'size': 5})]

        """
        @wraps(func)
        def wrapper(*args, **kwds):
            suboptions = copy(self.options)
            suboptions.update(kwds.pop(self.name+"options", {}))

            #Collect all the relevant keywords in kwds
            #and put them in suboptions
            for key, value in kwds.items():
                if key.startswith(self.name):
                    suboptions[key[len(self.name):]] = value
                    del kwds[key]

            kwds[self.name + "options"] = suboptions

            return func(*args, **kwds)

        from sage.misc.sageinspect import sage_getsource
        wrapper._sage_src_ = lambda: sage_getsource(func)

        return wrapper

class options(object):
    def __init__(self, **options):
        """
        A decorator for functions which allows for default options to be
        set and reset by the end user.  Additionally, if one needs to, one
        can get at the original keyword arguments passed into the
        decorator.

        TESTS::

            sage: from sage.plot.misc import options
            sage: o = options(rgbcolor=(0,0,1))
            sage: o.options
            {'rgbcolor': (0, 0, 1)}
            sage: o = options(rgbcolor=(0,0,1), __original_opts=True)
            sage: o.original_opts
            True
            sage: loads(dumps(o)).options
            {'rgbcolor': (0, 0, 1)}
        """
        self.options = options
        self.original_opts = options.pop('__original_opts', False)

    def __call__(self, func):
        """
        EXAMPLES::

            sage: from sage.plot.misc import options
            sage: o = options(rgbcolor=(0,0,1))
            sage: def f(*args, **kwds): print args, list(sorted(kwds.items()))
            sage: f1 = o(f)
            sage: f1()
            () [('rgbcolor', (0, 0, 1))]
            sage: f1(rgbcolor=1)
            () [('rgbcolor', 1)]
            sage: o = options(rgbcolor=(0,0,1), __original_opts=True)
            sage: f2 = o(f)
            sage: f2(alpha=1)
            () [('__original_opts', {'alpha': 1}), ('alpha', 1), ('rgbcolor', (0, 0, 1))]

        """
        @wraps(func)
        def wrapper(*args, **kwds):
            options = copy(wrapper.options)
            if self.original_opts:
                options['__original_opts'] = kwds
            options.update(kwds)
            return func(*args, **options)


        def defaults():
            """
            Return the default options.

            EXAMPLES::

                sage: from sage.plot.misc import options
                sage: o = options(rgbcolor=(0,0,1))
                sage: def f(*args, **kwds): print args, list(sorted(kwds.items()))
                sage: f = o(f)
                sage: f.options['rgbcolor']=(1,1,1)
                sage: f.defaults()
                {'rgbcolor': (0, 0, 1)}
            """
            return copy(self.options)

        def reset():
            """
            Reset the options to the defaults.

            EXAMPLES::

                sage: from sage.plot.misc import options
                sage: o = options(rgbcolor=(0,0,1))
                sage: def f(*args, **kwds): print args, list(sorted(kwds.items()))
                sage: f = o(f)
                sage: f.options
                {'rgbcolor': (0, 0, 1)}
                sage: f.options['rgbcolor']=(1,1,1)
                sage: f.options
                {'rgbcolor': (1, 1, 1)}
                sage: f()
                () [('rgbcolor', (1, 1, 1))]
                sage: f.reset()
                sage: f.options
                {'rgbcolor': (0, 0, 1)}
                sage: f()
                () [('rgbcolor', (0, 0, 1))]
            """
            wrapper.options = copy(self.options)

        wrapper.options = copy(self.options)
        wrapper.reset = reset
        wrapper.reset.__doc__ = """
        Reset the options to the defaults.

        Defaults:
        %s
        """%self.options

        wrapper.defaults = defaults
        wrapper.defaults.__doc__ = """
        Return the default options.

        Defaults:
        %s
        """%self.options


        from sage.misc.sageinspect import sage_getsource
        wrapper._sage_src_ = lambda: sage_getsource(func)

        return wrapper


class rename_keyword(object):
    def __init__(self, deprecated=None, **renames):
        """
        A decorator which renames keyword arguments and optionally
        deprecates the new keyword.

        INPUT:

        - ``deprecated`` - If the option being renamed is deprecated, this
          is the Sage version where the deprecation initially occurs.

        - the rest of the arguments is a list of keyword arguments in the
          form ``renamed_option='existing_option'``.  This will have the
          effect of renaming ``renamed_option`` so that the function only
          sees ``existing_option``.  If both ``renamed_option`` and
          ``existing_option`` are passed to the function, ``existing_option``
          will override the ``renamed_option`` value.

        EXAMPLES::

            sage: from sage.plot.misc import rename_keyword
            sage: r = rename_keyword(color='rgbcolor')
            sage: r.renames
            {'color': 'rgbcolor'}
            sage: loads(dumps(r)).renames
            {'color': 'rgbcolor'}

        To deprecate an old keyword::

            sage: r = rename_keyword(deprecated='Sage version 4.2', color='rgbcolor')
        """
        self.renames = renames
        self.deprecated=deprecated

    def __call__(self, func):
        """
        Rename keywords.

        EXAMPLES::

            sage: from sage.plot.misc import rename_keyword
            sage: r = rename_keyword(color='rgbcolor')
            sage: def f(*args, **kwds): print args, kwds
            sage: f = r(f)
            sage: f()
            () {}
            sage: f(alpha=1)
            () {'alpha': 1}
            sage: f(rgbcolor=1)
            () {'rgbcolor': 1}
            sage: f(color=1)
            () {'rgbcolor': 1}

        We can also deprecate the renamed keyword::

            sage: r = rename_keyword(deprecated='Sage version 4.2', deprecated_option='new_option')
            sage: def f(*args, **kwds): print args, kwds
            sage: f = r(f)
            sage: f()
            () {}
            sage: f(alpha=1)
            () {'alpha': 1}
            sage: f(new_option=1)
            () {'new_option': 1}
            sage: f(deprecated_option=1)
            doctest:...: DeprecationWarning: (Since Sage version 4.2) use the option 'new_option' instead of 'deprecated_option'
            () {'new_option': 1}
        """
        @wraps(func)
        def wrapper(*args, **kwds):
            for old_name, new_name in self.renames.items():
                if old_name in kwds and new_name not in kwds:
                    if self.deprecated is not None:
                        deprecation("use the option %r instead of %r"%(new_name,old_name),
                            version=self.deprecated)
                    kwds[new_name] = kwds[old_name]
                    del kwds[old_name]
            return func(*args, **kwds)

        from sage.misc.sageinspect import sage_getsource
        wrapper._sage_src_ = lambda: sage_getsource(func)

        return wrapper


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
        sage: sage.plot.misc.setup_for_eval_on_grid(x+y, [(x,1,-1),(y,-1,1)], return_vars=True)
        (<sage.ext...>, [(1.0, -1.0, 2.0), (-1.0, 1.0, 2.0)], [x, y])
        sage: sage.plot.misc.setup_for_eval_on_grid(x+y, [(y,1,-1),(x,-1,1)], return_vars=True)
        (<sage.ext...>, [(1.0, -1.0, 2.0), (-1.0, 1.0, 2.0)], [y, x])
    """
    if max(len(r) for r in ranges)!=min(len(r) for r in ranges):
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
