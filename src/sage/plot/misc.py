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

from sage.misc.misc import verbose

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
    def __init__(self, **renames):
        """
        A decorator which renames keyword arguments.

        EXAMPLES:
            sage: from sage.plot.misc import rename_keyword
            sage: r = rename_keyword(color='rgbcolor')
            sage: r.renames
            {'color': 'rgbcolor'}
            sage: loads(dumps(r)).renames
            {'color': 'rgbcolor'}

        """
        self.renames = renames

    def __call__(self, func):
        """
        EXAMPLES:
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
        """
        @wraps(func)
        def wrapper(*args, **kwds):
            for old_name, new_name in self.renames.items():
                if kwds.has_key(old_name) and not kwds.has_key(new_name):
                    kwds[new_name] = kwds[old_name]
                    del kwds[old_name]
            return func(*args, **kwds)

        from sage.misc.sageinspect import sage_getsource
        wrapper._sage_src_ = lambda: sage_getsource(func)

        return wrapper
