"""
Handling Superseded Functionality

The main mechanism in Sage to deal with superseded functionality is to
add a deprecation warning. This will be shown once, the first time
that the deprecated function is called.

Note that all doctests in the following use the trac ticket number
#13109, which is where this mandatory argument to :func:`deprecation`
was introduced.
"""



########################################################################
#       Copyright (C) 2012 William Stein <wstein@gmail.com>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#
#                  http://www.gnu.org/licenses/
########################################################################

from warnings import warn, resetwarnings
import inspect

from sage.misc.lazy_attribute import lazy_attribute


def _check_trac_number(trac_number):
    """
    Check that the argument is likely to be a valid trac issue number.

    INPUT:

    - ``trac_number`` -- anything.

    OUTPUT:

    This function returns nothing. A ``ValueError`` is raised if the
    argument can not be a valid trac number.

    EXAMPLES::

        sage: from sage.misc.superseded import _check_trac_number
        sage: _check_trac_number(1)
        sage: _check_trac_number(int(10))
        sage: _check_trac_number(long(1000))
        sage: _check_trac_number('10')
        Traceback (most recent call last):
        ...
        ValueError: The argument "10" is not a valid trac issue number.
    """
    try:
        trac_number.__index__()
        if trac_number >= 1:
            return
    except:
        pass
    raise ValueError('The argument "'+str(trac_number)+'" is not a valid trac issue number.')

def deprecation(trac_number, message):
    r"""
    Issue a deprecation warning.

    INPUT:

    - ``trac_number`` -- integer. The trac ticket number where the
      deprecation is introduced.

    - ``message`` -- string. an explanation why things are deprecated
      and by what it should be replaced.

    EXAMPLES::

        sage: def foo():
        ...    sage.misc.superseded.deprecation(13109, 'the function foo is replaced by bar')
        sage: foo()
        doctest:1: DeprecationWarning: the function foo is replaced by bar
        See http://trac.sagemath.org/13109 for details.
    """
    _check_trac_number(trac_number)
    message += '\n'
    message += 'See http://trac.sagemath.org/'+ str(trac_number) + ' for details.'
    resetwarnings()
    # Stack level 3 to get the line number of the code which called
    # the deprecated function which called this function.
    warn(message, DeprecationWarning, stacklevel=3)



class DeprecatedFunctionAlias(object):
    """
    A wrapper around methods or functions which automatically print
    the correct deprecation message. See
    :func:`deprecated_function_alias`.

    AUTHORS:

    - Florent Hivert (2009-11-23), with the help of Mike Hansen.
    - Luca De Feo (2011-07-11), printing the full module path when different from old path
    """
    def __init__(self, trac_number, func, module):
        r"""
        TESTS::

            sage: from sage.misc.superseded import deprecated_function_alias
            sage: g = deprecated_function_alias(13109, number_of_partitions)
            sage: g.__doc__
            'Deprecated: Use :func:`number_of_partitions` instead.\nSee :trac:`13109` for details.\n\n'
        """
        _check_trac_number(trac_number)
        try:
            self.__dict__.update(func.__dict__)
        except AttributeError:
            pass # Cython classes don't have __dict__
        self.func = func
        self.trac_number  = trac_number
        self.instance = None # for use with methods
        self.__module__ = module
        if type(func) == type(deprecation):
            sphinxrole = "func"
        else:
            sphinxrole = "meth"
        doc = 'Deprecated: '
        doc += 'Use :' + sphinxrole + ':`' + self.func.__name__ + '` instead.\n'
        doc += 'See :trac:`' + str(self.trac_number) + '` for details.\n\n'
        self.__doc__ = doc

    @lazy_attribute
    def __name__(self):
        """
        TESTS::

            sage: from sage.misc.superseded import deprecated_function_alias
            sage: g = deprecated_function_alias(13109, number_of_partitions)
            sage: g.__name__
            'g'

            sage: from sage.misc.superseded import deprecated_function_alias
            sage: class cls(object):
            ...      def new_meth(self): return 42
            ...      old_meth = deprecated_function_alias(13109, new_meth)
            ...
            sage: cls().old_meth.__name__
            'old_meth'

            sage: cython('\n'.join([
            ...       r"from sage.misc.superseded import deprecated_function_alias",
            ...       r"cdef class cython_cls(object):",
            ...       r"    def new_cython_meth(self):",
            ...       r"        return 1",
            ...       r"    old_cython_meth = deprecated_function_alias(13109, new_cython_meth)"
            ...   ]))
            ...
            sage: cython_cls().old_cython_meth.__name__
            'old_cython_meth'
        """
        # first look through variables in stack frames
        for frame in inspect.stack():
            for name, obj in frame[0].f_globals.iteritems():
                if obj is self:
                    return name
        # then search object that contains self as method
        import gc, copy
        gc.collect()
        def is_class(gc_ref):
            if not isinstance(gc_ref, dict):
                return False
            is_python_class = '__module__' in gc_ref
            is_cython_class = '__new__' in gc_ref
            return is_python_class or is_cython_class
        for ref in gc.get_referrers(self):
            if is_class(ref):
                ref_copy = copy.copy(ref)
                for key, val in ref_copy.iteritems():
                    if val is self:
                        return key
        raise AttributeError, "The name of this deprecated function can not be determined"

    def __call__(self, *args, **kwds):
        """
        TESTS::

            sage: from sage.misc.superseded import deprecated_function_alias
            sage: def bla(): return 42
            sage: blo = deprecated_function_alias(13109, bla)
            sage: blo()
            doctest:1: DeprecationWarning: blo is deprecated. Please use bla instead.
            See http://trac.sagemath.org/13109 for details.
            42
        """
        if self.instance is None and self.__module__ != self.func.__module__:
            other = self.func.__module__ + "." + self.func.__name__
        else:
            other = self.func.__name__

        deprecation(self.trac_number,
                    "%s is deprecated. Please use %s instead."%(self.__name__, other))
        if self.instance is None:
            return self.func(*args, **kwds)
        else:
            return self.func(self.instance, *args, **kwds)

    def __get__(self, inst, cls = None):
        """
        TESTS::

            sage: from sage.misc.superseded import deprecated_function_alias
            sage: class cls(object):
            ...      def new_meth(self): return 42
            ...      old_meth = deprecated_function_alias(13109, new_meth)
            sage: obj = cls()
            sage: obj.old_meth.instance is obj
            True
        """
        self.instance = inst
        return self


def deprecated_function_alias(trac_number, func):
    """
    Create an aliased version of a function or a method which raise a
    deprecation warning message.

    If f is a function or a method, write
    ``g = deprecated_function_alias(trac_number, f)``
    to make a deprecated aliased version of f.

    INPUT:

    - ``trac_number`` -- integer. The trac ticket number where the
      deprecation is introduced.

    - ``func`` -- the function or method to be aliased

    EXAMPLES::

        sage: from sage.misc.superseded import deprecated_function_alias
        sage: g = deprecated_function_alias(13109, number_of_partitions)
        sage: g(5)
        doctest:...: DeprecationWarning: g is deprecated. Please use sage.combinat.partition.number_of_partitions instead.
        See http://trac.sagemath.org/13109 for details.
        7

    This also works for methods::

        sage: class cls(object):
        ...      def new_meth(self): return 42
        ...      old_meth = deprecated_function_alias(13109, new_meth)
        sage: cls().old_meth()
        doctest:...: DeprecationWarning: old_meth is deprecated. Please use new_meth instead.
        See http://trac.sagemath.org/13109 for details.
        42

    Trac #11585::

        sage: def a(): pass
        sage: b = deprecated_function_alias(13109, a)
        sage: b()
        doctest:...: DeprecationWarning: b is deprecated. Please use a instead.
        See http://trac.sagemath.org/13109 for details.

    AUTHORS:

     - Florent Hivert (2009-11-23), with the help of Mike Hansen.
     - Luca De Feo (2011-07-11), printing the full module path when different from old path
    """
    _check_trac_number(trac_number)
    module_name = inspect.getmodulename(
        inspect.currentframe(1).f_code.co_filename)
    if module_name is None:
        module_name = '__main__'
    return DeprecatedFunctionAlias(trac_number, func, module_name)


def deprecated_callable_import(trac_number, module_name, globs, locs, fromlist, message=None):
    """
    Imports a list of callables into the namespace from
    which it is called.  These callables however give a deprecation
    warning whenever they are called.  This is primarily used from
    deprecating things from Sage's ``all.py`` files.

    INPUT:

    - ``trac_number`` -- integer. The trac ticket number where the
      deprecation is introduced.

    - ``param module_name`` -- string or ``None``. The name of the
      module from which to import the callables or ``None``.

    - ``globs`` -- dictionary. The ``globals()`` from where this is being called.

    - ``locs`` -- dictionary. The ``locals()`` from where this is being called.

    - ``param fromlist: -- list of strings. The list the names of the
      callables to deprecate

    - ``message`` --` string. Message to display when the deprecated functions are called.

    .. note::

       If ``module_name`` is ``None``, then no importing will be done, and
       it will be assumed that the functions have already been
       imported and are present in ``globs``

    .. warning::

       This should really only be used for functions.

    EXAMPLES::

       sage: from sage.misc.superseded import deprecated_callable_import
       sage: is_prime(3)
       True
       sage: message = "Using %(name)s from here is deprecated."
       sage: deprecated_callable_import(13109, None, globals(), locals(), ['is_prime'], message)
       sage: is_prime(3)
       doctest:...: DeprecationWarning:
       Using is_prime from here is deprecated.
       See http://trac.sagemath.org/13109 for details.
       True
       sage: del is_prime
       sage: deprecated_callable_import(13109, 'sage.rings.arith', globals(), locals(), ['is_prime'])
       sage: is_prime(3)
       doctest:...: DeprecationWarning:
       Using is_prime from here is deprecated.  If you need to use it, please import it directly from sage.rings.arith.
       See http://trac.sagemath.org/13109 for details.
       True
    """
    _check_trac_number(trac_number)
    if message is None:
        message = '\nUsing %(name)s from here is deprecated. ' + \
            'If you need to use it, please import it directly from %(module_name)s.'
    from functools import partial
    from sage.misc.misc import sage_wraps
    if module_name is None:
        mod_dict = globs
    else:
        mod_dict = __import__(module_name, globs, locs, fromlist).__dict__
    for name in fromlist:
        func = mod_dict[name]
        def wrapper(func, name, *args, **kwds):
            from sage.misc.superseded import deprecation
            deprecation(trac_number, message%{'name': name, 'module_name': module_name})
            return func(*args, **kwds)
        globs[name] = sage_wraps(func)(partial(wrapper, func, name))
    del name
