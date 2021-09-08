# Do not add anything to this file.
# It will be removed soon in order to turn 'sage' into a native namespace package.
# See https://trac.sagemath.org/ticket/29705


# IPython calls this when starting up
def load_ipython_extension(*args):
    import sage.repl.ipython_extension
    sage.repl.ipython_extension.load_ipython_extension(*args)


# Monkey-patch inspect.isfunction() to support Cython functions.
def isfunction(obj):
    """
    Check whether something is a function.

    We assume that anything which has a genuine ``__code__``
    attribute (not using ``__getattr__`` overrides) is a function.
    This is meant to support Cython functions.

    EXAMPLES::

        sage: from inspect import isfunction
        sage: def f(): pass
        sage: isfunction(f)
        True
        sage: isfunction(lambda x:x)
        True
        sage: from sage.categories.coercion_methods import _mul_parent
        sage: isfunction(_mul_parent)
        True
        sage: isfunction(Integer.digits)     # unbound method
        False
        sage: isfunction(Integer(1).digits)  # bound method
        False

    Verify that ipywidgets can correctly determine signatures of Cython
    functions::

        sage: from ipywidgets.widgets.interaction import signature
        sage: from sage.dynamics.complex_dynamics.mandel_julia_helper import fast_mandelbrot_plot
        sage: signature(fast_mandelbrot_plot)  # random
        <IPython.utils._signatures.Signature object at 0x7f3ec8274e10>
    """
    # We use type(obj) instead of just obj to avoid __getattr__().
    # Some types, like methods, will return the __code__ of the
    # underlying function in __getattr__() but we don't want to
    # detect those as functions.
    return hasattr(type(obj), "__code__")

import inspect
inspect.isfunction = isfunction
