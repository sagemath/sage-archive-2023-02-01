r"""
Fast Numerical Evaluation

For many applications such as numerical integration, differential
equation approximation, plotting a 3d surface, optimization problems,
monte-carlo simulations, etc., one wishes to pass around and evaluate
a single algebraic expression many, many times at various floating
point values. Doing this via recursive calls over a python
representation of the object (even if Maxima or other outside packages
are not involved) is extremely inefficient.

Up until now the solution has been to use lambda expressions, but this
is neither intuitive, Sage-like, nor efficient (compared to operating
on raw C doubles).  This module provides a representation of algebraic
expression in Reverse Polish Notation, and provides an efficient
interpreter on C double values as a callable python object. It does
what it can in C, and will call out to Python if necessary.

Essential to the understanding of this class is the distinction
between symbolic expressions and callable symbolic expressions (where
the latter binds argument names to argument positions). The
``*vars`` parameter passed around encapsulates this information.

See the function ``fast_float(f, *vars)`` to create a fast-callable
version of f.

.. NOTE::

    Sage temporarily has two implementations of this functionality ;
    one in this file, which will probably be deprecated soon, and one in
    fast_callable.pyx.  The following instructions are for the old
    implementation; you probably want to be looking at fast_callable.pyx
    instead.

To provide this interface for a class, implement ``fast_float_(self, *vars)``.  The basic building blocks are
provided by the functions ``fast_float_constant`` (returns a
constant function), ``fast_float_arg`` (selects the ``n``-th value
when called with ``\ge_n`` arguments), and ``fast_float_func`` which
wraps a callable Python function. These may be combined with the
standard Python arithmetic operators, and support many of the basic
math functions such ``sqrt``, ``exp``, and trig functions.

TESTS:

This used to segfault because of an assumption that assigning None to a
variable would raise a TypeError::

    sage: from sage.ext.fast_eval import fast_float_arg, fast_float
    sage: fast_float_arg(0)+None
    Traceback (most recent call last):
    ...
    TypeError

AUTHORS:

- Robert Bradshaw (2008-10): Initial version
"""

#*****************************************************************************
#       Copyright (C) 2008 Robert Bradshaw <robertwb@math.washington.edu>
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 2 of the License, or
# (at your option) any later version.
#                  http://www.gnu.org/licenses/
#*****************************************************************************

from sage.ext.fast_callable import fast_callable, Wrapper


def fast_float(f, *vars, old=None, expect_one_var=False):
    """
    Tries to create a function that evaluates f quickly using
    floating-point numbers, if possible.  There are two implementations
    of fast_float in Sage; by default we use the newer, which is
    slightly faster on most tests.

    On failure, returns the input unchanged.

    INPUT:

    - ``f``    -- an expression
    - ``vars`` -- the names of the arguments
    - ``old``  -- deprecated, do not use
    - ``expect_one_var`` -- don't give deprecation warning if ``vars`` is
      omitted, as long as expression has only one var

    EXAMPLES::

        sage: from sage.ext.fast_eval import fast_float
        sage: x,y = var('x,y')
        sage: f = fast_float(sqrt(x^2+y^2), 'x', 'y')
        sage: f(3,4)
        5.0

    Specifying the argument names is essential, as fast_float objects
    only distinguish between arguments by order. ::

        sage: f = fast_float(x-y, 'x','y')
        sage: f(1,2)
        -1.0
        sage: f = fast_float(x-y, 'y','x')
        sage: f(1,2)
        1.0
    """
    if old:
        raise ValueError("the old implementation of fast_float has been removed")
    if old is not None:
        from sage.misc.superseded import deprecation
        deprecation(32234, "passing old=False to fast_float is deprecated")

    if isinstance(f, (tuple, list)):
        return tuple([fast_float(x, *vars, expect_one_var=expect_one_var) for x in f])

    cdef int i
    for i from 0 <= i < len(vars):
        if not isinstance(vars[i], str):
            v = str(vars[i])
            # inexact generators display as 1.00..0*x
            if '*' in v:
                v = v[v.index('*')+1:]
            vars = vars[:i] + (v,) + vars[i+1:]

    try:
        return fast_callable(f, vars=vars, domain=float,
                             expect_one_var=expect_one_var)
    except AttributeError:
        pass

    try:
        from sage.symbolic.ring import SR
        return fast_float(SR(f), *vars)
    except TypeError:
        pass

    if f is None:
        raise TypeError("no way to make fast_float from None")

    return f


def is_fast_float(x):
    return isinstance(x, Wrapper)
