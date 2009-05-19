import calculus
from sage.symbolic.ring import var as new_var

def var(s, ns=True):
    r"""
    Create a symbolic variable with the name \emph{s}.

    INPUT:

    - ``s`` - a string, either a single variable name, or a space or
      comma separated list of variable names.

    .. note::

       The new variable is both returned and automatically injected
       into the global namespace.  If you use var in library code, it
       is better to use sage.calculus.calculus.var, since it won't
       touch the global namespace.

    EXAMPLES:

    We define some symbolic variables::

        sage: var('n xx yy zz')
        (n, xx, yy, zz)

    Then we make an algebraic expression out of them::

        sage: f = xx^n + yy^n + zz^n; f
        xx^n + yy^n + zz^n

    We can substitute a new variable name for n::

        sage: f(n = var('sigma'))
        xx^sigma + yy^sigma + zz^sigma

    If you make an important builtin variable into a symbolic variable,
    you can get back the original value using restore::

        sage: var('QQ RR')
        (QQ, RR)
        sage: QQ
        QQ
        sage: restore('QQ')
        sage: QQ
        Rational Field

    We make two new variables separated by commas::

        sage: var('theta, gamma')
        (theta, gamma)
        sage: theta^2 + gamma^3
        gamma^3 + theta^2

    The new variables are of type Expresion, and belong
    to the symbolic expression ring::

        sage: type(theta)
        <type 'sage.symbolic.expression.Expression'>
        sage: parent(theta)
        Symbolic Ring
    """
    G = globals()  # this is the reason the code must be in Cython.
    if ns:
        v = new_var(s)
    else:
        v = calculus.var(s)
    if isinstance(v, tuple):
        for x in v:
            G[repr(x)] = x
    else:
        G[repr(v)] = v
    return v

def function(s, *args):
    """
    Create a formal symbolic function with the name *s*.

    INPUT:

    - ``s`` - a string, either a single variable name, or a space or
      comma separated list of variable names.

    .. note::

       The new function is both returned and automatically injected
       into the global namespace.  If you use var in library code, it
       is better to use sage.calculus.calculus.function, since it
       won't touch the global namespace.

    EXAMPLES::

    We create a formal function called supersin::

        sage: f = function('supersin', x)
        sage: f
        supersin(x)

    We can immediately use supersin in symbolic expressions::

        sage: y, z, A = var('y z A')
        sage: supersin(y+z) + A^3
        A^3 + supersin(y + z)

    We can define other functions in terms of supersin::

        sage: g(x,y) = supersin(x)^2 + sin(y/2)
        sage: g
        (x, y) |--> supersin(x)^2 + sin(1/2*y)
        sage: g.diff(y)
        (x, y) |--> 1/2*cos(1/2*y)
        sage: k = g.diff(x); k
        (x, y) |--> 2*supersin(x)*D[0](supersin)(x)

    In Sage 4.0, you must now use the :meth:`substitute_function`
    method to replace functions::

        sage: k.substitute_function(supersin, sin)
        2*sin(x)*cos(x)
    """
    if len(args) > 0:
        return function(s)(*args)

    G = globals()  # this is the reason the code must be in Cython.
    v = calculus.function(s)
    if isinstance(v, tuple):
        for x in v:
            G[repr(x)] = x
    else:
        G[repr(v)] = v
    return v


def clear_vars():
    """
    Delete all 1-letter symbolic variables that are predefined at
    startup of Sage.  Any one-letter global variables that are not
    symbolic variables are not cleared.

    EXAMPLES::

        sage: var('x y z')
        (x, y, z)
        sage: (x+y)^z
        (x + y)^z
        sage: k = 15
        sage: clear_vars()
        sage: (x+y)^z
        Traceback (most recent call last):
        ...
        NameError: name 'x' is not defined
        sage: expand((e + i)^2)
        2*I*e + e^2 - 1
        sage: k
        15
    """
    G = globals()
    from sage.symbolic.ring import is_SymbolicVariable
    for i in range(65,65+26) + range(97,97+26):
        if G.has_key(chr(i)) and is_SymbolicVariable(G[chr(i)]):
            #We check to see if there is a corresponding pyobject
            #associated with the expression.  This will work for
            #constants which we want to keep by will fail for
            #variables that we want to delete.
            try:
                G[chr(i)].pyobject()
            except TypeError:
                del G[chr(i)]



