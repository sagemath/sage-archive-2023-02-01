import calculus

def var(s):
    """
    Create a symbolic variable with the name \emph{s}.

    INPUT:
        s -- a string, either a single variable name,
             or a space or comma separated list of
             variable names.

    EXAMPLES:
    We define three variables:
        sage: var('xx yy zz')
        (xx, yy, zz)

    Then we make an algebraic expression out of them.
        sage: f = xx^n + yy^n + zz^n; f
        zz^n + yy^n + xx^n

    We can substitute a new variable name for n.
        sage: f(n = var('sigma'))
        zz^sigma + yy^sigma + xx^sigma

    If you make an important builtin variable into a symbolic variable,
    you can get back the original value using restore:
        sage: var('QQ RR')
        (QQ, RR)
        sage: QQ
        QQ
        sage: restore('QQ')
        sage: QQ
        Rational Field

    We make two new variables separated by commas:
        sage: var('theta, gamma')
        (theta, gamma)
        sage: theta^2 + gamma^3
        gamma^3 + theta^2

    The new variables are of type SymbolicVariable, and belong
    to the symbolic expression ring:
        sage: type(theta)
        <class 'sage.calculus.calculus.SymbolicVariable'>
        sage: parent(theta)
        Ring of Symbolic Expressions

    NOTE: The new variable is both returned and automatically injected
    into the global namespace.  If you use var in library code, it is
    better to use sage.calculus.calculus.var, since it won't touch the
    global namespace.
    """
    G = globals()  # this is the reason the code must be in SageX.
    v = calculus.var(s)
    if isinstance(v, tuple):
        for x in v:
            G[str(x)] = x
    else:
        G[str(v)] = v
    return v

