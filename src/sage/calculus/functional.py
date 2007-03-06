from calculus import SER, SymbolicExpression, CallableFunction

def diff(f, *args):
    """
    Formally differentiate the function f.
    """
    if isinstance(f, CallableFunction):
        return f.derivative(*args)
    if not isinstance(f, SymbolicExpression):
        f = SER(f)
    return f.derivative(*args)

derivative = diff

def integrate(f, *args):
    """
    Returns the integral of f.

    EXAMPLES:
        sage: integrate(sin(x), x)
        -cos(x)
        sage: integrate(1/x, 1, 2)


    """
    if isinstance(f, CallableFunction):
        return f.derivative(*args)
    if not isinstance(f, SymbolicExpression):
        f = SER(f)
    return f.integral(*args)

def solve(f, *args):
    return f.solve(*args)
