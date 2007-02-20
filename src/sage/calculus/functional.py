from calculus import SER, SymbolicExpression, CallableFunction

def diff(f, *args):
    if isinstance(f, CallableFunction):
        return f.derivative(*args)
    if not isinstance(f, SymbolicExpression):
        f = SER(f)
    return f.derivative(*args)

derivative = diff

def solve(f, *args):
    return f.solve(*args)
