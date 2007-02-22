from calculus import SER, SymbolicExpression

def diff(f, *args):
    if not isinstance(f, SymbolicExpression):
        f = SER(f)
    return f.derivative(*args)

derivative = diff


