class Python:
    """
    Allows for evaluating a chunk of code without any preparsing.
    """
    def eval(self, x, globals={}, locals={}):
        x = x.strip()
        y = x.split('\n')
        if len(y) == 0:
            return ''
        s = '\n'.join(y[:-1]) + '\n'
        t = y[-1]
        try:
            z = compile(t + '\n', '', 'single')
        except SyntaxError:
            s += '\n' + t
            z = None
        #else:
        #    s += '\n' + t
        eval(compile(s, '', 'exec'), globals, locals)
        if not z is None:
            eval(z, globals, locals)
        return ''

python = Python()




