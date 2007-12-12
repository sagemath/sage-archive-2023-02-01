class Python:
    """
    Allows for evaluating a chunk of code without any preparsing.
    """
    def eval(self, x, globals={}, locals={}):
        """
        Evaluate x with given globals (locals is ignored).  This is
        specifically meant for evaluating code blocks with
        \code{%python} in the notebook.

        EXAMPLES:
            sage: from sage.misc.python import Python
            sage: python = Python()
            sage: python.eval('2+2')
            4
            ''
            sage: python.eval('a=5\nb=7\na+b')
            12
            ''

        The locals variable is ignored -- it is there only for
        completeness.  It is ignored since otherwise the following
        won't work:
            sage: python.eval("def foo():\n return 'foo'\nprint foo()\ndef mumble():\n print 'mumble',foo()\nmumble()")
            foo
            mumble foo
            ''
        """
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

        eval(compile(s, '', 'exec'), globals)

        if not z is None:
            eval(z, globals)
        return ''

python = Python()




