import os

class Sh:
    r"""
    Evaluates a shell script and returns the output.

    To use this from the notebook type \code{\%sh}
    at the beginning of the input cell.

    The working directory is preserved between calls.
    """
    def __init__(self):
        self._curdir = os.path.abspath('.')

    def eval(self, x, globals={}, locals={}):
        w = open('t','w')
        w.write('cd %s\n'%self._curdir)
        w.write(x)
        w.write('\npwd\n')
        w.close()
        os.system('chmod +x t; ./t > out')
        s = open('out').read()
        os.system('rm t out')
        t = s.split('\n')
        self._curdir = os.path.abspath(t[-1])
        return '\n'.join(t[:-1])

sh = Sh()
