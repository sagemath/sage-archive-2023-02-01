import os

from misc import tmp_filename

class Sh:
    r"""
    Evaluates a shell script and returns the output.

    To use this from the notebook type \code{\%sh}
    at the beginning of the input cell.

    The working directory is preserved between calls.
    """
    def __init__(self):
        self._curdir = None

    def chdir(self, dir):
        self._curdir = dir

    def eval(self, x, globals={}, locals={}):
        if self._curdir is None:
            self._curdir = os.path.abspath('.')
        t = tmp_filename()
        out = tmp_filename()
        w = open(t,'w')
        w.write('cd "%s"\n'%self._curdir)
        w.write(x)
        w.write('\npwd\n')
        w.close()
        os.system('chmod +x "%s"; "%s" > "%s"'%(t,t,out))
        os.unlink(t)
        s = open(out).read()
        os.unlink(out)
        t = s.split('\n')
        self._curdir = os.path.abspath(t[-1])
        return '\n'.join(t[:-1])

sh = Sh()
