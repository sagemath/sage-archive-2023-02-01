import os, shutil
from sage.all import tmp_dir, SAGE_ROOT

n = 0

def qasm_png(s, outfile=None, res=150, verbose=False):
    """
    INPUT:
        s -- input string that describes QASM diagram
        res -- output resolution
        outfile -- None or a filename
    """
    global n
    if outfile is None:
        outfile = "%s.png"%n
        n += 1

    E = "%s/local/lib/qasm"%SAGE_ROOT
    if not os.path.exists(E):
        raise RuntimeError, "the optional qasm package is not installed"

    D = tmp_dir()
    open('%s/input.qasm'%D,'w').write(s)
    shutil.copy('%s/xyqcirc.tex'%E, '%s/xyqcirc.tex'%D)
    cmd = 'cd "%s" && python "%s/qasm2tex.py" input.qasm > input.tex && latex \\\\nonstopmode \\\\input input.tex && dvipng -D%s input.dvi -o input.png'%(D, E, res)
    if verbose:
        print cmd
    r = os.popen(cmd).read()
    if verbose:
        print r
    if not os.path.exists(D + '/input.png'):
        raise RuntimeError, r + '\n Error running qasm'
    shutil.copy(D+'/input.png', outfile)
    shutil.rmtree(D, ignore_errors=True)

class Qasm:
    def __call__(self, s, outfile=None, res=150, verbose=False):
        qasm_png(s, outfile=outfile, res=res, verbose=verbose)

    def eval(self, s, **kwds):
        qasm_png(s)
        return ''

qasm = Qasm()
