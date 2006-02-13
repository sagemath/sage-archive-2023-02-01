r"""
Interface to Magma

\note{You must have \code{magma} installed on your computer
for this interface to work.   Magma is not free, so it is
not included with \sage, but you can obtain it from
\url{http://magma.maths.usyd.edu.au/}.}   You do not
have to install any optional \sage packages.


SAGE provides an interface to the Magma computational algebra
system.  This system provides extensive functionality for
number theory, group theory, combinatorics and algebra.

The Magma interface offers three pieces of functionality:
\begin{enumerate}

\item \code{magma_console()} -- A function that dumps you
into an interactive command-line Magma session.

\item \code{magma(expr)} -- Evaluation of arbitrary Magma
expressions, with the result returned as a string.

\item \code{magma.new(expr)} -- Creation of a SAGE object that wraps a
Magma object.  This provides a Pythonic interface to Magma.  For example,
if \code{f=magma.new(10)}, then \code{f.Factors()} returns the prime
factorization of $10$ computed using Magma.

\end{enumerate}



\subsection{Parameters}
Some Magma functions have optional ``parameters'', which
are arguments that in Magma go after a colon.  In SAGE,
you pass these using named function arguments.  For example,

    sage: E = magma.new('EllipticCurve([0,1,1,-1,0])')
    sage: E.Rank(Bound = 5)
    0

\subsection{Multiple Return Values}

Some Magma functions return more than one value.
You can control how many you get using the \code{nvals}
named parameter to a function call:

    sage: n = magma.new(100)
    sage: n.IsSquare(nvals = 1)
    true
    sage: n.IsSquare(nvals = 2)
    (true, 10)

\subsection{Long Input}
The Magma interface reads in even very long input (using files) in a
robust manner.

    sage: t = '"%s"'%10^10000   # ten thousand character string.
    sage: a = magma.eval(t)
    sage: a = magma(t)
"""

#*****************************************************************************
#       Copyright (C) 2005 William Stein <wstein@ucsd.edu>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#
#    This code is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
#    General Public License for more details.
#
#  The full text of the GPL is available at:
#
#                  http://www.gnu.org/licenses/
#*****************************************************************************

from sage.structure.element import RingElement
from expect import console, Expect, ExpectElement, FunctionElement
PROMPT = ">>>"

class Magma(Expect):
    """
    Interface to the Magma interpreter.
    """
    def __init__(self, maxread=10000, script_subdirectory="user", logfile=None, server=None):
        Expect.__init__(self,
                        name = "magma",
                        prompt = ">>SAGE>>",
                        command = "magma -b -n",
                        maxread = maxread,
                        server = server,
                        script_subdirectory = script_subdirectory,
                        restart_on_ctrlc = True,
                        logfile = logfile,
                        eval_using_file_cutoff=100)
        # We use "-n" above in the Magma startup command so
        # local user startup changes are not read.

        self.__seq = 0

    def __reduce__(self):
        return reduce_load_Magma, tuple([])

    def _read_in_file_command(self, filename):
        return 'load "%s";'%filename

    def _continuation_prompt(self):
        return self._prompt

    def eval(self, x):
        x = str(x).rstrip()
        if len(x) == 0 or x[len(x) - 1] != ';':
            x += ';'
        ans = Expect.eval(self, x).replace('\\\n','')
        if ans.find("Runtime error") != -1:
            raise RuntimeError, "Error evaluation Magma code.\nIN:%s\nOUT:%s"%(x, ans)
        return ans

    def _start(self):
        self._change_prompt('>')
        Expect._start(self)
        self.eval('SetPrompt("%s"); SetLineEditor(false); SetColumns(0);'%PROMPT)
        self._change_prompt(PROMPT)
        self.expect().expect(PROMPT)
        self.expect().expect(PROMPT)

    def set(self, var, value):
        """
        Set the variable var to the given value.
        """
        out = self.eval("%s := %s"%(var, value))
        if out.lower().find("error") != -1:
            raise TypeError, "Error executing Magma code:\n%s"%out

    def get(self, var):
        """
        Get the value of the variable var.
        """
        return self.eval("%s"%var)

    #def clear(self, var):
    #    """
    #    Clear the variable named var.
    #    """
        #self.eval("delete %s"%var)
    #    self.eval("%s:=0"%var)

    def attach(self, filename):
        self.eval('Attach("%s")'%filename)

    def _next_var_name(self):
        if self.__seq == 0:
            self.eval('sage := [* *];')
        else:
            self.eval('Append(~sage, 0);')
        self.__seq += 1
        return 'sage[%s]'%self.__seq

    def function_call(self, function, args=[], params={}, nvals=1):
        if not isinstance(args, list):
            args = [args]
        for i in range(len(args)):
            if not isinstance(args[i], ExpectElement):
                args[i] = self(args[i])
        nvals = int(nvals)
        if len(params) == 0:
            par = ''
        else:
            par = ' : ' + ','.join(['%s:=%s'%(a,b) for a,b in params.items()])

        fun = "%s(%s%s)"%(function, ",".join([s.name() for s in args]), par)
        if nvals <= 0:
            out = self.eval(fun)
            ans = None
        elif nvals == 1:
            return self(fun)
        else:
            v = [self._next_var_name() for _ in range(nvals)]
            vars = ", ".join(v)
            cmd = "%s := %s;"%(vars, fun)
            out = self.eval(cmd)
            ans = tuple([MagmaElement(self, x, is_name = True) for x in v])

        if out.lower().find("error") != -1:
            raise TypeError, "Error executing Magma code:\n%s"%out
        return ans

    #def new(self, x):
    #    if isinstance(x, MagmaElement) and x.parent() == self:
    #        return x
    #    return MagmaElement(self, x)

    def _object_class(self):
        return MagmaElement

    def _left_list_delim(self):
        return "[*"

    def _right_list_delim(self):
        return "*]"

    def console(self):
        magma_console()

    def version(self):
        return magma_version()

    def help(self, s):
        print magma.eval('? %s'%s)

class MagmaFunctionElement(FunctionElement):
    def __call__(self, *args, **kwds):
        nvals = 1
        if len(kwds) > 0:
            if kwds.has_key('nvals'):
                nvals = kwds['nvals']
                del kwds['nvals']
        M = self._obj.parent()
        return M.function_call(self._name,
                               [self._obj.name()] + list(args),
                               params = kwds,
                               nvals = nvals)

    def help(self):
        M = self._obj.parent()
        t = str(self.Type())
        print M(self._name)


class MagmaElement(ExpectElement):
    def __getattr__(self, attrname):
        return MagmaFunctionElement(self, attrname)

    def methods(self):
        return self.parent()('ListSignatures(%s)'%self.Type())




magma = Magma()

def reduce_load_Magma():
    return magma

def magma_console():
    console('magma')

def magma_version():
    t = tuple([int(n) for n in magma.eval('GetVersion()').split()])
    return t, 'V%s.%s-%s'%t
