r"""
Interface to Magma

\note{You must have \code{magma} installed on your computer
for this interface to work.   Magma is not free, so it is
not included with \sage, but you can obtain it from
\url{http://magma.maths.usyd.edu.au/}.}   You do not
have to install any optional \sage packages.


    Type \code{magma.[tab]} for a list of all the functions available
    from your Magma install.  Type \code{magma.[tab]?} for Magma's
    help about a given function.  Type \code{magma(...)} to create
    a new Magma object, and \code{magma.eval(...)} to run a string
    using Magma (and get the result back as a string).

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

    sage: E = magma('EllipticCurve([0,1,1,-1,0])')
    sage: E.Rank(Bound = 5)
    0

\subsection{Multiple Return Values}

Some Magma functions return more than one value.
You can control how many you get using the \code{nvals}
named parameter to a function call:

    sage: n = magma(100)
    sage: n.IsSquare(nvals = 1)
    true
    sage: n.IsSquare(nvals = 2)
    (true, 10)
    sage: n = magma(-2006)
    sage: n.Factorization()
    [ <2, 1>, <17, 1>, <59, 1> ]
    sage: n.Factorization(nvals=2)
    ([ <2, 1>, <17, 1>, <59, 1> ], -1)

\subsection{Long Input}
The Magma interface reads in even very long input (using files) in a
robust manner.

    sage: t = '"%s"'%10^10000   # ten thousand character string.
    sage: a = magma.eval(t)
    sage: a = magma(t)

AUTHOR:
    -- William Stein (2005): initial version
    -- William Stein (2006-02-28): added extensive tab completion and interactive
                                   IPython documentation support.
    -- William Stein (2006-03-09): added nvals argument for magma.functions...
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

import sys

from sage.structure.element import RingElement
from expect import console, Expect, ExpectElement, ExpectFunction, FunctionElement
PROMPT = ">>>"

import sage.misc.misc

INTRINSIC_CACHE = '%s/magma_intrinsic_cache.sobj'%sage.misc.misc.DOT_SAGE

class Magma(Expect):
    """
    Interface to the Magma interpreter.

    Type \code{magma.[tab]} for a list of all the functions available
    from your Magma install.  Type \code{magma.[tab]?} for Magma's
    help about a given function.  Type \code{magma(...)} to create
    a new Magma object, and \code{magma.eval(...)} to run a string
    using Magma (and get the result back as a string).

    EXAMPLES:

    You must use nvals = 0 to call a function that doesn't return
    anything, otherwise you'll get an error.  (nvals is the number
    of return values.)

        sage: magma.SetDefaultRealFieldPrecision(200, nvals=0)  # optional and requires MAGMA >= v2.12

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

    def _post_process_from_file(self, s):
        if not isinstance(s, str):
            raise RuntimeError, "Error evaluating object in %s:\n%s"%(self,s)
        i = s.find('\n')
        return s[i+1:]

    def _continuation_prompt(self):
        return self._prompt

    def __getattr__(self, attrname):
        if attrname[:1] == "_":
            raise AttributeError
        return MagmaFunction(self, attrname)

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

    def objgens(self, value, gens):
        var = self._next_var_name()
        value = self(value)
        out = self.eval("_z<%s> := %s; %s := _z"%(gens, value.name(), var))
        if out.lower().find("error") != -1:
            raise TypeError, "Error executing Magma code:\n%s"%out
        return self(var)

    def __call__(self, x, gens=None):
        """
        EXAMPLES:
            sage: magma(EllipticCurve('37a'))                   # optional
            Elliptic Curve defined by y^2 + y = x^3 - x over Rational Field
            sage: magma('EllipticCurve([GF(5)|1,2,3,4,1])')     # optional
            Elliptic Curve defined by y^2 + x*y + 3*y = x^3 + 2*x^2 + 4*x + 1 over GF(5)
            sage: magma('PowerSeriesRing(Rationals())', 't')    # optional
            Power series ring in t over Rational Field
            sage: magma('PolynomialRing(RationalField(), 3)', 'x,y,z')  # optional
            Polynomial ring of rank 3 over Rational Field
            Lexicographical Order
            Variables: x, y, z
        """
        if gens is None:
            return Expect.__call__(self, x)
        return self.objgens(x, gens)


    #def clear(self, var):
    #    """
    #    Clear the variable named var.
    #    """
        #self.eval("delete %s"%var)
    #    self.eval("%s:=0"%var)

    def attach(self, filename):
        self.eval('Attach("%s")'%filename)
    Attach = attach

    def attach_spec(self, filename):
        self.eval('AttachSpec("%s")'%filename)
    AttachSpec = attach_spec

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


    # Usually "Sequences" are what you want in MAGMA, not "lists".
    # It's very painful using the interface without this.
    def _left_list_delim(self):
        #return "[*"
        return "["

    def _right_list_delim(self):
        #return "*]"
        return "]"

    def _assign_symbol(self):
        return ":="

    def _equality_symbol(self):
        return 'eq'

    # For efficiency purposes, you should definitely override these
    # in your derived class.
    def _true_symbol(self):
        return 'true'

    def _false_symbol(self):
        return 'false'


    def console(self):
        magma_console()

    def version(self):
        return magma_version()

    def help(self, s):
        print self.eval('? %s'%s)

    def trait_names(self, verbose=True, use_disk_cache=True):
        try:
            return self.__trait_names
        except AttributeError:
            import sage.misc.persist
            if use_disk_cache:
                try:
                    self.__trait_names = sage.misc.persist.load(INTRINSIC_CACHE)
                    return self.__trait_names
                except IOError:
                    pass
            if verbose:
                print "\nCreating list of all MAGMA intrinsics for use in tab completion."
                print "This takes a few minutes the first time, but is saved to the"
                print "file '%s' for future instant use."%INTRINSIC_CACHE
                print "Delete that file to force recreation of this cache."
                print "Scanning MAGMA types ..."
                tm = sage.misc.misc.cputime()
            T = self.eval('ListTypes()').split()
            N = []
            for t in T:
                if verbose:
                    print t, " ",
                    sys.stdout.flush()
                try:
                    s = self.eval('ListSignatures(%s)'%t)
                    for x in s.split('\n'):
                        i = x.find('(')
                        N.append(x[:i])
                except RuntimeError, msg:  # weird internal problems in MAGMA type system
                    print 'Error -- %s'%msg
                    pass
            if verbose:
                print "Done! (%s seconds)"%sage.misc.misc.cputime(tm)
            N = list(set(N))
            N.sort()
            print "Saving cache to '%s' for future instant use."%INTRINSIC_CACHE
            print "Delete the above file to force re-creation of the cache."
            sage.misc.persist.save(N, INTRINSIC_CACHE)
            self.__trait_names = N
            return N


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

    def _sage_doc_(self):
        M = self._obj.parent()
        t = str(self._obj.Type())
        s = M.eval(self._name)
        Z = s.split('(<')[1:]
        W = []
        tt = '(<%s'%t
        for X in Z:
            X = '(<' + X
            if '(<All>' in X or tt in X:
                W.append(X)
        return '\n'.join(W)

    def __repr__(self):
        M = self._obj.parent()
        return M.eval('%s`%s'%(self._obj.name(), self._name))


class MagmaFunction(ExpectFunction):
    def __call__(self, *args, **kwds):
        nvals = 1
        if len(kwds) > 0:
            if kwds.has_key('nvals'):
                nvals = kwds['nvals']
                del kwds['nvals']
        M = self._parent
        return M.function_call(self._name,
                               list(args),
                               params = kwds,
                               nvals = nvals)
    def _sage_doc_(self):
        M = self._parent
        return M.eval(self._name)


class MagmaElement(ExpectElement):
    def __getattr__(self, attrname):
        return MagmaFunctionElement(self, attrname)

    def evaluate(self, *args):
        return ExpectElement.__call__(self, *args)

    def __call__(self, *args):
        """
        EXAMPLES:
            sage: M = magma.RMatrixSpace(magma.IntegerRing(), 2, 2)  # optional
            sage: A = M([1,2,3,4]); A        # optional
            [1 2]
            [3 4]
            sage: type(A)                    # optional
            <class 'sage.interfaces.magma.MagmaElement'>
            sage: A.Type()                   # optional
            ModMatRngElt
        """
        if len(args) > 1:
            return self.evaluate(*args)
        self._check_valid()
        P = self.parent()
        x = P(args[0])
        #try:
        return P('%s!%s'%(self.name(), x.name()))
        #except (RuntimeError, TypeError):
        #    return self.evaluate(*args)

    def set_magma_attribute(self, attrname, value):
        P = self.parent()   # instance of MAGMA that contains this element.
        if not (isinstance(value, MagmaElement) and value.parent() is P):
            value = P(value)
        P.eval('%s`%s := %s'%(self.name(), attrname, value.name()))

    def get_magma_attribute(self, attrname):
        P = self.parent()
        return P('%s`%s'%(self.name(), attrname))

    def list_attributes(self):
        return magma.eval('ListAttributes(Type(%s))'%\
                          self.name()).split()

    def trait_names(self):
        M = self.methods()
        N = []
        for x in M:
            i = x.find('(')
            N.append(x[:i])
        return N + self.list_attributes()

    def methods(self, any=False):
        """
        Return all MAGMA intrinsics that can take self as the first
        argument.

        INPUT:
            any -- (bool: default is False) if True, also include
                   signatures with <Any> as first argument.
        """
        t = str(self.Type())
        X = self.parent().eval('ListSignatures(%s)'%self.Type()).split('\n')
        tt = "(<"+t
        if any:
            Y = [x for x in X if tt in x or "(<Any>" in x]
        else:
            Y = [x for x in X if tt in x]
        return Y


magma = Magma()

def reduce_load_Magma():
    return magma

def magma_console():
    console('magma')

def magma_version():
    t = tuple([int(n) for n in magma.eval('GetVersion()').split()])
    return t, 'V%s.%s-%s'%t
