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

We verify that an obviously principal ideal is principal:

    sage: _ = magma.eval('R<x> := PolynomialRing(RationalField())')
    sage: O = magma.NumberField('x^2+23').MaximalOrder()
    sage: I = magma('ideal<%s|%s.1>'%(O.name(),O.name()))
    sage: I.IsPrincipal(nvals=2)
    (true, [1, 0])

\subsection{Long Input}
The Magma interface reads in even very long input (using files) in a
robust manner.

    sage: t = '"%s"'%10^10000   # ten thousand character string.
    sage: a = magma.eval(t)
    sage: a = magma(t)

\subsection{Other Examples}

We compute a space of modular forms with character.

    sage: N = 20
    sage: D = 20
    sage: eps_top = fundamental_discriminant(D)
    sage: eps = magma.KroneckerCharacter(eps_top, RationalField())
    sage: M2 = magma.ModularForms(eps)
    sage: print M2
    Space of modular forms on Gamma_1(5) with character $.1, weight 2, and dimension 2 over Integer Ring.
    sage: print M2.Basis()   # note -- this has been changed to be *wrong* as below in Magma 2.14!!
    [
    1 + 10*q^2 + 20*q^3 + 20*q^5 + 60*q^7 + 50*q^8 + 30*q^10 + O(q^12),
    q + q^2 + 2*q^3 + 3*q^4 + 5*q^5 + 2*q^6 + 6*q^7 + 5*q^8 + 7*q^9 + 5*q^10 + 12*q^11 + O(q^12)
    ]

In SAGE/Python (and sort of C++) coercion of an element x into a
structure S is denoted by S(x).  This also works for the MAGMA interface:

    sage: G = magma.DirichletGroup(20)
    sage: G.AssignNames(['a', 'b'])
    sage: (G.1).Modulus()
    20
    sage: e = magma.DirichletGroup(40)(G.1)
    sage: print e
    $.1
    sage: print e.Modulus()
    40

We coerce some polynomial rings into MAGMA:

    sage: R.<y> = PolynomialRing(QQ)
    sage: S = magma(R)
    sage: print S
    Univariate Polynomial Ring in y over Rational Field
    sage: S.1
    y

This example illustrates that SAGE doesn't magically extend how MAGMA
implicit coercion (what there is, at least) works:
    sage: R.<x> = ZZ[]
    sage: x * 5
    5*x
    sage: x * 1.0
    1.00000000000000*x
    sage: x * (2/3)
    2/3*x
    sage: y = magma(x)
    sage: y * 5
    5*x
    sage: y * 1.0
    Traceback (most recent call last):
    ...
    TypeError: unsupported operand parent(s) for '*': 'Magma' and 'Real Field with 53 bits of precision'
    sage: y * (2/3)
    Traceback (most recent call last):
    ...
    TypeError: unsupported operand parent(s) for '*': 'Magma' and 'Rational Field'


AUTHOR:
    -- William Stein (2005): initial version
    -- William Stein (2006-02-28): added extensive tab completion and interactive
                                   IPython documentation support.
    -- William Stein (2006-03-09): added nvals argument for magma.functions...
"""

#*****************************************************************************
#       Copyright (C) 2005 William Stein <wstein@gmail.com>
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

import os, sys

from sage.structure.element import RingElement
from expect import console, Expect, ExpectElement, ExpectFunction, FunctionElement
PROMPT = ">>>"

import sage.misc.misc
import sage.misc.sage_eval

INTRINSIC_CACHE = '%s/magma_intrinsic_cache.sobj'%sage.misc.misc.DOT_SAGE
MAGMA_SPEC = '%s/magma/spec'%sage.misc.misc.SAGE_EXTCODE



class Magma(Expect):
    r"""
    Interface to the Magma interpreter.

    Type \code{magma.[tab]} for a list of all the functions available
    from your Magma install.  Type \code{magma.[tab]?} for Magma's
    help about a given function.  Type \code{magma(...)} to create
    a new Magma object, and \code{magma.eval(...)} to run a string
    using Magma (and get the result back as a string).

    NOTE: If you do not own a local copy of MAGMA, try using the
    \code{magma\_free} command instead, which uses the free demo web
    interface to MAGMA.

    EXAMPLES:

    You must use nvals = 0 to call a function that doesn't return
    anything, otherwise you'll get an error.  (nvals is the number
    of return values.)

        sage: magma.SetDefaultRealFieldPrecision(200, nvals=0)  # optional and requires MAGMA >= v2.12
        sage: magma.eval('1.1')   # optional
        '1.1000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000'
        sage: magma.SetDefaultRealFieldPrecision(30, nvals=0)  # optional
    """
    def __init__(self, maxread=10000, script_subdirectory=None,
                 logfile=None, server=None, server_tmpdir=None, user_config=False):
        """
        INPUT:
            maxread -- affects buffering
            script_subdirectory -- directory where scripts are read from
            logfile -- output logged to this file
            server -- address of remote server
            user_config -- if True, then local user configuration files
                           will be read by MAGMA.  If False (the default),
                           then MAGMA is started with the -n option which
                           supresses user configuration files.
        """
        #command = 'magma -b '
        command = 'magma'
        if not user_config:
            command += ' -n'
        Expect.__init__(self,
                        name = "magma",
                        prompt = ">>SAGE>>",
                        command = command,
                        maxread = maxread,
                        server = server,
                        server_tmpdir = server_tmpdir,
                        script_subdirectory = script_subdirectory,
                        restart_on_ctrlc = True,
                        logfile = logfile,
                        eval_using_file_cutoff=100)
        # We use "-n" above in the Magma startup command so
        # local user startup configuration is not read.

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

    def chdir(self, dir):
        """
        Change the Magma interpreters current working directory.

        EXAMPLES:
            sage: magma.eval('System("pwd")')   # optional and random
            '/Users/was/s/devel/sage-main/sage'
            sage: magma.chdir('..')             # optional
            sage: magma.eval('System("pwd")')   # optional and random
            '/Users/was/s/devel/sage-main'
        """
        self.eval('ChangeDirectory("%s")'%dir, strip=False)

    def eval(self, x, strip=True):
        """
        INPUT:
            x -- string of code
            strip -- ignored
        """
        x = str(x).rstrip()
        if len(x) == 0 or x[len(x) - 1] != ';':
            x += ';'
        ans = Expect.eval(self, x).replace('\\\n','')
        if 'Runtime error' in ans or 'User error' in ans:
            raise RuntimeError, "Error evaluation Magma code.\nIN:%s\nOUT:%s"%(x, ans)
        return ans

    def _start(self):
        self._change_prompt('>')
        Expect._start(self)
        self.eval('SetPrompt("%s"); SetLineEditor(false); SetColumns(0);'%PROMPT)
        self._change_prompt(PROMPT)
        self.expect().expect(PROMPT)
        self.expect().expect(PROMPT)
        self.expect().expect(PROMPT)
        if os.path.exists(MAGMA_SPEC):
            self.attach_spec(MAGMA_SPEC)


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
            if isinstance(x, bool):
                return Expect.__call__(self, str(x).lower())
            return Expect.__call__(self, x)
        return self.objgens(x, gens)


    #def clear(self, var):
    #    """
    #    Clear the variable named var.
    #    """
        #self.eval("delete %s"%var)
    #    self.eval("%s:=0"%var)

    def cputime(self, t=None):
        if t:
            return float(self.eval('Cputime(%s)'%t))
        else:
            return float(self.eval('Cputime()'))

    def chdir(self, dir):
        """
        Change to the given directory.
        """
        magma.eval('ChangeDirectory("%s")'%dir)

    def attach(self, filename):
        r"""
        Attach the given file to the running instance of MAGMA.

        Attaching a file in MAGMA makes all intrinsics defined in the
        file available to the shell.  Moreover, if the file doesn't
        start with the \code{freeze;} command, then the file is
        reloaded whenever it is changed.  Note that functions and
        procedures defined in the file are \emph{not} available.
        For only those, use \code{magma.load(filename)}.
        """
        return self.eval('Attach("%s")'%filename)

    Attach = attach

    def attach_spec(self, filename):
        r"""
        Attach the given spec file to the running instance of MAGMA.

        This attaches numerous files to the running MAGMA (see the
        MAGMA documentation for more details).
        """
        return self.eval('AttachSpec("%s")'%filename)

    AttachSpec = attach_spec

    def load(self, filename):
        """
        Load the file with given filename using the 'load' command
        in the MAGMA shell.

        Loading a file in MAGMA makes all the functions and procedures
        in the file available. The file should not contain any
        intrinsics (or you'll get errors).
        """
        return self.eval('load "%s"'%filename)

    def _next_var_name(self):
        if self.__seq == 0:
            self.eval('_sage_ := [* *];')
        else:
            try:
                self.eval('Append(~_sage_, 0);')
            except:
                # this exception could happen if the MAGMA process
                # was interrupted during startup / initialization.
                self.eval('_sage_ := [* 0 : i in [1..%s] *];'%self.__seq)
        self.__seq += 1
        return '_sage_[%s]'%self.__seq

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
            par = ' : ' + ','.join(['%s:=%s'%(a,self(b)) for a,b in params.items()])

        fun = "%s(%s%s)"%(function, ",".join([s.name() for s in args]), par)

        return self._do_call(fun, nvals)

    def _do_call(self, code, nvals):
        if nvals <= 0:
            out = self.eval(code)
            ans = None
        elif nvals == 1:
            return self(code)
        else:
            v = [self._next_var_name() for _ in range(nvals)]
            vars = ", ".join(v)
            cmd = "%s := %s;"%(vars, code)
            out = self.eval(cmd)
            ans = tuple([MagmaElement(self, x, is_name = True) for x in v])

        if out.lower().find("error") != -1:
            raise TypeError, "Error executing Magma code:\n%s"%out
        return ans

    def bar_call(self, left, name, gens, nvals=1):
        """
        This is a wrapper around the Magma constructor

           name<left | gens>

        returning nvals.

        INPUT:
            left -- something coerceable to a magma object
            name -- name of the constructor, e.g., sub, quo, ideal, etc.
            gens -- if a list/tuple, each item is coerced to magma;
                    otherwise gens itself is converted to magma
            nvals -- positive integer; number of return values

        OUTPUT:
            a single magma object if nvals == 1; otherwise a tuple
            of nvals magma objects.

        EXAMPLES:
        The bar_call function is used by the sub, quo, and ideal
        methods of Magma elements.  Here we illustrate directly using
        bar_call to create quotients:

            sage: V = magma.RModule(ZZ,3)    # optional -- requires magma
            sage: V                          # optional
            RModule(IntegerRing(), 3)
            sage: magma.bar_call(V, 'quo', [[1,2,3]], nvals=1)  # optional
            RModule(IntegerRing(), 2)
            sage: magma.bar_call(V, 'quo', [[1,2,3]], nvals=2)  # optional
            (RModule(IntegerRing(), 2),
             Mapping from: RModule(IntegerRing(), 3) to RModule(IntegerRing(), 2))
            sage: magma.bar_call(V, 'quo', V, nvals=2)          # optional
            (RModule(IntegerRing(), 0),
             Mapping from: RModule(IntegerRing(), 3) to RModule(IntegerRing(), 0))
        """
        magma = self
        # coerce each arg to be a Magma element
        if isinstance(gens, (list, tuple)):
            gens = [magma(z) for z in gens]
            # make comma separated list of names (in Magma) of each of the gens
            v = ', '.join([w.name() for w in gens])
        else:
            gens = magma(gens)
            v = gens.name()
        # construct the string that evaluates in Magma to define the subobject,
        # and return it evaluated in Magma.
        s = '%s< %s | %s >'%(name, left.name(), v)
        return self._do_call(s, nvals)

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

    def _lessthan_symbol(self):
        return ' lt '

    def _greaterthan_symbol(self):
        return ' gt '

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
                print "MAGMA may produce errors during this process, which are safe to ignore."
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

    def ideal(self, L):
        """
        Return the MAGMA ideal defined by L.

        INPUT:
            L -- a list of elements of a SAGE multivariate polynomial ring.

        OUTPUT:
            The magma ideal generated by the elements of L.

        EXAMPLES:
            sage: R.<x,y> = QQ[]
            sage: magma.ideal([x^2, y^3*x])         # optional -- requires magma
            Ideal of Polynomial ring of rank 2 over Rational Field
            Graded Reverse Lexicographical Order
            Variables: x, y
            Basis:
            [
            x^2,
            x*y^3
            ]
        """
        P = iter(L).next().parent()
        Pn = P._magma_().name()
        k = P.base_ring()
        if k.degree() > 1:
            i = str(k.gen())
            o = self("BaseRing(%s).1"%Pn).name()
            self.eval("%s := %s"%(i,o))
        mlist = self(L)
        return self("ideal<%s|%s>"%(Pn,mlist.name()))

    def set_verbose(self, type, level):
        """
        Set the verbosity level for a given algorithm, class, etc. in
        MAGMA.

        INPUT:
            type -- string (e.g. 'Groebner')
            level -- integer >= 0

        """
        self.SetVerbose(type,level)

    def SetVerbose(self, type, level):
        """
        Set the verbosity level for a given algorithm class etc. in
        MAGMA.

        INPUT:
            type -- string (e.g. 'Groebner'), see MAGMA documentation
            level -- integer >= 0

        NOTE: This method is provided to be consistent with the MAGMA
        naming convention.
        """
        if level < 0:
            raise TypeError, "level must be >= 0"
        self.eval('SetVerbose("%s",%d)'%(type,level))

    def get_verbose(self, type):
        """
        Get the verbosity level of a given algorithm class etc. in
        MAGMA.

        INPUT:
            type -- string (e.g. 'Groebner'), see MAGMA documentation
        """
        return self.GetVerbose(type)

    def GetVerbose(self, type):
        """
        Get the verbosity level of a given algorithm class etc. in
        MAGMA.

        INPUT:
            type -- string (e.g. 'Groebner'), see MAGMA documentation

        NOTE: This method is provided to be consistent with the MAGMA
        naming convention.
        """
        return int(self.eval('GetVerbose("%s")'%type))

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
        s = '\n'.join(W)
        s = sage.misc.misc.word_wrap(s)
        return s

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
        s = M.eval(self._name)
        s = sage.misc.misc.word_wrap(s, 80)
        return s


def is_MagmaElement(x):
    return isinstance(x, MagmaElement)

class MagmaElement(ExpectElement):
    def __getattr__(self, attrname):
        if attrname[:1] == "_":
            raise AttributeError
        return MagmaFunctionElement(self, attrname)

    def _sage_(self):
        """
        Return Sage version of this object.  Use self.sage() to
        get the Sage version.

        EXAMPLES:
        Enumerated Sets:
            sage: a = magma('{1,2/3,-5/9}')       # optional
            sage: a.sage()                        # optional
            {1, -5/9, 2/3}
            sage: type(a.sage())                  # optional
            <class 'sage.sets.set.Set_object_enumerated'>
            sage: a = magma('{1,2/3,-5/9}'); a    # optional
            { -5/9, 2/3, 1 }
            sage: a.Type()                        # optional
            SetEnum
            sage: b = a.sage(); b             # optional
            {1, -5/9, 2/3}
            sage: type(b)                         # optional
            <class 'sage.sets.set.Set_object_enumerated'>
            sage: c = magma(b); c                 # optional
            { -5/9, 2/3, 1 }
            sage: c.Type()                        # optional
            SetEnum

        Multisets are converted to lists:
            sage: m = magma('{* 1,2,2,2,4^^2,3 *}')    # optional
            sage: z = m.sage(); z                      # optional
            [1, 2, 2, 2, 3, 4, 4]
            sage: type(z)                              # optional
            <type 'list'>

        Matrices:
            sage: a = matrix(ZZ,3,3,[1..9])
            sage: m = magma(a)                        # optional
            sage: b = m.sage(); b                     # optional
            [1 2 3]
            [4 5 6]
            [7 8 9]
            sage: b == a                             # optional
            True

        A nonsquare matrix:
            sage: a = matrix(ZZ,2,3,[1..6])
            sage: m = magma(a)                       # optional
            sage: m.sage()                           # optional
            [1 2 3]
            [4 5 6]
        """
        P = self._check_valid()
        cmd = "Sage(%s)"%self.name()
        s = str(P.eval(cmd))
        return sage.misc.sage_eval.sage_eval(s)

    def AssignNames(self, names):
        """
        EXAMPLES:
            sage: G = magma.DirichletGroup(20)   # optional
            sage: G.AssignNames(['a','b'])       # optional
            sage: G.1                            # optional
            a

            sage: G.Elements()                   # optional
            [
            1,
            a,
            b,
            a*b
            ]
        """
        P = self._check_valid()
        cmd = 'AssignNames(~%s, [%s])'%(self.name(),
                                        ','.join('"%s"'%x for x in names))
        P.eval(cmd)

    assign_names = AssignNames

    def gens(self):
        """
        Return generators for self.

        If self is named X is MAGMA, this function evaluates X.1, X.2,
        etc., in MAGMA until an error occurs.  It then returns a SAGE
        list of the resulting X.i.  Note -- I don't think there is a
        MAGMA command that returns the list of valid X.i.  There are
        numerous ad hoc functions for various classes but nothing
        systematic.  This function gets around that problem.

        AUTHOR:
            * William Stein -- 2006-07-02
        """
        try:
            return self.__gens
        except AttributeError:
            pass
        G = []
        i = 1
        P = self._check_valid()
        n = self.name()
        while True:
            try:
                G.append(P('%s.%s'%(n,i)))
            except (RuntimeError, TypeError):
                break
            i += 1
        self.__gens = G
        return G

    def evaluate(self, *args):
        P = self._check_valid()
        v = [P(a) for a in args]
        names = ','.join([str(x) for x in v])
        return P('%s(%s)'%(self.name(), names))
    eval = evaluate

    def __call__(self, *args):
        """
        Coerce something into the object (using the MAGMA ! notation).

        For function calls, use self.eval(...).

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
        P = self._check_valid()
        x = P(args[0])
        try:
            return P('%s!%s'%(self.name(), x.name()))
        except (RuntimeError, TypeError):
            return self.evaluate(*args)

    def __iter__(self):
        P = self._check_valid()
        i = 1
        while True:
            try:
                yield self[i]
            except (TypeError, RuntimeError):
                return
            i += 1

    def __len__(self):
        P = self._check_valid()
        return int(P.eval('#%s'%self.name()))

    def _polynomial_(self, R):
        return R(list(self.Eltseq()))

    def _latex_(self):
        r"""
        Return latex representation of self.

        AUTHOR:
            -- Jennifer Balakrishnan (jenb@mit.edu)

        Types that are nicely latex include:
        \begin{itemize}
          \item rationals
          \item matrices
          \item polynomials
          \item binary quadratic forms
           \item elements of quadratic, cyclotomic number fields, and general
           number fields
          \item points
          \item elliptic curves
          \item power series
        \end{itemize}

        IMPLEMENTATION:
            Calls latex.m, which is in SAGE_ROOT/data/extcode/magma/latex.m

        EXAMPLES:
            sage: latex(magma('-2/3'))                            # optional
            \frac{-2}{3}

            sage: magma.eval('R<x> := PolynomialRing(RationalField()); f := (x-17/2)^3;')     # optional
            ''
            sage: latex(magma('f'))                               # optional
            x^{3}-\frac{51}{2}x^{2}+\frac{867}{4}x-\frac{4913}{8}

            sage: latex(magma('(MatrixAlgebra(RationalField(),3)![0,2,3,4,5,6,7,8,9])^(-1)'))    # optional
            \left(\begin{array}{ccc}-1&2&-1\\2&-7&4\\-1&\frac{14}{3}&\frac{-8}{3}\end{array}\right)

            sage: magma.eval('K<a> := CyclotomicField(11)')       # optional
            ''
            sage: latex(magma('a^3 + a - 17/3'))                  # optional
            \frac{-17}{3}+\zeta_{11}+\zeta_{11}^{3}

            sage: latex(magma('EllipticCurve([1,2/3,3/4,4/5,-5/6])'))    # optional
            y^2+xy+\frac{3}{4}y=x^3+\frac{2}{3}x^2+\frac{4}{5}x-\frac{5}{6}


            sage: _=magma.eval('R<x> := PolynomialRing(RationalField())')    # optional
            sage: _=magma.eval('K<a> := NumberField(x^3+17*x+2)')            # optional
            sage: latex(magma('(1/3)*a^2 - 17/3*a + 2'))                     # optional
            2-\frac{17}{3}a+\frac{1}{3}a^{2}

        SAGE auto-detects the greek letters and puts backslashes in:
            sage: _=magma.eval('R<x> := PolynomialRing(RationalField())')    # optional
            sage: _=magma.eval('K<alpha> := NumberField(x^3+17*x+2)')        # optional
            sage: latex(magma('(1/3)*alpha^2 - 17/3*alpha + 2'))             # optional
            2-\frac{17}{3}\alpha+\frac{1}{3}\alpha^{2}

            sage: _=magma.eval('R<alpha> := PolynomialRing(RationalField())') # optional
            sage: latex(magma('alpha^3-1/7*alpha + 3'))                      # optional
            \alpha^{3}-\frac{1}{7}\alpha+3


        Finite field elements:
            sage: _=magma.eval('K<a> := GF(27)')                             # optional
            sage: latex(magma('a^2+2'))                                      # optional
            2+a^{2}

        Printing of unnamed (dollar sign) generators works correctly:
            sage: latex(magma('FiniteField(81).1^2+1'))                      # optional
            1+\$.1^{2}

        Finite fields:
            sage: latex(magma('FiniteField(3)'))                             # optional
            \mathbf{F}_{{3}}
            sage: latex(magma('FiniteField(27)'))                            # optional
            \mathbf{F}_{{3}^{3}}

        Power Series:
            sage: _=magma.eval('R<x> := PowerSeriesRing(RationalField())')   # optional
            sage: latex(magma('(1/(1+x))'))                                  # optional
            1-x+x^{2}-x^{3}+x^{4}-x^{5}+x^{6}-x^{7}+x^{8}-x^{9}+x^{10}-x^{11}+x^{12}-x^{13}+x^{14}-x^{15}+x^{16}-x^{17}+x^{18}-x^{19}+O(x^{20})
            sage: _=magma.eval('R<x> := PowerSeriesRing(RationalField())')   # optional
            sage: latex(magma('(-1/(2+x + O(x^3)))'))                        # optional
            \frac{-1}{2}+\frac{1}{4}x-\frac{1}{8}x^{2}+O(x^{3})

        p-adic Numbers:
            sage: latex(magma('pAdicField(7,4)!9333294394/49'))              # optional
            4\cdot{}7^{-2} + 5\cdot{}7^{-1} + 5+ 6\cdot{}7^{1} + O(7^{2})
        """
        P = self._check_valid()
        s = str(P.eval('Latex(%s)'%self.name()))
        v = '\\mbox{\\rm '
        if s[:len(v)] == v:
            raise AttributeError
        return s

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
        v = list(set(N + self.list_attributes()))
        v.sort()
        return v

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

    def __floordiv__(self, x):
        """
        Quotient of division of self by other.  This is denoted // ("div" in magma).

        EXAMPLE:
            sage: R.<x,y,z>=QQ[]
            sage: magma(5)//magma(2) # optional
            2
            sage: m=magma(x*z+x*y)   # optional
            sage: n=magma(x)         # optional
            sage: m//n               # optional
            y + z
        """
        return self.parent()('%s div %s'%(self.name(), x.name()))

    def __nonzero__(self):
        try:
            return not self.parent()("%s eq 0"%self.name()).bool()
        except TypeError:
            return self.bool()

    def sub(self, gens):
        """
        Return the sub-object of self with given gens.

        INPUT:
            gens -- object or list/tuple of generators

        EXAMPLES:
            sage: V = magma('VectorSpace(RationalField(),3)')       # optional -- requires magma
            sage: W = V.sub([ [1,2,3], [1,1,2] ]); W                # optional
            Vector space of degree 3, dimension 2 over Rational Field
            Generators:
            (1 2 3)
            (1 1 2)
            Echelonized basis:
            (1 0 1)
            (0 1 1)
        """
        return self.parent().bar_call(self, 'sub', gens)

    def quo(self, gens):
        """
        Return the quotient of self by the given object or list of generators.

        INPUT:
            gens -- object or list/tuple of generators
        OUTPUT:
            magma element -- the quotient object
            magma element -- mapping from self to the quotient object

        EXAMPLES:
            sage: V = magma('VectorSpace(RationalField(),3)')       # optional -- requires magma
            sage: V.quo([[1,2,3], [1,1,2]])                         # optional
            (Full Vector space of degree 1 over Rational Field, Mapping from: Full Vector space of degree 3 over Rational Field to Full Vector space of degree 1 over Rational Field)

        We illustrate quotienting out by an object instead of a list of generators:
            sage: W = V.sub([ [1,2,3], [1,1,2] ])                   # optional
            sage: V.quo(W)                                          # optional
            (Full Vector space of degree 1 over Rational Field, Mapping from: Full Vector space of degree 3 over Rational Field to Full Vector space of degree 1 over Rational Field)

        We quotient a ZZ module out by a submodule.
            sage: V = magma.RModule(ZZ,3); V
            RModule(IntegerRing(), 3)
            sage: W, phi = V.quo([[1,2,3]])
            sage: W
            RModule(IntegerRing(), 2)
            sage: phi
            Mapping from: RModule(IntegerRing(), 3) to RModule(IntegerRing(), 2)
        """
        return self.parent().bar_call(self, 'quo', gens, nvals=2)

    def ideal(self, gens):
        """
        Return the ideal of self with given list of generators.

        INPUT:
            gens -- object or list/tuple of generators
        OUTPUT:
            magma element -- a Magma ideal

        EXAMPLES:
            sage: R = magma('PolynomialRing(RationalField())')        # optional -- requires magma
            sage: R.assign_names(['x'])
            sage: x = R.1                                             # optional
            sage: R.ideal([x^2 - 1, x^3 - 1])                         # optional
            Ideal of Univariate Polynomial Ring in x over Rational Field generated by x - 1
        """
        return self.parent().bar_call(self, 'ideal', gens, nvals=1)

###########################################################################

magma = Magma()

def reduce_load_Magma():
    return magma

def magma_console():
    console('magma')

def magma_version():
    t = tuple([int(n) for n in magma.eval('GetVersion()').split()])
    return t, 'V%s.%s-%s'%t
