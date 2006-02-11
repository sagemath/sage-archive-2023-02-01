r"""
Interface to GAP

\sage provides an interface to the GAP system.  This system provides
extensive group theory, combinatorics, etc.

The GAP interface will only work if GAP is installed on your
computer; this should be the case, since GAP is included with
\sage.  The interface offers three pieces of functionality:
\begin{enumerate}
\item \code{gap_console()} -- A function that dumps you
into an interactive command-line GAP session.

\item \code{gap(expr)} -- Evaluation of arbitrary GAP
expressions, with the result returned as a string.

\item \code{gap.new(expr)} -- Creation of a \sage object that wraps a
GAP object.  This provides a Pythonic interface to GAP.  For example,
if \code{f=gap.new(10)}, then \code{f.Factors()} returns the prime
factorization of $10$ computed using GAP.

\end{enumerate}

\subsection{First Examples}

We factor an integer using GAP:

    sage: n = gap(20062006); n
    20062006
    sage: n.parent()
    Gap
    sage: fac = n.Factors(); fac
    [ 2, 17, 59, 73, 137 ]
    sage: fac.parent()
    Gap
    sage: fac[1]
    2

\subsection{GAP and Singular}
This example illustrates conversion between Singular and GAP via \sage
as an intermediate step.  First we create and factor a Singular polynomial.

    sage: singular(389)
    389
    sage: R1 = singular.ring(0, '(x,y)', 'dp')
    sage: R1 = singular.ring(0, '(x,y)', 'dp')
    sage: f = singular('9*x^16-18*x^13*y^2-9*x^12*y^3+9*x^10*y^4-18*x^11*y^2+36*x^8*y^4+18*x^7*y^5-18*x^5*y^6+9*x^6*y^4-18*x^3*y^6-9*x^2*y^7+9*y^8')
    sage: F = f.factorize()
    sage: print F
    [1]:
       _[1]=9
       _[2]=x^6-2*x^3*y^2-x^2*y^3+y^4
       _[3]=-x^5+y^2
    [2]:
       1,1,2

Next we convert the factor $-x^5+y^2$ to a \sage multivariate
polynomial.  Note that it is important to let $x$ and $y$ be the
generators of a polynomial ring, so the eval command works.

    sage: x, y = MPolynomialRing(RationalField(), 2).gens()
    sage: g = eval(F[1][3].sage_polystring()); g
    x_1^2 - x_0^5

Next we create a polynomial ring in GAP and obtain its indeterminates:

    sage: R = gap.PolynomialRing('Rationals', 2); R
    PolynomialRing(..., [ x_1, x_2 ])
    sage: I = R.IndeterminatesOfPolynomialRing(); I
    [ x_1, x_2 ]

In order to eval $g$ in GAP, we need to tell GAP to view the variables
\code{x_0} and \code{x_1} as the two generators of $R$.  This is the
one tricky part.  In the GAP interpreter the object \code{I} has its
    own name (which isn't \code{I}).  We can access its name using
    \code{I.name()}.

        sage: _ = gap.eval("x_0 := %s[1];; x_1 := %s[2];;"%(I.name(), I.name()))

    Now $x_0$ and $x_1$ are defined, so we can construct the GAP polynomial $f$
    corresponding to $g$:

        sage: f = gap(g); f
        -x_1^5+x_2^2

    We can call GAP functions on $f$.  For example, we evaluate
    the GAP \code{Value} function, which evaluates $f$ at the point $(1,2)$.

        sage: f.Value(I, [1,2])
        3
        sage: g(1,2)        # agrees
        3

\subsection{Saving and loading objects} Saving and loading GAP objects
(using the dumps method, etc.)  is \emph{not} supported, since the
output string representation of Gap objects is sometimes not valid
input to GAP.  Creating classes that wrap GAP objects \emph{is}
supported, via simply defining the a _gap_init_ member function that
returns a string that when evaluated in GAP constructs the object.
See \code{groups/permutation_group.py} for a nontrivial example of
this.

\subsection{Long Input}
The GAP interface reads in even very long input (using files) in a
robust manner, as long as you are creating a new object.
\note{Using \code{gap.eval} for long input
is much less robust, and is not recommended.}

    sage: t = '"%s"'%10^10000   # ten thousand character string.
    sage: a = gap(t)

AUTHORS:
    -- David Joyner and William Stein; initial version(s)
    -- William Stein (2006-02-01): modified gap_console command
       so it uses exactly the same startup command as Gap.__init__.

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

from expect import Expect, ExpectElement
from sage.misc.misc import SAGE_ROOT, is_64_bit
import os
DB_HOME = "%s/data/"%SAGE_ROOT
WORKSPACE = "%s/tmp/gap-workspace"%SAGE_ROOT

def gap_command(use_workspace_cache=True):
    if use_workspace_cache and os.path.exists(WORKSPACE):
        return "gap -L %s"%WORKSPACE, False
    else:
        return "gap ", True


class Gap(Expect):
    r"""
    Interface to the GAP interpreter.

    AUTHORS: William Stein and David Joyner
    """
    def __init__(self, max_workspace_size=None,
                 maxread=100000, script_subdirectory="user",
                 use_workspace_cache = True,
                 server=None,
                 logfile = None):

        self.__use_workspace_cache = use_workspace_cache
        cmd, self.__make_workspace = gap_command(use_workspace_cache)
        cmd += ' -T -n -b '
        if max_workspace_size != None:
            cmd += " -o %s"%int(max_workspace_size)
        else: # unlimited
            if is_64_bit:
                cmd += " -o 9999G"
            else:
                cmd += " -o 3900m"

        Expect.__init__(self,
                        name = 'gap',
                        prompt = 'gap> ',
                        command = cmd,
                        maxread = maxread,
                        server = server,
                        script_subdirectory = script_subdirectory,
                        restart_on_ctrlc = True,
                        verbose_start = False,
                        logfile = logfile,
                        eval_using_file_cutoff=100)

        self.__seq = 0

    def __reduce__(self):
        return reduce_load_GAP, tuple([])

    def _quit_string(self):
        return 'quit'

    def _next_var_name(self):
        if len(self._available_vars) != 0:
            v = self._available_vars[0]
            del self._available_vars[0]
            return v
        if self.__seq == 0:
            self.eval('sage := [ ];')
        self.__seq += 1
        return 'sage[%s]'%self.__seq

    #def _read_in_file_command(self, filename):
    #    return 'EvalString(ReadAll(InputTextFile("%s")));'%filename
        #return 'Read(InputTextFile("%s"));'%filename

    def _eval_line_using_file(self, line, tmp):
        F = open(tmp, 'w')
        F.write(line)
        F.close()
        return self._eval_line('Read(InputTextFile("%s"));'%tmp,
                               allow_use_file=False)

    # Change the default for Gap, since eval using a file doesn't
    # work except for setting variables.
    def _eval_line(self, line, allow_use_file=False, wait_for_prompt=True):
        return Expect._eval_line(self, line, allow_use_file=allow_use_file,
                                 wait_for_prompt=wait_for_prompt)

    def _start(self):
        Expect._start(self, "Failed to start GAP.  One possible reason for this is that your gap workspace may be corrupted.  Perhaps remove %s/tmp/gap-workspace"%SAGE_ROOT)
        if self.__use_workspace_cache and self.__make_workspace:
            self.eval('SaveWorkspace("%s");'%WORKSPACE)

    def load_package(self, pkg, verbose=False):
        """
        Load the Gap package with the given name.
        """
        if verbose:
            "Loading GAP package %s"%pkg
        self.eval('LoadPackage("%s")'%pkg)

    def save_workspace(self):
        self.eval('SaveWorkspace("%s");'%WORKSPACE)

    def eval(self, x, newlines=False):
        r"""
        Send the code in the string s to the GAP interpreter and return
        the output as a string.

        INPUT:
            s -- string containing GAP code.
            newlines -- bool (default: True); if False, remove all
                      backslash-newlines inserted by the GAP output formatter.
        """
        # newlines cause hang (i.e., error but no gap> prompt!)
        x = str(x).rstrip().replace('\n',' ')
        if len(x) == 0 or x[len(x) - 1] != ';':
            x += ';'
        s = Expect.eval(self, x)
        if newlines:
            return s
        else:
            return s.replace("\\\n","")

    # Todo -- this -- but there is a tricky "when does it end" issue!
    # Maybe do via a file somehow?
    #def help(self, s):
    #    """
    #    Print help on a given topic.#
    #    """
    #    print Expect.eval(self, "? %s"%s)

    def set(self, var, value):
        """
        Set the variable var to the given value.
        """
        cmd = ('%s:=%s;;'%(var,value)).replace('\n','')
        #out = self.eval(cmd)
        out = self._eval_line(cmd, allow_use_file=True)
        if out.lower().find('error') != -1:
            raise TypeError, "Error executing code in GAP\nCODE:\n\t%s\nGAP ERROR:\n\t%s"%(cmd, out)


    def get(self, var):
        """
        Get the value of the variable var.
        """
        # TODO: Steve Linton says -- use "Print()".
        return self.eval('%s;'%var, newlines=False)

    #def clear(self, var):
        #"""
        #Clear the variable named var.
        #"""
        #self.eval('Unbind(%s)'%var)
        #self._available_vars.append(var)

    def _contains(self, v1, v2):
        return self.eval('%s in %s'%(v1,v2))

    def _is_true_string(self, t):
        return t == "true"

    def _true_symbol(self):
        return "true"

    def _false_symbol(self):
        return "false"

    def _equality_symbol(self):
        return "="

    def console(self):
        gap_console()

    def version(self):
        return gap_version()

    def _object_class(self):
        return GapElement


############

def gap_reset_workspace(max_workspace_size=None):
    """
    Call this to completely reset the GAP workspace, which
    is used by default when SAGE first starts GAP.

    The first time you start GAP from SAGE, it saves the
    startup state of GAP in the file

        tmp/gap-workspace

    This is useful, since then subsequent startup of GAP
    is at least 10 times as fast.  Unfortunately, if you
    install any new code for GAP, it won't be noticed unless
    you explicitly load it, e.g., with
           gap.load_package("laguna")
    """
    g = Gap(use_workspace_cache=False, max_workspace_size=None)
    g.eval('SaveWorkspace("%s");'%WORKSPACE)


class GapElement(ExpectElement):
    def __getitem__(self, n):
        self._check_valid()
        if not isinstance(n, tuple):
            return self.parent().new('%s[%s]'%(self._name, n))
        else:
            return self.parent().new('%s%s'%(self._name, ''.join(['[%s]'%x for x in n])))

    def __reduce__(self):
        return reduce_load, ()  # default is an invalid object

    def __repr__(self):
        s = ExpectElement.__repr__(self)
        if s.find('must have a value') != -1:
            raise RuntimeError, "An error occured creating an object in %s from:\n'%s'\n%s"%(self.parent().name(), self._create, s)
        return s

    def __len__(self):
        return int(self.Length())

    def _matrix_(self, R):
        r"""
        Return matrix over the (\sage) ring R determined by self, where self
        should be a Gap matrix.

            sage: s = gap("Z(3)*[[1,2,3],[3,4,5]]"); s
            [ [ Z(3), Z(3)^0, 0*Z(3) ], [ 0*Z(3), Z(3), Z(3)^0 ] ]
            sage: matrix(s, GF(3))
            [2 0 1]
            [2 0 1]

            sage: s = gap("[[1,2], [3/4, 5/6]]"); s
            [ [ 1, 2 ], [ 3/4, 5/6 ] ]
            sage: m = matrix(s, QQ); m
            [  1 3/4]
            [  2 5/6]
            sage: parent(m)
            Full MatrixSpace of 2 by 2 dense matrices over Rational Field

            sage: s = gap('[[Z(16),Z(16)^2],[Z(16)^3,Z(16)]]')
            sage: matrix(s, GF(16))
            [  a a^3]
            [a^2   a]
        """
        P = self.parent()
        v = self.DimensionsMat()
        n = int(v[1])
        m = int(v[2])

        from sage.matrix.matrix_space import MatrixSpace
        M = MatrixSpace(R, n, m)
        entries = [[R(self[i,j]) for i in range(1,n+1)] for j in range(1,m+1)]

        return M(entries)

        mat = copy.copy(s)
        mat = mat.replace("\n","")
        mat = mat.replace("] ]","")
        mat = mat.replace("[ [","")
        mat = mat.split("],")
        rows_str = mat
        rowss = []
        for x in rows_str:
            rowss.append(x.split(","))
        rows_mat = [[] for i in range(len(rowss))]
        for i in range(len(rowss)):
            for x in rowss[i]:
                x = x.replace("]","")
                x = x.replace("[","")
                rows_mat[i].append(x)
        M = []
        for i in range(m):
            rows = []
            for x in rows_mat[i]:
                rows.append(gap2sage_finite_field(x,F))
            M.append(rows)
        MS = MatrixSpace(F,m,n)
        return MS(M)



def is_GapElement(x):
    return isinstance(x, GapElement)

###########

#############

gap = Gap()

def reduce_load_GAP():
    return gap

def reduce_load():
    return GapElement(None, None)

import os
def gap_console(use_workspace_cache=True):
    cmd, _ = gap_command(use_workspace_cache=use_workspace_cache)
    os.system(cmd)

def gap_version():
    return gap.eval('VERSION')[1:-1]


