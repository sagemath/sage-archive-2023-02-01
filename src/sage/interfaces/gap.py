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
AP object.  This provides a Pythonic interface to GAP.  For example,
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
    x1^2 - x0^5

Next we create a polynomial ring in GAP and obtain its indeterminates:

    sage: R = gap.PolynomialRing('Rationals', 2); R
    PolynomialRing(..., [ x_1, x_2 ])
    sage: I = R.IndeterminatesOfPolynomialRing(); I
    [ x_1, x_2 ]

In order to eval $g$ in GAP, we need to tell GAP to view the variables
\code{x0} and \code{x1} as the two generators of $R$.  This is the
one tricky part.  In the GAP interpreter the object \code{I} has its
    own name (which isn't \code{I}).  We can access its name using
    \code{I.name()}.

        sage: _ = gap.eval("x0 := %s[1];; x1 := %s[2];;"%(I.name(), I.name()))

    Now $x_0$ and $x_1$ are defined, so we can construct the GAP polynomial $f$
    corresponding to $g$:

        sage: f = gap(str(g)); f
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


\subsection{Changing which GAP is used}
Use this code to change which GAP interpreter is
run.   E.g.,
\begin{verbatim}
   import sage.interfaces.gap
   sage.interfaces.gap.gap_cmd = "/usr/local/bin/gap"
\end{verbatim}

AUTHORS:
    -- David Joyner and William Stein; initial version(s)
    -- William Stein (2006-02-01): modified gap_console command
       so it uses exactly the same startup command as Gap.__init__.
    -- William Stein (2006-03-02): added tab completions:
             gap.[tab], x = gap(...), x.[tab], and docs, e.g.,
                  gap.function?  and x.function?
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

from expect import Expect, ExpectElement, FunctionElement, ExpectFunction, tmp
from sage.misc.misc import SAGE_ROOT, DOT_SAGE, is_64_bit
from IPython.genutils import page
import re
import os
import pexpect
DB_HOME = "%s/data/"%SAGE_ROOT
WORKSPACE = "%s/gap-workspace"%DOT_SAGE

if not os.path.exists('%s/tmp'%SAGE_ROOT):
    os.makedirs('%s/tmp'%SAGE_ROOT)

gap_cmd = "gap"

def gap_command(use_workspace_cache=True):
    if use_workspace_cache and os.path.exists(WORKSPACE):
        return "%s -L %s"%(gap_cmd, WORKSPACE), False
    else:
        return gap_cmd, True


class Gap(Expect):
    r"""
    Interface to the GAP interpreter.

    AUTHORS: William Stein and David Joyner
    """
    def __init__(self, max_workspace_size=None,
                 maxread=100000, script_subdirectory=None,
                 use_workspace_cache = True,
                 server=None,
                 logfile = None):

        self.__use_workspace_cache = use_workspace_cache
        cmd, self.__make_workspace = gap_command(use_workspace_cache)
        cmd += " -b -p -T"
        if max_workspace_size != None:
            cmd += " -o %s"%int(max_workspace_size)
        else: # unlimited
            if is_64_bit:
                cmd += " -o 9999G"
            else:
                cmd += " -o 3900m"
        cmd += " %s/extcode/gap/sage.g"%DB_HOME
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
        #if self.__seq == 0:
        #    self.eval('sage := [ ];')
        self.__seq += 1
        return '$sage%s'%self.__seq

    def _read_in_file_command(self, filename):
        return 'Read("%s");'%filename



    def _start(self):
        Expect._start(self, "Failed to start GAP.  One possible reason for this is that your gap workspace may be corrupted.  Perhaps remove %s/tmp/gap-workspace"%SAGE_ROOT)
        if self.__use_workspace_cache and self.__make_workspace:
            self.eval('SaveWorkspace("%s");'%WORKSPACE)

    def _continuation_prompt(self):
        return '> '

    def __getattr__(self, attrname):
        if attrname[:1] == "_":
            raise AttributeError
        return GapFunction(self, attrname)

    def load_package(self, pkg, verbose=False):
        """
        Load the Gap package with the given name.

        If loading fails, raise a RuntimeError exception.
        """
        if verbose:
            "Loading GAP package %s"%pkg
        x = self.eval('LoadPackage("%s")'%pkg)
        if x == 'fail':
            raise RuntimeError, 'Error loading Gap package %s'%pkg

    def save_workspace(self):
        self.eval('SaveWorkspace("%s");'%WORKSPACE)

    def eval(self, x, newlines=False, strip=True):
        r"""
        Send the code in the string s to the GAP interpreter and return
        the output as a string.

        INPUT:
            s -- string containing GAP code.
            newlines -- bool (default: True); if False, remove all
                      backslash-newlines inserted by the GAP output formatter.
            strip -- ignored
        """
        # newlines cause hang (i.e., error but no gap> prompt!)
        x = str(x).rstrip().replace('\n','')
        if len(x) == 0 or x[len(x) - 1] != ';':
            x += ';'
        s = Expect.eval(self, x)
        if newlines:
            return s
        else:
            return s.replace("\\\n","")

    # Todo -- this -- but there is a tricky "when does it end" issue!
    # Maybe do via a file somehow?
    def help(self, s, pager=True):
        """
        Print help on a given topic.
        """
        self.eval('$SAGE.tempfile := "%s";'%tmp)
        line = Expect.eval(self, "? %s"%s)
        match = re.search("Page from (\d+)", line)
        if match == None:
            print line
        else:
            (sline,) = match.groups()
            F = open(tmp,"r")
            if pager:
                page(F.read(), start = int(sline)-1)
            else:
                return F.read()

    def set(self, var, value):
        """
        Set the variable var to the given value.
        """
        cmd = ('%s:=%s;;'%(var,value)).replace('\n','')
        #out = self.eval(cmd)
        out = self._eval_line(cmd, allow_use_file=True)
        #if out.lower().find('error') != -1:
        #    raise TypeError, "Error executing code in GAP\nCODE:\n\t%s\nGAP ERROR:\n\t%s"%(cmd, out)


    def get(self, var):
        """
        Get the value of the variable var.
        """
        return self.eval('%s;'%var, newlines=False)

    def __getattr__(self, attrname):
        if attrname[:1] == "_":
            raise AttributeError
        return GapFunction(self, attrname)

    def _pre_interact(self):
        self._eval_line("$SAGE.StartInteract();")

    def _post_interact(self):
        self._eval_line("$SAGE.StopInteract();")

    def _execute_line(self, line, wait_for_prompt=True, expect_eof=False):
        E = self._expect
        try:
            if len(line) > 4095:
                raise RuntimeError,"Passing commands this long to gap would hang"
            E.sendline(line)
        except OSError:
            return RuntimeError, "Error evaluating %s in %s"%(line, self)
        if wait_for_prompt == False:
            return ('','')
        if len(line)==0:
            return ('','')
        try:
            E.expect("\r\n") # seems to be necessary to skip TWO echoes
            E.expect("\r\n") # one from the pty and one from GAP, I guess
            normal_outputs = []
            error_outputs = []
            current_outputs = normal_outputs
            while True:
                x = E.expect(['@p\d+\.','@@','@[A-Z]','@[123456!"#$%&].*\+',
                              '@e','@c','@f','@h','@i','@m','@n','@r','@s\d','@w.*\+',
                              '@x','@z'])
                current_outputs.append(E.before)
                if x == 0:   # @p
                    if E.after != '@p1.':
                        print "Warning: possibly wrong version of GAP package interface\n"
                        print "Crossing fingers and continuing\n"
                elif x == 1: #@@
                    current_outputs.append('@')
                elif x == 2: #special char
                    current_outputs.append(chr(ord(E.after[1:2])-ord('A')+1))
                elif x == 3: # garbage collection info, ignore
                    pass
                elif x == 4: # @e -- break loop
                    E.sendline(" ")
                    #E.expect("\r\n")
                    #E.expect("\r\n")
                elif x == 5: # @c completion, doesn't seem to happen when -p is in use
                    print "I didn't think GAP could do this\n"
                elif x == 6: # @f GAP error message
                    current_outputs = error_outputs;
                elif x == 7: # @h help text, but this stopped happening with new help
                    print "I didn't think GAP could do this"
                elif x == 8: # @i awaiting normal input
                    break;
                elif x == 9: # @m finished running a child
                             # not generated in GAP 4
                    print "Warning: this should never happen"
                elif x==10: #@n normal output line
                    current_outputs = normal_outputs;
                elif x==11: #@r echoing input
                    E.expect('@J')
                elif x==12: #@sN shouldn't happen
                    print "Warning: this should never happen"
                elif x==13: #@w GAP is trying to send a Window command
                    print "Warning: this should never happen"
                elif x ==14: #@x seems to be safely ignorable
                    pass
                elif x == 15:#@z GAP starting a subprocess
                             # actually not used
                        print "Warning: this should never happen"
        except pexpect.EOF:
            if not expect_eof:
                raise RuntimeError, "Unexpected EOF from %s executing %s"%(self,line)
        except IOError:
            raise RuntimeError, "IO Error from %s executing %s"%(self,line)
        return ("".join(normal_outputs),"".join(error_outputs))

    def _keyboard_interrupt(self):
        print "Interrupting %s..."%self
        os.killpg(self._expect.pid, 2)
        raise KeyboardInterrupt, "Ctrl-c pressed while running %s"%self

    def _eval_line(self, line, allow_use_file=True, wait_for_prompt=True):
        #if line.find('\n') != -1:
        #    raise ValueError, "line must not contain any newlines"
        try:
            if self._expect is None:
                self._start()
            E = self._expect
            #import pdb; pdb.set_trace()
            if allow_use_file and len(line) > self._eval_using_file_cutoff:
                return self._eval_line_using_file(line, tmp)
            try:
                (normal, error) = self._execute_line(line, wait_for_prompt=wait_for_prompt,
                                                 expect_eof= (self._quit_string() in line))

                if len(error)> 0:
                    if 'completion files' in error:
                        error += "\nPossibly try the command gap_reset_workspace() from SAGE."
                    gap_reset_workspace()
                    raise RuntimeError, "%s produced error output\n%s\n   executing %s"%(self, error,line)
                if len(normal) == 0:
                    return ''

                if isinstance(wait_for_prompt, str):
                    n = len(wait_for_prompt)
                else:
                    n = len(self._prompt)
                out = normal[:-n]
                if len(out) > 0 and out[-1] == "\n":
                    out = out[:-1]
                return out

            except (RuntimeError,),message:
                if 'EOF' in message:
                    print "** %s crashed or quit executing '%s' **"%(self, line)
                    print "Restarting %s and trying again"%self
                    self._start()
                    if line != '':
                        return self._eval_line(line, allow_use_file=allow_use_file)
                    else:
                        return ''
                else:
                    raise RuntimeError, message

        except KeyboardInterrupt:
            self._keyboard_interrupt()
            raise KeyboardInterrupt, "Ctrl-c pressed while running %s"%self
#        i = out.find("\n")
#        j = out.rfind("\r")
#        return out[i+1:j].replace('\r\n','\n')

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

    def trait_names(self):
        try:
            return self.__trait_names
        except AttributeError:
            self.__trait_names = eval(self.eval('NamesSystemGVars()')) + \
                                 eval(self.eval('NamesUserGVars()'))
        return self.__trait_names



############

def gap_reset_workspace(max_workspace_size=None):
    r"""
    Call this to completely reset the GAP workspace, which
    is used by default when SAGE first starts GAP.

    The first time you start GAP from SAGE, it saves the
    startup state of GAP in the file
    \begin{verbatim}
        $HOME/.sage/gap-workspace
    \end{verbatim}
    This is useful, since then subsequent startup of GAP
    is at least 10 times as fast.  Unfortunately, if you
    install any new code for GAP, it won't be noticed unless
    you explicitly load it, e.g., with
           gap.load_package("my_package")

    The packages sonata, guava, factint, gapdoc, grape, design, toric,
    and laguna are loaded in all cases before the workspace is saved,
    if they are available.
    """
    if os.path.exists(WORKSPACE):
        os.unlink(WORKSPACE)
    g = Gap(use_workspace_cache=False, max_workspace_size=None)
    for pkg in ['sonata', 'guava', 'factint', \
                'gapdoc', 'grape', 'design', \
                'toric', 'laguna']:
        try:
            g.load_package(pkg, verbose=True)
        except RuntimeError:
            pass
    # end for
    g.eval('SaveWorkspace("%s");'%WORKSPACE)


class GapElement(ExpectElement):
    def __getattr__(self, attrname):
        if attrname[:1] == "_":
            raise AttributeError
        return GapFunctionElement(self, attrname)

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
            raise RuntimeError, "An error occured creating an object in %s from:\n'%s'\n%s"%(self.parent().name(), self._createu, s)
        return s

    def __len__(self):
        return int(self.Length())

    def _latex_(self):
        self._check_valid()
        P = self.parent()
        try:
            s = P.eval('LaTeXObj(%s)'%self.name())
            s = s.replace('\\\\','\\').replace('"','')
            s = s.replace('%\\n',' ')
            return s
        except RuntimeError:
            return str(self)

    def _matrix_(self, R):
        r"""
        Return matrix over the (\sage) ring R determined by self, where self
        should be a Gap matrix.

            sage: s = gap("(Z(7)^0)*[[1,2,3],[4,5,6]]"); s
            [ [ Z(7)^0, Z(7)^2, Z(7) ], [ Z(7)^4, Z(7)^5, Z(7)^3 ] ]
            sage: matrix(s, GF(7))
            [1 2 3]
            [4 5 6]

            sage: s = gap("[[1,2], [3/4, 5/6]]"); s
            [ [ 1, 2 ], [ 3/4, 5/6 ] ]
            sage: m = matrix(s, QQ); m
            [  1   2]
            [3/4 5/6]
            sage: parent(m)
            Full MatrixSpace of 2 by 2 dense matrices over Rational Field

            sage: s = gap('[[Z(16),Z(16)^2],[Z(16)^3,Z(16)]]')
            sage: matrix(s, GF(16))
            [  a a^2]
            [a^3   a]
        """
        P = self.parent()
        v = self.DimensionsMat()
        n = int(v[1])
        m = int(v[2])

        from sage.matrix.matrix_space import MatrixSpace
        M = MatrixSpace(R, n, m)
        entries = [[R(self[r,c]) for c in range(1,m+1)] for r in range(1,n+1)]
        return M(entries)

    def trait_names(self):
        if '__trait_names' in self.__dict__:
            return self.__trait_names
        P = self.parent()
        v = P.eval('$SAGE.OperationsAdmittingFirstArgument(%s)'%self.name())
        v = v.replace('Tester(','').replace('Setter(','').replace('<Operation ','').replace('>','').replace(')','')
        v = eval(v)
        v = list(set(v))
        v.sort()
        self.__trait_names = v
        return v


class GapFunctionElement(FunctionElement):
    def _sage_doc_(self):
        M = self._obj.parent()
        return M.help(self._name, pager=False)


class GapFunction(ExpectFunction):
    def _sage_doc_(self):
        M = self._parent
        return M.help(self._name, pager=False)


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


