r"""
Interface to Maxima

Maxima is a free GPL'd general purpose computer algebra system whose
development started in 1968 at MIT.  It contains symbolic manipulation
algorithms, as well as implementations of special functions, including
elliptic functions and generalized hypergeometric functions. Moreover,
Maxima has implementations of many functions relating to the invariant
theory of the symmetric group $S_n$.  (However, the commands for group
invariants, and the corresponding Maxima documenation, are in French.)
For many links to Maxima documentation see
         \url{http://maxima.sourceforge.net/docs.shtml/}.

AUTHORS OF THIS MODULE:
    - William Stein (2005-12): Initial version
    - David Joyner: Improved documentation
    - William Stein (2006-01-08): Fixed bug in parsing
    - William Stein (2006-02-22): comparisons (following suggestion of David Joyner)
    - William Stein (2006-02-24): *greatly* improved robustness by adding
                                  sequence numbers to IO bracketing in _eval_line

If the string "error" (case insensitive) occurs in the output of
anything from maxima, a RuntimeError exception is raised.

EXAMPLES:
We evaluate a very simple expression in maxima.
    sage: maxima('3 * 5')
    15

We factor $x^5 - y^5$ in Maxima in several different ways.
The first way yields a Maxima object.
    sage: F = maxima.factor('x^5 - y^5')
    sage: F
    -(y - x)*(y^4 + x*y^3 + x^2*y^2 + x^3*y + x^4)
    sage: type(F)
    <class 'sage.interfaces.maxima.MaximaElement'>

Note that Maxima objects can also be displayed using ``ASCII art'';
to see a normal linear representation of any Maxima object x,
use \code{str(x)}.
    sage: F.display2d()
                               4      3    2  2    3      4
                   - (y - x) (y  + x y  + x  y  + x  y + x )

We can make this the default:
    sage: maxima.display2d(True)
    sage: F
                               4      3    2  2    3      4
                   - (y - x) (y  + x y  + x  y  + x  y + x )

You can always use \code{x.str()} to obtain the linear representation
of an object, even without changing the display2d flag.  This can
be useful for moving maxima data to other systems.
    sage: F.str()
    '-(y - x)*(y^4 + x*y^3 + x^2*y^2 + x^3*y + x^4)'

    sage: maxima.display2d(False)
    sage: F
    -(y - x)*(y^4 + x*y^3 + x^2*y^2 + x^3*y + x^4)


The \code{maxima.eval} command evaluates an expression in maxima
and returns the result as a string.

    sage: print maxima.eval('factor(x^5 - y^5)')
    -(y - x)*(y^4 + x*y^3 + x^2*y^2 + x^3*y + x^4)

We can create the polynomial $f$ as a Maxima polynomial, then call
the factor method on it.  Notice that the notation \code{f.factor()}
is consistent with how the rest of \sage works.
    sage: f = maxima('x^5 - y^5')
    sage: f^2
    (x^5 - y^5)^2
    sage: f.factor()
    -(y - x)*(y^4 + x*y^3 + x^2*y^2 + x^3*y + x^4)

Control-C interruption works well with the maxima interface,
because of the excellent implementation of maxima.  For example,
try the following sum but with a much bigger range, and hit
control-C.
    sage: maxima('sum(1/x^2, x, 1, 10)')
    1968329/1270080

\subsection{Tutorial}
We follow the tutorial at
\url{http://maxima.sourceforge.net/docs/intromax/}.

    sage: maxima('1/100 + 1/101')
    201/10100

    sage: a = maxima('(1 + sqrt(2))^5'); a
    (sqrt(2) + 1)^5
    sage: a.expand()
    29*sqrt(2) + 41

    sage: a = maxima('(1 + sqrt(2))^5')
    sage: float(a)
    82.012193308819747
    sage: a.numer()
    82.01219330881975

    sage: maxima.eval('fpprec : 100')
    '100'
    sage: a.bfloat()
    8.20121933088197564152489730020812442785204843859314941221237124017312418754011041266612384955016056b1

    sage: maxima('100!')
    93326215443944152681699238856266700490715968264381621468592963895217599993229915608941463976156518286253697920827223758251185210916864000000000000000000000000

    sage: f = maxima('(x + 3*y + x^2*y)^3')
    sage: f.expand()
    x^6*y^3 + 9*x^4*y^3 + 27*x^2*y^3 + 27*y^3 + 3*x^5*y^2 + 18*x^3*y^2 + 27*x*y^2 + 3*x^4*y + 9*x^2*y + x^3
    sage: f.subst('x=5/z')
    (5/z + 25*y/z^2 + 3*y)^3
    sage: g = f.subst('x=5/z')
    sage: h = g.ratsimp(); h
    (27*y^3*z^6 + 135*y^2*z^5 + (675*y^3 + 225*y)*z^4 + (2250*y^2 + 125)*z^3 + (5625*y^3 + 1875*y)*z^2 + 9375*y^2*z + 15625*y^3)/z^6
    sage: h.factor()
    (3*y*z^2 + 5*z + 25*y)^3/z^6

    sage: eqn = maxima(['a+b*c=1', 'b-a*c=0', 'a+b=5'])
    sage: s = eqn.solve('[a,b,c]'); s
    [[a = (25*sqrt(79)*%i + 25)/(6*sqrt(79)*%i - 34),b = (5*sqrt(79)*%i + 5)/(sqrt(79)*%i + 11),c = (sqrt(79)*%i + 1)/10],[a = (25*sqrt(79)*%i - 25)/(6*sqrt(79)*%i + 34),b = (5*sqrt(79)*%i - 5)/(sqrt(79)*%i - 11),c =  - (sqrt(79)*%i - 1)/10]]

Here is an example of solving an algebraic equation:
    sage: maxima('x^2+y^2=1').solve('y')
    [y =  - sqrt(1 - x^2),y = sqrt(1 - x^2)]
    sage: maxima('x^2 + y^2 = (x^2 - y^2)/sqrt(x^2 + y^2)').solve('y')
    [y =  - sqrt(( - y^2 - x^2)*sqrt(y^2 + x^2) + x^2),y = sqrt(( - y^2 - x^2)*sqrt(y^2 + x^2) + x^2)]

You can even nicely typeset the solution in latex:
    sage: latex(s)
    \left[ \left[ a={{25\,\sqrt{79}\,i+25}\over{6\,\sqrt{79}\,i-34}} ,   b={{5\,\sqrt{79}\,i+5}\over{\sqrt{79}\,i+11}} , c={{\sqrt{79}\,i+1  }\over{10}} \right]  , \left[ a={{25\,\sqrt{79}\,i-25}\over{6\,  \sqrt{79}\,i+34}} , b={{5\,\sqrt{79}\,i-5}\over{\sqrt{79}\,i-11}} ,   c=-{{\sqrt{79}\,i-1}\over{10}} \right]  \right]

To have the above appear onscreen via \code{xdvi}, type \code{view(s)}.
(TODO: For OS X should create pdf output and use preview instead?)

    sage: e = maxima('sin(u + v) * cos(u)^3'); e
    cos(u)^3*sin(v + u)
    sage: f = e.trigexpand(); f
    cos(u)^3*(cos(u)*sin(v) + sin(u)*cos(v))
    sage: f.trigreduce()
    (sin(v + 4*u) + sin(v - 2*u))/8 + (3*sin(v + 2*u) + 3*sin(v))/8
    sage: w = maxima('3 + k*%i')
    sage: f = w^2 + maxima('%e')^w
    sage: f.realpart()
    %e^3*cos(k) - k^2 + 9

    sage: f = maxima('x^3 * %e^(k*x) * sin(w*x)'); f
    x^3*%e^(k*x)*sin(w*x)
    sage: f.diff('x')
    k*x^3*%e^(k*x)*sin(w*x) + 3*x^2*%e^(k*x)*sin(w*x) + w*x^3*%e^(k*x)*cos(w*x)
    sage: f.integrate('x')
    (((k*w^6 + 3*k^3*w^4 + 3*k^5*w^2 + k^7)*x^3 + (3*w^6 + 3*k^2*w^4 - 3*k^4*w^2 - 3*k^6)*x^2 + ( - 18*k*w^4 - 12*k^3*w^2 + 6*k^5)*x - 6*w^4 + 36*k^2*w^2 - 6*k^4)*%e^(k*x)*sin(w*x) + (( - w^7 - 3*k^2*w^5 - 3*k^4*w^3 - k^6*w)*x^3 + (6*k*w^5 + 12*k^3*w^3 + 6*k^5*w)*x^2 + (6*w^5 - 12*k^2*w^3 - 18*k^4*w)*x - 24*k*w^3 + 24*k^3*w)*%e^(k*x)*cos(w*x))/(w^8 + 4*k^2*w^6 + 6*k^4*w^4 + 4*k^6*w^2 + k^8)

    sage: f = maxima('1/x^2')
    sage: f.integrate('x', 1, 'inf')
    1
    sage: g = maxima('f/sinh(k*x)^4')
    sage: g.taylor('x', 0, 3)
    f/(k^4*x^4) - 2*f/(3*k^2*x^2) + 11*f/45 - 62*k^2*f*x^2/945

    sage: maxima.taylor('asin(x)','x',0, 10)
    x + x^3/6 + 3*x^5/40 + 5*x^7/112 + 35*x^9/1152

\subsection{Examples involving matrices}
We illustrate computing with the matrix whose $i,j$ entry
is $i/j$, for $i,j=1,\ldots,4$.

    sage: f = maxima.eval('f[i,j] := i/j')
    sage: A = maxima('genmatrix(f,4,4)'); A
    matrix([1,1/2,1/3,1/4],[2,1,2/3,1/2],[3,3/2,1,3/4],[4,2,4/3,1])
    sage: A.determinant()
    0
    sage: A.echelon()
    matrix([1,1/2,1/3,1/4],[0,0,0,0],[0,0,0,0],[0,0,0,0])
    sage: A.eigenvalues()
    [[0,4],[3,1]]
    sage: A.eigenvectors()
    [[[0,4],[3,1]],[1,0,0, - 4],[0,1,0, - 2],[0,0,1, - 4/3],[1,2,3,4]]

We can also compute the echelon form in \sage:
    sage: B = matrix(QQ, A)
    sage: B.echelon_form()
    [  1 1/2 1/3 1/4]
    [  0   0   0   0]
    [  0   0   0   0]
    [  0   0   0   0]
    sage: B.charpoly('x').factor()
    (x - 4) * x^3

\subsection{Laplace Transforms}
We illustrate Laplace transforms:
    sage: _ = maxima.eval("f(t) := t*sin(t)")
    sage: maxima("laplace(f(t),t,s)")
    2*s/(s^2 + 1)^2

    sage: maxima("laplace(delta(t-3),t,s)") #Dirac delta function
    %e^-(3*s)

    sage: _ = maxima.eval("f(t) := exp(t)*sin(t)")
    sage: maxima("laplace(f(t),t,s)")
    1/(s^2 - 2*s + 2)

    sage: _ = maxima.eval("f(t) := t^5*exp(t)*sin(t)")
    sage: maxima("laplace(f(t),t,s)")
    360*(2*s - 2)/(s^2 - 2*s + 2)^4 - 480*(2*s - 2)^3/(s^2 - 2*s + 2)^5 + 120*(2*s - 2)^5/(s^2 - 2*s + 2)^6
    sage: maxima("laplace(f(t),t,s)").display2d()
                                             3                 5
               360 (2 s - 2)    480 (2 s - 2)     120 (2 s - 2)
              --------------- - --------------- + ---------------
                2           4     2           5     2           6
              (s  - 2 s + 2)    (s  - 2 s + 2)    (s  - 2 s + 2)

    sage: maxima("laplace(diff(x(t),t),t,s)")
    s*?%laplace(x(t),t,s) - x(0)

    sage: maxima("laplace(diff(x(t),t,2),t,s)")
    -?%at('diff(x(t),t,1),t = 0) + s^2*?%laplace(x(t),t,s) - x(0)*s

It is difficult to read some of these without the 2d representation:
    sage.: maxima("laplace(diff(x(t),t,2),t,s)").display2d()
                         !
                d        !         2
              - -- (x(t))!      + s  laplace(x(t), t, s) - x(0) s
                dt       !
                         !t = 0

Even better, use \code{view(maxima("laplace(diff(x(t),t,2),t,s)"))} to see
a typeset version.

\subsection{Continued Fractions}

A continued fraction $a + 1/(b + 1/(c + \cdots))$ is
represented in maxima by the list $[a, b, c, \ldots]$.

    sage: maxima("cf((1 + sqrt(5))/2)")
    [1,1,1,1,2]
    sage: maxima("cf ((1 + sqrt(341))/2)")
    [9,1,2,1,2,1,17,1,2,1,2,1,17,1,2,1,2,1,17,2]

\subsection{Special examples}

In this section we illustrate calculations that would be awkward to do
(as far as I know) in non-symbolic computer algebra systems like MAGMA
or GAP.

We compute the gcd of $2x^{n+4} - x^{n+2}$ and $4x^{n+1} + 3x^n$
for arbitrary $n$.

    sage: f = maxima('2*x^(n+4) - x^(n+2)')
    sage: g = maxima('4*x^(n+1) + 3*x^n')
    sage: f.gcd(g)
    x^n

You can plot 3d graphs (via gnuplot):

    sage.: maxima('plot3d(x^2-y^2, [x,-2,2], [y,-2,2], [grid,12,12])')
    [displays a 3 dimensional graph]

You can formally evaluate sums (note the \code{nusum} command):

    sage: S = maxima('nusum(exp(1+2*i/n),i,1,n)')
    sage.: S.display2d()
                            2/n + 3                   2/n + 1
                          %e                        %e
                   ----------------------- - -----------------------
                      1/n         1/n           1/n         1/n
                   (%e    - 1) (%e    + 1)   (%e    - 1) (%e    + 1)

We formally compute the limit as $n\to\infty$ of $2S/n$ as follows:

    sage: T = S*maxima('2/n')
    sage: T.tlimit('n','inf')
    %e^3 - %e

\subsection{Miscellaneous}
Obtaining digits of $\pi$:
    sage: maxima.eval('fpprec : 100')
    '100'
    sage: maxima(pi).bfloat()
    3.141592653589793238462643383279502884197169399375105820974944592307816406286208998628034825342117068b0

Defining functions in maxima:
    sage: maxima.eval('fun[a] := a^2')
    'fun[a] := a^2'
    sage: maxima('fun[10]')
    100

\subsection{Interactivity}
Unfortunately maxima doesn't seem to have a non-interactive mode,
which is needed for the \sage interface.  If any \sage call leads
to maxima interactively answering questions, then the questions
can't be answered and the maxima session may hang.
See the discussion at \url{http://www.ma.utexas.edu/pipermail/maxima/2005/011061.html} for some ideas about how to fix this problem.  An
example that illustrates this problem is
\code{maxima.eval('integrate (exp(a*x), x, 0, inf)')}.

\subsection{Latex Output}
The latex output of Maxima is not perfect.  E.g.,

    sage: maxima.eval('tex(sin(u) + sinh(v^2))')
    '$$\\sinhv^2 + \\sinu$$false'

Notice the lack of space after the sin macro, which is a latex syntax
error.  In \sage this is automatically fixed via a substition for
trig functions, which may have potentially bad side effects:

    sage: latex(maxima('sin(u) + sinh(v^2)'))
    \sinh v^2+\sin u

It would be nice if somebody would fix this problem.  One way would
be to improve Maxima by making the fix to Maxima and giving this back
to the Maxima people.

Here's another example:

    sage: g = maxima('exp(3*%i*x)/(6*%i) + exp(%i*x)/(2*%i) + c')
    sage: latex(g)
    -{{i\,e^{3\,i\,x}}\over{6}}-{{i\,e^{i\,x}}\over{2}}+c

\subsection{Long Input}
The MAXIMA interface reads in even very long input (using files) in a
robust manner, as long as you are creating a new object.
\note{Using \code{maxima.eval} for long input
is much less robust, and is not recommended.}

    sage: t = '"%s"'%10^10000   # ten thousand character string.
    sage: a = maxima(t)
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

import os, re
import pexpect
cygwin = os.uname()[0][:6]=="CYGWIN"

from expect import Expect, ExpectElement, FunctionElement, ExpectFunction, tmp
from pexpect import EOF

import sage.rings.all
#import sage.rings.complex_number2 as complex_number

from sage.misc.misc import verbose, DOT_SAGE, SAGE_ROOT

from sage.misc.multireplace import multiple_replace

SAGE_START = '_s_start_'
SAGE_END = '_s_stop_'
cnt = 0
seq = 0

COMMANDS_CACHE = '%s/maxima_commandlist_cache.sobj'%DOT_SAGE

import sage.server.support

# The Maxima "apropos" command, e.g., apropos(det) gives a list
# of all identifiers that begin in a certain way.  This could
# maybe be useful somehow... (?)  Also maxima has a lot for getting
# documentation from the system -- this could also be useful.

class Maxima(Expect):
    """
    Interface to the Maxima interpreter.
    """
    def __call__(self, x):
        import sage.rings.all
        if sage.rings.all.is_Infinity(x):
            return Expect.__call__(self, 'inf')
        else:
            return Expect.__call__(self, x)

    def __init__(self, script_subdirectory=None, logfile=None, server=None):
        """
        Create an instance of the Maxima interpreter.
        """
        # TODO: Input and output prompts in maxima can be changed by
        # setting inchar and outchar..
        eval_using_file_cutoff = 200
        self.__eval_using_file_cutoff = eval_using_file_cutoff
        Expect.__init__(self,
                        name = 'maxima',
                        prompt = '\(\%i[0-9]+\)',
                        command = "maxima --disable-readline",
                        maxread = 1,    # CRUCIAL to use less buffering for maxima (or get all kinds of hangs on OS X and 64-bit machines, etc!
                        script_subdirectory = script_subdirectory,
                        restart_on_ctrlc = False,
                        verbose_start = False,
                        init_code = ['display2d : false',  # no ascii art output
                                     #'load("mactex-utilities")'   # latex instead of plain tex from tex command (broken in maxima-5.11.0!)
                                     ],
                        logfile = logfile,
                        eval_using_file_cutoff=eval_using_file_cutoff)
        self._display2d = False

    def __getattr__(self, attrname):
        if attrname[:1] == "_":
            raise AttributeError
        return MaximaExpectFunction(self, attrname)

    def _start(self):
        # For some reason sending a single input line at startup avoids
        # lots of weird timing issues when doing doctests.
        Expect._start(self)
        self(1)

    # this doesn't work.
    #def x_start(self):
    #    Expect._start(self)
    #    self._expect.sendline('inchar:"__SAGE__";')
    #    self._change_prompt('__SAGE__[0-9]+\)')
    #    self.expect().expect('__SAGE__[0-9]+\)')

    def _eval_line_using_file(self, line, tmp):
        F = open(tmp, 'w')
        F.write(line)
        F.close()
        if self._expect is None:
            self._start()
        # For some reason this trivial comp
        # keeps certain random freezes from occuring.  Do not remove this.
        # The space before the \n is also important.
        self._expect.sendline('0;batchload("%s"); \n'%tmp)
        self._expect.expect(self._prompt)
        return ''

    def __reduce__(self):
        return reduce_load_Maxima, tuple([])

    def _quit_string(self):
        return 'quit();'

    def _eval_line(self, line, reformat=True, allow_use_file=False,
                   wait_for_prompt=True):
        if not wait_for_prompt:
            return Expect._eval_line(self, line)
        line = line.rstrip().rstrip(';')
        if line == '':
            return ''
        global seq
        seq += 1
        start = SAGE_START + str(seq)
        end = SAGE_END + str(seq)
        line = '%s;\n%s; %s;'%(start, line, end)
        if self._expect is None:
            self._start()
        if allow_use_file and self.__eval_using_file_cutoff and \
                            len(line) > self.__eval_using_file_cutoff:
            return self._eval_line_using_file(line, tmp)
        try:
            E = self._expect
            #print "in = '%s'"%line
            E.sendline(line)
            self._expect.expect(end)
            # We have timeouts below, since getting the end above
            # means the computation completed, but on some systems
            # (Cygwin) the expect interface can sometimes hang getting
            # the final prompts.
            try:
                self._expect.expect(end, timeout=1)
                out = self._expect.before
                self._expect.expect(self._prompt, timeout=1)
                out += self._expect.before
            except pexpect.TIMEOUT:
                out = self._expect.before
            #print "out = '%s'"%out

        except EOF:
            if self._quit_string() in line:
                return ''
        except KeyboardInterrupt:
            self._keyboard_interrupt()
            return ''

        if 'Incorrect syntax:' in out:
            raise RuntimeError, out

        import os
        if cygwin:
            # for reasons I can't deduce yet, maxima behaves somewhat
            # differently under cygwin...
            out = out.lstrip(';')
        else:
	    i = out.rfind(start)
	    j = out.rfind(end)
            out = out[i+len(start):j]

        if not reformat:
            return out
        if 'error' in out:
            return out
        out = out.lstrip()
        i = out.find('(%o')
        out0 = out[:i].strip()
        i += out[i:].find(')')
        out1 = out[i+1:].strip()
        out = out0 + out1
        out = ''.join(out.split())    # no whitespace
        i = out.rfind(';;')
        if i != -1:
            out = out[i+2:]
        out = out.replace('-', ' - ').replace('+',' + ').replace('=',' = ').replace(': =',' :=').replace('^ - ','^-')
        if out[:3] == ' - ':
            out = '-' + out[3:]
        out = out.replace('E - ', 'E-')
        out = out.replace('%e - ', '%e-')
        i = out.rfind('(%o')
        return out[:i]


    ###########################################
    # Interactive help
    ###########################################

    # This doesn't work because of how weird interaction is.
    # System call method below is much more robust.
    #def help(self, s):
    #    return "Help on Maxima commands currently not implemented (see %s/devel/sage/interfaces/maxima.py if you want to try to implement it)."%SAGE_ROOT
##         if self._expect is None:
##             self._start()
##         E = self._expect
##         if E is None:
##             raise RuntimError, "unable to start maxima"
##         E.sendline('describe("%s");'%s)
##         #old_timeout = E.timeout
##         #E.timeout = 0.2
##         E.expect("`all' or `none': ")
##         E.sendline('all\n\n')
##         E.sendline('all\n\n')
##         #E.expect(self._prompt)
##         E.expect('%o')
##         print E.before
##     # override the builtin describe command

    def help(self, s):
        if sage.server.support.EMBEDDED_MODE:
            os.system('maxima -r "describe(%s); "< /dev/null'%s)
        else:
            os.system('maxima -r "describe(%s);"'%s)

    def example(self, s):
        if sage.server.support.EMBEDDED_MODE:
            os.system('maxima -r "example(%s);" < /dev/null'%s)
        else:
            os.system('maxima -r "example(%s);"'%s)

    describe = help

    def demo(self, s):
        if sage.server.support.EMBEDDED_MODE:
            os.system('maxima -r "demo(%s);" < /dev/null'%s)
        else:
            os.system('maxima -r "demo(%s);"'%s)

    def completions(self, s):
        """
        Return all commands that complete the command starting with the
        string s.   This is like typing s[tab] in the maple interpreter.
        """
        s = self.eval('apropos(%s)'%s).replace('\\ - ','-')
        return [x for x in s[1:-1].split(',') if x[0] != '?']

    def _commands(self):
        """
        Return list of all commands defined in Maxima.
        """
        try:
            return self.__commands
        except AttributeError:
            self.__commands = sum([self.completions(chr(97+n)) for n in range(26)], [])
        return self.__commands

    def trait_names(self, verbose=True, use_disk_cache=True):
        try:
            return self.__trait_names
        except AttributeError:
            import sage.misc.persist
            if use_disk_cache:
                try:
                    self.__trait_names = sage.misc.persist.load(COMMANDS_CACHE)
                    return self.__trait_names
                except IOError:
                    pass
            if verbose:
                print "\nBuilding Maxima command completion list (this takes"
                print "a few seconds only the first time you do it)."
                print "To force rebuild later, delete %s."%COMMANDS_CACHE
            v = self._commands()
            self.__trait_names = v
            sage.misc.persist.save(v, COMMANDS_CACHE)
            return v

    def _object_class(self):
        return MaximaElement

    def _true_symbol(self):
        return 'true'

    def _false_symbol(self):
        return 'false'

    def function(self, args, defn, repr=None, latex=None):
        """
        Return the Maxima function with given arguments
        and definition.

        INPUT:
            args -- a string with variable names separated by commas
            defn -- a string (or Maxima expression) that defines
                    a function of the arguments in Maxima.
            repr -- an optional string; if given, this is how the function will print.

        EXAMPLES:
            sage: f = maxima.function('x', 'sin(x)')
            sage: f(3.2)
            -.05837414342758009
            sage: f = maxima.function('x,y', 'sin(x)+cos(y)')
            sage: f(2,3.5)
            sin(2) - .9364566872907963
            sage: f
            sin(x)+cos(y)

            sage: g = f.integrate('z'); g
            (cos(y) + sin(x))*z
            sage: g(1,2,3)
            3*(cos(2) + sin(1))

        The function definition can be a maxima object:
            sage: an_expr = maxima('sin(x)*gamma(x)')
            sage: t = maxima.function('x', an_expr)
            sage: t
            gamma(x)*sin(x)
            sage: t(2)
            sin(2)
            sage: float(t(2))
            0.90929742682568171
            sage: loads(t.dumps())
            gamma(x)*sin(x)
        """
        name = self._next_var_name()
        defn = str(defn)
        args = str(args)
        maxima.eval('%s(%s) := %s'%(name, args, defn))
        if repr is None:
            repr = defn
        f = MaximaFunction(self, name, repr, args, latex)
        return f

    def set(self, var, value):
        """
        Set the variable var to the given value.
        """
        cmd = '%s : %s$"";'%(var, str(value).rstrip(';'))
        out = self._eval_line(cmd, reformat=False, allow_use_file=True)

        if out.find("error") != -1:
            raise TypeError, "Error executing code in Maxima\nCODE:\n\t%s\nMaxima ERROR:\n\t%s"%(cmd, out)


    def get(self, var):
        """
        Get the string value of the variable var.
        """
        s = self._eval_line('%s'%var)
        return s

    #def clear(self, var):
    #    """
    #    Clear the variable named var.
    #    """
    #    if self._expect is None:
    #        return
    #    self._expect.sendline('kill(%s);'%var)
    #    self._expect.expect(self._prompt)

    def console(self):
        maxima_console()

    def version(self):
        return maxima_version()

    def display2d(self, flag=True):
        """
        Set the flag that determines whether Maxima objects are
        printed using their 2-d ASCII art representation.  When the
        maxima interface starts the default is that objects are not
        represented in 2-d.

        INPUT:
            flag -- bool (default: True)

        EXAMPLES
            sage: maxima('1/2')
            1/2
            sage: maxima.display2d(True)
            sage: maxima('1/2')
                                           1
                                           -
                                           2
            sage: maxima.display2d(False)
        """
        self._display2d = bool(flag)

    def plot2d(self, *args):
        r"""
        Plot a 2d graph using Maxima / gnuplot.

        maxima.plot2d(f, '[var, min, max]', options)

        INPUT:
            f -- a string representing a function (such as f="sin(x)")
            [var, xmin, xmax]
            options -- an optional string representing plot2d options in gnuplot format

        EXAMPLES:
            sage.: maxima.plot2d('sin(x)','[x,-5,5]')
            sage.: opts = '[gnuplot_term, ps], [gnuplot_out_file, "sin-plot.eps"]'
            sage.: maxima.plot2d('sin(x)','[x,-5,5]',opts)

        The eps file is saved in the current directory.
        """
        self('plot2d(%s)'%(','.join([str(x) for x in args])))

    def plot2d_parametric(self, r, var, trange, nticks=50, options=None):
        r"""
        Plots r = [x(t), y(t)] for t = tmin...tmax using gnuplot with options

        INPUT:
            r -- a string representing a function (such as r="[x(t),y(t)]")
            var -- a string representing the variable (such as var = "t")
            trange -- [tmin, tmax] are numbers with tmin<tmax
            nticks -- int (default: 50)
            options -- an optional string representing plot2d options in gnuplot format

        EXAMPLES:
            sage.: maxima.plot2d_parametric(["sin(t)","cos(t)"], "t",[-3.1,3.1])

            sage.: opts = '[gnuplot_preamble, "set nokey"], [gnuplot_term, ps], [gnuplot_out_file, "circle-plot.eps"]'
            sage.: maxima.plot2d_parametric(["sin(t)","cos(t)"], "t", [-3.1,3.1], options=opts)

        The eps file is saved to the current working directory.

        Here is another fun plot:
            sage.: maxima.plot2d_parametric(["sin(5*t)","cos(11*t)"], "t", [0,2*pi()], nticks=400)
        """
        tmin = trange[0]
        tmax = trange[1]
        cmd = "plot2d([parametric, %s, %s, [%s, %s, %s], [nticks, %s]]"%( \
                   r[0], r[1], var, tmin, tmax, nticks)
        if options is None:
            cmd += ")"
        else:
            cmd += ", %s)"%options
        self(cmd)

    def plot3d(self, *args):
        r"""
        Plot a 3d graph using Maxima / gnuplot.

        maxima.plot3d(f, '[x, xmin, xmax]', '[y, ymin, ymax]', '[grid, nx, ny]', options)

        INPUT:
            f -- a string representing a function (such as f="sin(x)")
            [var, min, max]

        EXAMPLES:
            sage.: maxima.plot3d('1 + x^3 - y^2', '[x,-2,2]', '[y,-2,2]', '[grid,12,12]')
            sage.: maxima.plot3d('sin(x)*cos(y)', '[x,-2,2]', '[y,-2,2]', '[grid,30,30]')
            sage.: opts = '[gnuplot_term, ps], [gnuplot_out_file, "sin-plot.eps"]'
            sage.: maxima.plot3d('sin(x+y)', '[x,-5,5]', '[y,-1,1]', opts)

        The eps file is saved in the current working directory.
        """
        self('plot3d(%s)'%(','.join([str(x) for x in args])))

    def plot3d_parametric(self, r, vars, urange, vrange, options=None):
        r"""
        Plot a 3d parametric graph with r=(x,y,z), x = x(u,v), y = y(u,v), z = z(u,v),
        for u = umin...umax, v = vmin...vmax using gnuplot with options.

        INPUT:
            x, y, z -- a string representing a function (such as x="u^2+v^2", ...)
            vars is a list or two strings representing variables (such as vars = ["u","v"])
            urange -- [umin, umax]
            vrange -- [vmin, vmax] are lists of numbers with
            umin < umax, vmin < vmax
            options -- optional string representing plot2d options in gnuplot format

        OUTPUT:
            displays a plot on screen or saves to a file

        EXAMPLES:
            sage.: maxima.plot3d_parametric(["v*sin(u)","v*cos(u)","v"], ["u","v"],[-3.2,3.2],[0,3])
            sage.: opts = '[gnuplot_term, ps], [gnuplot_out_file, "sin-cos-plot.eps"]'
            sage.: maxima.plot3d_parametric(["v*sin(u)","v*cos(u)","v"], ["u","v"],[-3.2,3.2],[0,3],opts)

        The eps file is saved in the current working directory.

        Here is a torus:

            sage.: _ = maxima.eval("expr_1: cos(y)*(10.0+6*cos(x)); expr_2: sin(y)*(10.0+6*cos(x)); expr_3: -6*sin(x);")  # optional
            sage.: maxima.plot3d_parametric(["expr_1","expr_2","expr_3"], ["x","y"],[0,6],[0,6])

        Here is a Mobius strip:
            sage.: x = "cos(u)*(3 + v*cos(u/2))"
            sage.: y = "sin(u)*(3 + v*cos(u/2))"
            sage.: z = "v*sin(u/2)"
            sage.: maxima.plot3d_parametric([x,y,z],["u","v"],[-3.1,3.2],[-1/10,1/10])
        """
        umin = urange[0]
        umax = urange[1]
        vmin = vrange[0]
        vmax = vrange[1]
        cmd = 'plot3d([%s, %s, %s], [%s, %s, %s], [%s, %s, %s]'%(
            r[0], r[1], r[2], vars[0], umin, umax, vars[1], vmin, vmax)
        if options is None:
            cmd += ')'
        else:
            cmd += ', %s)'%options
        maxima(cmd)

    def de_solve(maxima, de, vars, ics=None):
        """
        Solves a 1st or 2nd order ordinary differential equation (ODE)
        in two variables, possibly with initial conditions.

        INPUT:
            de -- a string representing the ODE
            vars -- a list of strings representing the two variables.
            ics -- a triple of numbers [a,b1,b2] representing
                   y(a)=b1, y'(a)=b2

        EXAMPLES:
            sage.: maxima.de_solve('diff(y,x,2) + 3*x = y', ['x','y'], [1,1,1])
            y = 3*x - 2*%e^(x - 1)
            sage.: maxima.de_solve('diff(y,x,2) + 3*x = y', ['x','y'])
            y = %k1*%e^x + %k2*%e^ - x + 3*x
            sage.: maxima.de_solve('diff(y,x) + 3*x = y', ['x','y'])
            y = (%c - 3*( - x - 1)*%e^ - x)*%e^x
            sage.: maxima.de_solve('diff(y,x) + 3*x = y', ['x','y'],[1,1])
            y =  - %e^ - 1*(5*%e^x - 3*%e*x - 3*%e)
        """
        if not isinstance(vars, str):
            str_vars = '%s, %s'%(vars[1], vars[0])
        else:
            str_vars = vars
        maxima.eval('depends(%s)'%str_vars)
        m = maxima(de)
        a = 'ode2(%s, %s)'%(m.name(), str_vars)
        if ics != None:
            if len(ics) == 3:
                cmd = "ic2("+a+",%s=%s,%s=%s,diff(%s,%s)=%s);"%(vars[0],ics[0], vars[1],ics[1], vars[1], vars[0], ics[2])
                return maxima(cmd)
            if len(ics) == 2:
                return maxima("ic1("+a+",%s=%s,%s=%s);"%(vars[0],ics[0], vars[1],ics[1]))
        return maxima(a+";")

    def de_solve_laplace(self, de, vars, ics=None):
        """
        Solves an ordinary differential equation (ODE) using Laplace transforms.

        INPUT:
            de -- a string representing the ODE
                  (e.g., de = "diff(f(x),x,2)=diff(f(x),x)+sin(x)")
            vars -- a list of strings representing the variables
                  (e.g., vars = ["x","f"])
            ics -- a list of numbers representing initial conditions,
                   with symbols allowed which are represented by strings
                   (eg, f(0)=1, f'(0)=2 is ics = [0,1,2])

        EXAMPLES:
            sage.: maxima.clear('x'); maxima.clear('f')
            sage.: maxima.de_solve_laplace("diff(f(x),x,2) = 2*diff(f(x),x)-f(x)", ["x","f"], [0,1,2])
            f(x) = x*%e^x + %e^x

            sage.: maxima.clear('x'); maxima.clear('f')
            sage.: f = maxima.de_solve_laplace("diff(f(x),x,2) = 2*diff(f(x),x)-f(x)", ["x","f"])
            sage.: f
            f(x) = x*%e^x*(at('diff(f(x),x,1),x = 0)) - f(0)*x*%e^x + f(0)*%e^x
            sage.: f.display2d()
                                               !
                                   x  d        !                  x          x
                        f(x) = x %e  (-- (f(x))!     ) - f(0) x %e  + f(0) %e
                                      dx       !
                                               !x = 0


        \note{The second equation sets the values of $f(0)$ and
        $f'(0)$ in maxima, so subsequent ODEs involving these
        variables will have these initial conditions automatically
        imposed.}
        """
        if not (ics is None):
            d = len(ics)
            for i in range(0,d-1):
                ic = 'atvalue(diff(%s(%s), %s, %s), %s = %s, %s)'%(
                    vars[1], vars[0], vars[0], i, vars[0], ics[0], ics[1+i])
                maxima.eval(ic)
        return maxima('desolve(%s, %s(%s))'%(de, vars[1], vars[0]))

    def solve_linear(self, eqns,vars):
        """
        Wraps maxima's linsolve.

        INPUT:
        eqns is a list of m strings, each rperesenting a linear question
        in m <= n variables
        vars is a list of n strings, each representing a variable

        EXAMPLES:
            sage: eqns = ["x + z = y","2*a*x - y = 2*a^2","y - 2*z = 2"]
            sage: vars = ["x","y","z"]
            sage.: maxima.solve_linear(eqns, vars)
            [x = a + 1,y = 2*a,z = a - 1]
        """
        eqs = "["
        for i in range(len(eqns)):
            if i<len(eqns)-1:
                eqs = eqs + eqns[i]+","
            if  i==len(eqns)-1:
                eqs = eqs + eqns[i]+"]"
        vrs = "["
        for i in range(len(vars)):
            if i<len(vars)-1:
                vrs = vrs + vars[i]+","
            if  i==len(vars)-1:
                vrs = vrs + vars[i]+"]"
        return self('linsolve(%s, %s)'%(eqs, vrs))

    def unit_quadratic_integer(self, n):
        r"""
        Finds a unit of the ring of integers of the quadratic number
        field $\Q(\sqrt{n})$, $n>1$, using the qunit maxima command.

        EXAMPLE:
            sage: u = maxima.unit_quadratic_integer(101)
            sage: u.parent()
            Number Field in a with defining polynomial x^2 - 101
            sage: u
            a + 10
            sage: u = maxima.unit_quadratic_integer(13)
            sage: u
            5*a + 18
            sage: u.parent()
            Number Field in a with defining polynomial x^2 - 13
        """
        from sage.rings.all import QuadraticField, Integer
        # Take square-free part so sqrt(n) doesn't get simplified further by maxima
        # (The original version of this function would yield wrong answers if
        # n is not squarefree.)
        n = Integer(n).square_free_part()
        if n < 1:
            raise ValueError, "n (=%s) must be >= 1"%n
        s = str(self('qunit(%s)'%n)).lower()
        r = re.compile('sqrt\(.*\)')
        s = r.sub('a', s)
        a = QuadraticField(n, 'a').gen()
        return eval(s)

    def plot_list(self, ptsx, ptsy, options=None):
        r"""
        Plots a curve determined by a sequence of points.

        INPUT:
            ptsx -- [x1,...,xn], where the xi and yi are real,
            ptsy -- [y1,...,yn]
            options -- a string representing maxima plot2d options.

        The points are (x1,y1), (x2,y2), etc.

        This function requires maxima 5.9.2 or newer.

        \note{More that 150 points can sometimes lead to the program
        hanging.  Why?}

        EXAMPLES:
            sage.: zeta_ptsx = [ (pari(1/2 + i*I/10).zeta().real()).precision(1) for i in range (70,150)]
            sage.: zeta_ptsy = [ (pari(1/2 + i*I/10).zeta().imag()).precision(1) for i in range (70,150)]
            sage.: maxima.plot_list(zeta_ptsx, zeta_ptsy)
            sage.: opts='[gnuplot_preamble, "set nokey"], [gnuplot_term, ps], [gnuplot_out_file, "zeta.eps"]'
            sage.: maxima.plot_list(zeta_ptsx, zeta_ptsy, opts)
        """
        cmd = 'plot2d([discrete,%s, %s]'%(ptsx, ptsy)
        if options is None:
            cmd += ')'
        else:
            cmd += ', %s)'%options
        self(cmd)


    def plot_multilist(self, pts_list, options=None):
        r"""
        Plots a list of list of points pts_list=[pts1,pts2,...,ptsn],
        where each ptsi is of the form [[x1,y1],...,[xn,yn]]
        x's must be integers and y's reals
        options is a string representing maxima plot2d options.

        Requires maxima 5.9.2 at least.
        \note{More that 150 points can sometimes lead to the
        program hanging.}

        EXAMPLES:
            sage.: xx = [ i/10.0 for i in range (-10,10)]
            sage.: yy = [ i/10.0 for i in range (-10,10)]
            sage.: x0 = [ 0 for i in range (-10,10)]
            sage.: y0 = [ 0 for i in range (-10,10)]
            sage.: zeta_ptsx1 = [ (pari(1/2+i*I/10).zeta().real()).precision(1) for i in range (10)]
            sage.: zeta_ptsy1 = [ (pari(1/2+i*I/10).zeta().imag()).precision(1) for i in range (10)]
            sage.: maxima.plot_multilist([[zeta_ptsx1,zeta_ptsy1],[xx,y0],[x0,yy]])
            sage.: zeta_ptsx1 = [ (pari(1/2+i*I/10).zeta().real()).precision(1) for i in range (10,150)]
            sage.: zeta_ptsy1 = [ (pari(1/2+i*I/10).zeta().imag()).precision(1) for i in range (10,150)]
            sage.: maxima.plot_multilist([[zeta_ptsx1,zeta_ptsy1],[xx,y0],[x0,yy]])
            sage.: opts='[gnuplot_preamble, "set nokey"]'
            sage.: maxima.plot_multilist([[zeta_ptsx1,zeta_ptsy1],[xx,y0],[x0,yy]],opts)
        """
        n = len(pts_list)
        cmd = '['
        for i in range(n):
            if i < n-1:
                cmd = cmd+'[discrete,'+str(pts_list[i][0])+','+str(pts_list[i][1])+'],'
            if i==n-1:
                cmd = cmd+'[discrete,'+str(pts_list[i][0])+','+str(pts_list[i][1])+']]'
        #print cmd
        if options is None:
            self('plot2d('+cmd+')')
        else:
            self('plot2d('+cmd+','+options+')')


class MaximaElement(ExpectElement):
    def __call__(self, x):
        self._check_valid()
        P = self.parent()
        return P('%s(%s)'%(self.name(), x))

    def __cmp__(self, other):
        """
        EXAMPLES:
            sage: a = maxima(1); b = maxima(2)
            sage: a == b
            False
            sage: a < b
            True
            sage: a > b
            False
            sage: b < a
            False
            sage: b > a
            True

        We can also compare more complicated object such as functions:
            sage: f = maxima('sin(x)'); g = maxima('cos(x)')
            sage: -f == g.diff('x')
            True
        """

        # thanks to David Joyner for telling me about using "is".
        P = self.parent()
        if P.eval("is (%s < %s)"%(self.name(), other.name())) == P._true_symbol():
            return -1
        elif P.eval("is (%s > %s)"%(self.name(), other.name())) == P._true_symbol():
            return 1
        elif P.eval("is (%s = %s)"%(self.name(), other.name())) == P._true_symbol():
            return 0
        else:
            return -1  # everything is supposed to be comparable in Python, so we define
                       # the comparison thus when no comparable in interfaced system.
    def numer(self):
        return self.comma('numer')

    def real(self):
        return self.realpart()

    def imag(self):
        return self.imagpart()

    def _complex_mpfr_field_(self, CC):
        """
        EXAMPLES:
            sage: CC(maxima('1+%i'))
             1.00000000000000 + 1.00000000000000*I
            sage: CC(maxima('2342.23482943872+234*%i'))
             2342.23482943872 + 234.000000000000*I
            sage: ComplexField(10)(maxima('2342.23482943872+234*%i'))
             2300 + 230*I
        """
        return sage.rings.all.ComplexNumber( CC, self.real(), self.imag() )

    def str(self):
        self._check_valid()
        P = self.parent()
        return P.get(self._name)

    def __repr__(self):
        self._check_valid()
        P = self.parent()
        if P._display2d:
            return self.display2d(onscreen=False)
        else:
            return P.get(self._name)

    def display2d(self, onscreen=True):
        """
        EXAMPLES:
            sage: F = maxima('x^5 - y^5').factor()
            sage: F.display2d ()
                                   4      3    2  2    3      4
                       - (y - x) (y  + x y  + x  y  + x  y + x )
        """
        self._check_valid()
        P = self.parent()
        s = P._eval_line('display2d : true; %s'%self.name(), reformat=False)
        P._eval_line('display2d : false', reformat=False)
        if not cygwin:
            i = s.find('true')
            i += s[i:].find('\n')
            s = s[i+1:]
        i = s.find('true')
        i += s[i:].find('\n')
        j = s.rfind('(%o')
        s = s[i:j-2]
        i = s.find('(%o')
        j = i + s[i:].find(')')
        s = s[:i] + ' '*(j-i+1) + s[j+1:]
        s = s.lstrip('\n')
        if onscreen:
            print s
        else:
            return s

    def diff(self, var='x', n=1):
        """
        Return the n-th derivative of self.

        INPUT:
            var -- variable (default: 'x')
            n -- integer (default: 1)

        OUTPUT:
            n-th derivative of self with respect to the variable var

        EXAMPLES:
            sage: f = maxima('x^2')
            sage: f.diff()
            2*x
            sage: f.diff('x')
            2*x
            sage: f.diff('x', 2)
            2
            sage: maxima('sin(x^2)').diff('x',4)
            16*x^4*sin(x^2) - 12*sin(x^2) - 48*x^2*cos(x^2)

            sage: f = maxima('x^2 + 17*y^2')
            sage: f.diff('x')
            2*x
            sage: f.diff('y')
            34*y
        """
        return ExpectElement.__getattr__(self, 'diff')(var, n)

    derivative = diff

    def nintegral(self, var='x', a=0, b=1,
                  desired_relative_error='1e-8',
                  maximum_num_subintervals=200):
        r"""
        Return a numerical approximation to the integral
        of self from a to b.

        INPUT:
            var -- variable to integrate with respect to
            a -- lower endpoint of integration
            b -- upper endpoint of integration
            desired_relative_error -- (default: '1e-8') the desired
                 relative error
            maximum_num_subintervals -- (default: 200) maxima number
                 of subintervals

        OUTPUT:
            -- approximation to the integral
            -- estimated absolute error of the approximation
            -- the number of integrand evaluations
            -- an error code:
                  0 -- no problems were encountered
                  1 -- too many subintervals were done
                  2 -- excessive roundoff error
                  3 -- extremely bad integrand behavior
                  4 -- failed to converge
                  5 -- integral is probably divergent or slowly convergent
                  6 -- the input is invalid

        EXAMPLES:
            sage: maxima('exp(-sqrt(x))').nintegral('x',0,1)
            (.5284822353142306, 4.163314137883845E-11, 231, 0)

        Note that GP also does numerical integration, and can do
        so to very high precision very quickly:
            sage: gp('intnum(x=0,1,exp(-sqrt(x)))')
            0.5284822353142307136179049194             # 32-bit
            0.52848223531423071361790491935415653021   # 64-bit
            sage: _ = gp.set_precision(80)
            sage: gp('intnum(x=0,1,exp(-sqrt(x)))')
            0.52848223531423071361790491935415653021675547587292866196865279321015401702040079
        """
        from sage.rings.all import Integer
        v = self.quad_qags(var, a, b, desired_relative_error,
                           maximum_num_subintervals)
        return v[0], v[1], Integer(v[2]), Integer(v[3])

    def integral(self, var='x', min=None, max=None):
        r"""
        Return the integral of self with respect to the variable x.

        INPUT:
            var -- variable
            min -- default: None
            max -- default: None

        Returns the definite integral if xmin is not None,
        otherwise returns an indefinite integral.

        EXAMPLES:
            sage: maxima('x^2+1').integral()
            x^3/3 + x
            sage: maxima('x^2+ 1 + y^2').integral('y')
            y^3/3 + x^2*y + y
            sage: maxima('x / (x^2+1)').integral()
            log(x^2 + 1)/2
            sage: maxima('1/(x^2+1)').integral()
            atan(x)
            sage.: maxima('1/(x^2+1)').integral('x', 0, infinity)
            %pi/2
            sage: maxima('x/(x^2+1)').integral('x', -1, 1)
            0

            sage: f = maxima('exp(x^2)').integral('x',0,1); f
            -sqrt(%pi)*%i*erf(%i)/2
            sage: f.numer()         # I wonder how to get a real number (~1.463)??
            -.8862269254527579*%i*erf(%i)
        """
        I = ExpectElement.__getattr__(self, 'integrate')
        if min is None:
            return I(var)
        else:
            if max is None:
                raise ValueError, "neither or both of min/max must be specified."
            return I(var, min, max)

    integrate = integral




    def __float__(self):
        return float(str(self.numer()))

    def __len__(self):
        """
        Return the length of a list.

        EXAMPLES:
            sage: v = maxima('create_list(x^i,i,0,5)')
            sage: len(v)
            6
        """
        self._check_valid()
        return int(self.parent().eval('length(%s)'%self.name()))

    def __getattr__(self, attrname):
        if attrname[:1] == "_":
            raise AttributeError
        return MaximaFunctionElement(self, attrname)

    def __getitem__(self, n):
        r"""
        Return the n-th element of this list.

        \note{Lists are 0-based when accessed via the \sage interface,
        not 1-based as they are in the Maxima interpreter.}

        EXAMPLES:
            sage: v = maxima('create_list(i*x^i,i,0,5)'); v
            [0,x,2*x^2,3*x^3,4*x^4,5*x^5]
            sage: v[3]
            3*x^3
            sage: v[0]
            0
            sage: v[10]
            Traceback (most recent call last):
            ...
            IndexError: n = (10) must be between 0 and 5
        """
        n = int(n)
        if n < 0 or n >= len(self):
            raise IndexError, "n = (%s) must be between %s and %s"%(n, 0, len(self)-1)
        return ExpectElement.__getitem__(self, n+1)

    def subst(self, val):
        return self.comma(val)

    def comma(self, args):
        self._check_valid()
        P = self.parent()
        return P('%s, %s'%(self.name(), args))

    def _latex_(self):
        self._check_valid()
        P = self.parent()
        s = maxima._eval_line('tex(%s)'%self.name(), reformat=False)
        if not '$$' in s:
            raise RuntimeError, "Error texing maxima object."
        i = s.find('$$')
        j = s.rfind('$$')
        s = s[i+2:j]
        s = multiple_replace({'\r\n':' ',
                              '\\%':'',
                              '\\arcsin ':'\\sin^{-1} ',
                              '\\arccos ':'\\cos^{-1} ',
                              '\\arctan ':'\\tan^{-1} '}, s)

        # Fix a maxima bug, which gives a latex representation of multiplying
        # two numbers as a single space. This was really bad when 2*17^(1/3)
        # gets TeXed as '2 17^{\frac{1}{3}}'
        #
        # This regex matches a string of spaces preceeded by either a '}', a
        # decimal digit, or a ')', and followed by a decimal digit. The spaces
        # get replaced by a '\cdot'.
        s = re.sub(r'(?<=[})\d]) +(?=\d)', '\cdot', s)
        return s

    def trait_names(self):
        return self.parent().trait_names()

    def _matrix_(self, R):
        r"""
        If self is a Maxima matrix, return the corresponding \sage
        matrix over the \sage ring $R$.

        This may or may not work depending in how complicated the
        entries of self are!  It only works if the entries of self
        can be coerced as strings to produce meaningful elements
        of $R$.

        EXAMPLES:
            sage: _ = maxima.eval("f[i,j] := i/j")
            sage: A = maxima('genmatrix(f,4,4)'); A
            matrix([1,1/2,1/3,1/4],[2,1,2/3,1/2],[3,3/2,1,3/4],[4,2,4/3,1])
            sage: A._matrix_(QQ)
            [  1 1/2 1/3 1/4]
            [  2   1 2/3 1/2]
            [  3 3/2   1 3/4]
            [  4   2 4/3   1]

        You can also use the \code{matrix} command (which is defined
        in \code{sage.misc.functional}):
            sage: matrix(QQ, A)
            [  1 1/2 1/3 1/4]
            [  2   1 2/3 1/2]
            [  3 3/2   1 3/4]
            [  4   2 4/3   1]
        """
        from sage.matrix.all import MatrixSpace
        self._check_valid()
        P = self.parent()
        nrows = int(P.eval('length(%s)'%self.name()))
        if nrows == 0:
            return MatrixSpace(R, 0, 0)(0)
        ncols = int(P.eval('length(%s[1])'%self.name()))
        M = MatrixSpace(R, nrows, ncols)
        s = self.str().replace('matrix','').replace(',',"','").\
            replace("]','[","','").replace('([',"['").replace('])',"']")
        s = eval(s)
        return M([R(x) for x in s])

    def partial_fraction_decomposition(self, var='x'):
        """
        Return the partial fraction decomposition of self with respect to
        the variable var.

        EXAMPLES:
            sage: f = maxima('1/((1+x)*(x-1))')
            sage: f.partial_fraction_decomposition('x')
            1/(2*(x - 1)) - 1/(2*(x + 1))
            sage: f.partial_fraction_decomposition('x').display2d()
                                 1           1
                             --------- - ---------
                             2 (x - 1)   2 (x + 1)
        """
        return self.partfrac(var)


class MaximaFunctionElement(FunctionElement):
    def _sage_doc_(self):
        return self._obj.parent().help(self._name)

class MaximaExpectFunction(ExpectFunction):
    def _sage_doc_(self):
        M = self._parent
        return M.help(self._name)


class MaximaFunction(MaximaElement):
    def __init__(self, parent, name, defn, args, latex):
        MaximaElement.__init__(self, parent, name, is_name=True)
        self.__defn = defn
        self.__args = args
        self.__latex = latex

    def __call__(self, *x):
        self._check_valid()
        P = self.parent()
        if len(x) == 1:
            x = '(%s)'%x
        return P('%s%s'%(self.name(), x))

    def __repr__(self):
        return self.__defn

    def _latex_(self):
        if self.__latex is None:
            return '\\mbox{\\rm %s}'%self.__defn
        else:
            return self.__latex

    def integrate(self, var):
        return self.integral(var)

    def integral(self, var):
        self._check_valid()
        P = self.parent()
        f = P('integrate(%s(%s), %s)'%(self.name(), self.__args, var))
        if var in self.__args.split(','):
            args = self.__args
        else:
            args = self.__args + ',' + var
        return P.function(args, str(f))


def is_MaximaElement(x):
    return isinstance(x, MaximaElement)

# An instance
maxima = Maxima(script_subdirectory=None)

def reduce_load_Maxima():
    return maxima

import os
def maxima_console():
    os.system('maxima')

def maxima_version():
    return os.popen('maxima --version').read().split()[1]

def __doctest_cleanup():
    import sage.interfaces.quit
    sage.interfaces.quit.expect_quitall()
