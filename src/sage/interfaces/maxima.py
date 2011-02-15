r"""
Interface to Maxima

Maxima is a free GPL'd general purpose computer algebra system
whose development started in 1968 at MIT. It contains symbolic
manipulation algorithms, as well as implementations of special
functions, including elliptic functions and generalized
hypergeometric functions. Moreover, Maxima has implementations of
many functions relating to the invariant theory of the symmetric
group `S_n`. (However, the commands for group invariants,
and the corresponding Maxima documentation, are in French.) For many
links to Maxima documentation see
http://maxima.sourceforge.net/docs.shtml/.

AUTHORS:

- William Stein (2005-12): Initial version

- David Joyner: Improved documentation

- William Stein (2006-01-08): Fixed bug in parsing

- William Stein (2006-02-22): comparisons (following suggestion of
  David Joyner)

- William Stein (2006-02-24): *greatly* improved robustness by adding
  sequence numbers to IO bracketing in _eval_line

If the string "error" (case insensitive) occurs in the output of
anything from Maxima, a RuntimeError exception is raised.

EXAMPLES: We evaluate a very simple expression in Maxima.

::

    sage: maxima('3 * 5')
    15

We factor `x^5 - y^5` in Maxima in several different ways.
The first way yields a Maxima object.

::

    sage: F = maxima.factor('x^5 - y^5')
    sage: F
    -(y-x)*(y^4+x*y^3+x^2*y^2+x^3*y+x^4)
    sage: type(F)
    <class 'sage.interfaces.maxima_abstract.MaximaElement'>

Note that Maxima objects can also be displayed using "ASCII art";
to see a normal linear representation of any Maxima object x. Just
use the print command: use ``str(x)``.

::

    sage: print F
                               4      3    2  2    3      4
                   - (y - x) (y  + x y  + x  y  + x  y + x )

You can always use ``repr(x)`` to obtain the linear
representation of an object. This can be useful for moving maxima
data to other systems.

::

    sage: repr(F)
    '-(y-x)*(y^4+x*y^3+x^2*y^2+x^3*y+x^4)'
    sage: F.str()
    '-(y-x)*(y^4+x*y^3+x^2*y^2+x^3*y+x^4)'

The ``maxima.eval`` command evaluates an expression in
maxima and returns the result as a *string* not a maxima object.

::

    sage: print maxima.eval('factor(x^5 - y^5)')
    -(y-x)*(y^4+x*y^3+x^2*y^2+x^3*y+x^4)

We can create the polynomial `f` as a Maxima polynomial,
then call the factor method on it. Notice that the notation
``f.factor()`` is consistent with how the rest of Sage
works.

::

    sage: f = maxima('x^5 - y^5')
    sage: f^2
    (x^5-y^5)^2
    sage: f.factor()
    -(y-x)*(y^4+x*y^3+x^2*y^2+x^3*y+x^4)

Control-C interruption works well with the maxima interface,
because of the excellent implementation of maxima. For example, try
the following sum but with a much bigger range, and hit control-C.

::

    sage: maxima('sum(1/x^2, x, 1, 10)')
    1968329/1270080

Tutorial
--------

We follow the tutorial at
http://maxima.sourceforge.net/docs/intromax/.

::

    sage: maxima('1/100 + 1/101')
    201/10100

::

    sage: a = maxima('(1 + sqrt(2))^5'); a
    (sqrt(2)+1)^5
    sage: a.expand()
    3*2^(7/2)+5*sqrt(2)+41

::

    sage: a = maxima('(1 + sqrt(2))^5')
    sage: float(a)
    82.012193308819747
    sage: a.numer()
    82.01219330881975

::

    sage: maxima.eval('fpprec : 100')
    '100'
    sage: a.bfloat()
    8.20121933088197564152489730020812442785204843859314941221237124017312418754011041266612384955016056b1

::

    sage: maxima('100!')
    93326215443944152681699238856266700490715968264381621468592963895217599993229915608941463976156518286253697920827223758251185210916864000000000000000000000000

::

    sage: f = maxima('(x + 3*y + x^2*y)^3')
    sage: f.expand()
    x^6*y^3+9*x^4*y^3+27*x^2*y^3+27*y^3+3*x^5*y^2+18*x^3*y^2+27*x*y^2+3*x^4*y+9*x^2*y+x^3
    sage: f.subst('x=5/z')
    (5/z+25*y/z^2+3*y)^3
    sage: g = f.subst('x=5/z')
    sage: h = g.ratsimp(); h
    (27*y^3*z^6+135*y^2*z^5+(675*y^3+225*y)*z^4+(2250*y^2+125)*z^3+(5625*y^3+1875*y)*z^2+9375*y^2*z+15625*y^3)/z^6
    sage: h.factor()
    (3*y*z^2+5*z+25*y)^3/z^6

::

    sage: eqn = maxima(['a+b*c=1', 'b-a*c=0', 'a+b=5'])
    sage: s = eqn.solve('[a,b,c]'); s
    [[a=(25*sqrt(79)*%i+25)/(6*sqrt(79)*%i-34),b=(5*sqrt(79)*%i+5)/(sqrt(79)*%i+11),c=(sqrt(79)*%i+1)/10],[a=(25*sqrt(79)*%i-25)/(6*sqrt(79)*%i+34),b=(5*sqrt(79)*%i-5)/(sqrt(79)*%i-11),c=-(sqrt(79)*%i-1)/10]]

Here is an example of solving an algebraic equation::

    sage: maxima('x^2+y^2=1').solve('y')
    [y=-sqrt(1-x^2),y=sqrt(1-x^2)]
    sage: maxima('x^2 + y^2 = (x^2 - y^2)/sqrt(x^2 + y^2)').solve('y')
    [y=-sqrt((-y^2-x^2)*sqrt(y^2+x^2)+x^2),y=sqrt((-y^2-x^2)*sqrt(y^2+x^2)+x^2)]

You can even nicely typeset the solution in latex::

    sage: latex(s)
    \left[ \left[ a={{25\,\sqrt{79}\,i+25}\over{6\,\sqrt{79}\,i-34}} ,   b={{5\,\sqrt{79}\,i+5}\over{\sqrt{79}\,i+11}} , c={{\sqrt{79}\,i+1  }\over{10}} \right]  , \left[ a={{25\,\sqrt{79}\,i-25}\over{6\,  \sqrt{79}\,i+34}} , b={{5\,\sqrt{79}\,i-5}\over{\sqrt{79}\,i-11}} ,   c=-{{\sqrt{79}\,i-1}\over{10}} \right]  \right]

To have the above appear onscreen via ``xdvi``, type
``view(s)``. (TODO: For OS X should create pdf output
and use preview instead?)

::

    sage: e = maxima('sin(u + v) * cos(u)^3'); e
    cos(u)^3*sin(v+u)
    sage: f = e.trigexpand(); f
    cos(u)^3*(cos(u)*sin(v)+sin(u)*cos(v))
    sage: f.trigreduce()
    (sin(v+4*u)+sin(v-2*u))/8+(3*sin(v+2*u)+3*sin(v))/8
    sage: w = maxima('3 + k*%i')
    sage: f = w^2 + maxima('%e')^w
    sage: f.realpart()
    %e^3*cos(k)-k^2+9

::

    sage: f = maxima('x^3 * %e^(k*x) * sin(w*x)'); f
    x^3*%e^(k*x)*sin(w*x)
    sage: f.diff('x')
    k*x^3*%e^(k*x)*sin(w*x)+3*x^2*%e^(k*x)*sin(w*x)+w*x^3*%e^(k*x)*cos(w*x)
    sage: f.integrate('x')
    (((k*w^6+3*k^3*w^4+3*k^5*w^2+k^7)*x^3+(3*w^6+3*k^2*w^4-3*k^4*w^2-3*k^6)*x^2+(-18*k*w^4-12*k^3*w^2+6*k^5)*x-6*w^4+36*k^2*w^2-6*k^4)*%e^(k*x)*sin(w*x)+((-w^7-3*k^2*w^5-3*k^4*w^3-k^6*w)*x^3+(6*k*w^5+12*k^3*w^3+6*k^5*w)*x^2+(6*w^5-12*k^2*w^3-18*k^4*w)*x-24*k*w^3+24*k^3*w)*%e^(k*x)*cos(w*x))/(w^8+4*k^2*w^6+6*k^4*w^4+4*k^6*w^2+k^8)

::

    sage: f = maxima('1/x^2')
    sage: f.integrate('x', 1, 'inf')
    1
    sage: g = maxima('f/sinh(k*x)^4')
    sage: g.taylor('x', 0, 3)
    f/(k^4*x^4)-2*f/(3*k^2*x^2)+11*f/45-62*k^2*f*x^2/945

::

    sage: maxima.taylor('asin(x)','x',0, 10)
    x+x^3/6+3*x^5/40+5*x^7/112+35*x^9/1152

Examples involving matrices
---------------------------

We illustrate computing with the matrix whose `i,j` entry
is `i/j`, for `i,j=1,\ldots,4`.

::

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
    [[[0,4],[3,1]],[[[1,0,0,-4],[0,1,0,-2],[0,0,1,-4/3]],[[1,2,3,4]]]]

We can also compute the echelon form in Sage::

    sage: B = matrix(QQ, A)
    sage: B.echelon_form()
    [  1 1/2 1/3 1/4]
    [  0   0   0   0]
    [  0   0   0   0]
    [  0   0   0   0]
    sage: B.charpoly('x').factor()
    (x - 4) * x^3

Laplace Transforms
------------------

We illustrate Laplace transforms::

    sage: _ = maxima.eval("f(t) := t*sin(t)")
    sage: maxima("laplace(f(t),t,s)")
    2*s/(s^2+1)^2

::

    sage: maxima("laplace(delta(t-3),t,s)") #Dirac delta function
    %e^-(3*s)

::

    sage: _ = maxima.eval("f(t) := exp(t)*sin(t)")
    sage: maxima("laplace(f(t),t,s)")
    1/(s^2-2*s+2)

::

    sage: _ = maxima.eval("f(t) := t^5*exp(t)*sin(t)")
    sage: maxima("laplace(f(t),t,s)")
    360*(2*s-2)/(s^2-2*s+2)^4-480*(2*s-2)^3/(s^2-2*s+2)^5+120*(2*s-2)^5/(s^2-2*s+2)^6
    sage: print maxima("laplace(f(t),t,s)")
                                             3                 5
               360 (2 s - 2)    480 (2 s - 2)     120 (2 s - 2)
              --------------- - --------------- + ---------------
                2           4     2           5     2           6
              (s  - 2 s + 2)    (s  - 2 s + 2)    (s  - 2 s + 2)

::

    sage: maxima("laplace(diff(x(t),t),t,s)")
    s*'laplace(x(t),t,s)-x(0)

::

    sage: maxima("laplace(diff(x(t),t,2),t,s)")
    -?%at('diff(x(t),t,1),t=0)+s^2*'laplace(x(t),t,s)-x(0)*s

It is difficult to read some of these without the 2d
representation::

    sage: print maxima("laplace(diff(x(t),t,2),t,s)")
                         !
                d        !         2
              - -- (x(t))!      + s  laplace(x(t), t, s) - x(0) s
                dt       !
                         !t = 0

Even better, use
``view(maxima("laplace(diff(x(t),t,2),t,s)"))`` to see
a typeset version.

Continued Fractions
-------------------

A continued fraction `a + 1/(b + 1/(c + \cdots))` is
represented in maxima by the list `[a, b, c, \ldots]`.

::

    sage: maxima("cf((1 + sqrt(5))/2)")
    [1,1,1,1,2]
    sage: maxima("cf ((1 + sqrt(341))/2)")
    [9,1,2,1,2,1,17,1,2,1,2,1,17,1,2,1,2,1,17,2]

Special examples
----------------

In this section we illustrate calculations that would be awkward to
do (as far as I know) in non-symbolic computer algebra systems like
MAGMA or GAP.

We compute the gcd of `2x^{n+4} - x^{n+2}` and
`4x^{n+1} + 3x^n` for arbitrary `n`.

::

    sage: f = maxima('2*x^(n+4) - x^(n+2)')
    sage: g = maxima('4*x^(n+1) + 3*x^n')
    sage: f.gcd(g)
    x^n

You can plot 3d graphs (via gnuplot)::

    sage: maxima('plot3d(x^2-y^2, [x,-2,2], [y,-2,2], [grid,12,12])')  # not tested
    [displays a 3 dimensional graph]

You can formally evaluate sums (note the ``nusum``
command)::

    sage: S = maxima('nusum(exp(1+2*i/n),i,1,n)')
    sage: print S
                            2/n + 3                   2/n + 1
                          %e                        %e
                   ----------------------- - -----------------------
                      1/n         1/n           1/n         1/n
                   (%e    - 1) (%e    + 1)   (%e    - 1) (%e    + 1)

We formally compute the limit as `n\to\infty` of
`2S/n` as follows::

    sage: T = S*maxima('2/n')
    sage: T.tlimit('n','inf')
    %e^3-%e

Miscellaneous
-------------

Obtaining digits of `\pi`::

    sage: maxima.eval('fpprec : 100')
    '100'
    sage: maxima(pi).bfloat()
    3.141592653589793238462643383279502884197169399375105820974944592307816406286208998628034825342117068b0

Defining functions in maxima::

    sage: maxima.eval('fun[a] := a^2')
    'fun[a]:=a^2'
    sage: maxima('fun[10]')
    100

Interactivity
-------------

Unfortunately maxima doesn't seem to have a non-interactive mode,
which is needed for the Sage interface. If any Sage call leads to
maxima interactively answering questions, then the questions can't be
answered and the maxima session may hang. See the discussion at
http://www.ma.utexas.edu/pipermail/maxima/2005/011061.html for some
ideas about how to fix this problem. An example that illustrates this
problem is ``maxima.eval('integrate (exp(a*x), x, 0, inf)')``.

Latex Output
------------

To TeX a maxima object do this::

    sage: latex(maxima('sin(u) + sinh(v^2)'))
    \sinh v^2+\sin u

Here's another example::

    sage: g = maxima('exp(3*%i*x)/(6*%i) + exp(%i*x)/(2*%i) + c')
    sage: latex(g)
    -{{i\,e^{3\,i\,x}}\over{6}}-{{i\,e^{i\,x}}\over{2}}+c

Long Input
----------

The MAXIMA interface reads in even very long input (using files) in
a robust manner, as long as you are creating a new object.

.. note::

   Using ``maxima.eval`` for long input is much less robust, and is
   not recommended.

::

    sage: t = '"%s"'%10^10000   # ten thousand character string.
    sage: a = maxima(t)

TESTS: This working tests that a subtle bug has been fixed::

    sage: f = maxima.function('x','gamma(x)')
    sage: g = f(1/7)
    sage: g
    gamma(1/7)
    sage: del f
    sage: maxima(sin(x))
    sin(x)

This tests to make sure we handle the case where Maxima asks if an
expression is positive or zero.

::

    sage: var('Ax,Bx,By')
    (Ax, Bx, By)
    sage: t = -Ax*sin(sqrt(Ax^2)/2)/(sqrt(Ax^2)*sqrt(By^2 + Bx^2))
    sage: t.limit(Ax=0, dir='+')
    0

A long complicated input expression::

    sage: maxima._eval_line('((((((((((0) + ((1) / ((n0) ^ (0)))) + ((1) / ((n1) ^ (1)))) + ((1) / ((n2) ^ (2)))) + ((1) / ((n3) ^ (3)))) + ((1) / ((n4) ^ (4)))) + ((1) / ((n5) ^ (5)))) + ((1) / ((n6) ^ (6)))) + ((1) / ((n7) ^ (7)))) + ((1) / ((n8) ^ (8)))) + ((1) / ((n9) ^ (9)));')
    '1/n9^9+1/n8^8+1/n7^7+1/n6^6+1/n5^5+1/n4^4+1/n3^3+1/n2^2+1/n1+1'
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

from __future__ import with_statement

import os, re, sys, subprocess
import pexpect
cygwin = os.uname()[0][:6]=="CYGWIN"

from expect import Expect, ExpectElement, FunctionElement, ExpectFunction, gc_disabled, AsciiArtString
from pexpect import EOF

from random import randrange

##import sage.rings.all
import sage.rings.complex_number

from sage.misc.misc import verbose, DOT_SAGE, SAGE_ROOT

from sage.misc.multireplace import multiple_replace

COMMANDS_CACHE = '%s/maxima_commandlist_cache.sobj'%DOT_SAGE

import sage.server.support

import maxima_abstract
from maxima_abstract import MaximaFunctionElement, MaximaExpectFunction, MaximaElement, MaximaFunction, maxima_console

# The Maxima "apropos" command, e.g., apropos(det) gives a list
# of all identifiers that begin in a certain way.  This could
# maybe be useful somehow... (?)  Also maxima has a lot for getting
# documentation from the system -- this could also be useful.

class Maxima(maxima_abstract.Maxima):
    """
    Interface to the Maxima interpreter.
    """
    def __init__(self, script_subdirectory=None, logfile=None, server=None,
                 init_code = None):
        """
        Create an instance of the Maxima interpreter.

        TESTS::

            sage: maxima == loads(dumps(maxima))
            True

        We make sure labels are turned off (see trac 6816)::

            sage: 'nolabels:true' in maxima._Expect__init_code
            True
        """
        # TODO: Input and output prompts in maxima can be changed by
        # setting inchar and outchar..
        eval_using_file_cutoff = 256
        self.__eval_using_file_cutoff = eval_using_file_cutoff
        STARTUP = '%s/local/bin/sage-maxima.lisp'%SAGE_ROOT

        # We set maxima's configuration directory to $DOT_SAGE/maxima
        # This avoids that sage's maxima inadvertently loads
        # ~/.maxima/maxima-init.mac
        # If you absolutely want maxima instances that are started by
        # this interface to preload commands, put them in
        # $DOT_SAGE/maxima/maxima-init.mac
        # (we use the "--userdir" option in maxima for this)
        import sage.misc.misc
        SAGE_MAXIMA_DIR = os.path.join(sage.misc.misc.DOT_SAGE,"maxima")

        if not os.path.exists(STARTUP):
            raise RuntimeError, 'You must get the file local/bin/sage-maxima.lisp'
        if init_code is None:
            # display2d -- no ascii art output
            # keepfloat -- don't automatically convert floats to rationals
            init_code = ['display2d : false', 'keepfloat : true']

        # Turn off the prompt labels, since computing them *very
        # dramatically* slows down the maxima interpret after a while.
        # See the function makelabel in suprv1.lisp.
        # Many thanks to andrej.vodopivec@gmail.com and also
        # Robert Dodier for figuring this out!
        # See trac # 6818.
        init_code.append('nolabels:true')



        Expect.__init__(self,
                        name = 'maxima',
                        prompt = '\(\%i[0-9]+\)',
                        command = 'maxima-noreadline --userdir="%s" -p "%s"'%(SAGE_MAXIMA_DIR,STARTUP),
                        maxread = 10000,
                        script_subdirectory = script_subdirectory,
                        restart_on_ctrlc = False,
                        verbose_start = False,
                        init_code = init_code,
                        logfile = logfile,
                        eval_using_file_cutoff=eval_using_file_cutoff)
        self._display_prompt = '<sage-display>'  # must match what is in the file local/bin/sage-maxima.lisp!!
        self._output_prompt_re = re.compile('\(\%o[0-9]+\)')
        self._ask = ['zero or nonzero?', 'an integer?', 'positive, negative, or zero?',
                     'positive or negative?', 'positive or zero?']
        self._prompt_wait = [self._prompt] + [re.compile(x) for x in self._ask] + \
                            ['Break [0-9]+'] #note that you might need to change _expect_expr if you
                                             #change this
        self._error_re = re.compile('(Principal Value|debugmode|incorrect syntax|Maxima encountered a Lisp error)')
        self._display2d = False



    def _function_class(self):
        """
        EXAMPLES::

            sage: maxima._function_class()
            <class 'sage.interfaces.maxima_abstract.MaximaExpectFunction'>
        """
        return MaximaExpectFunction

    def _start(self):
        """
        Starts the Maxima interpreter.

        EXAMPLES::

            sage: m = Maxima()
            sage: m.is_running()
            False
            sage: m._start()
            sage: m.is_running()
            True

        Test that we can use more than 256MB RAM (see trac 6722):

            sage: a = maxima(10)^(10^5)
            sage: b = a^600              # long time -- about 10-15 seconds

        """
        Expect._start(self)
        self._sendline(r":lisp (defun tex-derivative (x l r) (tex (if $derivabbrev (tex-dabbrev x) (tex-d x '\\partial)) l r lop rop ))")

        # Remove limit on the max heapsize (since otherwise it defaults
        # to 256MB with ECL).
        self._sendline(":lisp (ext:set-limit 'ext:heap-size 0)")
        self._eval_line('0;')

    def __reduce__(self):
        """
        EXAMPLES::

            sage: maxima.__reduce__()
            (<function reduce_load_Maxima at 0x...>, ())
        """
        return reduce_load_Maxima, tuple([])

    def _sendline(self, str):
        self._sendstr(str)
        os.write(self._expect.child_fd, os.linesep)

    def _expect_expr(self, expr=None, timeout=None):
        """
        EXAMPLES:
            sage: a,b=var('a,b')
            sage: integrate(1/(x^3 *(a+b*x)^(1/3)),x)
            Traceback (most recent call last):
            ...
            TypeError: Computation failed since Maxima requested additional constraints (try the command 'assume(a>0)' before integral or limit evaluation, for example):
            Is  a  positive or negative?
            sage: assume(a>0)
            sage: integrate(1/(x^3 *(a+b*x)^(1/3)),x)
            2/9*sqrt(3)*b^2*arctan(1/3*(2*(b*x + a)^(1/3) + a^(1/3))*sqrt(3)/a^(1/3))/a^(7/3) + 2/9*b^2*log((b*x + a)^(1/3) - a^(1/3))/a^(7/3) - 1/9*b^2*log((b*x + a)^(2/3) + (b*x + a)^(1/3)*a^(1/3) + a^(2/3))/a^(7/3) + 1/6*(4*(b*x + a)^(5/3)*b^2 - 7*(b*x + a)^(2/3)*a*b^2)/((b*x + a)^2*a^2 - 2*(b*x + a)*a^3 + a^4)
            sage: var('x, n')
            (x, n)
            sage: integral(x^n,x)
            Traceback (most recent call last):
            ...
            TypeError: Computation failed since Maxima requested additional constraints (try the command 'assume(n+1>0)' before integral or limit evaluation, for example):
            Is  n+1  zero or nonzero?
            sage: assume(n+1>0)
            sage: integral(x^n,x)
            x^(n + 1)/(n + 1)
            sage: forget()
        """
        if expr is None:
            expr = self._prompt_wait
        if self._expect is None:
            self._start()
        try:
            if timeout:
                i = self._expect.expect(expr,timeout=timeout)
            else:
                i = self._expect.expect(expr)
            if i > 0:
                v = self._expect.before

                #We check to see if there is a "serious" error in Maxima.
                #Note that this depends on the order of self._prompt_wait
                if expr is self._prompt_wait and i > len(self._ask):
                    self.quit()
                    raise ValueError, "%s\nComputation failed due to a bug in Maxima -- NOTE: Maxima had to be restarted."%v

                j = v.find('Is ')
                v = v[j:]
                k = v.find(' ',4)
                msg = "Computation failed since Maxima requested additional constraints (try the command 'assume(" + v[4:k] +">0)' before integral or limit evaluation, for example):\n" + v + self._ask[i-1]
                self._sendline(";")
                self._expect_expr()
                raise ValueError, msg
        except KeyboardInterrupt, msg:
            #print self._expect.before
            i = 0
            while True:
                try:
                    print "Control-C pressed.  Interrupting Maxima. Please wait a few seconds..."
                    self._sendstr('quit;\n'+chr(3))
                    self._sendstr('quit;\n'+chr(3))
                    self.interrupt()
                    self.interrupt()
                except KeyboardInterrupt:
                    i += 1
                    if i > 10:
                        break
                    pass
                else:
                    break
            raise KeyboardInterrupt, msg

    def _eval_line(self, line, allow_use_file=False,
                   wait_for_prompt=True, reformat=True, error_check=True):
        """
        EXAMPLES:

        We check that errors are correctly checked::

            sage: maxima._eval_line('1+1;')
            '2'
            sage: maxima._eval_line('sage0: x == x;')
            Traceback (most recent call last):
            ...
            TypeError: Error executing code in Maxima...


        """
        if len(line) == 0:
            return ''
        line = line.rstrip()
        if line[-1] != '$' and line[-1] != ';':
            line += ';'

        self._synchronize()

        if len(line) > self.__eval_using_file_cutoff:
            # This implicitly uses the set method, then displays the result of the thing that was set.
            # This only works when the input line is an expression.   But this is our only choice, since
            # batchmode doesn't display expressions to screen.
            a = self(line)
            return repr(a)
        else:
            self._sendline(line)

        line_echo = self._expect.readline()
        if not wait_for_prompt:
            return
        assert line_echo.strip() == line.strip()

        # This broke in maxima-5.22.1 as discussed in http://trac.sagemath.org/sage_trac/ticket/10187
        #self._expect_expr(self._display_prompt)
        #pre_out = self._before()
        #self._expect_expr()
        #out = self._before()
        #
        # if error_check:
        #     self._error_check(line, pre_out)
        #     self._error_check(line, out)
        #
        # if not reformat:
        #     return out
        #
        # r = self._output_prompt_re
        # m = r.search(out)
        # if m is None:
        #     o = out[:-2]
        # else:
        #     o = out[m.end()+1:-2]
        # o = ''.join([x.strip() for x in o.split()])
        # return o
        #
        # i = o.rfind('(%o')
        # return o[:i]

        self._expect_expr(self._display_prompt)
        out = self._before()        # input echo + output prompt + output
        if error_check:
            self._error_check(line, out)

        if not reformat:
            return out

        self._expect_expr()
        assert len(self._before())==0, 'Maxima expect interface is confused!'

        r = self._output_prompt_re
        m = r.search(out)
        if m is None:
            o = out[:-2]
        else:
            o = out[m.end()+1:-2]
        o = ''.join([x.strip() for x in o.split()])
        return o


    def _synchronize(self):
        """
        Synchronize pexpect interface.

        This put a random integer (plus one!) into the output stream, then
        waits for it, thus resynchronizing the stream. If the random
        integer doesn't appear within 1 second, maxima is sent interrupt
        signals.

        This way, even if you somehow left maxima in a busy state
        computing, calling _synchronize gets everything fixed.

        EXAMPLES: This makes Maxima start a calculation::

            sage: maxima._sendstr('1/1'*500)

        When you type this command, this synchronize command is implicitly
        called, and the above running calculation is interrupted::

            sage: maxima('2+2')
            4
        """
        marker = '__SAGE_SYNCHRO_MARKER_'
        if self._expect is None: return
        r = randrange(2147483647)
        s = marker + str(r+1)

        # The 0; *is* necessary... it comes up in certain rare cases
        # that are revealed by extensive testing.  Don't delete it. -- william stein
        cmd = '''0;sconcat("%s",(%s+1));\n'''%(marker,r)
        self._sendstr(cmd)
        try:
            self._expect_expr(timeout=0.5)
            if not s in self._before():
                self._expect_expr(s,timeout=0.5)
                self._expect_expr(timeout=0.5)
        except pexpect.TIMEOUT, msg:
            self._interrupt()
        except pexpect.EOF:
            self._crash_msg()
            self.quit()

    ###########################################
    # Direct access to underlying lisp interpreter.
    ###########################################
    def lisp(self, cmd):
        """
        Send a lisp command to maxima.

        .. note::

           The output of this command is very raw - not pretty.

        EXAMPLES::

            sage: maxima.lisp("(+ 2 17)")   # random formatted output
             :lisp (+ 2 17)
            19
            (
        """
        self._eval_line(':lisp %s\n""'%cmd, allow_use_file=False, wait_for_prompt=False, reformat=False, error_check=False)
        self._expect_expr('(%i)')
        return self._before()

    ###########################################
    # Interactive help
    ###########################################
    def _command_runner(self, command, s, redirect=True):
        """
        Run ``command`` in a new Maxima session and return its
        output as an ``AsciiArtString``.

        If redirect is set to False, then the output of the command is not
        returned as a string. Instead, it behaves like os.system. This is
        used for interactive things like Maxima's demos. See maxima.demo?

        EXAMPLES::

            sage: maxima._command_runner('describe', 'gcd')
            -- Function: gcd (<p_1>, <p_2>, <x_1>, ...)
            ...
        """
        cmd = 'maxima --very-quiet -r "%s(%s);" '%(command, s)
        if sage.server.support.EMBEDDED_MODE:
            cmd += '< /dev/null'

        if redirect:
            p = subprocess.Popen(cmd, shell=True, stdin=subprocess.PIPE,
                             stdout=subprocess.PIPE, stderr=subprocess.PIPE)
            res = p.stdout.read()
            # ecl-10.2 : 3 lines
            # ecl-10.4 : 5 lines
            # ecl-11.1 : 4 lines fancy a tango?
            # We now get 4 lines of commented verbosity
            # every time Maxima starts, so we need to get rid of them
            for _ in range(4):
                res = res[res.find('\n')+1:]
            return AsciiArtString(res)
        else:
            subprocess.Popen(cmd, shell=True)

    def _object_class(self):
        """
        Return the Python class of Maxima elements.

        EXAMPLES::

            sage: maxima._object_class()
            <class 'sage.interfaces.maxima_abstract.MaximaElement'>
        """
        return MaximaElement

    def _function_element_class(self):
        """
        EXAMPLES::

            sage: maxima._function_element_class()
            <class 'sage.interfaces.maxima_abstract.MaximaFunctionElement'>
        """
        return MaximaFunctionElement

    def function(self, args, defn, rep=None, latex=None):
        """
        Return the Maxima function with given arguments and definition.

        INPUT:


        -  ``args`` - a string with variable names separated by
           commas

        -  ``defn`` - a string (or Maxima expression) that
           defines a function of the arguments in Maxima.

        -  ``rep`` - an optional string; if given, this is how
           the function will print.


        EXAMPLES::

            sage: f = maxima.function('x', 'sin(x)')
            sage: f(3.2)
            -.058374143427580...
            sage: f = maxima.function('x,y', 'sin(x)+cos(y)')
            sage: f(2,3.5)
            sin(2)-.9364566872907963
            sage: f
            sin(x)+cos(y)

        ::

            sage: g = f.integrate('z')
            sage: g
            (cos(y)+sin(x))*z
            sage: g(1,2,3)
            3*(cos(2)+sin(1))

        The function definition can be a maxima object::

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
        if isinstance(defn, MaximaElement):
            defn = defn.str()
        elif not isinstance(defn, str):
            defn = str(defn)
        if isinstance(args, MaximaElement):
            args = args.str()
        elif not isinstance(args, str):
            args = str(args)
        cmd = '%s(%s) := %s'%(name, args, defn)
        maxima._eval_line(cmd)
        if rep is None:
            rep = defn
        f = MaximaFunction(self, name, rep, args, latex)
        return f

    def set(self, var, value):
        """
        Set the variable var to the given value.

        INPUT:


        -  ``var`` - string

        -  ``value`` - string


        EXAMPLES::

            sage: maxima.set('xxxxx', '2')
            sage: maxima.get('xxxxx')
            '2'
        """
        if not isinstance(value, str):
            raise TypeError
        cmd = '%s : %s$'%(var, value.rstrip(';'))
        if len(cmd) > self.__eval_using_file_cutoff:
            self._batch(cmd, batchload=True)
        else:
            self._eval_line(cmd)
            #self._sendline(cmd)
            #self._expect_expr()
            #out = self._before()
            #self._error_check(cmd, out)

    def clear(self, var):
        """
        Clear the variable named var.

        EXAMPLES::

            sage: maxima.set('xxxxx', '2')
            sage: maxima.get('xxxxx')
            '2'
            sage: maxima.clear('xxxxx')
            sage: maxima.get('xxxxx')
            'xxxxx'
        """
        try:
            self._expect.send('kill(%s)$'%var)
        except (TypeError, AttributeError):
             pass

    def get(self, var):
        """
        Get the string value of the variable var.

        EXAMPLES::

            sage: maxima.set('xxxxx', '2')
            sage: maxima.get('xxxxx')
            '2'
        """
        s = self._eval_line('%s;'%var)
        return s

    def version(self):
        """
        Return the version of Maxima that Sage includes.

        EXAMPLES::

            sage: maxima.version()
            '5.23.2'
        """
        return maxima_version()


##     def display2d(self, flag=True):
##         """
##         Set the flag that determines whether Maxima objects are
##         printed using their 2-d ASCII art representation.  When the
##         maxima interface starts the default is that objects are not
##         represented in 2-d.

##         INPUT:
##             flag -- bool (default: True)

##         EXAMPLES
##             sage: maxima('1/2')
##             1/2
##             sage: maxima.display2d(True)
##             sage: maxima('1/2')
##                                            1
##                                            -
##                                            2
##             sage: maxima.display2d(False)
##         """
##         self._display2d = bool(flag)

def is_MaximaElement(x):
    """
    Returns True if x is of type MaximaElement.

    EXAMPLES::

        sage: from sage.interfaces.maxima import is_MaximaElement
        sage: m = maxima(1)
        sage: is_MaximaElement(m)
        True
        sage: is_MaximaElement(1)
        False
    """
    return isinstance(x, MaximaElement)

# An instance
maxima = Maxima(init_code = ['display2d:false; domain: complex; keepfloat: true'],
                script_subdirectory=None)

def reduce_load_Maxima():
    """
    EXAMPLES::

        sage: from sage.interfaces.maxima import reduce_load_Maxima
        sage: reduce_load_Maxima()
        Maxima
    """
    return maxima


def maxima_version():
    """
    EXAMPLES::

        sage: from sage.interfaces.maxima import maxima_version
        sage: maxima_version()
        '5.23.2'
    """
    return os.popen('maxima --version').read().split()[-1]

def __doctest_cleanup():
    import sage.interfaces.quit
    sage.interfaces.quit.expect_quitall()


