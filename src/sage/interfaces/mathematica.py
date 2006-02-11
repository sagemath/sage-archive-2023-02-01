r"""
Interface to Mathematica

The Mathematica interface will only work if Mathematica is
installed on your computer with a command line interface that runs
when you give the \code{math} command.  The interface offers
three pieces of functionality:
\begin{enumerate}

\item \code{mathematica_console()} -- A function that dumps you
into an interactive command-line Mathematica session.   This is
an enhanced version of the usual Mathematica command-line, in
that it provides readline editing and history (the usual one
doesn't!)

\item \code{mathematica(expr)} -- Creation of a SAGE object that
wraps a Mathematica object.  This provides a Pythonic
interface to Mathematica.    For example, if
\code{f=mathematica('x\^2-1')}, then \code{f.Factor()}
returns the factorization of $x^2-1$ computed using Mathematica.

\item \code{mathematica.eval(expr)} -- Evaluation of arbitrary Mathematica
expressions, with the result returned as a string.


\end{enumerate}


\subsection{Tutorial}
We follow some of the tutorial from
   \url{http://library.wolfram.com/conferences/devconf99/withoff/Basic1.html/}.
For any of this to work you must buy and install the optional mathematica
program, and it must be available as the command \code{math} in your
PATH.  You do not have to install any special SAGE packages.

\subsubsection{Syntax}
Now make 1 and add it to itself.  The result is a Mathematica object.
    sage: m = mathematica
    sage: a = m(1) + m(1); a
    2
    sage: a.parent()
    Mathematica
    sage: m('1+1')
    2
    sage: m(3)**m(50)
    717897987691852588770249

The following is equivalent to \code{Plus[2, 3]} in Mathematica:
    sage: m = mathematica
    sage: m(2).Plus(m(3))
    5

We can also compute $7(2+3)$.
    sage: m(7).Times(m(2).Plus(m(3)))
    35
    sage: m('7(2+3)')
    35

\subsubsection{Some typical input}
We solve an equation and a system of two equations:

    sage: eqn = mathematica('3x + 5 == 14')
    sage: eqn
    5 + 3 x == 14
    sage: eqn.Solve('x')
    {{x -> 3}}
    sage: sys = mathematica('{x^2 - 3y == 3, 2x - y == 1}')
    sage: print sys
              2
            {x  - 3 y == 3, 2 x - y == 1}
    sage: sys.Solve('{x, y}')
    {{y -> -1, x -> 0}, {y -> 11, x -> 6}}


\subsubsection{Assignments and definitions}

If you assign the mathematica $5$ to a variable $c$ in SAGE,
this does not affect the $c$ in Mathematica.

    sage: c = m(5)
    sage: m('b + c x')
    b + c x
    sage: m('b') + c*m('x')
    b + 5 x

The SAGE interfaces changes SAGE lists into Mathematica lists:
    sage: m = mathematica
    sage: eq1 = m('x^2 - 3y == 3')
    sage: eq2 = m('2x - y == 1')
    sage: v = m([eq1, eq2])
    sage: print v
              2
            {x  - 3 y == 3, 2 x - y == 1}
    sage: v.Solve(['x', 'y'])
    {{y -> -1, x -> 0}, {y -> 11, x -> 6}}

\subsubsection{Function definitions}

Define mathematica functions by simply sending the definition to the
interpreter.

    sage: m = mathematica
    sage: _ = mathematica('f[p_] = p^2');
    sage: m('f[9]')
    81

\subsubsection{Numerical Calculations}

We find the $x$ such that $e^x - 3x = 0$.
    sage: e = mathematica('Exp[x] - 3x == 0')
    sage: e.FindRoot(['x', 2])
    {x -> 1.51213}

Note that this agrees with what the PARI interpreter gp produces:
    sage: gp('solve(x=1,2,exp(x)-3*x)')
    1.512134551657842473896739678

Next we find the minimimum of a polynomial using the two
different ways of accessing Mathematica:

    sage: mathematica.eval('FindMinimum[x^3 - 6x^2 + 11x - 5, {x,3}]')
    '{0.6151, {x -> 2.57735}}'
    sage: f = mathematica('x^3 - 6x^2 + 11x - 5')
    sage: f.FindMinimum(['x', 3])
    {0.6151, {x -> 2.57735}}


\subsubsection{Polynomial and Integer Factorization}

We factor a polynomial of degree 200 over the integers.

    sage: x = PolynomialRing(IntegerRing()).gen()
    sage: f = (x**100+17*x+5)*(x**100-5*x+20)
    sage: f
    x^200 + 12*x^101 + 25*x^100 - 85*x^2 + 315*x + 100
    sage: g = mathematica(str(f))
    sage: g
                              2       100       101    200
            100 + 315 x - 85 x  + 25 x    + 12 x    + x
    sage: g.Factor()
                         100               100
            (20 - 5 x + x   ) (5 + 17 x + x   )

We can also factor a multivariate polynomial:
    sage: f = mathematica('x^6 + (-y - 2)*x^5 + (y^3 + 2*y)*x^4 - y^4*x^3')
    sage: f.Factor()
             3                  2    3
            x  (x - y) (-2 x + x  + y )

We factor an integer:
    sage: n = mathematica(2434500)
    sage: n.FactorInteger()
    {{2, 2}, {3, 2}, {5, 3}, {541, 1}}
    sage: n = mathematica(2434500)
    sage: F = n.FactorInteger(); F
    {{2, 2}, {3, 2}, {5, 3}, {541, 1}}
    sage: F[1]
    {2, 2}
    sage: F[4]
    {541, 1}

We can also load the ECM package and factoring using it:
    sage: _ = mathematica.eval("<<NumberTheory`FactorIntegerECM`");
    sage: mathematica.FactorIntegerECM('932901*939321')
    8396109

%\subsection{Module Documentation}

\subsection{Long Input}
The Mathematica interface reads in even very long input (using files)
in a robust manner.

    sage: t = '"%s"'%10^10000   # ten thousand character string.
    sage: a = mathematica(t)
    sage: a = mathematica.eval(t)
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

import os

from expect import Expect, ExpectElement

from sage.misc.misc import verbose

def clean_output(s):
    i = s.find('Out[')
    j = i + s[i:].find('=')
    s = s[:i] + '\n       ' + s[j+1:]
    p = s.split()
    if len(p) == 1:
        return p[0]
    s = s.strip('\n')
    if len(s.split('\n')) == 1:
        s = s.strip()
    else:
        s = s.replace('\n\n','\n')
    s = s.replace('\\\n \n>    ','')
    return s.rstrip()


class Mathematica(Expect):
    """
    Interface to the Mathematica interpreter.
    """
    def __init__(self, maxread=100, script_subdirectory="", logfile=None, server=None):
        Expect.__init__(self,
                        name = 'mathematica',
                        prompt = ':=',
                        command = "math",
                        maxread = maxread,
                        server = server,
                        script_subdirectory = script_subdirectory,
                        verbose_start = False,
                        logfile=logfile,
                        eval_using_file_cutoff=50)

    def _read_in_file_command(self, filename):
        return '<<"%s"'%filename

    def _keyboard_interrupt(self):
        print "Interrupting %s..."%self
        e = self._expect
        e.sendline(chr(3))  # send ctrl-c
        e.expect('Interrupt> ')
        e.sendline("a")  # a -- abort
        e.expect(self._prompt)
        return e.before

    def eval(self, code, strip=True):
        s = Expect.eval(self, code)
        if strip:
            return clean_output(s)
        else:
            return s

    #def _keyboard_interrupt(self):
    #    print "Keyboard interrupt pressed; trying to recover."
    #    E = self.expect()
    #    E.sendline(chr(3))
    #    E.sendline('a')
    #    E.expect(':= ')
    #    raise KeyboardInterrupt, "Ctrl-c pressed while running Mathematica command"


    def set(self, var, value):
        """
        Set the variable var to the given value.
        """
        cmd = '%s=%s;'%(var,value)
        #out = self.eval(cmd)
        out = self._eval_line(cmd, allow_use_file=True)
        if len(out) > 8:
            raise TypeError, "Error executing code in Mathematica\nCODE:\n\t%s\nMathematica ERROR:\n\t%s"%(cmd, out)

    def get(self, var):
        """
        Get the value of the variable var.
        """
        return self.eval(var, strip=True)

    #def clear(self, var):
    #    """
    #    Clear the variable named var.
    #    """
    #    self.eval('Clear[%s]'%var)

    def _eval_line(self, line,  allow_use_file=True):
        s = Expect._eval_line(self, line,  allow_use_file=allow_use_file)
        return str(s).strip('\n')

    def function_call(self, function, args=[]):
        if not isinstance(args, list):
            args = [args]
        for i in range(len(args)):
            if not isinstance(args[i], ExpectElement):
                args[i] = self.new(args[i])
        return self.new("%s[%s]"%(function, ",".join([s.name() for s in args])))

    def _left_list_delim(self):
        return "{"

    def _right_list_delim(self):
        return "}"

    def _object_class(self):
        return MathematicaElement

    def console(self, readline=True):
        mathematica_console(readline=readline)


class MathematicaElement(ExpectElement):
    def __getitem__(self, n):
        return self.parent().new('%s[[%s]]'%(self._name, n))

    def __float__(self):
        P = self.parent()
        # TODO: Is 16 enough?
        return float(P.eval('N[%s,16]'%self.name()))


# An instance
mathematica = Mathematica(script_subdirectory='user')

# Cleverly run Mathematica with the benefit of readline, which
# is something the usual commerical mathematica doesn't provide!
# See
#    http://aspn.activestate.com/ASPN/Cookbook/Python/Recipe/363500

import os, sys
def mathematica_console(readline=True):
    if not readline:
        os.system('math')
        return
    f1 = os.popen('math ', 'w')
    f1.flush()
    try:
        while True:
            sys.stdout.write('')
            try:
                line = raw_input('In[ ]:= ')
                f1.writelines(line+'\n')
                f1.flush()
            except KeyboardInterrupt:
                f1.close()
                break
    except EOFError:
        pass
    sys.stdout.write('\n')

#def mathematica_console():
#    os.system('mathematica')


