r"""
Interface to Mathematica

The Mathematica interface will only work if Mathematica is
installed on your computer with a command line interface that runs
when you give the ``math`` command. The interface
offers three pieces of functionality:


#. ``mathematica_console()`` - A function that dumps
   you into an interactive command-line Mathematica session. This is
   an enhanced version of the usual Mathematica command-line, in that
   it provides readline editing and history (the usual one doesn't!)

#. ``mathematica(expr)`` - Creation of a Sage object
   that wraps a Mathematica object. This provides a Pythonic interface
   to Mathematica. For example, if
   ``f=mathematica('x2-1')``, then
   ``f.Factor()`` returns the factorization of
   `x^2-1` computed using Mathematica.

#. ``mathematica.eval(expr)`` - Evaluation of arbitrary
   Mathematica expressions, with the result returned as a string.


Tutorial
--------

We follow some of the tutorial from
http://library.wolfram.com/conferences/devconf99/withoff/Basic1.html/.

For any of this to work you must buy and install the Mathematica
program, and it must be available as the command
``math`` in your PATH.

Syntax
~~~~~~

Now make 1 and add it to itself. The result is a Mathematica
object.

::

    sage: m = mathematica
    sage: a = m(1) + m(1); a                # optional - mathematica
    2
    sage: a.parent()                        # optional - mathematica
    Mathematica
    sage: m('1+1')                          # optional - mathematica
    2
    sage: m(3)**m(50)                       # optional - mathematica
    717897987691852588770249

The following is equivalent to ``Plus[2, 3]`` in
Mathematica::

    sage: m = mathematica
    sage: m(2).Plus(m(3))                   # optional - mathematica
    5

We can also compute `7(2+3)`.

::

    sage: m(7).Times(m(2).Plus(m(3)))       # optional - mathematica
    35
    sage: m('7(2+3)')                       # optional - mathematica
    35

Some typical input
~~~~~~~~~~~~~~~~~~

We solve an equation and a system of two equations::

    sage: eqn = mathematica('3x + 5 == 14') # optional - mathematica
    sage: eqn                               # optional - mathematica
    5 + 3*x == 14
    sage: eqn.Solve('x')                    # optional - mathematica
    {{x -> 3}}
    sage: sys = mathematica('{x^2 - 3y == 3, 2x - y == 1}')  # optional - mathematica
    sage: print sys                         # optional - mathematica
               2
             {x  - 3 y == 3, 2 x - y == 1}
    sage: sys.Solve('{x, y}')               # optional - mathematica
    {{y -> -1, x -> 0}, {y -> 11, x -> 6}}

Assignments and definitions
~~~~~~~~~~~~~~~~~~~~~~~~~~~

If you assign the mathematica `5` to a variable `c`
in Sage, this does not affect the `c` in Mathematica.

::

    sage: c = m(5)                          # optional - mathematica
    sage: print m('b + c x')                # optional - mathematica
                 b + c x
    sage: print m('b') + c*m('x')           # optional - mathematica
             b + 5 x

The Sage interfaces changes Sage lists into Mathematica lists::

    sage: m = mathematica
    sage: eq1 = m('x^2 - 3y == 3')          # optional - mathematica
    sage: eq2 = m('2x - y == 1')            # optional - mathematica
    sage: v = m([eq1, eq2]); v              # optional - mathematica
    {x^2 - 3*y == 3, 2*x - y == 1}
    sage: v.Solve(['x', 'y'])               # optional - mathematica
    {{y -> -1, x -> 0}, {y -> 11, x -> 6}}

Function definitions
~~~~~~~~~~~~~~~~~~~~

Define mathematica functions by simply sending the definition to
the interpreter.

::

    sage: m = mathematica
    sage: _ = mathematica('f[p_] = p^2');   # optional - mathematica
    sage: m('f[9]')                         # optional - mathematica
    81

Numerical Calculations
~~~~~~~~~~~~~~~~~~~~~~

We find the `x` such that `e^x - 3x = 0`.

::

    sage: e = mathematica('Exp[x] - 3x == 0') # optional - mathematica
    sage: e.FindRoot(['x', 2])                # optional - mathematica
    {x -> 1.512134551657842}

Note that this agrees with what the PARI interpreter gp produces::

    sage: gp('solve(x=1,2,exp(x)-3*x)')
    1.512134551657842473896739678              # 32-bit
    1.5121345516578424738967396780720387046    # 64-bit

Next we find the minimimum of a polynomial using the two different
ways of accessing Mathematica::

    sage: mathematica('FindMinimum[x^3 - 6x^2 + 11x - 5, {x,3}]')  # optional - mathematica
    {0.6150998205402516, {x -> 2.5773502699629733}}
    sage: f = mathematica('x^3 - 6x^2 + 11x - 5')  # optional - mathematica
    sage: f.FindMinimum(['x', 3])                  # optional - mathematica
    {0.6150998205402516, {x -> 2.5773502699629733}}

Polynomial and Integer Factorization
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

We factor a polynomial of degree 200 over the integers.

::

    sage: R.<x> = PolynomialRing(ZZ)
    sage: f = (x**100+17*x+5)*(x**100-5*x+20)
    sage: f
    x^200 + 12*x^101 + 25*x^100 - 85*x^2 + 315*x + 100
    sage: g = mathematica(str(f))            # optional - mathematica
    sage: print g                            # optional - mathematica
                               2       100       101    200
             100 + 315 x - 85 x  + 25 x    + 12 x    + x
    sage: g                                  # optional - mathematica
    100 + 315*x - 85*x^2 + 25*x^100 + 12*x^101 + x^200
    sage: print g.Factor()                   # optional - mathematica
                          100               100
             (20 - 5 x + x   ) (5 + 17 x + x   )

We can also factor a multivariate polynomial::

    sage: f = mathematica('x^6 + (-y - 2)*x^5 + (y^3 + 2*y)*x^4 - y^4*x^3')  # optional - mathematica
    sage: print f.Factor()                   # optional - mathematica
              3                  2    3
             x  (x - y) (-2 x + x  + y )

We factor an integer::

    sage: n = mathematica(2434500)           # optional - mathematica
    sage: n.FactorInteger()                  # optional - mathematica
    {{2, 2}, {3, 2}, {5, 3}, {541, 1}}       # optional - mathematica
    sage: n = mathematica(2434500)           # optional - mathematica
    sage: F = n.FactorInteger(); F           # optional - mathematica
    {{2, 2}, {3, 2}, {5, 3}, {541, 1}}
    sage: F[1]                               # optional - mathematica
    {2, 2}
    sage: F[4]                               # optional - mathematica
    {541, 1}

We can also load the ECM package and factoring using it::

    sage: _ = mathematica.eval("<<NumberTheory`FactorIntegerECM`");  # optional - mathematica
    sage: mathematica.FactorIntegerECM('932901*939321')              # optional - mathematica
    8396109

Long Input
----------

The Mathematica interface reads in even very long input (using
files) in a robust manner.

::

    sage: t = '"%s"'%10^10000   # ten thousand character string.
    sage: a = mathematica(t)        # optional - mathematica
    sage: a = mathematica.eval(t)   # optional - mathematica

Loading and saving
------------------

Mathematica has an excellent ``InputForm`` function,
which makes saving and loading Mathematica objects possible. The
first examples test saving and loading to strings.

::

    sage: x = mathematica(pi/2)     # optional - mathematica
    sage: print x                   # optional - mathematica
             Pi
             --
             2
    sage: loads(dumps(x)) == x      # optional - mathematica
    True
    sage: n = x.N(50)               # optional - mathematica
    sage: print n                   # optional - mathematica
                  1.5707963267948966192313216916397514420985846996876
    sage: loads(dumps(n)) == n      # optional - mathematica
    True

OTHER Examples::

    sage: def math_bessel_K(nu,x):
    ...       return mathematica(nu).BesselK(x).N(20).sage()
    ...
    sage: math_bessel_K(2,I)                      # optional - mathematica
    0.180489972066962*I - 2.592886175491197

AUTHORS:

- William Stein (2005): first version

- Doug Cutrell (2006-03-01): Instructions for use under Cygwin/Windows.
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

import os

from expect import (Expect, ExpectElement, ExpectFunction,
                    FunctionElement, AsciiArtString)

from sage.misc.misc import verbose, graphics_filename

def clean_output(s):
    if s is None:
        return ''
    i = s.find('Out[')
    j = i + s[i:].find('=')
    s = s[:i] + ' '*(j+1-i) + s[j+1:]
    s = s.replace('\\\n','')
    return s.strip('\n')

class Mathematica(Expect):
    """
    Interface to the Mathematica interpreter.
    """
    def __init__(self, maxread=100, script_subdirectory="", logfile=None, server=None, server_tmpdir=None):
        Expect.__init__(self,
                        name = 'mathematica',
                        prompt = 'In[[0-9]+]:=',
                        command = "math",
                        maxread = maxread,
                        server = server,
                        server_tmpdir = server_tmpdir,
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

    def _install_hints(self):
        """
        Hints for installing mathematica on your computer.

        AUTHORS:

        - William Stein and Justin Walker (2006-02-12)
        """
        return """
In order to use the Mathematica interface you need to have Mathematica
installed and have a script in your PATH called "math" that runs the
command-line version of Mathematica. Alternatively, you could use a
remote connection to a server running Mathematica -- for hints, type
    print mathematica._install_hints_ssh()


  (1) You might have to buy Mathematica (http://www.wolfram.com/).

  (2) * LINUX: The math script comes standard with your Mathematica install.

      * APPLE OS X:
          (a) create a file called math (in your PATH):
              #!/bin/sh
              /Applications/Mathematica\ 5.2.app/Contents/MacOS/MathKernel $@

          Note that the 5.2 part will depend on the version of
          Mathematica you have, and the above path could be different
          if you installed mathematica elsewhere.

          (b) Make the file executable.
                chmod +x math

      * WINDOWS:

        Install Mathematica for Linux into the VMware virtual machine (sorry,
        that's the only way at present).
"""

##         The following only works with SAGE for Cygwin (not colinux).
##         Note that SAGE colinux is the preferred way to run SAGE in Windows,
##         and I do not know how to use mathematica from colinux SAGE (unless
##         you install Mathematica-for-linux into the colinux machine, which
##         is possible).

##         Create a file named "math", which you place in the SAGE root
##         directory.  The file contained a single line, which was the
##         path to the mathematica math.exe file.  In my case, this might be:

##         C:/Program Files/Wolfram Research/Mathematica/4.0/math.exe

##         The key points are
##         1) there is a file named "math.exe", and it will generally be
##            located in a place analagous to the above (depending on where
##            Mathematica has been installed).  This file is used only for
##            launching the kernel with a text-based interface.
##         2) a cygwin batch file must be created which executes this file,
##            which means using forward slashes rather than back slashes,
##            and probably surrounding everything in quotes
##         3) this cygwin batch file must be on the path for SAGE (placing
##            it in <SAGE_ROOT>/local/bin/ is an easy way to ensure this).

    def eval(self, code, strip=True, **kwds):
        s = Expect.eval(self, code, **kwds)
        if strip:
            return AsciiArtString(clean_output(s))
        else:
            return AsciiArtString(s)

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

    def get(self, var, ascii_art=False):
        """
        Get the value of the variable var.

        AUTHORS:

        - William Stein

        - Kiran Kedlaya (2006-02-04): suggested using InputForm
        """
        if ascii_art:
            return self.eval(var, strip=True)
        else:
            return self.eval('InputForm[%s, NumberMarks->False]'%var, strip=True)

    #def clear(self, var):
    #    """
    #    Clear the variable named var.
    #    """
    #    self.eval('Clear[%s]'%var)

    def _eval_line(self, line,  allow_use_file=True, wait_for_prompt=True):
        s = Expect._eval_line(self, line,
             allow_use_file=allow_use_file, wait_for_prompt=wait_for_prompt)
        return str(s).strip('\n')

    def function_call(self, function, args=None, kwds=None):
        args, kwds = self._convert_args_kwds(args, kwds)
        return self.new("%s[%s]"%(function, ",".join([s.name() for s in args])))

    def _left_list_delim(self):
        return "{"

    def _right_list_delim(self):
        return "}"

    ###########################################
    # System -- change directory, etc
    ###########################################
    def chdir(self, dir):
        """
        Change Mathematica's current working directory.

        EXAMPLES::

            sage: mathematica.chdir('/')          # optional
            sage: mathematica('Directory[]')      # optional
            "/"
        """
        self.eval('SetDirectory["%s"]'%dir)

    def _true_symbol(self):
        return '         True'

    def _false_symbol(self):
        return '         False'

    def _equality_symbol(self):
        return '=='

    def _assign_symbol(self):
        return ":="

    def _object_class(self):
        return MathematicaElement

    def console(self, readline=True):
        mathematica_console(readline=readline)

    def trait_names(self):
        a = self.eval('Names["*"]')
        return a.replace('$','').replace('\n \n>','').replace(',','').replace('}','').replace('{','').split()


    def help(self, cmd):
        return self.eval('? %s'%cmd)

    def __getattr__(self, attrname):
        if attrname[:1] == "_":
            raise AttributeError
        return MathematicaFunction(self, attrname)

class MathematicaElement(ExpectElement):
    def __getitem__(self, n):
        return self.parent().new('%s[[%s]]'%(self._name, n))

    def __getattr__(self, attrname):
        self._check_valid()
        if attrname[:1] == "_":
            raise AttributeError
        return MathematicaFunctionElement(self, attrname)

    def __float__(self):
        P = self.parent()
        # TODO: Is 16 enough?
        return float(P.eval('N[%s,16]'%self.name()))

    def _reduce(self):
        return self.parent().eval('InputForm[%s]'%self.name())

    def __reduce__(self):
        return reduce_load, (self._reduce(), )

    def _latex_(self):
        z = self.parent().eval('TeXForm[%s]'%self.name())
        i = z.find('=')
        return z[i+1:]

    def __repr__(self):
        P = self._check_valid()
        return P.get(self._name, ascii_art=False).strip()

    def __str__(self):
        P = self._check_valid()
        return P.get(self._name, ascii_art=True)

    def show(self, filename=None, ImageSize=600):
        """
        Show a mathematica plot in the Sage notebook.

        EXAMPLES::

            sage: P = mathematica('Plot[Sin[x],{x,-2Pi,4Pi}]')   # optional - mathematica
            sage: show(P)                                        # optional - mathematica
            sage: P.show(ImageSize=800)                          # optional - mathematica
        """
        P = self._check_valid()
        if filename is None:
            filename = graphics_filename()
        orig_dir = P.eval('Directory[]').strip()
        P.chdir(os.path.abspath("."))
        s = 'Export["%s", %s, ImageSize->%s]'%(filename, self.name(), ImageSize)
        P.eval(s)
        P.chdir(orig_dir)

    def str(self):
        return str(self)

    def __cmp__(self, other):
        #if not (isinstance(other, ExpectElement) and other.parent() is self.parent()):
        #    return coerce.cmp(self, other)
        P = self.parent()
        if P.eval("%s < %s"%(self.name(), other.name())).strip() == 'True':
            return -1
        elif P.eval("%s > %s"%(self.name(), other.name())).strip() == 'True':
            return 1
        elif P.eval("%s == %s"%(self.name(), other.name())).strip() == 'True':
            return 0
        else:
            return -1  # everything is supposed to be comparable in Python, so we define
                       # the comparison thus when no comparable in interfaced system.

class MathematicaFunction(ExpectFunction):
    def _sage_doc_(self):
        M = self._parent
        return M.help(self._name)


class MathematicaFunctionElement(FunctionElement):
    def _sage_doc_(self):
        M = self._obj.parent()
        return M.help(self._name)


# An instance
mathematica = Mathematica(script_subdirectory='user')

def reduce_load(X):
    return mathematica(X)

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
                line = raw_input('        ')
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


