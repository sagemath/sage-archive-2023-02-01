r"""
Interface to Mathematica

The Mathematica interface will only work if Mathematica is installed on your
computer with a command line interface that runs when you give the ``math``
command. The interface lets you send certain Sage objects to Mathematica,
run Mathematica functions, import certain Mathematica expressions to Sage,
or any combination of the above.

To send a Sage object ``sobj`` to Mathematica, call ``mathematica(sobj)``.
This exports the Sage object to Mathematica and returns a new Sage object
wrapping the Mathematica expression/variable, so that you can use the
Mathematica variable from within Sage. You can then call Mathematica
functions on the new object; for example::

    sage: mobj = mathematica(x^2-1)             # optional - mathematica
    sage: mobj.Factor()                         # optional - mathematica
    (-1 + x)*(1 + x)

In the above example the factorization is done using Mathematica's
``Factor[]`` function.

To see Mathematica's output you can simply print the Mathematica wrapper
object. However if you want to import Mathematica's output back to Sage,
call the Mathematica wrapper object's ``sage()`` method. This method returns
a native Sage object::

    sage: mobj = mathematica(x^2-1)             # optional - mathematica
    sage: mobj2 = mobj.Factor(); mobj2          # optional - mathematica
    (-1 + x)*(1 + x)
    sage: mobj2.parent()                        # optional - mathematica
    Mathematica
    sage: sobj = mobj2.sage(); sobj             # optional - mathematica
    (x - 1)*(x + 1)
    sage: sobj.parent()                         # optional - mathematica
    Symbolic Ring


If you want to run a Mathematica function and don't already have the input
in the form of a Sage object, then it might be simpler to input a string to
``mathematica(expr)``. This string will be evaluated as if you had typed it
into Mathematica::

    sage: mathematica('Factor[x^2-1]')          # optional - mathematica
    (-1 + x)*(1 + x)
    sage: mathematica('Range[3]')               # optional - mathematica
    {1, 2, 3}

If you don't want Sage to go to the trouble of creating a wrapper for the
Mathematica expression, then you can call ``mathematica.eval(expr)``, which
returns the result as a Mathematica AsciiArtString formatted string. If you
want the result to be a string formatted like Mathematica's InputForm, call
``repr(mobj)`` on the wrapper object ``mobj``. If you want a string
formatted in Sage style, call ``mobj._sage_repr()``::

    sage: mathematica.eval('x^2 - 1')           # optional - mathematica
                   2
             -1 + x
    sage: repr(mathematica('Range[3]'))         # optional - mathematica
    '{1, 2, 3}'
    sage: mathematica('Range[3]')._sage_repr()  # optional - mathematica
    '[1, 2, 3]'

Finally, if you just want to use a Mathematica command line from within
Sage, the function ``mathematica_console()`` dumps you into an interactive
command-line Mathematica session. This is an enhanced version of the usual
Mathematica command-line, in that it provides readline editing and history
(the usual one doesn't!)

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
    {{x -> 0, y -> -1}, {x -> 6, y -> 11}}

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
    {{x -> 0, y -> -1}, {x -> 6, y -> 11}}

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

Next we find the minimum of a polynomial using the two different
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
    {{2, 2}, {3, 2}, {5, 3}, {541, 1}}
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

Complicated translations
------------------------

The ``mobj.sage()`` method tries to convert a Mathematica object to a Sage
object. In many cases, it will just work. In particular, it should be able to
convert expressions entirely consisting of:

- numbers, i.e. integers, floats, complex numbers;
- functions and named constants also present in Sage, where:

    - Sage knows how to translate the function or constant's name from
      Mathematica's, or
    - the Sage name for the function or constant is trivially related to
      Mathematica's;

- symbolic variables whose names don't pathologically overlap with
  objects already defined in Sage.

This method will not work when Mathematica's output includes:

- strings;
- functions unknown to Sage;
- Mathematica functions with different parameters/parameter order to
  the Sage equivalent.

If you want to convert more complicated Mathematica expressions, you can
instead call ``mobj._sage_()`` and supply a translation dictionary::

    sage: m = mathematica('NewFn[x]')       # optional - mathematica
    sage: m._sage_(locals={'NewFn': sin})   # optional - mathematica
    sin(x)

For more details, see the documentation for ``._sage_()``.


OTHER Examples::

    sage: def math_bessel_K(nu,x):
    ...       return mathematica(nu).BesselK(x).N(20)
    ...
    sage: math_bessel_K(2,I)                      # optional - mathematica
    0.180489972066962*I - 2.592886175491197             # 32-bit
    -2.59288617549119697817 + 0.18048997206696202663*I  # 64-bit

::

    sage: slist = [[1, 2], 3., 4 + I]
    sage: mlist = mathematica(slist); mlist     # optional - mathematica
    {{1, 2}, 3., 4 + I}
    sage: slist2 = list(mlist); slist2          # optional - mathematica
    [{1, 2}, 3., 4 + I]
    sage: slist2[0]                             # optional - mathematica
    {1, 2}
    sage: slist2[0].parent()                    # optional - mathematica
    Mathematica
    sage: slist3 = mlist.sage(); slist3         # optional - mathematica
    [[1, 2], 3.0, I + 4]

::

    sage: mathematica('10.^80')     # optional - mathematica
    1.*^80
    sage: mathematica('10.^80').sage()  # optional - mathematica
    1e+80

AUTHORS:

- William Stein (2005): first version

- Doug Cutrell (2006-03-01): Instructions for use under Cygwin/Windows.

- Felix Lawrence (2009-08-21): Added support for importing Mathematica lists
  and floats with exponents.
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
import re

from expect import (Expect, ExpectElement, ExpectFunction,
                    FunctionElement, AsciiArtString)

def clean_output(s):
    if s is None:
        return ''
    i = s.find('Out[')
    j = i + s[i:].find('=')
    s = s[:i] + ' '*(j+1-i) + s[j+1:]
    s = s.replace('\\\n','')
    return s.strip('\n')

def _un_camel(name):
    """
    Convert `CamelCase` to `camel_case`.

    EXAMPLES::

    sage: sage.interfaces.mathematica._un_camel('CamelCase')
    'camel_case'
    sage: sage.interfaces.mathematica._un_camel('EllipticE')
    'elliptic_e'
    sage: sage.interfaces.mathematica._un_camel('FindRoot')
    'find_root'
    sage: sage.interfaces.mathematica._un_camel('GCD')
    'gcd'
    """
    s1 = re.sub('(.)([A-Z][a-z]+)', r'\1_\2', name)
    return re.sub('([a-z0-9])([A-Z])', r'\1_\2', s1).lower()

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
              /Applications/Mathematica.app/Contents/MacOS/MathKernel $@

          The path in the above script must be modified if you installed
          Mathematica elsewhere or installed an old version of
          Mathematica that has the version in the .app name.

          (b) Make the file executable.
                chmod +x math

      * WINDOWS:

        Install Mathematica for Linux into the VMware virtual machine (sorry,
        that's the only way at present).
"""

##         The following only works with Sage for Cygwin (not colinux).
##         Note that Sage colinux is the preferred way to run Sage in Windows,
##         and I do not know how to use Mathematica from colinux Sage (unless
##         you install Mathematica-for-linux into the colinux machine, which
##         is possible).

##         Create a file named "math", which you place in the Sage root
##         directory.  The file contained a single line, which was the
##         path to the mathematica math.exe file.  In my case, this might be:

##         C:/Program Files/Wolfram Research/Mathematica/4.0/math.exe

##         The key points are
##         1) there is a file named "math.exe", and it will generally be
##            located in a place analogous to the above (depending on where
##            Mathematica has been installed).  This file is used only for
##            launching the kernel with a text-based interface.
##         2) a cygwin batch file must be created which executes this file,
##            which means using forward slashes rather than back slashes,
##            and probably surrounding everything in quotes
##         3) this cygwin batch file must be on the path for Sage (placing
##            it in <SAGE_LOCAL>/bin/ is an easy way to ensure this).

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
            raise TypeError("Error executing code in Mathematica\nCODE:\n\t%s\nMathematica ERROR:\n\t%s"%(cmd, out))

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

    def _eval_line(self, line,  allow_use_file=True, wait_for_prompt=True, restart_if_needed=False):
        s = Expect._eval_line(self, line,
             allow_use_file=allow_use_file, wait_for_prompt=wait_for_prompt)
        return str(s).strip('\n')

    def _function_call_string(self, function, args, kwds):
        """
        Returns the string used to make function calls.

        EXAMPLES::

            sage: mathematica._function_call_string('Sin', ['x'], [])
            'Sin[x]'
        """
        return "%s[%s]"%(function, ",".join(args))

    def _left_list_delim(self):
        return "{"

    def _right_list_delim(self):
        return "}"

    def _left_func_delim(self):
        return "["

    def _right_func_delim(self):
        return "]"

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

    def _exponent_symbol(self):
        """
        Returns the symbol used to denote the exponent of a number in
        Mathematica.

        EXAMPLES::

            sage: mathematica._exponent_symbol()        # optional - mathematica
            '*^'

        ::

            sage: bignum = mathematica('10.^80')  # optional - mathematica
            sage: repr(bignum)                          # optional - mathematica
            '1.*^80'
            sage: repr(bignum).replace(mathematica._exponent_symbol(), 'e').strip() # optional - mathematica
            '1.e80'
        """
        return "*^"

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
        return z[i+1:].strip()

    def __repr__(self):
        P = self._check_valid()
        return P.get(self._name, ascii_art=False).strip()

    def _sage_(self, locals={}):
        r"""
        Attempt to return a Sage version of this object.

        This method works successfully when Mathematica returns a result
        or list of results that consist only of:
        - numbers, i.e. integers, floats, complex numbers;
        - functions and named constants also present in Sage, where:
            - Sage knows how to translate the function or constant's name
            from Mathematica's naming scheme, or
            - you provide a translation dictionary `locals`, or
            - the Sage name for the function or constant is simply the
             Mathematica name in lower case;
        - symbolic variables whose names don't pathologically overlap with
          objects already defined in Sage.

        This method will not work when Mathematica's output includes:
        - strings;
        - functions unknown to Sage that are not specified in `locals`;
        - Mathematica functions with different parameters/parameter order to
          the Sage equivalent. In this case, define a function to do the
          parameter conversion, and pass it in via the locals dictionary.

        EXAMPLES:

        Mathematica lists of numbers/constants become Sage lists of
        numbers/constants::

            sage: m = mathematica('{{1., 4}, Pi, 3.2e100, I}')  # optional - mathematica
            sage: s = m.sage(); s       # optional - mathematica
            [[1.0, 4], pi, 3.2*e100, I]
            sage: s[1].n()              # optional - mathematica
            3.14159265358979
            sage: s[3]^2                # optional - mathematica
            -1

        ::

            sage: m = mathematica('x^2 + 5*y')      # optional - mathematica
            sage: m.sage()                          # optional - mathematica
            x^2 + 5*y

        ::

            sage: m = mathematica('Sin[Sqrt[1-x^2]] * (1 - Cos[1/x])^2')  # optional - mathematica
            sage: m.sage()                          # optional - mathematica
            (cos(1/x) - 1)^2*sin(sqrt(-x^2 + 1))

        ::

            sage: m = mathematica('NewFn[x]')       # optional - mathematica
            sage: m._sage_(locals={'NewFn': sin})   # optional - mathematica
            sin(x)

        ::

            sage: var('bla')                        # optional - mathematica
            bla
            sage: m = mathematica('bla^2')          # optional - mathematica
            sage: bla^2 - m.sage()                  # optional - mathematica
            0

        ::

            sage: m = mathematica('bla^2')          # optional - mathematica
            sage: mb = m.sage()                     # optional - mathematica
            sage: var('bla')                        # optional - mathematica
            bla
            sage: bla^2 - mb                        # optional - mathematica
            0


        AUTHORS:

        - Felix Lawrence (2010-11-03): Major rewrite to use ._sage_repr() and
          sage.calculus.calculus.symbolic_expression_from_string() for greater
          compatibility, while still supporting conversion of symbolic
          expressions.
        """
        from sage.symbolic.pynac import symbol_table
        from sage.symbolic.constants import constants_name_table as constants
        from sage.calculus.calculus import symbolic_expression_from_string
        from sage.calculus.calculus import _find_func as find_func

        # Get Mathematica's output and perform preliminary formatting
        res = self._sage_repr()
        if '"' in res:
            raise NotImplementedError("String conversion from Mathematica \
                does not work.  Mathematica's output was: %s" % res)

        # Find all the mathematica functions, constants and symbolic variables
        # present in `res`.  Convert MMA functions and constants to their
        # Sage equivalents (if possible), using `locals` and
        # `sage.symbolic.pynac.symbol_table['mathematica']` as translation
        # dictionaries.  If a MMA function or constant is not either
        # dictionary, then we use a variety of tactics listed in `autotrans`.
        # If a MMA variable is not in any dictionary, then create an
        # identically named Sage equivalent.

        # Merge the user-specified locals dictionary and the symbol_table
        # (locals takes priority)
        lsymbols = symbol_table['mathematica'].copy()
        lsymbols.update(locals)

        # Strategies for translating unknown functions/constants:
        autotrans = [   str.lower,      # Try it in lower case
                        _un_camel,    # Convert `CamelCase` to `camel_case`
                        lambda x: x     # Try the original name
                    ]

        # Find the MMA funcs/vars/constants - they start with a letter.
        # Exclude exponents (e.g. 'e8' from 4.e8)
        p = re.compile('(?<!\.)[a-zA-Z]\w*')
        for m in p.finditer(res):
            # If the function, variable or constant is already in the
            # translation dictionary, then just move on.
            if m.group() in lsymbols:
                pass
            # Now try to translate all other functions -- try each strategy
            # in `autotrans` and check if the function exists in Sage
            elif m.end() < len(res) and res[m.end()] == '(':
                for t in autotrans:
                    f = find_func(t(m.group()), create_when_missing = False)
                    if f is not None:
                        lsymbols[m.group()] = f
                        break
                else:
                    raise NotImplementedError("Don't know a Sage equivalent \
                        for Mathematica function '%s'.  Please specify one \
                        manually using the 'locals' dictionary" % m.group())
            # Check if Sage has an equivalent constant
            else:
                for t in autotrans:
                    if t(m.group()) in constants:
                        lsymbols[m.group()] = constants[t(m.group())]
                        break
            # If Sage has never heard of the variable, then
            # symbolic_expression_from_string will automatically create it
        try:
            return symbolic_expression_from_string(res, lsymbols,
                accept_sequence=True)
        except Exception:
            raise NotImplementedError("Unable to parse Mathematica \
                output: %s" % res)

    def __str__(self):
        P = self._check_valid()
        return P.get(self._name, ascii_art=True)

    def __len__(self):
        """
        Return the object's length, evaluated by mathematica.

        EXAMPLES::

            sage: len(mathematica([1,1.,2]))    # optional - mathematica
            3

        AUTHORS:
        - Felix Lawrence (2009-08-21)
        """
        return self.Length()

    def show(self, filename=None, ImageSize=600):
        r"""
        Show a mathematica expression or plot in the Sage notebook.

        EXAMPLES::

            sage: P = mathematica('Plot[Sin[x],{x,-2Pi,4Pi}]')   # optional - mathematica
            sage: show(P)                                        # optional - mathematica
            sage: P.show(ImageSize=800)                          # optional - mathematica
            sage: Q = mathematica('Sin[x Cos[y]]/Sqrt[1-x^2]')   # optional - mathematica
            sage: show(Q)                                        # optional - mathematica
            <html><div class="math">\frac{\sin (x \cos (y))}{\sqrt{1-x^2}}</div></html>
        """
        P = self._check_valid()
        if P.eval('InputForm[%s]' % self.name()).strip().startswith('Graphics['):
            if filename is None:
                from sage.misc.temporary_file import graphics_filename
                filename = graphics_filename()
            orig_dir = P.eval('Directory[]').strip()
            P.chdir(os.path.abspath("."))
            s = 'Export["%s", %s, ImageSize->%s]'%(filename, self.name(), ImageSize)
            P.eval(s)
            P.chdir(orig_dir)
        else:
            print '<html><div class="math">%s</div></html>' % self._latex_()

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

    def N(self, *args):
        """
        EXAMPLES::

            sage: mathematica('Pi').N(10)    # optional -- mathematica
            3.1415926536
            sage: mathematica('Pi').N(50)    # optional -- mathematica
            3.14159265358979323846264338327950288419716939937511
        """
        # The base class way up the hierarchy defines an "N" (modeled
        # after Mathematica's!)  which overwrites the Mathematica one,
        # and doesn't work at all. We restore it here.
        return self.parent().N(self, *args)


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
# is something the usual commercial mathematica doesn't provide!
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


