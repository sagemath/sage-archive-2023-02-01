r"""
Interface to Singular

AUTHORS:

- David Joyner and William Stein (2005): first version

- Martin Albrecht (2006-03-05): code so singular.[tab] and x =
  singular(...), x.[tab] includes all singular commands.

- Martin Albrecht (2006-03-06): This patch adds the equality symbol to
  singular. Also fix a problem in which " " as prompt means comparison
  will break all further communication with Singular.

- Martin Albrecht (2006-03-13): added current_ring() and
  current_ring_name()

- William Stein (2006-04-10): Fixed problems with ideal constructor

- Martin Albrecht (2006-05-18): added sage_poly.

- Simon King (2010-11-23): Reduce the overhead caused by waiting for
  the Singular prompt by doing garbage collection differently.

- Simon King (2011-06-06): Make conversion from Singular to Sage more flexible.

Introduction
------------

This interface is extremely flexible, since it's exactly like
typing into the Singular interpreter, and anything that works there
should work here.

The Singular interface will only work if Singular is installed on
your computer; this should be the case, since Singular is included
with Sage. The interface offers three pieces of functionality:


#. ``singular_console()`` - A function that dumps you
   into an interactive command-line Singular session.

#. ``singular(expr, type='def')`` - Creation of a
   Singular object. This provides a Pythonic interface to Singular.
   For example, if ``f=singular(10)``, then
   ``f.factorize()`` returns the factorization of
   `10` computed using Singular.

#. ``singular.eval(expr)`` - Evaluation of arbitrary
   Singular expressions, with the result returned as a string.

Of course, there are polynomial rings and ideals in Sage as well
(often based on a C-library interface to Singular). One can convert
an object in the Singular interpreter interface to Sage by the
method ``sage()``.


Tutorial
--------

EXAMPLES: First we illustrate multivariate polynomial
factorization::

    sage: R1 = singular.ring(0, '(x,y)', 'dp')
    sage: R1
    //   characteristic : 0
    //   number of vars : 2
    //        block   1 : ordering dp
    //                  : names    x y
    //        block   2 : ordering C
    sage: f = singular('9x16 - 18x13y2 - 9x12y3 + 9x10y4 - 18x11y2 + 36x8y4 + 18x7y5 - 18x5y6 + 9x6y4 - 18x3y6 - 9x2y7 + 9y8')
    sage: f
    9*x^16-18*x^13*y^2-9*x^12*y^3+9*x^10*y^4-18*x^11*y^2+36*x^8*y^4+18*x^7*y^5-18*x^5*y^6+9*x^6*y^4-18*x^3*y^6-9*x^2*y^7+9*y^8
    sage: f.parent()
    Singular

::

    sage: F = f.factorize(); F
    [1]:
       _[1]=9
       _[2]=x^6-2*x^3*y^2-x^2*y^3+y^4
       _[3]=-x^5+y^2
    [2]:
       1,1,2

::

    sage: F[1]
    9,
    x^6-2*x^3*y^2-x^2*y^3+y^4,
    -x^5+y^2
    sage: F[1][2]
    x^6-2*x^3*y^2-x^2*y^3+y^4

We can convert `f` and each exponent back to Sage objects
as well.

::

    sage: g = f.sage(); g
    9*x^16 - 18*x^13*y^2 - 9*x^12*y^3 + 9*x^10*y^4 - 18*x^11*y^2 + 36*x^8*y^4 + 18*x^7*y^5 - 18*x^5*y^6 + 9*x^6*y^4 - 18*x^3*y^6 - 9*x^2*y^7 + 9*y^8
    sage: F[1][2].sage()
    x^6 - 2*x^3*y^2 - x^2*y^3 + y^4
    sage: g.parent()
    Multivariate Polynomial Ring in x, y over Rational Field

This example illustrates polynomial GCD's::

    sage: R2 = singular.ring(0, '(x,y,z)', 'lp')
    sage: a = singular.new('3x2*(x+y)')
    sage: b = singular.new('9x*(y2-x2)')
    sage: g = a.gcd(b)
    sage: g
    x^2+x*y

This example illustrates computation of a Groebner basis::

    sage: R3 = singular.ring(0, '(a,b,c,d)', 'lp')
    sage: I = singular.ideal(['a + b + c + d', 'a*b + a*d + b*c + c*d', 'a*b*c + a*b*d + a*c*d + b*c*d', 'a*b*c*d - 1'])
    sage: I2 = I.groebner()
    sage: I2
    c^2*d^6-c^2*d^2-d^4+1,
    c^3*d^2+c^2*d^3-c-d,
    b*d^4-b+d^5-d,
    b*c-b*d^5+c^2*d^4+c*d-d^6-d^2,
    b^2+2*b*d+d^2,
    a+b+c+d

The following example is the same as the one in the Singular - Gap
interface documentation::

    sage: R  = singular.ring(0, '(x0,x1,x2)', 'lp')
    sage: I1 = singular.ideal(['x0*x1*x2 -x0^2*x2', 'x0^2*x1*x2-x0*x1^2*x2-x0*x1*x2^2', 'x0*x1-x0*x2-x1*x2'])
    sage: I2 = I1.groebner()
    sage: I2
    x1^2*x2^2,
    x0*x2^3-x1^2*x2^2+x1*x2^3,
    x0*x1-x0*x2-x1*x2,
    x0^2*x2-x0*x2^2-x1*x2^2
    sage: I2.sage()
    Ideal (x1^2*x2^2, x0*x2^3 - x1^2*x2^2 + x1*x2^3, x0*x1 - x0*x2 - x1*x2, x0^2*x2 - x0*x2^2 - x1*x2^2) of Multivariate Polynomial Ring in x0, x1, x2 over Rational Field


This example illustrates moving a polynomial from one ring to
another. It also illustrates calling a method of an object with an
argument.

::

    sage: R = singular.ring(0, '(x,y,z)', 'dp')
    sage: f = singular('x3+y3+(x-y)*x2y2+z2')
    sage: f
    x^3*y^2-x^2*y^3+x^3+y^3+z^2
    sage: R1 = singular.ring(0, '(x,y,z)', 'ds')
    sage: f = R.fetch(f)
    sage: f
    z^2+x^3+y^3+x^3*y^2-x^2*y^3

We can calculate the Milnor number of `f`::

    sage: _=singular.LIB('sing.lib')     # assign to _ to suppress printing
    sage: f.milnor()
    4

The Jacobian applied twice yields the Hessian matrix of
`f`, with which we can compute.

::

    sage: H = f.jacob().jacob()
    sage: H
    6*x+6*x*y^2-2*y^3,6*x^2*y-6*x*y^2,  0,
    6*x^2*y-6*x*y^2,  6*y+2*x^3-6*x^2*y,0,
    0,                0,                2
    sage: H.sage()
    [6*x + 6*x*y^2 - 2*y^3     6*x^2*y - 6*x*y^2                     0]
    [    6*x^2*y - 6*x*y^2 6*y + 2*x^3 - 6*x^2*y                     0]
    [                    0                     0                     2]
    sage: H.det()   # This is a polynomial in Singular
    72*x*y+24*x^4-72*x^3*y+72*x*y^3-24*y^4-48*x^4*y^2+64*x^3*y^3-48*x^2*y^4
    sage: H.det().sage()   # This is the corresponding polynomial in Sage
    72*x*y + 24*x^4 - 72*x^3*y + 72*x*y^3 - 24*y^4 - 48*x^4*y^2 + 64*x^3*y^3 - 48*x^2*y^4

The 1x1 and 2x2 minors::

    sage: H.minor(1)
    2,
    6*y+2*x^3-6*x^2*y,
    6*x^2*y-6*x*y^2,
    6*x^2*y-6*x*y^2,
    6*x+6*x*y^2-2*y^3
    sage: H.minor(2)
    12*y+4*x^3-12*x^2*y,
    12*x^2*y-12*x*y^2,
    12*x^2*y-12*x*y^2,
    12*x+12*x*y^2-4*y^3,
    -36*x*y-12*x^4+36*x^3*y-36*x*y^3+12*y^4+24*x^4*y^2-32*x^3*y^3+24*x^2*y^4

::

    sage: _=singular.eval('option(redSB)')
    sage: H.minor(1).groebner()
    1

Computing the Genus
-------------------

We compute the projective genus of ideals that define curves over
`\QQ`. It is *very important* to load the
``normal.lib`` library before calling the
``genus`` command, or you'll get an error message.

EXAMPLE::

    sage: singular.lib('normal.lib')
    sage: R = singular.ring(0,'(x,y)','dp')
    sage: i2 = singular.ideal('y9 - x2*(x-1)^9 + x')
    sage: i2.genus()
    40

Note that the genus can be much smaller than the degree::

    sage: i = singular.ideal('y9 - x2*(x-1)^9')
    sage: i.genus()
    0

An Important Concept
--------------------

AUTHORS:

- Neal Harris

The following illustrates an important concept: how Sage interacts
with the data being used and returned by Singular. Let's compute a
Groebner basis for some ideal, using Singular through Sage.

::

    sage: singular.lib('poly.lib')
    sage: singular.ring(32003, '(a,b,c,d,e,f)', 'lp')
            //   characteristic : 32003
            //   number of vars : 6
            //        block   1 : ordering lp
            //                        : names    a b c d e f
            //        block   2 : ordering C
    sage: I = singular.ideal('cyclic(6)')
    sage: g = singular('groebner(I)')
    Traceback (most recent call last):
    ...
    TypeError: Singular error:
    ...

We restart everything and try again, but correctly.

::

    sage: singular.quit()
    sage: singular.lib('poly.lib'); R = singular.ring(32003, '(a,b,c,d,e,f)', 'lp')
    sage: I = singular.ideal('cyclic(6)')
    sage: I.groebner()
    f^48-2554*f^42-15674*f^36+12326*f^30-12326*f^18+15674*f^12+2554*f^6-1,
    ...

It's important to understand why the first attempt at computing a
basis failed. The line where we gave singular the input
'groebner(I)' was useless because Singular has no idea what 'I' is!
Although 'I' is an object that we computed with calls to Singular
functions, it actually lives in Sage. As a consequence, the name
'I' means nothing to Singular. When we called
``I.groebner()``, Sage was able to call the groebner
function on'I' in Singular, since 'I' actually means something to
Sage.

Long Input
----------

The Singular interface reads in even very long input (using files)
in a robust manner, as long as you are creating a new object.

::

    sage: t = '"%s"'%10^15000   # 15 thousand character string (note that normal Singular input must be at most 10000)
    sage: a = singular.eval(t)
    sage: a = singular(t)

TESTS:

We test an automatic coercion::

    sage: a = 3*singular('2'); a
    6
    sage: type(a)
    <class 'sage.interfaces.singular.SingularElement'>
    sage: a = singular('2')*3; a
    6
    sage: type(a)
    <class 'sage.interfaces.singular.SingularElement'>

Create a ring over GF(9) to check that ``gftables`` has been installed,
see ticket #11645::

    sage: singular.eval("ring testgf9 = (9,x),(a,b,c,d,e,f),(M((1,2,3,0)),wp(2,3),lp);")
    'ring testgf9 = (9,x),(a,b,c,d,e,f),(M((1,2,3,0)),wp(2,3),lp);'
"""

#We could also do these calculations without using the singular
#interface (behind the scenes the interface is used by Sage):
#    sage: x, y = PolynomialRing(RationalField(), 2, names=['x','y']).gens()
#    sage: C = ProjectivePlaneCurve(y**9 - x**2*(x-1)**9)
#    sage: C.genus()
#    0
#    sage: C = ProjectivePlaneCurve(y**9 - x**2*(x-1)**9 + x)
#    sage: C.genus()
#    40

#*****************************************************************************
#       Copyright (C) 2005 David Joyner and William Stein
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

import os, re, sys

from expect import Expect, ExpectElement, FunctionElement, ExpectFunction

from sage.structure.sequence import Sequence

from sage.structure.element import RingElement

import sage.rings.integer

from sage.misc.misc import get_verbose

class SingularError(RuntimeError):
    """
    Raised if Singular printed an error message
    """
    pass


class Singular(Expect):
    r"""
    Interface to the Singular interpreter.

    EXAMPLES: A Groebner basis example.

    ::

        sage: R = singular.ring(0, '(x0,x1,x2)', 'lp')
        sage: I = singular.ideal([ 'x0*x1*x2 -x0^2*x2', 'x0^2*x1*x2-x0*x1^2*x2-x0*x1*x2^2', 'x0*x1-x0*x2-x1*x2'])
        sage: I.groebner()
        x1^2*x2^2,
        x0*x2^3-x1^2*x2^2+x1*x2^3,
        x0*x1-x0*x2-x1*x2,
        x0^2*x2-x0*x2^2-x1*x2^2

    AUTHORS:

    - David Joyner and William Stein
    """
    def __init__(self, maxread=1000, script_subdirectory=None,
                 logfile=None, server=None,server_tmpdir=None):
        """
        EXAMPLES::

            sage: singular == loads(dumps(singular))
            True
        """
        prompt = '> '
        Expect.__init__(self,
                        terminal_echo=False,
                        name = 'singular',
                        prompt = prompt,
                        command = "Singular -t --ticks-per-sec 1000", #no tty and fine grained cputime()
                        maxread = maxread,
                        server = server,
                        server_tmpdir = server_tmpdir,
                        script_subdirectory = script_subdirectory,
                        restart_on_ctrlc = True,
                        verbose_start = False,
                        logfile = logfile,
                        eval_using_file_cutoff=100 if os.uname()[0]=="SunOS" else 1000)
        self.__libs  = []
        self._prompt_wait = prompt
        self.__to_clear = []   # list of variable names that need to be cleared.

    def _start(self, alt_message=None):
        """
        EXAMPLES::

            sage: s = Singular()
            sage: s.is_running()
            False
            sage: s._start()
            sage: s.is_running()
            True
            sage: s.quit()
        """
        self.__libs = []
        Expect._start(self, alt_message)
        # Load some standard libraries.
        self.lib('general')   # assumed loaded by misc/constants.py

        # these options are required by the new coefficient rings
        # supported by Singular 3-1-0.
        self.option("redTail")
        self.option("redThrough")
        self.option("intStrategy")
        self._saved_options = self.option('get')

    def __reduce__(self):
        """
        EXAMPLES::

            sage: singular.__reduce__()
            (<function reduce_load_Singular at 0x...>, ())
        """
        return reduce_load_Singular, ()

    def _equality_symbol(self):
        """
        EXAMPLES::

            sage: singular._equality_symbol()
            '=='
        """
        return '=='

    def _true_symbol(self):
        """
        EXAMPLES::

            sage: singular._true_symbol()
            '1'
        """
        return '1'

    def _false_symbol(self):
        """
        EXAMPLES::

            sage: singular._false_symbol()
            '0'
        """
        return '0'

    def _quit_string(self):
        """
        EXAMPLES::

            sage: singular._quit_string()
            'quit'
        """
        return 'quit'

    def _read_in_file_command(self, filename):
        r"""
        EXAMPLES::

            sage: singular._read_in_file_command('test')
            '< "...";'

            sage: filename = tmp_filename()
            sage: f = open(filename, 'w')
            sage: f.write('int x = 2;\n')
            sage: f.close()
            sage: singular.read(filename)
            sage: singular.get('x')
            '2'
        """
        return '< "%s";'%filename

    def eval(self, x, allow_semicolon=True, strip=True, **kwds):
        r"""
        Send the code x to the Singular interpreter and return the output
        as a string.

        INPUT:


        -  ``x`` - string (of code)

        -  ``allow_semicolon`` - default: False; if False then
           raise a TypeError if the input line contains a semicolon.

        -  ``strip`` - ignored


        EXAMPLES::

            sage: singular.eval('2 > 1')
            '1'
            sage: singular.eval('2 + 2')
            '4'

        if the verbosity level is `> 1` comments are also printed
        and not only returned.

        ::

            sage: r = singular.ring(0,'(x,y,z)','dp')
            sage: i = singular.ideal(['x^2','y^2','z^2'])
            sage: s = i.std()
            sage: singular.eval('hilb(%s)'%(s.name()))
            '// 1 t^0\n// -3 t^2\n// 3 t^4\n// -1 t^6\n\n// 1 t^0\n//
            3 t^1\n// 3 t^2\n// 1 t^3\n// dimension (affine) = 0\n//
            degree (affine) = 8'

        ::

            sage: set_verbose(1)
            sage: o = singular.eval('hilb(%s)'%(s.name()))
            //         1 t^0
            //        -3 t^2
            //         3 t^4
            //        -1 t^6
            //         1 t^0
            //         3 t^1
            //         3 t^2
            //         1 t^3
            // dimension (affine) = 0
            // degree (affine)  = 8

        This is mainly useful if this method is called implicitly. Because
        then intermediate results, debugging outputs and printed statements
        are printed

        ::

            sage: o = s.hilb()
            //         1 t^0
            //        -3 t^2
            //         3 t^4
            //        -1 t^6
            //         1 t^0
            //         3 t^1
            //         3 t^2
            //         1 t^3
            // dimension (affine) = 0
            // degree (affine)  = 8
            // ** right side is not a datum, assignment ignored

        rather than ignored

        ::

            sage: set_verbose(0)
            sage: o = s.hilb()
        """
        # Simon King:
        # In previous versions, the interface was first synchronised and then
        # unused variables were killed. This created a considerable overhead.
        # By trac ticket #10296, killing unused variables is now done inside
        # singular.set(). Moreover, it is not done by calling a separate _eval_line.
        # In that way, the time spent by waiting for the singular prompt is reduced.

        # Before #10296, it was possible that garbage collection occured inside
        # of _eval_line. But collection of the garbage would launch another call
        # to _eval_line. The result would have been a dead lock, that could only
        # be avoided by synchronisation. Since garbage collection is now done
        # without an additional call to _eval_line, synchronisation is not
        # needed anymore, saving even more waiting time for the prompt.

        # Uncomment the print statements below for low-level debugging of
        # code that involves the singular interfaces.  Everything goes
        # through here.
        #print "input: %s"%x
        x = str(x).rstrip().rstrip(';')
        x = x.replace("> ",">\t") #don't send a prompt  (added by Martin Albrecht)
        if not allow_semicolon and x.find(";") != -1:
            raise TypeError("singular input must not contain any semicolons:\n%s"%x)
        if len(x) == 0 or x[len(x) - 1] != ';':
            x += ';'

        s = Expect.eval(self, x, **kwds)

        if s.find("error") != -1 or s.find("Segment fault") != -1:
            raise SingularError('Singular error:\n%s'%s)

        if get_verbose() > 0:
            ret = []
            for line in s.splitlines():
                if line.startswith("//"):
                    print line
            return s
        else:
            return s

    def set(self, type, name, value):
        """
        Set the variable with given name to the given value.

        REMARK:

        If a variable in the Singular interface was previously marked for
        deletion, the actual deletion is done here, before the new variable
        is created in Singular.

        EXAMPLES::

            sage: singular.set('int', 'x', '2')
            sage: singular.get('x')
            '2'

        We test that an unused variable is only actually deleted if this method
        is called::

            sage: a = singular(3)
            sage: n = a.name()
            sage: del a
            sage: singular.eval(n)
            '3'
            sage: singular.set('int', 'y', '5')
            sage: singular.eval('defined(%s)'%n)
            '0'

        """
        cmd = ''.join('if(defined(%s)){kill %s;};'%(v,v) for v in self.__to_clear)
        cmd += '%s %s=%s;'%(type, name, value)
        self.__to_clear = []
        self.eval(cmd)

    def get(self, var):
        """
        Get string representation of variable named var.

        EXAMPLES::

            sage: singular.set('int', 'x', '2')
            sage: singular.get('x')
            '2'
        """
        return self.eval('print(%s);'%var)

    def clear(self, var):
        """
        Clear the variable named ``var``.

        EXAMPLES::

            sage: singular.set('int', 'x', '2')
            sage: singular.get('x')
            '2'
            sage: singular.clear('x')

        "Clearing the variable" means to allow to free the memory
        that it uses in the Singular sub-process. However, the
        actual deletion of the variable is only committed when
        the next element in the Singular interface is created::

            sage: singular.get('x')
            '2'
            sage: a = singular(3)
            sage: singular.get('x')
            '`x`'

        """
        # We add the variable to the list of vars to clear when we do an eval.
        # We queue up all the clears and do them at once to avoid synchronizing
        # the interface at the same time we do garbage collection, which can
        # lead to subtle problems.    This was Willem Jan's ideas, implemented
        # by William Stein.
        self.__to_clear.append(var)

    def _create(self, value, type='def'):
        """
        Creates a new variable in the Singular session and returns the name
        of that variable.

        EXAMPLES::

            sage: singular._create('2', type='int')
            'sage...'
            sage: singular.get(_)
            '2'
        """
        name = self._next_var_name()
        self.set(type, name, value)
        return name

    def __call__(self, x, type='def'):
        """
        Create a singular object X with given type determined by the string
        x. This returns var, where var is built using the Singular
        statement type var = ... x ... Note that the actual name of var
        could be anything, and can be recovered using X.name().

        The object X returned can be used like any Sage object, and wraps
        an object in self. The standard arithmetic operators work. Moreover
        if foo is a function then X.foo(y,z,...) calls foo(X, y, z, ...)
        and returns the corresponding object.

        EXAMPLES::

            sage: R = singular.ring(0, '(x0,x1,x2)', 'lp')
            sage: I = singular.ideal([ 'x0*x1*x2 -x0^2*x2', 'x0^2*x1*x2-x0*x1^2*x2-x0*x1*x2^2', 'x0*x1-x0*x2-x1*x2'])
            sage: I
             -x0^2*x2+x0*x1*x2,
            x0^2*x1*x2-x0*x1^2*x2-x0*x1*x2^2,
            x0*x1-x0*x2-x1*x2
            sage: type(I)
            <class 'sage.interfaces.singular.SingularElement'>
            sage: I.parent()
            Singular
        """
        if isinstance(x, SingularElement) and x.parent() is self:
            return x
        elif isinstance(x, ExpectElement):
            return self(x.sage())
        elif not isinstance(x, ExpectElement) and hasattr(x, '_singular_'):
            return x._singular_(self)

        # some convenient conversions
        if type in ("module","list") and isinstance(x,(list,tuple,Sequence)):
            x = str(x)[1:-1]

        return SingularElement(self, type, x, False)

    def has_coerce_map_from_impl(self, S):
        # we want to implement this without coercing, since singular has state.
        if hasattr(S, 'an_element'):
            if hasattr(S.an_element(), '_singular_'):
                return True
            try:
                self._coerce_(S.an_element())
            except TypeError:
                return False
            return True
        elif S is int or S is long:
            return True
        raise NotImplementedError


    def cputime(self, t=None):
        r"""
        Returns the amount of CPU time that the Singular session has used.
        If ``t`` is not None, then it returns the difference
        between the current CPU time and ``t``.

        EXAMPLES::

            sage: t = singular.cputime()
            sage: R = singular.ring(0, '(x0,x1,x2)', 'lp')
            sage: I = singular.ideal([ 'x0*x1*x2 -x0^2*x2', 'x0^2*x1*x2-x0*x1^2*x2-x0*x1*x2^2', 'x0*x1-x0*x2-x1*x2'])
            sage: gb = I.groebner()
            sage: singular.cputime(t) #random
            0.02
        """
        if t:
            return float(self.eval('timer-(%d)'%(int(1000*t))))/1000.0
        else:
            return float(self.eval('timer'))/1000.0

    ###################################################################
    # Singular libraries
    ###################################################################
    def lib(self, lib, reload=False):
        """
        Load the Singular library named lib.

        Note that if the library was already loaded during this session it
        is not reloaded unless the optional reload argument is True (the
        default is False).

        EXAMPLES::

            sage: singular.lib('sing.lib')
            sage: singular.lib('sing.lib', reload=True)
        """
        if lib[-4:] != ".lib":
            lib += ".lib"
        if not reload and lib in self.__libs:
            return
        self.eval('LIB "%s"'%lib)
        self.__libs.append(lib)

    LIB = lib
    load = lib

    ###################################################################
    # constructors
    ###################################################################
    def ideal(self, *gens):
        """
        Return the ideal generated by gens.

        INPUT:


        -  ``gens`` - list or tuple of Singular objects (or
           objects that can be made into Singular objects via evaluation)


        OUTPUT: the Singular ideal generated by the given list of gens

        EXAMPLES: A Groebner basis example done in a different way.

        ::

            sage: _ = singular.eval("ring R=0,(x0,x1,x2),lp")
            sage: i1 = singular.ideal([ 'x0*x1*x2 -x0^2*x2', 'x0^2*x1*x2-x0*x1^2*x2-x0*x1*x2^2', 'x0*x1-x0*x2-x1*x2'])
            sage: i1
            -x0^2*x2+x0*x1*x2,
            x0^2*x1*x2-x0*x1^2*x2-x0*x1*x2^2,
            x0*x1-x0*x2-x1*x2

        ::

            sage: i2 = singular.ideal('groebner(%s);'%i1.name())
            sage: i2
            x1^2*x2^2,
            x0*x2^3-x1^2*x2^2+x1*x2^3,
            x0*x1-x0*x2-x1*x2,
            x0^2*x2-x0*x2^2-x1*x2^2
        """
        if isinstance(gens, str):
            gens = self(gens)

        if isinstance(gens, SingularElement):
            return self(gens.name(), 'ideal')

        if not isinstance(gens, (list, tuple)):
            raise TypeError("gens (=%s) must be a list, tuple, string, or Singular element"%gens)

        if len(gens) == 1 and isinstance(gens[0], (list, tuple)):
            gens = gens[0]
        gens2 = []
        for g in gens:
            if not isinstance(g, SingularElement):
                gens2.append(self.new(g))
            else:
                gens2.append(g)
        return self(",".join([g.name() for g in gens2]), 'ideal')

    def list(self, x):
        r"""
        Creates a list in Singular from a Sage list ``x``.

        EXAMPLES::

            sage: singular.list([1,2])
            [1]:
               1
            [2]:
               2
        """
        return self(x, 'list')

    def matrix(self, nrows, ncols, entries=None):
        """
        EXAMPLES::

            sage: singular.lib("matrix")
            sage: R = singular.ring(0, '(x,y,z)', 'dp')
            sage: A = singular.matrix(3,2,'1,2,3,4,5,6')
            sage: A
            1,2,
            3,4,
            5,6
            sage: A.gauss_col()
            2,-1,
            1,0,
            0,1

        AUTHORS:

        - Martin Albrecht (2006-01-14)
        """
        name = self._next_var_name()
        if entries is None:
            s = self.eval('matrix %s[%s][%s]'%(name, nrows, ncols))
        else:
            s = self.eval('matrix %s[%s][%s] = %s'%(name, nrows, ncols, entries))
        return SingularElement(self, None, name, True)

    def ring(self, char=0, vars='(x)', order='lp', check=True):
        r"""
        Create a Singular ring and makes it the current ring.

        INPUT:


        -  ``char`` - characteristic of the base ring (see
           examples below), which must be either 0, prime (!), or one of
           several special codes (see examples below).

        -  ``vars`` - a tuple or string that defines the
           variable names

        -  ``order`` - string - the monomial order (default:
           'lp')

        -  ``check`` - if True, check primality of the
           characteristic if it is an integer.


        OUTPUT: a Singular ring

        .. note::

           This function is *not* identical to calling the Singular
           ``ring`` function. In particular, it also attempts to
           "kill" the variable names, so they can actually be used
           without getting errors, and it sets printing of elements
           for this range to short (i.e., with \*'s and carets).

        EXAMPLES: We first declare `\QQ[x,y,z]` with degree reverse
        lexicographic ordering.

        ::

            sage: R = singular.ring(0, '(x,y,z)', 'dp')
            sage: R
            //   characteristic : 0
            //   number of vars : 3
            //        block   1 : ordering dp
            //                  : names    x y z
            //        block   2 : ordering C

        ::

            sage: R1 = singular.ring(32003, '(x,y,z)', 'dp')
            sage: R2 = singular.ring(32003, '(a,b,c,d)', 'lp')

        This is a ring in variables named x(1) through x(10) over the
        finite field of order `7`::

            sage: R3 = singular.ring(7, '(x(1..10))', 'ds')

        This is a polynomial ring over the transcendental extension
        `\QQ(a)` of `\QQ`::

            sage: R4 = singular.ring('(0,a)', '(mu,nu)', 'lp')

        This is a ring over the field of single-precision floats::

            sage: R5 = singular.ring('real', '(a,b)', 'lp')

        This is over 50-digit floats::

            sage: R6 = singular.ring('(real,50)', '(a,b)', 'lp')
            sage: R7 = singular.ring('(complex,50,i)', '(a,b)', 'lp')

        To use a ring that you've defined, use the set_ring() method on
        the ring. This sets the ring to be the "current ring". For
        example,

        ::

            sage: R = singular.ring(7, '(a,b)', 'ds')
            sage: S = singular.ring('real', '(a,b)', 'lp')
            sage: singular.new('10*a')
            1.000e+01*a
            sage: R.set_ring()
            sage: singular.new('10*a')
            3*a
        """
        if len(vars) > 2:
            s = '; '.join(['if(defined(%s)>0){kill %s;};'%(x,x)
                           for x in vars[1:-1].split(',')])
            self.eval(s)

        if check and isinstance(char, (int, long, sage.rings.integer.Integer)):
            if char != 0:
                n = sage.rings.integer.Integer(char)
                if not n.is_prime():
                    raise ValueError("the characteristic must be 0 or prime")
        R = self('%s,%s,%s'%(char, vars, order), 'ring')
        self.eval('short=0')  # make output include *'s for multiplication for *THIS* ring.
        return R

    def string(self, x):
        """
        Creates a Singular string from a Sage string. Note that the Sage
        string has to be "double-quoted".

        EXAMPLES::

            sage: singular.string('"Sage"')
            Sage
        """
        return self(x, 'string')

    def set_ring(self, R):
        """
        Sets the current Singular ring to R.

        EXAMPLES::

            sage: R = singular.ring(7, '(a,b)', 'ds')
            sage: S = singular.ring('real', '(a,b)', 'lp')
            sage: singular.current_ring()
            //   characteristic : 0 (real)
            //   number of vars : 2
            //        block   1 : ordering lp
            //                  : names    a b
            //        block   2 : ordering C
            sage: singular.set_ring(R)
            sage: singular.current_ring()
            //   characteristic : 7
            //   number of vars : 2
            //        block   1 : ordering ds
            //                  : names    a b
            //        block   2 : ordering C
        """
        if not isinstance(R, SingularElement):
            raise TypeError("R must be a singular ring")
        self.eval("setring %s; short=0"%R.name(), allow_semicolon=True)

    setring = set_ring

    def current_ring_name(self):
        """
        Returns the Singular name of the currently active ring in
        Singular.

        OUTPUT: currently active ring's name

        EXAMPLES::

            sage: r = PolynomialRing(GF(127),3,'xyz')
            sage: r._singular_().name() == singular.current_ring_name()
            True
        """
        ringlist = self.eval("listvar(ring)").splitlines()
        p = re.compile("// ([a-zA-Z0-9_]*).*\[.*\].*\*.*") #do this in constructor?
        for line in ringlist:
            m = p.match(line)
            if m:
                return m.group(int(1))
        return None

    def current_ring(self):
        """
        Returns the current ring of the running Singular session.

        EXAMPLES::

            sage: r = PolynomialRing(GF(127),3,'xyz', order='invlex')
            sage: r._singular_()
            //   characteristic : 127
            //   number of vars : 3
            //        block   1 : ordering rp
            //                  : names    x y z
            //        block   2 : ordering C
            sage: singular.current_ring()
            //   characteristic : 127
            //   number of vars : 3
            //        block   1 : ordering rp
            //                  : names    x y z
            //        block   2 : ordering C
        """
        name = self.current_ring_name()
        if name:
            return self(name)
        else:
            return None

    def trait_names(self):
        """
         Return a list of all Singular commands.

         EXAMPLES::

             sage: singular.trait_names()
             ['exteriorPower',
              ...
              'stdfglm']
         """
        p = re.compile("// *([a-z0-9A-Z_]*).*") #compiles regular expression
        proclist = self.eval("listvar(proc)").splitlines()
        return [p.match(line).group(int(1)) for line in proclist]

    def console(self):
        """
        EXAMPLES::

            sage: singular_console() #not tested
                                 SINGULAR                             /  Development
             A Computer Algebra System for Polynomial Computations   /   version 3-0-4
                                                                   0<
                 by: G.-M. Greuel, G. Pfister, H. Schoenemann        \   Nov 2007
            FB Mathematik der Universitaet, D-67653 Kaiserslautern    \
        """
        singular_console()

    def version(self):
        """
        EXAMPLES:
        """
        return singular_version()

    def _function_class(self):
        """
        EXAMPLES::

            sage: singular._function_class()
            <class 'sage.interfaces.singular.SingularFunction'>
        """
        return SingularFunction

    def _function_element_class(self):
        """
        EXAMPLES::

            sage: singular._function_element_class()
            <class 'sage.interfaces.singular.SingularFunctionElement'>
        """
        return SingularFunctionElement

    def option(self, cmd=None, val=None):
        """
        Access to Singular's options as follows:

        Syntax: option() Returns a string of all defined options.

        Syntax: option( 'option_name' ) Sets an option. Note to disable an
        option, use the prefix no.

        Syntax: option( 'get' ) Returns an intvec of the state of all
        options.

        Syntax: option( 'set', intvec_expression ) Restores the state of
        all options from an intvec (produced by option('get')).

        EXAMPLES::

            sage: singular.option()
            //options: redefine loadLib usage prompt
            sage: singular.option('get')
            0,
            10321
            sage: old_options = _
            sage: singular.option('noredefine')
            sage: singular.option()
            //options: loadLib usage prompt
            sage: singular.option('set', old_options)
            sage: singular.option('get')
            0,
            10321
        """
        if cmd is None:
            return SingularFunction(self,"option")()
        elif cmd == "get":
            #return SingularFunction(self,"option")("\"get\"")
            return self(self.eval("option(get)"),"intvec")
        elif cmd == "set":
            if not isinstance(val,SingularElement):
                raise TypeError("singular.option('set') needs SingularElement as second parameter")
            #SingularFunction(self,"option")("\"set\"",val)
            self.eval("option(set,%s)"%val.name())
        else:
            SingularFunction(self,"option")("\""+str(cmd)+"\"")

    def _keyboard_interrupt(self):
        print "Interrupting %s..."%self
        try:
            self._expect.sendline(chr(4))
        except pexpect.ExceptionPexpect as msg:
            raise pexcept.ExceptionPexpect("THIS IS A BUG -- PLEASE REPORT. This should never happen.\n" + msg)
        self._start()
        raise KeyboardInterrupt("Restarting %s (WARNING: all variables defined in previous session are now invalid)"%self)

class SingularElement(ExpectElement):
    def __init__(self, parent, type, value, is_name=False):
        """
        EXAMPLES::

            sage: a = singular(2)
            sage: loads(dumps(a))
            (invalid object -- defined in terms of closed session)
        """
        RingElement.__init__(self, parent)
        if parent is None: return
        if not is_name:
            try:
                self._name = parent._create( value, type)
            # Convert SingularError to TypeError for
            # coercion to work properly.
            except SingularError as x:
                self._session_number = -1
                raise TypeError, x, sys.exc_info()[2]
            except BaseException:
                self._session_number = -1
                raise
        else:
            self._name = value
        self._session_number = parent._session_number

    def __repr__(self):
        r"""
        Return string representation of ``self``.

        EXAMPLE::

            sage: r = singular.ring(0,'(x,y)','dp')
            sage: singular(0)
            0
            sage: singular('x') # indirect doctest
            x
            sage: singular.matrix(2,2)
            0,0,
            0,0
            sage: singular.matrix(2,2,"(25/47*x^2*y^4 + 63/127*x + 27)^3,y,0,1")
            15625/103823*x^6*y.., y,
            0,                    1

        Note that the output is truncated

        ::

            sage: M= singular.matrix(2,2,"(25/47*x^2*y^4 + 63/127*x + 27)^3,y,0,1")
            sage: M.rename('T')
            sage: M
            T[1,1],y,
            0,         1

        if ``self`` has a custom name, it is used to print the
        matrix, rather than abbreviating its contents
        """
        try:
            self._check_valid()
        except ValueError:
            return '(invalid object -- defined in terms of closed session)'
        try:
            if self._get_using_file:
                s = self.parent().get_using_file(self._name)
        except AttributeError:
            s = self.parent().get(self._name)
        if s.__contains__(self._name):
            if hasattr(self, '__custom_name'):
                s =  s.replace(self._name, self.__dict__['__custom_name'])
            elif self.type() == 'matrix':
                s = self.parent().eval('pmat(%s,20)'%(self.name()))
        return s

    def __copy__(self):
        r"""
        Returns a copy of ``self``.

        EXAMPLES::

            sage: R=singular.ring(0,'(x,y)','dp')
            sage: M=singular.matrix(3,3,'0,0,-x, 0,y,0, x*y,0,0')
            sage: N=copy(M)
            sage: N[1,1]=singular('x+y')
            sage: N
            x+y,0,-x,
            0,  y,0,
            x*y,0,0
            sage: M
            0,  0,-x,
            0,  y,0,
            x*y,0,0
            sage: L=R.ringlist()
            sage: L[4]=singular.ideal('x**2-5')
            sage: Q=L.ring()
            sage: otherR=singular.ring(5,'(x)','dp')
            sage: cpQ=copy(Q)
            sage: cpQ.set_ring()
            sage: cpQ
            //   characteristic : 0
            //   number of vars : 2
            //        block   1 : ordering dp
            //                  : names    x y
            //        block   2 : ordering C
            // quotient ring from ideal
            _[1]=x^2-5
            sage: R.fetch(M)
            0,  0,-x,
            0,  y,0,
            x*y,0,0
        """
        if (self.type()=='ring') or (self.type()=='qring'):
            # Problem: singular has no clean method to produce
            # a copy of a ring/qring. We use ringlist, but this
            # is only possible if we make self the active ring,
            # use ringlist, and switch back to the previous
            # base ring.
            br=self.parent().current_ring()
            self.set_ring()
            OUT = (self.ringlist()).ring()
            br.set_ring()
            return OUT
        else:
            return self.parent()(self.name())

    def __len__(self):
        """
        Returns the size of this Singular element.

        EXAMPLES::

            sage: R = singular.ring(0, '(x,y,z)', 'dp')
            sage: A = singular.matrix(2,2)
            sage: len(A)
            4
        """
        return int(self.size())

    def __reduce__(self):
        """
        Note that the result of the returned reduce_load is an invalid
        Singular object.

        EXAMPLES::

            sage: singular(2).__reduce__()
            (<function reduce_load at 0x...>, ())
        """
        return reduce_load, ()  # default is an invalid object

    def __setitem__(self, n, value):
        """
        Set the n-th element of self to x.

        INPUT:


        -  ``n`` - an integer *or* a 2-tuple (for setting
           matrix elements)

        -  ``value`` - anything (is coerced to a Singular
           object if it is not one already)


        OUTPUT: Changes elements of self.

        EXAMPLES::

            sage: R = singular.ring(0, '(x,y,z)', 'dp')
            sage: A = singular.matrix(2,2)
            sage: A
            0,0,
            0,0
            sage: A[1,1] = 5
            sage: A
            5,0,
            0,0
            sage: A[1,2] = '5*x + y + z3'
            sage: A
            5,z^3+5*x+y,
            0,0
        """
        P = self.parent()
        if not isinstance(value, SingularElement):
            value = P(value)
        if isinstance(n, tuple):
            if len(n) != 2:
                raise ValueError("If n (=%s) is a tuple, it must be a 2-tuple"%n)
            x, y = n
            P.eval('%s[%s,%s] = %s'%(self.name(), x, y, value.name()))
        else:
            P.eval('%s[%s] = %s'%(self.name(), n, value.name()))

    def __nonzero__(self):
        """
        Returns True if this Singular element is not zero.

        EXAMPLES::

            sage: singular(0).__nonzero__()
            False
            sage: singular(1).__nonzero__()
            True
        """
        P = self.parent()
        return P.eval('%s == 0'%self.name()) == '0'

    def sage_polystring(self):
        r"""
        If this Singular element is a polynomial, return a string
        representation of this polynomial that is suitable for evaluation
        in Python. Thus \* is used for multiplication and \*\* for
        exponentiation. This function is primarily used internally.

        The short=0 option *must* be set for the parent ring or this
        function will not work as expected. This option is set by default
        for rings created using ``singular.ring`` or set using
        ``ring_name.set_ring()``.

        EXAMPLES::

            sage: R = singular.ring(0,'(x,y)')
            sage: f = singular('x^3 + 3*y^11 + 5')
            sage: f
            x^3+3*y^11+5
            sage: f.sage_polystring()
            'x**3+3*y**11+5'
        """
        return str(self).replace('^','**')

    def sage_global_ring(self):
        """
        Return the current basering in Singular as a polynomial ring or quotient ring.

        EXAMPLE::

            sage: singular.eval('ring r1 = (9,x),(a,b,c,d,e,f),(M((1,2,3,0)),wp(2,3),lp)')
            'ring r1 = (9,x),(a,b,c,d,e,f),(M((1,2,3,0)),wp(2,3),lp);'
            sage: R = singular('r1').sage_global_ring()
            sage: R
            Multivariate Polynomial Ring in a, b, c, d, e, f over Finite Field in x of size 3^2
            sage: R.term_order()
            Block term order with blocks:
            (Matrix term order with matrix
            [1 2]
            [3 0],
             Weighted degree reverse lexicographic term order with weights (2, 3),
             Lexicographic term order of length 2)

        ::

            sage: singular.eval('ring r2 = (0,x),(a,b,c),dp')
            'ring r2 = (0,x),(a,b,c),dp;'
            sage: singular('r2').sage_global_ring()
            Multivariate Polynomial Ring in a, b, c over Fraction Field of Univariate Polynomial Ring in x over Rational Field

        ::

            sage: singular.eval('ring r3 = (3,z),(a,b,c),dp')
            'ring r3 = (3,z),(a,b,c),dp;'
            sage: singular.eval('minpoly = 1+z+z2+z3+z4')
            'minpoly = 1+z+z2+z3+z4;'
            sage: singular('r3').sage_global_ring()
            Multivariate Polynomial Ring in a, b, c over Finite Field in z of size 3^4

        Real and complex fields in both Singular and Sage are defined with a precision.
        The precision in Singular is given in terms of digits, but in Sage it is given
        in terms of bits. So, the digit precision is internally converted to a reasonable
        bit precision::

            sage: singular.eval('ring r4 = (real,20),(a,b,c),dp')
            'ring r4 = (real,20),(a,b,c),dp;'
            sage: singular('r4').sage_global_ring()
            Multivariate Polynomial Ring in a, b, c over Real Field with 70 bits of precision

        The case of complex coefficients is not fully supported, yet, since
        the generator of a complex field in Sage is always called "I"::

            sage: singular.eval('ring r5 = (complex,15,j),(a,b,c),dp')
            'ring r5 = (complex,15,j),(a,b,c),dp;'
            sage: R = singular('r5').sage_global_ring(); R
            Multivariate Polynomial Ring in a, b, c over Complex Field with 54 bits of precision
            sage: R.base_ring()('j')
            Traceback (most recent call last):
            ...
            NameError: name 'j' is not defined
            sage: R.base_ring()('I')
            1.00000000000000*I

        In our last example, the base ring is a quotient ring::

            sage: singular.eval('ring r6 = (9,a), (x,y,z),lp')
            'ring r6 = (9,a), (x,y,z),lp;'
            sage: Q = singular('std(ideal(x^2,x+y^2+z^3))', type='qring')
            sage: Q.sage_global_ring()
            Quotient of Multivariate Polynomial Ring in x, y, z over Finite Field in a of size 3^2 by the ideal (y^4 - y^2*z^3 + z^6, x + y^2 + z^3)

        AUTHOR:

        - Simon King (2011-06-06)

        """
        # extract the ring of coefficients
        singular = self.parent()
        charstr = singular.eval('charstr(basering)').split(',',1)
        from sage.all import ZZ
        is_extension = len(charstr)==2
        if charstr[0]=='integer':
            br = ZZ
            is_extension = False
        elif charstr[0]=='0':
            from sage.all import QQ
            br = QQ
        elif charstr[0]=='real':
            from sage.all import RealField, ceil, log
            prec = singular.eval('ringlist(basering)[1][2][1]')
            br = RealField(ceil((ZZ(prec)+1)/log(2,10)))
            is_extension = False
        elif charstr[0]=='complex':
            from sage.all import ComplexField, ceil, log
            prec = singular.eval('ringlist(basering)[1][2][1]')
            I = singular.eval('ringlist(basering)[1][3]')
            br = ComplexField(ceil((ZZ(prec)+1)/log(2,10)))
            is_extension = False
        else:
            # it ought to be a finite field
            q = ZZ(charstr[0])
            from sage.all import GF
            if q.is_prime():
                br = GF(q)
            else:
                br = GF(q,charstr[1])
                # Singular has no extension of a non-prime field
                is_extension = False

        # We have the base ring of the base ring. But is it
        # an extension?
        if is_extension:
            minpoly = singular.eval('minpoly')
            if minpoly == '0':
                from sage.all import Frac
                BR = Frac(br[charstr[1]])
            else:
                is_short = singular.eval('short')
                if is_short!='0':
                    singular.eval('short=0')
                    minpoly = ZZ[charstr[1]](singular.eval('minpoly'))
                    singular.eval('short=%s'%is_short)
                else:
                    minpoly = ZZ[charstr[1]](minpoly)
                BR = br.extension(minpoly,name=charstr[1])
        else:
            BR = br

        # Now, we form the polynomial ring over BR with the given variables,
        # using Singular's term order
        from sage.rings.polynomial.term_order import termorder_from_singular
        from sage.all import PolynomialRing
        if singular.eval('typeof(basering)')=='ring':
            return PolynomialRing(BR, names=singular.eval('varstr(basering)'), order=termorder_from_singular(singular))
        P = PolynomialRing(BR, names=singular.eval('varstr(basering)'), order=termorder_from_singular(singular))
        return P.quotient(singular('ringlist(basering)[4]')._sage_(P), names=singular.eval('varstr(basering)'))

    def sage_poly(self, R=None, kcache=None):
        """
        Returns a Sage polynomial in the ring r matching the provided poly
        which is a singular polynomial.

        INPUT:


        -  ``R`` - (default: None); an optional polynomial ring.
           If it is provided, then you have to make sure that it
           matches the current singular ring as, e.g., returned by
           singular.current_ring(). By default, the output of
           :meth:`sage_global_ring` is used.

        -  ``kcache`` - (default: None); an optional dictionary
           for faster finite field lookups, this is mainly useful for finite
           extension fields


        OUTPUT: MPolynomial

        EXAMPLES::

            sage: R = PolynomialRing(GF(2^8,'a'),2,'xy')
            sage: f=R('a^20*x^2*y+a^10+x')
            sage: f._singular_().sage_poly(R)==f
            True
            sage: R = PolynomialRing(GF(2^8,'a'),1,'x')
            sage: f=R('a^20*x^3+x^2+a^10')
            sage: f._singular_().sage_poly(R)==f
            True

        ::

            sage: P.<x,y> = PolynomialRing(QQ, 2)
            sage: f = x*y**3 - 1/9 * x + 1; f
            x*y^3 - 1/9*x + 1
            sage: singular(f)
            x*y^3-1/9*x+1
            sage: P(singular(f))
            x*y^3 - 1/9*x + 1

        TESTS::

            sage: singular.eval('ring r = (3,z),(a,b,c),dp')
            'ring r = (3,z),(a,b,c),dp;'
            sage: singular.eval('minpoly = 1+z+z2+z3+z4')
            'minpoly = 1+z+z2+z3+z4;'
            sage: p = singular('z^4*a^3+z^2*a*b*c')
            sage: p.sage_poly()
            (-z^3 - z^2 - z - 1)*a^3 + (z^2)*a*b*c
            sage: singular('z^4')
            (-z3-z2-z-1)

        AUTHORS:

        - Martin Albrecht (2006-05-18)
        - Simon King (2011-06-06): Deal with Singular's short polynomial representation,
          automatic construction of a polynomial ring, if it is not explicitly given.

        .. note::

           For very simple polynomials
           ``eval(SingularElement.sage_polystring())`` is faster than
           SingularElement.sage_poly(R), maybe we should detect the
           crossover point (in dependence of the string length) and
           choose an appropriate conversion strategy
        """
        # TODO: Refactor imports to move this to the top
        from sage.rings.polynomial.multi_polynomial_ring import MPolynomialRing_polydict
        from sage.rings.polynomial.multi_polynomial_libsingular import MPolynomialRing_libsingular
        from sage.rings.polynomial.multi_polynomial_element import MPolynomial_polydict
        from sage.rings.polynomial.polynomial_ring import is_PolynomialRing
        from sage.rings.polynomial.polydict import PolyDict,ETuple
        from sage.rings.polynomial.polynomial_singular_interface import can_convert_to_singular
        from sage.rings.quotient_ring import QuotientRing_generic
        from sage.rings.quotient_ring_element import QuotientRingElement

        ring_is_fine = False
        if R is None:
            ring_is_fine = True
            R = self.sage_global_ring()

        sage_repr = {}
        k = R.base_ring()

        variable_str = "*".join(R.variable_names())

        # This returns a string which looks like a list where the first
        # half of the list is filled with monomials occurring in the
        # Singular polynomial and the second half filled with the matching
        # coefficients.
        #
        # Our strategy is to split the monomials at "*" to get the powers
        # in the single variables and then to split the result to get
        # actual exponent.
        #
        # So e.g. ['x^3*y^3','a'] get's split to
        # [[['x','3'],['y','3']],'a']. We may do this quickly,
        # as we know what to expect.

        is_short = self.parent().eval('short')
        if is_short!='0':
            self.parent().eval('short=0')
            if isinstance(R, MPolynomialRing_libsingular):
                out = R(self)
                self.parent().eval('short=%s'%is_short)
                return out
            singular_poly_list = self.parent().eval("string(coef(%s,%s))"%(\
                    self.name(),variable_str)).split(",")
            self.parent().eval('short=%s'%is_short)
        else:
            if isinstance(R, MPolynomialRing_libsingular):
                return R(self)
            singular_poly_list = self.parent().eval("string(coef(%s,%s))"%(\
                    self.name(),variable_str)).split(",")

        if singular_poly_list == ['1','0'] :
            return R(0)

        coeff_start = len(singular_poly_list) // 2

        if isinstance(R,(MPolynomialRing_polydict,QuotientRing_generic)) and (ring_is_fine or can_convert_to_singular(R)):
            # we need to lookup the index of a given variable represented
            # through a string
            var_dict = dict(zip(R.variable_names(),range(R.ngens())))

            ngens = R.ngens()

            for i in range(coeff_start):
                exp = dict()
                monomial = singular_poly_list[i]

                if monomial!="1":
                    variables = [var.split("^") for var in monomial.split("*") ]
                    for e in variables:
                        var = e[0]
                        if len(e)==int(2):
                            power = int(e[1])
                        else:
                            power=1
                        exp[var_dict[var]]=power

                if kcache is None:
                    sage_repr[ETuple(exp,ngens)]=k(singular_poly_list[coeff_start+i])
                else:
                    elem = singular_poly_list[coeff_start+i]
                    if elem not in kcache:
                        kcache[elem] = k( elem )
                    sage_repr[ETuple(exp,ngens)]= kcache[elem]

            p = MPolynomial_polydict(R,PolyDict(sage_repr,force_int_exponents=False,force_etuples=False))
            if isinstance(R, MPolynomialRing_polydict):
                return p
            else:
                return QuotientRingElement(R,p,reduce=False)

        elif is_PolynomialRing(R) and (ring_is_fine or can_convert_to_singular(R)):

            sage_repr = [0]*int(self.deg()+1)

            for i in range(coeff_start):
                monomial = singular_poly_list[i]
                exp = int(0)

                if monomial!="1":
                    term =  monomial.split("^")
                    if len(term)==int(2):
                        exp = int(term[1])
                    else:
                        exp = int(1)

                if kcache is None:
                    sage_repr[exp]=k(singular_poly_list[coeff_start+i])
                else:
                    elem = singular_poly_list[coeff_start+i]
                    if elem not in kcache:
                        kcache[elem] = k( elem )
                    sage_repr[ exp ]= kcache[elem]

            return R(sage_repr)

        else:
            raise TypeError("Cannot coerce %s into %s"%(self,R))

    def sage_matrix(self, R, sparse=True):
        """
        Returns Sage matrix for self

        INPUT:

        -  ``R`` - (default: None); an optional ring, over which
           the resulting matrix is going to be defined.
           By default, the output of :meth:`sage_global_ring` is used.

        - ``sparse`` - (default: True); determines whether the
          resulting matrix is sparse or not.

        EXAMPLES::

            sage: R = singular.ring(0, '(x,y,z)', 'dp')
            sage: A = singular.matrix(2,2)
            sage: A.sage_matrix(ZZ)
            [0 0]
            [0 0]
            sage: A.sage_matrix(RDF)
            [0.0 0.0]
            [0.0 0.0]
        """
        from sage.matrix.constructor import Matrix
        nrows, ncols = int(self.nrows()),int(self.ncols())

        if R is None:
            R = self.sage_global_ring()
            A = Matrix(R, nrows, ncols, sparse=sparse)
            #this is slow
            for x in range(nrows):
                for y in range(ncols):
                    A[x,y]=self[x+1,y+1].sage_poly(R)
            return A

        A = Matrix(R, nrows, ncols, sparse=sparse)
        #this is slow
        for x in range(nrows):
            for y in range(ncols):
                A[x,y]=R(self[x+1,y+1])

        return A

    def _sage_(self, R=None):
        r"""
        Coerces self to Sage.

        EXAMPLES::

            sage: R = singular.ring(0, '(x,y,z)', 'dp')
            sage: A = singular.matrix(2,2)
            sage: A._sage_(ZZ)
            [0 0]
            [0 0]
            sage: A = random_matrix(ZZ,3,3); A
            [ -8   2   0]
            [  0   1  -1]
            [  2   1 -95]
            sage: As = singular(A); As
            -8     2     0
            0     1    -1
            2     1   -95
            sage: As._sage_()
            [ -8   2   0]
            [  0   1  -1]
            [  2   1 -95]

        ::

            sage: singular.eval('ring R = integer, (x,y,z),lp')
            '// ** redefining R **'
            sage: I = singular.ideal(['x^2','y*z','z+x'])
            sage: I.sage()  # indirect doctest
            Ideal (x^2, y*z, x + z) of Multivariate Polynomial Ring in x, y, z over Integer Ring

        ::

            sage: singular('ringlist(basering)').sage()
            [['integer'], ['x', 'y', 'z'], [['lp', (1, 1, 1)], ['C', (0)]], Ideal (0) of Multivariate Polynomial Ring in x, y, z over Integer Ring]

        ::

            sage: singular.eval('ring r10 = (9,a), (x,y,z),lp')
            'ring r10 = (9,a), (x,y,z),lp;'
            sage: singular.eval('setring R')
            'setring R;'
            sage: singular('r10').sage()
            Multivariate Polynomial Ring in x, y, z over Finite Field in a of size 3^2

        Note that the current base ring has not been changed by asking for another ring::

            sage: singular('basering')
            //   coeff. ring is : Integers
            //   number of vars : 3
            //        block   1 : ordering lp
            //                  : names    x y z
            //        block   2 : ordering C

        ::

            sage: singular.eval('setring r10')
            'setring r10;'
            sage: Q = singular('std(ideal(x^2,x+y^2+z^3))', type='qring')
            sage: Q.sage()
            Quotient of Multivariate Polynomial Ring in x, y, z over Finite Field in a of size 3^2 by the ideal (y^4 - y^2*z^3 + z^6, x + y^2 + z^3)
            sage: singular('x^2+y').sage()
            x^2 + y
            sage: singular('x^2+y').sage().parent()
            Quotient of Multivariate Polynomial Ring in x, y, z over Finite Field in a of size 3^2 by the ideal (y^4 - y^2*z^3 + z^6, x + y^2 + z^3)

        """
        typ = self.type()
        if typ=='poly':
            return self.sage_poly(R)
        elif typ == 'module':
            return self.sage_matrix(R,sparse=True)
        elif typ == 'matrix':
            return self.sage_matrix(R,sparse=False)
        elif typ == 'list':
            return [ f._sage_(R) for f in self ]
        elif typ == 'intvec':
            from sage.modules.free_module_element import vector
            return vector([sage.rings.integer.Integer(str(e)) for e in self])
        elif typ == 'intmat':
            from sage.matrix.constructor import matrix
            from sage.rings.integer_ring import ZZ
            A =  matrix(ZZ, int(self.nrows()), int(self.ncols()))
            for i in xrange(A.nrows()):
                for j in xrange(A.ncols()):
                    A[i,j] = sage.rings.integer.Integer(str(self[i+1,j+1]))
            return A
        elif typ == 'string':
            return repr(self)
        elif typ == 'ideal':
            R = R or self.sage_global_ring()
            return R.ideal([p.sage_poly(R) for p in self])
        elif typ in ['ring', 'qring']:
            br = singular('basering')
            self.set_ring()
            R = self.sage_global_ring()
            br.set_ring()
            return R
        raise NotImplementedError("Coercion of this datatype not implemented yet")

    def set_ring(self):
        """
        Sets the current ring in Singular to be self.

        EXAMPLES::

            sage: R = singular.ring(7, '(a,b)', 'ds')
            sage: S = singular.ring('real', '(a,b)', 'lp')
            sage: singular.current_ring()
            //   characteristic : 0 (real)
            //   number of vars : 2
            //        block   1 : ordering lp
            //                  : names    a b
            //        block   2 : ordering C
            sage: R.set_ring()
            sage: singular.current_ring()
            //   characteristic : 7
            //   number of vars : 2
            //        block   1 : ordering ds
            //                  : names    a b
            //        block   2 : ordering C
        """
        self.parent().set_ring(self)


    def sage_flattened_str_list(self):
        """
        EXAMPLES::

            sage: R=singular.ring(0,'(x,y)','dp')
            sage: RL = R.ringlist()
            sage: RL.sage_flattened_str_list()
            ['0', 'x', 'y', 'dp', '1,1', 'C', '0', '_[1]=0']
        """
        s = str(self)
        c = '\[[0-9]*\]:'
        r = re.compile(c)
        s = r.sub('',s).strip()
        return s.split()

    def sage_structured_str_list(self):
        r"""
        If self is a Singular list of lists of Singular elements, returns
        corresponding Sage list of lists of strings.

        EXAMPLES::

            sage: R=singular.ring(0,'(x,y)','dp')
            sage: RL=R.ringlist()
            sage: RL
            [1]:
               0
            [2]:
               [1]:
                  x
               [2]:
                  y
            [3]:
               [1]:
                  [1]:
                     dp
                  [2]:
                     1,1
               [2]:
                  [1]:
                     C
                  [2]:
                     0
            [4]:
               _[1]=0
            sage: RL.sage_structured_str_list()
            ['0', ['x', 'y'], [['dp', '1,\n1 '], ['C', '0 ']], '0']
        """
        if not (self.type()=='list'):
            return str(self)
        return [X.sage_structured_str_list() for X in self]

    def trait_names(self):
        """
        Returns the possible tab-completions for self. In this case, we
        just return all the tab completions for the Singular object.

        EXAMPLES::

            sage: R = singular.ring(0,'(x,y)','dp')
            sage: R.trait_names()
            ['exteriorPower',
             ...
             'stdfglm']
        """
        return self.parent().trait_names()

    def type(self):
        """
        Returns the internal type of this element.

        EXAMPLES::

            sage: R = PolynomialRing(GF(2^8,'a'),2,'x')
            sage: R._singular_().type()
            'ring'
            sage: fs = singular('x0^2','poly')
            sage: fs.type()
            'poly'
        """
        # singular reports // $varname $type $stuff
        p = re.compile("// [\w]+ (\w+) [\w]*")
        m = p.match(self.parent().eval("type(%s)"%self.name()))
        return m.group(1)

    def __iter__(self):
        """
        EXAMPLES::

            sage: R = singular.ring(0, '(x,y,z)', 'dp')
            sage: A = singular.matrix(2,2)
            sage: list(iter(A))
            [[0], [0]]
            sage: A[1,1] = 1; A[1,2] = 2
            sage: A[2,1] = 3; A[2,2] = 4
            sage: list(iter(A))
            [[1,3], [2,4]]
        """
        if self.type()=='matrix':
            l = self.ncols()
        else:
            l = len(self)
        for i in range(1, l+1):
            yield self[i]

    def _singular_(self):
        """
        EXAMPLES::

            sage: R = singular.ring(0, '(x,y,z)', 'dp')
            sage: A = singular.matrix(2,2)
            sage: A._singular_() is A
            True
        """
        return self

    def attrib(self, name, value=None):
        """
        Get and set attributes for self.

        INPUT:


        -  ``name`` - string to choose the attribute

        -  ``value`` - boolean value or None for reading,
           (default:None)


        VALUES: isSB - the standard basis property is set by all commands
        computing a standard basis like groebner, std, stdhilb etc.; used
        by lift, dim, degree, mult, hilb, vdim, kbase isHomog - the weight
        vector for homogeneous or quasihomogeneous ideals/modules isCI -
        complete intersection property isCM - Cohen-Macaulay property rank
        - set the rank of a module (see nrows) withSB - value of type
        ideal, resp. module, is std withHilb - value of type intvec is
        hilb(_,1) (see hilb) withRes - value of type list is a free
        resolution withDim - value of type int is the dimension (see dim)
        withMult - value of type int is the multiplicity (see mult)

        EXAMPLE::

            sage: P.<x,y,z> = PolynomialRing(QQ)
            sage: I = Ideal([z^2, y*z, y^2, x*z, x*y, x^2])
            sage: Ibar = I._singular_()
            sage: Ibar.attrib('isSB')
            0
            sage: singular.eval('vdim(%s)'%Ibar.name()) # sage7 name is random
            // ** sage7 is no standard basis
            4
            sage: Ibar.attrib('isSB',1)
            sage: singular.eval('vdim(%s)'%Ibar.name())
            '4'
        """
        if value is None:
            return int(self.parent().eval('attrib(%s,"%s")'%(self.name(),name)))
        else:
            self.parent().eval('attrib(%s,"%s",%d)'%(self.name(),name,value))

class SingularFunction(ExpectFunction):
    def _sage_doc_(self):
        """
        EXAMPLES::

            sage: 'groebner' in singular.groebner._sage_doc_()
            True
        """
        if not nodes:
            generate_docstring_dictionary()

        prefix = \
"""
This function is an automatically generated pexpect wrapper around the Singular
function '%s'.

EXAMPLE::

    sage: groebner = singular.groebner
    sage: P.<x, y> = PolynomialRing(QQ)
    sage: I = P.ideal(x^2-y, y+x)
    sage: groebner(singular(I))
    x+y,
    y^2-y
"""%(self._name,)
        prefix2 = \
"""

The Singular documentation for '%s' is given below.
"""%(self._name,)

        try:
            return prefix + prefix2 + nodes[node_names[self._name]]
        except KeyError:
            return prefix

class SingularFunctionElement(FunctionElement):
    def _sage_doc_(self):
        r"""
        EXAMPLES::

            sage: R = singular.ring(0, '(x,y,z)', 'dp')
            sage: A = singular.matrix(2,2)
            sage: 'matrix_expression' in A.nrows._sage_doc_()
            True
        """
        if not nodes:
            generate_docstring_dictionary()
        try:
            return nodes[node_names[self._name]]
        except KeyError:
            return ""

def is_SingularElement(x):
    r"""
    Returns True is x is of type ``SingularElement``.

    EXAMPLES::

        sage: from sage.interfaces.singular import is_SingularElement
        sage: is_SingularElement(singular(2))
        True
        sage: is_SingularElement(2)
        False
    """
    return isinstance(x, SingularElement)

def reduce_load():
    """
    Note that this returns an invalid Singular object!

    EXAMPLES::

        sage: from sage.interfaces.singular import reduce_load
        sage: reduce_load()
        (invalid object -- defined in terms of closed session)
    """
    return SingularElement(None, None, None)



nodes = {}
node_names = {}

def generate_docstring_dictionary():
    """
    Generate global dictionaries which hold the docstrings for
    Singular functions.

    EXAMPLE::

        sage: from sage.interfaces.singular import generate_docstring_dictionary
        sage: generate_docstring_dictionary()
    """
    global nodes
    global node_names

    nodes.clear()
    node_names.clear()

    singular_docdir = os.environ["SAGE_LOCAL"]+"/share/singular/"

    new_node = re.compile("File: singular\.hlp,  Node: ([^,]*),.*")
    new_lookup = re.compile("\* ([^:]*):*([^.]*)\..*")

    L, in_node, curr_node = [], False, None

    for line in open(singular_docdir + "singular.hlp"):
        m = re.match(new_node,line)
        if m:
            # a new node starts
            in_node = True
            nodes[curr_node] = "".join(L)
            L = []
            curr_node, = m.groups()
        elif in_node: # we are in a node
           L.append(line)
        else:
           m = re.match(new_lookup, line)
           if m:
               a,b = m.groups()
               node_names[a] = b.strip()

        if line == "6 Index\n":
            in_node = False

    nodes[curr_node] = "".join(L) # last node

def get_docstring(name):
    """
    Return the docstring for the function ``name``.

    INPUT:

    - ``name`` - a Singular function name

    EXAMPLE::

        sage: from sage.interfaces.singular import get_docstring
        sage: 'groebner' in get_docstring('groebner')
        True
        sage: 'standard.lib' in get_docstring('groebner')
        True

    """
    if not nodes:
        generate_docstring_dictionary()
    try:
        return nodes[node_names[name]]
    except KeyError:
        return ""

##################################

singular = Singular()

def reduce_load_Singular():
    """
    EXAMPLES::

        sage: from sage.interfaces.singular import reduce_load_Singular
        sage: reduce_load_Singular()
        Singular
    """
    return singular

import os
def singular_console():
    """
    Spawn a new Singular command-line session.

    EXAMPLES::

        sage: singular_console() #not tested
                             SINGULAR                             /  Development
         A Computer Algebra System for Polynomial Computations   /   version 3-0-4
                                                               0<
             by: G.-M. Greuel, G. Pfister, H. Schoenemann        \   Nov 2007
        FB Mathematik der Universitaet, D-67653 Kaiserslautern    \
    """
    os.system('Singular')


def singular_version():
    """
    Returns the version of Singular being used.

    EXAMPLES:
    """
    return singular.eval('system("--version");')



class SingularGBLogPrettyPrinter:
    """
    A device which prints Singular Groebner basis computation logs
    more verbatim.
    """
    rng_chng = re.compile("\[\d+:\d+\]")# [m:n] internal ring change to
                                        # poly representation with
                                        # exponent bound m and n words in
                                        # exponent vector
    new_elem = re.compile("s")          # found a new element of the standard basis
    red_zero = re.compile("-")          # reduced a pair/S-polynomial to 0
    red_post = re.compile("\.")         # postponed a reduction of a pair/S-polynomial
    cri_hilb = re.compile("h")          # used Hilbert series criterion
    hig_corn = re.compile("H\(\d+\)")   # found a 'highest corner' of degree d, no need to consider higher degrees
    num_crit = re.compile("\(\d+\)")    # n critical pairs are still to be reduced
    red_num =  re.compile("\(S:\d+\)")  # doing complete reduction of n elements
    deg_lead = re.compile("\d+")        # the degree of the leading terms is currently d

    # SlimGB
    red_para = re.compile("M\[(\d+),(\d+)\]") # parallel reduction of n elements with m non-zero output elements
    red_betr = re.compile("b")                # exchange of a reductor by a 'better' one
    non_mini = re.compile("e")                # a new reductor with non-minimal leading term

    crt_lne1 = re.compile("product criterion:(\d+) chain criterion:(\d+)")
    crt_lne2 = re.compile("NF:(\d+) product criterion:(\d+), ext_product criterion:(\d+)")

    pat_sync = re.compile("1\+(\d+);")

    global_pattern = re.compile("(\[\d+:\d+\]|s|-|\.|h|H\(\d+\)|\(\d+\)|\(S:\d+\)|\d+|M\[\d+,[b,e]*\d+\]|b|e).*")

    def __init__(self, verbosity=1):
        """
        Construct a new Singular Groebner Basis log pretty printer.

        INPUT:

        - ``verbosity`` - how much information should be printed
          (between 0 and 3)

        EXAMPLE::

            sage: from sage.interfaces.singular import SingularGBLogPrettyPrinter
            sage: s0 = SingularGBLogPrettyPrinter(verbosity=0)
            sage: s1 = SingularGBLogPrettyPrinter(verbosity=1)
            sage: s0.write("[1:2]12")

            sage: s1.write("[1:2]12")
            Leading term degree: 12.
        """
        self.verbosity = verbosity

        self.curr_deg = 0 # current degree
        self.max_deg = 0  # maximal degree in total

        self.nf = 0 # number of normal forms computed (SlimGB only)
        self.prod = 0 # number of S-polynomials discarded using product criterion
        self.ext_prod = 0 # number of S-polynomials discarded using extended product criterion
        self.chain = 0 # number of S-polynomials discarded using chain criterion

        self.storage = "" # stores incomplete strings
        self.sync = None # should we expect a sync integer?

    def write(self, s):
        """
        EXAMPLE::

            sage: from sage.interfaces.singular import SingularGBLogPrettyPrinter
            sage: s3 = SingularGBLogPrettyPrinter(verbosity=3)
            sage: s3.write("(S:1337)")
            Performing complete reduction of 1337 elements.
            sage: s3.write("M[389,12]")
            Parallel reduction of 389 elements with 12 non-zero output elements.
        """
        verbosity = self.verbosity

        if self.storage:
            s = self.storage + s
            self.storage = ""

        for line in s.splitlines():
            # deal with the Sage <-> Singular syncing code
            match = re.match(SingularGBLogPrettyPrinter.pat_sync,line)
            if match:
                self.sync = int(match.groups()[0])
                continue

            if self.sync and line == "%d"%(self.sync+1):
                self.sync = None
                continue

            if line.endswith(";"):
                continue
            if line.startswith(">"):
                continue

            if line.startswith("std") or line.startswith("slimgb"):
                continue

            # collect stats returned about avoided reductions to zero
            match = re.match(SingularGBLogPrettyPrinter.crt_lne1,line)
            if match:
                self.prod,self.chain = map(int,re.match(SingularGBLogPrettyPrinter.crt_lne1,line).groups())
                self.storage = ""
                continue
            match = re.match(SingularGBLogPrettyPrinter.crt_lne2,line)
            if match:
                self.nf,self.prod,self.ext_prod = map(int,re.match(SingularGBLogPrettyPrinter.crt_lne2,line).groups())
                self.storage = ""
                continue

            while line:
                match = re.match(SingularGBLogPrettyPrinter.global_pattern, line)
                if not match:
                    self.storage = line
                    line = None
                    continue

                token, = match.groups()
                line = line[len(token):]

                if re.match(SingularGBLogPrettyPrinter.rng_chng,token):
                    continue

                elif re.match(SingularGBLogPrettyPrinter.new_elem,token) and verbosity >= 3:
                    print "New element found."

                elif re.match(SingularGBLogPrettyPrinter.red_zero,token) and verbosity >= 2:
                    print "Reduction to zero."

                elif re.match(SingularGBLogPrettyPrinter.red_post, token) and verbosity >= 2:
                    print "Reduction postponed."

                elif re.match(SingularGBLogPrettyPrinter.cri_hilb, token) and verbosity >= 2:
                    print "Hilber series criterion applied."

                elif re.match(SingularGBLogPrettyPrinter.hig_corn, token) and verbosity >= 1:
                    print "Maximal degree found: %s"%token

                elif re.match(SingularGBLogPrettyPrinter.num_crit, token) and verbosity >= 1:
                    print "Leading term degree: %2d. Critical pairs: %s."%(self.curr_deg,token[1:-1])

                elif re.match(SingularGBLogPrettyPrinter.red_num, token) and verbosity >= 3:
                    print "Performing complete reduction of %s elements."%token[3:-1]

                elif re.match(SingularGBLogPrettyPrinter.deg_lead, token):
                    if verbosity >= 1:
                        print "Leading term degree: %2d."%int(token)
                    self.curr_deg = int(token)
                    if self.max_deg < self.curr_deg:
                        self.max_deg = self.curr_deg

                elif re.match(SingularGBLogPrettyPrinter.red_para, token) and verbosity >= 3:
                    m,n = re.match(SingularGBLogPrettyPrinter.red_para,token).groups()
                    print "Parallel reduction of %s elements with %s non-zero output elements."%(m,n)

                elif re.match(SingularGBLogPrettyPrinter.red_betr, token) and verbosity >= 3:
                    print "Replaced reductor by 'better' one."

                elif re.match(SingularGBLogPrettyPrinter.non_mini, token) and verbosity >= 2:
                    print "New reductor with non-minimal leading term found."

    def flush(self):
        """
        EXAMPLE::

            sage: from sage.interfaces.singular import SingularGBLogPrettyPrinter
            sage: s3 = SingularGBLogPrettyPrinter(verbosity=3)
            sage: s3.flush()
        """
        sys.stdout.flush()
