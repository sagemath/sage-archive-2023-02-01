r"""
Interface to Magma

.. note::

   You must have ``magma`` installed on your
   computer for this interface to work. Magma is not free, so it is
   not included with Sage, but you can obtain it from
   http://magma.maths.usyd.edu.au/.

Type ``magma.[tab]`` for a list of all the functions
available from your Magma install. Type
``magma.[tab]?`` for Magma's help about a given
function. Type ``magma(...)`` to create a new Magma
object, and ``magma.eval(...)`` to run a string using
Magma (and get the result back as a string).

Sage provides an interface to the Magma computational algebra
system. This system provides extensive functionality for number
theory, group theory, combinatorics and algebra.

The Magma interface offers three pieces of functionality:


#. ``magma_console()`` - A function that dumps you
   into an interactive command-line Magma session.

#. ``magma(expr)`` - Evaluation of arbitrary Magma
   expressions, with the result returned as a string.

#. ``magma.new(expr)`` - Creation of a Sage object that
   wraps a Magma object. This provides a Pythonic interface to Magma.
   For example, if ``f=magma.new(10)``, then
   ``f.Factors()`` returns the prime factorization of
   `10` computed using Magma.


Parameters
----------

Some Magma functions have optional "parameters", which are
arguments that in Magma go after a colon. In Sage, you pass these
using named function arguments. For example,

::

    sage: E = magma('EllipticCurve([0,1,1,-1,0])')                 # optional - magma
    sage: E.Rank(Bound = 5)                                        # optional - magma
    0

Multiple Return Values
----------------------

Some Magma functions return more than one value. You can control
how many you get using the ``nvals`` named parameter to
a function call::

    sage: n = magma(100)                                           # optional - magma
    sage: n.IsSquare(nvals = 1)                                    # optional - magma
    true
    sage: n.IsSquare(nvals = 2)                                    # optional - magma
    (true, 10)
    sage: n = magma(-2006)                                         # optional - magma
    sage: n.Factorization()                                        # optional - magma
    [ <2, 1>, <17, 1>, <59, 1> ]
    sage: n.Factorization(nvals=2)                                 # optional - magma
    ([ <2, 1>, <17, 1>, <59, 1> ], -1)

We verify that an obviously principal ideal is principal::

    sage: _ = magma.eval('R<x> := PolynomialRing(RationalField())')    # optional - magma
    sage: O = magma.NumberField('x^2+23').MaximalOrder()               # optional - magma
    sage: I = magma('ideal<%s|%s.1>'%(O.name(),O.name()))              # optional - magma
    sage: I.IsPrincipal(nvals=2)                                       # optional - magma
    (true, [1, 0])

Long Input
----------

The Magma interface reads in even very long input (using files) in
a robust manner.

::

    sage: t = '"%s"'%10^10000   # ten thousand character string.       # optional - magma
    sage: a = magma.eval(t)                                            # optional - magma
    sage: a = magma(t)                                                 # optional - magma

Other Examples
--------------

We compute a space of modular forms with character.

::

    sage: N = 20
    sage: D = 20
    sage: eps_top = fundamental_discriminant(D)
    sage: eps = magma.KroneckerCharacter(eps_top, RationalField())        # optional - magma
    sage: M2 = magma.ModularForms(eps)                                    # optional - magma
    sage: print M2                                                        # optional - magma
    Space of modular forms on Gamma_1(5) ...
    sage: print M2.Basis()                                                # optional - magma
    [
    1 + 10*q^2 + 20*q^3 + 20*q^5 + 60*q^7 + ...
    q + q^2 + 2*q^3 + 3*q^4 + 5*q^5 + 2*q^6 + ...
    ]

In Sage/Python (and sort of C++) coercion of an element x into a
structure S is denoted by S(x). This also works for the Magma
interface::

    sage: G = magma.DirichletGroup(20)                                    # optional - magma
    sage: G.AssignNames(['a', 'b'])                                       # optional - magma
    sage: (G.1).Modulus()                                                 # optional - magma
    20
    sage: e = magma.DirichletGroup(40)(G.1)                               # optional - magma
    sage: print e                                                         # optional - magma
    $.1
    sage: print e.Modulus()                                               # optional - magma
    40

We coerce some polynomial rings into Magma::

    sage: R.<y> = PolynomialRing(QQ)
    sage: S = magma(R)                                                    # optional - magma
    sage: print S                                                         # optional - magma
    Univariate Polynomial Ring in y over Rational Field
    sage: S.1                                                             # optional - magma
    y

This example illustrates that Sage doesn't magically extend how
Magma implicit coercion (what there is, at least) works. The errors
below are the result of Magma having a rather limited automatic
coercion system compared to Sage's::

    sage: R.<x> = ZZ[]
    sage: x * 5
    5*x
    sage: x * 1.0
    1.00000000000000*x
    sage: x * (2/3)
    2/3*x
    sage: y = magma(x)                                                    # optional - magma
    sage: y * 5                                                           # optional - magma
    5*x
    sage: y * 1.0                                                         # optional - magma
    Traceback (most recent call last):
    ...
    TypeError: unsupported operand parent(s) for '*': 'Magma' and 'Real Field with 53 bits of precision'
    sage: y * (2/3)                                                       # optional - magma
    Traceback (most recent call last):
    ...
    TypeError: unsupported operand parent(s) for '*': 'Magma' and 'Rational Field'

AUTHORS:

- William Stein (2005): initial version

- William Stein (2006-02-28): added extensive tab completion and
  interactive IPython documentation support.

- William Stein (2006-03-09): added nvals argument for
  magma.functions...
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

import os, re, sys

from sage.structure.parent import Parent
from sage.structure.element import RingElement
from expect import console, Expect, ExpectElement, ExpectFunction, FunctionElement
PROMPT = ">>>"

SAGE_REF = "_sage_ref"
SAGE_REF_RE = re.compile('%s\d+'%SAGE_REF)

import sage.misc.misc
import sage.misc.sage_eval

INTRINSIC_CACHE = '%s/magma_intrinsic_cache.sobj'%sage.misc.misc.DOT_SAGE


EXTCODE_DIR = None
def extcode_dir():
    """
    Return directory that contains all the Magma extcode.  This is put
    in a writable directory owned by the user, since when attached,
    Magma has to write sig and lck files.

    EXAMPLES:
        sage: sage.interfaces.magma.extcode_dir()
        '...dir_.../data/'
    """
    global EXTCODE_DIR
    if not EXTCODE_DIR:
        import shutil
        tmp = sage.misc.misc.tmp_dir()
        shutil.copytree('%s/magma/'%sage.misc.misc.SAGE_EXTCODE, tmp + '/data')
        EXTCODE_DIR = "%s/data/"%tmp
    return EXTCODE_DIR


class Magma(Expect):
    r"""
    Interface to the Magma interpreter.

    Type ``magma.[tab]`` for a list of all the functions
    available from your Magma install. Type
    ``magma.[tab]?`` for Magma's help about a given
    function. Type ``magma(...)`` to create a new Magma
    object, and ``magma.eval(...)`` to run a string using
    Magma (and get the result back as a string).

    .. note::

       If you do not own a local copy of Magma, try using the
       ``magma_free`` command instead, which uses the free demo web
       interface to Magma.

    EXAMPLES:

    You must use nvals = 0 to call a function that doesn't return
    anything, otherwise you'll get an error. (nvals is the number of
    return values.)

    ::

        sage: magma.SetDefaultRealFieldPrecision(200, nvals=0)  # magma >= v2.12; optional - magma
        sage: magma.eval('1.1')   # optional - magma
        '1.1000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000'
        sage: magma.SetDefaultRealFieldPrecision(30, nvals=0)  # optional - magma
    """
    def __init__(self, maxread=10000, script_subdirectory=None,
                 logfile=None, server=None, server_tmpdir=None, user_config=False):
        """
        INPUT:


        -  ``maxread`` - affects buffering

        -  ``script_subdirectory`` - directory where scripts
           are read from

        -  ``logfile`` - output logged to this file

        -  ``server`` - address of remote server

        -  ``user_config`` - if True, then local user
           configuration files will be read by Magma. If False (the default),
           then Magma is started with the -n option which supresses user
           configuration files.


        EXAMPLES::

            sage: Magma(maxread=1000, logfile=tmp_filename())
            Magma
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
                        restart_on_ctrlc = False,
                        logfile = logfile,
                        eval_using_file_cutoff=100)
        # We use "-n" above in the Magma startup command so
        # local user startup configuration is not read.

        self.__seq = 0
        self.__ref = 0
        self.__available_var = []
        self.__cache = {}

    def __reduce__(self):
        """
        Used to pickle a magma interface instance.

        Unpickling results in the default magma interpreter; this is a
        choice, and perhaps not the most logical one! It means that if you
        make two distinct magma interfaces, pickle both, then unpickle
        them, you get back exactly the same one. We illustrate this
        behavior below.

        OUTPUT: function, empty tuple

        EXAMPLES::

            sage: loads(dumps(magma)) is magma
            True

        Unpickling always gives the default global magma interpreter::

            sage: m1 = Magma(); m2 = Magma()
            sage: m1 is m2
            False
            sage: loads(dumps(m1)) is loads(dumps(m2))
            True
            sage: loads(dumps(m1)) is magma
            True
        """
        return reduce_load_Magma, tuple([])

    def _read_in_file_command(self, filename):
        """
        Return the command in Magma that reads in the contents of the given
        file.

        INPUT:


        -  ``filename`` - string


        OUTPUT:


        -  ``string`` - a magma command


        EXAMPLES::

            sage: magma._read_in_file_command('file.m')
            'load "file.m";'
        """
        return 'load "%s";'%filename

    def _post_process_from_file(self, s):
        r"""
        Used internally in the Magma interface to post-process the result
        of evaluating a string using a file. For Magma what this does is
        delete the first output line, since that is a verbose output line
        that Magma displays about loading a file.

        INPUT:


        -  ``s`` - a string


        OUTPUT: a string

        EXAMPLES::

            sage: magma._post_process_from_file("Loading ...\nHello")
            'Hello'
            sage: magma._post_process_from_file("Hello")
            ''
        """
        if not isinstance(s, str):
            raise RuntimeError, "Error evaluating object in %s:\n%s"%(self,s)
        # Chop off the annoying "Loading ... " message that Magma
        # always outputs no matter what.
        i = s.find('\n')
        if i == -1: # special case -- command produced no output, so no \n
            return ''
        return s[i+1:]

    def __getattr__(self, attrname):
        """
        Return a formal wrapper around a Magma function, or raise an
        AttributeError if attrname starts with an underscore.

        INPUT:


        -  ``attrname`` - a string


        OUTPUT: MagmaFunction instance

        EXAMPLES::

            sage: g = magma.__getattr__('EllipticCurve')
            sage: g
            EllipticCurve
            sage: type(g)
            <class 'sage.interfaces.magma.MagmaFunction'>

        In fact, __getattr__ is called implicitly in the following
        case::

            sage: f = magma.EllipticCurve
            sage: type(f)
            <class 'sage.interfaces.magma.MagmaFunction'>
            sage: f
            EllipticCurve
        """
        if attrname[:1] == "_":
            raise AttributeError
        return MagmaFunction(self, attrname)

    def chdir(self, dir):
        """
        Change the Magma interpreters current working directory.

        INPUT:


        -  ``dir`` - a string


        EXAMPLES::

            sage: magma.eval('System("pwd")')   # random, optional - magma
            '/Users/was/s/devel/sage-main/sage'
            sage: magma.chdir('..')             # optional - magma
            sage: magma.eval('System("pwd")')   # random, optional - magma
            '/Users/was/s/devel/sage-main'
        """
        self.eval('ChangeDirectory("%s")'%dir, strip=False)

    def eval(self, x, strip=True, **kwds):
        """
        Evaluate the given block x of code in Magma and return the output
        as a string.

        INPUT:


        -  ``x`` - string of code

        -  ``strip`` - ignored


        OUTPUT: string

        EXAMPLES: We evaluate a string that involves assigning to a
        variable and printing.

        ::

            sage: magma.eval("a := 10;print 2+a;")      # optional - magma
            '12'

        We evaluate a large input line (note that no weird output appears
        and that this works quickly).

        ::

            sage: magma.eval("a := %s;"%(10^10000))    # optional - magma
            "
        """
        x = str(x).rstrip()
        if len(x) == 0 or x[len(x) - 1] != ';':
            x += ';'
        ans = Expect.eval(self, x, **kwds).replace('\\\n','')
        if 'Runtime error' in ans or 'User error' in ans:
            raise RuntimeError, "Error evaluating Magma code.\nIN:%s\nOUT:%s"%(x, ans)
        return ans

    def _start(self):
        """
        Initialize a Magma interface instance. This involves (1) setting up
        an obfuscated prompt, and (2) attaching the MAGMA_SPEC file (see
        EXTCODE_DIR/spec file (see sage.interfaces.magma.EXTCODE_DIR/spec).

        EXAMPLES: This is not too exciting::

            sage: magma._start()          # optional - magma
        """
        self._change_prompt('>')
        Expect._start(self)
        self.eval('SetPrompt("%s"); SetLineEditor(false); SetColumns(0);'%PROMPT)
        self._change_prompt(PROMPT)
        self.expect().expect(PROMPT)
        self.expect().expect(PROMPT)
        self.expect().expect(PROMPT)
        self.attach_spec(extcode_dir() + '/spec')

    def set(self, var, value):
        """
        Set the variable var to the given value in the Magma interpreter.

        INPUT:


        -  ``var`` - string; a variable name

        -  ``value`` - string; what to set var equal to


        EXAMPLES::

            sage: magma.set('abc', '2 + 3/5')       # optional - magma
            sage: magma('abc')                      # optional - magma
            13/5
        """
        out = self.eval("%s:=%s"%(var, value))
        if out.lower().find("error") != -1:
            raise TypeError, "Error executing Magma code:\n%s"%out

    def get(self, var):
        """
        Get the value of the variable var.

        INPUT:


        -  ``var`` - string; name of a variable defined in the
           Magma session


        OUTPUT:


        -  ``string`` - string representation of the value of
           the variable.


        EXAMPLES::

            sage: magma.set('abc', '2 + 3/5')     # optional - magma
            sage: magma.get('abc')                # optional - magma
            '13/5'
        """
        return self.eval("%s"%var)

    def objgens(self, value, gens):
        """
        Create a new object with given value and gens.

        INPUT:


        -  ``value`` - something coercible to an element of
           this Magma interface

        -  ``gens`` - string; command separated list of
           variable names


        OUTPUT: new Magma element that is equal to value with given gens

        EXAMPLES::

            sage: R = magma.objgens('PolynomialRing(Rationals(),2)', 'alpha,beta')    # optional - magma
            sage: R.gens()          # optional - magma
            [alpha, beta]

        Because of how Magma works you can use this to change the variable
        names of the generators of an object::

            sage: S = magma.objgens(R, 'X,Y')          # optional - magma
            sage: R                                    # optional - magma
            Polynomial ring of rank 2 over Rational Field
            Lexicographical Order
            Variables: X, Y
            sage: S                                    # optional - magma
            Polynomial ring of rank 2 over Rational Field
            Lexicographical Order
            Variables: X, Y
        """
        var = self._next_var_name()
        value = self(value)
        out = self.eval("_zsage_<%s> := %s; %s := _zsage_"%(gens, value.name(), var))
        if out.lower().find("error") != -1:
            raise TypeError, "Error executing Magma code:\n%s"%out
        return self(var)

    def __call__(self, x, gens=None):
        """
        Coerce x into this Magma interpreter interface.

        INPUT:


        -  ``x`` - object

        -  ``gens`` - string; names of generators of self,
           separated by commas


        OUTPUT: MagmaElement

        EXAMPLES::

            sage: magma(EllipticCurve('37a'))                   # optional - magma
            Elliptic Curve defined by y^2 + y = x^3 - x over Rational Field
            sage: magma('EllipticCurve([GF(5)|1,2,3,4,1])')     # optional - magma
            Elliptic Curve defined by y^2 + x*y + 3*y = x^3 + 2*x^2 + 4*x + 1 over GF(5)
            sage: magma('PowerSeriesRing(Rationals())', 't')    # optional - magma
            Power series ring in t over Rational Field
            sage: magma('PolynomialRing(RationalField(), 3)', 'x,y,z')  # optional - magma
            Polynomial ring of rank 3 over Rational Field
            Lexicographical Order
            Variables: x, y, z

        We test a coercion between different Magma instances::

            sage: m = Magma()
            sage: n = Magma()
            sage: a = n(m(2))           # optional - magma
            sage: a.parent() is n       # optional - magma
            True
            sage: a.parent() is m       # optional - magma
            False

        We test caching::

            sage: R.<x> =  ZZ[]                     # optional - magma
            sage: magma(R) is magma(R)              # optional - magma
            True
            sage: m = Magma()                       # optional - magma
            sage: m(R)                              # optional - magma
            Univariate Polynomial Ring in x over Integer Ring
            sage: m(R) is magma(R)                  # optional - magma
            False
            sage: R._magma_cache                    # optional - magma
            {Magma: Univariate Polynomial Ring in x over Integer Ring,
             Magma: Univariate Polynomial Ring in x over Integer Ring}
        """
        if isinstance(x, bool):
            return Expect.__call__(self, 'true' if x else 'false')

        if gens is not None:  # get rid of this at some point -- it's weird
            return self.objgens(x, gens)

        # This is mostly about caching the Magma element in the object
        # itself below.  Note that it is *very* important that caching
        # happen on the object itself, and not in a dictionary that is
        # held by the Magma interface, since we want garbage collection
        # of the objects in the Magma interface to work correctly.
        has_cache = hasattr(x, '_magma_cache')
        try:
            if has_cache and x._magma_cache.has_key(self):
                A = x._magma_cache[self]
                if A._session_number == self._session_number:
                    return A
        except AttributeError:
            # This happens when x has _magma_cache as a cdef public object attribute.
            x._magma_cache = {}

        try:
            if self.__cache.has_key(x):
                A = self.__cache[x]
                if A._session_number == self._session_number:
                    return A
        except TypeError:  # if x isn't hashable
            pass

        A = Expect.__call__(self, x)
        if has_cache:
            x._magma_cache[self] = A
        else:
            try:  # use try/except here, because if x is cdef'd we won't be able to set this.
                x._magma_cache = {self:A}
            except AttributeError, msg:
                # Unfortunately, we *have* do have this __cache
                # attribute, which can lead to "leaks" in the working
                # Magma session.  This is because it is critical that
                # parent objects get cached, but sometimes they can't
                # be cached in the object itself, because the object
                # doesn't have a _magma_cache attribute.  So in such
                # cases when the object is a parent we cache it in
                # the magma interface.
                if isinstance(x, Parent):
                    self.__cache[x] = A
        return A


    def _coerce_from_special_method(self, x):
        """
        Tries to coerce to self by calling a special underscore method.

        If no such method is defined, raises an AttributeError instead of a
        TypeError.

        EXAMPLES::

            sage: magma._coerce_from_special_method(-3/5)     # optional - magma
            -3/5

        Note that AttributeError::

            sage: magma._coerce_from_special_method('2 + 3')  # optional - magma
            Traceback (most recent call last):
            ...
            AttributeError: 'str' object has no attribute '_magma_init_'
        """
        s = x._magma_init_(self)
        a = self(s)

        # dereference all _sage_ref's used in this string.
        while True:
            z = SAGE_REF_RE.search(s)
            if not z: break
            self.eval('delete %s;'%s[z.start():z.end()])
            s = s[z.end()+1:]
        return a

    def _with_names(self, s, names):
        """
        Return s but wrapped by a call to SageCreateWithNames. This is just
        a very simple convenience function so that code is cleaner.

        INPUT:


        -  ``s`` - string

        -  ``names`` - list of strings


        OUTPUT: string

        EXAMPLES::

            sage: magma._with_names('PolynomialRing(RationalField())', ['y'])     # optional - magma
            'SageCreateWithNames(PolynomialRing(RationalField()),["y"])'
        """
        return 'SageCreateWithNames(%s,[%s])'%(s, ','.join('"%s"'%x for x in names))

    def clear(self, var):
        """
        Clear the variable named var and make it available to be used
        again.

        INPUT:


        -  ``var`` - a string


        EXAMPLES::

            sage: magma = Magma()      # optional - magma
            sage: magma.clear('foo')   # sets foo to 0 in magma; optional - magma
            sage: magma.eval('foo')    # optional - magma
            '0'

        Because we cleared foo, it is set to be used as a variable name in
        the future::

            sage: a = magma('10')      # optional - magma
            sage: a.name()             # optional - magma
            'foo'

        The following tests that the whole variable clearing and freeing
        system is working correctly.

        ::

            sage: magma = Magma()      # optional - magma
            sage: a = magma('100')     # optional - magma
            sage: a.name()             # optional - magma
            '_sage_[1]'
            sage: del a                # optional - magma
            sage: b = magma('257')     # optional - magma
            sage: b.name()             # optional - magma
            '_sage_[1]'
            sage: del b                # optional - magma
            sage: magma('_sage_[1]')   # optional - magma
            0
        """
        self.__available_var.append(var)
        self.eval("%s:=0"%var)

    def cputime(self, t=None):
        """
        Return the CPU time in seconds that has elapsed since this Magma
        session started. This is a floating point number, computed by
        Magma.

        If t is given, then instead return the floating point time from
        when t seconds had elapsed. This is useful for computing elapsed
        times between two points in a running program.

        INPUT:


        -  ``t`` - float (default: None); if not None, return
           cputime since t


        OUTPUT:


        -  ``float`` - seconds


        EXAMPLES::

            sage: type(magma.cputime())         # optional - magma
            <type 'float'>
            sage: magma.cputime()                # random, optional - magma
            1.9399999999999999
            sage: t = magma.cputime()            # optional - magma
            sage: magma.cputime(t)               # random, optional - magma
            0.02
        """
        if t:
            return float(self.eval('Cputime(%s)'%t))
        else:
            return float(self.eval('Cputime()'))

    def chdir(self, dir):
        """
        Change to the given directory.

        INPUT:


        -  ``dir`` - string; name of a directory


        EXAMPLES::

            sage: magma.chdir('/')                 # optional - magma
            sage: magma.eval('System("pwd")')      # optional - magma
            '/'
        """
        magma.eval('ChangeDirectory("%s")'%dir)

    def attach(self, filename):
        r"""
        Attach the given file to the running instance of Magma.

        Attaching a file in Magma makes all intrinsics defined in the file
        available to the shell. Moreover, if the file doesn't start with
        the ``freeze;`` command, then the file is reloaded
        whenever it is changed. Note that functions and procedures defined
        in the file are *not* available. For only those, use
        ``magma.load(filename)``.

        INPUT:


        -  ``filename`` - a string


        EXAMPLES: Attaching a file that exists is fine::

            sage: magma.attach('%s/data/extcode/magma/sage/basic.m'%Sage_ROOT)    # optional - magma

        Attaching a file that doesn't exist raises an exception::

            sage: magma.attach('%s/data/extcode/magma/sage/basic2.m'%Sage_ROOT)     # optional - magma
            Traceback (most recent call last):
            ...
            RuntimeError: Error evaluating Magma code...
        """
        self.eval('Attach("%s")'%filename)

    Attach = attach

    def attach_spec(self, filename):
        r"""
        Attach the given spec file to the running instance of Magma.

        This can attach numerous other files to the running Magma (see the
        Magma documentation for more details).

        INPUT:


        -  ``filename`` - a string


        EXAMPLES::

            sage: magma.attach_spec('%s/data/extcode/magma/spec'%SAGE_ROOT)    # optional - magma
            sage: magma.attach_spec('%s/data/extcode/magma/spec2'%SAGE_ROOT)   # optional - magma
            Traceback (most recent call last):
            ...
            RuntimeError: Can't open package spec file .../data/extcode/magma/spec2 for reading (No such file or directory)
        """
        s = self.eval('AttachSpec("%s")'%filename)
        if s:
            raise RuntimeError, s.strip()

    AttachSpec = attach_spec

    def load(self, filename):
        r"""
        Load the file with given filename using the 'load' command in the
        Magma shell.

        Loading a file in Magma makes all the functions and procedures in
        the file available. The file should not contain any intrinsics (or
        you'll get errors). It also runs code in the file, which can
        produce output.

        INPUT:


        -  ``filename`` - string


        OUTPUT: output printed when loading the file

        EXAMPLES::

            sage: open(SAGE_TMP + 'a.m','w').write('function f(n) return n^2; end function;\nprint "hi";')
            sage: print magma.load(SAGE_TMP + 'a.m')      # optional - magma
            Loading ".../.sage//temp/.../a.m"
            hi
            sage: magma('f(12)')       # optional - magma
            144
        """
        return self.eval('load "%s"'%filename)

    def _next_var_name(self):
        """
        Return the next available variable name in Magma.

        OUTPUT: string

        EXAMPLES::

            sage: m = Magma()
            sage: m._next_var_name()     # optional - magma
            '_sage_[1]'
            sage: m._next_var_name()     # optional - magma
            '_sage_[2]'
            sage: a = m(3/8); a          # optional - magma
            3/8
            sage: a.name()               # optional - magma
            '_sage_[3]'
            sage: m._next_var_name()     # optional - magma
            '_sage_[4]'
        """
        if self.__seq == 0:
            self.eval('_sage_ := [* *];')
        else:
            try:
                self.eval('Append(~_sage_, 0);')
            except:
                # this exception could happen if the Magma process
                # was interrupted during startup / initialization.
                self.eval('_sage_ := [* 0 : i in [1..%s] *];'%self.__seq)
        if len(self.__available_var) > 0:
            n = self.__available_var[0]
            del self.__available_var[0]
            return n
        else:
            self.__seq += 1
            return '_sage_[%s]'%self.__seq

    def _next_ref_name(self):
        """
        Return the next reference name. This is used internally to deal
        with Magma objects that would be deallocated before they are used
        in constructing another object.

        OUTPUT: string

        EXAMPLES::

            sage: magma._next_ref_name()
            '_sage_ref...'
        """
        self.__ref += 1
        return '%s%s'%(SAGE_REF, self.__ref)

    def function_call(self, function, args=[], params={}, nvals=1):
        """
        Return result of evaluating a Magma function with given input,
        parameters, and asking for nvals as output.

        INPUT:


        -  ``function`` - string, a Magma function name

        -  ``args`` - list of objects coercible into this magma
           interface

        -  ``params`` - Magma parameters, passed in aftera
           colon

        -  ``nvals`` - number of return values from the
           function to ask Magma for


        OUTPUT: MagmaElement or tuple of nvals MagmaElement's

        EXAMPLES::

            sage: magma.function_call('Factorization', 100)    # optional - magma
            [ <2, 2>, <5, 2> ]
            sage: magma.function_call('NextPrime', 100, {'Proof':False})    # optional - magma
            101
            sage: magma.function_call('PolynomialRing', [QQ,2])      # optional - magma
            Polynomial ring of rank 2 over Rational Field
            Lexicographical Order
            Variables: $.1, $.2

        Next, we illustrate multiple return values::

            sage: magma.function_call('IsSquare', 100)         # optional - magma
            true
            sage: magma.function_call('IsSquare', 100, nvals=2)     # optional - magma
            (true, 10)
            sage: magma.function_call('IsSquare', 100, nvals=3)     # optional - magma
            Traceback (most recent call last):
            ...
            RuntimeError: Error evaluating Magma code...
            Runtime error in :=: Expected to assign 3 value(s) but only computed 2 value(s)
        """
        args, params = self._convert_args_kwds(args, params)
        nvals = int(nvals)
        if len(params) == 0:
            par = ''
        else:
            par = ' : ' + ','.join(['%s:=%s'%(a,b.name()) for a,b in params.items()])

        fun = "%s(%s%s)"%(function, ",".join([s.name() for s in args]), par)

        return self._do_call(fun, nvals)

    def _do_call(self, code, nvals):
        """
        Evaluate the given code expression assuming that it outputs nvals
        distinct values. Return the resulting values as a tuple if nvals =
        2.

        INPUT:


        -  ``code`` - a string; code to evaluate

        -  ``nvals`` - an integer; number of return values


        OUTPUT: nvals distinct values

        EXAMPLES::

            sage: magma._do_call('SetVerbose("Groebner",2)', 0)     # optional - magma
            sage: magma._do_call('Factorization(-5)', 1)            # optional - magma
            [ <5, 1> ]

        Here we get two outputs, as a tuple.

        ::

            sage: magma._do_call('Factorization(-5)', 2)            # optional - magma
            ([ <5, 1> ], -1)

        You can also do this::

            sage: F, sign = magma._do_call('Factorization(-5)', 2)  # optional - magma
            sage: F                                                 # optional - magma
            [ <5, 1> ]
            sage: sign                                              # optional - magma
            -1

        An expression that has one value.

        ::

            sage: magma._do_call('3^5', 1)                          # optional - magma
            243
        """
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

        nameleft gens

        returning nvals.

        INPUT:


        -  ``left`` - something coerceable to a magma object

        -  ``name`` - name of the constructor, e.g., sub, quo,
           ideal, etc.

        -  ``gens`` - if a list/tuple, each item is coerced to
           magma; otherwise gens itself is converted to magma

        -  ``nvals`` - positive integer; number of return
           values


        OUTPUT: a single magma object if nvals == 1; otherwise a tuple of
        nvals magma objects.

        EXAMPLES: The bar_call function is used by the sub, quo, and ideal
        methods of Magma elements. Here we illustrate directly using
        bar_call to create quotients::

            sage: V = magma.RModule(ZZ,3)    # optional - magma
            sage: V                          # optional - magma
            RModule(IntegerRing(), 3)
            sage: magma.bar_call(V, 'quo', [[1,2,3]], nvals=1)  # optional - magma
            RModule(IntegerRing(), 2)
            sage: magma.bar_call(V, 'quo', [[1,2,3]], nvals=2)  # optional - magma
            (RModule(IntegerRing(), 2),
             Mapping from: RModule(IntegerRing(), 3) to RModule(IntegerRing(), 2))
            sage: magma.bar_call(V, 'quo', V, nvals=2)          # optional - magma
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
        """
        Return the Python class of elements of the Magma interface.

        OUTPUT: a Python class

        EXAMPLES::

            sage: magma._object_class()
            <class 'sage.interfaces.magma.MagmaElement'>
        """
        return MagmaElement

    # Usually "Sequences" are what you want in Magma, not "lists".
    # It's very painful using the interface without this.
    def _left_list_delim(self):
        """
        Return the left sequence delimiter in Magma. Despite the name in
        this function, this is really the least painful choice.

        EXAMPLES::

            sage: magma._left_list_delim()
            '['
        """
        #return "[*"
        return "["

    def _right_list_delim(self):
        """
        Return the right sequence delimiter in Magma. Despite the name in
        this function, this is really the least painful choice.

        EXAMPLES::

            sage: magma._right_list_delim()
            ']'
        """
        #return "*]"
        return "]"

    def _assign_symbol(self):
        """
        Returns the assignment symbol in Magma.

        EXAMPLES::

            sage: magma._assign_symbol()
            ':='
        """
        return ":="

    def _equality_symbol(self):
        """
        Returns the equality testing logical symbol in Magma.

        EXAMPLES::

            sage: magma._equality_symbol()
            'eq'
        """
        return 'eq'

    def _lessthan_symbol(self):
        """
        Returns the less than testing logical symbol in Magma.

        EXAMPLES::

            sage: magma._lessthan_symbol()
            ' lt '
        """
        return ' lt '

    def _greaterthan_symbol(self):
        """
        Returns the greater than testing logical symbol in Magma.

        EXAMPLES::

            sage: magma._greaterthan_symbol()
            ' gt '
        """
        return ' gt '

    # For efficiency purposes, you should definitely override these
    # in your derived class.
    def _true_symbol(self):
        """
        Returns the string representation of "truth" in Magma.

        EXAMPLES::

            sage: magma._true_symbol()
            'true'
        """
        return 'true'

    def _false_symbol(self):
        """
        Returns the string representation of "false" in Magma.

        EXAMPLES::

            sage: magma._false_symbol()
            'false'
        """
        return 'false'

    def console(self):
        """
        Run a command line Magma session. This session is completely
        separate from this Magma interface.

        EXAMPLES::

            sage: magma.console()             # not tested
            Magma V2.14-9     Sat Oct 11 2008 06:36:41 on one      [Seed = 1157408761]
            Type ? for help.  Type <Ctrl>-D to quit.
            >
            Total time: 2.820 seconds, Total memory usage: 3.95MB
        """
        magma_console()

    def version(self):
        """
        Return the version of Magma that you have in your PATH on your
        computer.

        OUTPUT:


        -  ``numbers`` - 3-tuple: major, minor, etc.

        -  ``string`` - version as a string


        EXAMPLES::

            sage: magma.version()       # random, optional - magma
            ((2, 14, 9), 'V2.14-9')
        """
        return magma_version()

    def help(self, s):
        """
        Return Magma help on string s.

        This returns what typing ?s would return in Magma.

        INPUT:


        -  ``s`` - string


        OUTPUT: string

        EXAMPLES::

            sage: magma.help("NextPrime")       # optional - magma
            ===============================================================================
            PATH: /magma/ring-field-algebra/integer/prime/next-previous/NextPrime
            KIND: Intrinsic
            ===============================================================================
            NextPrime(n) : RngIntElt -> RngIntElt
            NextPrime(n: parameter) : RngIntElt -> RngIntElt
            ...
        """
        print self.eval('? %s'%s)

    def trait_names(self, verbose=True, use_disk_cache=True):
        """
        Return a list of all Magma commands.

        This is used as a hook to enable custom command completion.

        Magma doesn't provide any fast way to make a list of all commands,
        which is why caching is done by default. Note that an adverse
        impact of caching is that *new* commands are not picked up, e.g.,
        user defined variables or functions.

        INPUT:


        -  ``verbose`` - bool (default: True); whether to
           verbosely output status info the first time the command list is
           built

        -  ``use_disk_cache`` - bool (default: True); use
           cached command list, which is saved to disk.


        OUTPUT: list of strings

        EXAMPLES::

            sage: len(magma.trait_names(verbose=False))    # random, optional - magma
            7261
        """
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
                print "\nCreating list of all Magma intrinsics for use in tab completion."
                print "This takes a few minutes the first time, but is saved to the"
                print "file '%s' for future instant use."%INTRINSIC_CACHE
                print "Magma may produce errors during this process, which are safe to ignore."
                print "Delete that file to force recreation of this cache."
                print "Scanning Magma types ..."
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
                except RuntimeError, msg:  # weird internal problems in Magma type system
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
        Return the Magma ideal defined by L.

        INPUT:


        -  ``L`` - a list of elements of a Sage multivariate
           polynomial ring.


        OUTPUT: The magma ideal generated by the elements of L.

        EXAMPLES::

            sage: R.<x,y> = QQ[]
            sage: magma.ideal([x^2, y^3*x])         # optional - magma
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
        Pn = self(P).name()
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
        Magma.

        INPUT:


        -  ``type`` - string (e.g. 'Groebner')

        -  ``level`` - integer = 0


        EXAMPLES::

            sage: magma.set_verbose("Groebner", 2)      # optional - magma
            sage: magma.get_verbose("Groebner")         # optional - magma
            2
        """
        self.SetVerbose(type,level)

    def SetVerbose(self, type, level):
        """
        Set the verbosity level for a given algorithm class etc. in Magma.

        INPUT:


        -  ``type`` - string (e.g. 'Groebner'), see Magma
           documentation

        -  ``level`` - integer = 0


        .. note::

           This method is provided to be consistent with the Magma
           naming convention.

        ::

            sage: magma.SetVerbose("Groebner", 2)      # optional - magma
            sage: magma.GetVerbose("Groebner")         # optional - magma
            2
        """
        if level < 0:
            raise TypeError, "level must be >= 0"
        self.eval('SetVerbose("%s",%d)'%(type,level))

    def get_verbose(self, type):
        """
        Get the verbosity level of a given algorithm class etc. in Magma.

        INPUT:


        -  ``type`` - string (e.g. 'Groebner'), see Magma
           documentation


        EXAMPLES::

            sage: magma.set_verbose("Groebner", 2)        # optional - magma
            sage: magma.get_verbose("Groebner")           # optional - magma
            2
        """
        return self.GetVerbose(type)

    def GetVerbose(self, type):
        """
        Get the verbosity level of a given algorithm class etc. in Magma.

        INPUT:


        -  ``type`` - string (e.g. 'Groebner'), see Magma
           documentation


        .. note::

           This method is provided to be consistent with the Magma
           naming convention.

        EXAMPLES::

            sage: magma.SetVerbose("Groebner", 2)      # optional - magma
            sage: magma.GetVerbose("Groebner")         # optional - magma
            2
        """
        return int(self.eval('GetVerbose("%s")'%type))

class MagmaFunctionElement(FunctionElement):
    def __call__(self, *args, **kwds):
        """
        Return the result of calling this Magma function at given inputs.

        Use the optional nvals keyword argument to specify that there are
        multiple return values.

        EXAMPLES: We create a MagmaFunctionElement::

            sage: n = magma(-15)                        # optional - magma
            sage: f = n.Factorisation                   # optional - magma
            sage: type(f)                               # optional - magma
            <class 'sage.interfaces.magma.MagmaFunctionElement'>
            sage: f()                                   # optional - magma
            [ <3, 1>, <5, 1> ]

        We verify that the nvals argument works.

        ::

            sage: f(nvals=2)                            # optional - magma
            ([ <3, 1>, <5, 1> ], -1)

        This illustrates the more conventional way of calling a method on
        an object. It's equivalent to the above, but done in all in one
        step.

        ::

            sage: n.Factorization(nvals = 2)            # optional - magma
            ([ <3, 1>, <5, 1> ], -1)
        """
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
        """
        Return the docstring for this function of an element.

        OUTPUT: string

        EXAMPLES::

            sage: n = magma(-15)             # optional - magma
            sage: f = n.Factorisation        # optional - magma
            sage: print f._sage_doc_()       # optional - magma
            (<RngIntElt> n) -> RngIntEltFact, RngIntElt, SeqEnum
            ...
            sage: print n.Factorisation._sage_doc_()    # optional - magma
            (<RngIntElt> n) -> RngIntEltFact, RngIntElt, SeqEnum
            ...
        """
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
        """
        Return string representation of this partially evaluated function.

        This is basically the docstring (as returned by _sage_doc_)
        unless self._name is the name of an attribute of the object, in
        which case this returns the value of that attribute.

        EXAMPLES::

            sage: magma(-15).Factorisation           # optional - magma
            Partially evaluated Magma function or intrinsic 'Factorisation'
            ...

        We create a vector space, set its M attribute to a number, then
        display/get the attribute as a string.

        ::

            sage: V = magma('VectorSpace(RationalField(),2)')    # optional - magma
            sage: V.set_magma_attribute('M', 290398)             # optional - magma
            sage: V.M                                            # optional - magma
            290398
            sage: type(V.M)                                      # optional - magma
            <class 'sage.interfaces.magma.MagmaFunctionElement'>
            sage: type(V.M.__repr__())                           # optional - magma
            <type 'str'>

        Displaying a non-attribute function works as above.

        ::

            sage: V.Dimension                                    # optional - magma
            Partially evaluated Magma function or intrinsic 'Dimension'
            ...
        """
        M = self._obj.parent()
        try:
            return M.eval('%s`%s'%(self._obj.name(), self._name))
        except RuntimeError:
            return "Partially evaluated Magma function or intrinsic '%s'\n\nSignature:\n\n%s"%(
                self._name, self._sage_doc_())


class MagmaFunction(ExpectFunction):
    def __call__(self, *args, **kwds):
        """
        Return the result of calling this Magma function at given inputs.

        Use the optional nvals keyword argument to specify that there are
        multiple return values.

        EXAMPLES: We create a MagmaFunction::

            sage: f = magma.Factorisation                   # optional - magma
            sage: type(f)                                   # optional - magma
            <class 'sage.interfaces.magma.MagmaFunction'>
            sage: f(-15)                                    # optional - magma
            [ <3, 1>, <5, 1> ]

        We verify that the nvals argument works.

        ::

            sage: f(-15, nvals=2)                           # optional - magma
            ([ <3, 1>, <5, 1> ], -1)
            sage: f.__call__(-15, nvals=2)                  # optional - magma
            ([ <3, 1>, <5, 1> ], -1)
        """
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
        """
        Return docstring about this function.

        OUTPUT: string

        EXAMPLES::

            sage: f = magma.Factorisation
            sage: type(f)
            <class 'sage.interfaces.magma.MagmaFunction'>
            sage: print f._sage_doc_()                          # optional - magma
            Intrinsic 'Factorisation'
            ...
        """
        M = self._parent
        s = M.eval(self._name)
        s = sage.misc.misc.word_wrap(s, 80)
        return s


def is_MagmaElement(x):
    """
    Return True if x is of type MagmaElement, and False otherwise.

    INPUT:


    -  ``x`` - any object


    OUTPUT: bool

    EXAMPLES::

        sage: from sage.interfaces.magma import is_MagmaElement
        sage: is_MagmaElement(2)
        False
        sage: is_MagmaElement(magma(2))                    # optional - magma
        True
    """
    return isinstance(x, MagmaElement)

class MagmaElement(ExpectElement):
    def _ref(self):
        """
        Return a variable name that is a new reference to this particular
        MagmaElement in Magma. This keeps this object from being garbage
        collected by Magma, even if all the Sage references to it are
        freed.

        Important special behavior: When _ref is used during an implicit
        call to _magma_init_, then the reference disappears after the
        coercion is done. More precisely, if the output of _ref() appears
        as part of the output of a call to _magma_init_ that is then
        going to be input to magma(...), then it is deleted in the Magma
        interface. The main use for this behavior is that in
        _magma_init_ it allows you to get a reference to one object, and
        use it exactly once in constructing a string to evaluate in Magma,
        without having to worry about that object being deallocated. There
        are more sophisticated ways that the same problem (with
        _magma_init_ and references) could have been solved, but this
        solution is much simpler and easier to understand than all others I
        came up with. If this doesn't make sense, read the source code to
        _coerce_from_special_method, which is much shorter than this
        paragraph.

        .. warning::

           Use _ref sparingly, since it involves a full call to Magma,
           which can be slow.

        OUTPUT: string

        EXAMPLES::

            sage: a = magma('-2/3')                          # optional - magma
            sage: s = a._ref(); s                            # optional - magma
            '_sage_ref...'
            sage: magma(s)                                   # optional - magma
            -2/3
        """
        P = self._check_valid()
        n = P._next_ref_name()
        P.set(n, self.name())
        return n

    def __getattr__(self, attrname):
        """
        INPUT:


        -  ``attrname`` - string


        OUTPUT: a Magma function partially evaluated with self as the first
        input.

        .. note::

           If the input ``attrname`` starts with an underscore, an
           AttributeError is raised so that the actual Python _ method/value
           can be accessed.

        EXAMPLES::

            sage: n = magma(-15)                                     # optional - magma
            sage: type(n)                                            # optional - magma
            <class 'sage.interfaces.magma.MagmaElement'>
            sage: f = n.__getattr__('Factorization')                 # optional - magma
            sage: type(f)                                            # optional - magma
            <class 'sage.interfaces.magma.MagmaFunctionElement'>
            sage: f                                                  # optional - magma
            Partially evaluated Magma function or intrinsic 'Factorization'
            ...
        """
        if attrname[:1] == "_":
            raise AttributeError
        return MagmaFunctionElement(self, attrname)

    def _sage_(self):
        """
        Return Sage version of this object. Use self.sage() to get the Sage
        version.

        EXAMPLES: Enumerated Sets::

            sage: a = magma('{1,2/3,-5/9}')       # optional - magma
            sage: a.sage()                        # optional - magma
            {1, -5/9, 2/3}
            sage: a._sage_()                      # optional - magma
            {1, -5/9, 2/3}
            sage: type(a.sage())                  # optional - magma
            <class 'sage.sets.set.Set_object_enumerated'>
            sage: a = magma('{1,2/3,-5/9}'); a    # optional - magma
            { -5/9, 2/3, 1 }
            sage: a.Type()                        # optional - magma
            SetEnum
            sage: b = a.sage(); b             # optional - magma
            {1, -5/9, 2/3}
            sage: type(b)                         # optional - magma
            <class 'sage.sets.set.Set_object_enumerated'>
            sage: c = magma(b); c                 # optional - magma
            { -5/9, 2/3, 1 }
            sage: c.Type()                        # optional - magma
            SetEnum

        Multisets are converted to lists::

            sage: m = magma('{* 1,2,2,2,4^^2,3 *}')    # optional - magma
            sage: z = m.sage(); z                      # optional - magma
            [1, 2, 2, 2, 3, 4, 4]
            sage: type(z)                              # optional - magma
            <type 'list'>

        Matrices::

            sage: a = matrix(ZZ,3,3,[1..9])
            sage: m = magma(a)                        # optional - magma
            sage: b = m.sage(); b                     # optional - magma
            [1 2 3]
            [4 5 6]
            [7 8 9]
            sage: b == a                             # optional - magma
            True

        A nonsquare matrix::

            sage: a = matrix(ZZ,2,3,[1..6])
            sage: m = magma(a)                       # optional - magma
            sage: m.sage()                           # optional - magma
            [1 2 3]
            [4 5 6]
        """
        z, preparse = self.Sage(nvals = 2)
        s = str(z)
        preparse =  str(preparse) == 'true'
        return sage.misc.sage_eval.sage_eval(s, preparse=preparse)

    def AssignNames(self, names):
        """
        EXAMPLES::

            sage: G = magma.DirichletGroup(20)   # optional - magma
            sage: G.AssignNames(['a','b'])       # optional - magma
            sage: G.1                            # optional - magma
            a

        ::

            sage: G.Elements()                   # optional - magma
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

    def gen(self, n):
        """
        Return the n-th generator of this Magma element. Note that
        generators are 1-based in Magma rather than 0 based!

        INPUT:


        -  ``n`` - a *positive* integer


        OUTPUT: MagmaElement

        EXAMPLES::

            sage: k.<a> = GF(9)
            sage: magma(k).gen(1)         # optional -- magma
            a
            sage: R.<s,t,w> = k[]
            sage: m = magma(R)            # optional -- magma
            sage: m.gen(1)                # optional -- magma
            s
            sage: m.gen(2)                # optional -- magma
            t
            sage: m.gen(3)                # optional -- magma
            w
            sage: m.gen(0)                # optional -- magma
            Traceback (most recent call last):
            ...
            IndexError: index must be positive since Magma indexes are 1-based
            sage: m.gen(4)                # optional -- magma
            Traceback (most recent call last):
            ...
            IndexError: list index out of range
        """
        if n <= 0:
            raise IndexError, "index must be positive since Magma indexes are 1-based"
        return self.gens()[n-1]

    def gens(self):
        """
        Return generators for self.

        If self is named X is Magma, this function evaluates X.1, X.2,
        etc., in Magma until an error occurs. It then returns a Sage list
        of the resulting X.i. Note - I don't think there is a Magma command
        that returns the list of valid X.i. There are numerous ad hoc
        functions for various classes but nothing systematic. This function
        gets around that problem. Again, this is something that should
        probably be reported to the Magma group and fixed there.

        AUTHORS:

        - William Stein (2006-07-02)

        EXAMPLES::

            sage: magma("VectorSpace(RationalField(),3)").gens()         # optional - magma
            [(1 0 0), (0 1 0), (0 0 1)]
            sage: magma("AbelianGroup(EllipticCurve([1..5]))").gens()    # optional - magma
            [$.1]
        """
        try:
            return self._magma_gens
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
        self._magma_gens = G
        return G

    def gen_names(self):
        """
        Return list of Magma variable names of the generators of self.

        .. note::

           As illustrated below, these are not the print names of the
           the generators of the Magma object, but special variable
           names in the Magma session that reference the generators.

        EXAMPLES::

            sage: R.<x,zw> = QQ[]
            sage: S = magma(R)               # optional - magma
            sage: S.gen_names()              # optional - magma
            ('_sage_[...]', '_sage_[...]')
            sage: magma(S.gen_names()[1])    # optional - magma
            zw
        """
        try:
            return self.__gen_names
        except AttributeError:
            self.__gen_names = tuple([x.name() for x in self.gens()])
        return self.__gen_names

    def evaluate(self, *args):
        """
        Evalute self at the inputs.

        INPUT:


        -  ``*args`` - import arguments


        OUTPUT: self(\*args)

        EXAMPLES::

            sage: f = magma('Factorization')            # optional - magma
            sage: f.evaluate(15)                        # optional - magma
            [ <3, 1>, <5, 1> ]
            sage: f(15)                                 # optional - magma
            [ <3, 1>, <5, 1> ]
            sage: f = magma('GCD')                      # optional - magma
            sage: f.evaluate(15,20)                     # optional - magma
            5
        """
        P = self._check_valid()
        v = [P(a) for a in args]
        names = ','.join([str(x) for x in v])
        return P('%s(%s)'%(self.name(), names))

    eval = evaluate

    def __call__(self, *args):
        """
        Coerce something into the object (using the Magma ! notation).

        For function calls, use self.eval(...).

        EXAMPLES::

            sage: M = magma.RMatrixSpace(magma.IntegerRing(), 2, 2)  # optional - magma
            sage: A = M([1,2,3,4]); A        # optional - magma
            [1 2]
            [3 4]
            sage: type(A)                    # optional - magma
            <class 'sage.interfaces.magma.MagmaElement'>
            sage: A.Type()                   # optional - magma
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
        """
        Return iterator over this Magma element.

        OUTPUT: generator object

        .. warning::

           Internally this constructs the list of elements in self in
           Magma, which is not a lazy operation. This is because Magma
           doesn't have a notion of lazy iterators, unfortunately.

        EXAMPLES::

            sage: V = magma('VectorSpace(GF(3),2)')             # optional - magma
            sage: V                                             # optional - magma
            Full Vector space of degree 2 over GF(3)
            sage: w = V.__iter__(); w                           # optional - magma
            <generator object at ...>
            sage: w.next()                                      # optional - magma
            (0 0)
            sage: w.next()                                      # optional - magma
            (1 0)
            sage: list(w)                                       # optional - magma
            [(2 0), (0 1), (1 1), (2 1), (0 2), (1 2), (2 2)]
        """
        P = self._check_valid()
        z = P('[_a : _a in %s]'%self.name())
        for i in range(1,len(z)+1):
            yield z[i]

    def __len__(self):
        r"""
        Return cardinality of this Magma element.

        This is the same as ``#self`` in Magma.

        EXAMPLES::

            sage: V = magma('VectorSpace(GF(3),2)')           # optional - magma
            sage: V                                           # optional - magma
            Full Vector space of degree 2 over GF(3)
            sage: len(V)                                      # optional - magma
            9
            sage: V.__len__()                                 # optional - magma
            9
        """
        P = self._check_valid()
        return int(P.eval('#%s'%self.name()))

    def _polynomial_(self, R):
        """
        Try to convert self into a polynomial in the univariate polynomial
        ring R.

        EXAMPLES::

            sage: R.<x> = QQ[]
            sage: f = magma(x^2 + 2/3*x + 5)                 # optional - magma
            sage: f                                          # optional - magma
            x^2 + 2/3*x + 5
            sage: f.Type()                                   # optional - magma
            RngUPolElt
            sage: f._polynomial_(R)                          # optional - magma
            x^2 + 2/3*x + 5
        """
        return R(list(self.Eltseq()))

    def _latex_(self):
        r"""
        Return latex representation of self.

        AUTHORS:

        - Jennifer Balakrishnan

        Types that are nicely latex include:


        -  rationals

        -  matrices

        -  polynomials

        -  binary quadratic forms

        -  elements of quadratic, cyclotomic number fields, and general
           number fields

        -  points

        -  elliptic curves

        -  power series


        IMPLEMENTATION: Calls latex.m, which is in
        Sage_ROOT/data/extcode/magma/latex.m

        EXAMPLES::

            sage: latex(magma('-2/3'))                            # optional - magma
            \frac{-2}{3}
            sage: magma('-2/3')._latex_()                         # optional - magma
            '\\frac{-2}{3}'

        ::

            sage: magma.eval('R<x> := PolynomialRing(RationalField()); f := (x-17/2)^3;')     # optional - magma
            "
            sage: latex(magma('f'))                               # optional - magma
            x^{3}-\frac{51}{2}x^{2}+\frac{867}{4}x-\frac{4913}{8}

        ::

            sage: latex(magma('(MatrixAlgebra(RationalField(),3)![0,2,3,4,5,6,7,8,9])^(-1)'))    # optional - magma
            \left(\begin{array}{ccc}-1&2&-1\\2&-7&4\\-1&\frac{14}{3}&\frac{-8}{3}\end{array}\right)

        ::

            sage: magma.eval('K<a> := CyclotomicField(11)')       # optional - magma
            "
            sage: latex(magma('a^3 + a - 17/3'))                  # optional - magma
            \frac{-17}{3}+\zeta_{11}+\zeta_{11}^{3}

        ::

            sage: latex(magma('EllipticCurve([1,2/3,3/4,4/5,-5/6])'))    # optional - magma
            y^2+xy+\frac{3}{4}y=x^3+\frac{2}{3}x^2+\frac{4}{5}x-\frac{5}{6}

        ::

            sage: _=magma.eval('R<x> := PolynomialRing(RationalField())')    # optional - magma
            sage: _=magma.eval('K<a> := NumberField(x^3+17*x+2)')            # optional - magma
            sage: latex(magma('(1/3)*a^2 - 17/3*a + 2'))                     # optional - magma
            2-\frac{17}{3}a+\frac{1}{3}a^{2}

        Sage auto-detects the greek letters and puts backslashes in::

            sage: _=magma.eval('R<x> := PolynomialRing(RationalField())')    # optional - magma
            sage: _=magma.eval('K<alpha> := NumberField(x^3+17*x+2)')        # optional - magma
            sage: latex(magma('(1/3)*alpha^2 - 17/3*alpha + 2'))             # optional - magma
            2-\frac{17}{3}\alpha+\frac{1}{3}\alpha^{2}

        ::

            sage: _=magma.eval('R<alpha> := PolynomialRing(RationalField())') # optional - magma
            sage: latex(magma('alpha^3-1/7*alpha + 3'))                      # optional - magma
            \alpha^{3}-\frac{1}{7}\alpha+3

        Finite field elements::

            sage: _=magma.eval('K<a> := GF(27)')                             # optional - magma
            sage: latex(magma('a^2+2'))                                      # optional - magma
            2+a^{2}

        Printing of unnamed (dollar sign) generators works correctly::

            sage: latex(magma('FiniteField(81).1^2+1'))                      # optional - magma
            1+\$.1^{2}

        Finite fields::

            sage: latex(magma('FiniteField(3)'))                             # optional - magma
            \mathbf{F}_{{3}}
            sage: latex(magma('FiniteField(27)'))                            # optional - magma
            \mathbf{F}_{{3}^{3}}

        Power Series::

            sage: _=magma.eval('R<x> := PowerSeriesRing(RationalField())')   # optional - magma
            sage: latex(magma('(1/(1+x))'))                                  # optional - magma
            1-x+x^{2}-x^{3}+x^{4}-x^{5}+x^{6}-x^{7}+x^{8}-x^{9}+x^{10}-x^{11}+x^{12}-x^{13}+x^{14}-x^{15}+x^{16}-x^{17}+x^{18}-x^{19}+O(x^{20})
            sage: _=magma.eval('R<x> := PowerSeriesRing(RationalField())')   # optional - magma
            sage: latex(magma('(-1/(2+x + O(x^3)))'))                        # optional - magma
            \frac{-1}{2}+\frac{1}{4}x-\frac{1}{8}x^{2}+O(x^{3})

        p-adic Numbers::

            sage: latex(magma('pAdicField(7,4)!9333294394/49'))              # optional - magma
            4\cdot{}7^{-2} + 5\cdot{}7^{-1} + 5+ 6\cdot{}7^{1} + O(7^{2})
        """
        P = self._check_valid()
        s = str(P.eval('Latex(%s)'%self.name()))
        v = '\\mbox{\\rm '
        if s[:len(v)] == v:
            raise AttributeError
        return s

    def set_magma_attribute(self, attrname, value):
        """
        INPUTS: attrname - string value - something coercible to a
        MagmaElement

        EXAMPLES::

            sage: V = magma("VectorSpace(RationalField(),2)")   # optional - magma
            sage: V.set_magma_attribute('M',10)                 # optional - magma
            sage: V.get_magma_attribute('M')                    # optional - magma
            10
            sage: V.M                                           # optional - magma
            10
        """
        P = self.parent()   # instance of Magma that contains this element.
        if not (isinstance(value, MagmaElement) and value.parent() is P):
            value = P(value)
        P.eval('%s`%s := %s'%(self.name(), attrname, value.name()))

    def get_magma_attribute(self, attrname):
        """
        Return value of a given Magma attribute. This is like selfattrname
        in Magma.

        OUTPUT: MagmaElement

        EXAMPLES::

            sage: V = magma("VectorSpace(RationalField(),10)")   # optional - magma
            sage: V.set_magma_attribute('M','"hello"')           # optional - magma
            sage: V.get_magma_attribute('M')                     # optional - magma
            hello
            sage: V.M                                            # optional - magma
            hello
        """
        P = self.parent()
        return P('%s`%s'%(self.name(), attrname))

    def list_attributes(self):
        """
        Return the attributes of self, obtained by calling the
        ListAttributes function in Magma.

        OUTPUT: list of strings

        EXAMPLES: We observe that vector spaces in Magma have numerous
        funny and mysterious attributes.

        ::

            sage: V = magma("VectorSpace(RationalField(),2)")        # optional - magma
            sage: V.list_attributes()                                # optional - magma
            ['Coroots', 'Roots', 'decomp', 'ssbasis', 'M', 'StrLocalData', 'eisen', 'weights', 'RootDatum', 'T', 'p']
        """
        return magma.eval('ListAttributes(Type(%s))'%\
                          self.name()).split()

    def trait_names(self):
        """
        Return all Magma functions that have this Magma element as first
        input. This is used for tab completion.

        .. note::

           This function can unfortunately be slow if there are a very
           large number of functions, e.g., when self is an
           integer. (This could be fixed by the addition of an
           appropriate function to the Magma kernel, which is
           something that can only be done by the Magma developers.)

        OUTPUT:


        -  ``list`` - sorted list of distinct strings


        EXAMPLES::

            sage: v = magma('2/3').trait_names()       # optional - magma
            sage: type(v[0])                           # optional - magma
            <type 'str'>
        """
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
        Return signatures of all Magma intrinsics that can take self as the
        first argument, as strings.

        INPUT:


        -  ``any`` - (bool: default is False) if True, also
           include signatures with Any as first argument.


        OUTPUT: list of strings

        EXAMPLES::

            sage: v = magma('2/3').methods()          # optional - magma
            sage: v[0]                                # optional - magma
            "'*'..."
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
        Quotient of division of self by other. This is denoted // ("div" in
        magma).

        EXAMPLE::

            sage: R.<x,y,z> = QQ[]
            sage: magma(5)//magma(2)     # optional - magma
            2
            sage: m = magma(x*z + x*y)   # optional - magma
            sage: n = magma(x)           # optional - magma
            sage: m//n                   # optional - magma
            y + z
        """
        return self.parent()('%s div %s'%(self.name(), x.name()))

    def __nonzero__(self):
        """
        Return True if self is nonzero according to Magma. If Magma can't
        decide, i.e., raising an error then also return True.

        EXAMPLES: We define a Magma vector space.

        ::

            sage: V = magma('VectorSpace(GF(3),2)'); V            # optional - magma
            Full Vector space of degree 2 over GF(3)

        The first generator is nonzero.

        ::

            sage: V.gen(1).__nonzero__()                          # optional - magma
            True

        The zero element is zero.

        ::

            sage: V(0).__nonzero__()                              # optional - magma
            False

        The space itself is nonzero (the default - in Magma no comparison
        to 0 is possible).

        ::

            sage: V.__nonzero__()                                 # optional - magma
            True

        We can also call bool which is the opposite of __nonzero__.

        ::

            sage: bool(V(0))                                      # optional - magma
            False
            sage: bool(V.gen(1))                                  # optional - magma
            True
            sage: bool(V)                                         # optional - magma
            True
        """
        try:
            return not self.parent()("%s eq 0"%self.name()).bool()
        except TypeError:
            return True

    def sub(self, gens):
        """
        Return the sub-object of self with given gens.

        INPUT:


        -  ``gens`` - object or list/tuple of generators


        EXAMPLES::

            sage: V = magma('VectorSpace(RationalField(),3)')       # optional - magma
            sage: W = V.sub([ [1,2,3], [1,1,2] ]); W                # optional - magma
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
        Return the quotient of self by the given object or list of
        generators.

        INPUT:


        -  ``gens`` - object or list/tuple of generators


        OUTPUT:


        -  ``magma element`` - the quotient object

        -  ``magma element`` - mapping from self to the
           quotient object


        EXAMPLES::

            sage: V = magma('VectorSpace(RationalField(),3)')       # optional - magma
            sage: V.quo([[1,2,3], [1,1,2]])                         # optional - magma
            (Full Vector space of degree 1 over Rational Field, Mapping from: Full Vector space of degree 3 over Rational Field to Full Vector space of degree 1 over Rational Field)

        We illustrate quotienting out by an object instead of a list of
        generators::

            sage: W = V.sub([ [1,2,3], [1,1,2] ])                   # optional - magma
            sage: V.quo(W)                                          # optional - magma
            (Full Vector space of degree 1 over Rational Field, Mapping from: Full Vector space of degree 3 over Rational Field to Full Vector space of degree 1 over Rational Field)

        We quotient a ZZ module out by a submodule.

        ::

            sage: V = magma.RModule(ZZ,3); V # optional - magma
            RModule(IntegerRing(), 3)
            sage: W, phi = V.quo([[1,2,3]])  # optional - magma
            sage: W                          # optional - magma
            RModule(IntegerRing(), 2)
            sage: phi                        # optional - magma
            Mapping from: RModule(IntegerRing(), 3) to RModule(IntegerRing(), 2)
        """
        return self.parent().bar_call(self, 'quo', gens, nvals=2)

    def ideal(self, gens):
        """
        Return the ideal of self with given list of generators.

        INPUT:


        -  ``gens`` - object or list/tuple of generators


        OUTPUT:


        -  ``magma element`` - a Magma ideal


        EXAMPLES::

            sage: R = magma('PolynomialRing(RationalField())')        # optional - magma
            sage: R.assign_names(['x'])                               # optional - magma
            sage: x = R.1                                             # optional - magma
            sage: R.ideal([x^2 - 1, x^3 - 1])                         # optional - magma
            Ideal of Univariate Polynomial Ring in x over Rational Field generated by x - 1
        """
        return self.parent().bar_call(self, 'ideal', gens, nvals=1)

###########################################################################

magma = Magma()

def reduce_load_Magma():
    """
    Used in unpickling a Magma interface.

    This functions just returns the global default Magma interface.

    EXAMPLES::

        sage: sage.interfaces.magma.reduce_load_Magma()
        Magma
    """
    return magma

def magma_console():
    """
    Run a command line Magma session.

    EXAMPLES::

        sage: magma_console()             # not tested
        Magma V2.14-9     Sat Oct 11 2008 06:36:41 on one      [Seed = 1157408761]
        Type ? for help.  Type <Ctrl>-D to quit.
        >
        Total time: 2.820 seconds, Total memory usage: 3.95MB
    """
    console('magma')

def magma_version():
    """
    Return the version of Magma that you have in your PATH on your
    computer.

    OUTPUT:


    -  ``numbers`` - 3-tuple: major, minor, etc.

    -  ``string`` - version as a string


    EXAMPLES::

        sage: magma_version()       # random, optional - magma
        ((2, 14, 9), 'V2.14-9')
    """
    t = tuple([int(n) for n in magma.eval('GetVersion()').split()])
    return t, 'V%s.%s-%s'%t
