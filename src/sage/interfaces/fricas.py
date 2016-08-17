r"""
Interface to FriCAS

.. TODO::

    - Evaluation using a file is not done. Any input line with more
      than a few thousand characters would hang the system, so currently it
      automatically raises an exception.

    - All completions of a given command.

    - Interactive help.

FriCAS is a free GPL-compatible (modified BSD license) general
purpose computer algebra system based on Axiom.  The FriCAS
website can be found at http://fricas.sourceforge.net/.

AUTHORS:

- Mike Hansen (2009-02): Split off the FriCAS interface from
  the Axiom interface.

If the string "error" (case insensitive) occurs in the output of
anything from FriCAS, a RuntimeError exception is raised.

EXAMPLES: We evaluate a very simple expression in FriCAS.

::

    sage: fricas('3 * 5')                                                       # optional - fricas
    15
    sage: a = fricas(3) * fricas(5); a                                          # optional - fricas
    15

The type of a is FriCASElement, i.e., an element of the FriCAS
interpreter.

::

    sage: type(a)                                                               # optional - fricas
    <class 'sage.interfaces.fricas.FriCASElement'>
    sage: a.parent()                                                            # optional - fricas
    FriCAS

The underlying FriCAS type of a is also available, via the type
method::

    sage: a.typeOf()                                                            # optional - fricas
    PositiveInteger

We factor `x^5 - y^5` in FriCAS in several different ways.
The first way yields a FriCAS object.

::

    sage: F = fricas.factor('x^5 - y^5'); F                                     # optional - fricas
               4      3    2 2    3     4
    - (y - x)(y  + x y  + x y  + x y + x )
    sage: type(F)                                                               # optional - fricas
    <class 'sage.interfaces.fricas.FriCASElement'>
    sage: F.typeOf()                                                            # optional - fricas
    Factored(Polynomial(Integer))

Note that FriCAS objects are normally displayed using "ASCII art".

::

    sage: a = fricas(2/3); a                                                    # optional - fricas
      2
      -
      3
    sage: a = fricas('x^2 + 3/7'); a                                            # optional - fricas
       2   3
      x  + -
           7

The ``fricas.eval`` command evaluates an expression in FriCAS and
returns the result as a string. This is exact as if we typed in the
given line of code to FriCAS; the return value is what FriCAS would
print out, except that type information and line number are stripped
away.

::

    sage: print(fricas.eval('factor(x^5 - y^5)'))                               # optional - fricas
    |startAlgebraOutput|
               4      3    2 2    3     4
    - (y - x)(y  + x y  + x y  + x y + x )
    |endOfAlgebraOutput|

We can create the polynomial `f` as a FriCAS polynomial,
then call the factor method on it. Notice that the notation
``f.factor()`` is consistent with how the rest of Sage
works.

::

    sage: f = fricas('x^5 - y^5')                                               # optional - fricas
    sage: f^2                                                                   # optional - fricas
       10     5 5    10
      y   - 2x y  + x
    sage: f.factor()                                                            # optional - fricas
               4      3    2 2    3     4
    - (y - x)(y  + x y  + x y  + x y + x )

Control-C interruption works well with the FriCAS interface. For
example, try the following sum but with a much bigger range, and hit
control-C.

::

    sage:  f = fricas('(x^5 - y^5)^10000')                                      # not tested - fricas
    Interrupting FriCAS...
    ...
    <type 'exceptions.TypeError'>: Ctrl-c pressed while running FriCAS

::

    sage: fricas('1/100 + 1/101')                                               # optional - fricas
       201
      -----
      10100
    sage: a = fricas('(1 + sqrt(2))^5'); a                                      # optional - fricas
         +-+
      29\|2  + 41

Let us demonstrate some features of FriCAS.  FriCAS can guess a
differential equation for the generating function for integer
partitions::

    sage: fricas("guessADE([partition n for n in 0..40], homogeneous==4)")      # optional - fricas
     [
       [
           n
         [x ]f(x):
              2    3 (iv)          2    2 ,             3  ,,,         2    2 ,,   2
             x f(x) f    (x) + (20x f(x) f (x) + 5x f(x) )f   (x) - 39x f(x) f  (x)
    <BLANKLINE>
           +
                 2     ,   2           2 ,           3  ,,        2 ,   4
             (12x f(x)f (x)  - 15x f(x) f (x) + 4f(x) )f  (x) + 6x f (x)
    <BLANKLINE>
           +
                      ,   3         2 ,   2
             10x f(x)f (x)  - 16f(x) f (x)
    <BLANKLINE>
             =
             0
         ,
                        2     3      4
        f(x)= 1 + x + 2x  + 3x  + O(x )]
       ]

FriCAS can solve linear ordinary differential equations::

    sage: fricas.set("y", "operator y")                                         # optional - fricas
    sage: fricas.set("deq", "x^3*D(y x, x, 3) + x^2*D(y x, x, 2) - 2*x*D(y x, x) + 2*y x - 2*x^4")  # optional - fricas
    sage: fricas.set("sol", "solve(deq, y, x)"); fricas("sol")                  # optional - fricas
                  5      3      2               3     2      3      3     2
                 x  - 10x  + 20x  + 4         2x  - 3x  + 1 x  - 1 x  - 3x  - 1
    [particular= --------------------,basis= [-------------,------,------------]]
                          15x                       x          x         x

    sage: fricas("sol.particular").sage()                                       # optional - fricas
    1/15*(x^5 - 10*x^3 + 20*x^2 + 4)/x
    sage: fricas("sol.basis").sage()                                            # optional - fricas
    [(2*x^3 - 3*x^2 + 1)/x, (x^3 - 1)/x, (x^3 - 3*x^2 - 1)/x]
    sage: fricas.eval(")clear values y deq sol")                                # optional - fricas
    ''

FriCAS can expand expressions into series::

    sage: a = fricas("series(sqrt(cos(x)),x=0)"); a                             # optional - fricas
        1  2    1  4    19   6     559   8     29161    10      11
    1 - - x  - -- x  - ---- x  - ------ x  - --------- x   + O(x  )
        4      96      5760      645120      116121600

    sage: a.coefficients()[38].sage()                                           # optional - fricas
    -29472026335337227150423659490832640468979/274214482066329363682430667508979749984665600000000

    sage: a = fricas("series(sqrt(cos(x)/x),x=0)"); a                           # optional - fricas
       1      3       7
     - -      -       -
       2   1  2    1  2      5
    x    - - x  - -- x  + O(x )
           4      96

    sage: a.coefficient(3/2).sage()                                             # optional - fricas
    -1/4

FriCAS does some limits right:

    sage: fricas("limit(x^2*exp(-x)*Ei(x) - x, x=%plusInfinity)")               # optional - fricas
    1

"""

###########################################################################
#       Copyright (C) 2008 Mike Hansen <mhansen@gmail.com>
#                     2007 Bill Page
#                     2006 William Stein <wstein@gmail.com>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#  The full text of the GPL is available at:
#
#                  http://www.gnu.org/licenses/
###########################################################################
from __future__ import print_function
# from __future__ import absolute_import

from sage.interfaces.tab_completion import ExtraTabCompletion
from sage.interfaces.expect import Expect, ExpectElement, FunctionElement, ExpectFunction
from sage.env import DOT_SAGE
import re
import six

FRICAS_SINGLE_LINE_START = 3 # where the output starts when it fits next to the line number
FRICAS_MULTI_LINE_START = 2  # and when it doesn't
FRICAS_LINE_LENGTH = 80      # length of a line, should match the line length in sage
FRICAS_LINE_BREAK = "\r\n"   # line ending character
FRICAS_WHAT_OPERATIONS_STRING = "Operations whose names satisfy the above pattern\(s\):"
# only the last command should be necessary to make the interface
# work, the other are optimizations.  Beware that lisp distinguishes
# between ' and ".
FRICAS_INIT_CODE = (
")set functions compile on",
")set message autoload off",
")set message type off",
")set output length " + str(FRICAS_LINE_LENGTH),
")lisp (setf |$ioHook|"
"            (lambda (x &optional args)"
"              (when (member x '(|startAlgebraOutput| |endOfAlgebraOutput|"
"                                |startKeyedMsg|      |endOfKeyedMsg|))"
"               (prin1 x)"
"               (princ #\\Newline))))")

FRICAS_LINENUMBER_OFF_CODE = ")lisp (setf |$IOindex| NIL)"
FRICAS_FIRST_PROMPT = "\(1\) -> "
FRICAS_LINENUMBER_OFF_PROMPT = "\(NIL\) -> "

class FriCAS(ExtraTabCompletion, Expect):
>>>>>>> 5b7e6ab... initial version of new fricas interface
    """
    Interface to a FriCAS interpreter.
    """
    def __init__(self, name='fricas', command='fricas -nox -noclef',
                 script_subdirectory=None, logfile=None,
                 server=None, server_tmpdir=None):
        """
        Create an instance of the FriCAS interpreter.

        TESTS::

            sage: fricas == loads(dumps(fricas))
            True
        """
        eval_using_file_cutoff = 600
        assert max(len(c) for c in FRICAS_INIT_CODE) < eval_using_file_cutoff
        self.__eval_using_file_cutoff = eval_using_file_cutoff
        self._COMMANDS_CACHE = '%s/%s_commandlist_cache.sobj'%(DOT_SAGE, name)
        # we run the init code in _start to avoid spurious output
        Expect.__init__(self,
                        name = name,
                        prompt = FRICAS_FIRST_PROMPT,
                        command = command,
                        script_subdirectory = script_subdirectory,
                        server=server,
                        server_tmpdir=server_tmpdir,
                        restart_on_ctrlc = False,
                        verbose_start = False,
                        init_code = [],
                        logfile = logfile,
                        eval_using_file_cutoff=eval_using_file_cutoff)

    def _start(self):
        """
        Start the FriCAS interpreter and switch off the linenumbers.

        EXAMPLES::

            sage: a = FriCAS()                                                  # optional - fricas
            sage: a.is_running()                                                # optional - fricas
            False
            sage: a._start()                                                    # optional - fricas
            sage: a.is_running()                                                # optional - fricas
            True
            sage: a.quit()                                                      # optional - fricas
        """
        # setting the prompt properly is necessary for restarting FriCAS
        self._prompt = FRICAS_FIRST_PROMPT
        Expect._start(self)
        for line in FRICAS_INIT_CODE:
            self.eval(line)
        # switching off the line numbers also modified the prompt
        self._prompt = FRICAS_LINENUMBER_OFF_PROMPT
        self.eval(FRICAS_LINENUMBER_OFF_CODE)

    def _quit_string(self):
        """
        Returns the string used to quit FriCAS.

        EXAMPLES::

            sage: fricas._quit_string()                                         # optional - fricas
            ')quit'
            sage: a = FriCAS()                                                  # optional - fricas
            sage: a.is_running()                                                # optional - fricas
            False
            sage: a._start()                                                    # optional - fricas
            sage: a.is_running()                                                # optional - fricas
            True
            sage: a.quit()                                                      # optional - fricas
            sage: a.is_running()                                                # optional - fricas
            False
        """
        return ')quit'

    def _commands(self):
        """
        Returns a list of commands available. This is done by parsing the
        result of the first section of the output of ')what things'.

        EXAMPLES::

            sage: cmds = fricas._commands()                                     # optional - fricas
            sage: len(cmds) > 100                                               # optional - fricas
            True
            sage: '<' in cmds                                                   # optional - fricas
            True
            sage: 'factor' in cmds                                              # optional - fricas
            True
        """
        output = self.eval(")what operations")
        m = re.search(FRICAS_WHAT_OPERATIONS_STRING + "\r\n(.*)\r\n\|startKeyedMsg\|", output, flags = re.DOTALL)
        l = m.groups()[0].split()
        return l

    def _tab_completion(self, verbose=True, use_disk_cache=True):
        """
        Returns a list of all the commands defined in Fricas and optionally
        (per default) store them to disk.

        EXAMPLES::

            sage: c = fricas._tab_completion(use_disk_cache=False, verbose=False)         # optional - fricas
            sage: len(c) > 100                                                  # optional - fricas
            True
            sage: 'factor' in c                                                 # optional - fricas
            True
            sage: '**' in c                                                     # optional - fricas
            False
            sage: 'upperCase?' in c                                             # optional - fricas
            False
            sage: 'upperCase_q' in c                                            # optional - fricas
            True
            sage: 'upperCase_e' in c                                            # optional - fricas
            True
        """
        try:
            return self.__tab_completion
        except AttributeError:
            import sage.misc.persist
            if use_disk_cache:
                try:
                    self.__tab_completion = sage.misc.persist.load(self._COMMANDS_CACHE)
                    return self.__tab_completion
                except IOError:
                    pass
            if verbose:
                print("\nBuilding %s command completion list (this takes" % self)
                print("a few seconds only the first time you do it).")
                print("To force rebuild later, delete %s." % self._COMMANDS_CACHE)
            v = self._commands()

            #Process we now need process the commands to strip out things which
            #are not valid Python identifiers.
            valid = re.compile('[^a-zA-Z0-9_]+')
            names = [x for x in v if valid.search(x) is None]

            #Change everything that ends with ? to _q and
            #everything that ends with ! to _e
            names += [x[:-1]+"_q" for x in v if x.endswith("?")]
            names += [x[:-1]+"_e" for x in v if x.endswith("!")]

            self.__tab_completion = names
            if len(v) > 200:
                # Fricas is actually installed.
                sage.misc.persist.save(v, self._COMMANDS_CACHE)
            return names

# what I expect from FriCAS:

# 1.) in set(self, var, value)
#
# no markers:
# there could be some "debugging" output, as in fricas("guessADE([1,1,1,1], debug==true)")
#
# startKeyedMsg: an error happened
#
# 2.) in get(self, var)
#
# |startAlgebraOutput\|...|endOfAlgebraOutput\|
#
# 3.) I also need a routine to send a system command and get its output.

    def set(self, var, value):
        """Set the variable var to the given value.

        INPUT:

        - var, value: strings, the first representing a valid FriCAS
          variable identifier, the second a FriCAS expression.

        OUTPUT: None

        EXAMPLES::

            sage: fricas.set('xx', '2')                                         # optional - fricas
            sage: fricas.get('xx')                                              # optional - fricas
            '2'

        """
        # print("fricas.set %s := %s" %(var, value))
        cmd = '%s%s%s;'%(var,self._assign_symbol(), value)
        output = self.eval(cmd)
        m = re.search("\r\n\|startKeyedMsg\|\r\n(.*)\r\n\|endOfKeyedMsg\|\r", output, flags = re.DOTALL)
        if m:
            raise RuntimeError("An error occurred when FriCAS evaluated '%s': %s" % (value, output))

    def get(self, var):
        r"""
        Get the string representation of the value of a variable.

        EXAMPLES::

            sage: fricas.set('xx', '2')                                         # optional - fricas
            sage: fricas.get('xx')                                              # optional - fricas
            '2'
            sage: a = fricas('(1 + sqrt(2))^5')                                 # optional - fricas
            sage: fricas.get(a.name())                                          # optional - fricas
            '   +-+\r\n29\\|2  + 41'
        """
        # print("fricas.get %s" %var)
        output = self.eval(str(var))
        m = re.search("\r\n\|startAlgebraOutput\|\r\n(.*)\r\n\|endOfAlgebraOutput\|\r", output, flags = re.DOTALL)
        if m:
            lines = m.groups()[0].split(FRICAS_LINE_BREAK)
            if max(len(line) for line in lines) < FRICAS_LINE_LENGTH:
                return "\r\n".join(line[FRICAS_SINGLE_LINE_START:] for line in lines)
            else:
                return "\r\n".join(line[FRICAS_MULTI_LINE_START:] for line in lines)

    def _assign_symbol(self):
        """Returns the symbol used for setting a variable.

        This would be used if `self.set` from Interface would not be
        overloaded.
        """
        return ":="

    def _equality_symbol(self):
        """Returns the equality testing logical symbol in FriCAS.

        EXAMPLES::

            sage: a = fricas(x==6); a                                           # optional - fricas
            x= 6
        """
        return "="

    def __repr__(self):
        """
        EXAMPLES::

            sage: fricas                                                        # optional - fricas
            FriCAS
        """
        return "FriCAS"

    def __reduce__(self):
        """
        EXAMPLES::

            sage: fricas.__reduce__()                                           # optional - fricas
            (<function reduce_load_fricas at 0x...>, ())
            sage: f, args = _                                                   # optional - fricas
            sage: f(*args)                                                      # optional - fricas
            FriCAS
        """
        return reduce_load_fricas, tuple([])

    def _function_class(self):
        """
        Return the FriCASExpectFunction class.

        EXAMPLES::

            sage: fricas._function_class()                                      # optional - fricas
            <class 'sage.interfaces.fricas.FriCASExpectFunction'>
            sage: type(fricas.gcd)                                              # optional - fricas
            <class 'sage.interfaces.fricas.FriCASExpectFunction'>
        """
        return FriCASExpectFunction

    def _object_class(self):
        """
        EXAMPLES::

            sage: fricas._object_class()                                        # optional - fricas
            <class 'sage.interfaces.fricas.FriCASElement'>
            sage: type(fricas(2))                                               # optional - fricas
            <class 'sage.interfaces.fricas.FriCASElement'>
        """
        return FriCASElement

    def _function_element_class(self):
        """
        Returns the FriCAS function element class.

        EXAMPLES::

            sage: fricas._function_element_class()                              # optional - fricas
            <class 'sage.interfaces.fricas.FriCASFunctionElement'>
            sage: type(fricas(2).gcd) #optional - fricas                        # optional - fricas
            <class 'sage.interfaces.fricas.FriCASFunctionElement'>
        """
        return FriCASFunctionElement

    def console(self):
        """
        Spawn a new FriCAS command-line session.

        EXAMPLES::

            sage: fricas.console()                                              # not tested
                             FriCAS (AXIOM fork) Computer Algebra System
                                    Version: FriCAS 1.0.5
                     Timestamp: Thursday February 19, 2009 at 06:57:33
            -----------------------------------------------------------------------------
               Issue )copyright to view copyright notices.
               Issue )summary for a summary of useful system commands.
               Issue )quit to leave AXIOM and return to shell.
            -----------------------------------------------------------------------------
        """
        fricas_console()

class FriCASElement(ExpectElement):

    def __len__(self):
        """
        Return the length of a list.

        EXAMPLES::

            sage: v = fricas('[x^i for i in 0..5]')                             # optional - fricas
            sage: len(v)                                                        # optional - fricas
            6
        """
        P = self._check_valid()
        l = P('# %s ' %self._name)
        return l.sage()

    def __getitem__(self, n):
        """We implement the sage conventions here, translating to 0-based iterables.

        We do not check validity, since many objects in FriCAS are
        iterable, in particular Streams

        TODO: we should also implement negative arguments and tuples

        """
        n = int(n)
        if n < 0:
            raise IndexError("index out of range")
        return self.elt(n+1)

    def __int__(self):
        return int(self.sage())

    def bool(self):
        raise NotImplementedError

    def __nonzero__(self):
        raise NotImplementedError

    def __long__(self):
        raise NotImplementedError

    def __float__(self):
        """
        TEST::
        sage: float(fricas(2))                                                  # optional - fricas
        2.0
        """
        return float(self.sage())

    def _integer_(self, ZZ=None):
        raise NotImplementedError

    def _rational_(self):
        raise NotImplementedError

    def gen(self, n):
        raise NotImplementedError

    def _unparsed_InputForm(self):
        """Return the output from FriCAS as a string without the quotes.

        TODO:

        - catch errors, especially when InputForm is not available:

        -- for example when integration returns "failed"
        -- UnivariatePolynomial

        should we provide workarounds, too?

        TESTS:

        We test that strings are returned properly:

        sage: r = fricas('concat([concat(string(i)," ") for i in 0..299])')._unparsed_InputForm()   # optional - fricas
        sage: r == " ".join([str(i) for i in range(300)]) + ' '                                     # optional - fricas
        True

        sage: fricas('concat([string(1) for i in 1..5])')._unparsed_InputForm() == "1"*5            # optional - fricas
        True

        sage: fricas('concat([string(1) for i in 1..10000])')._unparsed_InputForm() == "1"*10000    # optional - fricas
        True
        """
        P = self._check_valid()
        return P.get('unparse(%s::InputForm)' %self._name).replace("\r\n", "")[1:-1]


    def _get_sage_type(self, domain):
        """
        INPUT:

        - domain, a FriCAS SExpression

        OUTPUT:

        - a corresponding Sage type
        """
        from sage.rings.all import ZZ, QQ, QQbar, PolynomialRing, RDF
        from sage.rings.fraction_field import FractionField
        from sage.rings.finite_rings.integer_mod_ring import Integers
        from sage.rings.real_mpfr import RealField
        from sage.symbolic.ring import SR
        from sage.matrix.constructor import matrix

        # first implement domains without arguments
        head = str(domain.car())
        if head in ["Integer", "NonNegativeInteger", "PositiveInteger"]:
            return ZZ
        if head == "String":
            return str
        if head == "Float":
            prec = max(self.mantissa().length()._sage_(), 53)
            return RealField(prec)
        if head == "DoubleFloat":
            return RDF
        if head == "AlgebraicNumber":
            return QQbar

        # now implement "functorial" types
        if head == "OrderedCompletion":
            # this is a workaround, I don't know how translate this
            return SR

        if head == "IntegerMod":
            return Integers(domain[1].integer().sage())

        if head == "Fraction":
            return FractionField(self._get_sage_type(domain[1]))

        if head == "Expression":
            return SR

        if head == "Polynomial":
            # this is a workaround, since in sage we always have to specify the variables
            return SR

        raise NotImplementedError

    def _sage_(self):
        """
        Convert self to a Sage object.

        TESTS:

        Floats:

        sage: fricas(2.1234).sage()                                             # optional - fricas
        2.12340000000000
        sage: _.parent()                                                        # optional - fricas
        Real Field with 53 bits of precision
        sage: a = RealField(100)(pi)
        sage: fricas(a).sage()                                                  # optional - fricas
        3.1415926535897932384626433833
        sage: _.parent()                                                        # optional - fricas
        Real Field with 100 bits of precision
        sage: fricas(a).sage() == a                                             # optional - fricas
        True
        sage: fricas(2.0).sage()                                                # optional - fricas
        2.00000000000000
        sage: _.parent()                                                        # optional  - fricas
        Real Field with 53 bits of precision

        Algebraic numbers:

        sage: a = fricas('(1 + sqrt(2))^5'); a                                  # optional - fricas
           +-+
        29\|2  + 41
        sage: b = a.sage(); b                                                   # optional - fricas
        82.0121933088198?
        sage: b.radical_expression()                                            # optional - fricas
        29*sqrt(2) + 41

        Integers modulo n:

        sage: fricas("((42^17)^1783)::IntegerMod(5^(5^5))").sage() == Integers(5^(5^5))((42^17)^1783) # optional - fricas
        True

        We can also convert FriCAS's polynomials to Sage polynomials.

        sage: a = fricas(x^2 + 1); a.typeOf()                                   # optional - fricas
        Polynomial(Integer)
        sage: a.sage()                                                          # optional - fricas
        x^2 + 1
        sage: _.parent()                                                        # optional - fricas
        Univariate Polynomial Ring in x over Integer Ring
        sage: fricas('x^2 + y^2 + 1/2').sage()                                  # optional - fricas
        y^2 + x^2 + 1/2
        sage: _.parent()                                                        # optional - fricas
        Multivariate Polynomial Ring in y, x over Rational Field

        sage: fricas("1$Polynomial Integer").sage()                             # optional - fricas
        1

        sage: fricas("x^2/2").sage()                                            # optional - fricas
        1/2*x^2

        Rational functions:

        sage: fricas("x^2 + 1/z").sage()                                        # optional - fricas
        x^2 + 1/z

        Expressions:

        sage: fricas("sin(x+y)/exp(z)*log(1+%e)").sage()                        # optional - fricas
        e^(-z)*log(e + 1)*sin(x + y)

        sage: fricas("factorial(n)").sage()                                     # optional - fricas
        factorial(n)

        sage: fricas("integrate(sin(x+y), x=0..1)").sage()                      # optional - fricas
        -cos(y + 1) + cos(y)

        sage: fricas("integrate(x*sin(1/x), x=0..1)").sage()                    # optional - fricas
        'failed'

        Matrices:

        sage: fricas("matrix [[x^n/2^m for n in 0..5] for m in 0..3]").sage()   # optional - fricas
        [      1       x     x^2     x^3     x^4     x^5]
        [    1/2   1/2*x 1/2*x^2 1/2*x^3 1/2*x^4 1/2*x^5]
        [    1/4   1/4*x 1/4*x^2 1/4*x^3 1/4*x^4 1/4*x^5]
        [    1/8   1/8*x 1/8*x^2 1/8*x^3 1/8*x^4 1/8*x^5]

        Lists:

        sage: fricas("[2^n/x^n for n in 0..5]").sage()                          # optional - fricas
        [1, 2/x, 4/x^2, 8/x^3, 16/x^4, 32/x^5]

        sage: fricas("[matrix [[i for i in 1..n]] for n in 0..5]").sage()       # optional - fricas
        [[], [1], [1 2], [1 2 3], [1 2 3 4], [1 2 3 4 5]]

        """
        from sage.rings.all import ZZ, QQ, QQbar, PolynomialRing, RDF
        from sage.rings.fraction_field import FractionField
        from sage.rings.finite_rings.integer_mod_ring import Integers
        from sage.rings.real_mpfr import RealField
        from sage.symbolic.ring import SR
        from sage.matrix.constructor import matrix
        from sage.misc.sage_eval import sage_eval

        # TODO: perhaps we should translate the type first?
        # TODO: perhaps we should get the InputForm as SExpression?

        # remember: fricas.get gives a str and should only used if
        # you know what you are doing, whereas fricas.new gives a
        # FriCASElement
        
        # the coercion to Any gets rid of the Union domain
        P = self._check_valid()
        domain = P.new("dom(%s::Any)" % self._name) # domain is now a fricas SExpression

        # first translate dummy domains such as "failed"
        # we must not recurse here!
        if P.get("string?(%s)" % domain._name) == "true":
            return P.get("string(%s)" % domain._name)[1:-1]

        # now translate domains without arguments
        head = str(domain.car())
        if head == "Boolean":
            return self._unparsed_InputForm() == "true"

        if head in ["Integer", "NonNegativeInteger", "PositiveInteger"]:
            return ZZ(self._unparsed_InputForm())

        if head == "String":
            return self._unparsed_InputForm()

        if head == "Float":
            # Warning: precision$Float gives the current precision,
            # whereas length(mantissa(self)) gives the precision of
            # self.
            prec = max(self.mantissa().length().sage(), 53)
            R = RealField(prec)
            x, e, b = self._unparsed_InputForm().lstrip('float(').rstrip(')').split(',')
            return R(ZZ(x)*ZZ(b)**ZZ(e))

        if head == "DoubleFloat":
            return RDF(self._unparsed_InputForm())

        if head == "AlgebraicNumber":
            s = self._unparsed_InputForm()[:-len("::AlgebraicNumber()")]
            return sage_eval("QQbar(" + s + ")")

        # finally implement "functorial" types
        if head == "List":
            n = int(P.get('# %s'%self._name))
            return [P.new('%s(%s)'%(self._name, k)).sage() for k in range(1, n+1)]

        if head == "Matrix":
            base_ring = self._get_sage_type(domain[1])
            rows = P.new('listOfLists(%s)' %self._name).sage()
            return matrix(base_ring, rows)

        if head == "IntegerMod":
            # one might be tempted not to go via InputForm here, but
            # it turns out to be safer to do it.
            n = self._unparsed_InputForm()[len("index("):]
            n = n[:n.find(")")]
            return self._get_sage_type(domain)(n)

        if head == "Fraction":
            return self.numer().sage()/self.denom().sage()

        if head == "Polynomial":
            base_ring = self._get_sage_type(domain[1])
            vars = self.variables()._unparsed_InputForm()[1:-1]
            if vars == "":
                return base_ring(self._unparsed_InputForm())
            else:
                R = PolynomialRing(base_ring, vars)
                return R(self._unparsed_InputForm())

        if head == "OrderedCompletion":
            # this is a workaround, I don't know how translate this
            if str(domain[1].car()) == "Expression":
                s = self._unparsed_InputForm()
                return SR(s.replace("pi()", "pi"))

        if head == "Expression":
            # TODO: we also have Expression Complex Integer and the like
            if str(domain[1].car()) == "Integer":
                s = self._unparsed_InputForm()
                return SR(s.replace("pi()", "pi"))

        if head == 'DistributedMultivariatePolynomial':
            base_ring = self._get_sage_type(domain[2])
            vars = domain[1].car()
            R = PolynomialRing(base_ring, vars)
            return R(self._unparsed_InputForm())



class FriCASFunctionElement(FunctionElement):
    def __init__(self, object, name):
        """Make FriCAS operation names valid python function identifiers.

        TESTS::

        sage: a = fricas('"Hello"')                                             # optional - fricas
        sage: a.upperCase_q                                                     # optional - fricas
        upperCase?
        sage: a.upperCase_e                                                     # optional - fricas
        upperCase!
        sage: a.upperCase_e()                                                   # optional - fricas
        "HELLO"

        """
        if name.endswith("_q"):
            name = name[:-2] + "?"
        elif name.endswith("_e"):
            name = name[:-2] + "!"
        FunctionElement.__init__(self, object, name)

    pass

class FriCASExpectFunction(ExpectFunction):
    def __init__(self, parent, name):
        """Translate the pythonized function identifier back to a FriCAS
        operation name.

        TESTS::

        sage: fricas.upperCase_q                                                # optional - fricas
        upperCase?
        sage: fricas.upperCase_e                                                # optional - fricas
        upperCase!

        """
        if name.endswith("_q"):
            name = name[:-2] + "?"
        elif name.endswith("_e"):
            name = name[:-2] + "!"
        ExpectFunction.__init__(self, parent, name)

    pass

def is_FriCASElement(x):
    """
    Return True if x is of type FriCASElement.

    EXAMPLES::

        sage: from sage.interfaces.fricas import is_FriCASElement               # optional - fricas
        sage: is_FriCASElement(fricas(2))                                       # optional - fricas
        True
        sage: is_FriCASElement(2)                                               # optional - fricas
        False
    """
    return isinstance(x, FriCASElement)

fricas = FriCAS()

def reduce_load_fricas():
    """
    Returns the FriCAS interface object defined in
    sage.interfaces.fricas.

    EXAMPLES::

        sage: from sage.interfaces.fricas import reduce_load_fricas             # optional - fricas
        sage: reduce_load_fricas()                                              # optional - fricas
        FriCAS
    """
    return fricas

import os

def fricas_console():
    """
    Spawn a new FriCAS command-line session.

    EXAMPLES::

        sage: fricas_console()                                                  # not tested
                         FriCAS (AXIOM fork) Computer Algebra System
                                    Version: FriCAS 1.0.5
                     Timestamp: Thursday February 19, 2009 at 06:57:33
        -----------------------------------------------------------------------------
           Issue )copyright to view copyright notices.
           Issue )summary for a summary of useful system commands.
           Issue )quit to leave AXIOM and return to shell.
        -----------------------------------------------------------------------------
    """
    from sage.repl.rich_output.display_manager import get_display_manager
    if not get_display_manager().is_in_terminal():
        raise RuntimeError('Can use the console only in the terminal. Try %%fricas magics instead.')
    os.system('fricas -nox')

def __doctest_cleanup():
    """
    EXAMPLES::

        sage: from sage.interfaces.fricas import __doctest_cleanup              # optional - fricas
        sage: a = FriCAS()                                                      # optional - fricas
        sage: two = a(2)                                                        # optional - fricas
        sage: a.is_running()                                                    # optional - fricas
        True
        sage: __doctest_cleanup()                                               # optional - fricas
        sage: a.is_running()                                                    # optional - fricas
        False
    """
    import sage.interfaces.quit
    sage.interfaces.quit.expect_quitall()
