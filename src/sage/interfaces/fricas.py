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

    sage: fricas('3 * 5')                     # optional - fricas
    15
    sage: a = fricas(3) * fricas(5); a        # optional - fricas
    15

The type of a is FriCASElement, i.e., an element of the FriCAS
interpreter.

::

    sage: type(a)                            # optional - fricas
    <class 'sage.interfaces.fricas.FriCASElement'>
    sage: a.parent()                          # optional - fricas
    FriCAS

The underlying FriCAS type of a is also available, via the type
method::

    sage: a.typeOf()                             # optional - fricas
    PositiveInteger

We factor `x^5 - y^5` in FriCAS in several different ways.
The first way yields a FriCAS object.

::

    sage: F = fricas.factor('x^5 - y^5'); F      # optional - fricas
               4      3    2 2    3     4
    - (y - x)(y  + x y  + x y  + x y + x )
    sage: type(F)                               # optional - fricas
    <class 'sage.interfaces.fricas.FriCASElement'>
    sage: F.typeOf()                              # optional - fricas
    Factored(Polynomial(Integer))

Note that FriCAS objects are normally displayed using "ASCII art".

::

    sage: a = fricas(2/3); a          # optional - fricas
      2
      -
      3
    sage: a = fricas('x^2 + 3/7'); a      # optional - fricas
       2   3
      x  + -
           7

The ``fricas.eval`` command evaluates an expression in
FriCAS and returns the result as a string. This is exact as if we
typed in the given line of code to FriCAS; the return value is what
FriCAS would print out.

::

    sage: print(fricas.eval('factor(x^5 - y^5)'))   # optional - fricas
               4      3    2 2    3     4
    - (y - x)(y  + x y  + x y  + x y + x )

We can create the polynomial `f` as a FriCAS polynomial,
then call the factor method on it. Notice that the notation
``f.factor()`` is consistent with how the rest of Sage
works.

::

    sage: f = fricas('x^5 - y^5')                  # optional - fricas
    sage: f^2                                      # optional - fricas
       10     5 5    10
      y   - 2x y  + x
    sage: f.factor()                               # optional - fricas
               4      3    2 2    3     4
    - (y - x)(y  + x y  + x y  + x y + x )

Control-C interruption works well with the FriCAS interface. For
example, try the following sum but with a much bigger range, and hit
control-C.

::

    sage:  f = fricas('(x^5 - y^5)^10000')       # not tested - fricas
    Interrupting FriCAS...
    ...
    <type 'exceptions.TypeError'>: Ctrl-c pressed while running FriCAS

::

    sage: fricas('1/100 + 1/101')                  # optional - fricas
       201
      -----
      10100
    sage: a = fricas('(1 + sqrt(2))^5'); a         # optional - fricas
         +-+
      29\|2  + 41

TESTS:

We check to make sure the subst method works with keyword
arguments.

::

    sage: a = fricas(x+2); a  #optional - fricas
    x + 2
    sage: a.subst(x=3)       #optional - fricas
    5

We verify that FriCAS floating point numbers can be converted to
Python floats.

::

    sage: float(fricas(2))     #optional - fricas
    2.0
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
# from __future__ import print_function
# from __future__ import absolute_import

from sage.interfaces.axiom import PanAxiomElement, PanAxiomFunctionElement, PanAxiomExpectFunction
from sage.interfaces.tab_completion import ExtraTabCompletion
from sage.interfaces.expect import Expect, ExpectElement, FunctionElement, ExpectFunction
from sage.env import DOT_SAGE
import re
import six

FRICAS_SINGLE_LINE_START = 3      # where the output starts
FRICAS_MULTI_LINE_START = 2
FRICAS_LINE_LENGTH = 245   # length of a line
FRICAS_LINE_BREAK = "\r\n" # line ending
# beware that lisp distinguishes between ' and ".
FRICAS_INIT_CODE = (
")lisp (princ \"running FriCAS intialisation code\")",
")set functions compile on",
")set output length " + str(FRICAS_LINE_LENGTH),
")set message autoload off",
")lisp (setf |$ioHook|"
"            (lambda (x &optional args)"
"              (when (member x '(|startAlgebraOutput| |endOfAlgebraOutput|"
"                                |startKeyedMsg|      |endOfKeyedMsg|))"
"               (let ((typetime (member (car args) '(S2GL0012 S2GL0013 S2GL0014))))"
"                 (cond"
"                   ((and (eq x '|startKeyedMsg|) typetime) (princ \"|startTypeTime|\"))"
"                   ((and (eq x '|endOfKeyedMsg|) typetime) (princ \"|endOfTypeTime|\"))"
"                   (t (prin1 x))))"
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
        Expect.__init__(self,
                        name = name,
                        prompt = FRICAS_FIRST_PROMPT,
                        command = command,
                        script_subdirectory = script_subdirectory,
                        server=server,
                        server_tmpdir=server_tmpdir,
                        restart_on_ctrlc = False,
                        verbose_start = False,
                        init_code = FRICAS_INIT_CODE,
                        logfile = logfile,
                        eval_using_file_cutoff=eval_using_file_cutoff)

    def _start(self):
        """
        Start the FriCAS interpreter and switch off the linenumbers.

        EXAMPLES::

            sage: a = FriCAS()
            sage: a.is_running()
            False
            sage: a._start()     #optional - axiom
            sage: a.is_running() #optional - axiom
            True
            sage: a.quit()       #optional - axiom
        """
        Expect._start(self)
        # switching off the line numbers also modified the prompt
        self._prompt = FRICAS_LINENUMBER_OFF_PROMPT
        self.eval(FRICAS_LINENUMBER_OFF_CODE)

    def eval(self, code, strip=True, synchronize=False, locals=None, allow_use_file=True,
             split_lines="nofile", **kwds):
        """Send the code to the FriCAS interpreter and return the pretty
        printed output as a string.

        INPUT:

        -  ``code``       -- text to evaluate

        """
        output = Expect.eval(self, code, strip=strip, synchronize=synchronize,
                             locals=locals, allow_use_file=allow_use_file, split_lines=split_lines,
                             **kwds)
#        print "running eval", code
#        print output
        m = re.match("\r\n\|startAlgebraOutput\|\r\n(.*)\r\n\|endOfAlgebraOutput\|", output, flags = re.DOTALL)
        if m:
            lines = m.groups()[0].split(FRICAS_LINE_BREAK)
            if max(len(line) for line in lines) < FRICAS_LINE_LENGTH:
                return "\r\n".join(line[FRICAS_SINGLE_LINE_START:] for line in lines)
            else:
                return "\r\n".join(line[FRICAS_MULTI_LINE_START:] for line in lines)

        else:
#            print "TODO"
            return output

    def set(self, var, value):
        """Set the variable var to the given value.

        INPUT:

        - var, value: strings, the first representing a valid FriCAS
          variable identifier, the second a FriCAS expression.

        OUTPUT: None

        EXAMPLES::

            sage: fricas.set('xx', '2')    #optional - fricas
            sage: fricas.get('xx')         #optional - fricas
            '2'

        """
#        print "running set", var, value
        if not isinstance(var, str) or not isinstance(value, str):
            raise TypeError

        self._eval_line('%s := %s;'%(var, value))

    def get(self, var):
        r"""
        Get the string representation of the value of a variable.

        This is only used for display.

        EXAMPLES::

            sage: fricas.set('xx', '2')    #optional - fricas
            sage: fricas.get('xx')         #optional - fricas
            '2'
            sage: a = fricas('(1 + sqrt(2))^5') #optional - fricas
            sage: fricas.get(a.name())          #optional - fricas
            '   +-+\r\n29\\|2  + 41'
        """
#        print "running get", var
        return self.eval(str(var))

    def _assign_symbol(self):
        """Returns the symbol used for setting a variable.

        This would be used if `self.set` from Interface would not be
        overloaded.
        """
        return ":="

    def _equality_symbol(self):
        """Returns the equality testing logical symbol in FriCAS.

        EXAMPLES::

            sage: a = fricas(x==6); a   #optional fricas
            x= 6
        """
        return "="

    def __repr__(self):
        """
        EXAMPLES::

            sage: fricas
            FriCAS
        """
        return "FriCAS"

    def __reduce__(self):
        """
        EXAMPLES::

            sage: fricas.__reduce__()
            (<function reduce_load_fricas at 0x...>, ())
            sage: f, args = _
            sage: f(*args)
            FriCAS
        """
        return reduce_load_fricas, tuple([])

    def _function_class(self):
        """
        Return the FriCASExpectFunction class.

        EXAMPLES::

            sage: fricas._function_class()
            <class 'sage.interfaces.fricas.FriCASExpectFunction'>
            sage: type(fricas.gcd)
            <class 'sage.interfaces.fricas.FriCASExpectFunction'>
        """
        return FriCASExpectFunction

    def _object_class(self):
        """
        EXAMPLES::

            sage: fricas._object_class()
            <class 'sage.interfaces.fricas.FriCASElement'>
            sage: type(fricas(2)) #optional - fricas
            <class 'sage.interfaces.fricas.FriCASElement'>
        """
        return FriCASElement

    def _function_element_class(self):
        """
        Returns the FriCAS function element class.

        EXAMPLES::

            sage: fricas._function_element_class()
            <class 'sage.interfaces.fricas.FriCASFunctionElement'>
            sage: type(fricas(2).gcd) #optional - fricas
            <class 'sage.interfaces.fricas.FriCASFunctionElement'>
        """
        return FriCASFunctionElement

    def console(self):
        """
        Spawn a new FriCAS command-line session.

        EXAMPLES::

            sage: fricas.console() #not tested
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

            sage: v = axiom('[x^i for i in 0..5]')            # optional - axiom
            sage: len(v)                                      # optional - axiom
            6
        """
        P = self._check_valid()
        l = P('# %s ' %self._name)
        return l.sage()

    def __getitem__(self, n):
        """
        We implement the sage conventions here!

        TODO: should we really check the length?

        TODO: we should also implement negative arguments and tuples
        """
        n = int(n)
        if n < 0 or n >= len(self):
            raise IndexError("index out of range")
        P = self._check_valid()
        # in FriCAS, retrieving an item is the same as calling it with an integer
        return P('%s(%s)'%(self._name, n+1))

    def __int__(self):
        return int(self.sage())

    def bool(self):
        raise NotImplementedError

    def __nonzero__(self):
        raise NotImplementedError

    def __long__(self):
        raise NotImplementedError

    def __float__(self):
        return float(self.sage())

    def _integer_(self, ZZ=None):
        raise NotImplementedError

    def _rational_(self):
        raise NotImplementedError

    def gen(self, n):
        raise NotImplementedError

    def _get_1d_output(self):
        """Return the output from FriCAS as a string without the quotes.

        TODO:

        - catch errors, especially when InputForm is not available:

        -- for example when integration returns "failed"
        -- UnivariatePolynomial

        should we provide workarounds, too?

        TESTS:

        We test that strings are returned properly:

        sage: r = fricas('concat([concat(string(i)," ") for i in 0..299])')._get_1d_output()
        sage: r == " ".join([str(i) for i in range(300)]) + ' '
        True

        sage: fricas('concat([string(1) for i in 1..5])')._get_1d_output() == "1"*5
        True

        sage: fricas('concat([string(1) for i in 1..10000])')._get_1d_output() == "1"*10000
        True
        """
        P = self._check_valid()
        return P.eval('unparse(%s::InputForm)' %self._name).replace("\r\n", "")[1:-1]

    def _get_sage_type(self, type):
        """
        INPUT:

        - type, a FriCAS SExpression

        OUTPUT:

        - a corresponding Sage type
        """
        from sage.rings.all import ZZ, QQ, QQbar, PolynomialRing, RDF
        from sage.rings.fraction_field import FractionField
        from sage.rings.finite_rings.integer_mod_ring import Integers
        from sage.rings.real_mpfr import RealField
        from sage.symbolic.ring import SR

        # first implement domains without arguments
        head = str(type.car())
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
        if head == "Union":
            return self._get_sage_type(type[1][2])

        if head == "OrderedCompletion":
            # this is a workaround, I don't know how translate this
            return SR
        
        if head == "IntegerMod":
            return Integers(type[1].integer().sage())

        if head == "Fraction":
            return FractionField(self._get_sage_type(type[1]))

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

        sage: fricas(2.1234).sage()                               # optional - fricas
        2.12340000000000
        sage: _.parent()                                          # optional - fricas
        Real Field with 53 bits of precision
        sage: a = RealField(100)(pi)
        sage: fricas(a).sage()                                    # optional - fricas
        3.1415926535897932384626433833
        sage: _.parent()                                          # optional - fricas
        Real Field with 100 bits of precision
        sage: fricas(a).sage() == a                               # optional - fricas
        True
        sage: fricas(2.0).sage()                                  # optional - fricas
        2.00000000000000
        sage: _.parent()                                          # optional  - fricas
        Real Field with 53 bits of precision

        Algebraic numbers:

        sage: a = fricas('(1 + sqrt(2))^5'); a                    # optional - fricas
           +-+
        29\|2  + 41
        sage: b = a.sage(); b                                     # optional - fricas
        82.0121933088198?
        sage: b.radical_expression()                              # optional - fricas
        29*sqrt(2) + 41

        Integers modulo n:

        sage: fricas("((42^17)^1783)::IntegerMod(5^(5^5))").sage() == Integers(5^(5^5))((42^17)^1783) # optional - fricas
        True

        We can also convert Fricas's polynomials to Sage polynomials.

        sage: a = fricas(x^2 + 1); a.typeOf()                     # optional - fricas
        Polynomial(Integer)
        sage: a.sage()                                            # optional - fricas
        x^2 + 1
        sage: _.parent()                                          # optional - fricas
        Univariate Polynomial Ring in x over Integer Ring
        sage: fricas('x^2 + y^2 + 1/2').sage()                    # optional - fricas
        y^2 + x^2 + 1/2
        sage: _.parent()                                          # optional - fricas
        Multivariate Polynomial Ring in y, x over Rational Field

        sage: fricas("1$Polynomial Integer").sage()               # optional - fricas
        1

        sage: fricas("x^2/2").sage()                              # optional - fricas
        1/2*x^2

        Rational functions:

        sage: fricas("x^2 + 1/z").sage()                          # optional - fricas
        x^2 + 1/z

        Expressions:

        sage: fricas("sin(x+y)/exp(z)*log(1+%e)").sage()          # optional - fricas
        e^(-z)*log(e + 1)*sin(x + y)

        sage: fricas("factorial(n)").sage()                       # optional - fricas
        factorial(n)

        sage: fricas("integrate(sin(x+y), x=0..1)").sage()        # optional - fricas
        -cos(y + 1) + cos(y)

        sage: fricas("integrate(x*sin(1/x), x=0..1)").sage()      # optional - fricas

        sage: fricas("guessPade [1,1,2,3,5,8]")                   # optional - fricas
            n        1
        [[[x ]- ----------]]
                 2
                x  + x - 1

        """
        from sage.rings.all import ZZ, QQ, QQbar, PolynomialRing, RDF
        from sage.rings.fraction_field import FractionField
        from sage.rings.finite_rings.integer_mod_ring import Integers
        from sage.rings.real_mpfr import RealField
        from sage.symbolic.ring import SR
        from sage.misc.sage_eval import sage_eval

        domain = self.dom() # type is now a fricas SExpression
        head = str(domain.car())
        # the special domain Union needs special treatment
        if head == "Union":
            domain = domain[1][2]
            head = str(domain.car())

        # first implement domains without arguments
        if head in ["Integer", "NonNegativeInteger", "PositiveInteger"]:
            return ZZ(self._get_1d_output())

        if head == "String":
            return self._get_1d_output()

        if head == "Float":
            prec = max(self.mantissa().length().sage(), 53)
            R = RealField(prec)
            x, e, b = self._get_1d_output().lstrip('float(').rstrip(')').split(',')
            return R(ZZ(x)*ZZ(b)**ZZ(e))

        if head == "DoubleFloat":
            return RDF(self._get_1d_output())

        if head == "AlgebraicNumber":
            s = self._get_1d_output()[:-len("::AlgebraicNumber()")]
            return sage_eval("QQbar(" + s + ")")

        # now implement "functorial" types
        
        if head == "List":
            P = self._check_valid()
            n = P.new('# %s'%self.name()).sage()
            return [P.new('%s(%s)'%(self._name, k)).sage() for k in range(1, n+1)]

        if head == "IntegerMod":
            # one might be tempted not to go via unparse and
            # InputForm here, but it turns out to be safer to do
            # it.
            n = self._get_1d_output()[len("index("):]
            n = n[:n.find(")")]
            return self._get_sage_type(domain)(n)

        if head == "Fraction":
            return self.numer().sage()/self.denom().sage()

        if head == "Polynomial":
            base_ring = self._get_sage_type(domain[1])
            vars = self.variables()._get_1d_output()[1:-1]
            if vars == "":
                return base_ring(self._get_1d_output())
            else:
                R = PolynomialRing(base_ring, vars)
                return R(self._get_1d_output())

        if head == "OrderedCompletion":
            # this is a workaround, I don't know how translate this
            if str(domain[1].car()) == "Expression":
                s = self._get_1d_output()
                return SR(s)
            
        if head == "Expression":
            # TODO: we also have Expression Complex Integer and the like
            if str(domain[1].car()) == "Integer":
                s = self._get_1d_output()
                return SR(s)

        if head == 'DistributedMultivariatePolynomial':
            base_ring = self._get_sage_type(domain[2])
            vars = domain[1].car()
            R = PolynomialRing(base_ring, vars)
            return R(self._get_1d_output())



class FriCASFunctionElement(PanAxiomFunctionElement):
    pass

class FriCASExpectFunction(PanAxiomExpectFunction):
    pass

def is_FriCASElement(x):
    """
    Return True if x is of type FriCASElement.

    EXAMPLES::

        sage: from sage.interfaces.fricas import is_FriCASElement #optional - fricas
        sage: is_FriCASElement(fricas(2))  #optional - fricas
        True
        sage: is_FriCASElement(2)  #optional - fricas
        False
    """
    return isinstance(x, FriCASElement)

fricas = FriCAS()

def reduce_load_fricas():
    """
    Returns the FriCAS interface object defined in
    sage.interfaces.fricas.

    EXAMPLES::

        sage: from sage.interfaces.fricas import reduce_load_fricas
        sage: reduce_load_fricas()
        FriCAS
    """
    return fricas

import os

def fricas_console():
    """
    Spawn a new FriCAS command-line session.

    EXAMPLES::

        sage: fricas_console() #not tested
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

        sage: from sage.interfaces.fricas import __doctest_cleanup #optional - fricas
        sage: a = FriCAS()   #optional - fricas
        sage: two = a(2)     #optional - fricas
        sage: a.is_running() #optional - fricas
        True
        sage: __doctest_cleanup()   #optional - fricas
        sage: a.is_running()   #optional - fricas
        False
    """
    import sage.interfaces.quit
    sage.interfaces.quit.expect_quitall()
