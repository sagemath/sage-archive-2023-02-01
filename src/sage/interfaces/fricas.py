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

    sage: a.type()                           # optional - fricas
    PositiveInteger

We factor `x^5 - y^5` in FriCAS in several different ways.
The first way yields a FriCAS object.

::

    sage: F = fricas.factor('x^5 - y^5'); F      # optional - fricas
               4      3    2 2    3     4
    - (y - x)(y  + x y  + x y  + x y + x )
    sage: type(F)                               # optional - fricas
    <class 'sage.interfaces.fricas.FriCASElement'>
    sage: F.type()                              # optional - fricas
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
    - (y - x)(y  + x y  + x y  + x y + x )    Type: Factored(Polynomial(Integer))

We can create the polynomial `f` as a FriCAS polynomial,
then call the factor method on it. Notice that the notation
``f.factor()`` is consistent with how the rest of Sage
works.

::

    sage: f = fricas('x^5 - y^5')                  # optional - fricas
    sage: f^2                                     # optional - fricas
       10     5 5    10
      y   - 2x y  + x
    sage: f.factor()                              # optional - fricas
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
from axiom import PanAxiomElement, PanAxiomFunctionElement, PanAxiomExpectFunction
from sage.interfaces.tab_completion import ExtraTabCompletion
from expect import Expect, ExpectElement, FunctionElement, ExpectFunction
import re
import six
from sage.rings.all import RealField, ZZ, QQ

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

    def __getitem__(self, n):
        raise NotImplementedError

    def __int__(self):
        raise NotImplementedError

    def bool(self):
        raise NotImplementedError

    def __nonzero__(self):
        raise NotImplementedError

    def __long__(self):
        raise NotImplementedError

    def __float__(self):
        raise NotImplementedError

    def _integer_(self, ZZ=None):
        raise NotImplementedError

    def _rational_(self):
        raise NotImplementedError

    def gen(self, n):
        raise NotImplementedError

    def _get_1d_output(self):
        """Return the output from FriCAS as a string without the quotes.

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

    def _get_type_str(self):
        P = self._check_valid()
        output = P._eval_line(self._name + ";")

        m = re.match("\r\n\|startTypeTime\|\r\n *Type: (.*)\r\n\|endOfTypeTime\|", output, flags = re.DOTALL)
        return m.groups()[0]

    def _get_type(self):
        def to_list_of_lists(node):
            if isinstance(node, ast.Name):
                return node.id
            elif isinstance(node, ast.Num):
                return node.n
            elif isinstance(node, ast.List):
                return [to_list_of_lists(arg) for arg in node.elts]
            elif isinstance(node, ast.Call):
                return [node.func.id] + [to_list_of_lists(arg) for arg in node.args]

        # the fricas union type displays as Union(domain, ...), or Union(fail: failed, ...)
        # for ast.parse this would be a syntax error
        A = ast.parse(self._get_type_str().replace("...", "_").replace(": ", "_"))
        assert len(A.body) == 1
        return to_list_of_lists(A.body[0].value)

    def _get_sage_type(self, type):
        if type in ["Integer", "NonNegativeInteger", "PositiveInteger"]:
            return ZZ
        elif type == "String":
            return str
        elif type == "Float":
            prec = max(self.mantissa().length()._sage_(), 53)
            return RealField(prec)
        elif type == "DoubleFloat":
            return RDF
        elif type == "AlgebraicNumber":
            return QQbar
        elif isinstance(type, list):
            if type[0] == "Union":
                return self._get_sage_type(type[1])

            elif type[0] == "IntegerMod":
                return Integers(ZZ(type[1]))

            elif type[0] == "Fraction":
                return FractionField(self._get_sage_type(type[1]))

            elif type[0] == "Expression":
                return SR

            elif type[0] == "Polynomial":
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

        sage: a = fricas(x^2 + 1); a._get_type()                  # optional - fricas
        ['Polynomial', 'Integer']
        sage: a.sage()                                            # optional - fricas
        x^2 + 1
        sage: _.parent()                                          # optional - fricas
        Univariate Polynomial Ring in x over Integer Ring
        sage: fricas('x^2 + y^2 + 1/2').sage()                    # optional - fricas
        y^2 + x^2 + 1/2
        sage: _.parent()                                          # optional - fricas
        Multivariate Polynomial Ring in y, x over Rational Field

        Rational functions:

        sage: fricas("x^2 + 1/z").sage()                          # optional - fricas
        x^2 + 1/z

        Expressions:

        sage: fricas("sin(x+y)/exp(z)*log(1+%e)").sage()          # optional - fricas
        e^(-z)*log(e + 1)*sin(x + y)

        sage: fricas("factorial(n)").sage()                       # optional - fricas
        factorial(n)

        sage: fricas("guessPade [1,1,2,3,5,8]")                   # optional - fricas
            n        1
        [[[x ]- ----------]]
                 2
                x  + x - 1

        """
        type = self._get_type()
        if type in ["Integer", "NonNegativeInteger", "PositiveInteger"]:
            return ZZ(self._get_1d_output())

        elif type == "String":
            return self._get_1d_output()

        elif type == "Float":
            prec = max(self.mantissa().length().sage(), 53)
            R = RealField(prec)
            x, e, b = self._get_1d_output().lstrip('float(').rstrip(')').split(',')
            return R(ZZ(x)*ZZ(b)**ZZ(e))

        elif type == "DoubleFloat":
            from sage.rings.all import RDF
            return RDF(self._get_1d_output())

        elif type == "AlgebraicNumber":
            from sage.rings.all import QQbar
            s = self._get_1d_output()[:-len("::AlgebraicNumber()")]
            return sage.misc.sage_eval.sage_eval("QQbar(" + s + ")")

        elif isinstance(type, list):
            if type[0] == "List":
                P = self._check_valid()
                n = P.new('# %s'%self.name()).sage()
                return [P.new('%s(%s)'%(self._name, k)).sage() for k in range(1, n+1)]

            elif type[0] == "IntegerMod":
                # one might be tempted not to go via unparse and
                # InputForm here, but it turns out to be safer to do
                # it.
                n = self._get_1d_output()[len("index("):]
                n = n[:n.find(")")]
                return self._get_sage_type(type)(n)

            elif type[0] == "Fraction":
                return self.numer().sage()/self.denom().sage()

            elif type[0] == "Polynomial":
                from sage.rings.all import PolynomialRing
                base_ring = self._get_sage_type(type[1])
                vars = self.variables()._get_1d_output()[1:-1]
                R = PolynomialRing(base_ring, vars)
                return R(self._get_1d_output())

            elif type[0] == "Expression":
                if type[1] == "Integer":
                    s = self._get_1d_output()
                    return SR(s)

            elif type[0] == 'DistributedMultivariatePolynomial':
                from sage.rings.all import PolynomialRing
                base_ring = self._get_sage_type(type[2])
                vars = type[1]
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
