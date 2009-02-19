r"""
Interface to FriCAS

TODO:

- Evaluation using a file is not done. Any input line with more than a
  few thousand characters would hang the system, so currently it
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

EXAMPLES: We evaluate a very simple expression in FriCas.

::

    sage: fricas('3 * 5')                     # optional - fricas
    15
    sage: a = fricas(3) * fricas(5); a        # optional - fricas
    15

The type of a is FriCASElement, i.e., an element of the axiom
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
axiom and returns the result as a string. This is exact as if we
typed in the given line of code to axiom; the return value is what
FriCAS would print out.

::

    sage: print fricas.eval('factor(x^5 - y^5)')   # optional - fricas
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

TESTS: We check to make sure the subst method works with keyword
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
from axiom import PanAxiom, PanAxiomElement, PanAxiomFunctionElement, PanAxiomExpectFunction
from pexpect import EOF
from sage.misc.multireplace import multiple_replace
from expect import Expect, ExpectElement, FunctionElement, ExpectFunction
from sage.misc.misc import verbose, DOT_SAGE, SAGE_ROOT


class FriCAS(PanAxiom):
    """
    The FriCAS interpreter interfaces.

    TESTS::

        sage: fricas == loads(dumps(fricas))
        True
    """
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

class FriCASElement(PanAxiomElement):
    pass

class FriCASFunctionElement(PanAxiomFunctionElement):
    pass

class FriCASExpectFunction(PanAxiomExpectFunction):
    pass

def is_FriCASElement(x):
    """
    Returns True of x is of type FriCASElement.

    EXAMPLES::

        sage: from sage.interfaces.fricas import is_FriCASElement
        sage: is_FriCASElement(fricas(2)) #optional - fricas
        True
        sage: is_FriCASElement(2)
        False
    """
    return isinstance(x, FriCASElement)

fricas = FriCAS(name='fricas', command='fricas -nox -noclef')

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
    os.system('fricas -nox')

def __doctest_cleanup():
    """
    EXAMPLES::

        sage: from sage.interfaces.fricas import __doctest_cleanup
        sage: a = FriCAS()
        sage: two = a(2)     #optional - fricas
        sage: a.is_running() #optional - fricas
        True
        sage: __doctest_cleanup()
        sage: a.is_running()
        False
    """
    import sage.interfaces.quit
    sage.interfaces.quit.expect_quitall()
