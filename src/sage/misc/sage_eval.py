r"""
Evaluating a String in SAGE
"""

#*****************************************************************************
#       Copyright (C) 2006 William Stein <wstein@gmail.com>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#
#                  http://www.gnu.org/licenses/
#*****************************************************************************

from preparser import preparse

def sage_eval(source, locals={}):
    r"""
    Obtain a \SAGE object from the input string by evaluate it using
    SAGE.  This means calling eval after preparsing and with
    globals equal to everything included in the scope of
    \code{from sage.all import *}.).

    If the object has an _sage_ method it is called and the value is
    returned.  Otherwise str is called on the object, and all
    preparsing is applied and the resulting expression is evaluated in the
    context of \code{from sage.all import *}.  To evaluate the
    expression with certain variables set, use the locals argument,
    which should be a dictionary.

    EXAMPLES:
    This example illustrates that preparsing is applied.
        sage: eval('2^3')
        1
        sage: sage_eval('2^3')
        8

    Note that you have explicitly define variables and pass
    them as the second option:

        sage: x = PolynomialRing(RationalField(),"x").gen()
        sage: sage_eval('x^2+1', locals={'x':x})
        x^2 + 1

    This illustrates interfaces:
        sage: f = gp('2/3')
        sage: type(f)
        <class 'sage.interfaces.gp.GpElement'>
        sage: f._sage_()
        2/3
        sage: type(f._sage_())
        <type 'sage.rings.rational.Rational'>
        sage: a = gap(939393/2433)
        sage: a._sage_()
        313131/811
        sage: type(a._sage_())
        <type 'sage.rings.rational.Rational'>

    This example illustrates that evaluation occurs in the context of
    \code{from sage.all import *}.  Even though bernoulli has been
    redefined in the local scope, when calling \code{sage_eval} the
    default value meaning of bernoulli is used.  Likewise for QQ
    below.

        sage: bernoulli = lambda x : x^2
        sage: bernoulli(6)
        36
        sage: eval('bernoulli(6)')
        36
        sage: sage_eval('bernoulli(6)')
        1/42

        sage: QQ = lambda x : x^2
        sage: QQ(2)
        4
        sage: sage_eval('QQ(2)')
        2
        sage: parent(sage_eval('QQ(2)'))
        Rational Field

    This example illustrates setting a variable for use in evaluation.
        sage: x = 5
        sage: eval('4/3 + x', {'x':25})
        26
        sage: sage_eval('4/3 + x',  locals={'x':25})
        79/3

    This example illustrates how \code{sage_eval} can be useful
    when evaluating the output of other computer algebra systems.

        sage: R.<x> = PolynomialRing(RationalField())
        sage: gap.eval('R:=PolynomialRing(Rationals,["x"]);')
        'PolynomialRing(..., [ x ])'
        sage: ff = gap.eval('x:=IndeterminatesOfPolynomialRing(R);; f:=x^2+1;'); ff
        'x^2+1'
        sage: sage_eval(ff, locals={'x':x})
        x^2 + 1
        sage: eval(ff)
        Traceback (most recent call last):
        ...
        RuntimeError: Use ** for exponentiation, not '^', which means xor
        in Python, and has the wrong precedence.

    Here you can see eval simply will not work but \code{sage_eval} will.
    """
    if not isinstance(source, str):
        raise TypeError, "source must be a string."

    import sage.all
    p = preparse(source)
    try:
        return eval(p, sage.all.__dict__, locals)
    except SyntaxError, msg:
        raise SyntaxError, "%s\nError using SAGE to evaluate '%s'"%(msg, p)



def sageobj(x, vars=None):
    """
    Return a native SAGE object associated to x, if possible
    and implemented.

    If x is a string it is evaluated with SAGE preparsing.

    EXAMPLES:
        sage: type(sageobj(gp('34/56')))
        <type 'sage.rings.rational.Rational'>
        sage: n = 5/2
        sage: sageobj(n) is n
        True
        sage: k = sageobj('Z(8^3/1)', {'Z':ZZ}); k
        512
        sage: type(k)
        <type 'sage.rings.integer.Integer'>
    """
    try:
       return x._sage_()
    except (TypeError, NotImplementedError, AttributeError):
       return sage_eval(str(x), vars)
