"""
Evaluating a string in \sage
"""

#*****************************************************************************
#       Copyright (C) 2006 William Stein <wstein@ucsd.edu>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#
#                  http://www.gnu.org/licenses/
#*****************************************************************************

from preparser import preparse

def sage_eval(_obj_, extra_locals={}):
    r"""
    Obtain a \SAGE object from the input string by evaluate it using
    SAGE (including preparsing and in the context of the \SAGE
    library, etc.).

    All preparsing is applied and the expression is evaluated in the
    context of \code{from sage.all import *}.  To evaluate the
    expression with certain variables set, use the extra_locals
    argument, which should be a dictionary.

    EXAMPLES:
    This example illustrates that preparsing is applied.
        sage: eval('2^3')
        1
        sage: sage_eval('2^3')
        8

    This illustrates interfaces:
        sage: f = gp('2/3')
        sage: type(f)
        <class 'sage.interfaces.gp.GpElement'>
        sage: f._sage_()
        2/3
        sage: type(f._sage_())
        <type 'rational.Rational'>
        sage: a = gap(939393/2433)
        sage: a._sage_()
        313131/811
        sage: type(a._sage_())
        <type 'rational.Rational'>

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
        sage: eval('4/3 + x',  {'x':25})
        26
        sage: sage_eval('4/3 + x',  {'x':25})
        79/3

    This example illustrates how \code{sage_eval} can be useful
    when evaluating the output of other computer algebra systems.

        sage: x = PolynomialRing(Q).gen()
        sage: gap.eval('R:=PolynomialRing(Rationals,["x"]);')
        'PolynomialRing(..., [ x ])'
        sage: ff = gap.eval('x:=IndeterminatesOfPolynomialRing(R);; f:=x^2+1;')
        sage: ff; sage_eval(ff, {'x':x})
        'x^2+1'
        x^2 + 1
        sage: eval(ff)
        Traceback (most recent call last):
        ...
        RuntimeError: Use ** for exponentiation, not '^', which means xor
        in Python, and has the wrong precedence.

    Here you can see eval simply will not work but \code{sage_eval} will.
    """
    if not isinstance(_obj_, str):
        try:
            return _obj_._sage_()
        except (RuntimeError, NotImplementedError, AttributeError):
            _obj_ = str(_obj_)

    # The following line generates a warning when the file is compiled.
    # However, I can't think of any way around having it here.  I could
    # move sage_eval to sage/all.py, but then no part of sage could
    # use it, which would make it useless.
    from sage.all import *
    L = locals()
    for k, x in extra_locals.items():
        L[k] = x
    return eval(preparse(_obj_), globals(), L)
