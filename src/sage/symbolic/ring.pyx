include "../ext/interrupt.pxi"
include "../ext/stdsage.pxi"
include "../ext/cdefs.pxi"

from sage.rings.integer cimport Integer
from sage.rings.real_mpfr import RealNumber

from expression cimport Expression, new_Expression_from_GEx

from sage.structure.element import RingElement

cdef class NSymbolicRing(Ring):
    """
    Symbolic Ring, parent object for all symbolic expressions.
    """
    def __init__(self):
        pass

    def _repr_(self):
        """
        Return a string representation of self.

        EXAMPLES:
            sage: from sage.symbolic.ring import NSR
            sage: NSR
            New Symbolic Ring
        """
        return "New Symbolic Ring"

    cdef _coerce_c_impl(self, other):
        """

        EXAMPLES:
        sage: from sage.symbolic.ring import NSR
        sage: NSR._coerce_(int(5))
        5
        sage: NSR._coerce_(5)
        5
        sage: NSR._coerce_(float(5))
        5.0
        sage: NSR._coerce_(5.0)
        5.0
        """
        cdef GEx exp

        if isinstance(other, (int, long)):
            GEx_construct_long(&exp, other)
            #GEx_construct_pyobject(exp, Integer(other))
        elif isinstance(other, float):
            GEx_construct_double(&exp, other)
            #GEx_construct_pyobject(exp, other)
        elif isinstance(other, Integer):
            GEx_construct_pyobject(exp, other)
            #GEx_construct_long(&exp,mpz_get_si((<Integer>other).value))
        elif isinstance(other, RealNumber):
            #GEx_construct_double(&exp, float(other))
            GEx_construct_pyobject(exp, other)
        elif isinstance(other, RingElement):
            GEx_construct_pyobject(exp, other)
        else:
            raise TypeError
        return new_Expression_from_GEx(exp)

    def __call__(self, other):
#TODO:
        try:
            return self._coerce_(other)
        except TypeError:
            raise TypeError, "conversion to %s from %s not defined."%(self, other.parent())


NSR = NSymbolicRing()

def var(name):
    """
    EXAMPLES:
        sage: from sage.symbolic.ring import var
        sage: var("x y z")
        (x, y, z)
        sage: var("x,y,z")
        (x, y, z)
        sage: var("x , y , z")
        (x, y, z)
        sage: var("z")
        z
    """
    if ',' in name:
        return tuple([new_symbol(s.strip()) for s in name.split(',')])
    if ' ' in name:
        return tuple([new_symbol(s.strip()) for s in name.split(' ')])
    return new_symbol(name)


def new_symbol(name):
    """
    EXAMPLES:
        sage: from sage.symbolic.ring import new_symbol
        sage: new_symbol("asdfasdfasdf")
        asdfasdfasdf
    """
    cdef GSymbol symb = get_symbol(name)
    cdef Expression e
    global NSR
    e = <Expression>PY_NEW(Expression)
    GEx_construct_symbol(&e._gobj, symb)
    e._parent = NSR
    return e


pi = new_Expression_from_GEx(g_Pi)
catalan = new_Expression_from_GEx(g_Catalan)
euler = new_Expression_from_GEx(g_Euler)

import sage.rings.integer
ginac_pyinit_Integer(sage.rings.integer.Integer)

import sage.rings.real_double
ginac_pyinit_Float(sage.rings.real_double.RDF)

import sage.rings.arith
ginac_pyinit_gcd(sage.rings.arith.gcd)

import sage.rings.arith
def f(x,y):
    print "x,y = ", (x,y)
    ans = sage.rings.arith.lcm(x,y)
    print "ans = ", ans
    return ans

#ginac_pyinit_lcm(f)
ginac_pyinit_lcm(sage.rings.arith.lcm)

import sage.rings.arith
ginac_pyinit_binomial(sage.rings.arith.binomial)
