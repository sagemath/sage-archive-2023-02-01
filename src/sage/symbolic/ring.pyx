include "../ext/interrupt.pxi"
include "../ext/stdsage.pxi"
include "../ext/cdefs.pxi"

#################################################################
# Initialize the library
#################################################################

#initialize_ginac()

from sage.rings.integer cimport Integer
from sage.rings.real_mpfr import RealNumber

from expression cimport Expression, new_Expression_from_GEx

from sage.structure.element import RingElement

cdef class NSymbolicRing(Ring):
    """
    Symbolic Ring, parent object for all symbolic expressions.
    """
    def __init__(self):
        """
        Initialize the New Symbolic Ring.

        EXAMPLES:
            sage: sage.symbolic.ring.NSymbolicRing()
            New Symbolic Ring
        """

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
            5.00000000000000

        An interval arithmetic number:
            sage: NSR._coerce_(RIF(pi))
            3.141592653589794?

        A number modulo 7:
            sage: a = NSR._coerce_(Mod(3,7)); a
            3
            sage: a^2
            2
        """
        cdef GEx exp

        if isinstance(other, int):
            GEx_construct_long(&exp, other)
        elif isinstance(other, float):
            GEx_construct_double(&exp, other)
        elif isinstance(other, (Integer, long, complex)):
            GEx_construct_pyobject(exp, other)
        elif isinstance(other, RealNumber):
            GEx_construct_pyobject(exp, other)
        elif isinstance(other, RingElement):
            GEx_construct_pyobject(exp, other)
        else:
            raise TypeError
        return new_Expression_from_GEx(exp)

    def __call__(self, other):
        """
        INPUT:
            other -- python object

        EXAMPLES:
            sage: from sage.symbolic.ring import NSR
            sage: NSR(2/5)
            2/5
            sage: NSR.__call__(2/5)
            2/5
            sage: NSR.__call__('foo')
            Traceback (most recent call last):
            ...
            TypeError: conversion not defined
        """
        try:
            return self._coerce_(other)
        except TypeError:
            raise TypeError, "conversion not defined"


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

from sage.rings.complex_double import CDF
ginac_pyinit_I(CDF.gen())
