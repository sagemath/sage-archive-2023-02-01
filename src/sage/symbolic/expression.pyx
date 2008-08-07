include "../ext/interrupt.pxi"
include "../ext/stdsage.pxi"
include "../ext/cdefs.pxi"

import ring

from sage.structure.element cimport ModuleElement, RingElement, Element

cdef class Expression(CommutativeRingElement):
    def __dealloc__(self):
        GEx_destruct(&self._gobj)

    def _repr_(self):
        """
            sage: var("x y", ns=1)
            (x, y)
            sage: x+y
            x+y
        """
        return GEx_to_str(&self._gobj)

    def __hash__(self):
        """
            sage: var("x y", ns=1)
            (x, y)
            sage: hash(x+y)
            46142460
            sage: d = {x+y: 5}
            sage: d
            {x+y: 5}
        """
        return self._gobj.gethash()

    cdef ModuleElement _add_c_impl(left, ModuleElement right):
        """
            sage.: var("x y", ns=1)
            (x, y)
            sage.: x+y+y+x
            2*x+2*y
        """
        _sig_on
        cdef GEx e = gadd(left._gobj, (<Expression>right)._gobj)
        _sig_off
        return new_Expression_from_GEx(e)

    cdef RingElement _mul_c_impl(left, RingElement right):
        """
            sage: var("x y", ns=1)
            (x, y)
            sage: x*y*y
            x*y^2
        """
        _sig_on
        cdef GEx e = gmul(left._gobj, (<Expression>right)._gobj)
        _sig_off
        return new_Expression_from_GEx(e)

    def __richcmp__(left, right, int op):
        #boilerplate from sage.structure.element
        return (<Element>left)._richcmp(right, op)

    cdef int _cmp_c_impl(left, Element right) except -2:
        return left._gobj.compare((<Expression>right)._gobj)

    def __pow__(Expression self, exp, ignored):
        cdef Expression nexp = self._parent._coerce_c(exp)
        _sig_on
        return new_Expression_from_GEx(gpow(self._gobj, nexp._gobj))
        _sig_off

    def expand(Expression self):
        _sig_on
        return new_Expression_from_GEx(self._gobj.expand(0))
        _sig_off

    def collect(Expression self, Expression s):
#FIXME convert second argument if necessary
        _sig_on
        return new_Expression_from_GEx(self._gobj.collect(s._gobj, False))
        _sig_off

cdef Expression new_Expression_from_GEx(GEx juice):
    cdef Expression nex
    nex = <Expression>PY_NEW(Expression)
    GEx_construct_ex(&nex._gobj, juice)
    nex._parent = ring.NSR
    return nex
