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

include "../../ext/interrupt.pxi"
include "../../ext/stdsage.pxi"
include 'misc.pxi'
include 'decl.pxi'

ZZ_pContextDict = {}

from sage.libs.ntl.ntl_ZZ cimport ntl_ZZ


cdef class ntl_ZZ_pContext_class:
    def __init__(self, ntl_ZZ v):
        """
        EXAMPLES:
            # You can construct contexts manually.
            sage: c=ntl.ZZ_pContext(ntl.ZZ(11))
            sage: n1=c.ZZ_p(12)
            sage: n1
            1

            # or You can construct contexts implicitly.
            sage: n2=ntl.ZZ_p(12, 7)
            sage: n2
            5
            sage: ntl.ZZ_p(2,3)+ntl.ZZ_p(1,3)
            0
            sage: n2+n1  # Mismatched moduli:  It will go BOOM!
            Traceback (most recent call last):
            ...
            ValueError: You can not perform arithmetic with elements of different moduli.
        """
        pass

    def __new__(self, ntl_ZZ v):
        ZZ_pContext_construct_ZZ(&self.x, &v.x)
        ZZ_pContextDict[repr(v)] = self
        self.p = v
        self.p_bits = self.p._integer_().bits()

    def __dealloc__(self):
        ZZ_pContext_destruct(&self.x)

    def __reduce__(self):
        """
        sage: c=ntl.ZZ_pContext(ntl.ZZ(13))
        sage: loads(dumps(c)) is c
        True
        """
        return ntl_ZZ_pContext, (self.p,)

    def restore(self):
        self.restore_c()

    cdef void restore_c(self):
        self.x.restore()

    def ZZ_p(self,v = None):
        from ntl_ZZ_p import ntl_ZZ_p
        return ntl_ZZ_p(v,modulus=self)

    def ZZ_pX(self,v = None):
        from ntl_ZZ_pX import ntl_ZZ_pX
        return ntl_ZZ_pX(v,modulus=self)

def ntl_ZZ_pContext( ntl_ZZ v ):
    try:
        return ZZ_pContextDict[repr(v)]
    except KeyError:
        return ntl_ZZ_pContext_class(v)
