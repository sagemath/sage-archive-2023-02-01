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

zz_pContextDict = {}

cdef class ntl_zz_pContext_class:
    def __init__(self, long v):
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

    def __new__(self, long v):
        ## TODO: v smaller than NTL_SP_BOUND
        zz_pContext_construct_long(&self.x, v)
        zz_pContextDict[repr(v)] = self
        self.p = v

    def __dealloc__(self):
        zz_pContext_destruct(&self.x)

    def __reduce__(self):
        """
        sage: c=ntl.ZZ_pContext(ntl.ZZ(13))
        sage: loads(dumps(c)) is c
        True
        """
        return ntl_zz_pContext, (self.p,)

    def restore(self):
        self.restore_c()

    cdef void restore_c(self):
        zz_pContext_restore(&self.x)

##    def zz_p(self,v = None):
##        from ntl_lzz_p import ntl_zz_p
##        return ntl_zz_p(v,modulus=self)

##    def zz_pX(self,v = None):
##        from ntl_lzz_pX import ntl_zz_pX
##        return ntl_zz_pX(v,modulus=self)

def ntl_zz_pContext( long v ):
    try:
        return zz_pContextDict[repr(v)]
    except KeyError:
        return ntl_zz_pContext_class(v)
