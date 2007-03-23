"""
Elements of p-Adic Fields

AUTHOR:
    -- David Roe
    -- Genya Zaytman: documentation
    -- David Harvey: doctests
"""

#*****************************************************************************
#       Copyright (C) 2007 David Roe <roed@math.harvard.edu>
#                          William Stein <wstein@gmail.com>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#
#                  http://www.gnu.org/licenses/
#*****************************************************************************

#import sage.rings.padics.precision_error
import sage.rings.padics.padic_generic_element

class pAdicFieldGenericElement(sage.rings.padics.padic_generic_element.pAdicGenericElement):
    def __getitem__(self, n):
        r"""
        Returns the coefficient of $p^n$ in the series expansion of self, as an integer in the range $0$ to $p-1$.

        EXAMPLES:
            sage: R = Qp(7, 4, 'capped-rel', 'series'); a = R(1/21); a
            5*7^-1 + 4 + 4*7 + 4*7^2 + O(7^3)
            sage: a[-1]
            5
            sage: a[0]
            4
        """
        v = self.valuation()
        if n >= self.precision_absolute():
            raise sage.rings.padics.precision_error.PrecisionError, "element not known to enough precision."
        elif n < v:
            return 0
        else:
            return self.list()[n - v]
