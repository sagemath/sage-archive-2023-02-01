"""
Elliptic curves over a general field
"""

#*****************************************************************************
#       Copyright (C) 2006 William Stein <wstein@gmail.com>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#
#                  http://www.gnu.org/licenses/
#*****************************************************************************

import ell_generic

class EllipticCurve_field(ell_generic.EllipticCurve_generic):
    def _check_satisfies_equations(self, v):
        """
        Verify that the coordinates of v = [a,b,c] define a point on
        this scheme, or raise a TypeError.  Note that c is assumed
        to equal either 0 or 1.

        EXAMPLES:
            sage: E = EllipticCurve('37a1')
            sage: E._check_satisfies_equations([0,-1,0])
            sage: E._check_satisfies_equations([0,1,0])
            sage: E._check_satisfies_equations([0,0,0])
            Traceback (most recent call last):
            ...
            TypeError: coordinates [0, 0, 0] do not define a point on Elliptic Curve defined by y^2 + y = x^3 - x over Rational Field

        """
        if v[2] == 0:
            if not (v[0] == 0 and v[1] != 0):
                self._error_bad_coords(v)
            return
        x, y, a = v[0], v[1], self.ainvs()
        if y**2 + a[0]*x*y + a[2]*y != x**3 + a[1]*x**2 + a[3]*x + a[4]:
            self._error_bad_coords(v)


