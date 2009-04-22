r"""
Elliptic curves over a general field

This module defines the class ``EllipticCurve_field``, based on
``EllipticCurve_generic``, for elliptic curves over general fields.

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
        r"""
        Verify that a coordinate triple satisfies the equations of
        this elliptic curve.

        INPUT:

        - ``v`` (list or tuple) -- homogeneous coordinates of a point in `\mathbb{P}^2`.

        .. warning::

           ``v[2]`` must be 0 or 1.

        OUTPUT:

        Nothing if the equations are satisfied; a TypeError is raised
        if not.

        EXAMPLES::

            sage: E = EllipticCurve('37a1')
            sage: E._check_satisfies_equations([0,-1,0])
            sage: E._check_satisfies_equations([0,1,0])
            sage: E._check_satisfies_equations([0,0,0])
            Traceback (most recent call last):
            ...
            TypeError: coordinates [0, 0, 0] do not define a point on Elliptic Curve defined by y^2 + y = x^3 - x over Rational Field
        """
        x, y, z = v
        if z == 0:
            if not (x == 0 and y != 0):
                self._error_bad_coords(v)
            return
        a1, a2, a3, a4, a6 = self.ainvs()
        if y**2 + (a1*x + a3)*y != x*(x*(x+a2)+a4)+a6:
#        if y**2 + (a[0]*x + a[2])*y != x*(x*(x+a[1])+a[3])+a[4]:
            self._error_bad_coords(v)


