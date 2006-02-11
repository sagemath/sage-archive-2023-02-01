"""
Constructor for Jacobian of a hyperelliptic curve
"""

#*****************************************************************************
#  Copyright (C) 2006 David Kohel <kohel@maths.usyd.edu>
#  Distributed under the terms of the GNU General Public License (GPL)
#                  http://www.gnu.org/licenses/
#*****************************************************************************

from hyperelliptic_generic import is_HyperellipticCurve
import jacobian_g2
import jacobian_generic

def Jacobian(C):
    if not is_HyperellipticCurve(C):
        raise TypeError, "Argument C (= %s) must be a hyperelliptic curve."
    if C.genus() == 2:
        return jacobian_g2.HyperellipticJacobian_g2(C)
    else:
        return jacobian_generic.HyperellipticJacobian_generic(C)
