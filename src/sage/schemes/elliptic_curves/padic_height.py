"""
Access to MAGMA code for computing p-adic height pairings.

(deprecated)
"""

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

from sage.interfaces.magma import Magma
from sage.rings.all import O, pAdicField

_KeyboardInterrupt = KeyboardInterrupt   # bug in Python requires this???

magma = None
def init():
    global magma
    if magma is None:
        magma = Magma(script_subdirectory='padic_height')
        magma.attach('kedlaya.m')
        magma.attach('padic_height.m')

def padic_eval(cmd, p, prec):
    init()
    global magma
    bad = False
    try:
        x = magma.eval(cmd)
    except _KeyboardInterrupt:
        bad = True
        magma = None
        raise _KeyboardInterrupt

    if not bad:
        try:
            v, n = eval(x)
            if x.find("error") != -1:
                bad = True
        except SyntaxError, KeyboardInterrupt:
            bad = True

    if bad or not (isinstance(v, (int, long)) and isinstance(n, (int, long))):
        magma = None
        raise RuntimeError, "error computing p-adic height: %s"%x

    z = (n + O(p**(prec)))
    if v < 0:
        z = z / (p**(-v))
    elif v > 0:
        z = z * (p**v)
    return z

def padic_height(E, p, point, prec=20):
    """
    INPUT:
        E -- five-tuple of integers that define a weierstrass equation
        p -- a prime number
        point -- point on E
        prec -- precision parameter
    OUTPUT:
        p-adic number
    """
    # TODO: add type checking for better error messages.
    cmd = "E := EllipticCurve(%s); r := height_function(E,%s,%s)(E!%s); v:=Valuation(r); print [v, Integers()!(r*%s^(-v))];"%(
                E, p, prec, list(point), p)
    return padic_eval(cmd, p, prec)

def padic_regulator(E, p, points, prec=20):
    """
    Returns the p-adic regulator, where points is a basis for the Mordell-Weil group
    modulo torsion.

    INPUT:
        E -- five-tuple of integers that define a weierstrass equation
        p -- a prime number
        points -- list of points on E
        prec -- precision parameter
    OUTPUT:
        p-adic number
    """
    points = [list(x) for x in points]
    cmd = "E := EllipticCurve(%s); r := regulator(E,%s,%s,%s); v:=Valuation(r); print [v, Integers()!(r*%s^(-v))];"%(
        E, p, prec, points, p)
    return padic_eval(cmd, p, prec)

def padic_E2(E, p, prec=20):
    """
    INPUT:
        E -- five-tuple of integers that define a weierstrass equation
        p -- a prime number
        prec -- precision parameter
    OUTPUT:
        p-adic number
    """
    cmd = "E := EllipticCurve(%s); r := E2(E,%s,%s); v:=Valuation(r); print [v, Integers()!(r*%s^(-v))];"%(
         E, p, prec, p)
    return padic_eval(cmd, p, prec)

def padic_E2_of_c4c6(c4, c6, p, prec=20):
    K = pAdicField(p)
    c4 = K(c4); c6 = K(c6)
    cmd = "r := E2_c4c6(%s,%s,%s,%s); v:=Valuation(r); print [v, Integers()!(r*%s^(-v))];"%(
         c4.lift(), c6.lift(), p, prec, p)
    return padic_eval(cmd, p, prec)
