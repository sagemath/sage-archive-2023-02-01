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
from sage.rings.all import O

_KeyboardInterrupt = KeyboardInterrupt   # bug in Python requires this???

magma = None
def init():
    """
    Start up Magma (if necessary) and load the relevant files
    """
    global magma
    if magma is None:
        magma = Magma(script_subdirectory='ell_padic')
        magma('Attach("kedlaya.m");')
        magma('Attach("padic_height.m");')
        magma('Attach("gl2.m");')
        magma('Attach("myl.m");')
        magma('Attach("myl.m");')
        magma('Attach("shark.m");')

def padic_height(E, p, point, prec=20):
    """
    INPUT:
        E -- five-tuple of integers that define a Weierstrass equation
        p -- a prime number
        point -- point on E
        prec -- precision parameter
    OUTPUT:
        p-adic number
    """
    # TODO: add type checking for better error messages.

    init()
    global magma
    cmd = "E := EllipticCurve(%s); r := height_function(E,%s,%s)(E!%s); v:=Valuation(r); print [v, Integers()!(r*%s^(-v))];"%(
        E, p, prec, list(point), p)
    bad = False
    try:
        x = magma(cmd)
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
    init()
    points = [list(x) for x in points]
    cmd = "E := EllipticCurve(%s); r := regulator(E,%s,%s,%s); v:=Valuation(r); print [v, Integers()!(r*%s^(-v))];"%(
        E, p, prec, points, p)
    global magma
    bad = False
    try:
        x = magma(cmd)
    except KeyboardInterrupt:
        bad = True
        magma = None
        raise KeyboardInterrupt

    if not bad:
        try:
            v, n = eval(x)
            if x.find("error") != -1:
                bad = True
        except SyntaxError:
            bad = True

    if bad or not (isinstance(v, (int, long)) and isinstance(n, (int, long))):
        raise RuntimeError, "error computing p-adic regulator: %s"%x

    z = (n + O(p**(prec)))
    if v < 0:
        z = z / (p**(-v))
    elif v > 0:
        z = z * (p**v)
    return z
