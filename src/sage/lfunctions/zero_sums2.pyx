r"""
Class file for computing sums over zeros of motivic L-functions.
All computations are done to double precision.

AUTHORS:

- Simon Spicer (2014-06): first version

"""

##############################################################################
#       Copyright (C) 2014 Simon Spicer <mlungu@uw.edu>
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
##############################################################################

#from sage.structure.sage_object import SageObject
#from scipy.special import erfcx, spence, psi
#from sage.rings.integer_ring import ZZ
#from sage.rings.real_double import RDF
#from sage.rings.complex_double import CDF
#from sage.rings.infinity import PlusInfinity
#from sage.rings.arith import prime_powers
#from sage.functions.log import log, exp
#from sage.functions.other import real, imag
#from sage.symbolic.constants import pi, euler_gamma
#from sage.libs.pari.all import pari

cdef class TestClass:
    r"""
    Test Class
    """

    cdef _E
    cdef _N

    def __init__(self,E,N=None):
        r"""
        Initializes self.

        INPUT:

        - ``E`` -- An elliptic curve defined over the rational numbers

        - ``N`` -- (default: None) If not None, a positive integer equal to
          the conductor of E. This is passable so that rank estimation
          can be done for curves whose (large) conductor has been precomputed.

        EXAMPLES:

            sage: from sage.lfunctions.zero_sums import LFunctionZeroSum_EllipticCurve
            sage: E = EllipticCurve([1,0,0,3,-4])
            sage: Z = LFunctionZeroSum_EllipticCurve(E); Z
            Zero sum estimator for L-function attached to Elliptic Curve defined by y^2 + x*y = x^3 + 3*x - 4 over Rational Field
            sage: E = EllipticCurve('5077a')
            sage: Z = LFunctionZeroSum_EllipticCurve(E); Z
            Zero sum estimator for L-function attached to Elliptic Curve defined by y^2 + y = x^3 - 7*x + 6 over Rational Field
        """
        self._E = E
        if N is not None:
            self._N = N
        else:
            self._N = E.conductor()

    def __repr__(self):
        r"""
        Representation of self.

        EXAMPLES::

            sage: Z = LFunctionZeroSum(EllipticCurve('37a')); Z
            Zero sum estimator for L-function attached to Elliptic Curve defined by y^2 + y = x^3 - x over Rational Field
        """
        s = "Zero sum estimator for L-function attached to "
        return s+str(self._E)

    def elliptic_curve(self):
        """
        Return the elliptic curve associated with self.

        EXAMPLES::

            sage: E = EllipticCurve([23,100])
            sage: Z = LFunctionZeroSum(E)
            sage: Z.elliptic_curve()
            Elliptic Curve defined by y^2 = x^3 + 23*x + 100 over Rational Field

        """
        return self._E
