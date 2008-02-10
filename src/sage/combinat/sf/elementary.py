"""
Elementary symmetric functions
"""
#*****************************************************************************
#       Copyright (C) 2007 Mike Hansen <mhansen@gmail.com>,
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
import sfa, multiplicative, classical, dual
from sage.rings.all import PolynomialRing
import sage.libs.symmetrica.all as symmetrica

###################################
#                                 #
# Elementary Symmetric Functions  #
#                                 #
###################################
class SymmetricFunctionAlgebra_elementary(multiplicative.SymmetricFunctionAlgebra_multiplicative):
    def __init__(self, R):
        """
        TESTS:
            sage: e = SFAElementary(QQ)
            sage: e == loads(dumps(e))
            True
        """
        classical.SymmetricFunctionAlgebra_classical.__init__(self, R, "elementary", SymmetricFunctionAlgebraElement_elementary, 'e')

    def dual_basis(self, scalar=None, prefix=None):
        """
        Returns the dual basis of the elementary basis with
        respect to the scalar product scalar.  If scalar is None,
        then the standard scalar product for the classical
        symmetric functions is used.

        EXAMPLES:

        """
        if scalar is None:
            scalar = zee

        return dual.SymmetricFunctionAlgebra_dual(self, scalar, prefix=prefix)


class SymmetricFunctionAlgebraElement_elementary(classical.SymmetricFunctionAlgebraElement_classical):
    def frobenius(self):
        """
        Returns the image of self under the Frobenius / omega automorphism.

        EXAMPLES:
            sage: e = SFAElementary(QQ)
            sage: a = e([2,1]); a
            e[2, 1]
            sage: a.frobenius()
            e[1, 1, 1] - e[2, 1]
            sage: a.omega()
            e[1, 1, 1] - e[2, 1]

            sage: h = SFAHomogeneous(QQ)
            sage: h(e([2,1]).omega())
            h[2, 1]
        """
        base_ring = self.parent().base_ring()
        h = sfa.SFAHomogeneous(base_ring)
        mcs = self.monomial_coefficients()
        res = h(0)
        res._monomial_coefficients = mcs
        return self.parent()(res)

    def expand(self, n, alphabet='x'):
        """
        Expands the symmetric function as a symmetric polynomial in n variables.

        EXAMPLES:
            sage: e = SFAElementary(QQ)
            sage: e([2,1]).expand(3)
            x0^2*x1 + x0*x1^2 + x0^2*x2 + 3*x0*x1*x2 + x1^2*x2 + x0*x2^2 + x1*x2^2
            sage: e([1,1,1]).expand(2)
            x0^3 + 3*x0^2*x1 + 3*x0*x1^2 + x1^3
            sage: e([3]).expand(2)
            0
            sage: e([2]).expand(3)
            x0*x1 + x0*x2 + x1*x2
        """
        e = eval('symmetrica.compute_' + str(classical.translate[self.parent().basis_name()]).lower() + '_with_alphabet')
        resPR = PolynomialRing(self.parent().base_ring(), n, alphabet)
        res = resPR(0)
        self_mc = self._monomial_coefficients
        for part in self_mc:
            if max(part) > n:
                continue
            res += self_mc[part] * e(part, n, alphabet)
        return res
