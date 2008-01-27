"""
Power-sum symmetric functions
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
import sfa, multiplicative, classical

class SymmetricFunctionAlgebra_power(multiplicative.SymmetricFunctionAlgebra_multiplicative):
    def __init__(self, R):
        """
        TESTS:
            sage: p = SFAPower(QQ)
            sage: p == loads(dumps(p))
            True
        """
        classical.SymmetricFunctionAlgebra_classical.__init__(self, R, "power", SymmetricFunctionAlgebraElement_power, 'p')


    def dual_basis(self, scalar=None, prefix=None):
        """
        Returns the dual basis of the power-sum basis with
        respect to the scalar product scalar.  If scalar is None,
        then the standard scalar product for the classical
        symmetric functions is used.

        EXAMPLES:

        """
        if scalar is None:
            scalar = sfa.zee

        return dual.SymmetricFunctionAlgebra_dual(self, scalar, prefix=prefix)

class SymmetricFunctionAlgebraElement_power(classical.SymmetricFunctionAlgebraElement_classical):
    def frobenius(self):
        """
        Returns the image of self under the Frobenius / omega automorphism.

        EXAMPLES:
            sage: p = SFAPower(QQ)
            sage: a = p([2,1]); a
            p[2, 1]
            sage: a.frobenius()
            -p[2, 1]
            sage: a.omega()
            -p[2, 1]
        """
        parent = self.parent()
        base_ring = parent.base_ring()
        z = {}
        mcs = self.monomial_coefficients()
        for part in mcs:
            z[part] = (-1)**(sum(part)-len(part))*mcs[part]
        res = parent(0)
        res._monomial_coefficients = z
        return res

    def scalar(self, x):
        """
        Returns the standard scalar product of self and x.

        Note that the power-sum symmetric functions are orthogonal
        under this scalar product.  The value of <p_lambda, p_lambda>
        is given by the size of the centralizer in S_n of a permutation
        of cycle type lambda.

        EXAMPLES:
            sage: p = SFAPower(QQ)
            sage: p4 = Partitions(4)
            sage: matrix([ [p(a).scalar(p(b)) for a in p4] for b in p4])
            [ 4  0  0  0  0]
            [ 0  3  0  0  0]
            [ 0  0  8  0  0]
            [ 0  0  0  4  0]
            [ 0  0  0  0 24]
        """

        R = self.parent().base_ring()

        if self.parent() != x.parent():
            try:
                x = self.parent()( x )
            except:
                raise TypeError, "cannot compute the scalar product of self and x (= %s)"%x

        if len(self) < len(x):
            smaller = self
            greater = x
        else:
            smaller = x
            greater = self

        res = R(0)
        smcs = smaller._monomial_coefficients
        gmcs = greater._monomial_coefficients
        for s_part in smcs :
            if s_part in gmcs:
                res += smcs[s_part]*gmcs[s_part]*sfa.zee(s_part)

        return res
