"""
Multiplicative symmetric functions
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
import classical
import sage.combinat.partition

class SymmetricFunctionAlgebra_multiplicative(classical.SymmetricFunctionAlgebra_classical):
    def _multiply_basis(self, left, right):
        """
        TESTS:
            sage: e = SFAElementary(QQ)
            sage: e([2,1])^2
            e[2, 2, 1, 1]

            sage: h = SFAHomogeneous(QQ)
            sage: h([2,1])^2
            h[2, 2, 1, 1]

            sage: p = SFAPower(QQ)
            sage: p([2,1])^2
            p[2, 2, 1, 1]

            sage: QQx.<x> = QQ[]
            sage: p = SFAPower(QQx) # indirect doctest
            sage: (x*p([2]))^2
            x^2*p[2, 2]
        """
        m = list(left)+list(right)
        m.sort(reverse=True)
        return sage.combinat.partition.Partition_class(m)


