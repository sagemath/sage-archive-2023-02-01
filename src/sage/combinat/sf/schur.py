"""
Schur symmetric functions
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
import classical, dual
import sage.libs.symmetrica.all as symmetrica
from sage.rings.all import ZZ, QQ, Integer

class SymmetricFunctionAlgebra_schur(classical.SymmetricFunctionAlgebra_classical):
    def __init__(self, R):
        """
        TESTS::

            sage: s = SFASchur(QQ)
            sage: s == loads(dumps(s))
            True
        """
        classical.SymmetricFunctionAlgebra_classical.__init__(self, R, "schur", 's')

    def is_schur_basis(self):
        """
        EXAMPLES::

            sage: s = SFASchur(QQ)
            sage: s.is_schur_basis()
            True
        """
        return True

    def dual_basis(self, scalar=None, scalar_name="",  prefix=None):
        """
        The dual basis to the Schur basis with respect to the standard
        scalar product is the Schur basis since it is self-dual.

        EXAMPLES::

            sage: s = SFASchur(QQ)
            sage: ds = s.dual_basis()
            sage: s is ds
            True
        """
        if scalar is None:
            return self
        else:
            return dual.SymmetricFunctionAlgebra_dual(self, scalar, scalar_name, prefix)

    def _multiply(self, left, right): # TODO: factor out this code for all bases (as is done for coercions)
        """
        TESTS::

            sage: s = SFASchur(QQ)
            sage: a = s([2,1]) + 1; a
            s[] + s[2, 1]
            sage: a^2   # indirect doctest
            s[] + 2*s[2, 1] + s[2, 2, 1, 1] + s[2, 2, 2] + s[3, 1, 1, 1] + 2*s[3, 2, 1] + s[3, 3] + s[4, 1, 1] + s[4, 2]

        ::

            sage: QQx.<x> = QQ[]
            sage: s = SFASchur(QQx)
            sage: a = x^2*s([2,1]) + 2*x; a
            2*x*s[] + x^2*s[2, 1]
            sage: a^2
            4*x^2*s[] + 4*x^3*s[2, 1] + x^4*s[2, 2, 1, 1] + x^4*s[2, 2, 2] + x^4*s[3, 1, 1, 1] + 2*x^4*s[3, 2, 1] + x^4*s[3, 3] + x^4*s[4, 1, 1] + x^4*s[4, 2]

        ::

        ::

            sage: 0*s([2,1])
            0
        """
        #Use symmetrica to do the multiplication
        A = left.parent()
        R = A.base_ring()

        if  R is ZZ or R is QQ:
            if left and right:
                return symmetrica.mult_schur_schur(left, right)
            else:
                return A._from_dict({})

        z_elt = {}
        for (left_m, left_c) in left._monomial_coefficients.iteritems():
            for (right_m, right_c) in right._monomial_coefficients.iteritems():
                c = left_c * right_c
                d = symmetrica.mult_schur_schur({left_m:Integer(1)}, {right_m:Integer(1)})._monomial_coefficients
                for m in d:
                    if m in z_elt:
                        z_elt[ m ] = z_elt[m] + c * d[m]
                    else:
                        z_elt[ m ] = c * d[m]
        return A._from_dict(z_elt)


    class Element(classical.SymmetricFunctionAlgebra_classical.Element):
        # TODO: fix indentation


      def __pow__(self, n):
          """
          Naive powering

          See ``Monoids.Element.__pow__`` and ``Monoids.Element._pow_naive``.

          EXAMPLES::

              sage: s = SFASchur(QQ[x])
              sage: len(s([2,1])^8) # long time (~ 4 s)
              1485
              sage: len(s([2,1])^9) # long time (~10 s)
              2876

          Binary exponentiation does not seem to bring any speedup for
          schur functions. This most likely is because of the
          explosion of the number of terms.

          #    sage: s = SFASchur(QQ); y = s([1])
          #    sage: n = 24
          #    sage: %timeit y**n    # using binary exponentiation
          #    10 loops, best of 3: 1.22 s per loop
          #    sage: %timeit prod(y for i in range(n))
          #    10 loops, best of 3: 1.06 s per loop

          With polynomial coefficients, this is actually much *slower*
          (although this should be profiled further; there seems to
          be an unreasonable number of polynomial multiplication involved,
          besides the fact that 1 * QQ[x].one() currently involves a
          polynomial multiplication)

          #    sage: sage: s = SFASchur(QQ[x])
          #    sage: y = s([2,1])
          #    sage: %timeit y**7
          #    10 loops, best of 3: 18.9 s per loop
          #    sage: %timeit y*y*y*y*y*y*y
          #    10 loops, best of 3: 1.73 s per loop

          Todo: do the same for the other non multiplicative bases?

          """
          return self._pow_naive(n)

      def omega(self):
        """
        Returns the image of self under the Frobenius / omega
        automorphism.

        EXAMPLES::

            sage: s = SFASchur(QQ)
            sage: s([2,1]).omega()
            s[2, 1]
            sage: s([2,1,1]).omega()
            s[3, 1]
        """
        conj = lambda part: part.conjugate()
        return self.map_support(conj)


      def scalar(self, x):
        """
        Returns the standard scalar product between self and x.

        Note that the Schur functions are self-dual with respect to this
        scalar product. They are also lower-triangularly related to the
        monomial symmetric functions with respect to this scalar product.

        EXAMPLES::

            sage: s = SFASchur(ZZ)
            sage: a = s([2,1])
            sage: b = s([1,1,1])
            sage: c = 2*s([1,1,1])
            sage: d = a + b
            sage: a.scalar(a)
            1
            sage: b.scalar(b)
            1
            sage: b.scalar(a)
            0
            sage: b.scalar(c)
            2
            sage: c.scalar(c)
            4
            sage: d.scalar(a)
            1
            sage: d.scalar(b)
            1
            sage: d.scalar(c)
            2

        ::

            sage: m = SFAMonomial(ZZ)
            sage: p4 = Partitions(4)
            sage: l = [ [s(p).scalar(m(q)) for q in p4] for p in p4]
            sage: matrix(l)
            [ 1  0  0  0  0]
            [-1  1  0  0  0]
            [ 0 -1  1  0  0]
            [ 1 -1 -1  1  0]
            [-1  2  1 -3  1]
        """
        s = self.parent()
        R = s.base_ring()
        one = R(1)
        f = lambda p1, p2: one
        x = s(x)
        return s._apply_multi_module_morphism(self, x, f, orthogonal=True)

      def expand(self, n, alphabet='x'):
        """
        Expands the symmetric function as a symmetric polynomial in n
        variables.

        EXAMPLES::

            sage: s = SFASchur(QQ)
            sage: a = s([2,1])
            sage: a.expand(2)
            x0^2*x1 + x0*x1^2
            sage: a.expand(3)
            x0^2*x1 + x0*x1^2 + x0^2*x2 + 2*x0*x1*x2 + x1^2*x2 + x0*x2^2 + x1*x2^2
            sage: a.expand(4)
            x0^2*x1 + x0*x1^2 + x0^2*x2 + 2*x0*x1*x2 + x1^2*x2 + x0*x2^2 + x1*x2^2 + x0^2*x3 + 2*x0*x1*x3 + x1^2*x3 + 2*x0*x2*x3 + 2*x1*x2*x3 + x2^2*x3 + x0*x3^2 + x1*x3^2 + x2*x3^2
            sage: a.expand(2, alphabet='y')
            y0^2*y1 + y0*y1^2
            sage: a.expand(2, alphabet=['a','b'])
            a^2*b + a*b^2
            sage: s([1,1,1,1]).expand(3)
            0
        """
        condition = lambda part: len(part) > n
        return self._expand(condition, n, alphabet)

# Backward compatibility for unpickling
from sage.structure.sage_object import register_unpickle_override
register_unpickle_override('sage.combinat.sf.schur', 'SymmetricFunctionAlgebraElement_schur',  SymmetricFunctionAlgebra_schur.Element)
