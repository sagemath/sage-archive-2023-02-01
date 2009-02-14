"""nodoctest
Matrices over a domain
"""

########################################################################
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
########################################################################


cimport matrix
import  matrix

import sage.structure.sequence


cdef class Matrix_domain_dense(matrix.Matrix):


    def charpoly(self, *args, **kwds):
        r"""
        Return the characteristic polynomial of self, as a polynomial
        over the base ring.

        ALGORITHM: Use \code{self.matrix_over_field()} to obtain the
        corresponding matrix over a field, then call the
        \code{charpoly} function on it.

        EXAMPLES:
        First a matrix over $\Z$:
            sage: A = MatrixSpace(IntegerRing(),2)( [[1,2], [3,4]] )
            sage: f = A.charpoly('x')
            sage: f
            x^2 - 5*x - 2
            sage: f.parent()
            Univariate Polynomial Ring in x over Integer Ring

        We compute the characteristic polynomial of a matrix over
        the polynomial ring $\Z[a]$:
            sage: R = PolynomialRing(IntegerRing(),'a'); a = R.gen()
            sage: M = MatrixSpace(R,2)([[a,1], [a,a+1]])
            sage: M
            [    a     1]
            [    a a + 1]
            sage: f = M.charpoly('x')
            sage: f
            x^2 + (-2*a - 1)*x + a^2
            sage: f.parent()
            Univariate Polynomial Ring in x over Univariate Polynomial Ring in a over Integer Ring
            sage: M.trace()
            2*a + 1
            sage: M.determinant()
            a^2

        We compute the characteristic polynomial of a matrix over the
        multi-variate polynomial ring $\Z[x,y]$:
            sage: R = PolynomialRing(IntegerRing(),2); x,y = R.gens()
            sage: A = MatrixSpace(R,2)([x, y, x^2, y^2])
            sage: f = A.charpoly('x')
            sage: f
            x^2 + (-1*x1^2 - x0)*x + x0*x1^2 - x0^2*x1

        It's a little difficult to distinguish the variables.  To fix this,
        we rename the indeterminate $Z$:
            sage: f.parent()._assign_names("Z")
            sage: f
            Z^2 + (-1*x1^2 - x0)*Z + x0*x1^2 - x0^2*x1

        We can pass parameters in, which are passed on to the charpoly
        function for matrices over a field.
            sage: A = 1000*MatrixSpace(ZZ,10)(range(100))
            sage: A.charpoly(bound=2)
            x^10 + 14707*x^9 - 21509*x^8
            sage: A = 1000*MatrixSpace(ZZ,10)(range(100))
            sage: A.charpoly('x')
            x^10 - 495000*x^9 - 8250000000*x^8
        """
        f = self.matrix_over_field(copy=True).charpoly(*args, **kwds)
        return f.base_extend(self.base_ring())

    def determinant(self):
        """
        Return the determinant of this matrix.

        INPUT:
            -- a square matrix

        ALGORITHM: Find the characteristic polynomial and take its
        constant term (up to sign).

        EXAMPLES:
        We create a matrix over $\Z[x,y]$ and compute its determinant.
            sage: R = PolynomialRing(IntegerRing(),2); x,y = R.gens()
            sage: A = MatrixSpace(R,2)([x, y, x**2, y**2])
            sage: A.determinant()
            x0*x1^2 - x0^2*x1
        """
        if self.nrows() != self.ncols():
            raise ArithmeticError, "Matrix must be square, but is %sx%s"%(
                self.nrows(), self.ncols())
        # Use stupid slow but completely general method.
        d = (-1)**self.nrows() * self.charpoly('x')[0]
        return self.base_ring()(d)

    def is_invertible(self):
        r"""
        Return True if this matrix is invertible.

        EXAMPLES:
        The following matrix is invertible over $\Q$ but not over $\Z$.
            sage: A = MatrixSpace(IntegerRing(), 2)(range(4))
            sage: A.is_invertible()
            False
            sage: A.matrix_over_field().is_invertible()
            True

        The inverse function is a constructor for matrices over the
        fraction field, so it can work even if A is not invertible.
            sage: ~A   # inverse of A
            [-3/2  1/2]
            [   1    0]

        The next matrix is invertible over $\Z$.
            sage: A = MatrixSpace(IntegerRing(),2)([1,10,0,-1])
            sage: A.is_invertible()
            True
            sage: ~A                # compute the inverse
            [ 1 10]
            [ 0 -1]

        The following nontrivial matrix is invertible over $\Z[x]$.
            sage: R = PolynomialRing(IntegerRing())
            sage: x = R.gen()
            sage: A = MatrixSpace(R,2)([1,x,0,-1])
            sage: A.is_invertible()
            True
            sage: ~A
            [ 1  x]
            [ 0 -1]
        """
        return self.is_square() and self.determinant().is_unit()

    def __invert__(self):
        r"""
        Return this inverse of this matrix, as a matrix over the fraction field.

        Raises a \code{ZeroDivisionError} if the matrix has zero
        determinant, and raises an \code{ArithmeticError}, if the
        inverse doesn't exist because the matrix is nonsquare.

        EXAMPLES:
            sage: A = MatrixSpace(IntegerRing(), 2)([1,1,3,5])
            sage: ~A
            [ 5/2 -1/2]
            [-3/2  1/2]

        Even if the inverse lies in the base field, the result is still a matrix
        over the fraction field.
            sage: I = MatrixSpace(IntegerRing(),2)( 1 )  # identity matrix
            sage: ~I
            [1 0]
            [0 1]
            sage: (~I).parent()
            Full MatrixSpace of 2 by 2 dense matrices over Rational Field

        TESTS:
            sage: MatrixSpace(IntegerRing(), 0)().inverse()
            []
        """
        return ~self.matrix_over_field()


    def matrix_over_field(self, copy=False):
        """
        Return this matrix, but with entries viewed as elements
        of the fraction field of the base ring.

        EXAMPLES:
            sage: A = MatrixSpace(IntegerRing(),2)([1,2,3,4])
            sage: B = A.matrix_over_field()
            sage: B
            [1 2]
            [3 4]
            sage: B.parent()
            Full MatrixSpace of 2 by 2 dense matrices over Rational Field
        """
        return self.change_ring(self.base_ring().fraction_field(), copy=copy)

    def _singular_(self, singular=None):
        """
        Tries to coerce this matrix to a singular matrix.
        """
        if singular is None:
            from sage.interfaces.all import singular as singular_default
            singular = singular_default
        try:
            self.base_ring()._singular_(singular)
        except (NotImplementedError, AttributeError):
            raise TypeError, "Cannot coerce to Singular"

        return singular.matrix(self.nrows(),self.ncols(),singular(self.list()))


