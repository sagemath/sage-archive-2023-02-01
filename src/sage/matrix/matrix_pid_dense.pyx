"""nodoctest
Matrices over a PID
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


cimport matrix_domain_dense
import matrix_domain_dense

import sage.modules.free_module


cdef class Matrix_pid_dense(matrix_domain_dense.Matrix_domain_dense):

    def column_module(self):
        """
        Return the free module over the base ring spanned by the
        columns of this matrix.
        """
        return self.transpose().row_space()

    def decomposition(self, is_diagonalizable=False, dual=False):
        """
        Returns the decomposition of the free module on which this
        matrix acts from the right, along with whether this matrix
        acts irreducibly on each factor.  The factors are guaranteed
        to be sorted in the same way as the corresponding factors of
        the characteristic polynomial.

        Let A be the matrix acting from the on the vector space V of
        column vectors.  Assume that A is square.  This function
        computes maximal subspaces W_1, ..., W_n corresponding to
        Galois conjugacy classes of eigenvalues of A.  More precisely,
        let f(X) be the characteristic polynomial of A.  This function
        computes the subspace $W_i = ker(g_(A)^n)$, where g_i(X) is an
        irreducible factor of f(X) and g_i(X) exactly divides f(X).
        If the optional parameter is_diagonalizable is True, then we
        let W_i = ker(g(A)), since then we know that ker(g(A)) =
        $ker(g(A)^n)$.

        If dual is True, also returns the corresponding decomposition
        of V under the action of the transpose of A.  The factors are
        guarenteed to correspond.

        OUTPUT:
            list -- list of pairs (V,t), where V is a vector spaces
                    and t is a bool, and t is True exactly when the
                    charpoly of self on V is irreducible.

            (optional) list -- list of pairs (W,t), where W is a vector
                    space and t is a bool, and t is True exactly
                    when the charpoly of the transpose of self on W
                    is irreducible.

        EXAMPLES:
            sage: MS1 = MatrixSpace(ZZ,4)
            sage: MS2 = MatrixSpace(QQ,6)
            sage: A = MS1.matrix([3,4,5,6,7,3,8,10,14,5,6,7,2,2,10,9])
            sage: B = MS2(range(36))
            sage: B*11   # random output
            [-11  22 -11 -11 -11 -11]
            [ 11 -22 -11 -22  11  11]
            [-11 -11 -11 -22 -22 -11]
            [-22  22 -22  22 -11  11]
            [ 22 -11  11 -22  11  22]
            [ 11  11  11 -22  22  22]
            sage: decomposition(A)
            [(Ambient free module of rank 4 over the principal ideal domain Integer Ring, True)]
            sage: decomposition(B)
            [(Vector space of degree 6 and dimension 4 over Rational Field
            Basis matrix:
            [ 1  0  0  0 -5  4]
            [ 0  1  0  0 -4  3]
            [ 0  0  1  0 -3  2]
            [ 0  0  0  1 -2  1],
              False),
             (Vector space of degree 6 and dimension 2 over Rational Field
            Basis matrix:
            [ 1  0 -1 -2 -3 -4]
            [ 0  1  2  3  4  5],
              True)]
        """
        if not self.is_square():
            raise ArithmeticError, "self must be a square matrix"

        if self.nrows() == 0:
            return []

        f = self.charpoly('x')
        E = []

        # Idea: For optimization, could compute powers of self
        #       up to max degree of any factor.  Then get g(self)
        #       by taking a linear combination.   ??????

        if dual:
            Edual = []
        F = f.factor()
        if len(F) == 1:
            V = sage.modules.free_module.FreeModule(
                              self.base_ring(), self.nrows())
            m = F[0][1]
            if dual:
                return [(V,bool(m==1))], [(V,bool(m==1))]
            else:
                return [(V,bool(m==1))]
        F.sort()
        for g, m in f.factor():
            if is_diagonalizable:
                B = g(self)
            else:
                B = g(self) ** m
            E.append((B.kernel(), bool(m==1)))
            if dual:
                Edual.append((B.transpose().kernel(), bool(m==1)))
        if dual:
            return E, Edual
        return E

    def echelon_form(self, include_zero_rows=True):
        """
        Return the echelon form of this matrix over the integers.

        This is a matrix over the base ring (a PID) which is, \emph{by
        definition}, what is also called the Hermite normal form.
        """
        raise NotImplementedError


    def image(self):
        """
        Return the image of the homomorphism on rows defined by this matrix.

        EXAMPLES:
            sage: MS1 = MatrixSpace(ZZ,4)
            sage: MS2 = MatrixSpace(QQ,6)
            sage: A = MS1.matrix([3,4,5,6,7,3,8,10,14,5,6,7,2,2,10,9])
            sage: B = MS2.random_element()

            sage: image(A)
            Free module of degree 4 and rank 4 over Integer Ring
            Echelon basis matrix:
            [  1   0   0 426]
            [  0   1   0 518]
            [  0   0   1 293]
            [  0   0   0 687]

            sage: image(B) == B.row_module()
            True
        """
        return self.row_module()

    def row_module(self):
        """
        Return the free module over the base ring spanned by the rows
        of self.

        EXAMPLES:
            sage: A = MatrixSpace(IntegerRing(), 2)([1,2,3,4])
            sage: A.row_module()
            Free module of degree 2 and rank 2 over Integer Ring
            Echelon basis matrix:
            [1 0]
            [0 2]
        """
        M = sage.modules.free_module.FreeModule(self.base_ring(), self.ncols())
        return M.span(self.rows())

    def kernel_on(self, V, poly=None, check=False):
        """
        Return the kernel of self restricted to the invariant subspace V.
        The result is a vector subspace of V, which is also a subspace
        of the ambient space.

        INPUT:
            V -- vector subspace
            check -- (optional) default: False
            poly -- (optional) default: None; if not None, compute instead
                    the kernel of poly(self) on V.

        OUTPUT:
            a subspace

        WARNING: This function does \emph{not} check that V is in fact
        invariant under self, unless check is True (not the default).
        """
        A = self.restrict(V, check=check)
        if not poly is None:
            A = poly(A)
        W = A.kernel()
        if V.is_ambient():
            return W
        else:
            A = V.basis_matrix()
            B = W.basis_matrix()
            C = B*A
            return C.row_module()
