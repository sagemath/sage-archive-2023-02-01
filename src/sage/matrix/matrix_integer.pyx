import matrix_pid
cimport matrix_pid

import sage.rings.integer_ring

def _parimatrix_to_strlist(A):
    s = str(A)
    s = s.replace('Mat(','').replace(')','')
    s = s.replace(';',',').replace(' ','')
    s = s.replace(",", "','")
    s = s.replace("[", "['")
    s = s.replace("]", "']")
    return eval(s)

def _parimatrix_to_reversed_strlist(A):
    s = str(A)
    if s.find('Mat') != -1:
        return _parimatrix_to_strlist(A)
    s = s.replace('[','').replace(']','').replace(' ','')
    v = s.split(';')
    v.reverse()
    s = "['" + (','.join(v)) + "']"
    s = s.replace(",", "','")
    return eval(s)

def convert_parimatrix(z):
    n = z.ncols();
    r = []
    for i from 0 <= i < n:
        r.append(n-i)
    z = z.vecextract(r)
    z = z.mattranspose()
    n = z.ncols();
    r = []
    for i from 0 <= i < n:
        r.append(n-i)
    z = z.vecextract(r)
    return _parimatrix_to_strlist(z)


cdef class Matrix_integer(matrix_pid.Matrix_pid):

    def __init__(self, parent):
        matrix_pid.Matrix_pid.__init__(self, parent)

    def echelon_form(self, include_zero_rows=True):
        r"""
        Return the echelon form of this matrix over the integers.

        EXAMPLES:
            sage: A = MatrixSpace(IntegerRing(),2)([1,2,3,4])
            sage: A.echelon_form()
            [1 0]
            [0 2]

            sage: A = MatrixSpace(IntegerRing(),5)(range(25))
            sage: A.echelon_form()
            [  5   0  -5 -10 -15]
            [  0   1   2   3   4]
            [  0   0   0   0   0]
            [  0   0   0   0   0]
            [  0   0   0   0   0]
        """
        if self.nrows() == 0 or self.ncols() == 0:
            self.__rank = 0
            return self
        # The following complicated sequence of column reversals
        # and transposes is needed since PARI's Hermite Normal Form
        # does column operations instead of row operations.
        n = self.ncols()
        r = []
        for i from 0 <= i < n:
            r.append(n-i)
        v = self._pari_()
        v = v.vecextract(r) # this reverses the order of columns
        v = v.mattranspose()
        w = v.mathnf(1)

        cdef Matrix_integer H_m
        H = convert_parimatrix(w[0])
        if self.ncols() == 1:
            H = [H]

        # We can do a 'fast' change of the above into a list of ints,
        # since we know all entries are ints:
        (nr,nc) = (self.nrows(), self.ncols())
        num_missing_rows = (nr*nc - len(H)) / nc
        rank = nr - num_missing_rows
        if include_zero_rows:
            H = H + ['0']*(num_missing_rows*nc)
            H_m = self.new_matrix(nrows=nr, ncols=nc, entries=H, coerce_entries=True)
        else:
            H_m = self.new_matrix(nrows=rank, ncols=nc, entries=H, coerce_entries=True)
        H_m.__rank = rank
        H_m.set_immutable()
        return H_m

    def rank(self):
        """
        Return the rank of self, which is the rank of the space
        spanned by the rows of self.
        """
        if self.__rank is not None:
            return self.__rank
        else:
            r = self.change_ring(rational_field.RationalField()).rank()
            if self.is_immutable():
                self.__rank = r
            return r

    def elementary_divisors(self):
        """
        Return the elementary divisors of self, in order.

        The elementary divisors are the invariants of the finite
        abelian group that is the cokernel of this matrix.  They are
        ordered in reverse by divisibility.

        INPUT:
            matrix
        OUTPUT:
            list of int's

        EXAMPLES:
            sage: A = MatrixSpace(IntegerRing(), 3)(range(9))
            sage: A.elementary_divisors()
            [0, 3, 1]
            sage: C = MatrixSpace(ZZ,4)([3,4,5,6,7,3,8,10,14,5,6,7,2,2,10,9])
            sage: C.elementary_divisors()
            [687, 1, 1, 1]

        SEE ALSO: smith_form
        """
        try:
            return self.__elementary_divisors
        except AttributeError:
            if self.nrows() == 0 or self.ncols() == 0:
                return []
            d = self._pari_().matsnf(0).python()
            if self.is_immutable():
                self.__elementary_divisors = d
            return d

    def smith_form(self, transformation=False):
        """
        Returns matrices S, U, and V such that S = U*self*V, and S
        is in Smith normal form.  Thus S is diagonal with diagonal
        entries the ordered elementary divisors of S.

        The elementary divisors are the invariants of the finite
        abelian group that is the cokernel of this matrix.  They are
        ordered in reverse by divisibility.

        EXAMPLES:
            sage: A = MatrixSpace(IntegerRing(), 3)(range(9))
            sage: D, U, V = A.smith_form()
            sage: D
            [0 0 0]
            [0 3 0]
            [0 0 1]
            sage: U
            [-1  2 -1]
            [ 0 -1  1]
            [ 0  1  0]
            sage: V
            [ 1  4 -1]
            [-2 -3  1]
            [ 1  0  0]
            sage: U*A*V
            [0 0 0]
            [0 3 0]
            [0 0 1]

        It also makes sense for nonsquare matrices:

            sage: A = Matrix(ZZ,3,2,range(6))
            sage: D, U, V = A.smith_form()
            sage: D
            [0 0]
            [2 0]
            [0 1]
            sage: U
            [-1  2 -1]
            [ 0 -1  1]
            [ 0  1  0]
            sage: V
            [ 3 -1]
            [-2  1]
            sage: U * A * V
            [0 0]
            [2 0]
            [0 1]

        SEE ALSO: elementary_divisors
        """
        v = self._pari_().matsnf(1).python()
        print "v =", v
        D = self.matrix_space()(v[2])
        U = self.matrix_space(ncols = self.nrows())(v[0])
        V = self.matrix_space(nrows = self.ncols())(v[1])
        return D, U, V

    def frobenius(self,flag=0):
        """
        Return the Frobenius form (rational canonical form) of this matrix.

        If flag is 1, return only the elementary divisors.  If flag is
        2, return a two-components vector [F,B] where F is the
        Frobenius form and B is the basis change so that $M=B^{-1}FB$.

        INPUT:
           flag -- 0 (default), 1 or 2 as described above

        ALGORITHM: uses pari's matfrobenius()

        EXAMPLE:
           sage: A = MatrixSpace(IntegerRing(), 3)(range(9))
           sage: A.frobenius(0)
           [ 0  0  0]
           [ 1  0 18]
           [ 0  1 12]
           sage: A.frobenius(1)
           [x^3 - 12*x^2 - 18*x]
           sage: A.frobenius(2)
           ([ 0  0  0]
           [ 1  0 18]
           [ 0  1 12],
           [    -1      2     -1]
           [     0  23/15 -14/15]
           [     0  -2/15   1/15])

        AUTHOR:
           -- 2006-04-02: Martin Albrecht

        TODO:
           -- move this to work for more general matrices than just over Z.
              This will require fixing how PARI polynomials are coerced
              to SAGE polynomials.
        """
        if self.nrows()!=self.ncols():
            raise ArithmeticError, "frobenius matrix of non-square matrix not defined."

        v = self._pari_().matfrobenius(flag)
        if flag==0:
            return self.matrix_space()(v.python())
        elif flag==1:
            r = polynomial_ring.PolynomialRing(self.base_ring())
            #TODO: this should be handled in PolynomialRing not here
            retr = []
            for x in v:
                retr.append(eval(str(x).replace("^","**"), {}, r.gens_dict()))
            return retr
        elif flag==2:
            F = matrix_space.MatrixSpace(rational_field.RationalField(),
                                         self.nrows())(v[0].python())
            B = matrix_space.MatrixSpace(rational_field.RationalField(),
                                         self.nrows())(v[1].python())
            return F,B

    def kernel(self, LLL=False):
        r"""
        Return the kernel of this matrix, as a module over the integers.

        INPUT:
           LLL -- bool (default: False); if True the basis is an LLL
                  reduced basis; otherwise, it is an echelon basis.

        EXAMPLES:
            sage: M = MatrixSpace(IntegerRing(),4,2)(range(8))
            sage: M.kernel()
            Free module of degree 4 and rank 2 over Integer Ring
            Echelon basis matrix:
            [ 1  0 -3  2]
            [ 0  1 -2  1]
        """

        Z = self.base_ring()

        if self.nrows() == 0:    # from a 0 space
            M = sage.modules.free_module.FreeModule(Z, self.nrows())
            return M.zero_subspace()

        elif self.ncols() == 0:  # to a 0 space
            return sage.modules.free_module.FreeModule(Z, self.nrows())

        A = self._pari_().mattranspose()
        B = A.matkerint()


        n = self.nrows()
        Z = sage.rings.integer_ring.IntegerRing()
        M = sage.modules.free_module.FreeModule(Z, n)

        if B.ncols() == 0:
            return M.zero_submodule()

        # Now B is a basis for the LLL-reduced integer kernel as a
        # PARI object.  The basis vectors or B[0], ..., B[n-1],
        # where n is the dimension of the kernel.
        B = []
        for b in B:
          tmp = []
          for x in b:
              tmp.append(Z(x))
          B.append(M(tmp))
#        B = [M([Z(x) for x in b]) for b in B]
        if LLL:
            return M.span_of_basis(B)
        else:
            return M.span(B)

    def _adjoint(self):
        """assumes self is a square matrix (checked in adjoint)"""
        return self.parent()(self._pari_().matadjoint().python())

    def _lllgram(self):
        """assumes self is a square matrix (checked in lllgram)"""
        Z = sage.rings.integer_ring.IntegerRing()
        n = self.nrows()
        # pari does not like negative definite forms
        if n > 0 and self[0,0] < 0:
            self = -self
        # maybe should be /unimodular/ matrices ?
        MS = matrix_space.MatrixSpace(Z,n,n)
        try:
            U = MS(self._pari_().lllgramint().python())
        except (RuntimeError, ArithmeticError):
            raise ValueError, "not a definite matrix"
        # Fix last column so that det = +1
        if U.det() == -1:
            for i in range(n):
                U[i,n-1] = - U[i,n-1]
        return U

    def _ntl_(self):
        """
        ntl.mat_ZZ representation of self.

        \note{NTL only knows dense matrices, so if you provide a
        sparse matrix NTL will allocate memory for every zero entry.}
        """
        return mat_ZZ(self.nrows(),self.ncols(),self.list())

