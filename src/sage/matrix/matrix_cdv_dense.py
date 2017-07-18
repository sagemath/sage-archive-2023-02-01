from sage.matrix.matrix_generic_dense import Matrix_generic_dense
from sage.rings.infinity import Infinity

class Matrix_cdv_dense(Matrix_generic_dense):
    def _elementary_divisors_unimodular(self, transformation):
        n = self.nrows()
        m = self.ncols()
        S = self.parent()(self.list())
        R = self.base_ring()
        if transformation:
            from sage.matrix.special import identity_matrix
            left = identity_matrix(R,n)
            right = identity_matrix(R,m)
        else:
            left = right = None
        val = -Infinity
        divisors = [ ]
        for piv in range(min(n,m)):
            curval = Infinity
            for i in range(piv,n):
                for j in range(piv,m):
                    v = S[i,j].valuation()
                    if v < curval:
                        pivi = i; pivj = j
                        curval = v
                        if v == val: break
                else:
                    continue
                break
            if S[piv,piv] == 0:
                divisors += (min(n,m)-piv) * [ S[piv,piv] ]
                break
            if pivi > piv:
                S.swap_rows(pivi,piv)
                S.rescale_row(piv,-1)
                if transformation:
                    left.swap_rows(pivi,piv)
                    left.rescale_row(piv,-1)
            if pivj > piv:
                S.swap_columns(pivj,piv)
                S.rescale_col(piv,-1)
                if transformation:
                    right.swap_columns(pivj,piv)
                    right.rescale_col(piv,-1)
            divisors.append(S[piv,piv])
            inv = ~(S[piv,piv] >> curval)
            for i in range(piv+1,n):
                scalar = -inv * (S[i,piv] >> curval)
                scalar = scalar.lift_to_precision()
                S.add_multiple_of_row(i,piv,scalar,piv+1)
                if transformation:
                    left.add_multiple_of_row(i,piv,scalar)
            if transformation:
                for j in range(piv+1,m):
                    scalar = -inv * (S[piv,j] >> curval)
                    scalar = scalar.lift_to_precision()
                    right.add_multiple_of_column(j,piv,scalar)
        return divisors, left, right

    def elementary_divisors(self, transformation=False):
        """
        Return the elementary divisors of this matrix.

        INPUT:

        - ``transformation`` -- a boolean (default: False)
          Indicates whether the transformation matrices are returned

        NOTE:

        We recall that a Smith decomposition of a matrix `M`
        defined over a complete discrete valuation ring/field
        is a writing of the form `L*M*R = S` where:

        - `L` and `R` are invertible matrices in the ring of
          integers

        - the only non-vanishing entries of `S` are located on
          the diagonal (through `S` might be not a square matrix)

        - if `d_i` denotes the `(i,i)` entry of `S`, then `d_i`
          divides `d_{i+1}` for all `i`.

        A Smith decomposition is unique if the `d_i` are normalized
        so that they are all either `0` or a power of the distinguished
        uniformizer of the base ring.
        Normalized this way, they are called the elementary divisors of `M`.

        EXAMPLES::

            sage: A = Zp(5, prec=10, print_mode="digits")
            sage: M = matrix(A, 2, 2, [2, 7, 1, 6])
            sage: M.elementary_divisors()
            [...1, ...10]

            sage: divisors, L, R = M.elementary_divisors(transformation=True)
            sage: divisors
            [...1, ...10]
            sage: L
            [...2222222223             0]
            [ ...444444444          ...2]
            sage: R
            [         ...1 ...2222222214]
            [            0          ...1]

        This method works for rectangular matrices as well::

            sage: M = matrix(A, 3, 2, [2, 7, 1, 6, 3, 8])
            sage: divisors, L, R = M.elementary_divisors(transformation=True)
            sage: divisors
            [...1, ...10]
            sage: L
            [...2222222223             0             0]
            [ ...444444444          ...2             0]
            [...4444444443          ...1          ...1]
            sage: R
            [         ...1 ...2222222214]
            [            0          ...1]

        It may happen that elementary divisors cannot be determined due 
        to a lack of precision, in which case non-exact zeros are used::

            sage: M = matrix(R, 2, 2, [1,1,1,1])
            sage: M.elementary_divisors()
            [...1, ...]

        TESTS::

            sage: M = random_matrix(R, 4)
            sage: divisors, L, R = M.elementary_divisors(transformation=True)
            sage: L*M*R == diagonal_matrix(divisors)
            True
        """
        R = self.base_ring()
        divisors, left, right = self._elementary_divisors_unimodular(transformation)
        for i in range(len(divisors)):
            div = divisors[i]
            if div != 0:
                v = div.valuation()
                if transformation:
                    scalar = ~(div >> v)
                    left.rescale_row(i,scalar)
                divisors[i] = R(1) << v
        if transformation:
            return divisors, left, right
        else:
            return divisors

    def smith_form(self, transformation=True):
        """
        Return the Smith normal form of this matrix.

        INPUT:

        - ``transformation`` -- a boolean (default: True)
          Indicates whether the transformation matrices are returned

        NOTE:

        We recall that a Smith decomposition of a matrix `M`
        defined over a complete discrete valuation ring/field
        is a writing of the form `L*M*R = S` where:

        - `L` and `R` are invertible matrices in the ring of
          integers

        - the only non-vanishing entries of `S` are located on
          the diagonal (through `S` might be not a square matrix)

        - if `d_i` denotes the `(i,i)` entry of `D`, then `d_i`
          divides `d_{i+1}` for all `i`.

        A Smith decomposition is unique if the `d_i` are normalized
        so that they are all either `0` or a power of the distinguished
        uniformizer of the base ring.
        Normalized this way, the writing `L*M*R = S` is called the
        Smith normal form of `M`.

        EXAMPLES::

            sage: A = Zp(5, prec=10, print_mode="digits")
            sage: M = matrix(A, 2, 2, [2, 7, 1, 6])

            sage: S, L, R = M.smith_form()
            sage: S
            [ ...1     0]
            [    0 ...10]
            sage: L
            [...2222222223             0]
            [ ...444444444          ...2]
            sage: R
            [         ...1 ...2222222214]
            [            0          ...1]

        This method works for rectangular matrices as well::

            sage: M = matrix(A, 3, 2, [2, 7, 1, 6, 3, 8])
            sage: S, L, R = M.smith_form()
            [ ...1     0]
            [    0 ...10]
            [    0     0]
            sage: L
            [...2222222223             0             0]
            [ ...444444444          ...2             0]
            [...4444444443          ...1          ...1]
            sage: R
            [         ...1 ...2222222214]
            [            0          ...1]

        TESTS::

            sage: M = random_matrix(R, 4)
            sage: S, L, R = M.smith_form()
            sage: L*M*R == S
            True
        """
        R = self.base_ring()
        divisors, left, right = self._elementary_divisors_unimodular(transformation)
        S = self.parent()(0)
        for i in range(len(divisors)):
            div = divisors[i]
            if div == 0:
                S[i,i] = div
            else:
                v = div.valuation()
                if transformation:
                    scalar = ~(div >> v)
                    left.rescale_row(i,scalar)
                S[i,i] = R(1) << v
        if transformation:
            return S, left, right
        else:
            return S
        
    def determinant(self):
        """
        Return the determinant of this matrix.

        ALGORITHM:

        We compute the unimodular Smith normal form of the 
        input matrix `M`, that is we decompose `M` as a product
        `M = L*D*R` where `D` is diagonal and `L` are `R` are 
        unimodular matrices with coefficients in the ring of
        integers.
        Then, we compute `det(M)` as `det(D)` which is the
        product of its diagonal entries.

        This algorithm has a good behaviour regarding precision.
        (The result it outputs has optimal precision when each
        entry of `M` is given at the same absolute precision.)

        """
        if self.nrows() != self.ncols():
            raise ValueError("self must be a square matrix")
        divisors, _, _ = self._elementary_divisors_unimodular(False)
        det = self.base_ring()(1)
        for div in divisors:
            det *= div
        return det
