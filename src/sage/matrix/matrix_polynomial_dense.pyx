"""
Dense matrices over univariate polynomials over fields

The implementation inherits from Matrix_generic_dense but some algorithms
are optimized for polynomial matrices.

AUTHORS:

- Kwankyu Lee (2016-12-15): initial version with code moved from other files.

- Johan Rosenkilde (2017-02-07): added weak_popov_form()

"""
#*****************************************************************************
#       Copyright (C) 2016 Kwankyu Lee <ekwankyu@gmail.com>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#  as published by the Free Software Foundation; either version 2 of
#  the License, or (at your option) any later version.
#                  http://www.gnu.org/licenses/
#*****************************************************************************

from sage.matrix.matrix_generic_dense cimport Matrix_generic_dense
from sage.matrix.matrix2 cimport Matrix
from sage.rings.integer_ring import ZZ
from sage.misc.superseded import deprecated_function_alias

cdef class Matrix_polynomial_dense(Matrix_generic_dense):
    r"""
    Dense matrix over a univariate polynomial ring over a field.

    For a field $\Bold{K}$, we consider matrices over the univariate
    polynomial ring $\Bold{K}[x]$.
    
    They are often used to represent bases of some $\Bold{K}[x]$-modules. In
    this context, there are two possible representations which are both
    commonly used in the literature.

    - Working column-wise: each column of the matrix is a vector in the basis;
      then, a $\\Bold{K}[x]$-submodule of $\\Bold{K}[x]^{m}$ of rank $n$ is
      represented by an $m \\times n$ matrix, whose columns span the module
      (via $\\Bold{K}[x]$-linear combinations). This matrix has full rank,
      and $n \\leq m$.

    - Working row-wise: each row of the matrix is a vector in the basis; then,
      a $\\Bold{K}[x]$-submodule of $\\Bold{K}[x]^{n}$ of rank $m$ is
      represented by an $m \\times n$ matrix, whose rows span the module (via
      $\\Bold{K}[x]$-linear combinations). This matrix has full rank, and $m
      \\leq n$.

    For the rest of this class description, we assume that one is working
    row-wise. For a given such module, all its bases are equivalent under
    left-multiplication by a unimodular matrix, that is, a matrix which has
    determinant in $\Bold{K}\setminus\{0\}$.

    There are bases which are called reduced or minimal: their rows have the
    minimal degree possible among all bases of this module. The degree of a row
    is the maximum of the degrees of the entries of the row. An equivalent
    condition is that the leading matrix of this basis has full rank (see the
    description of :meth:`leading_matrix`). There is a unique minimal basis,
    called the Popov basis of the module, which satisfies some additional
    normalization condition (see the description of :meth:`row_degree`).

    These notions can be extended via a more general degree measure, involving
    a tuple of integers which is called shift and acts as column degree shifts
    in the definition of row degree. Precisely, for given $s_1,\ldots,s_n \in
    \ZZ$ and a row vector $[p_1 \; \cdots \; p_n] \in \Bold{K}[x]^{1 \times
    n}$, its shifted row degree is the maximum of $\deg(p_j) + s_j$ for $1 \leq
    j \leq n$. Then, minimal bases and Popov bases are defined similarly, with
    respect to this notion of degree.

    Another important canonical basis is the Hermite basis, which is a lower
    triangular matrix satisfying a normalization condition similar to that for
    the Popov basis. In fact, if $d$ is the largest degree appearing in the
    Hermite basis, then the Hermite basis coincide with the shifted Popov basis
    with the shift $(0,d,2d,\ldots,(n-1)d)$.
    """

    def _check_shift_dimension(self, shifts=None, row_wise=True):
        r"""
        Return a boolean indicating if the ``shifts`` argument has the
        right length.

        For an $m \times n$ polynomial matrix, if working row-wise then
        ``shifts`` should have $n$ entries; if working column-wise, it should
        have $m$ entries.

        INPUT:

        - ``shifts`` -- (optional, default: ``None``) list of integers;
          ``None`` is interpreted as ``shifts=[0,...,0]``.

        - ``row_wise`` -- (optional, default: ``True``) boolean, if ``True``
          then shifts apply to the columns of the matrix and otherwise to its
          rows (see the class description for more details).

        EXAMPLES::

            sage: pR.<x> = GF(7)[]
            sage: M = Matrix( pR, [[3*x+1, 0, 1], [x^3+3, 0, 0]])
            sage: M._check_shift_dimension()

            sage: M._check_shift_dimension(shifts=[1,3,2])

            sage: M._check_shift_dimension(shifts=[1,3,2], row_wise=False)
            Traceback (most recent call last):
            ...
            ValueError: Shifts length should be the row dimension.
        """
        if shifts!=None and (not row_wise) and len(shifts) != self.nrows():
            raise ValueError('Shifts length should be the row dimension.')
        if shifts!=None and (row_wise and len(shifts) != self.ncols()):
            raise ValueError('Shifts length should be the column dimension.')

    def degree(self):
        r"""
        Return the degree of this matrix.

        For a given polynomial matrix, its degree is the maximum of the degrees
        of all its entries. If the matrix is nonzero, this is a nonnegative
        integer; here, the degree of the zero matrix is -1.

        OUTPUT: an integer.

        EXAMPLES::

            sage: pR.<x> = GF(7)[]
            sage: M = Matrix( pR, [[3*x+1, 0, 1], [x^3+3, 0, 0]])
            sage: M.degree()
            3

        The zero matrix has degree ``-1``::

            sage: M = Matrix( pR, 2, 3 )
            sage: M.degree()
            -1

        Any "empty" matrix has degree ``None``::

            sage: M = Matrix( pR, 3, 0 )
            sage: M.degree()

        """
        if self.nrows()==0 or self.ncols()==0:
            raise ValueError('Empty matrix does not have a degree.')
        return max([ self[i,j].degree()
            for i in range(self.nrows()) for j in range(self.ncols()) ])

    def degree_matrix(self, shifts=None, row_wise=True):
        r"""
        Return the matrix of the (shifted) degrees in this matrix.

        For a given polynomial matrix $M = (M_{i,j})_{i,j}$, its degree matrix
        is the matrix $(\deg(M_{i,j}))_{i,j}$ formed by the degrees of its
        entries. Here, the degree of the zero polynomial is $-1$.

        For given shifts $s_1,\ldots,s_m \in \ZZ$, the shifted degree
        matrix of $M$ is either $(\deg(M_{i,j})+s_j)_{i,j}$ if working
        row-wise, or $(\deg(M_{i,j})+s_i)_{i,j}$ if working column-wise. In the
        former case, $m$ has to be the number of columns of $M$; in the latter
        case, the number of its rows. Here, if $M_{i,j}=0$ then the
        corresponding entry in the shifted degree matrix is
        $\min(s_1,\ldots,s_m)-1$. For more on shifts and working row-wise
        versus column-wise, see the class documentation.

        INPUT:

        - ``shifts`` -- (optional, default: ``None``) list of integers;
          ``None`` is interpreted as ``shifts=[0,...,0]``.

        - ``row_wise`` -- (optional, default: ``True``) boolean, if ``True``
          then shifts apply to the columns of the matrix and otherwise to its
          rows (see the class description for more details).

        OUTPUT: an integer matrix.

        EXAMPLES::

            sage: pR.<x> = GF(7)[]
            sage: M = Matrix( pR, [[3*x+1, 0, 1], [x^3+3, 0, 0]])
            sage: M.degree_matrix()
            [ 1 -1  0]
            [ 3 -1 -1]

            sage: M.degree_matrix(shifts=[0,1,2])
            [ 1 -1  2]
            [ 3 -1 -1]

        The zero entries in the polynomial matrix can be identified in the
        (shifted) degree matrix as the entries equal to ``min(shifts)-1``::

            sage: M.degree_matrix(shifts=[-2,1,2])
            [-1 -3  2]
            [ 1 -3 -3]

        Using ``row_wise=False``, the function supports shifts applied to the
        rows of the matrix (which, in terms of modules, means that we are
        working column-wise, see the class documentation)::

            sage: M.degree_matrix(shifts=[-1,2], row_wise=False)
            [ 0 -2 -1]
            [ 5 -2 -2]
        """
        self._check_shift_dimension(shifts,row_wise)
        if shifts is None:
            return self.apply_map(lambda x: x.degree())
        from sage.matrix.constructor import Matrix
        zero_degree = min(shifts)-1
        if row_wise: 
            return Matrix( ZZ, [[ self[i,j].degree()+shifts[j]
                if self[i,j]!=0 else zero_degree
                for j in range(self.ncols()) ] for i in range(self.nrows())] )
        else:
            return Matrix( ZZ, [[ self[i,j].degree()+shifts[i]
                if self[i,j]!=0 else zero_degree
                for j in range(self.ncols()) ] for i in range(self.nrows())] )

    def row_degree(self, shifts=None):
        r"""
        Return the (shifted) row degree of this matrix.

        For a given polynomial matrix $M = (M_{i,j})_{i,j}$ with $m$ rows and
        $n$ columns, its row degree is the tuple $(d_1,\ldots,d_m)$ where $d_i
        = \max_j(\deg(M_{i,j}))$ for $1\leq i \leq m$. Thus, $d_i=-1$ if
        the $i$-th row of $M$ is zero, and $d_i \geq 0$ otherwise.

        For given shifts $s_1,\ldots,s_n \in \ZZ$, the shifted row degree of
        $M$ is $(d_1,\ldots,d_m)$ where $d_i = \max_j(\deg(M_{i,j})+s_j)$.
        Here, if the $i$-th row of $M$ is zero then $d_i
        =\min(s_1,\ldots,s_n)-1$; otherwise, $d_i$ is larger than this value.

        INPUT:

        - ``shifts`` -- (optional, default: ``None``) list of integers;
          ``None`` is interpreted as ``shifts=[0,...,0]``.

        OUTPUT: a list of integers.

        REFERENCES:

        - [Wol1974]_ (Section 2.5, without shifts), and [VBB1992]_ (Section 3).

        - Up to changes of signs, shifted row degrees coincide with the notion
          of *defect* commonly used in the rational approximation literature
          (see for example [Bec1992]_ ).

        EXAMPLES::

            sage: pR.<x> = GF(7)[]
            sage: M = Matrix(pR, [ [3*x+1, 0, 1], [x^3+3, 0, 0] ])
            sage: M.row_degree()
            [1, 3]

            sage: M.row_degree(shifts=[0,1,2])
            [2, 3]

        A zero row in a polynomial matrix can be identified in the (shifted)
        row degree as the entries equal to ``min(shifts)-1``::

            sage: M = Matrix(pR, [[3*x+1, 0, 1], [x^3+3, 0, 0], [0, 0, 0]])
            sage: M.row_degree()
            [1, 3, -1]

            sage: M.row_degree(shifts=[-2,1,2])
            [2, 1, -3]

        The row degree of an empty matrix ($0\times n$ or $m\times 0$) is
        not defined::
            
            sage: M = Matrix( pR, 0, 3 )
            sage: M.row_degree()
            []

            sage: M = Matrix( pR, 3, 0 )
            sage: M.row_degree()
            [None, None, None]
        """
        self._check_shift_dimension(shifts,row_wise=True)
        if self.ncols()==0 or self.nrows()==0:
            raise ValueError('Empty matrix does not have a degree.')
        if shifts is None:
            return [ max([ self[i,j].degree() for j in range(self.ncols()) ])
                    for i in range(self.nrows()) ]
        zero_degree = min(shifts)-1
        return [ max([ self[i,j].degree() + shifts[j]
            if self[i,j]!=0 else zero_degree
            for j in range(self.ncols()) ]) for i in range(self.nrows()) ]

    def column_degree(self, shifts=None):
        r"""
        Return the (shifted) column degree of this matrix.

        For a given polynomial matrix $M = (M_{i,j})_{i,j}$ with $m$ rows and
        $n$ columns, its column degree is the tuple $(d_1,\ldots,d_n)$ where
        $d_j = \max_i(\deg(M_{i,j}))$ for $1\leq j \leq n$. Thus, $d_j=-1$ if
        the $j$-th column of $M$ is zero, and $d_j \geq 0$ otherwise.

        For given shifts $s_1,\ldots,s_m \in \ZZ$, the shifted column degree of
        $M$ is $(d_1,\ldots,d_n)$ where $d_j = \max_i(\deg(M_{i,j})+s_i)$.
        Here, if the $j$-th column of $M$ is zero then $d_j =
        \min(s_1,\ldots,s_m)-1$; otherwise $d_j$ is larger than this value.

        INPUT:

        - ``shifts`` -- (optional, default: ``None``) list of integers;
          ``None`` is interpreted as ``shifts=[0,...,0]``.

        OUTPUT: a list of integers.

        REFERENCES:
        
        see the documentation of :meth:`row_degree`.

        EXAMPLES::

            sage: pR.<x> = GF(7)[]
            sage: M = Matrix(pR, [ [3*x+1, 0, 1], [x^3+3, 0, 0] ])
            sage: M.column_degree()
            [3, -1, 0]

            sage: M.column_degree(shifts=[0,2])
            [5, -1, 0]

        A zero column in a polynomial matrix can be identified in the (shifted)
        column degree as the entries equal to ``min(shifts)-1``::

            sage: M.column_degree(shifts=[-2,1])
            [4, -3, -2]

        The column degree of a $0\times n$ matrix is the list ``[None]*n``,
        while the column degree of a $m\times 0$ matrix is the empty list::

            sage: M = Matrix( pR, 0, 3 )
            sage: M.column_degree()
            [None, None, None]

            sage: M = Matrix( pR, 3, 0 )
            sage: M.column_degree()
            []
        """
        if self.ncols()==0 or self.nrows()==0:
            raise ValueError('Empty matrix does not have a degree.')
        self._check_shift_dimension(shifts,row_wise=False)
        if self.nrows()==0:
            return [None]*(self.ncols())
        if shifts is None:
            return [ max([ self[i,j].degree() for i in range(self.nrows()) ])
                    for j in range(self.ncols()) ]
        zero_degree = min(shifts)-1
        return [ max([ self[i,j].degree() + shifts[i]
            if self[i,j]!=0 else zero_degree
            for i in range(self.nrows()) ]) for j in range(self.ncols()) ]

    def leading_matrix(self, shifts=None, row_wise=True):
        r"""
        Return the (shifted) leading matrix of this matrix.

        Let $M$ be a univariate polynomial matrix in $\Bold{K}[x]^{m \times
        n}$. Working row-wise and without shifts, its leading matrix is the
        matrix in $\Bold{K}^{m \times n}$ formed by the leading coefficients of
        the entries of $M$ which reach the degree of the corresponding row.
  
        More precisely, if working row-wise, let $s_1,\ldots,s_n \in \ZZ$
        be a shift, and let $(d_1,\ldots,d_m)$ denote the shifted row degree of
        $M$. Then, the shifted leading matrix of $M$ is the matrix in
        $\Bold{K}^{m \times n}$ whose entry $i,j$ is the coefficient of degree
        $d_i-s_j$ of the entry $i,j$ of $M$. Going over the Laurent
        polynomials, the shifted leading matrix of $M$ can also be described as
        the coefficient of degree $0$ of the polynomial matrix
        $\mathrm{diag}(x^{-d_1},\ldots,x^{-d_m}) M
        \mathrm{diag}(x^{s_1},\ldots,x^{s_m})$ (which only has entries of
        nonpositive degree).

        If working column-wise, let $s_1,\ldots,s_m \in \ZZ$ be a shift,
        and let $(d_1,\ldots,d_n)$ denote the shifted column degree of $M$.
        Then, the shifted leading matrix of $M$ is the matrix in $\Bold{K}^{m
        \times n}$ whose entry $i,j$ is the coefficient of degree $d_j-s_i$ of
        the entry $i,j$ of $M$. Going over the Laurent polynomials, the shifted
        leading matrix of $M$ can also be described as the coefficient of
        degree $0$ of the polynomial matrix
        $\mathrm{diag}(x^{s_1},\ldots,x^{s_m}) M
        \mathrm{diag}(x^{-d_1},\ldots,x^{-d_m})$ (which only has entries of
        nonpositive degree).

        INPUT:

        - ``shifts`` -- (optional, default: ``None``) list of integers;
          ``None`` is interpreted as ``shifts=[0,...,0]``.

        - ``row_wise`` -- (optional, default: ``True``) boolean, ``True`` if
          working row-wise (see the class description).

        OUTPUT: a matrix over the base field.

        REFERENCES:
        
        [Wol1974]_ (Section 2.5, without shifts) and [VBB1992]_ (Section 3).

        EXAMPLES::

            sage: pR.<x> = GF(7)[]
            sage: M = Matrix(pR, [ [3*x+1, 0, 1], [x^3+3, 0, 0] ])
            sage: M.leading_matrix()
            [3 0 0]
            [1 0 0]

            sage: M.leading_matrix(shifts=[0,1,2])
            [0 0 1]
            [1 0 0]

            sage: M.leading_matrix(row_wise=False)
            [0 0 1]
            [1 0 0]

            sage: M.leading_matrix(shifts=[-2,1], row_wise=False)
            [0 0 1]
            [1 0 0]

            sage: M.leading_matrix(shifts=[2,0], row_wise=False)
            [3 0 1]
            [1 0 0]
        """
        self._check_shift_dimension(shifts,row_wise)
        from sage.matrix.constructor import Matrix
        if row_wise:
            row_degree = self.row_degree(shifts)
            if shifts is None:
                return Matrix([ [ self[i,j].leading_coefficient()
                    if self[i,j].degree() == row_degree[i] else 0
                    for j in range(self.ncols()) ]
                    for i in range(self.nrows()) ])
            return Matrix([ [ self[i,j].leading_coefficient()
                if self[i,j].degree() + shifts[j] == row_degree[i] else 0
                for j in range(self.ncols()) ]
                for i in range(self.nrows()) ])
        else:
            column_degree = self.column_degree(shifts)
            if shifts is None:
                return Matrix([ [ self[i,j].leading_coefficient()
                    if self[i,j].degree() == column_degree[j] else 0
                    for j in range(self.ncols()) ]
                    for i in range(self.nrows()) ])
            return Matrix([ [ self[i,j].leading_coefficient()
                if self[i,j].degree() + shifts[i] == column_degree[j] else 0
                for j in range(self.ncols()) ]
                for i in range(self.nrows()) ])

    def is_reduced(self, shifts=None, row_wise=True):
        #TODO should use def with zero rows
        r"""
        Return ``True`` if and only if this matrix is in (shifted) reduced
        form.

        An $m \times n$ univariate polynomial matrix $M$ is said to be in
        shifted row (resp.  column) reduced form if its shifted leading matrix
        has rank $m$, with $m \leq n$ (resp. $n$, with $n \leq m$).
        
        Equivalently, when considering all the matrices obtained by
        left-multiplying (resp. right-multiplying) $M$ by a unimodular matrix,
        then the shifted row (resp. column) degree of $M$ -- once sorted in
        nondecreasing order -- is lexicographically minimal.

        INPUT:

        - ``shifts`` -- (optional, default: ``None``) list of integers;
          ``None`` is interpreted as ``shifts=[0,...,0]``.

        - ``row_wise`` -- (optional, default: ``True``) boolean, ``True`` if
          working row-wise (see the class description)

        OUTPUT: a boolean value.

        REFERENCES:
        
        [Wol1974]_ (Section 2.5, without shifts) and [VBB1992]_ (Section 3).

        EXAMPLES::

            sage: pR.<x> = GF(7)[]
            sage: M = Matrix(pR, [ [3*x+1, 0, 1], [x^3+3, 0, 0] ])
            sage: M.is_reduced()
            False

            sage: M.is_reduced(shifts=[0,1,2])
            True

            sage: M.is_reduced(shifts=[2,0], row_wise=False)
            False

            sage: M = Matrix(pR, [ [3*x+1, 0, 1], [x^3+3, 0, 0], [0, 1, 0] ])
            sage: M.is_reduced(shifts=[2,0,0], row_wise=False)
            True
        """
        self._check_shift_dimension(shifts,row_wise)
        number_generators = self.nrows() if row_wise else self.ncols()
        return number_generators == self.leading_matrix(shifts, row_wise).rank()

    def leading_positions(self, shifts=None, row_wise=True, return_degree=False):
        r"""
        Return the (shifted) leading positions (also known as the pivot index),
        and optionally the (shifted) pivot degree of this matrix.

        If working row-wise, for a given shift $s_1,\ldots,s_n \in
        \ZZ$, taken as $(0,\ldots,0)$ by default, and a row vector of
        univariate polynomials $[p_1,\ldots,p_n]$, the leading positions of
        this vector is the index $j$ of the rightmost nonzero entry $p_j$ such
        that $\deg(p_j) + s_j$ is equal to the row degree of the vector. Then,
        for this index $j$, the pivot degree of the vector is the degree
        $\deg(p_j)$.
        
        For the zero row, both the leading positions and degree are $-1$.  For
        a $m \times n$ polynomial matrix, the leading positions and pivot
        degree are the two lists containing the leading positions and the pivot
        degree of its rows.

        The definition is similar if working column-wise.

        INPUT:

        - ``shifts`` -- (optional, default: ``None``) list of integers;
          ``None`` is interpreted as ``shifts=[0,...,0]``.

        - ``row_wise`` -- (optional, default: ``True``) boolean, ``True`` if
          working row-wise (see the class description).

        - ``return_degree`` -- (optional, default: ``False``) boolean, ``True``
          implies that the pivot degrees are returned.

        OUTPUT: a list of integers if ``return_degree=False``; a pair of lists
        of integers otherwise.

        REFERENCES:
        
        [Kai1980]_ (Section 6.7.2, without shifts).

        EXAMPLES::

            sage: pR.<x> = GF(7)[]
            sage: M = Matrix(pR, [ [3*x+1, 0, 1], [x^3+3, 0, 0] ])
            sage: M.leading_positions(return_degree=True)
            ([0, 0], [1, 3])

            sage: M.leading_positions(shifts=[0,5,2], return_degree=True)
            ([2, 0], [0, 3])

            sage: M.leading_positions(row_wise=False, return_degree=True)
            ([1, -1, 0], [3, -1, 0])

            sage: M.leading_positions(shifts=[1,2], row_wise=False, \
                    return_degree=True)
            ([1, -1, 0], [3, -1, 0])

        In case several entries in the row (resp. column) reach the shifted row
        (resp. column) degree, the leading position is chosen as the rightmost
        (resp. bottommost) such entry::

            sage: M.leading_positions(shifts=[0,5,1],return_degree=True)
            ([2, 0], [0, 3])

            sage: M.leading_positions(shifts=[2,0], row_wise=False,return_degree=True)
            ([1, -1, 0], [3, -1, 0])

        If working row-wise, both the leading positions index and pivot degree
        of a $0\times n$ matrix are the empty list ``[]``, while for a $m\times
        0$ matrix both are the list ``[None]*m``. A similar property holds when
        working column-wise::
            
            sage: M = Matrix( pR, 0, 3 )
            sage: M.leading_positions(return_degree=True)
            ([], [])

            sage: M.leading_positions(row_wise=False,return_degree=True)
            ([None, None, None], [None, None, None])

            sage: M = Matrix( pR, 3, 0 )
            sage: M.leading_positions(return_degree=True)
            ([None, None, None], [None, None, None])

            sage: M.leading_positions(row_wise=False,return_degree=True)
            ([], [])
        """
        self._check_shift_dimension(shifts,row_wise)

        if self.ncols()==0 or self.nrows()==0:
            raise ValueError('Empty matrix does not have leading positions.')

        if row_wise:
            row_degree = self.row_degree(shifts)
            if shifts is None:
                pivot_index = [ (-1 if row_degree[i]==-1 else
                    max( [ j for j in range(self.ncols()) if
                    (self[i,j].degree() == row_degree[i]) ] ))
                    for i in range(self.nrows()) ]
            else:
                zero_degree=min(shifts)-1
                pivot_index = [ (-1 if row_degree[i]==zero_degree else
                    max( [ j for j in range(self.ncols()) if
                    (self[i,j]!=0 and
                    self[i,j].degree()+shifts[j] == row_degree[i]) ] ))
                    for i in range(self.nrows()) ]
            pivot_degree = [ (-1 if pivot_index[i]==-1 else
                self[i,pivot_index[i]].degree())
                for i in range(self.nrows()) ]
            return (pivot_index,pivot_degree) if return_degree else pivot_index
                    
        # now in the column-wise case
        column_degree = self.column_degree(shifts)
        if shifts is None:
            pivot_index = [ (-1 if column_degree[j]==-1 else
                max( [ i for i in range(self.nrows()) if
                (self[i,j].degree() == column_degree[j]) ] ))
                for j in range(self.ncols()) ]
        else:
            zero_degree=min(shifts)-1
            pivot_index = [ (-1 if column_degree[j]==zero_degree else
                max( [ i for i in range(self.nrows()) if
                (self[i,j]!=0 and
                self[i,j].degree()+shifts[i] == column_degree[j]) ] ))
                for j in range(self.ncols()) ]
        pivot_degree = [ (-1 if pivot_index[j]==-1 else
            self[pivot_index[j],j].degree())
            for j in range(self.ncols()) ]
        return (pivot_index,pivot_degree) if return_degree else pivot_index

    def is_weak_popov(self, shifts=None, row_wise=True, ordered=False):
        #TODO should use def with zero rows
        r"""
        Return a boolean indicating whether this matrix is in (shifted)
        (ordered) weak Popov form.

        If working row-wise (resp. column-wise), a polynomial matrix is said to
        be in weak Popov form if it has no zero row (resp. column) and its
        pivot index has pairwise distinct entries (for the ordered weak Popov
        form, this pivot index must be strictly increasing).

        Concerning square matrices, the ordered weak Popov form is sometimes
        also called the quasi-Popov form.

        INPUT:

        - ``shifts`` -- (optional, default: ``None``) list of integers;
          ``None`` is interpreted as ``shifts=[0,...,0]``.

        - ``row_wise`` -- (optional, default: ``True``) boolean, ``True`` if
          working row-wise (see the class description).

        - ``ordered`` -- (optional, default: ``False``) boolean, ``True`` if
          checking for an ordered weak Popov form.

        OUTPUT: a boolean.

        REFERENCES:
        
        [Kai1980]_ (Section 6.7.2, square case without shifts), [MS2003]_
        (without shifts), [BLV1999]_ .

        SEEALSO:  :meth:`weak_popov_form` .

        EXAMPLES::

            sage: pR.<x> = GF(7)[]
            sage: M = Matrix([ [x^3+3*x^2+6*x+6, 3*x^2+3*x+6, 4*x^2+x+3], \
                               [5,               1,           0        ], \
                               [2*x^2+2,         2*x+5,       x^2+4*x+6] ])
            sage: M.is_weak_popov()
            True

        One can check whether the pivot index, in addition to being pairwise
        distinct, are actually in increasing order::

            sage: M.is_weak_popov(ordered=True)
            True

            sage: N = M.with_swapped_rows(1,2)
            sage: N.is_weak_popov(ordered=True)
            False

        Shifts and orientation (row-wise or column-wise) are supported::

            sage: M.is_weak_popov(shifts=[2,3,1])
            False

            sage: M.is_weak_popov(shifts=[0,2,0],row_wise=False,ordered=True)
            True

        Rectangular matrices are supported::

            sage: M = Matrix([ \
                    [  x^3+5*x^2+5*x+1,       5,       6*x+4,         0], \
                    [      6*x^2+3*x+1,       1,           2,         0], \
                    [2*x^3+4*x^2+6*x+4, 5*x + 1, 2*x^2+5*x+5, x^2+5*x+6]  \
                    ])
            sage: M.is_weak_popov(shifts=[0,2,1,3])
            True

            sage: M.is_weak_popov(shifts=[0,2,1,3],ordered=True)
            True
        """
        self._check_shift_dimension(shifts,row_wise)
        pivot_index = self.leading_positions(shifts, row_wise)
        # pivot index should not have duplicates, which is equivalent to:
        # once sorted, it doesn't contain a pair of equal successive entries
        if not ordered:
            pivot_index.sort()
        # there should be no zero row (or column if not row_wise)
        # and the matrix should not be m x 0 (or 0 x n if not row_wise)
        if len(pivot_index) > 0 and \
            (pivot_index[0] == None or pivot_index[0] < 0):
            return False
        # pivot_index should be strictly increasing
        for index,next_pivot_index in enumerate(pivot_index[1:]):
            if next_pivot_index <= pivot_index[index]:
                return False
        return True

    def is_popov(self, shifts=None, row_wise=True, ordered_by_degree=False):
        #TODO should use def with zero rows
        #TODO order by degree version
        r"""
        Return a boolean indicating whether this matrix is in (shifted) Popov
        form.

        If working row-wise (resp. column-wise), a polynomial matrix is said to
        be in Popov form if it has no zero row (resp. column), its pivot index
        is strictly increasing, and for each row (resp. column) the pivot entry
        is monic and has degree strictly larger than the other entries in its
        column (resp. row).

        TODO. There is another convention, which replaces "pivot index strictly
        increasing" by "row (resp. column) degree nondecreasing, and for rows
        (resp. columns) of same degree, pivot indices increasing".

        INPUT:

        - ``shifts`` -- (optional, default: ``None``) list of integers;
          ``None`` is interpreted as ``shifts=[0,...,0]``.

        - ``row_wise`` -- (optional, default: ``True``) boolean, ``True`` if
          working row-wise (see the class description).

        OUTPUT: a boolean.

        REFERENCES:
        
        For the square case, without shifts: [Pop1972]_ and [Kai1980]_ (Section
        6.7.2). For the general case: [BLV2006]_ .

        EXAMPLES::

            sage: pR.<x> = GF(7)[]
            sage: M = Matrix(pR, [ [x^4+6*x^3+4*x+4, 3*x+6,     3  ], \
                                   [x^2+6*x+6,       x^2+5*x+5, 2  ], \
                                   [3*x,             6*x+5,     x+5] ])
            sage: M.is_popov()
            True

            sage: M.is_popov(shifts=[0,1,2])
            True

            sage: M[:,:2].is_popov()
            False

            sage: M[:2,:].is_popov(shifts=[0,1,2])
            True

            sage: M = Matrix(pR, [ [x^4+3*x^3+x^2+2*x+6, x^3+5*x^2+5*x+1], \
                                   [6*x+1,               x^2 + 4*x + 1  ], \
                                   [6,                   6              ] ])
            sage: M.is_popov(row_wise=False)
            False

            sage: M.is_popov(shifts=[0,2,3], row_wise=False)
            True
        """
        self._check_shift_dimension(shifts,row_wise)
        pivot_index,pivot_degree = self.leading_positions(shifts, row_wise,
                return_degree=True)
        # there should be no zero row (or column if not row_wise)
        # and the matrix should not be m x 0 (or 0 x n if not row_wise)
        if len(pivot_index) > 0 and \
            (pivot_index[0] == None or pivot_index[0] < 0):
            return False
        # pivot_index should be strictly increasing
        for index,next_pivot_index in enumerate(pivot_index[1:]):
            if next_pivot_index <= pivot_index[index]:
                return False
        # pivot entries should be monic, and pivot degrees must be the greatest
        # in their column (or in their row if column-wise)
        for i,index in enumerate(pivot_index):
            if row_wise:
                if not self[i,index].is_monic():
                    return False
                for k in range(self.nrows()):
                    if k==i:
                        continue
                    if self[k,index].degree() >= pivot_degree[i]:
                        return False
            else: # now column-wise
                if not self[index,i].is_monic():
                    return False
                for k in range(self.ncols()):
                    if k==i:
                        continue
                    if self[index,k].degree() >= pivot_degree[i]:
                        return False
        return True

    #def is_hermite(self, row_wise=True, lower_tri=True):
    # TODO

    #def is_weak_popov(self):
    #    r"""
    #    Return ``True`` if this matrix is in weak Popov form.

    #    OUTPUT:

    #    A matrix over a polynomial ring is in weak Popov form if all
    #    leading positions are different [MS2003]_. A leading position
    #    is the position `i` in a row with the highest degree; in case of tie,
    #    the maximal `i` is used (i.e. furthest to the right).

    #    EXAMPLES:

    #    A matrix with the same leading position in two rows is not in weak
    #    Popov form::

    #        sage: PF = PolynomialRing(GF(2^12,'a'),'x')
    #        sage: A = matrix(PF,3,[x,   x^2, x^3,\
    #                               x^2, x^2, x^2,\
    #                               x^3, x^2, x    ])
    #        sage: A.is_weak_popov()
    #        False

    #    If a matrix has different leading positions, it is in weak Popov
    #    form::

    #        sage: B = matrix(PF,3,[1,    1,  x^3,\
    #                               x^2,  1,  1,\
    #                               1,x^  2,  1  ])
    #        sage: B.is_weak_popov()
    #        True

    #    Weak Popov form is not restricted to square matrices::

    #        sage: PF = PolynomialRing(GF(7),'x')
    #        sage: D = matrix(PF,2,4,[x^2+1, 1, 2, x,\
    #                                 3*x+2, 0, 0, 0 ])
    #        sage: D.is_weak_popov()
    #        False

    #    Even a matrix with more rows than columns can still be in weak Popov
    #    form::

    #        sage: E = matrix(PF,4,2,[4*x^3+x, x^2+5*x+2,\
    #                                 0,       0,\
    #                                 4,       x,\
    #                                 0,       0         ])
    #        sage: E.is_weak_popov()
    #        True

    #    A matrix with fewer columns than non-zero rows is never in weak
    #    Popov form::

    #        sage: F = matrix(PF,3,2,[x^2,   x,\
    #                                 x^3+2, x,\
    #                                 4,     5])
    #        sage: F.is_weak_popov()
    #        False

    #    TESTS:

    #    Verify tie breaking by selecting right-most index::

    #        sage: F = matrix(PF,2,2,[x^2, x^2,\
    #                                 x,   5   ])
    #        sage: F.is_weak_popov()
    #        True

    #    .. SEEALSO::

    #        - :meth:`weak_popov_form <sage.matrix.matrix_polynomial_dense.weak_popov_form>`

    #    AUTHOR:

    #    - David Moedinger (2014-07-30)
    #    """
    #    t = set()
    #    for r in range(self.nrows()):
    #        max = -1
    #        for c in range(self.ncols()):
    #            if self[r, c].degree() >= max:
    #                max = self[r, c].degree()
    #                p = c
    #        if not max == -1:
    #            if p in t:
    #                return False
    #            t.add(p)
    #    return True

    #def order_by_degree(self, shifts=None, row_wise=True):
    # TODO

    #def order_by_leading_position(self, shifts=None, row_wise=True):
    # TODO

    def weak_popov_form(self, transformation=False, shifts=None):
        r"""
        Return a weak Popov form of this matrix.

        A matrix is in weak Popov form if the leading positions of the nonzero
        rows are all different. The leading position of a row is the right-most
        position whose entry has the maximal degree in the row.

        The weak Popov form is non-canonical, so an input matrix have many weak
        Popov forms.

        INPUT:

        - ``transformation`` -- boolean (default: ``False``) If ``True``, the
          transformation matrix is returned together with the weak Popov form.

        - ``shifts`` -- (default: ``None``) A tuple or list of integers
          `s_1, \ldots, s_n`, where `n` is the number of columns of the matrix.
          If given, a "shifted weak Popov form" is computed, i.e. such that the
          matrix `A\,\mathrm{diag}(x^{s_1}, \ldots, x^{s_n})` is in weak Popov
          form, where `\mathrm{diag}` denotes a diagonal matrix.

        ALGORITHM:

        This method implements the Mulders-Storjohann algorithm of [MS2003]_.

        EXAMPLES::

            sage: F.<a> = GF(2^4,'a')
            sage: PF.<x> = F[]
            sage: A = matrix(PF,[[1,  a*x^17 + 1 ],\
                                 [0,  a*x^11 + a^2*x^7 + 1 ]])
            sage: M, U = A.weak_popov_form(transformation=True)
            sage: U * A == M
            True
            sage: M.is_weak_popov()
            True
            sage: U.is_invertible()
            True

        A zero matrix will return itself::

            sage: Z = matrix(PF,5,3)
            sage: Z.weak_popov_form()
            [0 0 0]
            [0 0 0]
            [0 0 0]
            [0 0 0]
            [0 0 0]

        Shifted weak popov form is computed if ``shifts`` is given::

            sage: PF.<x> = QQ[]
            sage: A = matrix(PF,3,[x,   x^2, x^3,\
                                   x^2, x^1, 0,\
                                   x^3, x^3, x^3])
            sage: A.weak_popov_form()
            [        x       x^2       x^3]
            [      x^2         x         0]
            [  x^3 - x x^3 - x^2         0]
            sage: H,U = A.weak_popov_form(transformation=True, shifts=[16,8,0])
            sage: H
            [               x              x^2              x^3]
            [               0         -x^2 + x       -x^4 + x^3]
            [               0                0 -x^5 + x^4 + x^3]
            sage: U * A == H
            True

        .. SEEALSO::

            :meth:`is_weak_popov <sage.matrix.matrix_polynomial_dense.is_weak_popov>`
        """
        self._check_shift_dimension(shifts,row_wise=True)
        M = self.__copy__()
        U = M._weak_popov_form(transformation=transformation, shifts=shifts)
        M.set_immutable()
        if transformation:
            U.set_immutable()
        return (M,U) if transformation else M

    def _weak_popov_form(self, transformation=False, shifts=None):
        """
        Transform this matrix in place into weak Popov form.

        EXAMPLES::

            sage: F.<a> = GF(2^4,'a')
            sage: PF.<x> = F[]
            sage: A = matrix(PF,[[1,  a*x^17 + 1 ],\
                                 [0,  a*x^11 + a^2*x^7 + 1 ]])
            sage: M = A.__copy__()
            sage: U = M._weak_popov_form(transformation=True)
            sage: U * A == M
            True
            sage: M.is_weak_popov()
            True
            sage: U.is_invertible()
            True

            sage: PF.<x> = QQ[]
            sage: A = matrix(PF,3,[x,   x^2, x^3,\
                                   x^2, x^1, 0,\
                                   x^3, x^3, x^3])
            sage: A.weak_popov_form()
            [        x       x^2       x^3]
            [      x^2         x         0]
            [  x^3 - x x^3 - x^2         0]
            sage: M = A.__copy__()
            sage: U = M._weak_popov_form(transformation=True, shifts=[16,8,0])
            sage: M
            [               x              x^2              x^3]
            [               0         -x^2 + x       -x^4 + x^3]
            [               0                0 -x^5 + x^4 + x^3]
            sage: U * A == M
            True
        """
        cdef Py_ssize_t i, j
        cdef Py_ssize_t c, d, best, bestp

        cdef Py_ssize_t m = self.nrows()
        cdef Py_ssize_t n = self.ncols()

        cdef Matrix M = self
        cdef Matrix U

        cdef list to_row, conflicts

        R = self.base_ring()
        one = R.one()

        if transformation:
            from sage.matrix.constructor import identity_matrix
            U = identity_matrix(R, m)

        if shifts and len(shifts) != M.ncols():
            raise ValueError("the number of shifts must equal the number of columns")

        # initialise to_row and conflicts list
        to_row = [[] for i in range(n)]
        conflicts = []
        for i in range(m):
            bestp = -1
            best = -1
            for c in range(n):
                d = M.get_unsafe(i,c).degree()

                if shifts and d >= 0 :
                    d += shifts[c]

                if d >= best:
                    bestp = c
                    best = d

            if best >= 0:
                to_row[bestp].append((i,best))
                if len(to_row[bestp]) > 1:
                    conflicts.append(bestp)

        # while there is a conflict, do a simple transformation
        while conflicts:
            c = conflicts.pop()
            row = to_row[c]
            i,ideg = row.pop()
            j,jdeg = row.pop()

            if jdeg > ideg:
                i,j = j,i
                ideg,jdeg = jdeg,ideg

            coeff = - M.get_unsafe(i,c).lc() / M.get_unsafe(j,c).lc()
            s = coeff * one.shift(ideg - jdeg)

            M.add_multiple_of_row_c(i, j, s, 0)
            if transformation:
                U.add_multiple_of_row_c(i, j, s, 0)

            row.append((j,jdeg))

            bestp = -1
            best = -1
            for c in range(n):
                d = M.get_unsafe(i,c).degree()

                if shifts and d >= 0:
                    d += shifts[c]

                if d >= best:
                    bestp = c
                    best = d

            if best >= 0:
                to_row[bestp].append((i,best))
                if len(to_row[bestp]) > 1:
                    conflicts.append(bestp)

        if transformation:
            return U

    def reduced_form(self, transformation=None, shifts=None, row_wise=True):
        r"""
        Return a row reduced form of this matrix.

        A matrix `M` is row reduced if the (row-wise) leading term matrix has
        the same rank as `M`. The (row-wise) leading term matrix of a polynomial
        matrix `M` is the matrix over `k` whose `(i,j)`'th entry is the
        `x^{d_i}` coefficient of `M[i,j]`, where `d_i` is the greatest degree
        among polynomials in the `i`'th row of `M_0`.

        A row reduced form is non-canonical so a given matrix has many row
        reduced forms.

        INPUT:

        - ``transformation`` -- (default: ``False``). If this ``True``, the
          transformation matrix `U` will be returned as well: this is an
          invertible matrix over `k[x]` such that ``self`` equals `UW`, where
          `W` is the output matrix.

        - ``shifts`` -- (default: ``None``) A tuple or list of integers
          `s_1, \ldots, s_n`, where `n` is the number of columns of the matrix.
          If given, a "shifted row reduced form" is computed, i.e. such that the
          matrix `A\,\mathrm{diag}(x^{s_1}, \ldots, x^{s_n})` is row reduced, where
          `\mathrm{diag}` denotes a diagonal matrix.

        OUTPUT:

        - `W` -- a row reduced form of this matrix.

        EXAMPLES::

            sage: R.<t> = GF(3)['t']
            sage: M = matrix([[(t-1)^2],[(t-1)]])
            sage: M.reduced_form()
            [    0]
            [t + 2]

        If the matrix is an `n \times 1` matrix with at least one non-zero entry,
        `W` has a single non-zero entry and that entry is a scalar multiple of
        the greatest-common-divisor of the entries of the matrix::

            sage: M1 = matrix([[t*(t-1)*(t+1)],[t*(t-2)*(t+2)],[t]])
            sage: output1 = M1.reduced_form()
            sage: output1
            [0]
            [0]
            [t]

        The following is the first half of example 5 in [Hes2002]_ *except* that we
        have transposed the matrix; [Hes2002]_ uses column operations and we use row::

            sage: R.<t> = QQ['t']
            sage: M = matrix([[t^3 - t,t^2 - 2],[0,t]]).transpose()
            sage: M.reduced_form()
            [      t    -t^2]
            [t^2 - 2       t]

        The next example demonstrates what happens when the matrix is a zero matrix::

            sage: R.<t> = GF(5)['t']
            sage: M = matrix(R, 2, [0,0,0,0])
            sage: M.reduced_form()
            [0 0]
            [0 0]

        In the following example, the original matrix is already row reduced, but
        the output is a different matrix. This is because currently this method
        simply computes a weak Popov form, which is always also a row reduced matrix
        (see :meth:`weak_popov_form`). This behavior is likely to change when a faster
        algorithm designed specifically for row reduced form is implemented in Sage::

            sage: R.<t> = QQ['t']
            sage: M = matrix([[t,t,t],[0,0,t]]); M
            [t t t]
            [0 0 t]
            sage: M.reduced_form()
            [ t  t  t]
            [-t -t  0]

        The last example shows the usage of the transformation parameter::

            sage: Fq.<a> = GF(2^3)
            sage: Fx.<x> = Fq[]
            sage: A = matrix(Fx,[[x^2+a,x^4+a],[x^3,a*x^4]])
            sage: W,U = A.reduced_form(transformation=True);
            sage: W,U
            (
            [          x^2 + a           x^4 + a]  [1 0]
            [x^3 + a*x^2 + a^2               a^2], [a 1]
            )
            sage: U*W == A
            True
            sage: U.is_invertible()
            True

        """
        self._check_shift_dimension(shifts,row_wise)
        if not row_wise:
            return self.T.reduced_form(transformation, shifts, row_wise=True).T
        return self.weak_popov_form(transformation, shifts)

    row_reduced_form = deprecated_function_alias(23619, reduced_form)

    def hermite_form(self, include_zero_rows=True, transformation=False):
        """
        Return the Hermite form of this matrix.

        The Hermite form is also normalized, i.e., the pivot polynomials
        are monic.

        INPUT:

        - ``include_zero_rows`` -- boolean (default: ``True``); if ``False``,
          the zero rows in the output matrix are deleted

        - ``transformation`` -- boolean (default: ``False``); if ``True``,
          return the transformation matrix

        OUTPUT:

        - the Hermite normal form `H` of this matrix `A`

        - (optional) transformation matrix `U` such that `UA = H`

        EXAMPLES::

            sage: M.<x> = GF(7)[]
            sage: A = matrix(M, 2, 3, [x, 1, 2*x, x, 1+x, 2])
            sage: A.hermite_form()
            [      x       1     2*x]
            [      0       x 5*x + 2]
            sage: A.hermite_form(transformation=True)
            (
            [      x       1     2*x]  [1 0]
            [      0       x 5*x + 2], [6 1]
            )
            sage: A = matrix(M, 2, 3, [x, 1, 2*x, 2*x, 2, 4*x])
            sage: A.hermite_form(transformation=True, include_zero_rows=False)
            ([  x   1 2*x], [0 4])
            sage: H, U = A.hermite_form(transformation=True, include_zero_rows=True); H, U
            (
            [  x   1 2*x]  [0 4]
            [  0   0   0], [5 1]
            )
            sage: U * A == H
            True
            sage: H, U = A.hermite_form(transformation=True, include_zero_rows=False)
            sage: U * A
            [  x   1 2*x]
            sage: U * A == H
            True
        """
        A = self.__copy__()
        U = A._hermite_form_euclidean(transformation=transformation,
                                      normalization=lambda p: ~p.lc())
        if not include_zero_rows:
            i = A.nrows() - 1
            while i >= 0 and A.row(i) == 0:
                i -= 1
            A = A[:i+1]
            if transformation:
                U = U[:i+1]

        A.set_immutable()
        if transformation:
            U.set_immutable()

        return (A, U) if transformation else A
