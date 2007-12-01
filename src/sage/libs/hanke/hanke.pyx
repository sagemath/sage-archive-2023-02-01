"""nodoctest
Interface to the Hanke Quadratic Forms Library
"""

doc=r"""
Interface to the Hanke Quadratic Forms Library

AUTHOR: Jonathan Hanke (Python integration: William Stein)

This package provides an interface to some of the C++ code that Jon
Hanke wrote for certain quadratic forms computations related to the
290 problem.
"""



include 'interrupt.pxi'

cdef extern from "Python.h":
    int PyErr_CheckSignals()

cdef extern from "stdlib.h":
    void free(void *ptr)

cdef extern from "wrap.h":
    struct Matrix_mpz
    Matrix_mpz* Matrix_mpz_new(int r, int s)
    void Matrix_mpz_del(Matrix_mpz* m)
    char* Matrix_mpz_repr(Matrix_mpz* m)
    int Matrix_mpz_nrows(Matrix_mpz* m)
    int Matrix_mpz_ncols(Matrix_mpz* m)
    void Matrix_mpz_setitem(Matrix_mpz* x, int i, int j, char* z)
    char* Matrix_mpz_getitem(Matrix_mpz* x, int i, int j)
    char* Matrix_mpz_determinant(Matrix_mpz* x)
    Matrix_mpz* Matrix_mpz_adjoint(Matrix_mpz* m)
    char* Matrix_mpz_Local_Density(Matrix_mpz* m, char* p, char* m)
    char* Matrix_mpz_Local_Primitive_Density(Matrix_mpz* m, char* p, char* m)
    char* Matrix_mpz_level(Matrix_mpz* x)

    void Matrix_mpz_symmetric_swap(Matrix_mpz* x, int i, int j)
    void Matrix_mpz_symmetric_multiply(Matrix_mpz* x, int i, char* y)
    void Matrix_mpz_symmetric_divide(Matrix_mpz* x, int i, char* y)
    void Matrix_mpz_symmetric_add(Matrix_mpz* x, int i, int j, char* y)

    Matrix_mpz* Matrix_mpz_local_normal_form(Matrix_mpz* x, char* p)

    Matrix_mpz* Matrix_mpz_local_diagonal_form(Matrix_mpz* x, char* p)

    long Matrix_mpz_hasse_invariant(Matrix_mpz* x, char* p)

    int Matrix_mpz_is_anisotropic(Matrix_mpz* x, char* p)

    int Matrix_mpz_is_isotropic(Matrix_mpz* x, char* p)

    int Matrix_mpz_is_quadratic_form(Matrix_mpz* x)

    int Matrix_mpz_is_symmetric(Matrix_mpz* x)

    char* Matrix_mpz_anisotropic_primes(Matrix_mpz* x)

    char* Matrix_mpz_local_constants(Matrix_mpz* x, char* p, char* T)

    int Matrix_mpz_is_stable(Matrix_mpz* x, char* p, char* T)

    int Matrix_mpz_cmp(Matrix_mpz* x, Matrix_mpz* y)

cdef object string(char* s):
    _sig_off
    # Makes a python string and deletes what is pointed to by s.
    t = str(s)
    free(s)
    return t

class __init:
    pass
_INIT = __init()


cdef class hanke_Matrix:
    """
    A matrix over the integers, with functions supporting computing
    with the quadratic form it defines, when it is symmetric.
    """
    cdef Matrix_mpz* x
    cdef int _nrows, _ncols

    def __init__(self, r, s=0, v=[]):
        """
        The following examples illustrates three ways to create
        a hanke_Matrix.

        EXAMPLES:
            sage: hanke_Matrix(2,3, [1,2,3, 4,5,6])
            [ 1, 2, 3 ]
            [ 4, 5, 6 ]
            sage: hanke_Matrix(2, [1,2,3,4,])
            [ 1, 2 ]
            [ 3, 4 ]
            sage: hanke_Matrix(2,2)
            [ 0, 0 ]
            [ 0, 0 ]
        """
        if r is _INIT:
            return
        if s == 0:
            s = r
        if isinstance(s, list):
            v = s
            s = r
        _sig_on
        self.x = Matrix_mpz_new(r,s)
        _sig_off
        self._nrows = r
        self._ncols = s
        if len(v) > 0:
            self.set(v)

    def __dealloc__(self):
        Matrix_mpz_del(self.x)

    def __repr__(self):
        """
        EXAMPLE:
            sage: A = hanke_Matrix(3, range(9))
            sage: A
            [ 0, 1, 2 ]
            [ 3, 4, 5 ]
            [ 6, 7, 8 ]
        """
        _sig_on
        s = string(Matrix_mpz_repr(self.x))
        return s[:-1].replace(" [","[")

    def __setitem__(self, ij, x):
        """
        Set the i,j and j,i entries of self to value (so matrix stays symmetric).

        EXAMPLES:
            sage: Q = hanke_Matrix(2); Q
            [ 0, 0 ]
            [ 0, 0 ]
            sage: Q[0,1]=5
            sage: Q
            [ 0, 5 ]
            [ 0, 0 ]
        """
        cdef int i, j
        if not isinstance(ij, tuple) or len(ij) != 2:
            raise TypeError, 'ij must be a 2-tuple'
        i, j = ij
        if i < 0 or i >= self._nrows or j < 0 or j >= self._ncols:
            raise IndexError, "array index out of range"
        t = str(x)
        Matrix_mpz_setitem(self.x, i, j, t)

    def __getitem__(self, ij):
        """
        Access a matrix by Q[i, j] to get entry i, j of matrix corresponding to Q.

        EXAMPLES:
            sage: A = hanke_Matrix(2)
            sage: A[1,1] = 10
            sage: A[1,1]
            10
        """
        cdef int i, j
        if not isinstance(ij, tuple) or len(ij) != 2:
            raise TypeError, 'ij must be a 2-tuple'
        i, j = ij
        if i < 0 or i >= self._nrows or j < 0 or j >= self._ncols:
            raise IndexError, "array index out of range"
        return int(Matrix_mpz_getitem(self.x, i, j))

    ######################################################
    # Functions for setting the entries
    ######################################################
    def set_diagonal(self, entries):
        """
        Set the diagonal elements of the matrix to the entries in the list.

        EXAMPLES:
            sage: A = hanke_Matrix(2)
            sage: A.set_diagonal([7, 9])
            sage: A
            [ 7, 0 ]
            [ 0, 9 ]
            sage: A.set_diagonal([7, 9, 5])
            Traceback (most recent call last):
                ...
            ArithmeticError: entries must have length 2
        """
        cdef int i, m, n
        m = self._nrows
        n = self._ncols
        if n != m:
            raise IndexError, "matrix must be square"
        if not isinstance(entries, list):
            entries = list(entries)
        if len(entries) != n:
            raise ArithmeticError, "entries must have length %s"%n
        for i from 0 <= i < n:
            if PyErr_CheckSignals(): raise KeyboardInterrupt
            t = str(entries[i])
            Matrix_mpz_setitem(self.x, i, i, t)


    def set(self, entries):
        """
        Set the entries of self to the given list of entries.

        INPUT:
            entries -- list of n*m integers, which fill in the
                       matrix row by row.

        EXAMPLES:
            sage: Q = hanke_Matrix(2,[2,5,5,2])
            sage: Q
            [ 2, 5 ]
            [ 5, 2 ]
            sage: Q.set([4,1,1,4]); Q
            [ 4, 1 ]
            [ 1, 4 ]
            sage: Q.set([1,2,3])
            Traceback (most recent call last):
            ...
            ArithmeticError: entries must have length 4
        """
        cdef int i, j, k, m, n
        m = self._nrows
        n = self._ncols
        if not isinstance(entries, list):
            entries = list(entries)
        if len(entries) != m*n:
            raise ArithmeticError, "entries must have length %s"%(m*n)
        k = 0
        for i from 0 <= i < m:
            if PyErr_CheckSignals(): raise KeyboardInterrupt
            for j from 0 <= j < n:
                t = str(entries[k])
                Matrix_mpz_setitem(self.x, i, j, t)
                k = k + 1

    def list(self):
        """
        Return the matrix as a single list of int's.  This is useful
        for object serialization (loading and saving to disk).

        EXAMPLES:
            sage: a = hanke_Matrix(2)
            sage: print a.list()
            [0, 0, 0, 0]
            sage: a = hanke_Matrix(2); a.set(range(4));
            sage: print a.list()
            [0, 1, 2, 3]
        """
        cdef int i, j
        v = []
        for i from 0 <= i < self._nrows:
            for j from 0 <= j < self._ncols:
                v.append(int(Matrix_mpz_getitem(self.x, i, j)))
        return v

    ######################################################
    # Very simple invariants
    ######################################################

    def nrows(self):
        """
        The number of rows of the matrix.

        EXAMPLES:
            sage: Q = hanke_Matrix(10);
            sage: Q.nrows()
            10
        """
        return self._nrows

    def ncols(self):
        """
        The number of columns of the matrix.

        EXAMPLES:
            sage: Q = hanke_Matrix(10);
            sage: Q.ncols()
            10
        """
        return self._ncols

    ######################################################
    # Required for object persistence (e.g., via Pickle)
    ######################################################
    def __reduce__(self):
        import hanke_py
        return hanke_py.init, (self.nrows(), self.ncols(), self.list())

    def __getstate__(self):
        """
        EXAMPLES:
            sage: Q = hanke_Matrix(2,[2,5,5,2]); Q
            [ 2, 5 ]
            [ 5, 2 ]
            sage: import cPickle
            sage: s = cPickle.dumps(Q)
            sage: R = cPickle.loads(s); R
            [ 2, 5 ]
            [ 5, 2 ]
            sage: R.determinant()
            -21
            sage: Q = hanke_Matrix(50,[2094084390852 for _ in xrange(50*50)])
            sage: s = cPickle.dumps(Q)
            sage: R = cPickle.loads(s)
            sage: Q == R
            True
        """
        return {'m':self._m, 'n':self._n, 'entries':self.list()}

    def __setstate__(self, v):
        self._nrows = v['m']
        self._ncols = v['n']
        self.x = Matrix_mpz_new(self._nrows, self._ncols)
        self.set(v['entries'])

    ######################################################
    # Determinant, adjoint, and level
    ######################################################
    def determinant(self):
        """
        The determinant of the matrix.

        EXAMPLES:
            sage: Q = hanke_Matrix(2,[2,5,5,2])
            sage: Q
            [ 2, 5 ]
            [ 5, 2 ]
            sage: Q.determinant()
            -21
        """
        if self._nrows != self._ncols:
            raise ArithmeticError, "matrix must be square.";
        if self._nrows == 0:
            return 1
        _sig_on
        d = Matrix_mpz_determinant(self.x)
        _sig_off
        return int(d)


    def adjoint(self):
        r"""
        The adjoint of this matrix.

        \begin{notice}

        This function is very slow for matrices that are about
        $20\times 20$ or larger.

        \end{notice}

        EXAMPLES:
            sage: Q = hanke_Matrix(2,[2,5,5,2])
            sage: Q.adjoint()
            [ 2, -5 ]
            [ -5, 2 ]
        """
        _sig_on
        return make_hanke_Matrix(Matrix_mpz_adjoint(self.x))

    def local_density(self, p, m):
        """
        Returns the local density at p.

        INPUT:
            p -- integer, a prime number
            m -- integer

        OUTPUT:
            A rational number (represented as a string)

        EXAMPLES:
            sage: Q = hanke_Matrix(2, [1, 2, 2, 1])
            sage: Q.local_density(2,3)
            '51/50'

            sage: D = hanke_diagonal_form([2,6,10,14])
            sage: D.local_density(19,19)
            '117648/117649'
        """
        _p = str(p)
        _m = str(m)
        _sig_on
        return string(Matrix_mpz_Local_Density(self.x, _p, _m))

    def local_primitive_density(self, p, m):
        """
        Returns the local primitive density at p.

        INPUT:
            p -- integer, a prime number
            m -- integer
        OUTPUT:
            A rational number, represented as a string.

        EXAMPLES:
            sage: Q = hanke_Matrix(2); Q[0,0]=10; Q[0,1] = 20;
            sage: Q.local_primitive_density(2,3)
            '1'
        """
        _p = str(p)
        _m = str(m)
        _sig_on
        return string(Matrix_mpz_Local_Primitive_Density(self.x, _p, _m))

    def is_quadratic_form(self):
        """
        Return True if this matrix defines an integer valued quadratic form, i.e.,
        the matrix is symmetric and the diagonal entries are all even.

        EXAMPLES:

            sage: Q = hanke_Matrix(2); Q[0,0]=10; Q[0,1] = 20;
            sage: Q.is_quadratic_form()
            False
            sage: hanke_Matrix(2, [2, 1, 1, 2]).is_quadratic_form()
            True
            sage: hanke_Matrix(2, [1, 2, 2, 1]).is_quadratic_form()
            False
            sage: hanke_Matrix(2, [1, 1, 1, 1]).is_quadratic_form()
            False
            sage: hanke_Matrix(2, [2, 3, 1, 2]).is_quadratic_form()
            False
        """
        return bool(Matrix_mpz_is_quadratic_form(self.x))


    def __cmp__(self, other):
        """
        Compare self to other.

        EXAMPLES:
            sage: h1 = hanke_Matrix(2, [2, 3, 1, 2])
            sage: h2 = hanke_Matrix(2, [2, 3, 1, 0])
            sage: h1 == h2
            False
            sage: h1 == h1
            True
        """
        cdef hanke_Matrix _other
        if not isinstance(other, hanke_Matrix):
            return -1
        _other = other
        if Matrix_mpz_cmp(self.x, _other.x):
            return 0
        return 1

##     def __call__(self, v):
##         """
##         Evaluate the quadratic form defined by this matrix on the given
##         list of values for the variables.

##         INPUT:
##             v -- iterable (e.g., list) with at least n arithmetic
##                  items, where n is the number of variables.

##         EXAMPLES:
##             sage: Q = hanke_Matrix(2,[0,5,5,0])
##             sage: Q([1,2])
##             10
##             sage: Q = hanke_Matrix(2,[2,5,5,2])
##             sage: Q([3,1])
##             25

##         Even iterators are allowed for v:
##             sage: Q(xrange(3,5))
##             85
##             sage: Q([3,4])
##             85
##         """
##         raise NotImplementedError


    ######################################################
    # The rest of the functional below are mainly related
    # to the quadratic form associated to this matrix.
    ######################################################

    def level(self):
        """
        The level of the quadratic form associated to this matrix.

        EXAMPLES:
            sage: Q = hanke_Matrix(2,[2,5,5,2])
            sage: Q.level()
            42
            sage: Q = hanke_Matrix(2,[1,0,0,0])
            sage: Q.level()
            Traceback (most recent call last):
            ...
            ArithmeticError: Singular matrix, so level not defined.
        """
        d = self.determinant()
        if d == 0:
            raise ArithmeticError, "Singular matrix, so level not defined."
        _sig_on
        N = int(string(Matrix_mpz_level(self.x)))
        _sig_off
        return N


    ######################################################
    # Symmetric row and column operations
    ######################################################

    def is_symmetric(self):
        """
        Return True if this is a symmetric matrix.

        EXAMPLES:
            sage: Q = hanke_Matrix(3); Q[0,0]=10; Q[0,1]=3; Q[1,2]=2; Q[1,0]=3; Q[2,1]=2;
            sage: Q.is_symmetric()
            True
            sage: Q = hanke_Matrix(2,3,range(6))
            sage: Q.is_symmetric()
            False
            sage: Q = hanke_Matrix(2,2,[1,2,3,4])
            sage: Q.is_symmetric()
            False
        """
        return bool(Matrix_mpz_is_symmetric(self.x))


    ######################################################
    # Local Normal Form
    ######################################################
    def local_normal_form(self, p):
        r"""
        Returns the local normal form at the prime p.

        INPUT:
            p -- integer, a prime number
        OUTPUT:
            ???

        EXAMPLES:
            sage: Q = hanke_Matrix(2, [10,20, 20,10])
            sage: Q
            [ 10, 20 ]
            [ 20, 10 ]
            sage: Q.determinant()
            -300
            sage: Q.local_normal_form(2)
            [ 5, 0 ]
            [ 0, -15 ]
            sage: Q.local_normal_form(3)
            [ 5, 0 ]
            [ 0, -15 ]
            sage: Q.local_normal_form(5)
            [ 5, 0 ]
            [ 0, -15 ]

        METHOD:
            Uses the algorithm described in Conway and Sloane's book
            to produce a matrix in local normal form at p, which is
            block-diagonal and ordered in increasing p-divisibility.
            When p=2 there are 1x1 and 2x2 blocks, while when p>2 this
            is diagonal.  Our representation of this is as an upper
            triangular matrix with integer entries.  (JH, 5/3/05)

        \begin{notice}

            This returns an upper triangular matrix with integer
            coefficients, though when p>2 it is diagonal.  Thus in
            general this is \emph{not} a symmetric matrix!  This was
            very convenient for computing local densities (since it
            avoided dividing by 2 all p>2), but is inconvenient for
            using in this class.  Because of this, it should be
            returned only as a matrix_mpz class, or perhaps not
            returned at all...  (JH, 5/3/05)

        \end{notice}

        """
        if not self.is_symmetric():
            raise RuntimeError, "Matrix must be symmetric."
        p = str(p)
        _sig_on
        return make_hanke_Matrix(Matrix_mpz_local_normal_form(self.x, p))


    ######################################################
    # Local Invariants
    ######################################################
    def local_diagonal(self, p):
        """
        Returns a diagonal form equivalent to Q over the p-adic numbers Q_p.

        EXAMPLES:
            sage: Q = hanke_Matrix(2, [10,20,20,0])
            sage: Q
            [ 10, 20 ]
            [ 20, 0 ]
            sage: Q.local_diagonal(2)
            [ 5, 0 ]
            [ 0, -20 ]
        """
        if not self.is_symmetric():
            raise RuntimeError, "Matrix must be symmetric."

        p = str(p)
        _sig_on
        return make_hanke_Matrix(Matrix_mpz_local_diagonal_form(self.x, p))

    def hasse_invariant(self, p):
        """
        Return the Hasse invariant at p.

        EXAMPLES:
            sage: Q = hanke_Matrix(2, [10,20,20,0])
            sage: Q
            [ 10, 20 ]
            [ 20, 0 ]
            sage: Q.hasse_invariant(2)
            1
        """
        if not self.is_symmetric():
            raise RuntimeError, "Matrix must be symmetric."
        p = str(p)
        _sig_on
        n = Matrix_mpz_hasse_invariant(self.x, p)
        _sig_off
        return n


    def is_anisotropic(self, p):
        """
        Returns True if and only if this quadratic form is anisotropic
        locally at p.

        EXAMPLES:
            sage: Q = hanke_Matrix(2, [10,20,20,0])
            sage: Q
            [ 10, 20 ]
            [ 20, 0 ]
            sage: Q.is_anisotropic(2)
            True
            sage: Q.is_isotropic(2)
            False
        """
        if not self.is_symmetric():
            raise RuntimeError, "Matrix must be symmetric."

        p = str(p)
        _sig_on
        b = Matrix_mpz_is_anisotropic(self.x, p)
        _sig_off
        return bool(b)

    def is_isotropic(self, p):
        """
        Returns True if and only if this quadratic form is isotropic
        locally at p.
        EXAMPLES:
            sage: Q = hanke_Matrix(2, [10,20,20,0])
            sage: Q
            [ 10, 20 ]
            [ 20, 0 ]
            sage: Q.is_isotropic(2)
            False
        """
        if not self.is_symmetric():
            raise RuntimeError, "Matrix must be symmetric."

        p = str(p)
        _sig_on
        b = Matrix_mpz_is_isotropic(self.x, p)
        _sig_off
        return bool(b)

    def anisotropic_primes(self):
        r"""
        Returns all primes at which this quadratic form is anisotropic.

        TODO: This is broken.  Here's debuging info, which says that you're
        trying to take the sqrt of a negative number:

        \begin{verbatim}
            was@form:~$ gdb python
            ...
            Python 2.4.1 (#2, Jul 12 2005, 13:09:38)
            [GCC 4.0.1 (Debian 4.0.1-1)] on linux2
            Type "help", "copyright", "credits" or "license" for more information.
            >>.> from hanke import hanke_Matrix;  Q = hanke_Matrix(2, [10,20,20,5]);  Q.anisotropic_primes()

            Program received signal SIGFPE, Arithmetic exception.
            [Switching to Thread -1209675008 (LWP 29424)]
            0xb7ca797b in __gmp_exception (error_bit=4) at errno.c:40
            40        __gmp_junk = 10 / __gmp_0;
            (gdb) bt
            #0  0xb7ca797b in __gmp_exception (error_bit=4) at errno.c:40
            #1  0xb7ca79ab in __gmp_sqrt_of_negative () at errno.c:51
            #2  0xb7cbb690 in __gmpz_sqrt (root=0xbfffeff8, op=Variable "op" is not available.
            ) at sqrt.c:40
            #3  0xb7d32167 in PrimeDivisors (m=@0xbffff080) at hanke/GMP_class_extras/mpz_class_extras.cc:368
            #4  0xb7d31500 in Matrix_mpz::AnisotropicPrimes (this=0x819cd58) at hanke/Matrix_mpz/Local_Invariants.cc:199
            #5  0xb7d08e35 in Matrix_mpz_anisotropic_primes (x=0x819cd58) at hanke/wrap.cc:178
        \end{verbatim}

        """
        if not self.is_symmetric():
            raise RuntimeError, "Matrix must be symmetric."
        _sig_on
        return string(Matrix_mpz_anisotropic_primes(self.x))


    ######################################################
    # Local Constants
    ######################################################
    def local_constant(self, p, T):
        """
        This computes the local constants C_p(T) required for finding
        an explicit lower bound for the growth of local densities
        within a p-squareclass $T*p^(2i)$.  They are precisely defined
        in [Ha, top of p370].

        TODO: What is Ha???

        EXAMPLES:
            sage: D = hanke_diagonal_form([2,6,10,14])
            sage: print D
            [ 2, 0, 0, 0 ]
            [ 0, 6, 0, 0 ]
            [ 0, 0, 10, 0 ]
            [ 0, 0, 0, 14 ]

            sage: D.local_constant(3,10)         # not tested
            '15/2'

            sage: D.local_constant(3,30)         # not tested
            '153/20'

        TODO: How long does this take?  I tried running the above
              commands and similar ones and didn't get any output (at least
              not in a minute).

        """
        if not self.is_symmetric():
            raise RuntimeError, "Matrix must be symmetric."

        p = str(p)
        T = str(T)
        _sig_on
        return string(Matrix_mpz_local_constants(self.x, p, T))

    def is_stable(self, p, T):
        r"""
        Checks if T is p-stable for this quadratic form, meaning that
        the growth of local densities as grow by powers of $p^2$ from
        $T*p^2$ onwards is independent of the level.  This is
        precisely defined in [Ha, Definitons 3.6 and 5.2].

        TODO: first sentence of above description needs editing -- it
        sounds weird.

        \begin{notice} This function depends on the local density
        routines, which only have been tested when n=3 and n=4.  Also,
        the local density routines depended on being passed a local
        normal form of Q at that prime, which was causing major issues
        and has now been fixed.  So with all of this, the n=3 and n=4
        cases should work, and the others should work better. =) (JH,
        5/3/05) \end{notice}

        It's still very slow when p=2, even for the $x^2 + 3y^2 + 5z^2
        + 7w^2$ form...  so possibly something is wrong, however $p=3$
        it works very well. =) There's also some spurious output to get
        rid of...  (JH, 5/3/05)

        EXAMPLES:
            sage: Q=hanke_diagonal_form([2,2])
            sage: Q.is_stable(2,3)
            1

            sage: D = hanke_diagonal_form([2,6,10,14])
            sage: print D
            [ 2, 0, 0, 0 ]
            [ 0, 6, 0, 0 ]
            [ 0, 0, 10, 0 ]
            [ 0, 0, 0, 14 ]
            sage: D.is_stable(3,90)      # TODO: the program used to give False??
            True
            sage: D.is_stable(3,270)     # TODO: the program used to give False??
            True
            sage: D.is_stable(3,810)
            True
            sage: D.is_stable(3,2430)
            True

        """
        if not self.is_symmetric():
            raise RuntimeError, "Matrix must be symmetric."

        p = str(p)
        T = str(T)
        _sig_on
        b = Matrix_mpz_is_stable(self.x, p, T)
        _sig_off
        return bool(b)

cdef make_hanke_Matrix(Matrix_mpz* m):
    cdef hanke_Matrix h
    _sig_off
    h = hanke_Matrix(_INIT)
    h._nrows = Matrix_mpz_nrows(m)
    h._ncols = Matrix_mpz_ncols(m)
    h.x = m
    return h

def hanke_diagonal_form(v):
    n = len(v)
    A = hanke_Matrix(n)
    A.set_diagonal(v)
    return A


