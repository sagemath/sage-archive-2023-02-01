"""
Utility classes for multi-modular algorithms.
"""

######################################################################
#       Copyright (C) 2006 William Stein
#
#  Distributed under the terms of the GNU General Public License (GPL)
#
#  The full text of the GPL is available at:
#                  http://www.gnu.org/licenses/
######################################################################


include "../ext/interrupt.pxi"
include "../ext/stdsage.pxi"
include "../ext/gmp.pxi"


from sage.rings.integer_ring import ZZ
from sage.rings.arith import next_prime # should I just use probable primes?

# should I have mod_int versions of these functions?
# c_inverse_mod_longlong modular inverse used exactly once in _refresh_precomputations
from sage.ext.arith cimport arith_llong
cdef arith_llong ai
ai = arith_llong()

# TODO: have one global instance for sharing, copy for MutableMultiModularBasis

cdef class MultiModularBasis_base:
    r"""
    This class stores a list of machine-sized prime numbers,
    and can do reduction and Chinese Remainder Theorem lifting
    modulo these primes.

    Lifting implemented via Garner's algorithm, which has the advantage
    that all reductions are word-sized. For each $i$ precompute

       $\prod_j=1^{i-1} m_j$ and $\prod_j=1^{i-1} m_j^{-1} (mod m_i)$
    """

    def __new__(self, height):
        r"""
        Allocate the space for the moduli and precomputation lists
        and initalize the first element of that list.
        """
        cdef mod_int p
        p = next_prime(START_PRIME_MAX-1)

        self.n = 1

        self.moduli = <mod_int*>sage_malloc(sizeof(mod_int))
        self.partial_products = <mpz_t*>sage_malloc(sizeof(mpz_t))
        self.C = <mod_int*>sage_malloc(sizeof(mod_int))

        if self.moduli == NULL or self.partial_products == NULL or self.C == NULL:
            raise MemoryError, "out of memory allocating multi-modular prime list"

        self.moduli[0] = p
        mpz_init_set_ui(self.partial_products[0], p)
        self.C[0] = 1

    def __dealloc__(self):
        sage_free(self.moduli)
        sage_free(self.partial_products)
        sage_free(self.C)

    def __init__(self, height):
        r"""
        Initialize a multi-modular basis and perform precomputations.

        INPUT:
            height -- as integer
                          determines how many primes are computed
                          (their product must be at least height)
                      as list
                          a list of prime moduli to start with
        """
        cdef int i
        if isinstance(height, list):
            m = height
            self.n = len(m)
            self.moduli = <mod_int*>sage_realloc(self.moduli, sizeof(mod_int) * self.n)
            self.partial_products = <mpz_t*>sage_realloc(self.partial_products, sizeof(mpz_t) * self.n)
            self.C = <mod_int*>sage_realloc(self.C, sizeof(mod_int) * self.n)
            if self.moduli == NULL or self.partial_products == NULL or self.C == NULL:
                raise MemoryError, "out of memory allocating multi-modular prime list"
            for i from 0 <= i < len(m):
                if m[i] > MOD_INT_MAX:
                    raise ValueError, "given modulus cannot be manipulated in a single machine word"
                self.moduli[i] = m[i]
            for i from 1 <= i < len(m):
                mpz_init(self.partial_products[i])
            self._refresh_products(0)
            self._refresh_precomputations(0)
        else:
            self._extend_moduli_to_height(height)


    def _extend_moduli_to_height(self, height):
        cdef Integer h
        h = ZZ(height)
        self._extend_moduli_to_height_c(h.value)

    cdef int _extend_moduli_to_height_c(self, mpz_t height) except -1:
        r"""
        Expand the list of primes and perform precomputations.

        INPUT:
            height -- determines how many primes are computed
                      (their product must be at least height)
        """
        if mpz_cmp(height, self.partial_products[self.n-1]) <= 0:
            return self.n
        cdef int i
        new_moduli = []
        new_partial_products = []
        cdef Integer p, M
        M = PY_NEW(Integer)
        mpz_set(M.value, self.partial_products[self.n-1])
        p = ZZ(self.last_prime())
        while mpz_cmp(height, M.value) > 0:
            p = next_prime(p)
            new_moduli.append(p)
            M *= p
            new_partial_products.append(M)

        cdef int new_count, old_count
        old_count = self.n
        new_count = self.n + len(new_moduli)

        self.moduli = <mod_int*>sage_realloc(self.moduli, sizeof(mod_int)*new_count)
        self.partial_products = <mpz_t*>sage_realloc(self.partial_products, sizeof(mpz_t)*new_count)
        self.C = <mod_int*>sage_realloc(self.C, sizeof(mod_int)*new_count)
        if self.moduli == NULL or self.partial_products == NULL or self.C == NULL:
            raise MemoryError, "out of memory allocating multi-modular prime list"

        for i from self.n <= i < new_count:
            self.moduli[i] = new_moduli[i-self.n]
            mpz_init_set(self.partial_products[i], (<Integer>new_partial_products[i-self.n]).value)

        self.n = new_count
        self._refresh_precomputations(old_count)
        return new_count

    def _extend_moduli_to_count(self, int count):
        r"""
        Expand the list of primes and perform precomputations.

        INPUT:
            count -- the minimum number of moduli in the resulting list
        """
        if count <= self.n:
            return self.n
        self.moduli = <mod_int*>sage_realloc(self.moduli, sizeof(mod_int) * count)
        self.partial_products = <mpz_t*>sage_realloc(self.partial_products, sizeof(mpz_t) * count)
        self.C = <mod_int*>sage_realloc(self.C, sizeof(mod_int) * count)
        if self.moduli == NULL or self.partial_products == NULL or self.C == NULL:
            raise MemoryError, "out of memory allocating multi-modular prime list"

        cdef int i
        cdef Integer p
        p = ZZ(self.last_prime())
        for i from self.n <= i < count:
            p = next_prime(p)
            self.moduli[i] = p
            mpz_init(self.partial_products[i])
            mpz_mul(self.partial_products[i], self.partial_products[i-1], p.value)
        old_count = self.n
        self.n = count
        self._refresh_precomputations(old_count)
        return count

    def _extend_moduli(self, count):
        self._extend_moduli_to_count(self.n + count)

    cdef void _refresh_products(self, int start):
        r"""
        Compute and store $\prod_j=1^{i-1} m_j$ of i > start.
        """
        cdef mpz_t z
        mpz_init(z)
        if start == 0:
            mpz_set_ui(self.partial_products[0], self.moduli[0])
            start += 1
        for i from start <= i < self.n:
            mpz_set_ui(z, self.moduli[i])
            mpz_mul(self.partial_products[i], self.partial_products[i-1], z)
        mpz_clear(z)


    cdef void _refresh_precomputations(self, int start):
        r"""
        Compute and store $\prod_j=1^{i-1} m_j^{-1} (mod m_i)$ of i >= start.
        """
        if start == 0:
            start = 1 # first one is trivial, never used
        for i from start <= i < self.n:
            self.C[i] = ai.c_inverse_mod_longlong(mpz_fdiv_ui(self.partial_products[i-1], self.moduli[i]), self.moduli[i])

    cdef int min_moduli_count(self, mpz_t height) except -1:
        r"""
        Compute the minimum number of primes needed to uniquely determin
        an integer mod height.
        """
        self._extend_moduli_to_height_c(height)

        cdef int count
        count = self.n * mpz_sizeinbase(height, 2) / mpz_sizeinbase(self.partial_products[self.n-1], 2) # an estimate
        count = max(min(count, self.n), 1)
        while count > 1 and mpz_cmp(height, self.partial_products[count-1]) < 0:
           count -= 1
        while mpz_cmp(height, self.partial_products[count-1]) > 0:
           count += 1

        return count

    cdef int moduli_list_c(self, mod_int** moduli):
        memcpy(moduli[0], self.moduli, sizeof(mod_int)*self.n)
        return self.n

    cdef mod_int last_prime(self):
        return self.moduli[self.n-1]

    cdef int mpz_reduce_tail(self, mpz_t z, mod_int* b, int offset, int len) except -1:
        r"""
        Performs reduction mod $m_i$ for offset <= i < len

        b[i] = z mod $m_{i+offset}$ for 0 <= i < len

        INPUT:
            z -- the integer being reduced
            b -- array to hold the reductions mod each m_i.
                 It MUST be allocated and have length at least len
            offset -- first prime in list to reduce against
            len    -- number of primes in list to reduce against
        """
        cdef int i
        cdef mod_int* m
        m = self.moduli + offset
        for i from 0 <= i < len:
            b[i] = mpz_fdiv_ui(z, m[i])
        return 0

    cdef int mpz_reduce_vec_tail(self, mpz_t* z, mod_int** b, int vn, int offset, int len) except -1:
        r"""
        Performs reduction mod $m_i$ for offset <= i < len

        b[i][j] = z[j] mod $m_{i+offset}$ for 0 <= i < len

        INPUT:
            z  -- an array of integers being reduced
            b  -- array to hold the reductions mod each m_i.
                 It MUST be fully allocated and each
                 have length at least len
            vn -- length of z and each b[i]
            offset -- first prime in list to reduce against
            len    -- number of primes in list to reduce against
        """
        cdef int i, j
        cdef mod_int* m
        m = self.moduli + offset
        for i from 0 <= i < len:
            mi = m[i]
            for j from 0 <= j < vn:
                b[i][j] = mpz_fdiv_ui(z[j], mi)
        return 0

    cdef int mpz_crt_tail(self, mpz_t z, mod_int* b, int offset, int len) except -1:
        r"""
        Calculate lift mod $\prod_{i=0}^{offset+len-1} m_i$.

        z = b[i] mod $m_{i+offset}$ for 0 <= i < len

        In the case that offset > 0,
        z remains unchanged mod $\prod_{i=0}^{offset-1} m_i$

        INPUT:
            z  -- a placeholder for the constructed integer
                  z MUST be initalized IF and ONLY IF offset > 0
            b  -- array holding the reductions mod each m_i.
                  It MUST have length at least len
            offset -- first prime in list to reduce against
            len    -- number of primes in list to reduce against
        """
        cdef int i, s
        cdef mpz_t u
        cdef mod_int* m
        m = self.moduli + offset
        mpz_init(u)
        if offset == 0:
            s = 1
            mpz_init_set_ui(z, b[0])
            if b[0] == 0:
                while s < len and b[s] == 0: # fast forward to first non-zero
                    s += 1
        else:
            s = 0
        for i from s <= i < len:
            mpz_set_ui(u, ((b[i] + m[i] - mpz_fdiv_ui(z, m[i])) * self.C[i]) % m[i])
            mpz_mul(u, u, self.partial_products[i-1])
            mpz_add(z, z, u)
        mpz_clear(u)
        return 0

    cdef int mpz_crt_vec_tail(self, mpz_t* z, mod_int** b, int vc, int offset, int len) except -1:
        r"""
        Calculate lift mod $\prod_{i=0}^{offset+len-1} m_i$.

        z[j] = b[i][j] mod $m_{i+offset}$ for 0 <= i < len

        In the case that offset > 0,
        z[j] remains unchanged mod $\prod_{i=0}^{offset-1} m_i$

        INPUT:
            z  -- a placeholder for the constructed integers
                  z MUST be allocated and have length at least vc
                  z[j] MUST be initalized IF and ONLY IF offset > 0
            b  -- array holding the reductions mod each m_i.
                  MUST have length at least len
            vn -- length of z and each b[i]
            offset -- first prime in list to reduce against
            len    -- number of primes in list to reduce against
        """
        cdef int i, j
        cdef mpz_t u
        cdef mod_int* m
        m = self.moduli + offset
        mpz_init(u)
        if offset == 0:
            s = 1
        else:
            s = 0
        for j from 0 <= j < vc:
            i = s
            if offset == 0:
                mpz_init_set_ui(z[j], b[0][j])
                if b[0][j] == 0:
                    while i < len and b[i][j] == 0: # fast forward to first non-zero
                        i += 1
            while i < len:
                mpz_set_ui(u, ((b[i][j] + m[i] - mpz_fdiv_ui(z[j], m[i])) * self.C[i]) % m[i]) # u = ((b_i - z) * C_i) % m_i
                mpz_mul(u, u, self.partial_products[i-1])
                mpz_add(z[j], z[j], u)
                i += 1
        mpz_clear(u)
        return 0

    def crt(self, b):
        r"""
        Calculate lift mod $\prod_{i=0}^{len(b)-1} m_i$.

        In the case that offset > 0,
        z[j] remains unchanged mod $\prod_{i=0}^{offset-1} m_i$

        INPUT:
            b  -- a list of length at most self.n
        OUTPUT:
            Integer z where z = b[i] mod $m_i$ for 0 <= i < len(b)
        """
        cdef int i, n
        n = len(b)
        if n > self.n:
            raise IndexError, "beyond bound for multi-modular prime list"
        cdef mod_int* bs
        bs = <mod_int*>sage_malloc(sizeof(mod_int)*n)
        if bs == NULL:
            raise MemoryError, "out of memory allocating multi-modular prime list"
        for i from 0 <= i < n:
            bs[i] = b[i]
        cdef Integer z
        z = PY_NEW(Integer)
        self.mpz_crt_tail(z.value, bs, 0, n)
        sage_free(bs)
        return z

    def precomputation_list(self):
        cdef int i
        return [ZZ(self.C[i]) for i from 0 <= i < self.n]

    def partial_product(self, n):
        if n > self.n:
            raise IndexError, "beyond bound for multi-modular prime list"
        cdef Integer z
        z = PY_NEW(Integer)
        mpz_init_set(z.value, self.partial_products[n-1])
        return z

    def prod(self):
        cdef Integer z
        z = PY_NEW(Integer)
        mpz_init_set(z.value, self.partial_products[self.n-1])
        return z

    def list(self):
        cdef int i
        cdef mod_int* moduli
        moduli = <mod_int*>sage_malloc(sizeof(mod_int) * self.n)
        n = self.moduli_list_c(&moduli)
        m = [ZZ(moduli[i]) for i from 0 <= i < n]
        sage_free(moduli)
        return m

    def __len__(self):
        return self.n

    def __iter__(self):
        return self.list().__iter__()

    def __getitem__(self, ix):

        if isinstance(ix, slice):
            return type(self)(self.list()[ix])

        cdef int i
        i = ix
        if i < 0 or i >= self.n:
            raise IndexError, "index out of range"
        return self.moduli[i]

    def __repr__(self):
        return str(type(self))+str(self.list())


cdef class MultiModularBasis(MultiModularBasis_base):
    """
    Class used for storing a MultiModular bases of a fixed length.
    """


    cdef int mpz_reduce(self, mpz_t z, mod_int* b) except -1:
        r"""
        Performs reduction mod $m_i$ for each modulus $m_i$

        b[i] = z mod $m_i$ for 0 <= i < len(self)

        INPUT:
            z -- the integer being reduced
            b -- array to hold the reductions mod each m_i.
                 It MUST be allocated and have length at least len
        """
        self.mpz_reduce_tail(z, b, 0, self.n)

    cdef int mpz_reduce_vec(self, mpz_t* z, mod_int** b, int vn) except -1:
        r"""
        Performs reduction mod $m_i$ for each modulus $m_i$

        b[i][j] = z[j] mod $m_i$ for 0 <= i < len(self)

        INPUT:
            z  -- an array of integers being reduced
            b  -- array to hold the reductions mod each m_i.
                 It MUST be fully allocated and each
                 have length at least len
            vn -- length of z and each b[i]
        """
        self.mpz_reduce_vec_tail(z, b, vn, 0, self.n)

    cdef int mpz_crt(self, mpz_t z, mod_int* b) except -1:
        r"""
        Calculate lift mod $\prod m_i$.

        z = b[i] mod $m_{i+offset}$ for 0 <= i < len(self)

        INPUT:
            z  -- a placeholder for the constructed integer
                  z MUST NOT be initalized
            b  -- array holding the reductions mod each $m_i$.
                  It MUST have length at least len(self)
        """
        self.mpz_crt_tail(z, b, 0, self.n)

    cdef int mpz_crt_vec(self, mpz_t* z, mod_int** b, int vn) except -1:
        r"""
        Calculate lift mod $\prod m_i$.

        z[j] = b[i][j] mod $m_i$ for 0 <= i < len(self)

        INPUT:
            z  -- a placeholder for the constructed integers
                  z MUST be allocated and have length at least vn,
                  but each z[j] MUST NOT be initalized
            b  -- array holding the reductions mod each $m_i$.
                  It MUST have length at least len(self)
            vn -- length of z and each b[i]
        """
        self.mpz_crt_vec_tail(z, b, vn, 0, self.n)


cdef class MutableMultiModularBasis(MultiModularBasis):
    """
    Class used for performing multi-modular methods,
    with the possiblity of removing bad primes.
    """

    cdef mod_int last_prime(self):
        if self.moduli[self.n-1] > self.__last_prime:
            self.__last_prime = self.moduli[self.n-1]
        return self.__last_prime

    def next_prime(self):
        return self.next_prime_c()

    cdef mod_int next_prime_c(self) except -1:
        self._extend_moduli(1)
        return self.moduli[self.n-1]

    def replace_prime(self, ix):
        return self.replace_prime_c(ix)

    cdef mod_int replace_prime_c(self, int ix) except -1:
        r"""
        Replace $m_{ix} in the list of moduli with a new
        prime number greater than all others in the list,
        and recompute all precomputations.

        INPUT:
            ix -- index into list of moduli

        OUTPUT:
            p -- the new prime modulus
        """
        cdef int i
        cdef mod_int new_p

        if ix < 0 or ix >= self.n:
            raise IndexError, "index out of range"

        new_p = next_prime(self.last_prime())
        self.__last_prime = new_p
        self.moduli[ix] = new_p

        self._refresh_products(ix)
        self._refresh_precomputations(ix)
        return new_p
