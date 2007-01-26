include "../ext/interrupt.pxi"
include "../ext/stdsage.pxi"
include "../ext/gmp.pxi"


from sage.rings.integer_ring import ZZ
from sage.rings.arith import next_prime

# should I have mod_int versions of these functions?
# c_inverse_mod_longlong modular inverse used exactly once in _extend_moduli_to_height
from sage.ext.arith cimport arith_llong
cdef arith_llong ai
ai = arith_llong()

cdef class MultiModularBasis_base:

    def __new__(self, height):
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


    cdef int _extend_moduli_to_height(self, mpz_t height) except -1:
        if mpz_cmp(height, self.partial_products[self.n-1]) <= 0: return self.n
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

    cdef int _extend_moduli_to_count(self, int count) except -1:
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

    cdef void _refresh_products(self, int start):
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
        if start == 0:
            start = 1 # first one is trivial, never used
        for i from start <= i < self.n:
            self.C[i] = ai.c_inverse_mod_longlong(mpz_fdiv_ui(self.partial_products[i-1], self.moduli[i]), self.moduli[i])

    cdef int min_moduli_count(self, mpz_t height) except -1:
        self._extend_moduli_to_height(height)

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
        cdef int i
        cdef mod_int* m
        m = self.moduli + offset
        for i from 0 <= i < len:
            b[i] = mpz_fdiv_ui(z, m[i])
        return 0

    cdef int mpz_reduce_vec_tail(self, mpz_t* z, mod_int** b, int vn, int offset, int len) except -1:
        cdef int i, j
        cdef mod_int* m
        m = self.moduli + offset
        for i from 0 <= i < len:
            mi = m[i]
            for j from 0 <= j < vn:
                b[i][j] = mpz_fdiv_ui(z[j], mi)
        return 0

    cdef int mpz_crt_tail(self, mpz_t z, mod_int* b, int offset, int len) except -1:
        # Garner's Algorithm
        # z is not yet initalized
        cdef int i, s
        cdef mpz_t u
        cdef mod_int* m
        m = self.moduli + offset
        mpz_init(u)
        if offset == 0:
            s = 1
            mpz_init_set_ui(z, b[0])
        else:
            s = 0
        for i from s <= i < len:
            mpz_set_ui(u, ((b[i] + m[i] - mpz_fdiv_ui(z, m[i])) * self.C[i]) % m[i])
            mpz_mul(u, u, self.partial_products[i-1])
            mpz_add(z, z, u)
        mpz_clear(u)
        return 0

    cdef int mpz_crt_vec_tail(self, mpz_t* z, mod_int** b, int vc, int offset, int len) except -1:
        # Garner's Algorithm
        # z allocated but not initalized
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
            if offset == 0:
                mpz_init_set_ui(z[j], b[0][j])
            for i from s <= i < len:
                mpz_set_ui(u, ((b[i][j] + m[i] - mpz_fdiv_ui(z[j], m[i])) * self.C[i]) % m[i]) # u = ((b_i - z) * C_i) % m_i
                mpz_mul(u, u, self.partial_products[i-1])
                mpz_add(z[j], z[j], u)
        return 0

    def crt(self, b):
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
        cdef int i
        i = ix
        if i < 0 or i >= self.n:
            raise IndexError, "index out of range"
        return self.moduli[i]

    def __repr__(self):
        return str(type(self))+str(self.list())


cdef class MultiModularBasis(MultiModularBasis_base):

    def __init__(self, height):
        cdef Integer h
        h = ZZ(height)
        self._extend_moduli_to_height(h.value)

    cdef int mpz_reduce(self, mpz_t z, mod_int* b) except -1:
        self.mpz_reduce_tail(z, b, 0, self.n)

    cdef int mpz_reduce_vec(self, mpz_t* z, mod_int** b, int vn) except -1:
        self.mpz_reduce_vec_tail(z, b, vn, 0, self.n)

    cdef int mpz_crt(self, mpz_t z, mod_int* b) except -1:
        self.mpz_crt_tail(z, b, 0, self.n)

    cdef int mpz_crt_vec(self, mpz_t* z, mod_int** b, int vn) except -1:
        self.mpz_crt_vec_tail(z, b, vn, 0, self.n)


cdef class MutableMultiModularBasis(MultiModularBasis):

    def __init__(self, height):
        MultiModularBasis.__init__(self, height)

    cdef mod_int last_prime(self):
        if self.moduli[self.n-1] > self.__last_prime:
            self.__last_prime = self.moduli[self.n-1]
        return self.__last_prime

    def next_prime(self):
        return self.next_prime_c()

    cdef mod_int next_prime_c(self) except -1:
        self._extend_moduli_to_count(self.n+1)
        return self.moduli[self.n-1]

    def replace_prime(self, ix):
        return self.replace_prime_c(ix)

    cdef mod_int replace_prime_c(self, int ix) except -1:
        cdef int i
        cdef Integer z
        cdef mod_int bad_p, new_p

        if ix < 0 or ix >= self.n:
            raise IndexError, "index out of range"

        bad_b = self.moduli[ix]
        z = next_prime(self.last_prime())
        new_p = z
        self.__last_prime = new_p
        self.moduli[ix] = new_p

        if ix > 0:
            mpz_mul(self.partial_products[ix], self.partial_products[ix-1], z.value)
        else:
            mpz_set(self.partial_products[ix], z.value)

        for i from ix < i < self.n:
            mpz_set_ui(z.value, self.moduli[i])
            mpz_mul(self.partial_products[i], self.partial_products[i-1], z.value)

        self._refresh_products(ix)
        self._refresh_precomputations(ix)
        return new_p
