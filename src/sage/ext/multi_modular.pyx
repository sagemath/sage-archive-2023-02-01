include "../ext/interrupt.pxi"
include "../ext/stdsage.pxi"
include "../ext/gmp.pxi"


from sage.rings.integer_ring import ZZ
from sage.rings.arith import next_prime

from sage.rings.integer import Integer
from sage.rings.integer cimport Integer

# should I have mod_int versions of these functions?
# c_inverse_mod_longlong modular inverse used exactly once in _extend_moduli_list
from sage.ext.arith cimport arith_llong
cdef arith_llong ai
ai = arith_llong()

cdef class MultiModularBasis:

    def __init__(self, start_prime=START_PRIME_MAX):
        cdef mod_int p
        p = next_prime(start_prime-1)

        self.moduli_count = 1

        self.moduli = <mod_int*>sage_malloc(sizeof(mod_int))
        self.moduli_partial_product = <mpz_t*>sage_malloc(sizeof(mpz_t))
        self.C = <mod_int*>sage_malloc(sizeof(mod_int))

        if self.moduli == NULL or self.moduli_partial_product == NULL or self.C == NULL:
            raise MemoryError, "out of memory allocating multi-modular prime list"

        self.moduli[0] = p
        mpz_init_set_ui(self.moduli_partial_product[0], p)
        self.C[0] = 1

    def __dealloc__(self):
        sage_free(self.moduli)
        sage_free(self.moduli_partial_product)
        sage_free(self.C)


    cdef int _extend_moduli_list(self, mpz_t height) except -1:
        if mpz_cmp(height, self.moduli_partial_product[self.moduli_count-1]) <= 0: return self.moduli_count
        cdef int i
        new_moduli = []
        new_partial_products = []
        new_C = []
        cdef Integer p, M
        p = PY_NEW(Integer)
        mpz_set_ui(p.value, self.moduli[self.moduli_count-1])
        M = PY_NEW(Integer)
        mpz_set(M.value, self.moduli_partial_product[self.moduli_count-1])
        while mpz_cmp(height, M.value) > 0:
            p = next_prime(p)
            new_moduli.append(p)
            new_C.append(ai.c_inverse_mod_longlong(mpz_fdiv_ui(M.value, p), p))
#            new_C.append((mod(M, p)**(-1)).lift())
            M *= p
            new_partial_products.append(M)

        cdef int new_count
        new_count = self.moduli_count + len(new_moduli)
        self.moduli = <mod_int*>sage_realloc(self.moduli, sizeof(mod_int)*new_count)
        self.moduli_partial_product = <mpz_t*>sage_realloc(self.moduli_partial_product, sizeof(mpz_t)*new_count)
        self.C = <mod_int*>sage_realloc(self.C, sizeof(mod_int)*new_count)
        if self.moduli == NULL or self.moduli_partial_product == NULL or self.C == NULL:
            raise MemoryError, "out of memory allocating multi-modular prime list"

        for i from self.moduli_count <= i < new_count:
            self.moduli[i] = new_moduli[i-self.moduli_count]
            self.C[i] = new_C[i-self.moduli_count]
            mpz_init_set(self.moduli_partial_product[i], (<Integer>new_partial_products[i-self.moduli_count]).value)
        self.moduli_count = new_count
        return new_count


    cdef int moduli_list_c(self, mod_int** moduli, mpz_t height) except -1:
        self._extend_moduli_list(height)

        cdef int count
        count = self.moduli_count * mpz_sizeinbase(height, 2) / mpz_sizeinbase(self.moduli_partial_product[self.moduli_count-1], 2) # an estimate
        count = max(min(count, self.moduli_count), 1)
        while count > 1 and mpz_cmp(height, self.moduli_partial_product[count-1]) < 0:
           count -= 1
        while mpz_cmp(height, self.moduli_partial_product[count-1]) > 0:
           count += 1

        moduli[0] = <mod_int*>sage_malloc(sizeof(mod_int)*count)
        if moduli[0] == NULL:
            raise MemoryError, "out of memory allocating multi-modular prime list"
        memcpy(moduli[0], self.moduli, sizeof(mod_int)*count)
        return count

    cdef int mpz_crt(self, mpz_t z, mod_int* b, int n) except -1:
        # Garner's Algorithm
        # z is not yet initalized
        cdef int i
        cdef mpz_t u
        mpz_init(u)
        mpz_init_set_ui(z, b[0])
        for i from 1 <= i < n:
            mpz_set_ui(u, ((b[i] + self.moduli[i] - mpz_fdiv_ui(z, self.moduli[i])) * self.C[i]) % self.moduli[i])
            mpz_mul(u, u, self.moduli_partial_product[i-1])
            mpz_add(z, z, u)
        return 0

    cdef int mpz_crt_vec(self, mpz_t* z, mod_int** b, int n, int vc) except -1:
        # Garner's Algorithm
        # z allocated but not initalized
        cdef int i, j
        cdef mpz_t u
        mpz_init(u)
        for j from 0 <= j < vc:
            mpz_init_set_ui(z[j], b[0][j])
            for i from 1 <= i < n:
                mpz_set_ui(u, ((b[i][j] + self.moduli[i] - mpz_fdiv_ui(z[j], self.moduli[i])) * self.C[i]) % self.moduli[i]) # u = ((b_i - z) * C_i) % m_i
                mpz_mul(u, u, self.moduli_partial_product[i-1])
                mpz_add(z[j], z[j], u)
        return 0

    def moduli_list(self, height):
        cdef Integer h
        h = ZZ(height)
        cdef int i, n
        cdef mod_int* moduli
        n = self.moduli_list_c(&moduli, h.value)
        m = [ZZ(moduli[i]) for i from 0 <= i < n]
        sage_free(moduli)
        return m

    def crt(self, b):
        cdef int i, n
        n = len(b)
        if n > self.moduli_count:
            raise IndexError, "beyond bound for multi-modular prime list"
        cdef mod_int* bs
        bs = <mod_int*>sage_malloc(sizeof(mod_int)*n)
        if bs == NULL:
            raise MemoryError, "out of memory allocating multi-modular prime list"
        for i from 0 <= i < n:
            bs[i] = b[i]
        cdef Integer z
        z = PY_NEW(Integer)
        self.mpz_crt(z.value, bs, n)
        sage_free(bs)
        return z

    def precomputations(self):
        cdef int i
        return [ZZ(self.C[i]) for i from 0 <= i < self.moduli_count]

    def partial_product(self, n):
        if n > self.moduli_count:
            raise IndexError, "beyond bound for multi-modular prime list"
        cdef Integer z
        z = PY_NEW(Integer)
        mpz_init_set(z.value, self.moduli_partial_product[self.moduli_count-1])
        return z
