"""
Fast computation of combinatorial functions (Cython + mpz).

Currently implemented:

- Stirling numbers of the second kind
- iterators for set partitions

AUTHORS:

- Fredrik Johansson (2010-10): Stirling numbers of second kind
- Martin Rubey and Travis Scrimshaw (2018): iterators for set partitions
"""

cimport cython

from cysignals.memory cimport check_allocarray, sig_free

from sage.libs.gmp.all cimport *
from sage.rings.integer cimport Integer

cdef void mpz_addmul_alt(mpz_t s, mpz_t t, mpz_t u, unsigned long parity):
    """
    Set s = s + t*u * (-1)^parity
    """
    if parity & 1:
        mpz_submul(s, t, u)
    else:
        mpz_addmul(s, t, u)


cdef mpz_stirling_s2(mpz_t s, unsigned long n, unsigned long k):
    """
    Set s = S(n,k) where S(n,k) denotes a Stirling number of the
    second kind.

    Algorithm: S(n,k) = (sum_{j=0}^k (-1)^(k-j) C(k,j) j^n) / k!

    TODO: compute S(n,k) efficiently for large n when n-k is small
    (e.g. when k > 20 and n-k < 20)
    """
    cdef mpz_t t, u
    cdef mpz_t *bc
    cdef unsigned long j, max_bc
    # Some important special cases
    if k+1 >= n:
        # Upper triangle of n\k table
        if k > n:
            mpz_set_ui(s, 0)
        elif n == k:
            mpz_set_ui(s, 1)
        elif k+1 == n:
            # S(n,n-1) = C(n,2)
            mpz_set_ui(s, n)
            mpz_mul_ui(s, s, n-1)
            mpz_tdiv_q_2exp(s, s, 1)
    elif k <= 2:
        # Leftmost three columns of n\k table
        if k == 0:
            mpz_set_ui(s, 0)
        elif k == 1:
            mpz_set_ui(s, 1)
        elif k == 2:
            # 2^(n-1)-1
            mpz_set_ui(s, 1)
            mpz_mul_2exp(s, s, n-1)
            mpz_sub_ui(s, s, 1)
    # Direct sequential evaluation of the sum
    elif n < 200:
        mpz_init(t)
        mpz_init(u)
        mpz_set_ui(t, 1)
        mpz_set_ui(s, 0)
        for j in range(1, k//2+1):
            mpz_mul_ui(t, t, k+1-j)
            mpz_tdiv_q_ui(t, t, j)
            mpz_set_ui(u, j)
            mpz_pow_ui(u, u, n)
            mpz_addmul_alt(s, t, u, k+j)
            if 2*j != k:
                # Use the fact that C(k,j) = C(k,k-j)
                mpz_set_ui(u, k-j)
                mpz_pow_ui(u, u, n)
                mpz_addmul_alt(s, t, u, j)
        # Last term not included because loop starts from 1
        mpz_set_ui(u, k)
        mpz_pow_ui(u, u, n)
        mpz_add(s, s, u)
        mpz_fac_ui(t, k)
        mpz_tdiv_q(s, s, t)
        mpz_clear(t)
        mpz_clear(u)
    # Only compute odd powers, saving about half of the time for large n.
    # We need to precompute binomial coefficients since they will be accessed
    # out of order, adding overhead that makes this slower for small n.
    else:
        mpz_init(t)
        mpz_init(u)
        max_bc = (k+1)//2
        bc = <mpz_t*> check_allocarray(max_bc+1, sizeof(mpz_t))
        mpz_init_set_ui(bc[0], 1)
        for j in range(1, max_bc+1):
            mpz_init_set(bc[j], bc[j-1])
            mpz_mul_ui(bc[j], bc[j], k+1-j)
            mpz_tdiv_q_ui(bc[j], bc[j], j)
        mpz_set_ui(s, 0)
        for j in range(1, k+1, 2):
            mpz_set_ui(u, j)
            mpz_pow_ui(u, u, n)
            # Process each 2^p * j, where j is odd
            while True:
                if j > max_bc:
                    mpz_addmul_alt(s, bc[k-j], u, k+j)
                else:
                    mpz_addmul_alt(s, bc[j], u, k+j)
                j *= 2
                if j > k:
                    break
                mpz_mul_2exp(u, u, n)
        for j in range(max_bc+1):   # careful: 0 ... max_bc
            mpz_clear(bc[j])
        sig_free(bc)
        mpz_fac_ui(t, k)
        mpz_tdiv_q(s, s, t)
        mpz_clear(t)
        mpz_clear(u)

def _stirling_number2(n, k):
    """
    Python wrapper of mpz_stirling_s2.

        sage: from sage.combinat.combinat_cython import _stirling_number2
        sage: _stirling_number2(3, 2)
        3

    This is wrapped again by stirling_number2 in combinat.py.
    """
    cdef Integer s = Integer.__new__(Integer)
    mpz_stirling_s2(s.value, n, k)
    return s

#####################################################################
## Set partition iterators

@cython.wraparound(False)
@cython.boundscheck(False)
cdef list from_word(list w, list base_set):
    cdef list sp = []
    cdef Py_ssize_t i
    cdef Py_ssize_t b
    for i in range(len(w)):
        b = <Py_ssize_t> (w[i])
        x = base_set[i]
        if len(sp) <= b:
            sp.append([x])
        else:
            sp[b].append(x)
    return sp

@cython.wraparound(False)
@cython.boundscheck(False)
def set_partition_iterator(base_set):
    """
    A fast iterator for the set partitions of the base set, which
    returns lists of lists instead of set partitions types.

    EXAMPLES::

        sage: from sage.combinat.combinat_cython import set_partition_iterator
        sage: list(set_partition_iterator([1,-1,x]))
        [[[1, -1, x]],
         [[1, -1], [x]],
         [[1, x], [-1]],
         [[1], [-1, x]],
         [[1], [-1], [x]]]
    """
    cdef list base = list(base_set)

    # Knuth, TAOCP 4A 7.2.1.5, Algorithm H
    cdef Py_ssize_t N = len(base)
    # H1: initialize
    cdef list a = [0] * N
    if N <= 1:
        yield from_word(a, base)
        return

    cdef list b = [1] * N
    cdef Py_ssize_t j
    cdef Py_ssize_t last = N - 1
    while True:
        # H2: visit
        yield from_word(a, base)
        if a[last] == b[last]:
            # H4: find j
            j = N - 2
            while a[j] == b[j]:
                j -= 1
            # H5: increase a_j
            if j == 0:
                break
            a[j] += 1
            # H6: zero out a_{j+1},...,a_{n-1}
            b[last] = b[j] + int(a[j] == b[j])
            j += 1
            while j < N - 1:
                a[j] = 0
                b[j] = b[last]
                j += 1
            a[last] = 0
        else:
            # H3: increase a_{n-1}
            a[last] += 1

@cython.wraparound(False)
@cython.boundscheck(False)
def _set_partition_block_gen(Py_ssize_t n, Py_ssize_t k, list a):
    r"""
    Recursively generate set partitions of ``n`` with fixed block
    size ``k`` using Algorithm 4.23 from *Combinatorial Generation*
    by Ruskey. ``a`` is a list of size ``n``.

    EXAMPLES::

        sage: from sage.combinat.combinat_cython import _set_partition_block_gen
        sage: a = list(range(3))
        sage: for p in _set_partition_block_gen(3, 2, a):
        ....:     print(p)
        [0, 1, 0]
        [0, 1, 1]
        [0, 0, 1]
    """
    cdef Py_ssize_t i
    if n == k:
        yield a
        return

    for i in range(k):
        a[n-1] = i
        for P in _set_partition_block_gen(n-1, k, a):
            yield P
        a[n-1] = n-1
    if k > 1:
        a[n-1] = k-1
        for P in _set_partition_block_gen(n-1, k-1, a):
            yield P
        a[n-1] = n-1

@cython.wraparound(False)
@cython.boundscheck(False)
def set_partition_iterator_blocks(base_set, Py_ssize_t k):
    """
    A fast iterator for the set partitions of the base set into the
    specified number of blocks, which returns lists of lists
    instead of set partitions types.

    EXAMPLES::

        sage: from sage.combinat.combinat_cython import set_partition_iterator_blocks
        sage: list(set_partition_iterator_blocks([1,-1,x], 2))
        [[[1, x], [-1]], [[1], [-1, x]], [[1, -1], [x]]]
    """
    cdef list base = list(base_set)
    cdef Py_ssize_t n = len(base)
    cdef list a = list(range(n))
    # TODO: implement _set_partition_block_gen as an iterative algorithm
    for P in _set_partition_block_gen(n, k, a):
        yield from_word(<list> P, base)

#####################################################################
## Perfect matchings iterator

cdef class PerfectMatchingsIterator(object):

    def __init__(self, Py_ssize_t n):
        cdef Py_ssize_t i
        self.n = n
        self.available = <bint*> check_allocarray(n, sizeof(bint))
        for i in range(n):
            self.available[i] = False
        self.current = [(i,i+1) for i in range(0,n,2)]
        self.at_end = (n <= 1) or (n % 2 != 0)

    def __dealloc__(self):
        sig_free(self.available)

    cpdef list next(self):
        if self.at_end:
            return None
        cdef list ret = list(self.current)
        cdef Py_ssize_t n = self.n
        cdef Py_ssize_t a, b
        a,b = self.current.pop()
        self.available[a] = True
        self.available[b] = True
        while True:

            a += 1 # Move a to the next position
            while a < n and not self.available[a]:
                a += 1
            if a == n:
                # Nothing more we can do at this level
                if not self.current:
                    break
                a,b = self.current.pop()
                self.available[a] = True
                self.available[b] = True
                continue
            b = a + 1
            while b < n and not self.available[b]:
                b += 1
            if b != n:
                self.available[a] = False
                self.available[b] = False
                self.current.append((a,b))
                a = 1
                if 2 * len(self.current) == n:
                    return ret

        self.at_end = True
        return ret

