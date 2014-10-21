from sage.libs.gmp.types cimport *
from sage.libs.gmp.types cimport mp_bits_per_limb
from sage.libs.gmp.mpn cimport mpn_rshift, mpn_lshift, mpn_copyi, mpn_ior_n, mpn_zero, mpn_copyd, mpn_cmp
from sage.misc.bitset cimport *

ctypedef struct biseq_s:       # bounded integer sequence
    # Use bitsets to store the data, which in turn is based on GMP integers
    bitset_t data
    # Bitsize of one element of this sequence. Note: We do not store the exact
    # bound for the items of this sequence, but store the bitlength that is
    # sufficient to store one item.
    mp_bitcnt_t itembitsize
    # If `L` is a limb, then `L&mask_item` extracts the first `itembitsize`
    # bits of it, i.e., it extracts one item.
    mp_limb_t mask_item
    # Number of items in this sequence.
    mp_size_t length

ctypedef biseq_s biseq_t[1]

ctypedef mp_limb_t biseq_item_t

cdef bint biseq_init(biseq_t R, mp_size_t l, mp_size_t itemsize) except -1
   # Allocate memory (filled with zero) for a bounded integer sequence
   # of length l with items fitting in itemsize bits.

cdef void biseq_dealloc(biseq_t S)
   # Free the data stored in S.

cdef bint biseq_copy(biseq_t R, biseq_t S) except -1
   # Replace the content of R by a copy of S

cdef bint list_to_biseq(biseq_t R, list data, mp_limb_t bound) except -1
   # Convert a list to a bounded integer sequence

cdef list biseq_to_list(biseq_t S)
   # Convert a bounded integer sequence to a list

cdef bint biseq_concat(biseq_t R, biseq_t S1, biseq_t S2) except -1
   # Does not test whether the sequences have the same bound!

cdef bint biseq_first_bits_equal(mp_limb_t* b1, mp_limb_t* b2, mp_bitcnt_t d) except -1
   # Boilerplate function for comparison of the first bits of two gmp bitsets,
   # mimmicking mpz_congruent_2exp_p

cdef inline bint biseq_startswith(biseq_t S1, biseq_t S2) except -1
   # Is S1=S2+something? Does not check whether the sequences have the same
   # bound!

cdef int biseq_contains(biseq_t S1, biseq_t S2, mp_size_t start) except -2
   # Returns the position *in S1* of S2 as a subsequence of S1[start:], or -1
   # if S2 is not a subsequence. Does not check whether the sequences have the
   # same bound!

cdef int biseq_max_overlap(biseq_t S1, biseq_t S2) except 0
   # Returns the minimal *positive* integer i such that S2 starts with S1[i:].
   # This function will *not* test whether S2 starts with S1!

cdef int biseq_index(biseq_t S, mp_limb_t item, mp_size_t start) except -2
   # Returns the position *in S* of the item in S[start:], or -1 if S[start:]
   # does not contain the item.

cdef biseq_item_t biseq_getitem(biseq_t S, mp_size_t index) except -1
   # Returns S[index], without checking margins

cdef bint biseq_slice(biseq_t R, biseq_t S, mp_size_t start, mp_size_t stop, int step) except -1
   # Fills R with S[start:stop:step]

cdef class BoundedIntegerSequence:
    cdef biseq_t data
    cdef str str(self)
    cpdef bint startswith(self, BoundedIntegerSequence other)
    cpdef list list(self)
    cpdef BoundedIntegerSequence maximal_overlap(self, BoundedIntegerSequence other)

cpdef BoundedIntegerSequence NewBISEQ(tuple data, mp_bitcnt_t itembitsize, mp_size_t length)
