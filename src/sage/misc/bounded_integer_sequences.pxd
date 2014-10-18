from sage.libs.gmp.types cimport *
from sage.misc.bitset cimport *

ctypedef struct biseq:       # bounded integer sequence
    bitset_t data            # Use GMP integers as bitarrays
    mp_bitcnt_t itembitsize  # Bitsize of one element of this sequence. Note:
                             # We do not store the exact bound for the items
                             # of this sequence, but store the bitlength that
                             # is sufficient to store one item.
    mp_limb_t mask_item      # "bla&mask_item" greps onethe item bla starts with
    mp_size_t length         # Number of items in this sequence.
                             # bitsize=length*itembitsize, but we do not want
                             # to repeat this multiplication and thus store
                             # the result.

ctypedef biseq biseq_t[1]

cdef bint allocate_biseq(biseq_t R, mp_size_t l, mp_size_t itemsize) except -1
   # Allocate memory (filled with zero) for a bounded integer sequence
   # of length l with items fitting in itemsize bits.

cdef void dealloc_biseq(biseq_t S)
   # Free the data stored in S.

cdef bint copy_biseq(biseq_t R, biseq_t S) except -1
   # Replace the content of R by a copy of S

cdef bint list_to_biseq(biseq_t R, list data, unsigned int bound) except -1
   # Convert a list to a bounded integer sequence

cdef list biseq_to_list(biseq_t S)
   # Convert a bounded integer sequence to a list

cdef bint concat_biseq(biseq_t R, biseq_t S1, biseq_t S2) except -1
   # Does not test whether the sequences have the same bound!

cdef bint first_bits_equal(mp_limb_t* b1, mp_limb_t* b2, mp_bitcnt_t d) except -1
   # Boilerplate function for comparison of the first bits of two gmp bitsets,
   # mimmicking mpz_congruent_2exp_p

cdef inline bint startswith_biseq(biseq_t S1, biseq_t S2) except -1
   # Is S1=S2+something? Does not check whether the sequences have the same
   # bound!

cdef int contains_biseq(biseq_t S1, biseq_t S2, mp_size_t start) except -2
   # Returns the position *in S1* of S2 as a subsequence of S1[start:], or -1
   # if S2 is not a subsequence. Does not check whether the sequences have the
   # same bound!

cdef int max_overlap_biseq(biseq_t S1, biseq_t S2) except 0
   # Returns the minimal *positive* integer i such that S2 starts with S1[i:].
   # This function will *not* test whether S2 starts with S1!

cdef int index_biseq(biseq_t S, mp_limb_t item, mp_size_t start) except -2
   # Returns the position *in S* of the item in S[start:], or -1 if S[start:]
   # does not contain the item.

cdef int getitem_biseq(biseq_t S, mp_size_t index) except -1
   # Returns S[index], without checking margins

cdef bint slice_biseq(biseq_t R, biseq_t S, mp_size_t start, mp_size_t stop, int step) except -1
   # Fills R with S[start:stop:step]

cdef class BoundedIntegerSequence:
    cdef biseq_t data
    cdef str str(self)
    cpdef bint startswith(self, BoundedIntegerSequence other)
    cpdef list list(self)
    cpdef BoundedIntegerSequence maximal_overlap(self, BoundedIntegerSequence other)

cpdef BoundedIntegerSequence NewBISEQ(tuple data, mp_bitcnt_t itembitsize, mp_size_t length)
