include "sage/libs/ntl/decl.pxi"

ctypedef struct biseq:  # bounded integer sequence
    mpz_t data                # Use GMP integers as bitarrays
    unsigned long int bitsize # Bitsize of "data"
    unsigned int itembitsize  # Bitsize of one element of this sequence.
                              # Note: We do not store the exact bound for the
                              # items of this sequence, but store the
                              # bitlength that is sufficient to store one
                              # item.
    unsigned int mask_item    # "bla&mask_item" greps onethe item bla starts with
    size_t length             # Number of items in this sequence.
                              # bitsize=length*itembitsize, but we do not want
                              # to repeat this multiplication and thus store
                              # the result.

ctypedef biseq *biseq_t

cdef biseq_t allocate_biseq(size_t l, unsigned long int itemsize) except NULL
   # Allocate memory (filled with zero) for a bounded integer sequence
   # of length l with items fitting in itemsize bits.

cdef void dealloc_biseq(biseq_t S)

cdef biseq_t list_to_biseq(list data, unsigned int bound) except NULL
   # Convert a list to a bounded integer sequence

cdef list biseq_to_list(biseq_t S)
   # Convert a bounded integer sequence to a list

cdef biseq_t concat_biseq(biseq_t S1, biseq_t S2) except NULL
   # Does not test whether the sequences have the same bound!

cdef inline bint startswith_biseq(biseq_t S1, biseq_t S2)
   # Is S1=S2+something? Does not check whether the sequences have the same
   # bound!

cdef int contains_biseq(biseq_t S1, biseq_t S2, size_t start)
   # Returns the position *in S1* of S2 as a subsequence of S1[start:], or -1
   # if S2 is not a subsequence. Does not check whether the sequences have the
   # same bound!

cdef int max_overlap_biseq(biseq_t S1, biseq_t S2) except 0
   # Returns the minimal *positive* integer i such that S2 starts with S1[i:].
   # This function will *not* test whether S2 starts with S1!

cdef int index_biseq(biseq_t S, int item, size_t start) except -2
   # Returns the position *in S* of the item in S[start:], or -1 if S[start:]
   # does not contain the item.

cdef int getitem_biseq(biseq_t S, unsigned long int index) except -1
   # Returns S[index], without checking margins

cdef biseq_t slice_biseq(biseq_t S, int start, int stop, int step) except NULL
   # Returns the biseq S[start:stop:step]

cdef class BoundedIntegerSequence:
    cdef biseq_t data
    cdef str str(self)
    cpdef bint startswith(self, BoundedIntegerSequence other)
    cpdef list list(self)
    cpdef BoundedIntegerSequence maximal_overlap(self, BoundedIntegerSequence other)

cpdef BoundedIntegerSequence NewBISEQ(data, unsigned long int bitsize, unsigned int itembitsize, size_t length)
