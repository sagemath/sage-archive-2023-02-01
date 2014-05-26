include "sage/libs/ntl/decl.pxi"

ctypedef struct biseq_t:  # bounded integer sequence type
    mpz_t data                # Use GMP integers as bitarrays
    unsigned long int bitsize # Bitsize of "data"
    unsigned int itembitsize  # Bitsize of one element of this sequence.
                              # Note: We do not store the exact bound for the
                              # items of this sequence, but store the
                              # bitlength that is sufficient to store one
                              # item.
    size_t length             # Number of items in this sequence.
                              # bitsize=length*itembitsize, but we do not want
                              # to repeat this multiplication and thus store
                              # the result.

cdef biseq_t* allocate_biseq(size_t l, unsigned long int itemsize) except NULL

#cdef inline void dealloc_biseq(biseq_t S)

cdef biseq_t* list2biseq(biseq_t S, list data) except NULL
   # Assumes that S is allocated

cdef list biseq2list(biseq_t S)

cdef str biseq2str(biseq_t S)

cdef biseq_t* concat_biseq(biseq_t S1, biseq_t S2) except NULL
   # Does not test whether the sequences have the same bound!

#cdef inline int cmp_biseq(biseq_t S1, biseq_t S2)

cdef inline bint startswith_biseq(biseq_t S1, biseq_t S2)
   # Is S1=S2+something? Does not check whether the sequences have the same
   # bound!

cdef int contains_biseq(biseq_t S1, biseq_t S2, size_t start)
   # Returns the position *in S1* of S2 as a subsequence of S1[start:], or -1
   # if S2 is not a subsequence. Does not check whether the sequences have the
   # same bound!

cdef int index_biseq(biseq_t S, int item, size_t start)
   # Returns the position *in S* of the item in S[start:], or -1 if S does not
   # contain the item.

cdef int getitem_biseq(biseq_t S, unsigned long int index)
   # Returns S[index], without checking margins

cdef biseq_t* slice_biseq(biseq_t S, int start, int stop, int step) except NULL
   # Returns the biseq S[start:stop:step]

cdef class BoundedIntegerSequence:
    cdef biseq_t data
    cdef str str(self)
    cpdef bint startswith(self, BoundedIntegerSequence other)
