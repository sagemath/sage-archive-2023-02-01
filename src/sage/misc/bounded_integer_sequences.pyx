r"""
Sequences of bounded integers

AUTHORS:

- Simon King (2014): initial version

This module provides :class:`BoundedIntegerSequence`, which implements
sequences of bounded integers and is for many (but not all) operations faster
than representing the same sequence as a Python :class:`tuple`.

It also provides some boilerplate functions that can be cimported in Cython
modules::

    cdef biseq_t* allocate_biseq(size_t l, unsigned long int itemsize) except NULL
       # Allocate memory (filled with zero) for a bounded integer sequence
       # of length l with items fitting in itemsize bits.

    cdef biseq_t* list_to_biseq(biseq_t S, list data) except NULL
       # Assumes that S is allocated

    cdef list biseq_to_list(biseq_t S)
       # Convert a bounded integer sequence to a list

    cdef str biseq_to_str(biseq_t S)
       # String representation of bounded integer sequence as comma
       # separated list

    cdef biseq_t* concat_biseq(biseq_t S1, biseq_t S2) except NULL
       # Does not test whether the sequences have the same bound!

    cdef inline bint startswith_biseq(biseq_t S1, biseq_t S2)
       # Is S1=S2+something? Does not check whether the sequences have the same
       # bound!

    cdef int contains_biseq(biseq_t S1, biseq_t S2, size_t start)
       # Returns the position *in S1* of S2 as a subsequence of S1[start:], or -1
       # if S2 is not a subsequence. Does not check whether the sequences have the
       # same bound!

    cdef int index_biseq(biseq_t S, int item, size_t start)
       # Returns the position *in S* of the item in S[start:], or -1 if S[start:]
       # does not contain the item.

    cdef int getitem_biseq(biseq_t S, unsigned long int index) except -1
       # Returns S[index], without checking margins

    cdef biseq_t* slice_biseq(biseq_t S, int start, int stop, int step) except NULL
       # Returns the biseq S[start:stop:step]

"""
#*****************************************************************************
#       Copyright (C) 2014 Simon King <simon.king@uni-jena.de>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#  as published by the Free Software Foundation; either version 2 of
#  the License, or (at your option) any later version.
#                  http://www.gnu.org/licenses/
#*****************************************************************************

include "sage/ext/cdefs.pxi"
include "sage/ext/interrupt.pxi"
include "sage/ext/stdsage.pxi"
include "sage/ext/python.pxi"

cdef extern from "mpz_pylong.h":
    cdef long mpz_pythonhash(mpz_t src)

cdef extern from "gmp.h":
    cdef int mp_bits_per_limb
    mp_limb_t mpn_rshift(mp_ptr res, mp_ptr src, mp_size_t n, unsigned int count)
    mp_limb_t mpn_lshift(mp_ptr res, mp_ptr src, mp_size_t n, unsigned int count)
    void mpn_copyi(mp_ptr res, mp_ptr src, mp_size_t n)
    void mpn_ior_n (mp_limb_t *rp, mp_limb_t *s1p, mp_limb_t *s2p, mp_size_t n)
    void mpn_zero (mp_limb_t *rp, mp_size_t n)
    void mpn_copyd (mp_limb_t *rp, mp_limb_t *s1p, mp_size_t n)

cdef extern from "Python.h":
    bint PySlice_Check(PyObject* ob)

cdef size_t times_size_of_limb = [1,2,4,8,16,32,64].index(sizeof(mp_limb_t))
cdef size_t times_mp_bits_per_limb = [1,2,4,8,16,32,64,128,256].index(mp_bits_per_limb)
cdef size_t mod_mp_bits_per_limb = ((<size_t>1)<<times_mp_bits_per_limb)-1

cdef tuple ZeroNone = (0,None)
cdef PyObject* zeroNone = <PyObject*>ZeroNone
cdef dict EmptyDict = {}
cdef PyObject* emptyDict = <PyObject*>EmptyDict

from cython.operator import dereference as deref, preincrement as preinc

###################
#  Boilerplate
#   cdef functions
###################

#
# (De)allocation, copying
#

cdef biseq_t* allocate_biseq(size_t l, unsigned long int itemsize) except NULL:
    """
    Allocate memory for a bounded integer sequence of length l with items
    fitting in itemsize bits. Returns a pointer to the bounded integer
    sequence, or NULL on error.
    """
    cdef biseq_t out
    out.bitsize = l*itemsize
    out.length = l
    out.itembitsize = itemsize
    out.mask_item = ((<unsigned int>1)<<itemsize)-1
    sig_on()
    mpz_init2(out.data, out.bitsize+mp_bits_per_limb)
    sig_off()
    return &out

#
# Conversion
#

cdef biseq_t* list_to_biseq(biseq_t S, list data) except NULL:
    """
    Fill the content of a list into a bounded integer sequence.

    The bounded integer sequence is supposed to have enough memory allocated
    and initialised to zero. This function returns a pointer to the originally
    given bounded integer sequence, or NULL on error.

    """
    cdef unsigned long int item
    cdef __mpz_struct *S_boil
    S_boil = <__mpz_struct*>S.data
    cdef mp_limb_t item_limb, *tmp_limb
    S_boil._mp_size = (((len(data)*S.itembitsize)-1)>>times_mp_bits_per_limb)+1
    tmp_limb = <mp_limb_t*>sage_malloc(2<<times_size_of_limb)
    if tmp_limb==NULL:
        raise MemoryError("Cannot even allocate two long integers")
    cdef int offset_mod, offset_div
    cdef int offset = 0
    if not data:
        return &S
    for item in data:
        item_limb = <mp_limb_t>(item&S.mask_item)
        offset_mod = offset&mod_mp_bits_per_limb
        offset_div = offset>>times_mp_bits_per_limb
        if offset_mod:
            tmp_limb[1] = mpn_lshift(tmp_limb, &item_limb, 1, offset_mod)
            if offset_div < S_boil._mp_size:
                S_boil._mp_d[offset_div] |= tmp_limb[0]
                S_boil._mp_d[offset_div+1] = tmp_limb[1]
            else:
                S_boil._mp_d[offset_div] |= tmp_limb[0]
        else:
            tmp_limb[0] = item_limb
            S_boil._mp_d[offset_div] = tmp_limb[0]
        offset += S.itembitsize
    sage_free(tmp_limb)
    return &S

cdef list biseq_to_list(biseq_t S):
    """
    Convert a bounded integer sequence to a list of integers.
    """
    cdef __mpz_struct seq
    seq = deref(<__mpz_struct*>S.data)
    cdef unsigned int index, limb_index, bit_index
    index = 0
    cdef mp_limb_t *tmp_limb
    tmp_limb = <mp_limb_t*>sage_malloc(2<<times_size_of_limb)
    if tmp_limb==NULL:
        raise MemoryError("Cannot even allocate two long integers")
    cdef list L = []
    cdef size_t n
    for n from S.length>=n>0:
        limb_index = index>>times_mp_bits_per_limb
        bit_index  = index&mod_mp_bits_per_limb
        if bit_index:
            if bit_index+S.itembitsize >= mp_bits_per_limb:
                mpn_rshift(tmp_limb, seq._mp_d+limb_index, 2, bit_index)
                L.append(tmp_limb[0]&S.mask_item)
            else:
                L.append(((seq._mp_d[limb_index])>>bit_index)&S.mask_item)
        else:
            L.append(seq._mp_d[limb_index]&S.mask_item)
        index += S.itembitsize
    sage_free(tmp_limb)
    return L

cdef str biseq_to_str(biseq_t S):
    """
    String representation of bounded integer sequence as comma
    separated list of integers
    """
    cdef __mpz_struct seq
    seq = deref(<__mpz_struct*>S.data)
    cdef unsigned int index, limb_index, bit_index
    index = 0
    cdef mp_limb_t *tmp_limb
    tmp_limb = <mp_limb_t*>sage_malloc(2<<times_size_of_limb)
    if tmp_limb==NULL:
        raise MemoryError("Cannot even allocate two long integers")
    cdef list L = []
    cdef size_t n
    for n from S.length>=n>0:
        limb_index = index>>times_mp_bits_per_limb
        bit_index  = index&mod_mp_bits_per_limb
        if bit_index:
            if bit_index+S.itembitsize >= mp_bits_per_limb:
                mpn_rshift(tmp_limb, seq._mp_d+limb_index, 2, bit_index)
                L.append(repr(<int>(tmp_limb[0]&S.mask_item)))
            else:
                L.append(repr(<int>(((seq._mp_d[limb_index])>>bit_index)&S.mask_item)))
        else:
            L.append(repr(<int>(seq._mp_d[limb_index]&S.mask_item)))
        index += S.itembitsize
    sage_free(tmp_limb)
    return ', '.join(L)

#
# Arithmetics
#

cdef biseq_t *concat_biseq(biseq_t S1, biseq_t S2) except NULL:
    """
    Concatenate two bounded integer sequences.

    ASSUMPTION:

    - The two sequences must have equivalent bounds, i.e., the items on the
      sequences must fit into the same number of bits. This condition is not
      tested.

    OUTPUT:

    - A pointer to the concatenated sequence, or NULL on error.

    """
    cdef biseq_t out
    out.bitsize = S1.bitsize+S2.bitsize
    out.length = S1.length+S2.length
    out.itembitsize = S1.itembitsize # do not test == S2.itembitsize
    out.mask_item = S1.mask_item
    sig_on()
    mpz_init2(out.data, out.bitsize+mp_bits_per_limb)
    sig_off()
    mpz_mul_2exp(out.data, S2.data, S1.bitsize)
    mpz_ior(out.data, out.data, S1.data)
    return &out

cdef inline bint startswith_biseq(biseq_t S1, biseq_t S2):
    """
    Tests if bounded integer sequence S1 starts with bounded integer sequence S2

    ASSUMPTION:

    - The two sequences must have equivalent bounds, i.e., the items on the
      sequences must fit into the same number of bits. This condition is not
      tested.

    """
    return mpz_congruent_2exp_p(S1.data, S2.data, S2.bitsize)

cdef int contains_biseq(biseq_t S1, biseq_t S2, size_t start):
    """
    Tests if bounded integer sequence S1[start:] contains a sub-sequence S2

    INPUT:

    - ``S1``, ``S2`` -- two bounded integer sequences
    - ``start`` -- integer, start index

    ASSUMPTION:

    - The two sequences must have equivalent bounds, i.e., the items on the
      sequences must fit into the same number of bits. This condition is not
      tested.

    """
    if S1.length<S2.length+start:
        return -1
    cdef mpz_t tmp
    sig_on()
    mpz_init_set(tmp, S1.data)
    sig_off()
    mpz_fdiv_q_2exp(tmp, tmp, start*S1.itembitsize)
    cdef size_t i
    for i from start<=i<=S1.length-S2.length:
        if mpz_congruent_2exp_p(tmp, S2.data, S2.bitsize):
            mpz_clear(tmp)
            return i
        mpz_fdiv_q_2exp(tmp, tmp, S1.itembitsize)
    mpz_clear(tmp)
    return -1

cdef int index_biseq(biseq_t S, int item, size_t start):
    """
    Returns the position in S of an item in S[start:], or -1 if S[start:] does
    not contain the item.

    """
    if start>=S.length:
        return -1
    cdef mpz_t tmp, mpz_item
    sig_on()
    mpz_init_set(tmp, S.data)
    sig_off()
    mpz_fdiv_q_2exp(tmp, tmp, start*S.itembitsize)
    mpz_init_set_ui(mpz_item, item)
    cdef size_t i
    for i from start<=i<S.length:
        if mpz_congruent_2exp_p(tmp, mpz_item, S.itembitsize):
            mpz_clear(tmp)
            mpz_clear(mpz_item)
            return i
        mpz_fdiv_q_2exp(tmp, tmp, S.itembitsize)
    mpz_clear(tmp)
    mpz_clear(mpz_item)
    return -1

cdef int getitem_biseq(biseq_t S, unsigned long int index) except -1:
    """
    Get item S[index], without checking margins.

    """
    cdef __mpz_struct seq
    seq = deref(<__mpz_struct*>S.data)
    cdef unsigned long int limb_index, bit_index, limb_size
    index *= S.itembitsize
    limb_index = index>>times_mp_bits_per_limb
    bit_index  = index&mod_mp_bits_per_limb
    cdef mp_limb_t* tmp_limb
    cdef int out
    if bit_index:
        limb_size = ((bit_index+S.itembitsize-1)>>times_mp_bits_per_limb)+1
        if limb_size!=1:
            tmp_limb = <mp_limb_t*>sage_malloc(limb_size<<times_size_of_limb)
            if tmp_limb==NULL:
                raise MemoryError
            mpn_rshift(tmp_limb, seq._mp_d+limb_index, limb_size, bit_index)
            out = <int>deref(tmp_limb)
            sage_free(tmp_limb)
        else:
            out = (<int>(seq._mp_d[limb_index]))>>bit_index
    else:
        out = <int>(seq._mp_d[limb_index])
    return out&S.mask_item

cdef biseq_t* slice_biseq(biseq_t S, int start, int stop, int step) except NULL:
    """
    Create the slice S[start:stop:step] as bounded integer sequence.

    Return:

    - A pointer to the resulting bounded integer sequence, or NULL on error.

    """
    cdef unsigned long int length, length1
    if step>0:
        if stop>start:
            length = ((stop-start-1)//step)+1
        else:
            length = 0
    else:
        if stop>=start:
            length = 0
        else:
            length = ((stop-start+1)//step)+1
    cdef biseq_t out
    out = deref(allocate_biseq(length, S.itembitsize))
    cdef mpz_t tmp
    if step==1:
        sig_on()
        mpz_init_set(tmp, S.data)
        sig_off()
        mpz_fdiv_q_2exp(tmp, tmp, start*S.itembitsize)
        mpz_fdiv_r_2exp(out.data, tmp, out.bitsize)
        mpz_clear(tmp)
        return &out
    cdef mpz_t item
    sig_on()
    mpz_init_set(tmp, S.data)
    mpz_init2(item, S.itembitsize)
    sig_off()
    cdef unsigned long int bitstep
    if step>0:
        mpz_fdiv_q_2exp(tmp, tmp, start*S.itembitsize)
        bitstep = out.itembitsize*step
        for i from 0<=i<length:
            mpz_fdiv_r_2exp(item, tmp, out.itembitsize)
            mpz_mul_2exp(item, item, i*out.itembitsize)
            mpz_ior(out.data, out.data, item)
            mpz_fdiv_q_2exp(tmp, tmp, bitstep)
    else:
        bitstep = out.itembitsize*(-step)
        length1 = length-1
        mpz_fdiv_q_2exp(tmp, tmp, (start+length1*step)*S.itembitsize)
        for i from 0<=i<length:
            mpz_fdiv_r_2exp(item, tmp, out.itembitsize)
            mpz_mul_2exp(item, item, (length1-i)*out.itembitsize)
            mpz_ior(out.data, out.data, item)
            mpz_fdiv_q_2exp(tmp, tmp, bitstep)
    mpz_clear(tmp)
    mpz_clear(item)
    return &out

###########################################
# A cdef class that wraps the above, and
# behaves like a tuple

from sage.rings.integer import Integer
cdef class BoundedIntegerSequence:
    """
    A sequence of non-negative uniformely bounded integers

    INPUT:

    - ``bound``, non-negative integer. When zero, a :class:`ValueError`
      will be raised. Otherwise, the given bound is replaced by the next
      power of two that is greater than the given bound.
    - ``data``, a list of integers. The given integers will be truncated
      to be less than the bound.

    EXAMPLES:

    We showcase the similarities and differences between bounded integer
    sequences and lists respectively tuples.

    To distinguish from tuples or lists, we use pointed brackets for the
    string representation of bounded integer sequences::

        sage: from sage.misc.bounded_integer_sequences import BoundedIntegerSequence
        sage: S = BoundedIntegerSequence(21, [2, 7, 20]); S
        <2, 7, 20>

    Each bounded integer sequence has a bound that is a power of two, such
    that all its item are less than this bound (they are in fact truncated)::

        sage: S.bound()
        32
        sage: BoundedIntegerSequence(16, [2, 7, 20])
        <2, 7, 4>

    Bounded integer sequences are iterable, and we see that we can recover the
    originally given list, modulo the bound::

        sage: L = [randint(0,31) for i in range(5000)]
        sage: S = BoundedIntegerSequence(32, L)
        sage: list(L) == L
        True
        sage: S16 = BoundedIntegerSequence(16, L)
        sage: list(S16) == [x%S16.bound() for x in L]
        True

    Getting items and slicing works in the same way as for lists. Note,
    however, that slicing is an operation that is relatively slow for bounded
    integer sequences.  ::

        sage: n = randint(0,4999)
        sage: S[n] == L[n]
        True
        sage: m = randint(0,1000)
        sage: n = randint(3000,4500)
        sage: s = randint(1, 7)
        sage: list(S[m:n:s]) == L[m:n:s]
        True
        sage: list(S[n:m:-s]) == L[n:m:-s]
        True

    The :meth:`index` method works different for bounded integer sequences and
    tuples or lists. If one asks for the index of an item, the behaviour is
    the same. But we can also ask for the index of a sub-sequence::

        sage: L.index(L[200]) == S.index(L[200])
        True
        sage: S.index(S[100:2000])    # random
        100

    Similarly, containment tests work for both items and sub-sequences::

        sage: S[200] in S
        True
        sage: S[200:400] in S
        True

    Note, however, that containment of items will test for congruence modulo
    the bound of the sequence. Thus, we have::

        sage: S[200]+S.bound() in S
        True
        sage: L[200]+S.bound() in L
        False

    Bounded integer sequences are immutable, and thus copies are
    identical. This is the same for tuples, but of course not for lists::

        sage: T = tuple(S)
        sage: copy(T) is T
        True
        sage: copy(S) is S
        True
        sage: copy(L) is L
        False

    Concatenation works in the same way for list, tuples and bounded integer
    sequences::

        sage: M = [randint(0,31) for i in range(5000)]
        sage: T = BoundedIntegerSequence(32, M)
        sage: list(S+T)==L+M
        True
        sage: list(T+S)==M+L
        True
        sage: (T+S == S+T) == (M+L == L+M)
        True

    However, comparison works different for lists and bounded integer
    sequences. Bounded integer sequences are first compared by bound, then by
    length, and eventually by *reverse* lexicographical ordering::

        sage: S = BoundedIntegerSequence(21, [4,1,6,2,7,20,9])
        sage: T = BoundedIntegerSequence(51, [4,1,6,2,7,20])
        sage: S < T   # compare by bound, not length
        True
        sage: T < S
        False
        sage: S.bound() < T.bound()
        True
        sage: len(S) > len(T)
        True

    ::

        sage: T = BoundedIntegerSequence(21, [0,0,0,0,0,0,0,0])
        sage: S < T    # compare by length, not lexicographically
        True
        sage: T < S
        False
        sage: list(T) < list(S)
        True
        sage: len(T) > len(S)
        True

    ::

        sage: T = BoundedIntegerSequence(21, [4,1,5,2,8,20,9])
        sage: T > S   # compare by reverse lexicographic ordering...
        True
        sage: S > T
        False
        sage: len(S) == len(T)
        True
        sage: list(S)> list(T) # direct lexicographic ordering is different
        True

    TESTS:

    We test against various corner cases::

        sage: BoundedIntegerSequence(16, [2, 7, -20])
        Traceback (most recent call last):
        ...
        OverflowError: can't convert negative value to unsigned long
        sage: BoundedIntegerSequence(1, [2, 7, 0])
        <0, 1, 0>
        sage: BoundedIntegerSequence(0, [2, 7, 0])
        Traceback (most recent call last):
        ...
        ValueError: Positive bound expected
        sage: BoundedIntegerSequence(2, [])
        <>
        sage: BoundedIntegerSequence(2, []) == BoundedIntegerSequence(4, []) # The bounds differ
        False
        sage: BoundedIntegerSequence(16, [2, 7, 20])[1:1]
        <>

    """
    def __cinit__(self, unsigned long int bound, list data):
        """
        Allocate memory for underlying data

        INPUT:

        - ``bound``, non-negative integer
        - ``data``, ignored

        .. WARNING::

            If ``bound=0`` then no allocation is done.  Hence, this should
            only be done internally, when calling :meth:`__new__` without :meth:`__init__`.

        EXAMPLES::

            sage: from sage.misc.bounded_integer_sequences import BoundedIntegerSequence
            sage: BoundedIntegerSequence(21, [4,1,6,2,7,20,9])  # indirect doctest
            <4, 1, 6, 2, 7, 20, 9>

        """
        # In __init__, we'll raise an error if the bound is 0.
        cdef mpz_t tmp
        if bound!=0:
            mpz_init_set_ui(tmp, bound-1)
            self.data = deref(allocate_biseq(len(data), mpz_sizeinbase(tmp, 2)))
            mpz_clear(tmp)

    def __dealloc__(self):
        """
        Free the memory from underlying data

        EXAMPLES::

            sage: from sage.misc.bounded_integer_sequences import BoundedIntegerSequence
            sage: S = BoundedIntegerSequence(21, [4,1,6,2,7,20,9])
            sage: del S     # indirect doctest

        """
        mpz_clear(self.data.data)

    def __init__(self, unsigned long int bound, list data):
        """
        INPUT:

        - ``bound``, non-negative integer. When zero, a :class:`ValueError`
          will be raised. Otherwise, the given bound is replaced by the next
          power of two that is greater than the given bound.
        - ``data``, a list of integers. The given integers will be truncated
          to be less than the bound.

        EXAMPLES::

            sage: from sage.misc.bounded_integer_sequences import BoundedIntegerSequence
            sage: L = [randint(0,26) for i in range(5000)]
            sage: S = BoundedIntegerSequence(57, L)   # indirect doctest
            sage: list(S) == L
            True

        The given data are truncated according to the bound::

            sage: S = BoundedIntegerSequence(11, [4,1,6,2,7,20,9]); S
            <4, 1, 6, 2, 7, 4, 9>
            sage: S.bound()
            16
            sage: S = BoundedIntegerSequence(11, L)
            sage: [x%S.bound() for x in L] == list(S)
            True

        Non-positive bounds result in errors::

            sage: BoundedIntegerSequence(-1, L)
            Traceback (most recent call last):
            ...
            OverflowError: can't convert negative value to unsigned long
            sage: BoundedIntegerSequence(0, L)
            Traceback (most recent call last):
            ...
            ValueError: Positive bound expected

        """
        if bound==0:
            raise ValueError("Positive bound expected")
        self.data = deref(list_to_biseq(self.data, data))

    def __copy__(self):
        """
        :class:`BoundedIntegerSequence` is immutable, copying returns ``self``.

        EXAMPLES::

            sage: from sage.misc.bounded_integer_sequences import BoundedIntegerSequence
            sage: S = BoundedIntegerSequence(11, [4,1,6,2,7,20,9])
            sage: copy(S) is S
            True

        """
        return self

    def __reduce__(self):
        """
        Pickling of :class:`BoundedIntegerSequence`


        EXAMPLES::

            sage: from sage.misc.bounded_integer_sequences import BoundedIntegerSequence
            sage: L = [randint(0,26) for i in range(5000)]
            sage: S = BoundedIntegerSequence(32, L)
            sage: loads(dumps(S)) == S    # indirect doctest
            True

        """
        cdef size_t n
        cdef char *s
        n = mpz_sizeinbase(self.data.data, 32) + 2
        s = <char *>PyMem_Malloc(n)
        if s == NULL:
            raise MemoryError, "Unable to allocate enough memory for the string defining a bounded integer sequence."
        sig_on()
        mpz_get_str(s, 32, self.data.data)
        sig_off()
        data_str = <object> PyString_FromString(s)
        PyMem_Free(s)
        return NewBISEQ, (data_str, self.data.bitsize, self.data.itembitsize, self.data.length)

    def __len__(self):
        """
        EXAMPLES::

            sage: from sage.misc.bounded_integer_sequences import BoundedIntegerSequence
            sage: L = [randint(0,26) for i in range(5000)]
            sage: S = BoundedIntegerSequence(57, L)   # indirect doctest
            sage: len(S) == len(L)
            True

        """
        return self.data.length

    def __repr__(self):
        """
        String representation.

        To distinguish it from Python tuples or lists, we use pointed brackets
        as delimiters.

        EXAMPLES::

            sage: from sage.misc.bounded_integer_sequences import BoundedIntegerSequence
            sage: BoundedIntegerSequence(21, [4,1,6,2,7,20,9])   # indirect doctest
            <4, 1, 6, 2, 7, 20, 9>

        """
        return '<'+self.str()+'>'

    cdef str str(self):
        """
        A cdef helper function, returns the string representation without the brackets.

        Used in :meth:`__repr__`.

        EXAMPLES::

            sage: from sage.misc.bounded_integer_sequences import BoundedIntegerSequence
            sage: BoundedIntegerSequence(21, [4,1,6,2,7,20,9])   # indirect doctest
            <4, 1, 6, 2, 7, 20, 9>

        """
        return biseq_to_str(self.data)

    def bound(self):
        """
        The bound of this bounded integer sequence

        All items of this sequence are non-negative integers less than the
        returned bound. The bound is a power of two.

        EXAMPLES::

            sage: from sage.misc.bounded_integer_sequences import BoundedIntegerSequence
            sage: S = BoundedIntegerSequence(21, [4,1,6,2,7,20,9])
            sage: T = BoundedIntegerSequence(51, [4,1,6,2,7,20,9])
            sage: S.bound()
            32
            sage: T.bound()
            64

        """
        cdef long b = 1
        return (b<<self.data.itembitsize)

    def __iter__(self):
        """
        EXAMPLES::

            sage: from sage.misc.bounded_integer_sequences import BoundedIntegerSequence
            sage: L = [randint(0,26) for i in range(5000)]
            sage: S = BoundedIntegerSequence(27, L)
            sage: list(S) == L   # indirect doctest
            True

        """
        cdef __mpz_struct seq
        seq = deref(<__mpz_struct*>self.data.data)
        cdef unsigned int index, limb_index, bit_index
        index = 0
        cdef mp_limb_t *tmp_limb
        tmp_limb = <mp_limb_t*>sage_malloc(2<<times_size_of_limb)
        if tmp_limb==NULL:
            raise MemoryError("Cannot even allocate two long integers")
        cdef size_t n
        for n from self.data.length>=n>0:
            limb_index = index>>times_mp_bits_per_limb
            bit_index  = index&mod_mp_bits_per_limb
            if bit_index:
                if bit_index+self.data.itembitsize >= mp_bits_per_limb:
                    mpn_rshift(tmp_limb, seq._mp_d+limb_index, 2, bit_index)
                    yield <int>(tmp_limb[0]&self.data.mask_item)
                else:
                    yield <int>(((seq._mp_d[limb_index])>>bit_index)&self.data.mask_item)
            else:
                yield <int>(seq._mp_d[limb_index]&self.data.mask_item)
            index += self.data.itembitsize
        sage_free(tmp_limb)

    def __getitem__(self, index):
        """
        Get single items or slices.

        EXAMPLES::

            sage: from sage.misc.bounded_integer_sequences import BoundedIntegerSequence
            sage: S = BoundedIntegerSequence(21, [4,1,6,2,7,20,9])
            sage: S[2]
            6
            sage: S[1::2]
            <1, 2, 20>
            sage: S[-1::-2]
            <9, 7, 6, 4>

        TESTS::

            sage: L = [randint(0,26) for i in range(5000)]
            sage: S = BoundedIntegerSequence(27, L)
            sage: S[1234] == L[1234]
            True
            sage: list(S[100:2000:3]) == L[100:2000:3]
            True
            sage: list(S[3000:10:-7]) == L[3000:10:-7]
            True
            sage: S[:] == S
            True
            sage: S[:] is S
            True

        """
        cdef BoundedIntegerSequence out
        cdef int start,stop,step
        if PySlice_Check(<PyObject *>index):
            start,stop,step = index.indices(self.data.length)
            if start==0 and stop==self.data.length and step==1:
                return self
            out = BoundedIntegerSequence.__new__(BoundedIntegerSequence, 0, None)
            out.data = deref(slice_biseq(self.data, start, stop, step))
            return out
        cdef long Index
        try:
            Index = index
        except TypeError:
            raise TypeError("Sequence index must be integer or slice")
        if Index<0:
            Index = <long>(self.data.length)+Index
        if Index<0 or Index>=self.data.length:
            raise IndexError("Index out of range")
        return getitem_biseq(self.data, <size_t>Index)

    def __contains__(self, other):
        """
        Tells whether this bounded integer sequence contains an item or a sub-sequence

        EXAMPLES::

            sage: from sage.misc.bounded_integer_sequences import BoundedIntegerSequence
            sage: S = BoundedIntegerSequence(21, [4,1,6,2,7,20,9])
            sage: 6 in S
            True
            sage: BoundedIntegerSequence(21, [2, 7, 20]) in S
            True

        The bound of the sequences matters::

            sage: BoundedIntegerSequence(51, [2, 7, 20]) in S
            False

        Note that the items are compared up to congruence modulo the bound of
        the sequence. Thus we have::

            sage: 6+S.bound() in S
            True
            sage: S.index(6) == S.index(6+S.bound())
            True

        """
        if not isinstance(other, BoundedIntegerSequence):
            return index_biseq(self.data, other, 0)>=0
        cdef BoundedIntegerSequence right = other
        if self.data.itembitsize!=right.data.itembitsize:
            return False
        return contains_biseq(self.data, right.data, 0)>=0

    cpdef list list(self):
        """
        Converts this bounded integer sequence to a list

        NOTE:

        A conversion to a list is also possible by iterating over the
        sequence, which actually is faster.

        EXAMPLES::

            sage: from sage.misc.bounded_integer_sequences import BoundedIntegerSequence
            sage: L = [randint(0,26) for i in range(5000)]
            sage: S = BoundedIntegerSequence(32, L)
            sage: S.list() == list(S) == L
            True

        """
        return biseq_to_list(self.data)

    cpdef bint startswith(self, BoundedIntegerSequence other):
        """
        Tells whether ``self`` starts with a given bounded integer sequence

        EXAMPLES::

            sage: from sage.misc.bounded_integer_sequences import BoundedIntegerSequence
            sage: L = [randint(0,26) for i in range(5000)]
            sage: S = BoundedIntegerSequence(27, L)
            sage: L0 = L[:1000]
            sage: T = BoundedIntegerSequence(27, L0)
            sage: S.startswith(T)
            True
            sage: L0[-1] += 1
            sage: T = BoundedIntegerSequence(27, L0)
            sage: S.startswith(T)
            False
            sage: L0[-1] -= 1
            sage: L0[0] += 1
            sage: T = BoundedIntegerSequence(27, L0)
            sage: S.startswith(T)
            False
            sage: L0[0] -= 1

        The bounds of the sequences must be compatible, or :meth:`startswith`
        returns ``False``::

            sage: T = BoundedIntegerSequence(51, L0)
            sage: S.startswith(T)
            False

        """
        if self.data.itembitsize!=other.data.itembitsize:
            return False
        return startswith_biseq(self.data, other.data)

    def index(self, other):
        """
        The index of a given item or sub-sequence of ``self``

        EXAMPLES::

            sage: from sage.misc.bounded_integer_sequences import BoundedIntegerSequence
            sage: S = BoundedIntegerSequence(21, [4,1,6,2,6,20,9])
            sage: S.index(6)
            2
            sage: S.index(5)
            Traceback (most recent call last):
            ...
            ValueError: BoundedIntegerSequence.index(x): x(=5) not in sequence
            sage: S.index(-3)
            Traceback (most recent call last):
            ...
            ValueError: BoundedIntegerSequence.index(x): x(=-3) not in sequence
            sage: S.index(BoundedIntegerSequence(21, [6, 2, 6]))
            2
            sage: S.index(BoundedIntegerSequence(21, [6, 2, 7]))
            Traceback (most recent call last):
            ...
            ValueError: Not a sub-sequence

        The bound of (sub-)sequences matters::

            sage: S.index(BoundedIntegerSequence(51, [6, 2, 6]))
            Traceback (most recent call last):
            ...
            ValueError: Not a sub-sequence

        Note that items are compared up to congruence modulo the bound of the
        sequence::

            sage: S.index(6) == S.index(6+S.bound())
            True

        """
        cdef int out
        if not isinstance(other, BoundedIntegerSequence):
            if other<0:
                raise ValueError("BoundedIntegerSequence.index(x): x(={}) not in sequence".format(other))
            try:
                out = index_biseq(self.data, <unsigned int>other, 0)
            except TypeError:
                raise ValueError("BoundedIntegerSequence.index(x): x(={}) not in sequence".format(other))
            if out>=0:
                return out
            raise ValueError("BoundedIntegerSequence.index(x): x(={}) not in sequence".format(other))
        cdef BoundedIntegerSequence right = other
        if self.data.itembitsize!=right.data.itembitsize:
            raise ValueError("Not a sub-sequence")
        out = contains_biseq(self.data, right.data, 0)
        if out>=0:
            return out
        raise ValueError("Not a sub-sequence")

    def __add__(self, other):
        """
        Concatenation of bounded integer sequences.

        NOTE:

        There is no coercion happening, as bounded integer sequences are not
        considered to be elements of an object.

        EXAMPLES::

            sage: from sage.misc.bounded_integer_sequences import BoundedIntegerSequence
            sage: S = BoundedIntegerSequence(21, [4,1,6,2,7,20,9])
            sage: T = BoundedIntegerSequence(21, [4,1,6,2,8,15])
            sage: S+T
            <4, 1, 6, 2, 7, 20, 9, 4, 1, 6, 2, 8, 15>
            sage: T+S
            <4, 1, 6, 2, 8, 15, 4, 1, 6, 2, 7, 20, 9>
            sage: S in S+T
            True
            sage: T in S+T
            True
            sage: T+list(S)
            Traceback (most recent call last):
            ...
            TypeError:  Cannot convert list to sage.misc.bounded_integer_sequences.BoundedIntegerSequence
            sage: T+None
            Traceback (most recent call last):
            ...
            TypeError: Cannot concatenate bounded integer sequence and None

        """
        cdef BoundedIntegerSequence myself, right, out
        if other is None or self is None:
            raise TypeError('Cannot concatenate bounded integer sequence and None')
        myself = self  # may result in a type error
        right = other  #  --"--
        if right.data.itembitsize!=myself.data.itembitsize:
            raise ValueError("can only concatenate bounded integer sequences of compatible bounds")
        out = BoundedIntegerSequence.__new__(BoundedIntegerSequence, 0, None)
        out.data = deref(concat_biseq(myself.data, right.data))
        return out

    def __cmp__(self, other):
        """
        Comparison of bounded integer sequences

        We compare, in this order:

        - The bound of ``self`` and ``other``

        - The length of ``self`` and ``other``

        - Reverse lexicographical ordering, i.e., the sequences' items
          are compared starting with the last item.

        EXAMPLES:

        Comparison by bound::

            sage: from sage.misc.bounded_integer_sequences import BoundedIntegerSequence
            sage: S = BoundedIntegerSequence(21, [4,1,6,2,7,20,9])
            sage: T = BoundedIntegerSequence(51, [4,1,6,2,7,20,9])
            sage: S < T
            True
            sage: T < S
            False
            sage: list(T) == list(S)
            True

        Comparison by length::

            sage: T = BoundedIntegerSequence(21, [0,0,0,0,0,0,0,0])
            sage: S < T
            True
            sage: T < S
            False
            sage: list(T) < list(S)
            True
            sage: len(T) > len(S)
            True

        Comparison by *reverse* lexicographical ordering::

            sage: T = BoundedIntegerSequence(21, [4,1,5,2,8,20,9])
            sage: T > S
            True
            sage: S > T
            False
            sage: list(S)> list(T)
            True

        """
        cdef BoundedIntegerSequence right
        if other is None:
            return 1
        try:
            right = other
        except TypeError:
            return -1
        cdef int c = cmp(self.data.itembitsize, right.data.itembitsize)
        if c:
            return c
        c = cmp(self.data.length, right.data.length)
        if c:
            return c
        return mpz_cmp(self.data.data, right.data.data)

    def __hash__(self):
        """
        The hash takes into account the content and the bound of the sequence.

        EXAMPLES::

            sage: from sage.misc.bounded_integer_sequences import BoundedIntegerSequence
            sage: S = BoundedIntegerSequence(21, [4,1,6,2,7,20,9])
            sage: T = BoundedIntegerSequence(51, [4,1,6,2,7,20,9])
            sage: S == T
            False
            sage: list(S) == list(T)
            True
            sage: S.bound() == T.bound()
            False
            sage: hash(S) == hash(T)
            False
            sage: T = BoundedIntegerSequence(31, [4,1,6,2,7,20,9])
            sage: T.bound() == S.bound()
            True
            sage: hash(S) == hash(T)
            True

        """
        return mpz_pythonhash(self.data.data)

cpdef BoundedIntegerSequence NewBISEQ(data, unsigned long int bitsize, unsigned int itembitsize, size_t length):
    """
    Helper function for unpickling of :class:`BoundedIntegerSequence`.

    EXAMPLES::

        sage: from sage.misc.bounded_integer_sequences import BoundedIntegerSequence
        sage: L = [randint(0,26) for i in range(5000)]
        sage: S = BoundedIntegerSequence(32, L)
        sage: loads(dumps(S)) == S    # indirect doctest
        True

    """
    cdef BoundedIntegerSequence out = BoundedIntegerSequence.__new__(BoundedIntegerSequence, 0, None)
    mpz_init2(out.data.data, bitsize+mp_bits_per_limb)
    out.data.bitsize = bitsize
    out.data.itembitsize = itembitsize
    out.data.mask_item = ((<unsigned int>1)<<itembitsize)-1
    out.data.length = length
    mpz_set_str(out.data.data, data, 32)
    return out
