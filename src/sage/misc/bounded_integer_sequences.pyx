r"""
Sequences of bounded integers

AUTHORS:

- Simon King (2014): initial version

This module provides :class:`BoundedIntegerSequence`, which implements
sequences of bounded integers and is for many (but not all) operations faster
than representing the same sequence as a Python :class:`tuple`.

The underlying data structure is similar to :class:`~sage.misc.bitset.Bitset`,
which means that certain operations are implemented by using fast shift
operations from GMP.  The following boilerplate functions can be cimported in
Cython modules:

- ``cdef bint allocate_biseq(biseq_t R, mp_size_t l, mp_size_t itemsize) except -1``

  Allocate memory (filled with zero) for a bounded integer sequence
  of length l with items fitting in itemsize bits.

- ``cdef void dealloc_biseq(biseq_t S)``

   Free the data stored in ``S``.

- ``cdef bint copy_biseq(biseq_t R, biseq_t S)``

  Replace the content of ``R`` by a copy of ``S``.

- ``cdef bint list_to_biseq(biseq_t R, list data, unsigned int bound) except -1``

  Convert a list to a bounded integer sequence.

- ``cdef list biseq_to_list(biseq_t S)``

  Convert a bounded integer sequence to a list.

- ``cdef bint concat_biseq(biseq_t R, biseq_t S1, biseq_t S2) except -1``

  Concatenate ``S1`` and ``S2`` and write the result to ``R``. Does not test
  whether the sequences have the same bound!

- ``cdef bint first_bits_equal(mp_limb_t* b1, mp_limb_t* b2, mp_bitcnt_t d) except -1``

  Boilerplate function for comparison of the first bits of two gmp bitsets,
  mimmicking ``mpz_congruent_2exp_p``.

- ``cdef inline bint startswith_biseq(biseq_t S1, biseq_t S2)``

  Is ``S1=S2+something``? Does not check whether the sequences have the same
  bound!

- ``cdef int contains_biseq(biseq_t S1, biseq_t S2, mp_size_t start) except -2``

  Returns the position *in ``S1``* of ``S2`` as a subsequence of
  ``S1[start:]``, or ``-1`` if ``S2`` is not a subsequence. Does not check
  whether the sequences have the same bound!

- ``cdef int max_overlap_biseq(biseq_t S1, biseq_t S2) except 0``

  Returns the minimal *positive* integer ``i`` such that ``S2`` starts with
  ``S1[i:]``.  This function will *not* test whether ``S2`` starts with ``S1``!

- ``cdef int index_biseq(biseq_t S, mp_limb_t item, mp_size_t start) except -2``

  Returns the position *in S* of the item in ``S[start:]``, or ``-1`` if
  ``S[start:]`` does not contain the item.

- ``cdef int getitem_biseq(biseq_t S, mp_size_t index) except -1``

  Returns ``S[index]``, without checking margins.

- ``cdef bint slice_biseq(biseq_t R, biseq_t S, mp_size_t start, mp_size_t stop, int step) except -1``

  Fills ``R`` with ``S[start:stop:step]``

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
include "sage/libs/ntl/decl.pxi"
include 'sage/misc/bitset.pxi'

#cdef extern from "mpz_pylong.h":
#    cdef long mpz_pythonhash(mpz_t src)

cdef extern from "gmp.h":
    cdef int mp_bits_per_limb
    mp_limb_t mpn_rshift(mp_ptr res, mp_ptr src, mp_size_t n, unsigned int count)
    mp_limb_t mpn_lshift(mp_ptr res, mp_ptr src, mp_size_t n, unsigned int count)
    void mpn_copyi(mp_ptr res, mp_ptr src, mp_size_t n)
    void mpn_ior_n (mp_limb_t *rp, mp_limb_t *s1p, mp_limb_t *s2p, mp_size_t n)
    void mpn_zero (mp_limb_t *rp, mp_size_t n)
    void mpn_copyd (mp_limb_t *rp, mp_limb_t *s1p, mp_size_t n)
    int mpn_cmp (mp_limb_t *s1p, mp_limb_t *s2p, mp_size_t n)

cdef extern from "Python.h":
    bint PySlice_Check(PyObject* ob)

cdef mp_size_t times_size_of_limb = [1,2,4,8,16,32,64].index(sizeof(mp_limb_t))
cdef mp_size_t times_mp_bits_per_limb = [1,2,4,8,16,32,64,128,256].index(mp_bits_per_limb)
cdef mp_size_t mod_mp_bits_per_limb = ((<mp_size_t>1)<<times_mp_bits_per_limb)-1

cdef tuple ZeroNone = (0,None)
cdef PyObject* zeroNone = <PyObject*>ZeroNone
cdef dict EmptyDict = {}
cdef PyObject* emptyDict = <PyObject*>EmptyDict

from cython.operator import dereference as deref, preincrement as preinc, predecrement as predec, postincrement as postinc

###################
#  Boilerplate
#   cdef functions
###################

#
# (De)allocation, copying
#

cdef bint allocate_biseq(biseq_t R, mp_size_t l, mp_size_t itemsize) except -1:
    """
    Allocate memory for a bounded integer sequence of length l with items
    fitting in itemsize bits.
    """
    R.length = l
    R.itembitsize = itemsize
    R.mask_item = ((<mp_limb_t>1)<<itemsize)-1
    if l:
        bitset_init(R.data, l*itemsize)

cdef void dealloc_biseq(biseq_t S):
    if S.length:
        bitset_free(S.data)

cdef bint copy_biseq(biseq_t R, biseq_t S) except -1:
    """
    Create a copy of S and store it in R, overriding R's previous content.
    """
    if S.length:
        if R.length:
            bitset_realloc(R.data, S.data.size)
        else:
            bitset_init(R.data, S.data.size)
        mpn_copyi(R.data.bits, S.data.bits, S.data.limbs)
    elif R.length:
        bitset_free(R.data)
    R.length = S.length
    R.itembitsize = S.itembitsize
    R.mask_item = S.mask_item

#
# Conversion
#

cdef bint list_to_biseq(biseq_t R, list data, unsigned int bound) except -1:
    """
    Convert a list into a bounded integer sequence.

    """
    cdef mpz_t tmp
    cdef mp_size_t ln
    if bound!=0:
        mpz_init_set_ui(tmp, bound-1)
        ln = mpz_sizeinbase(tmp, 2)
        if ln>mp_bits_per_limb:
            raise ValueError("The integer bound {} does not fit into {}".format(bound, mp_bits_per_limb))
        allocate_biseq(R, len(data), ln)
        mpz_clear(tmp)
    else:
        raise ValueError("The bound for the items of bounded integer sequences must be positive")
    cdef unsigned long int item
    cdef mp_limb_t item_limb
    cdef mp_limb_t tmp_limb1, tmp_limb2
    cdef int offset_mod, offset_div
    cdef int offset = 0
    if not data:
        return True
    for item in data:
        item_limb = <mp_limb_t>(item&R.mask_item)
        offset_mod = offset&mod_mp_bits_per_limb
        offset_div = offset>>times_mp_bits_per_limb
        if offset_mod:
            tmp_limb2 = mpn_lshift(&tmp_limb1, &item_limb, 1, offset_mod)
            if tmp_limb2:
                # This can only happen, if offset_div is small enough to not
                # write out of bounds, since initially we have allocated
                # enough memory.
                R.data.bits[offset_div+1] = tmp_limb2
            R.data.bits[offset_div] |= tmp_limb1
        else:
            R.data.bits[offset_div] = item_limb
        offset += R.itembitsize

cdef list biseq_to_list(biseq_t S):
    """
    Convert a bounded integer sequence to a list of integers.
    """
    cdef mp_size_t limb_index, max_limbs
    cdef mp_bitcnt_t bit_index, local_index
    cdef list L = []
    if S.length==0:
        return L
    max_limbs = S.data.limbs
    predec(max_limbs)
    cdef mp_size_t n
    cdef mp_limb_t tmp_limb[2]
    limb_index = 0
    bit_index = 0   # This is the index mod mp_bits_per_limb in S.data
    local_index = 0 # This is the index in tmp_limb
    tmp_limb[0] = S.data.bits[0]
    for n from S.length>=n>0:
        #print local_index, limb_index, bit_index
        if local_index+S.itembitsize > mp_bits_per_limb:
            #print "need more stuff"
            local_index = 0
            # tmp_limb[0] can not provide the next item, we need to get new stuff
            if bit_index and limb_index < max_limbs:
                # We move two limbs now, so that tmp_limb[0] will be completely full
                #print "shifting two limbs"
                mpn_rshift(<mp_limb_t*>tmp_limb, S.data.bits+limb_index, 2, bit_index)
            else:
                #print "one limb is enough"
                # We shift within one limb
                tmp_limb[0] = S.data.bits[limb_index]>>bit_index
        #else: tmp_limb[0] provides the next item
        local_index += S.itembitsize
        #print "next item from",Integer(tmp_limb[0]).bits(),"is",<int>(tmp_limb[0]&S.mask_item)
        L.append(<int>(tmp_limb[0]&S.mask_item))
        tmp_limb[0]>>=S.itembitsize
        bit_index += S.itembitsize
        if bit_index>=mp_bits_per_limb:
            bit_index -= mp_bits_per_limb
            preinc(limb_index)
    return L

#
# Arithmetics
#

cdef bint concat_biseq(biseq_t R, biseq_t S1, biseq_t S2) except -1:
    """
    Concatenate two bounded integer sequences ``S1`` and ``S2``.

    ASSUMPTION:

    - The two sequences must have equivalent bounds, i.e., the items on the
      sequences must fit into the same number of bits. This condition is not
      tested.

    OUTPUT:

    The result is written into ``R``, which must not be initialised

    """
    R.itembitsize = S1.itembitsize # do not test == S2.itembitsize
    R.mask_item = S1.mask_item
    if S1.length==0:
        if S2.length==0:
            R.length = 0
            return True
        else:
            R.length = S2.length
            bitset_init(R.data, S2.data.size)
            mpn_copyi(R.data.bits, S2.data.bits, S2.data.limbs)
    else:
        if S2.length==0:
            R.length = S1.length
            bitset_init(R.data, S1.data.size)
            mpn_copyi(R.data.bits, S1.data.bits, S1.data.limbs)
        else:
            R.length = S1.length+S2.length
            bitset_init(R.data, S1.data.size+S2.data.size)
            mpn_lshift(R.data.bits+(S1.data.size>>times_mp_bits_per_limb), S2.data.bits,
                       S2.data.limbs, S1.data.size&mod_mp_bits_per_limb)
            mpn_ior_n(R.data.bits, R.data.bits, S1.data.bits, S1.data.limbs)

cdef inline bint first_bits_equal(mp_limb_t* b1, mp_limb_t* b2, mp_bitcnt_t d) except -1:
    """
    Compare the first bits of two gmp bitsets

    INPUT:

    - b1,b2: Pointers to bitsets
    - d: The number of bits to be compared

    ASSUMPTION:

    d must not exceed the length of either bitset.

    """
    cdef mp_size_t i, dlimbs
    cdef unsigned long dbits
    dlimbs = d>>times_mp_bits_per_limb
    if mpn_cmp(b1,b2,dlimbs)!=0:
        return False
    dbits  = d&mod_mp_bits_per_limb
    if dbits==0:
        return True
    return (((b1[dlimbs])^(b2[dlimbs]))&(((<mp_limb_t>1)<<dbits)-1))==0


cdef inline bint startswith_biseq(biseq_t S1, biseq_t S2) except -1:
    """
    Tests if bounded integer sequence S1 starts with bounded integer sequence S2

    ASSUMPTION:

    - The two sequences must have equivalent bounds, i.e., the items on the
      sequences must fit into the same number of bits. This condition is not
      tested.

    """
    if S2.length == 0:
        return True
    if S2.length>S1.length:
        return False
    return first_bits_equal(S1.data.bits, S2.data.bits, S2.data.size)

cdef int contains_biseq(biseq_t S1, biseq_t S2, mp_size_t start) except -2:
    """
    Tests if the bounded integer sequence ``S1[start:]`` contains a sub-sequence ``S2``

    INPUT:

    - ``S1``, ``S2`` -- two bounded integer sequences
    - ``start`` -- integer, start index

    OUTPUT:

    Index ``i>=start`` such that ``S1[i:]`` starts with ``S2``, or ``-1`` if
    ``S1[start:]`` does not contain ``S2``.

    ASSUMPTION:

    - The two sequences must have equivalent bounds, i.e., the items on the
      sequences must fit into the same number of bits. This condition is not
      tested.

    """
    if S1.length<S2.length+start:
        return -1
    if S2.length==0:
        return start
    cdef bitset_t tmp
    bitset_init(tmp, S2.data.size+mp_bits_per_limb)
    cdef mp_size_t index = 0
    cdef mp_bitcnt_t n = 0
    cdef mp_size_t limb_index, bit_index
    cdef mp_bitcnt_t offset = S2.data.size&mod_mp_bits_per_limb
    for index from start <= index <= S1.length-S2.length:
        limb_index = n>>times_mp_bits_per_limb
        bit_index = n&mod_mp_bits_per_limb
        # We shift a part of S1 (with length S2) by bit_index, and compare
        # with S2.
        if bit_index:
            if offset+bit_index>mp_bits_per_limb:
                mpn_rshift(tmp.bits, S1.data.bits+limb_index, S2.data.limbs+1, bit_index)
            else:
                mpn_rshift(tmp.bits, S1.data.bits+limb_index, S2.data.limbs, bit_index)
        else:
            mpn_copyi(tmp.bits, S1.data.bits+limb_index, S2.data.limbs)
        if first_bits_equal(tmp.bits, S2.data.bits, S2.data.size):
            bitset_free(tmp)
            return index
        n += S1.itembitsize
    bitset_free(tmp)
    return -1

cdef int index_biseq(biseq_t S, mp_limb_t item, mp_size_t start) except -2:
    """
    Returns the position in S of an item in S[start:], or -1 if S[start:] does
    not contain the item.

    """
    cdef mp_size_t limb_index, max_limbs
    cdef mp_bitcnt_t bit_index, local_index
    if S.length==0:
        return -1
    max_limbs = S.data.limbs
    predec(max_limbs)
    cdef mp_size_t n
    cdef mp_limb_t tmp_limb[2]
    limb_index = 0
    bit_index = 0   # This is the index mod mp_bits_per_limb in S.data
    local_index = 0 # This is the index in tmp_limb
    tmp_limb[0] = S.data.bits[0]
    item &= S.mask_item
    for n from 0<=n<S.length:
        if local_index+S.itembitsize > mp_bits_per_limb:
            local_index = 0
            # tmp_limb[0] can not provide the next item, we need to get new stuff
            if bit_index and limb_index < max_limbs:
                # We move two limbs now, so that tmp_limb[0] will be completely full
                mpn_rshift(<mp_limb_t*>tmp_limb, S.data.bits+limb_index, 2, bit_index)
            else:
                # We shift within one limb
                tmp_limb[0] = S.data.bits[limb_index]>>bit_index
        #else: tmp_limb[0] provides the next item
        local_index += S.itembitsize
        if item == tmp_limb[0]&S.mask_item:
            return n
        tmp_limb[0]>>=S.itembitsize
        bit_index += S.itembitsize
        if bit_index>=mp_bits_per_limb:
            bit_index -= mp_bits_per_limb
            preinc(limb_index)
    return -1

cdef int getitem_biseq(biseq_t S, mp_size_t index) except -1:
    """
    Get item S[index], without checking margins.

    """
    cdef mp_size_t limb_index, bit_index
    index *= S.itembitsize
    limb_index = index>>times_mp_bits_per_limb
    bit_index  = index&mod_mp_bits_per_limb
    cdef mp_limb_t tmp_limb[2]
    cdef int out
    if bit_index:
        # limb_size is 1 or 2
        if bit_index+S.itembitsize>mp_bits_per_limb:
            mpn_rshift(tmp_limb, S.data.bits+limb_index, 2, bit_index)
            out = <int>(tmp_limb[0])
        else:
            out = <int>((S.data.bits[limb_index])>>bit_index)
    else:
        out = <int>(S.data.bits[limb_index])
    return out&S.mask_item

cdef bint slice_biseq(biseq_t R, biseq_t S, mp_size_t start, mp_size_t stop, int step) except -1:
    """
    Create the slice S[start:stop:step] as bounded integer sequence.

    """
    cdef mp_size_t length, total_shift, n
    if step>0:
        if stop>start:
            length = ((stop-start-1)//step)+1
        else:
            allocate_biseq(R, 0, S.itembitsize)
            return True
    else:
        if stop>=start:
            allocate_biseq(R, 0, S.itembitsize)
            return True
        else:
            length = ((stop-start+1)//step)+1
    allocate_biseq(R, length, S.itembitsize)
    cdef mp_size_t offset_mod, offset_div
    total_shift = start*S.itembitsize
    if step==1:
        # Slicing essentially boils down to a shift operation.
        offset_mod = total_shift&mod_mp_bits_per_limb
        offset_div = total_shift>>times_mp_bits_per_limb
        if offset_mod:
            #print "mpn_rshift",R.data.limbs,"by",offset_mod,"bits"
            mpn_rshift(R.data.bits, S.data.bits+offset_div, R.data.limbs, offset_mod)
            # Perhaps the last to-be-moved item partially belongs to the next
            # limb of S?
            #print "->",bitset_string(R.data)
            if (R.data.size&mod_mp_bits_per_limb)+offset_mod>mp_bits_per_limb:
                #print "adjust the final limb"
                R.data.bits[R.data.limbs-1] |= (S.data.bits[offset_div+R.data.limbs]<<(mp_bits_per_limb-offset_mod))
                #print "->", bitset_string(R.data)
        else:
            #print "just copying",R.data.limbs,"limbs"
            mpn_copyi(R.data.bits, S.data.bits+offset_div, R.data.limbs)
            #print "->", bitset_string(R.data)
        # Perhaps the last few bits of the last limb have to be set to zero?
        offset_mod = (length*S.itembitsize)&mod_mp_bits_per_limb
        if offset_mod:
            #print "focus on",offset_mod,"bits"
            R.data.bits[R.data.limbs-1] &= ((<mp_limb_t>1)<<offset_mod) - 1
            #print "->",bitset_string(R.data)
        return True
    # Now for the difficult case: 
    #
    # Convention for variable names: We move data from *_src to *_tgt
    # tmp_limb is used to temporarily store the data that we move/shift
    cdef mp_limb_t tmp_limb[2]
    cdef mp_limb_t tmp_limb2
    cdef mp_size_t offset_div_src, max_limb
    cdef mp_bitcnt_t offset_mod_src
    cdef mp_size_t offset_src = total_shift
    cdef mp_bitcnt_t bitstep = step*S.itembitsize
    cdef mp_bitcnt_t offset_mod_tgt
    cdef mp_size_t offset_div_tgt
    cdef mp_size_t offset_tgt = 0
    for n from length>=n>0:
        offset_div_src = offset_src>>times_mp_bits_per_limb
        offset_mod_src = offset_src&mod_mp_bits_per_limb
        offset_div_tgt = offset_tgt>>times_mp_bits_per_limb
        offset_mod_tgt = offset_tgt&mod_mp_bits_per_limb
        # put data from src to tmp_limb[0].
        if offset_mod_src:
            if (offset_mod_src+S.itembitsize > mp_bits_per_limb):
                mpn_rshift(tmp_limb, S.data.bits+offset_div_src, 2, offset_mod_src)
                tmp_limb[0] &= S.mask_item
            else:
                tmp_limb[0] = ((S.data.bits[offset_div_src])>>offset_mod_src) & S.mask_item
        else:
            tmp_limb[0] = S.data.bits[offset_div_src] & S.mask_item
        # put data from tmp_limb[0] to tgt
        if offset_mod_tgt:
            tmp_limb2 = mpn_lshift(tmp_limb, tmp_limb, 1, offset_mod_tgt)
            if tmp_limb2:
                R.data.bits[offset_div_tgt+1] = tmp_limb2
            R.data.bits[offset_div_tgt]  |= tmp_limb[0]
        else:
            R.data.bits[offset_div_tgt] = tmp_limb[0]
        offset_tgt += S.itembitsize
        offset_src += bitstep

cdef int max_overlap_biseq(biseq_t S1, biseq_t S2) except 0:
    """
    Returns the smallest **positive** integer ``i`` such that ``S2`` starts with ``S1[i:]``.

    Returns ``-1`` if there is no overlap. Note that ``i==0`` (``S2`` starts with ``S1``)
    will not be considered!

    INPUT:

    - ``S1``, ``S2`` -- two bounded integer sequences

    ASSUMPTION:

    The two sequences must have equivalent bounds, i.e., the items on the
    sequences must fit into the same number of bits. This condition is not
    tested.

    """
    # Idea: We store the maximal possible tail of S1 in a temporary bitset,
    # and then rightshift it (taking care that the number of to-be-shifted
    # limbs will decrease over time), comparing with the initial part of S2
    # after each step.
    if S1.length == 0:
        raise ValueError("First argument must be of positive length")
    cdef bitset_t tmp
    cdef mp_size_t limb_index
    cdef mp_bitcnt_t bit_index, n
    cdef mp_size_t index, i, start_index
    if S2.length>=S1.length:
        start_index = 1
    else:
        start_index = S1.length-S2.length
    if start_index == S1.length:
        return -1
    n = S1.itembitsize*start_index
    ## Initilise the temporary bitset
    # Allocate an additional limb, to cope with shifting accross limb borders
    bitset_init(tmp, (mp_bits_per_limb+S1.data.size)-n)
    cdef mp_size_t nlimbs = tmp.limbs # number of to-be-shifted limbs
    limb_index = n>>times_mp_bits_per_limb
    bit_index  = n&mod_mp_bits_per_limb
    if bit_index==0:
        mpn_copyi(tmp.bits, S1.data.bits+limb_index, predec(nlimbs))
    else:
        if ((S1.data.size-n)&mod_mp_bits_per_limb)+bit_index>mp_bits_per_limb:
            # Need to shift an additional limb
            mpn_rshift(tmp.bits, S1.data.bits+limb_index, nlimbs, bit_index)
        else:
            mpn_rshift(tmp.bits, S1.data.bits+limb_index, predec(nlimbs), bit_index)
    cdef mp_bitcnt_t remaining_size = S1.data.size - n
    n = 0  # This now counts how many bits of the highest limb of tmp we have shifted.
    for index from start_index<=index<S1.length:
        if first_bits_equal(tmp.bits, S2.data.bits, remaining_size):
            bitset_free(tmp)
            return index
        remaining_size -= S1.itembitsize
        if n>=mp_bits_per_limb:
            mpn_rshift(tmp.bits, tmp.bits, predec(nlimbs), S1.itembitsize)
            n = 0
        else:
            mpn_rshift(tmp.bits, tmp.bits, nlimbs, S1.itembitsize)
            n += S1.itembitsize
    bitset_free(tmp)
    return -1


###########################################
# A cdef class that wraps the above, and
# behaves like a tuple

from sage.rings.integer import Integer
cdef class BoundedIntegerSequence:
    """
    A sequence of non-negative uniformely bounded integers.

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
        sage: list(S) > list(T) # direct lexicographic ordering is different
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
    def __cinit__(self, *args, **kwds):
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
        self.data.length = 0

    def __dealloc__(self):
        """
        Free the memory from underlying data

        EXAMPLES::

            sage: from sage.misc.bounded_integer_sequences import BoundedIntegerSequence
            sage: S = BoundedIntegerSequence(21, [4,1,6,2,7,20,9])
            sage: del S     # indirect doctest

        """
        dealloc_biseq(self.data)

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
        list_to_biseq(self.data, data, bound)

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

        TESTS:

        The discussion at :trac:`15820` explains why the following is a good test::

            sage: X = BoundedIntegerSequence(21, [4,1,6,2,7,2,3])
            sage: S = BoundedIntegerSequence(21, [0,0,0,0,0,0,0])
            sage: loads(dumps(X+S))
            <4, 1, 6, 2, 7, 2, 3, 0, 0, 0, 0, 0, 0, 0>
            sage: loads(dumps(X+S)) == X+S
            True
            sage: T = BoundedIntegerSequence(21, [0,4,0,1,0,6,0,2,0,7,0,2,0,3,0,0,0,0,0,0,0,0,0,0,0,0,0,0])
            sage: T[1::2]
            <4, 1, 6, 2, 7, 2, 3, 0, 0, 0, 0, 0, 0, 0>
            sage: T[1::2] == X+S
            True
            sage: loads(dumps(X[1::2])) == X[1::2]
            True

        """
        return NewBISEQ, (bitset_pickle(self.data.data), self.data.itembitsize, self.data.length)

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

    def __nonzero__(self):
        """
        A bounded integer sequence is nonzero if and only if its length is nonzero.

        EXAMPLES::

            sage: from sage.misc.bounded_integer_sequences import BoundedIntegerSequence
            sage: S = BoundedIntegerSequence(13, [0,0,0])
            sage: bool(S)
            True
            sage: bool(S[1:1])
            False

        """
        return self.data.length!=0

    def __repr__(self):
        """
        String representation.

        To distinguish it from Python tuples or lists, we use pointed brackets
        as delimiters.

        EXAMPLES::

            sage: from sage.misc.bounded_integer_sequences import BoundedIntegerSequence
            sage: BoundedIntegerSequence(21, [4,1,6,2,7,20,9])   # indirect doctest
            <4, 1, 6, 2, 7, 20, 9>
            sage: BoundedIntegerSequence(21, [0,0]) + BoundedIntegerSequence(21, [0,0])
            <0, 0, 0, 0>

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
        return ', '.join([repr(n) for n in biseq_to_list(self.data)])

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
        return Integer(1)<<self.data.itembitsize

    def __iter__(self):
        """
        EXAMPLES::

            sage: from sage.misc.bounded_integer_sequences import BoundedIntegerSequence
            sage: L = [randint(0,26) for i in range(5000)]
            sage: S = BoundedIntegerSequence(27, L)
            sage: list(S) == L   # indirect doctest
            True

        TESTS:

        The discussion at :trac:`15820` explains why this is a good test::

            sage: S = BoundedIntegerSequence(21, [0,0,0,0,0,0,0])
            sage: X = BoundedIntegerSequence(21, [4,1,6,2,7,2,3])
            sage: list(X)
            [4, 1, 6, 2, 7, 2, 3]
            sage: list(X+S)
            [4, 1, 6, 2, 7, 2, 3, 0, 0, 0, 0, 0, 0, 0]
            sage: list(BoundedIntegerSequence(21, [0,0]) + BoundedIntegerSequence(21, [0,0]))
            [0, 0, 0, 0]

        """
        if self.data.length==0:
            return
        cdef mp_size_t limb_index, max_limbs
        cdef mp_bitcnt_t bit_index, local_index
        max_limbs = self.data.data.limbs
        predec(max_limbs)
        cdef mp_size_t n
        cdef mp_limb_t tmp_limb[2]
        limb_index = 0
        bit_index = 0   # This is the index mod mp_bits_per_limb in self.data.data
        local_index = 0 # This is the index in tmp_limb
        tmp_limb[0] = self.data.data.bits[0]
        for n from self.data.length>=n>0:
            if local_index+self.data.itembitsize > mp_bits_per_limb:
                local_index = 0
                # tmp_limb[0] can not provide the next item, we need to get new stuff
                if bit_index and limb_index < max_limbs:
                    # We move two limbs now, so that tmp_limb[0] will be completely full
                    mpn_rshift(<mp_limb_t*>tmp_limb, self.data.data.bits+limb_index, 2, bit_index)
                else:
                    # We shift within one limb
                    tmp_limb[0] = self.data.data.bits[limb_index]>>bit_index
            #else: tmp_limb[0] provides the next item
            local_index += self.data.itembitsize
            yield <int>(tmp_limb[0]&self.data.mask_item)
            tmp_limb[0]>>=self.data.itembitsize
            bit_index += self.data.itembitsize
            if bit_index>=mp_bits_per_limb:
                bit_index -= mp_bits_per_limb
                preinc(limb_index)


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

            sage: S = BoundedIntegerSequence(8, [4,1,6,2,7,2,5,5,2])
            sage: S[-1::-2]
            <2, 5, 7, 6, 4>
            sage: S[1::2]
            <1, 2, 2, 5>

        ::

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

        ::

            sage: S = BoundedIntegerSequence(21, [0,0,0,0,0,0,0])
            sage: X = BoundedIntegerSequence(21, [4,1,6,2,7,2,3])
            sage: (X+S)[6]
            3
            sage: (X+S)[10]
            0
            sage: (X+S)[12:]
            <0, 0>

        ::

            sage: S[2:2] == X[4:2]
            True

        ::

            sage: S = BoundedIntegerSequence(6, [3, 5, 3, 1, 5, 2, 2, 5, 3, 3, 4])
            sage: S[10]
            4

        ::

            sage: B = BoundedIntegerSequence(27, [8, 8, 26, 18, 18, 8, 22, 4, 17, 22, 22, 7, 12, 4, 1, 7, 21, 7, 10, 10])
            sage: B[8:]
            <17, 22, 22, 7, 12, 4, 1, 7, 21, 7, 10, 10>

        """
        cdef BoundedIntegerSequence out
        cdef int start,stop,step
        if PySlice_Check(<PyObject *>index):
            start,stop,step = index.indices(self.data.length)
            if start==0 and stop==self.data.length and step==1:
                return self
            out = BoundedIntegerSequence.__new__(BoundedIntegerSequence, 0, None)
            slice_biseq(out.data, self.data, start, stop, step)
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
        return getitem_biseq(self.data, <mp_size_t>Index)

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

        TESTS:

        The discussion at :trac:`15820` explains why the following are good tests::

            sage: X = BoundedIntegerSequence(21, [4,1,6,2,7,2,3])
            sage: S = BoundedIntegerSequence(21, [0,0,0,0,0,0,0])
            sage: loads(dumps(X+S))
            <4, 1, 6, 2, 7, 2, 3, 0, 0, 0, 0, 0, 0, 0>
            sage: loads(dumps(X+S)) == X+S
            True
            sage: T = BoundedIntegerSequence(21, [0,4,0,1,0,6,0,2,0,7,0,2,0,3,0,0,0,0,0,0,0,0,0,0,0,0,0,0])
            sage: T[3::2]==(X+S)[1:]
            True
            sage: T[3::2] in X+S
            True
            sage: T = BoundedIntegerSequence(21, [0,4,0,1,0,6,0,2,0,7,0,2,0,3,0,0,0,16,0,0,0,0,0,0,0,0,0,0])
            sage: T[3::2] in (X+S)
            False

        ::

            sage: S1 = BoundedIntegerSequence(3,[1,3])
            sage: S2 = BoundedIntegerSequence(3,[0])
            sage: S2 in S1
            False

        ::

            sage: B = BoundedIntegerSequence(27, [8, 8, 26, 18, 18, 8, 22, 4, 17, 22, 22, 7, 12, 4, 1, 7, 21, 7, 10, 10])
            sage: B.index(B[8:])
            8

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

        The discussion at :trac:`15820` explains why the following is a good test::

            sage: (BoundedIntegerSequence(21, [0,0]) + BoundedIntegerSequence(21, [0,0])).list()
            [0, 0, 0, 0]

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

        TESTS::

            sage: S = BoundedIntegerSequence(8, [2, 2, 2, 1, 2, 4, 3, 3, 3, 2, 2, 0])
            sage: S[11]
            0
            sage: S.index(0)
            11

        """
        cdef int out
        if not isinstance(other, BoundedIntegerSequence):
            if other<0:
                raise ValueError("BoundedIntegerSequence.index(x): x(={}) not in sequence".format(other))
            try:
                out = index_biseq(self.data, <mp_limb_t>other, 0)
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
            sage: BoundedIntegerSequence(21, [4,1,6,2,7,20,9,4]) in S+T
            True
            sage: T+list(S)
            Traceback (most recent call last):
            ...
            TypeError:  Cannot convert list to sage.misc.bounded_integer_sequences.BoundedIntegerSequence
            sage: T+None
            Traceback (most recent call last):
            ...
            TypeError: Cannot concatenate bounded integer sequence and None

        TESTS:

        The discussion at :trac:`15820` explains why the following is a good test::

            sage: BoundedIntegerSequence(21, [0,0]) + BoundedIntegerSequence(21, [0,0])
            <0, 0, 0, 0>

        """
        cdef BoundedIntegerSequence myself, right, out
        if other is None or self is None:
            raise TypeError('Cannot concatenate bounded integer sequence and None')
        myself = self  # may result in a type error
        right = other  #  --"--
        if right.data.itembitsize!=myself.data.itembitsize:
            raise ValueError("can only concatenate bounded integer sequences of compatible bounds")
        out = BoundedIntegerSequence.__new__(BoundedIntegerSequence, 0, None)
        concat_biseq(out.data, myself.data, right.data)
        return out

    cpdef BoundedIntegerSequence maximal_overlap(self, BoundedIntegerSequence other):
        """
        Returns ``self``'s maximal trailing sub-sequence that ``other`` starts with.

        Returns ``None`` if there is no overlap

        EXAMPLES::

            sage: from sage.misc.bounded_integer_sequences import BoundedIntegerSequence
            sage: X = BoundedIntegerSequence(21, [4,1,6,2,7,2,3])
            sage: S = BoundedIntegerSequence(21, [0,0,0,0,0,0,0])
            sage: T = BoundedIntegerSequence(21, [2,7,2,3,0,0,0,0,0,0,0,1])
            sage: (X+S).maximal_overlap(T)
            <2, 7, 2, 3, 0, 0, 0, 0, 0, 0, 0>
            sage: print (X+S).maximal_overlap(BoundedIntegerSequence(21, [2,7,2,3,0,0,0,0,0,1]))
            None
            sage: (X+S).maximal_overlap(BoundedIntegerSequence(21, [0,0]))
            <0, 0>

        """
        if other.startswith(self):
            return self
        cdef int i = max_overlap_biseq(self.data, other.data)
        if i==-1:
            return None
        return self[i:]

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
        cdef BoundedIntegerSequence Self
        if other is None:
            return 1
        if self is None:
            return -1
        try:
            right = other
        except TypeError:
            return -1
        try:
            Self = self
        except TypeError:
            return 1
        cdef int c = cmp(Self.data.itembitsize, right.data.itembitsize)
        if c:
            return c
        c = cmp(Self.data.length, right.data.length)
        if c:
            return c
        return mpn_cmp(Self.data.data.bits, right.data.data.bits, Self.data.data.limbs)

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
        # This is similar to mpz_pythonhash, which is a very bad additive hash:
        cdef mp_limb_t h = 0
        cdef mp_limb_t h0
        cdef mp_size_t i
        cdef mp_limb_t* p = self.data.data.bits
        for i from self.data.data.limbs>i>=0:
            h0 = h
            h += deref(postinc(p))
            if h<h0: # overflow
                preinc(h)
        if h==-1:
            return -2
        return h
        ## bitset_hash is not a good hash either
        ## We should consider using FNV-1a hash, see http://www.isthe.com/chongo/tech/comp/fnv/,
        ## Or the hash defined in http://burtleburtle.net/bob/hash/doobs.html
        ## Or http://www.azillionmonkeys.com/qed/hash.html

cpdef BoundedIntegerSequence NewBISEQ(tuple bitset_data, mp_bitcnt_t itembitsize, mp_size_t length):
    """
    Helper function for unpickling of :class:`BoundedIntegerSequence`.

    EXAMPLES::

        sage: from sage.misc.bounded_integer_sequences import BoundedIntegerSequence
        sage: L = [randint(0,26) for i in range(5000)]
        sage: S = BoundedIntegerSequence(32, L)
        sage: loads(dumps(S)) == S    # indirect doctest
        True

    """
    cdef BoundedIntegerSequence out = BoundedIntegerSequence.__new__(BoundedIntegerSequence)
    out.data.itembitsize = itembitsize
    out.data.mask_item = ((<mp_limb_t>1)<<itembitsize)-1
    out.data.length = length
    # bitset_unpickle assumes that out.data.data is initialised.
    bitset_init(out.data.data, mp_bits_per_limb)
    bitset_unpickle(out.data.data, bitset_data)
    return out

def _biseq_stresstest():
    """
    This function creates many bounded integer sequences and manipulates them
    in various ways, in order to try to detect random memory corruptions.

    TESTS::

        sage: from sage.misc.bounded_integer_sequences import _biseq_stresstest
        sage: _biseq_stresstest()

    """
    cdef int i
    from sage.misc.prandom import randint
    cdef list L = [BoundedIntegerSequence(6, [randint(0,5) for x in range(randint(4,10))]) for y in range(100)]
    cdef int branch
    cdef BoundedIntegerSequence S, T
    for i from 0<=i<10000:
        branch = randint(0,4)
        if branch == 0:
            L[randint(0,99)] = L[randint(0,99)]+L[randint(0,99)]
        elif branch == 1:
            x = randint(0,99)
            if len(L[x]):
                y = randint(0,len(L[x])-1)
                z = randint(y,len(L[x])-1)
                L[randint(0,99)] = L[x][y:z]
            else:
                L[x] = BoundedIntegerSequence(6, [randint(0,5) for x in range(randint(4,10))])
        elif branch == 2:
            t = list(L[randint(0,99)])
            t = repr(L[randint(0,99)])
            t = L[randint(0,99)].list()
        elif branch == 3:
            x = randint(0,99)
            if len(L[x]):
                y = randint(0,len(L[x])-1)
                t = L[x][y]
                try:
                    t = L[x].index(t)
                except ValueError:
                    raise ValueError("{} should be in {} (bound {}) at position {}".format(t,L[x],L[x].bound(),y))
            else:
                L[x] = BoundedIntegerSequence(6, [randint(0,5) for x in range(randint(4,10))])
        elif branch == 4:
            S = L[randint(0,99)]
            T = L[randint(0,99)]
            startswith_biseq(S.data,T.data)
            contains_biseq(S.data, T.data, 0)
            if S.data.length>0:
                max_overlap_biseq(S.data, T.data)
