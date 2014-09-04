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

    cdef biseq_t* list_to_biseq(list data, unsigned int bound) except NULL
       # Assumes that S is allocated

    cdef list biseq_to_list(biseq_t S)
       # Convert a bounded integer sequence to a list

    cdef biseq_t* concat_biseq(biseq_t S1, biseq_t S2) except NULL
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
    int mpn_cmp (mp_limb_t *s1p, mp_limb_t *s2p, mp_size_t n)

cdef extern from "Python.h":
    bint PySlice_Check(PyObject* ob)

cdef size_t times_size_of_limb = [1,2,4,8,16,32,64].index(sizeof(mp_limb_t))
cdef size_t times_mp_bits_per_limb = [1,2,4,8,16,32,64,128,256].index(mp_bits_per_limb)
cdef size_t mod_mp_bits_per_limb = ((<size_t>1)<<times_mp_bits_per_limb)-1

cdef tuple ZeroNone = (0,None)
cdef PyObject* zeroNone = <PyObject*>ZeroNone
cdef dict EmptyDict = {}
cdef PyObject* emptyDict = <PyObject*>EmptyDict

from cython.operator import dereference as deref, preincrement as preinc, predecrement as predec

###################
#  Boilerplate
#   cdef functions
###################

#
# (De)allocation, copying
#

cdef biseq_t allocate_biseq(size_t l, unsigned long int itemsize) except NULL:
    """
    Allocate memory for a bounded integer sequence of length l with items
    fitting in itemsize bits. Returns a pointer to the bounded integer
    sequence, or NULL on error.
    """
    cdef biseq_t out = <biseq_t>sage_malloc(sizeof(biseq))
    if out==NULL:
        raise MemoryError("Can not allocate bounded integer sequence")
    out.bitsize = l*itemsize
    out.length = l
    out.itembitsize = itemsize
    out.mask_item = ((<unsigned int>1)<<itemsize)-1
    sig_on()
    mpz_init2(out.data, out.bitsize+mp_bits_per_limb)
    sig_off()
    return out

cdef void dealloc_biseq(biseq_t S):
    mpz_clear(S.data)
    sage_free(S)

#
# Conversion
#

cdef biseq_t list_to_biseq(list data, unsigned int bound) except NULL:
    """
    Convert a list into a bounded integer sequence.

    """
    cdef mpz_t tmp
    cdef biseq_t S
    if bound!=0:
        mpz_init_set_ui(tmp, bound-1)
        S = allocate_biseq(len(data), mpz_sizeinbase(tmp, 2))
        mpz_clear(tmp)
    else:
        raise ValueError("The bound for the items of bounded integer sequences must be positive")
    cdef unsigned long int item
    cdef mp_limb_t item_limb
    cdef mp_limb_t tmp_limb
    cdef int offset_mod, offset_div
    cdef int offset = 0
    if not data:
        return S
    for item in data:
        item_limb = <mp_limb_t>(item&S.mask_item)
        offset_mod = offset&mod_mp_bits_per_limb
        offset_div = offset>>times_mp_bits_per_limb
        if offset_mod:
            assert offset_div+1<(<__mpz_struct*>S.data)._mp_alloc
            (<__mpz_struct*>S.data)._mp_d[offset_div+1] = mpn_lshift(&tmp_limb, &item_limb, 1, offset_mod)
            (<__mpz_struct*>S.data)._mp_d[offset_div] |= tmp_limb
            if (<__mpz_struct*>S.data)._mp_d[offset_div+1]!=0 and offset_div+1>=(<__mpz_struct*>S.data)._mp_size:
                (<__mpz_struct*>S.data)._mp_size = offset_div+2
            elif tmp_limb!=0 and offset_div>=(<__mpz_struct*>S.data)._mp_size:
                (<__mpz_struct*>S.data)._mp_size = offset_div+1
        else:
            assert offset_div<(<__mpz_struct*>S.data)._mp_alloc
            (<__mpz_struct*>S.data)._mp_d[offset_div] = item_limb
            if item_limb!=0 and offset_div>=(<__mpz_struct*>S.data)._mp_size:
                (<__mpz_struct*>S.data)._mp_size = offset_div+1
        offset += S.itembitsize
    return S

cdef list biseq_to_list(biseq_t S):
    """
    Convert a bounded integer sequence to a list of integers.
    """
    cdef int index, limb_index, bit_index, max_limb
    index = 0
    cdef list L = []
    cdef size_t n
    # If limb_index<max_limb, then we safely are in allocated memory.
    # If limb_index==max_limb, then we must only consider *one* limb (as it is the last)
    # Otherwise, GMP was clever and has removed trailing zeroes.
    max_limb = (<__mpz_struct*>S.data)._mp_size - 1
    if max_limb<0:
        return [int(0)]*S.length
    cdef mp_limb_t tmp_limb[2]
    for n from S.length>=n>0:
        limb_index = index>>times_mp_bits_per_limb
        bit_index  = index&mod_mp_bits_per_limb
        # Problem: GMP tries to be clever, and sometimes removes trailing
        # zeroes. Hence, we must take care to not read bits from non-allocated
        # memory.
        if limb_index < max_limb:
            if bit_index:
                if bit_index+S.itembitsize >= mp_bits_per_limb:
                    mpn_rshift(<mp_limb_t*>tmp_limb, (<__mpz_struct*>S.data)._mp_d+limb_index, 2, bit_index)
                    L.append(<int>(tmp_limb[0]&S.mask_item))
                else:
                    L.append(<int>((((<__mpz_struct*>S.data)._mp_d[limb_index])>>bit_index)&S.mask_item))
            else:
                L.append(<int>((<__mpz_struct*>S.data)._mp_d[limb_index]&S.mask_item))
        elif limb_index == max_limb:
            if bit_index:
                L.append(<int>((((<__mpz_struct*>S.data)._mp_d[limb_index])>>bit_index)&S.mask_item))
            else:
                L.append(<int>((<__mpz_struct*>S.data)._mp_d[limb_index]&S.mask_item))
        else:
            break
        index += S.itembitsize
    if n:
        L.extend([int(0)]*n)
    return L

#
# Arithmetics
#

cdef biseq_t concat_biseq(biseq_t S1, biseq_t S2) except NULL:
    """
    Concatenate two bounded integer sequences.

    ASSUMPTION:

    - The two sequences must have equivalent bounds, i.e., the items on the
      sequences must fit into the same number of bits. This condition is not
      tested.

    OUTPUT:

    - A pointer to the concatenated sequence, or NULL on error.

    """
    cdef biseq_t out = <biseq_t>sage_malloc(sizeof(biseq))
    if out==NULL:
        raise MemoryError("Can not allocate bounded integer sequence")
    out.bitsize = S1.bitsize+S2.bitsize
    out.length = S1.length+S2.length
    out.itembitsize = S1.itembitsize # do not test == S2.itembitsize
    out.mask_item = S1.mask_item
    sig_on()
    mpz_init2(out.data, out.bitsize)
    sig_off()
    mpz_mul_2exp(out.data, S2.data, S1.bitsize)
    mpz_ior(out.data, out.data, S1.data)
    return out

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
    cdef __mpz_struct seq
    cdef int limb_size, limb_size_orig
    cdef mpz_t tmp
    # We will use tmp to store enough bits to be able to compare with S2.
    # Hence, we say that its _mp_size is equal to that of S2. This may not be
    # correct, as we move bits from S1 (hence, there might be different
    # trailing zeroes), but this is not critical for the mpz_* functions we
    # are using here: _mp_size must be big enough, but can be bigger than
    # needed.
    #
    # Problem: When bitshifting data of S1, it may happen that S1 has some
    # cut-off trailing zeroes. In this case, the bit-shift will not fill all
    # allocated memory, and thus the highest limbs of tmp may still be
    # non-zero from previous bitshifts. If the highest limb of S2 happens to
    # be equal to this non-zero artifact of previous bitshift operations on
    # S1, then the congruence tests will give the wrong result (telling that
    # S2 in S1 even though the last items of S2 are non-zero and the
    # corresponding items of S1 are zero.
    #
    # Solution: We only shift the limbs of S1 that are guaranteed to be
    # correct. Then, we compare S2 with the shifted version of S1 at most up
    # to the number of valid limbs. If S2's _mp_size is not larger than that,
    # then it is a sub-sequence; otherwise (provided that we are keeping track
    # of _mp_size correctly!!!) it has non-zero bits where S1 has (cut-off)
    # trailing zeroes, and can thus be no sub-sequence (also not for higher
    # index).
    seq = deref(<__mpz_struct*>S2.data)
    # The number of limbs required to store data of S2's bitsize (including
    # trailing zeroes)
    limb_size = (S2.bitsize>>times_mp_bits_per_limb)+1
    # Size of S2 without trailing zeroes:
    limb_size_orig = seq._mp_size

    seq = deref(<__mpz_struct*>S1.data)
    cdef int n, limb_index, bit_index
    # Idea: We shift-copy enough limbs from S1 to tmp and then compare with
    # S2, for each shift.
    sig_on()
    mpz_init2(tmp, S2.bitsize+mp_bits_per_limb)
    sig_off()
    # The maximal number of limbs of S1 that is safe to use after the current
    # position, but not more than what fits into tmp:
    cdef int max_S1_limbs = seq._mp_size-((start*S1.itembitsize)>>times_mp_bits_per_limb)
    cdef int max_limb_size = min((<__mpz_struct*>tmp)._mp_alloc, max_S1_limbs)
    if ((start*S1.itembitsize)&mod_mp_bits_per_limb)==0:
        # We increase it, since in the loop it is decreased when the index is
        # zero modulo bits per limb
        preinc(max_S1_limbs)
    n = 0
    cdef int index = 0
    # tmp may have trailing zeroes, and GMP can NOT cope with that!
    # Hence, we need to adjust the size later.
    (<__mpz_struct*>tmp)._mp_size = limb_size
    for index from 0<=index<=S1.length-S2.length:
        limb_index = n>>times_mp_bits_per_limb
        # Adjust the number of limbs we are shifting
        bit_index  = n&mod_mp_bits_per_limb
        if bit_index:
            if limb_size < max_limb_size:
                assert limb_size<(<__mpz_struct*>tmp)._mp_alloc
                mpn_rshift((<__mpz_struct*>tmp)._mp_d, seq._mp_d+limb_index, limb_size+1, bit_index)
            else:
                assert max_limb_size<=(<__mpz_struct*>tmp)._mp_alloc
                (<__mpz_struct*>tmp)._mp_size = max_limb_size
                mpn_rshift((<__mpz_struct*>tmp)._mp_d, seq._mp_d+limb_index, max_limb_size, bit_index)
        else:
            # The bit_index is zero, and hence one limb less is remaining. We
            # thus decrement max_S1_limbs.
            predec(max_S1_limbs)
            max_limb_size = min((<__mpz_struct*>tmp)._mp_alloc, max_S1_limbs)
            if limb_size < max_limb_size:
                assert limb_size<=(<__mpz_struct*>tmp)._mp_alloc
                mpn_copyi((<__mpz_struct*>tmp)._mp_d, seq._mp_d+limb_index, limb_size)
            else:
                assert max_limb_size<=(<__mpz_struct*>tmp)._mp_alloc
                (<__mpz_struct*>tmp)._mp_size = max_limb_size
                mpn_copyi((<__mpz_struct*>tmp)._mp_d, seq._mp_d+limb_index, max_limb_size)
        # Now, the first min(limb_size,max_limb_size) limbs of tmp are
        # guaranteed to be valid. The rest may be junk.

        # If the number of valid limbs of tmp is smaller than limb_size_orig
        # (this can only happen if max_limb_size became small), then we were
        # running into trailing zeroes of S1. Hence, as explained above, we
        # can decide now whether S2 is a sub-sequence of S1 in the case that
        # the non-cutoff bits match.
        if max_limb_size<limb_size_orig:
            # We have exactly max_limb_size limbs to compare
            if mpn_cmp((<__mpz_struct*>tmp)._mp_d, (<__mpz_struct*>(S2.data))._mp_d, max_limb_size)==0:
                mpz_clear(tmp)
                if limb_size_orig > max_limb_size:
                    # S2 has further non-zero entries
                    return -1
                # Both S2 and S1 have enough trailing zeroes
                return index
        else:
            if mpz_congruent_2exp_p(S2.data, tmp, S2.bitsize):
                mpz_clear(tmp)
                return index
        n += S1.itembitsize
    mpz_clear(tmp)
    return -1

cdef int max_overlap_biseq(biseq_t S1, biseq_t S2) except 0:
    """
    Returns the smallest **positive** integer ``i`` such that ``S2`` starts with ``S1[i:]``.

    Returns ``-1`` if there is no overlap. Note that ``i==0`` will not be considered!

    INPUT:

    - ``S1``, ``S2`` -- two bounded integer sequences

    ASSUMPTION:

    - The two sequences must have equivalent bounds, i.e., the items on the
      sequences must fit into the same number of bits. This condition is not
      tested.

    """
    cdef __mpz_struct seq
    seq = deref(<__mpz_struct*>S1.data)
    cdef mpz_t tmp
    # Idea: We shift-copy enough limbs from S1 to tmp and then compare with
    # the initial part of S2, for each shift.
    sig_on()
    mpz_init2(tmp, S1.bitsize+mp_bits_per_limb)
    sig_off()
    # We will use tmp to store enough bits to be able to compare with the tail
    # of S1.
    cdef int n, limb_index, bit_index
    cdef int index, i, start_index
    if S2.length>=S1.length:
        start_index = 1
    else:
        start_index = S1.length-S2.length
    n = S1.itembitsize*start_index
    for index from start_index<=index<S1.length:
        limb_index = n>>times_mp_bits_per_limb
        if limb_index>=seq._mp_size:
            break
        bit_index  = n&mod_mp_bits_per_limb
        if bit_index:
            assert seq._mp_size-limb_index<=(<__mpz_struct*>tmp)._mp_alloc
            mpn_rshift((<__mpz_struct*>tmp)._mp_d, seq._mp_d+limb_index, seq._mp_size-limb_index, bit_index)
        else:
            assert seq._mp_size-limb_index<=(<__mpz_struct*>tmp)._mp_alloc
            mpn_copyi((<__mpz_struct*>tmp)._mp_d, seq._mp_d+limb_index, seq._mp_size-limb_index)
        (<__mpz_struct*>tmp)._mp_size = seq._mp_size-limb_index
        if mpz_congruent_2exp_p(tmp, S2.data, S1.bitsize-n):
            mpz_clear(tmp)
            return index
        n += S1.itembitsize
    # now, it could be that we have to check for leading zeroes on S2.
    if index<S1.length:
        mpz_clear(tmp)
        mpz_init_set_ui(tmp, 0)
        for i from index <= i < S1.length:
            if mpz_congruent_2exp_p(tmp, S2.data, S1.bitsize-n):
                mpz_clear(tmp)
                return i
            n += S1.itembitsize
    mpz_clear(tmp)
    return -1

cdef int index_biseq(biseq_t S, int item, size_t start) except -2:
    """
    Returns the position in S of an item in S[start:], or -1 if S[start:] does
    not contain the item.

    """
    cdef __mpz_struct seq
    seq = deref(<__mpz_struct*>S.data)

    cdef int n, limb_index, bit_index, max_limb
    cdef mp_limb_t tmp_limb[2]
    item &= S.mask_item
    n = 0
    # If limb_index<max_limb, then we safely are in allocated memory.
    # If limb_index==max_limb, then we must only consider *one* limb (as it is the last)
    # Otherwise, GMP was clever and has removed trailing zeroes.
    max_limb = seq._mp_size - 1
    cdef int index = 0
    for index from 0<=index<S.length:
        limb_index = n>>times_mp_bits_per_limb
        bit_index  = n&mod_mp_bits_per_limb
        if limb_index < max_limb:
            if bit_index:
                if bit_index+S.itembitsize >= mp_bits_per_limb:
                    mpn_rshift(tmp_limb, seq._mp_d+limb_index, 2, bit_index)
                    if item==(tmp_limb[0]&S.mask_item):
                        return index
                elif item==(((seq._mp_d[limb_index])>>bit_index)&S.mask_item):
                    return index
            elif item==(seq._mp_d[limb_index]&S.mask_item):
                return index
        elif limb_index == max_limb:
            if bit_index:
                if item==(((seq._mp_d[limb_index])>>bit_index)&S.mask_item):
                    return index
            elif item==(seq._mp_d[limb_index]&S.mask_item):
                return index
        else:
            if item==0 and index<S.length:
                return index
            else:
                return -1
        n += S.itembitsize
    return -1

cdef int getitem_biseq(biseq_t S, unsigned long int index) except -1:
    """
    Get item S[index], without checking margins.

    """
    cdef __mpz_struct seq
    seq = deref(<__mpz_struct*>S.data)
    cdef long int limb_index, bit_index, limb_size
    index *= S.itembitsize
    limb_index = index>>times_mp_bits_per_limb
    bit_index  = index&mod_mp_bits_per_limb
    cdef mp_limb_t tmp_limb[2]
    cdef int out
    if limb_index>=seq._mp_size:
        return 0
    if bit_index:
        # limb_size is 1 or 2
        limb_size = ((bit_index+S.itembitsize-1)>>times_mp_bits_per_limb)+1
        if limb_index+limb_size>seq._mp_size:
            limb_size = 1 # the second limb may be non-allocated memory and
                          # should be treated as zero.
        if limb_size!=1:
            assert limb_size<=2
            mpn_rshift(tmp_limb, seq._mp_d+limb_index, limb_size, bit_index)
            out = <int>(tmp_limb[0])
        else:
            out = <int>((seq._mp_d[limb_index])>>bit_index)
    else:
        out = <int>(seq._mp_d[limb_index])
    return out&S.mask_item

cdef biseq_t slice_biseq(biseq_t S, int start, int stop, int step) except NULL:
    """
    Create the slice S[start:stop:step] as bounded integer sequence.

    Return:

    - A pointer to the resulting bounded integer sequence, or NULL on error.

    """
    cdef long int length, total_shift, n
    if step>0:
        if stop>start:
            length = ((stop-start-1)//step)+1
        else:
            return allocate_biseq(0, S.itembitsize)
    else:
        if stop>=start:
            return allocate_biseq(0, S.itembitsize)
        else:
            length = ((stop-start+1)//step)+1
    cdef biseq_t out = allocate_biseq(length, S.itembitsize)
    cdef int offset_mod, offset_div
    cdef int limb_size
    total_shift = start*S.itembitsize
    if step==1:
        # allocate_biseq allocates one limb more than strictly needed. This is
        # enough to shift a limb* array directly to its ._mp_d field.
        offset_mod = total_shift&mod_mp_bits_per_limb
        offset_div = total_shift>>times_mp_bits_per_limb
        limb_size = ((offset_mod+out.bitsize-1)>>times_mp_bits_per_limb)+1
        # GMP tries to be clever and thus may cut trailing zeroes. Hence,
        # limb_size may be so large that it would access non-allocated
        # memory. We adjust this now:
        limb_size = min(limb_size, (<__mpz_struct*>(S.data))._mp_size-offset_div)
        if limb_size<=0:
            (<__mpz_struct*>(out.data))._mp_size = 0
            return out
        assert limb_size<=(<__mpz_struct*>(out.data))._mp_alloc
        mpn_rshift((<__mpz_struct*>(out.data))._mp_d, (<__mpz_struct*>(S.data))._mp_d+offset_div, limb_size, offset_mod)
        for n from limb_size>n>=0:
            if (<__mpz_struct*>(out.data))._mp_d[n]:
                break
        (<__mpz_struct*>(out.data))._mp_size = n+1
        mpz_fdiv_r_2exp(out.data, out.data, out.bitsize)
        return out

    # Convention for variable names: We move data from *_src to *_tgt
    # tmp_limb is used to temporarily store the data that we move/shift
    cdef mp_limb_t tmp_limb[2]
    cdef int offset_mod_src, offset_div_src, max_limb
    cdef int offset_src = total_shift
    cdef __mpz_struct *seq_src
    seq_src = <__mpz_struct*>S.data
    # We need to take into account that GMP tries to be clever, and removes
    # trailing zeroes from the allocated memory.
    max_limb = seq_src._mp_size - 1
    cdef int bitstep = step*S.itembitsize
    cdef int offset_mod_tgt, offset_div_tgt
    cdef int offset_tgt = 0
    cdef __mpz_struct *seq_tgt
    seq_tgt = <__mpz_struct*>out.data
    for n from length>=n>0:
        offset_div_src = offset_src>>times_mp_bits_per_limb
        offset_mod_src = offset_src&mod_mp_bits_per_limb
        offset_div_tgt = offset_tgt>>times_mp_bits_per_limb
        offset_mod_tgt = offset_tgt&mod_mp_bits_per_limb
        # put data from src to tmp_limb[0].
        if offset_div_src < max_limb:
            if offset_mod_src:
                if (offset_mod_src+S.itembitsize >= mp_bits_per_limb):
                    mpn_rshift(tmp_limb, seq_src._mp_d+offset_div_src, 2, offset_mod_src)
                    tmp_limb[0] &= S.mask_item
                else:
                    tmp_limb[0] = ((seq_src._mp_d[offset_div_src])>>offset_mod_src) & S.mask_item
            else:
                tmp_limb[0] = seq_src._mp_d[offset_div_src] & S.mask_item
        elif offset_div_src == max_limb:
            if offset_mod_src:
                tmp_limb[0] = ((seq_src._mp_d[offset_div_src])>>offset_mod_src) & S.mask_item
            else:
                tmp_limb[0] = seq_src._mp_d[offset_div_src] & S.mask_item
        else:
            tmp_limb[0] = 0
        # put data from tmp_limb[0] to tgt
        if offset_mod_tgt:
            assert offset_div_tgt+1<(<__mpz_struct*>(out.data))._mp_alloc
            seq_tgt._mp_d[offset_div_tgt+1] = mpn_lshift(tmp_limb, tmp_limb, 1, offset_mod_tgt)
            seq_tgt._mp_d[offset_div_tgt]  |= tmp_limb[0]
            if seq_tgt._mp_d[offset_div_tgt+1]!=0 and offset_div_tgt+1>=seq_tgt._mp_size:
                seq_tgt._mp_size = offset_div_tgt+2
            elif tmp_limb[0]!=0 and offset_div_tgt>=seq_tgt._mp_size:
                seq_tgt._mp_size = offset_div_tgt+1
        else:
            assert offset_div_tgt<(<__mpz_struct*>(out.data))._mp_alloc
            seq_tgt._mp_d[offset_div_tgt] = tmp_limb[0]
            if tmp_limb[0]!=0:
                seq_tgt._mp_size = offset_div_tgt+1
        offset_tgt += S.itembitsize
        offset_src += bitstep
    return out

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
        self.data = NULL

    def __dealloc__(self):
        """
        Free the memory from underlying data

        EXAMPLES::

            sage: from sage.misc.bounded_integer_sequences import BoundedIntegerSequence
            sage: S = BoundedIntegerSequence(21, [4,1,6,2,7,20,9])
            sage: del S     # indirect doctest

        """
        if self.data!=NULL:
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
        self.data = list_to_biseq(data, bound)

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
        cdef __mpz_struct seq
        seq = deref(<__mpz_struct*>self.data.data)
        cdef int index, limb_index, bit_index
        cdef int max_limb
        index = 0
        cdef mp_limb_t tmp_limb[2]
        cdef size_t n, n2
        # If limb_index<max_limb, then we safely are in allocated memory.
        # If limb_index==max_limb, then we must only consider *one* limb (as it is the last)
        # Otherwise, GMP was clever and has removed trailing zeroes.
        max_limb = seq._mp_size - 1
        if max_limb<0:
            for n from self.data.length>=n>0:
                yield <int>0
            return
        for n from self.data.length>=n>0:
            limb_index = index>>times_mp_bits_per_limb
            bit_index  = index&mod_mp_bits_per_limb
            if limb_index < max_limb:
                if bit_index:
                    if bit_index+self.data.itembitsize >= mp_bits_per_limb:
                        mpn_rshift(tmp_limb, seq._mp_d+limb_index, 2, bit_index)
                        yield <int>(tmp_limb[0]&self.data.mask_item)
                    else:
                        yield <int>(((seq._mp_d[limb_index])>>bit_index)&self.data.mask_item)
                else:
                    yield <int>(seq._mp_d[limb_index]&self.data.mask_item)
            elif limb_index == max_limb:
                yield <int>(((seq._mp_d[limb_index])>>bit_index)&self.data.mask_item)
            else:
                for n2 from n>=n2>0:
                    yield <int>0
                break
            index += self.data.itembitsize

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

        """
        cdef BoundedIntegerSequence out
        cdef int start,stop,step
        if PySlice_Check(<PyObject *>index):
            start,stop,step = index.indices(self.data.length)
            if start==0 and stop==self.data.length and step==1:
                return self
            out = BoundedIntegerSequence.__new__(BoundedIntegerSequence, 0, None)
            out.data = slice_biseq(self.data, start, stop, step)
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
        out.data = concat_biseq(myself.data, right.data)
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
    out.data = <biseq_t>sage_malloc(sizeof(biseq))
    if out.data==NULL:
        raise MemoryError("Can not allocate bounded integer sequence")
    mpz_init2(out.data.data, bitsize+mp_bits_per_limb)
    out.data.bitsize = bitsize
    out.data.itembitsize = itembitsize
    out.data.mask_item = ((<unsigned int>1)<<itembitsize)-1
    out.data.length = length
    mpz_set_str(out.data.data, data, 32)
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
            max_overlap_biseq(S.data, T.data)
