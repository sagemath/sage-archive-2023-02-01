"""
Boilerplate functions for a cython implementation of elements of path algebras.

AUTHORS:

- Simon King (2014-12-04)

"""

#*****************************************************************************
#     Copyright (C) 2014 Simon King <simon.king@uni-jena.de>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#  as published by the Free Software Foundation; either version 2 of
#  the License, or (at your option) any later version.
#                  http://www.gnu.org/licenses/
#*****************************************************************************

include "sage/ext/interrupt.pxi"
include "sage/ext/stdsage.pxi"
include "sage/data_structures/bitset.pxi"

from cpython.ref cimport *
from cython.operator cimport preincrement as preinc, predecrement as predec, postincrement as postinc
from sage.libs.gmp.mpn cimport mpn_cmp

########################################
##
## Allocation and Deallocation of monomials
#
# Monomials are expensive, hence, copying will just be done by increasing a
# reference counter.

# Create a monomial by copying the given bounded integer sequence
cdef inline path_mon_t *mon_create(biseq_t Mon, unsigned int Pos, int Mid, mp_size_t S_len) except NULL:
    if Mid == -1:
        assert Pos==S_len==0
    cdef path_mon_t *out = <path_mon_t*>sage_malloc(sizeof(path_mon_t))
    if out==NULL:
        raise MemoryError("Out of memory while allocating a path term")
    biseq_init_copy(out.path, Mon)
    out.pos = Pos
    out.mid = Mid
    out.s_len = S_len
    out.ref = 1
    return out

# the following is only used in the free-list for terms.
cdef inline bint mon_realloc(path_mon_t *out, biseq_t Mon, unsigned int Pos, int Mid, mp_size_t S_len) except -1:
    if Mid == -1:
        assert Pos==S_len==0
    bitset_realloc(out.path.data, Mon.data.size)
    bitset_copy(out.path.data, Mon.data)
    out.path.length = Mon.length
    out.path.itembitsize = Mon.itembitsize
    out.path.mask_item = Mon.mask_item
    out.pos = Pos
    out.mid = Mid
    out.s_len = S_len
    out.ref = 1
    return True

# Create a monomial without copying the given bounded integer sequence
cdef inline path_mon_t *mon_create_keep(biseq_t Mon, unsigned int Pos, int Mid, mp_size_t S_len) except NULL:
    if Mid == -1:
        assert Pos==S_len==0
    cdef path_mon_t *out = <path_mon_t*>sage_malloc(sizeof(path_mon_t))
    if out==NULL:
        raise MemoryError("Out of memory while allocating a path term")
    out.path[0] = Mon[0]
    out.pos = Pos
    out.mid = Mid
    out.s_len = S_len
    out.ref = 1
    return out

cdef inline bint mon_realloc_keep(path_mon_t *out, biseq_t Mon, unsigned int Pos, int Mid, mp_size_t S_len) except -1:
    if Mid == -1:
        assert Pos==S_len==0
    biseq_dealloc(out.path)
    out.path[0] = Mon[0]
    out.pos = Pos
    out.mid = Mid
    out.s_len = S_len
    out.ref = 1
    return True

# Copying a monomial just means to create another reference
cdef inline path_mon_t *mon_copy(path_mon_t *M):
    preinc(M.ref)
    return M

# Deallocate the monomial
cdef inline void mon_free(path_mon_t *M):
    if predec(M.ref):
        return
    biseq_dealloc(M.path)
    sage_free(M)

# Linearisation
cdef inline tuple mon_pickle(path_mon_t *M):
    return (bitset_pickle(M.path.data) if M.path.length>0 else (),
            M.path.itembitsize, M.path.length, M.pos, M.mid, M.s_len)

cdef path_mon_t *mon_unpickle(tuple data) except NULL:
    cdef tuple bitset_data
    cdef mp_bitcnt_t itembitsize
    cdef mp_size_t length
    cdef unsigned int Pos
    cdef int Mid
    cdef mp_size_t S_len
    bitset_data, itembitsize, length, Pos, Mid, S_len = data
    cdef path_mon_t *out = <path_mon_t*>sage_malloc(sizeof(path_mon_t))
    if out==NULL:
        raise MemoryError("Out of memory while allocating a path term")
    out.path.itembitsize = itembitsize
    out.path.mask_item = limb_lower_bits_up(itembitsize)
    out.path.length = length

    # bitset_unpickle assumes that out.path.data is initialised.
    bitset_init(out.path.data, GMP_LIMB_BITS)
    if bitset_data: bitset_unpickle(out.path.data, bitset_data)
    out.pos = Pos
    out.mid = Mid
    out.s_len = S_len
    out.ref = 1
    return out


########################################
##
## Monomial orders

cdef inline int negdegrevlex(path_mon_t *M1, path_mon_t *M2):
    # a*s_i*b<c*s_j*d <=>
    # 1. deg(a*b) > deg(c*d), otherwise
    # 2. deg(a) > deg(c) (note that one of them may be -1), otherwise
    # 3. deg(s_i) < deg(s_j), otherwise
    # 4. a*s_i*b <_revlex c*s_j*d, otherwise
    # 5. i<j
    cdef mp_size_t l1 = M1.path.length + M2.s_len # sic!
    cdef mp_size_t l2 = M2.path.length + M1.s_len
    if l2<l1:
        return -1
    if l1<l2:
        return 1
    if M2.mid < M1.mid:
        return -1
    if M1.mid < M2.mid:
        return 1
    if M1.s_len < M2.s_len:
        return -1
    if M2.s_len < M1.s_len:
        return 1
    # mpn_cmp does comparison of long integers. If the two long integers have
    # the same number of digits (this is the case her), it is the same as
    # lexicographic comparison of the numbers. The highest digit corresponds
    # to the right-most item in the path. Hence, it becomes
    # reverse-lexicographic order.
    cdef int c = mpn_cmp(M1.path.data.bits, M2.path.data.bits, M1.path.data.limbs)
    if c!=0:
        return c
    if M1.pos < M2.pos:
        return -1
    if M2.pos < M1.pos:
        return 1
    return 0

cdef inline int degrevlex(path_mon_t *M1, path_mon_t *M2):
    # a*s_i*b<c*s_j*d <=>
    # 1. deg(a*b) < deg(c*d), otherwise
    # 2. deg(a) < deg(c) (note that one of them may be -1), otherwise
    # 3. deg(s_i) > deg(s_j), otherwise
    # 4. a*s_i*b <_revlex c*s_j*d, otherwise
    # 5. i<j
    cdef mp_size_t l1 = M1.path.length + M2.s_len # sic!
    cdef mp_size_t l2 = M2.path.length + M1.s_len
    if l2<l1:
        return 1
    if l1<l2:
        return -1
    if M2.mid < M1.mid:
        return 1
    if M1.mid < M2.mid:
        return -1
    if M1.s_len < M2.s_len:
        return 1
    if M2.s_len < M1.s_len:
        return -1
    # mpn_cmp does comparison of long integers. If the two long integers have
    # the same number of digits (this is the case her), it is the same as
    # lexicographic comparison of the numbers. The highest digit corresponds
    # to the right-most item in the path. Hence, it becomes
    # reverse-lexicographic order.
    cdef int c = mpn_cmp(M1.path.data.bits, M2.path.data.bits, M1.path.data.limbs)
    if c!=0:
        return c
    if M1.pos < M2.pos:
        return -1
    if M2.pos < M1.pos:
        return 1
    return 0

cdef inline int negdeglex(path_mon_t *M1, path_mon_t *M2):
    # a*s_i*b<c*s_j*d <=>
    # 1. deg(a*b) > deg(c*d), otherwise
    # 2. deg(a) > deg(c) (note that one of them may be -1), otherwise
    # 3. deg(s_i) < deg(s_j), otherwise
    # 4. a*s_i*b <_lex c*s_j*d, otherwise
    # 5. i<j
    cdef mp_size_t l1 = M1.path.length + M2.s_len # sic!
    cdef mp_size_t l2 = M2.path.length + M1.s_len
    cdef size_t item1, item2
    if l2<l1:
        return -1
    if l1<l2:
        return 1
    if M2.mid < M1.mid:
        return -1
    if M1.mid < M2.mid:
        return 1
    if M1.s_len < M2.s_len:
        return -1
    if M2.s_len < M1.s_len:
        return 1
    for index from 0 <= index < M1.path.length:
        item1 = biseq_getitem(M1.path, index)
        item2 = biseq_getitem(M2.path, index)
        if item1<item2:
            return -1
        if item2<item1:
            return 1
    if M1.pos < M2.pos:
        return -1
    if M2.pos < M1.pos:
        return 1
    return 0

cdef inline int deglex(path_mon_t *M1, path_mon_t *M2):
    # a*s_i*b<c*s_j*d <=>
    # 1. deg(a*b) < deg(c*d), otherwise
    # 2. deg(a) < deg(c) (note that one of them may be -1), otherwise
    # 3. deg(s_i) > deg(s_j), otherwise
    # 4. a*s_i*b <_lex c*s_j*d, otherwise
    # 5. i<j
    cdef mp_size_t l1 = M1.path.length + M2.s_len # sic!
    cdef mp_size_t l2 = M2.path.length + M1.s_len
    cdef size_t item1, item2
    if l2<l1:
        return 1
    if l1<l2:
        return -1
    if M2.mid < M1.mid:
        return 1
    if M1.mid < M2.mid:
        return -1
    if M1.s_len < M2.s_len:
        return 1
    if M2.s_len < M1.s_len:
        return -1
    for index from 0 <= index < M1.path.length:
        item1 = biseq_getitem(M1.path, index)
        item2 = biseq_getitem(M2.path, index)
        if item1<item2:
            return -1
        if item2<item1:
            return 1
    if M1.pos < M2.pos:
        return -1
    if M2.pos < M1.pos:
        return 1
    return 0

########################################
##
## Allocation and Deallocation of terms
###########################
# We use a kill list.
cdef path_poly_t *kill_list = <path_poly_t*>sage_malloc(sizeof(path_poly_t))
kill_list.nterms = 0
kill_list.lead = NULL

# Deallocate the term, and return the pointer .nxt, without using kill list
cdef inline path_term_t *term_free_force(path_term_t *T):
    mon_free(T.mon)
    cdef path_term_t *out = T.nxt
    sage_free(T)
    return out

cdef class _KillListProtector:
    def __dealloc__(self):
        cdef path_term_t *T = kill_list.lead
        while T != NULL:
            T = term_free_force(T)
        sage_free(kill_list)

cdef size_t poolsize = 1000   # The kill list contains at most that many terms.
_kill_list_protector = _KillListProtector()

# Put the term on the free list, and return the pointer .nxt
cdef inline path_term_t *term_free(path_term_t *T):
    if T.coef!=NULL:
        Py_DECREF(<object>(T.coef))
    cdef path_term_t *out = T.nxt
    if kill_list.nterms < poolsize:
        preinc(kill_list.nterms)
        T.nxt = kill_list.lead
        kill_list.lead = T
        return out
    mon_free(T.mon)
    sage_free(T)
    return out

# Create a term by copying the given bounded integer sequence
cdef inline path_term_t *term_create(object coef, biseq_t Mon, unsigned int Pos, int Mid, mp_size_t S_len) except NULL:
    cdef path_term_t *out
    if kill_list.nterms > 0:
        predec(kill_list.nterms)
        out = kill_list.lead
        kill_list.lead = out.nxt
        if out.mon.ref > 1:
            # This monomial is in use by other terms. It must be preserved!
            predec(out.mon.ref)
            out.mon = mon_create(Mon, Pos, Mid, S_len)
        else:
            mon_realloc(out.mon, Mon, Pos, Mid, S_len)
    else:
        out = <path_term_t*>sage_malloc(sizeof(path_term_t))
        if out==NULL:
            raise MemoryError("Out of memory while allocating a path term")
        out.mon = mon_create(Mon, Pos, Mid, S_len)
    Py_INCREF(coef)
    out.coef = <PyObject*>coef
    out.nxt = NULL
    return out

# Create a term without copying the given bounded integer sequence
cdef inline path_term_t *term_create_keep(object coef, biseq_t Mon, unsigned int Pos, int Mid, mp_size_t S_len) except NULL:
    cdef path_term_t *out
    if kill_list.nterms > 0:
        predec(kill_list.nterms)
        out = kill_list.lead
        kill_list.lead = out.nxt
        if out.mon.ref > 1:
            # This monomial is in use by other terms. It must be preserved!
            predec(out.mon.ref)
            out.mon = mon_create_keep(Mon, Pos, Mid, S_len)
        else:
            mon_realloc_keep(out.mon, Mon, Pos, Mid, S_len)
    else:
        out = <path_term_t*>sage_malloc(sizeof(path_term_t))
        if out==NULL:
            raise MemoryError("Out of memory while allocating a path term")
        out.mon = mon_create_keep(Mon, Pos, Mid, S_len)
    Py_INCREF(coef)
    out.coef = <PyObject*>coef
    #out.nxt = NULL  # to be taken care of externally
    return out

# Create a term from a path_mon_t, which is not copied
cdef inline path_term_t *term_create_keep_mon(object coef, path_mon_t *Mon) except NULL:
    cdef path_term_t *out
    if kill_list.nterms > 0:
        predec(kill_list.nterms)
        out = kill_list.lead
        kill_list.lead = out.nxt
        mon_free(out.mon)
    else:
        out = <path_term_t*>sage_malloc(sizeof(path_term_t))
        if out==NULL:
            raise MemoryError("Out of memory while allocating a path term")
    Py_INCREF(coef)
    out.coef = <PyObject*>coef
    out.mon = Mon
    #out.nxt = NULL  # to be taken care of externally
    return out

######################################################################
######################################################################

cdef inline path_term_t *term_copy(path_term_t *T) except NULL:
    cdef path_term_t *out
    if kill_list.nterms > 0:
        predec(kill_list.nterms)
        out = kill_list.lead
        kill_list.lead = out.nxt
        mon_free(out.mon)
    else:
        out = <path_term_t*>sage_malloc(sizeof(path_term_t))
        if out==NULL:
            raise MemoryError("Out of memory while allocating a path term")
    out.mon = mon_copy(T.mon)
    Py_INCREF(<object>T.coef)
    out.coef = T.coef
    # out.nxt is supposed to be taken care of externally
    return out

cdef inline path_term_t *term_copy_scale(path_term_t *T, object coef) except? NULL:
    cdef path_term_t *out
    cdef outcoef = coef*<object>(T.coef)
    if not outcoef:
        return NULL
    if kill_list.nterms > 0:
        predec(kill_list.nterms)
        out = kill_list.lead
        kill_list.lead = out.nxt
        mon_free(out.mon)
    else:
        out = <path_term_t*>sage_malloc(sizeof(path_term_t))
        if out==NULL:
            raise MemoryError("Out of memory while allocating a path term")
    out.mon = mon_copy(T.mon)
    Py_INCREF(outcoef)
    out.coef = <PyObject*>outcoef
    # out.nxt is supposed to be taken care of externally
    return out

cdef inline path_term_t *term_copy_recursive(path_term_t *T) except NULL:
    cdef path_term_t *out = term_copy(T)
    cdef path_term_t *first = out
    T = T.nxt
    while T!=NULL:
        out.nxt = term_copy(T)
        out = out.nxt
        T = T.nxt
    out.nxt = NULL
    return first

cdef inline long term_hash(path_term_t *T):
    return (<long>hash(<object>T.coef)+(T.mon.mid<<5)+(T.mon.pos<<10))^bitset_hash(T.mon.path.data)

cdef inline mp_size_t term_total_degree(path_term_t *T):
    return T.mon.path.length-T.mon.s_len

cdef inline tuple term_pickle(path_term_t *T):
    return (<object>T.coef, mon_pickle(T.mon))

cdef inline path_term_t *term_unpickle(object coef, tuple mon_data) except NULL:
    return term_create_keep_mon(coef, mon_unpickle(mon_data))

########################################
##
## Multiplication of monomials

cdef inline path_mon_t *mon_mul_path(path_mon_t *T, biseq_t p) except NULL:
    if p.length == 0:
        return mon_copy(T)
    cdef path_mon_t *out = <path_mon_t*>sage_malloc(sizeof(path_mon_t))
    if out==NULL:
        raise MemoryError("Out of memory while allocating a path monomial")
    biseq_init_concat(out.path, T.path, p)
    out.pos = T.pos
    out.mid = T.mid
    out.s_len = T.s_len
    out.ref = 1
    return out

cdef inline path_mon_t *path_mul_mon(biseq_t p, path_mon_t *T) except NULL:
    if p.length == 0:
        return mon_copy(T)
    cdef path_mon_t *out = <path_mon_t*>sage_malloc(sizeof(path_mon_t))
    if out==NULL:
        raise MemoryError("Out of memory while allocating a path monomial")
    biseq_init_concat(out.path, p, T.path)
    out.pos = T.pos
    out.mid = -1 if T.mid==-1 else T.mid+p.length
    out.s_len = T.s_len
    out.ref = 1
    return out

cdef path_mon_t *path_mul_mon_mul_path(biseq_t p, path_mon_t *T, biseq_t q) except NULL:
    # .mid and .s_len are taken care of externally!
    if p.length==0 and q.length==0:
        return mon_copy(T)
    cdef path_mon_t *out
    if p.length == 0:
        return mon_mul_path(T, q)
    if q.length == 0:
        return path_mul_mon(p, T)
    out = <path_mon_t*>sage_malloc(sizeof(path_mon_t))
    if out==NULL:
        raise MemoryError("Out of memory while allocating a path monomial")
    out.ref = 1
    out.pos = T.pos
    if T.path.length == 0:
        biseq_init_concat(out.path, p, q)
        return out
    cdef mp_size_t pTlength = p.length + T.path.length
    cdef mp_size_t res_length = pTlength + q.length
    biseq_init(out.path, res_length, p.itembitsize)
    if res_length == 0:
        return out
    cdef mp_bitcnt_t pTsize = p.data.size+T.path.data.size
    bitset_lshift(out.path.data, q.data, pTsize)
    cdef mp_bitcnt_t p_offset = p.data.size%GMP_LIMB_BITS
    # p_limbs gives the index of the limb that will store the first bit of the
    # shifted version of T.
    cdef mp_bitcnt_t p_limbs = (p.data.limbs - 1) if p_offset>0 else p.data.limbs

    # pT_limbs gives the index of the last limb used to store p+T
    cdef mp_bitcnt_t pT_limbs = (pTsize-1)//GMP_LIMB_BITS
    if ((T.path.data.size-1)%GMP_LIMB_BITS)+p_offset >= GMP_LIMB_BITS:
        # We shift all limbs of T. The highest bits of the highest limbs are
        # pushed out and returned by mpn_lshift. We need to assign them to the
        # beginning of the last limb that is (partially) occupied by p+T
        out.path.data.bits[pT_limbs] |= mpn_lshift(out.path.data.bits+p_limbs,
                                              T.path.data.bits, T.path.data.limbs, p_offset)
    else:
        if T.path.data.limbs>1:
            # If we would move all limbs of T, then the result would override
            # the lowest limb of the shifted copy of q. We thus only move all
            # but the last limb of T, assigning to the beginning of the last
            # limb of p+T the bits that have been pushed out.
            out.path.data.bits[pT_limbs] |= mpn_lshift(out.path.data.bits+p_limbs,
                                                  T.path.data.bits, T.path.data.limbs-1, p_offset)
            # Last, we need to move the last limb of T (which is only
            # partially occupied), namely into the spot between the previously
            # moved parts of T and the beginning of the shifted copy of q.
            out.path.data.bits[pT_limbs] |= (T.path.data.bits[T.path.data.limbs-1]<<p_offset)
        else:
            out.path.data.bits[p_limbs] |= (T.path.data.bits[T.path.data.limbs-1]<<p_offset)
    bitset_or(out.path.data, out.path.data, p.data)
    return out

########################################
## Addition and scaling of terms

cdef inline path_term_t *term_neg(path_term_t *T) except NULL:
    cdef path_term_t *out
    if kill_list.nterms > 0:
        predec(kill_list.nterms)
        out = kill_list.lead
        kill_list.lead = out.nxt
        mon_free(out.mon)
    else:
        out = <path_term_t*>sage_malloc(sizeof(path_term_t))
        if out==NULL:
            raise MemoryError("Out of memory while allocating a path term")
    cdef object coef = -<object>T.coef
    out.coef = <PyObject*>coef
    Py_INCREF(coef)
    out.mon = mon_copy(T.mon)
    # out.nxt is supposed to be taken care of externally
    return out

cdef inline path_term_t *term_neg_recursive(path_term_t *T) except NULL:
    cdef path_term_t *out = term_neg(T)
    cdef path_term_t *first = out
    T = T.nxt
    while T!=NULL:
        out.nxt = term_neg(T)
        out = out.nxt
        T = T.nxt
    out.nxt = NULL
    return first

cdef inline path_term_t *term_scale(path_term_t *T, object coef) except NULL:
    cdef path_term_t *out
    if kill_list.nterms > 0:
        predec(kill_list.nterms)
        out = kill_list.lead
        kill_list.lead = out.nxt
        mon_free(out.mon)
    else:
        out = <path_term_t*>sage_malloc(sizeof(path_term_t))
        if out==NULL:
            raise MemoryError("Out of memory while allocating a path term")
    cdef object new_coef = coef*<object>T.coef
    if new_coef:
        out.coef = <PyObject*>new_coef
        Py_INCREF(new_coef)
        out.mon = mon_copy(T.mon)
    else:
        out.coef = NULL
    # out.nxt is supposed to be taken care of externally
    return out

cdef inline path_term_t *term_scale_recursive(path_term_t *T, object coef) except NULL:
    cdef path_term_t *out = term_scale(T,coef)
    cdef path_term_t *first = out
    T = T.nxt
    while T!=NULL:
        out.nxt = term_scale(T, coef)
        if out.nxt.coef == NULL:
            #print "zwischendurch frei"
            sage_free(out.nxt)
            out.nxt = NULL
        else:
            out = out.nxt
        T = T.nxt
    out.nxt = NULL
    return first

## Multiplication of terms

cdef path_term_t *term_mul_coef(path_term_t *T, object coef) except NULL:
    cdef path_term_t *out
    if kill_list.nterms > 0:
        predec(kill_list.nterms)
        out = kill_list.lead
        kill_list.lead = out.nxt
        mon_free(out.mon)
    else:
        out = <path_term_t*>sage_malloc(sizeof(path_term_t))
        if out==NULL:
            raise MemoryError("Out of memory while allocating a path term")
    cdef object new_coef = (<object>T.coef)*coef
    if new_coef:
        out.coef = <PyObject*>(new_coef)
        Py_INCREF(new_coef)
        out.mon = mon_copy(T.mon)
    else:
        out.coef = NULL
    out.nxt = NULL
    return out

cdef path_term_t *term_mul_term(path_term_t *T1, path_term_t *T2) except NULL:
    cdef int new_mid
    cdef int new_pos
    cdef mp_size_t new_s_len
    if T1.mon.mid!=-1:
        if T2.mon.mid!=-1:
            raise ValueError("We cannot multiply two module elements")
        new_mid = T1.mon.mid
        new_pos = T1.mon.pos
        new_s_len = T1.mon.s_len
    elif T2.mon.mid!=-1:
        new_mid = T2.mon.mid+T1.mon.path.length
        new_pos = T2.mon.pos
        new_s_len = T2.mon.s_len
    else:
        new_mid = -1
        new_pos = 0
        new_s_len = 0
    cdef object new_coef = (<object>T1.coef)*(<object>T2.coef)

    cdef path_term_t *out
    if kill_list.nterms > 0:
        predec(kill_list.nterms)
        out = kill_list.lead
        kill_list.lead = out.nxt
        if new_coef:
            out.coef = <PyObject*>(new_coef)
            Py_INCREF(new_coef)
            if out.mon.ref>1:
                predec(out.mon.ref)
                out.mon = <path_mon_t*>sage_malloc(sizeof(path_mon_t))
                if out.mon == NULL:
                    raise MemoryError("Out of memory while allocating a path monomial")
                out.mon.ref = 1
                biseq_init_concat(out.mon.path, T1.mon.path, T2.mon.path)
            else:
                biseq_realloc_concat(out.mon.path, T1.mon.path, T2.mon.path)
        else:
            out.coef = NULL
    else:
        out = <path_term_t*>sage_malloc(sizeof(path_term_t))
        if out==NULL:
            raise MemoryError("Out of memory while allocating a path term")
        if new_coef:
            out.coef = <PyObject*>(new_coef)
            Py_INCREF(new_coef)
            out.mon = <path_mon_t*>sage_malloc(sizeof(path_mon_t))
            if out.mon == NULL:
                raise MemoryError("Out of memory while allocating a path monomial")
            out.mon.ref = 1
            biseq_init_concat(out.mon.path, T1.mon.path, T2.mon.path)
        else:
            out.coef = NULL
    out.mon.pos = new_pos
    out.mon.mid = new_mid
    out.mon.s_len = new_s_len
    out.nxt = NULL
    return out

########################################
##
## Basics for polynomials

cdef inline path_poly_t *poly_create() except NULL:
    cdef path_poly_t *out = <path_poly_t*>sage_malloc(sizeof(path_poly_t))
    if out==NULL:
        raise MemoryError("Out of memory while allocating a path polynomial")
    out.lead = NULL
    out.nterms = 0
    return out

cdef inline void poly_dealloc(path_poly_t *P):
    cdef path_term_t *T = P.lead
    while T!=NULL:
        T = term_free(T)

cdef inline void poly_free(path_poly_t *P):
    poly_dealloc(P)
    sage_free(P)

cdef inline bint poly_icopy(path_poly_t *out, path_poly_t *P) except -1:
    cdef path_term_t *T = P.lead
    out.nterms = P.nterms
    out.lead = term_copy_recursive(T)
    return True

cdef inline bint poly_icopy_neg(path_poly_t *out, path_poly_t *P) except -1:
    cdef path_term_t *T = P.lead
    out.nterms = P.nterms
    out.lead = term_neg_recursive(T)
    return True

cdef inline bint poly_icopy_scale(path_poly_t *out, path_poly_t *P, object coef) except -1:
    cdef path_term_t *T = P.lead
    cdef path_term_t *res = term_scale(T, coef)
    out.nterms = 0
    out.lead = NULL
    while res.coef == NULL:
        sage_free(res)
        T = T.nxt
        if T == NULL:
            return True
        res = term_scale(T, coef)
    out.lead = res
    preinc(out.nterms)
    T = T.nxt
    while T != NULL:
        res.nxt = term_scale(T, coef)
        if res.nxt.coef == NULL:
            sage_free(res.nxt)
        else:
            res = res.nxt
            preinc(out.nterms)
        T = T.nxt
    if res != NULL:
        res.nxt = NULL
    return True

cdef bint poly_is_sane(path_poly_t *P):
    cdef int count = 0
    cdef path_term_t *T = P.lead
    while T != NULL:
        preinc(count)
        T = T.nxt
    return count == P.nterms

cdef bint poly_is_ordered(path_poly_t *P, path_order_t cmp_terms):
    cdef path_term_t *T = P.lead
    while T != NULL:
        if T.nxt != NULL and cmp_terms(T.mon, T.nxt.mon)!=1:
            return False
        T = T.nxt
    return True

cdef inline list poly_pickle(path_poly_t *P):
    cdef list L = []
    cdef path_term_t *T = P.lead
    while T != NULL:
        L.append(term_pickle(T))
        T = T.nxt
    return L

cdef bint poly_inplace_unpickle(path_poly_t *P, list data) except -1:
    cdef tuple term_data
    cdef object coef
    cdef path_term_t *T
    if not data:
        P.nterms = 0
        P.lead = NULL
        return True
    P.nterms = len(data)
    coef, term_data = data.pop(0)
    P.lead = term_unpickle(coef, term_data)
    T = P.lead
    for coef, term_data in data:
        T.nxt = term_unpickle(coef, term_data)
        T = T.nxt
    T.nxt = NULL
    if not poly_is_sane(P):
        raise RuntimeError("Unpickling doesn't work")
    return True

############################################
##
## Polynomial arithmetics

# Comparison

cdef inline int poly_cmp(path_poly_t *P1, path_poly_t *P2, path_order_t cmp_terms):
    cdef path_term_t *T1 = P1.lead
    cdef path_term_t *T2 = P2.lead
    cdef int c
    while T1 != NULL and T2 != NULL:
        c = cmp_terms(T1.mon, T2.mon)
        if c != 0:
            return c
        c = cmp(<object>T1.coef, <object>T2.coef)
        if c != 0:
            return c
        T1 = T1.nxt
        T2 = T2.nxt
    if T1 == NULL:
        if T2 == NULL:
            return 0
        return -1
    return 1

cdef inline long poly_hash(path_poly_t *P):
    cdef path_term_t *T = P.lead
    cdef long out = 0
    while T != NULL:
        out = out<<7 | (out>>(sizeof(long)-7))
        out += term_hash(T)
        T = T.nxt
    return out

# Addition of a term
cdef inline void term_iadd(path_term_t *T1, path_term_t *T2):
    cdef object coef = <object>(T1.coef) + <object>(T2.coef)
    Py_DECREF(<object>(T1.coef))
    if coef:
        Py_INCREF(coef)
        T1.coef = <PyObject*>coef
    else:
        T1.coef = NULL

cdef bint poly_iadd_term_d(path_poly_t *P, path_term_t *T, path_order_t cmp_terms) except -1:
    if P.lead == NULL:
        preinc(P.nterms)
        T.nxt = NULL
        P.lead = T
        return True
    cdef path_term_t *tmp = P.lead
    cdef int c
    cdef object coef
    # First,
    c = cmp_terms(tmp.mon, T.mon)
    if c==-1:
        # The poly's lead term is smaller than T. Hence, we need to prepend
        # it.
        preinc(P.nterms)
        T.nxt = tmp
        P.lead = T
        return True
    elif c==0:
        term_iadd(tmp, T)
        term_free(T)
        if tmp.coef==NULL:
            predec(P.nterms)
            P.lead = term_free(tmp)
        elif <object>(tmp.coef)==0:
            raise RuntimeError("This should never happen")
        return True
    while True:
        # At this point, we have tmp>T.
        #
        # We need to append the term, or continue until we can
        # insert/append
        if tmp.nxt == NULL:
            preinc(P.nterms)
            T.nxt = NULL
            tmp.nxt = T
            return True
        c = cmp_terms(tmp.nxt.mon, T.mon)
        if c==-1:
            preinc(P.nterms)
            T.nxt = tmp.nxt
            tmp.nxt = T
            return True
        elif c==0:
            term_iadd(tmp.nxt, T)
            term_free(T)
            if tmp.nxt.coef==NULL:
                predec(P.nterms)
                tmp.nxt = term_free(tmp.nxt)
            elif <object>(tmp.coef)==0:
                raise RuntimeError("This should never happen")
            return True
        # otherwise, tmp is still larger than T. Hence, move to the next term
        # of P.
        tmp = tmp.nxt

# Addition of a poly - inplace for the first and destructive for the second argument
cdef inline bint poly_iadd_d(path_poly_t *P1, path_poly_t *P2, path_order_t cmp_terms) except -1:
    # Terms of P2 will be moved to P1, so that deallocation of P2 will not be
    # needed.
    if P1.lead == NULL:
        P1.lead = P2.lead
        P1.nterms = P2.nterms
        P2.nterms = 0
        P2.lead = NULL
        return 1
    if P2.lead == NULL:
        return 1
    cdef path_term_t *T1 = P1.lead
    cdef path_term_t *T2 = P2.lead
    cdef path_term_t *prev = NULL
    cdef int c
    cdef object new_coef
    while True:
        # Is one of the summands consumed already? Then we can use an easier
        # method.
        if T1 == NULL:
            if prev==NULL:
                P1.lead = T2
            else:
                prev.nxt = T2
            P1.nterms += P2.nterms
            P2.nterms = 0
            P2.lead = NULL
            return 1
        elif T2 == NULL:
            if P2.nterms != 0:
                print "term counting of second summand was wrong!",P2.nterms
            P2.lead = NULL
            return 1
        c = cmp_terms(T1.mon, T2.mon)
        if c == 0:
            # T1==T2 --- We need to add
            new_coef = <object>(T1.coef)+<object>(T2.coef)
            if new_coef:
                Py_INCREF(new_coef)
                Py_DECREF(<object>(T1.coef))
                T1.coef = <PyObject*>new_coef
                prev = T1
                T1 = T1.nxt
            else:
                T1 = term_free(T1)
                if prev==NULL:
                    P1.lead = T1
                else:
                    prev.nxt = T1
                predec(P1.nterms)
            predec(P2.nterms)
            T2 = term_free(T2)
        elif c == 1:
            # We move the prev/T1 through P1, until prev>T2>=T1. But now,
            # T1>T2. So, we move prev/T1 further down.
            prev = T1
            T1 = T1.nxt
        else:
            # prev > T2 > T1. Hence, we insert T2, without copying
            if prev==NULL:
                P1.lead = T2
                T2 = T2.nxt
                prev = P1.lead
            else:
                prev.nxt = T2
                T2 = T2.nxt
                prev = prev.nxt
            prev.nxt = T1
            preinc(P1.nterms)
            predec(P2.nterms)

# Addition of a poly, yielding a new poly and preserving the second argument
cdef path_poly_t *poly_add(path_poly_t *P1, path_poly_t *P2, path_order_t cmp_terms) except NULL:
    cdef path_poly_t *out = poly_create()
    cdef path_term_t *T1 = P1.lead
    cdef path_term_t *T2 = P2.lead
    cdef path_term_t *T = NULL
    cdef path_term_t *res
    cdef size_t count1, count2 # How many terms of P1/P2 have been considered?
    count1 = 0
    count2 = 0
    cdef object coef
    cdef int c
    while True:
        if T1 == NULL:
            out.nterms += (P2.nterms-count2)
            if T == NULL:
                if T2 == NULL:
                    out.lead = NULL
                else:
                    out.lead = term_copy_recursive(T2)
            else:
                if T2 == NULL:
                    T.nxt = NULL
                else:
                    T.nxt = term_copy_recursive(T2)
            return out
        if T2 == NULL:
            out.nterms += (P1.nterms-count1)
            if T == NULL:
                out.lead = term_copy_recursive(T1)
            else:
                T.nxt = term_copy_recursive(T1)
            return out

        c = cmp_terms(T1.mon,T2.mon)
        if c == 1:
            if T == NULL:
                out.lead = term_copy(T1)
                T = out.lead
            else:
                T.nxt = term_copy(T1)
                T = T.nxt
            T1 = T1.nxt
            preinc(count1)
            preinc(out.nterms)
        elif c == -1:
            if T == NULL:
                out.lead = term_copy(T2)
                T = out.lead
            else:
                T.nxt = term_copy(T2)
                T = T.nxt
            T2 = T2.nxt
            preinc(count2)
            preinc(out.nterms)
        else:
            coef = (<object>T1.coef)+(<object>T2.coef)
            if coef:
                preinc(out.nterms)
                if T == NULL:
                    out.lead = term_create(coef, T1.mon.path, T1.mon.pos, T1.mon.mid, T1.mon.s_len)
                    T = out.lead
                else:
                    T.nxt = term_create(coef, T1.mon.path, T1.mon.pos, T1.mon.mid, T1.mon.s_len)
                    T = T.nxt
            preinc(count1)
            preinc(count2)
            T1 = T1.nxt
            T2 = T2.nxt

# Subtraction of a poly, yielding a new poly and preserving the second argument
cdef path_poly_t *poly_sub(path_poly_t *P1, path_poly_t *P2, path_order_t cmp_terms) except NULL:
    cdef path_poly_t *out = poly_create()
    cdef path_term_t *T1 = P1.lead
    cdef path_term_t *T2 = P2.lead
    cdef path_term_t *T = NULL
    cdef path_term_t *res
    cdef size_t count1, count2 # How many terms of P1/P2 have been considered?
    count1 = 0
    count2 = 0
    cdef object coef
    cdef int c
    while True:
        if T1 == NULL:
            out.nterms += (P2.nterms-count2)
            if T == NULL:
                if T2 == NULL:
                    out.lead = NULL
                else:
                    out.lead = term_neg_recursive(T2)
            else:
                if T2 == NULL:
                    T.nxt = NULL
                else:
                    T.nxt = term_neg_recursive(T2)
            return out
        if T2 == NULL:
            out.nterms += (P1.nterms-count1)
            if T == NULL:
                out.lead = term_copy_recursive(T1)
            else:
                T.nxt = term_copy_recursive(T1)
            return out

        c = cmp_terms(T1.mon,T2.mon)
        if c == 1:
            if T == NULL:
                out.lead = term_copy(T1)
                T = out.lead
            else:
                T.nxt = term_copy(T1)
                T = T.nxt
            T1 = T1.nxt
            preinc(count1)
            preinc(out.nterms)
        elif c == -1:
            if T == NULL:
                out.lead = term_neg(T2)
                T = out.lead
            else:
                T.nxt = term_neg(T2)
                T = T.nxt
            T2 = T2.nxt
            preinc(count2)
            preinc(out.nterms)
        else:
            coef = (<object>T1.coef)-(<object>T2.coef)
            if coef:
                preinc(out.nterms)
                if T == NULL:
                    out.lead = term_create(coef, T1.mon.path, T1.mon.pos, T1.mon.mid, T1.mon.s_len)
                    T = out.lead
                else:
                    T.nxt = term_create(coef, T1.mon.path, T1.mon.pos, T1.mon.mid, T1.mon.s_len)
                    T = T.nxt
            preinc(count1)
            preinc(count2)
            T1 = T1.nxt
            T2 = T2.nxt

##
## In-place addition of a multiple of a polynomial

cdef int poly_iadd_lrmul(path_poly_t *P1, object coef, biseq_t L, path_poly_t *P2, biseq_t R, path_order_t cmp_terms, int mid, int pos, mp_size_t s_len) except -3:
    # Replace P1 by P1+coef*L*P2*R.
    #
    # Let m be a monomial of P2. If it is of module-type (i.e., m.mid!=-1),
    # then it is clear what L*m*R is (.mid and .s_len are obtained from
    # m.mid+L.length and m.s_len). Otherwise, however, it could be that we
    # want m2=L*m*R to denote a module-type monomial.
    #
    # In that case, "mid" encodes m2.mid.
    # - If 0<mid<=1+L.length, then m2.mid = mid-1.
    # - If -R.length-1<=mid<0, then m2.mid = R.length+mid+1-s_len
    # - If mid==0, the m2.mid=-1 and m2.s_len=0 and m2.pos=0.
    #
    # Also in that case, "pos" and "s_len" provide m2.pos and m2.s_len (unless mid==0).
    #
    # Return the largest degree of a monomial of the result, or -1 if it vanishes.
    #print ("start addmul poly")
    cdef int deg = -1
    cdef path_term_t *T1 = P1.lead
    if not coef or P2.lead==NULL:
        while T1 != NULL:
            deg = max(deg, term_total_degree(T1))
            T1 = T1.nxt
        return deg
    cdef path_mon_t *new_mon
    cdef object new_coef
    cdef path_term_t *prev = NULL
    cdef path_term_t *T2 = P2.lead
    cdef int c
    cdef size_t bla
    while T2!=NULL:
        new_coef = coef*<object>(T2.coef)
        if not new_coef:
            T2 = T2.nxt
            continue
        if T2.mon.mid!=-1:
            new_mon = path_mul_mon_mul_path(L, T2.mon, R)
            new_mon.mid = T2.mon.mid+L.length
            new_mon.pos = T2.mon.pos
            new_mon.s_len = T2.mon.s_len
        elif mid==0:
            new_mon = path_mul_mon_mul_path(L, T2.mon, R)
            new_mon.mid = -1
            new_mon.pos = 0
            new_mon.s_len = 0
        elif mid>0:
            #print ("new mid ==",mid-1)
            new_mon = path_mul_mon_mul_path(L, T2.mon, R)
            new_mon.mid = mid-1
            new_mon.pos = pos
            new_mon.s_len = s_len
        else:
            #print ("other new mid ==",L.length+T2.mon.path.length+R.length+mid+1)
            new_mon = path_mul_mon_mul_path(L, T2.mon, R)
            new_mon.mid = (L.length+T2.mon.path.length+R.length+mid+1)-s_len
            new_mon.pos = pos
            new_mon.s_len = s_len
        # Now new_term is L*T2*R
        # We go down in P1 until we may append, insert or add
        while T1!=NULL:
            c = cmp_terms(T1.mon, new_mon)
            if c!=1:
                break
            deg = max(deg, term_total_degree(T1))
            prev = T1
            T1 = prev.nxt
        if T1==NULL:
            # We need to append to P1
            if prev==NULL:
                P1.lead = term_create_keep_mon(new_coef, new_mon)
                prev = P1.lead
            else:
                prev.nxt = term_create_keep_mon(new_coef, new_mon)
                prev = prev.nxt
            deg = max(deg, term_total_degree(prev))
            prev.nxt = NULL
            preinc(P1.nterms)
        elif c==-1:
            # We need to insert between prev and T1
            if prev==NULL:
                P1.lead = term_create_keep_mon(new_coef, new_mon)
                prev = P1.lead
            else:
                prev.nxt = term_create_keep_mon(new_coef, new_mon)
                prev = prev.nxt
            deg = max(deg, term_total_degree(prev))
            prev.nxt = T1
            preinc(P1.nterms)
        else:
            # we add the coefficients and see what happens
            new_coef+=<object>(T1.coef)
            mon_free(new_mon)
            if new_coef:
                Py_INCREF(new_coef)
                Py_DECREF(<object>(T1.coef))
                T1.coef = <PyObject*>new_coef
                prev = T1
                T1 = T1.nxt
                deg = max(deg, term_total_degree(prev))
            else:
                predec(P1.nterms)
                T1 = term_free(T1)
                if prev==NULL:
                    P1.lead = T1
                else:
                    prev.nxt = T1
        T2 = T2.nxt
    while T1 != NULL:
        deg = max(deg, term_total_degree(T1))
        T1 = T1.nxt
    return deg

cdef path_term_t *poly_iadd_lmul(path_poly_t *P1, object coef, path_poly_t *P2, biseq_t R, path_order_t cmp_terms, int mid, int pos, mp_size_t s_len, path_term_t *P1start) except NULL:
    # Replace P1 by P1+coef*P2*R. Return a pointer to the first term of P1
    # that may be involved in a change when calling the function again with
    # P1, P2 and a "smaller" cofactor R. The return value should then be
    # provided as argument "P1start" of the next function call.
    #
    # We return P1start if P2.lead is NULL, otherwise if P1.lead becomes NULL during
    # addition, then we return P2.lead.
    #
    # Let m be a monomial of P2. If it is of module-type (i.e., m.mid!=-1),
    # then it is clear what m*R is (.mid and .s_len are obtained from
    # m.mid and m.s_len). Otherwise, however, it could be that we
    # want m2=m*R to denote a module-type monomial.
    #
    # In that case, "mid" encodes m2.mid.
    # - If 0<mid, then m2.mid = mid-1.
    # - If -R.length-1<=mid<0, then m2.mid = R.length+mid+1-s_len
    # - If mid==0, the m2.mid=-1 and m2.s_len=0 and m2.pos=0.
    #
    # Also in that case, "pos" and "s_len" provide m2.pos and m2.s_len (unless mid==0).
    if not coef or P2.lead==NULL:
        return P1start
    cdef path_mon_t *new_mon
    cdef object new_coef
    cdef path_term_t *prev = NULL
    cdef path_term_t *T1
    if P1start == NULL:
        T1 = P1.lead
    else:
        T1 = P1start
    cdef path_term_t *T2 = P2.lead
    cdef int c
    cdef path_term_t *out = P1start
    while T2!=NULL:
        new_coef = coef*<object>(T2.coef)
        if not new_coef:
            T2 = T2.nxt
            continue

        if T2.mon.mid!=-1:
            new_mon = mon_mul_path(T2.mon, R)
            new_mon.mid = T2.mon.mid
            new_mon.pos = T2.mon.pos
            new_mon.s_len = T2.mon.s_len
        elif mid==0:
            new_mon = mon_mul_path(T2.mon, R)
            new_mon.mid = -1
            new_mon.pos = 0
            new_mon.s_len = 0
        elif mid>0:
            #print ("new mid ==",mid-1)
            new_mon = mon_mul_path(T2.mon, R)
            new_mon.mid = mid-1
            new_mon.pos = pos
            new_mon.s_len = s_len
        else:
            #print ("other new mid ==",T2.mon.path.length+R.length+mid+1)
            new_mon = mon_mul_path(T2.mon, R)
            new_mon.mid = (T2.mon.path.length+R.length+mid+1)-s_len
            new_mon.pos = pos
            new_mon.s_len = s_len
        # Now new_term is T2*R
        # We go down in P1 until we may append, insert or add
        while T1!=NULL:
            c = cmp_terms(T1.mon, new_mon)
            if c!=1:
                break
            prev = T1
            T1 = prev.nxt
        if T1==NULL:
            # We need to append to P1
            if prev==NULL:
                P1.lead = term_create_keep_mon(new_coef, new_mon)
                prev = P1.lead
            else:
                prev.nxt = term_create_keep_mon(new_coef, new_mon)
                prev = prev.nxt
            if T2 == P2.lead:
                out = prev
            prev.nxt = NULL
            preinc(P1.nterms)
        elif c==-1:
            # We need to insert between prev and T1
            if prev==NULL:
                P1.lead = term_create_keep_mon(new_coef, new_mon)
                prev = P1.lead
            else:
                prev.nxt = term_create_keep_mon(new_coef, new_mon)
                prev = prev.nxt
            if T2 == P2.lead:
                out = prev
            prev.nxt = T1
            preinc(P1.nterms)
        else:
            # we add the coefficients and see what happens
            new_coef+=<object>(T1.coef)
            mon_free(new_mon)
            if new_coef:
                Py_INCREF(new_coef)
                Py_DECREF(<object>(T1.coef))
                T1.coef = <PyObject*>new_coef
                if T2 == P2.lead:
                    out = T1
                prev = T1
                T1 = T1.nxt
            else:
                predec(P1.nterms)
                T1 = term_free(T1)
                if prev==NULL:
                    P1.lead = T1
                else:
                    prev.nxt = T1
                if T2 == P2.lead:
                    out = prev
        T2 = T2.nxt
    #poly_is_sane(P1)
    if out == NULL:
        return P2.lead
    return out

########################################
##
## Basics homogeneous polynomials

cdef inline path_homog_poly_t *homog_poly_create(int start, int end) except NULL:
    cdef path_homog_poly_t *out = <path_homog_poly_t*>sage_malloc(sizeof(path_homog_poly_t))
    if out==NULL:
        raise MemoryError("Out of memory while allocating a path polynomial")
    out.poly = poly_create()
    out.start = start
    out.end = end
    out.nxt = NULL
    return out

cdef inline path_homog_poly_t *homog_poly_init_poly(int start, int end, path_poly_t *P) except NULL:
    # The polynomial P will not be copied!
    cdef path_homog_poly_t *out = <path_homog_poly_t*>sage_malloc(sizeof(path_homog_poly_t))
    if out==NULL:
        raise MemoryError("Out of memory while allocating a path polynomial")
    out.poly = P
    out.start = start
    out.end = end
    out.nxt = NULL
    return out

cdef path_homog_poly_t *homog_poly_init_list(int start, int end, list L, path_order_t cmp_terms, int mid, unsigned int pos, mp_size_t s_len) except NULL:
    # mid==-2 means that we have a left module element. I.e., the new mid is pathlength-s_len
    cdef path_homog_poly_t * out = homog_poly_create(start, end)
    cdef QuiverPath P
    for P,coef in L:
        poly_iadd_term_d(out.poly, term_create(coef, P._path, pos, mid if mid!=-2 else P._path.length-s_len, s_len), cmp_terms)
    return out

cdef inline void homog_poly_free(path_homog_poly_t *P):
    cdef path_homog_poly_t *nxt
    while P!=NULL:
        nxt = P.nxt
        poly_free(P.poly)
        sage_free(P)
        P = nxt

cdef inline path_homog_poly_t *homog_poly_copy(path_homog_poly_t *H) except NULL:
    cdef path_homog_poly_t *out
    cdef path_homog_poly_t *tmp
    if H == NULL:
        raise ValueError("The polynomial to be copied is the NULL pointer")
    out = homog_poly_create(H.start, H.end)
    poly_icopy(out.poly, H.poly)
    tmp = out
    H = H.nxt
    while H != NULL:
        tmp.nxt = homog_poly_create(H.start, H.end)
        tmp = tmp.nxt
        poly_icopy(tmp.poly, H.poly)
        H = H.nxt
    return out

cdef inline list homog_poly_pickle(path_homog_poly_t *H):
    cdef list L = []
    while H != NULL:
        L.append((H.start, H.end, poly_pickle(H.poly)))
        H = H.nxt
    return L

cdef path_homog_poly_t *homog_poly_unpickle(list data) except NULL:
    #ASSUMPTION: data is not empty
    cdef int start, end
    cdef list poly_data
    cdef path_homog_poly_t *out
    start, end, poly_data = data.pop(0)
    out = homog_poly_create(start, end)
    poly_inplace_unpickle(out.poly, poly_data)
    cdef path_homog_poly_t *tmp = out
    for start, end, poly_data in data:
        tmp.nxt = homog_poly_create(start, end)
        tmp = tmp.nxt
        poly_inplace_unpickle(tmp.poly, poly_data)
    return out

cdef inline path_homog_poly_t *homog_poly_neg(path_homog_poly_t *H) except NULL:
    cdef path_homog_poly_t *out
    cdef path_homog_poly_t *tmp
    if H == NULL:
        raise ValueError("The polynomial to be copied is the NULL pointer")
    out = homog_poly_create(H.start, H.end)
    poly_icopy_neg(out.poly, H.poly)
    tmp = out
    H = H.nxt
    while H != NULL:
        tmp.nxt = homog_poly_create(H.start, H.end)
        tmp = tmp.nxt
        poly_icopy_neg(tmp.poly, H.poly)
        H = H.nxt
    return out

cdef inline path_homog_poly_t *homog_poly_scale(path_homog_poly_t *H, object coef) except NULL:
    # The first component may be zero, all other zero components are removed.
    cdef path_homog_poly_t *out
    cdef path_homog_poly_t *tmp
    if H == NULL:
        raise ValueError("The polynomial to be copied is the NULL pointer")
    out = homog_poly_create(H.start, H.end)
    poly_icopy_scale(out.poly, H.poly, coef)
    tmp = out
    H = H.nxt
    while H != NULL:
        tmp.nxt = homog_poly_create(H.start, H.end)
        poly_icopy_scale(tmp.nxt.poly, H.poly, coef)
        if tmp.nxt.poly.nterms == 0:
            homog_poly_free(tmp.nxt)
            tmp.nxt = NULL
        else:
            tmp = tmp.nxt
        H = H.nxt
    return out

cdef inline path_homog_poly_t *homog_poly_get_predecessor_of_component(path_homog_poly_t *H, int s, int e):
    # Search through H.nxt.nxt... and return the pointer C to a component of H
    # such that either C.nxt.start==s and C.nxt.end==e, or the component for
    # (s,e) should be inserted between C and C.nxt. Return NULL if H==NULL or
    # (s,e) should be inserted in front of H.
    if H == NULL:
        return NULL
    if H.start > s:
        return NULL
    elif H.start == s and H.end >= e:
        return NULL
    while True:
        if H.nxt == NULL:
            return H
        if H.nxt.start == s:
            if H.nxt.end >= e:
                return H
        elif H.nxt.start > s:
            return H
        H = H.nxt
