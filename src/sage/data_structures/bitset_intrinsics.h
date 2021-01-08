#include "gmp.h"

#if __POPCNT__ || __BMI__
    #include <immintrin.h>
#endif


// This file contains functions of ``bitset_base.pxd``
// that can be optimized using intrinsics.

/*
#############################################################################
# Bitset Comparison
#############################################################################
*/

const mp_bitcnt_t LIMB_SIZE = sizeof(mp_limb_t);
const mp_bitcnt_t ALIGNMENT = sizeof(mp_limb_t);

static inline int _bitset_isempty(mp_limb_t* bits, mp_bitcnt_t limbs){
    /*
    Test whether bits is empty.  Return True (i.e., 1) if the set is
    empty, False (i.e., 0) otherwise.
    */
    // First check lowest limb
    if (bits[0])
        return 0;
    if (limbs == 1)
        return 1;
    // Compare bits to itself shifted by 1 limb. If these compare equal,
    // all limbs must be 0.
    return mpn_cmp(bits+1, bits, limbs-1) == 0;
}

static inline int _bitset_eq(mp_limb_t* a, mp_limb_t* b, mp_bitcnt_t limbs){
    /*
    Compare bitset a and b.  Return True (i.e., 1) if the sets are
    equal, and False (i.e., 0) otherwise.
    */
    return mpn_cmp(a, b, limbs) == 0;
}

static inline int _bitset_issubset(mp_limb_t* a, mp_limb_t* b, mp_bitcnt_t limbs){
    /*
    Test whether a is a subset of b (i.e., every element in a is also
    in b).
    */
    mp_bitcnt_t i;
    for(i = 0; i < limbs; i++){
        if ((a[i] & ~b[i]) != 0)
            return 0;
    }
    return 1;
}

static inline int _bitset_are_disjoint(mp_limb_t* a, mp_limb_t* b, mp_bitcnt_t limbs){
    /*
    Tests whether ``a`` and ``b`` have an empty intersection.
    */
    mp_bitcnt_t i;
    for(i = 0; i < limbs; i++){
        if (a[i] & b[i])
            return 0;
    }
    return 1;
}

/*
#############################################################################
# Bitset Searching
#############################################################################
*/

static inline long _bitset_first_in_limb_nonzero(mp_limb_t limb){
    /*
    Given a non-zero limb of a bitset, return the index of the first
    nonzero bit.
    */
#if (__BMI__) && (GMP_LIMB_BITS == 64) && (INTPTR_MAX == INT64_MAX)
    // If available we use intrinsics trailing zero count.
    // Also we check that ``mp_bitcnt_t`` is in fact 8 bytes long
    // and the architecture is 64-bit.
    return _tzcnt_u64(limb);
#else
    return mpn_scan1(&limb, 0);
#endif
}

static inline long _bitset_first_in_limb(mp_limb_t limb){
    /*
    Given a limb of a bitset, return the index of the first nonzero
    bit. If there are no bits set in the limb, return -1.
    */
    if (limb == 0)
        return -1;
    return _bitset_first_in_limb_nonzero(limb);
}

static inline long _bitset_len(mp_limb_t* bits, mp_bitcnt_t limbs){
    /*
    Calculate the number of items in the set (i.e., the number of nonzero bits).
    */
#if (__POPCNT__) && (GMP_LIMB_BITS == 64) && (INTPTR_MAX == INT64_MAX)
    // If available we use intrinsics popcount.
    // Also we check that ``mp_bitcnt_t`` is in fact 8 bytes long
    // and the architecture is 64-bit.
    mp_bitcnt_t i;
    uint64_t count = 0;
    for (i=0; i<limbs; i++){
        count += _mm_popcnt_u64(bits[i]);
    }
    return count;
#else
    return mpn_popcount(bits, limbs);
#endif
}

/*
#############################################################################
# Bitset Arithmetic
#############################################################################
*/

static inline void _bitset_intersection(mp_limb_t* dst, mp_limb_t* a, mp_limb_t* b, mp_bitcnt_t limbs){
    /*
    Set dst to the intersection of a and b, overwriting dst.
    */
    mpn_and_n(dst, a, b, limbs);
}

static inline void _bitset_union(mp_limb_t* dst, mp_limb_t* a, mp_limb_t* b, mp_bitcnt_t limbs){
    /*
    Set dst to the union of a and b, overwriting dst.
    */
    mpn_ior_n(dst, a, b, limbs);
}

static inline void _bitset_difference(mp_limb_t* dst, mp_limb_t* a, mp_limb_t* b, mp_bitcnt_t limbs){
    /*
    Set dst to the difference of a and b (i.e., things in a that are not
    in b), overwriting dst.
    */
    mpn_andn_n(dst, a, b, limbs);
}

static inline void _bitset_symmetric_difference(mp_limb_t* dst, mp_limb_t* a, mp_limb_t* b, mp_bitcnt_t limbs){
    /*
    Set dst to the symmetric difference of a and b, overwriting dst.
    */
    mpn_xor_n(dst, a, b, limbs);
}
