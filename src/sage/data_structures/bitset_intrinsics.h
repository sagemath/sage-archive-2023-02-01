#ifndef SAGE_BITSET_INTRINSICS
#define SAGE_BITSET_INTRINSICS

#include "gmp.h"

#if __POPCNT__ || __BMI__
    #include <immintrin.h>
#endif

#if __AVX2__
    #include <immintrin.h>
#elif __AVX__
    #include <immintrin.h>
    #include <emmintrin.h>
#elif __SSE4_1__
    #include <emmintrin.h>
    #include <smmintrin.h>
#endif


// This file contains functions of ``bitset_base.pxd``
// that can be optimized using intrinsics.

/*
#############################################################################
# Bitset Comparison
#############################################################################
*/

const mp_bitcnt_t LIMB_SIZE = sizeof(mp_limb_t);

#if __AVX__
    const mp_bitcnt_t ALIGNMENT = 32;
#elif __SSE4_1__
    const mp_bitcnt_t ALIGNMENT = 16;
#else
    const mp_bitcnt_t ALIGNMENT = 8;
#endif

// We assume that ``GMP_LIMB_BITS`` is a divisor of 64
#define LIMBS_PER_64 64/GMP_LIMB_BITS

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
    mp_bitcnt_t i;
#if __AVX__
    for(i = 0; i + 3*LIMBS_PER_64 < limbs; i +=4*LIMBS_PER_64){
        __m256i A = _mm256_loadu_si256((const __m256i*)&a[i]);
        __m256i B = _mm256_loadu_si256((const __m256i*)&b[i]);
        if (!_mm256_testc_si256(A, B) || !_mm256_testc_si256(B, A))
            return 0;
    }
    for(i; i < limbs; i++){
        if (a[i] != b[i])
            return 0;
    }
    return 1;
#else
    return mpn_cmp(a, b, limbs) == 0;
#endif
}

static inline int _bitset_issubset(mp_limb_t* a, mp_limb_t* b, mp_bitcnt_t limbs){
    /*
    Test whether a is a subset of b (i.e., every element in a is also
    in b).
    */
    mp_bitcnt_t i = 0;
#if __AVX__
    for(i; i + 3*LIMBS_PER_64 < limbs; i +=4*LIMBS_PER_64){
        __m256i A = _mm256_loadu_si256((const __m256i*)&a[i]);
        __m256i B = _mm256_loadu_si256((const __m256i*)&b[i]);
        if (!_mm256_testc_si256(B, A)) // Need to be opposite order!
            return 0;
    }
#elif __SSE4_1__
    for(i; i + LIMBS_PER_64 < limbs; i += 2*LIMBS_PER_64){
        __m128i A = _mm_loadu_si128((const __m128i*)&a[i]);
        __m128i B = _mm_loadu_si128((const __m128i*)&b[i]);
        if (!_mm_testc_si128(B, A)) // Need to be opposite order!
            return 0;
    }
#endif
    for(i; i < limbs; i++){
        if ((a[i] & ~b[i]) != 0)
            return 0;
    }
    return 1;
}

static inline int _bitset_are_disjoint(mp_limb_t* a, mp_limb_t* b, mp_bitcnt_t limbs){
    /*
    Tests whether ``a`` and ``b`` have an empty intersection.
    */
    mp_bitcnt_t i = 0;
#if __AVX__
    for(i; i + 3*LIMBS_PER_64 < limbs; i +=4*LIMBS_PER_64){
        __m256i A = _mm256_loadu_si256((const __m256i*)&a[i]);
        __m256i B = _mm256_loadu_si256((const __m256i*)&b[i]);
        if (!_mm256_testz_si256(A, B))
            return 0;
    }
#elif __SSE4_1__
    for(i; i + LIMBS_PER_64 < limbs; i += 2*LIMBS_PER_64){
        __m128i A = _mm_loadu_si128((const __m128i*)&a[i]);
        __m128i B = _mm_loadu_si128((const __m128i*)&b[i]);
        if (!_mm_testz_si128(A, B))
            return 0;
    }
#endif
    for(i; i < limbs; i++){
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
    mp_bitcnt_t i;
#if __AVX2__
    for(i = 0; i + 3*LIMBS_PER_64 < limbs; i +=4*LIMBS_PER_64){
        __m256i A = _mm256_loadu_si256((const __m256i*)&a[i]);
        __m256i B = _mm256_loadu_si256((const __m256i*)&b[i]);
        __m256i D = _mm256_and_si256(A, B);
        _mm256_storeu_si256((__m256i*)&dst[i], D);
    }
    for(i; i < limbs; i++){
        dst[i] = a[i] & b[i];
    }
#elif __SSE4_1__
    for(i = 0; i + LIMBS_PER_64 < limbs; i += 2*LIMBS_PER_64){
        __m128i A = _mm_loadu_si128((const __m128i*)&a[i]);
        __m128i B = _mm_loadu_si128((const __m128i*)&b[i]);
        __m128i D = _mm_and_si128(A, B);
        _mm_storeu_si128((__m128i*)&dst[i], D);
    }
    for(i; i < limbs; i++){
        dst[i] = a[i] & b[i];
    }
#else
    mpn_and_n(dst, a, b, limbs);
#endif
}

static inline void _bitset_union(mp_limb_t* dst, mp_limb_t* a, mp_limb_t* b, mp_bitcnt_t limbs){
    /*
    Set dst to the union of a and b, overwriting dst.
    */
    mp_bitcnt_t i;
#if __AVX2__
    for(i = 0; i + 3*LIMBS_PER_64 < limbs; i +=4*LIMBS_PER_64){
        __m256i A = _mm256_loadu_si256((const __m256i*)&a[i]);
        __m256i B = _mm256_loadu_si256((const __m256i*)&b[i]);
        __m256i D = _mm256_or_si256(A, B);
        _mm256_storeu_si256((__m256i*)&dst[i], D);
    }
    for(i; i < limbs; i++){
        dst[i] = a[i] | b[i];
    }
#elif __SSE4_1__
    for(i = 0; i + LIMBS_PER_64 < limbs; i += 2*LIMBS_PER_64){
        __m128i A = _mm_loadu_si128((const __m128i*)&a[i]);
        __m128i B = _mm_loadu_si128((const __m128i*)&b[i]);
        __m128i D = _mm_or_si128(A, B);
        _mm_storeu_si128((__m128i*)&dst[i], D);
    }
    for(i; i < limbs; i++){
        dst[i] = a[i] | b[i];
    }
#else
    mpn_ior_n(dst, a, b, limbs);
#endif
}

static inline void _bitset_difference(mp_limb_t* dst, mp_limb_t* a, mp_limb_t* b, mp_bitcnt_t limbs){
    /*
    Set dst to the difference of a and b (i.e., things in a that are not
    in b), overwriting dst.
    */
    mp_bitcnt_t i;
#if __AVX2__
    for(i = 0; i + 3*LIMBS_PER_64 < limbs; i +=4*LIMBS_PER_64){
        __m256i A = _mm256_loadu_si256((const __m256i*)&a[i]);
        __m256i B = _mm256_loadu_si256((const __m256i*)&b[i]);
        __m256i D = _mm256_andnot_si256(B, A);  // Need to be opposite order.
        _mm256_storeu_si256((__m256i*)&dst[i], D);
    }
    for(i; i < limbs; i++){
        dst[i] = a[i] & ~b[i];
    }
#elif __SSE4_1__
    for(i = 0; i + LIMBS_PER_64 < limbs; i += 2*LIMBS_PER_64){
        __m128i A = _mm_loadu_si128((const __m128i*)&a[i]);
        __m128i B = _mm_loadu_si128((const __m128i*)&b[i]);
        __m128i D = _mm_andnot_si128(B, A);  // Need to be opposite order.
        _mm_storeu_si128((__m128i*)&dst[i], D);
    }
    for(i; i < limbs; i++){
        dst[i] = a[i] & ~b[i];
    }
#else
    mpn_andn_n(dst, a, b, limbs);
#endif
}

static inline void _bitset_symmetric_difference(mp_limb_t* dst, mp_limb_t* a, mp_limb_t* b, mp_bitcnt_t limbs){
    /*
    Set dst to the symmetric difference of a and b, overwriting dst.
    */
    mp_bitcnt_t i;
#if __AVX2__
    for(i = 0; i + 3*LIMBS_PER_64 < limbs; i +=4*LIMBS_PER_64){
        __m256i A = _mm256_loadu_si256((const __m256i*)&a[i]);
        __m256i B = _mm256_loadu_si256((const __m256i*)&b[i]);
        __m256i D = _mm256_xor_si256(A, B);
        _mm256_storeu_si256((__m256i*)&dst[i], D);
    }
    for(i; i < limbs; i++){
        dst[i] = a[i] ^ b[i];
    }
#elif __SSE4_1__
    for(i = 0; i + LIMBS_PER_64 < limbs; i += 2*LIMBS_PER_64){
        __m128i A = _mm_loadu_si128((const __m128i*)&a[i]);
        __m128i B = _mm_loadu_si128((const __m128i*)&b[i]);
        __m128i D = _mm_xor_si128(A, B);
        _mm_storeu_si128((__m128i*)&dst[i], D);
    }
    for(i; i < limbs; i++){
        dst[i] = a[i] ^ b[i];
    }
#else
    mpn_xor_n(dst, a, b, limbs);
#endif
}

#endif
