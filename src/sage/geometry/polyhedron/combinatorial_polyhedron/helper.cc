/*
#*****************************************************************************
#       Copyright (C) 2019 Jonathan Kliem <jonathan.kliem@fu-berlin.de>
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 2 of the License, or
# (at your option) any later version.
#                  http://www.gnu.org/licenses/
#*****************************************************************************
*/

// Copyright: see base.pyx

#include <cstdint>
#include <cstdio>


// ************** intrinsics prepared, but not enabled for now ***************
/*
#if __AVX__
    #include <immintrin.h>
#elif __SSE4_1__
    #include <emmintrin.h>
    #include <smmintrin.h>
#endif

#if __POPCNT__
    #include <immintrin.h>
#endif
*/

// as of now, 512bit does not have something like _mm256_testc_si256,
// which is the bottle neck of this function,
// so it does not make sense to implement it

// inline int is_subset(uint64_t *A, uint64_t *B, size_t face_length)
// the bottlen-neck is checking for subsets, which requires something as
// _mm256_testc_si256, trying to determine, what is the best way of doing it:
#if 0//__AVX__
    // 256-bit commands, those operations are equivalent to the operations
    // defined in `#else`
    // intrics defined in immintrin.h
    const size_t chunksize = 256;
    inline int is_subset(uint64_t *A, uint64_t *B, size_t face_length){
        // A & ~B == 0
        // returns 1 if A is a subset of B, otherwise returns 0
        // this is done by checking if there is an element in A,
        // which is not in B
        // `face_length` is the length of A and B in terms of uint64_t
        // note that A,B need to be 32-byte-aligned
        size_t i;
        for (i = 0; i < face_length; i += 4){
            __m256i a = _mm256_load_si256((const __m256i*)&A[i]);
            __m256i b = _mm256_load_si256((const __m256i*)&B[i]);
            if (!_mm256_testc_si256(b, a)){ //need to be opposite order !!
                return 0;
            }
        }
        return 1;
    }

#elif 0//__SSE4_1__
    // 128-bit commands, those operations are equivalent to the operations
    // defined in `#else`
    // intrics defined in smmintrin.h and emmintrin.h
    const size_t chunksize = 128;
    inline int is_subset(uint64_t *A, uint64_t *B, size_t face_length){
        // A & ~B == 0
        // returns 1 if A is a subset of B, otherwise returns 0
        // this is done by checking if there is an element in A,
        // which is not in B
        // `face_length` is the length of A and B in terms of uint64_t
        // note that A,B need to be 16-byte-aligned
        size_t i;
        for (i = 0; i < face_length; i += 2){
            __m128i a = _mm_load_si128((const __m128i*)&A[i]);
            __m128i b = _mm_load_si128((const __m128i*)&B[i]);
            if (!_mm_testc_si128(b, a)){ //need to be opposite order !!
                return 0;
            }
        }
        return 1;
    }

#else
    // no intrinsics
    const size_t chunksize = 64;
    inline int is_subset(uint64_t *A, uint64_t *B, size_t face_length){
        // A & ~B == 0
        // returns 1 if A is a subset of B, otherwise returns 0
        // this is done by checking if there is an element in A,
        // which is not in B
        // `face_length` is the length of A and B in terms of uint64_t
        size_t i;
        for (i = 0; i < face_length; i++){
            if (A[i] & ~B[i]){
                return 0;
            }
        }
        return 1;
    }

#endif

// inline void intersection(uint64_t *A, uint64_t *B, uint64_t *C,
//                          size_t face_length)
// now determining, how to do insersection
#if 0//__AVX2__
    // 256-bit commands, those operations are equivalent to the operations
    // defined in `#else`
    // intrics defined in immintrin.h
    inline void intersection(uint64_t *A, uint64_t *B, uint64_t *C, \
                             size_t face_length){
        // C = A & B
        // will set C to be the intersection of A and B
        // `face_length` is the length of A, B and C in terms of uint64_t
        // note that A,B,C need to be 32-byte-aligned
        size_t i;
        for (i = 0; i < face_length; i += 4){
            __m256i a = _mm256_load_si256((const __m256i*)&A[i]);
            __m256i b = _mm256_load_si256((const __m256i*)&B[i]);
            __m256i c = _mm256_and_si256(a, b);
            _mm256_store_si256((__m256i*)&C[i],c);
        }
    }

#elif 0//__SSE4_1__
    // actually SSE2 would be fine, but we don't want to force greater chunks,
    // because of intersection, which is not the bottleneck
    // 128-bit commands, those operations are equivalent to the operations
    // defined in `#else`
    // intrinsics defined in emmintrin.h
    inline void intersection(uint64_t *A, uint64_t *B, uint64_t *C, \
                             size_t face_length){
        // C = A & B
        // will set C to be the intersection of A and B
        // `face_length` is the length of A, B and C in terms of uint64_t
        // note that A,B,C need to be 16-byte-aligned
        size_t i;
        for (i = 0; i < face_length; i += 2){
            __m128i a = _mm_load_si128((const __m128i*)&A[i]);
            __m128i b = _mm_load_si128((const __m128i*)&B[i]);
            __m128i c = _mm_and_si128(a, b);
            _mm_store_si128((__m128i*)&C[i],c);
        }
    }

#else
    // commands, without intrinsics
    inline void intersection(uint64_t *A, uint64_t *B, uint64_t *C, \
                             size_t face_length){
        // C = A & B
        // will set C to be the intersection of A and B
        // `face_length` is the length of A, B and C in terms of uint64_t
        size_t i;
        for (i = 0; i < face_length; i++){
            C[i] = A[i] & B[i];
        }
    }

#endif

// inline size_t CountFaceBits(uint64_t* A, size_t face_length)
// determine the best way to count the set bits in uint64_t *
#if 0//(__POPCNT__) && (INTPTR_MAX == INT64_MAX) // 64-bit and popcnt
    inline size_t CountFaceBits(uint64_t* A, size_t face_length) {
        // counts the number of vertices in a face by counting bits set to one
        // `face_length` is the length of A in terms of uint64_t
        size_t i;
        unsigned int count = 0;
        for (i=0; i<face_length; i++){
            count += (size_t) _mm_popcnt_u64(A[i]);
        }
        return count;
    }

#else // popcount without intrinsics
    inline size_t CountFaceBits(uint64_t* A, size_t face_length) {
        // counts the number of vertices in a face by counting bits set to one
        // `face_length` is the length of A in terms of uint64_t
        size_t i;
        unsigned int count = 0;
        for (i=0; i<face_length; i++){
            uint64_t a = A[i];
            while (a){
                count += a & 1;
                a >>= 1;
            }
        }
        return count;
    }

#endif



size_t get_next_level(\
        uint64_t **faces, size_t lenfaces, uint64_t **nextfaces, \
        uint64_t **nextfaces2, uint64_t **forbidden, \
        size_t nr_forbidden, size_t face_length){
    // intersects the first `lenfaces - 1` faces of `faces`
    // with faces[lenfaces-1]`
    // determines which ones are exactly of one dimension less
    // by considering containment
    // newfaces2 will point at those of exactly one dimension less
    // which are not contained in any of the faces in `forbidden`
    // returns the number of those faces
    const size_t constlenfaces = lenfaces;
    int addfacearray[constlenfaces - 1] = { };
    size_t j,k, addthisface;
    size_t newfacescounter = 0;
    for (j = 0; j < lenfaces - 1; j++){
        intersection(faces[j], faces[lenfaces - 1], nextfaces[j], face_length);
        addfacearray[j] = 1;
    }
    // we have create all possible intersection with the i_th-face,
    // but some of them might not be of exactly one dimension less
    for (j = 0; j < lenfaces-1; j++){
        for(k = 0; k < j; k++){
            // testing if nextfaces[j] is contained in different nextface
            if(is_subset(nextfaces[j], nextfaces[k],face_length)){
                addfacearray[j] = 0;
                // nextfaces[j] is a subset of nextfaces[k] -> do not add it
                break;
                }
            }
        if (!addfacearray[j]) {
            continue;
        }

        for(k = j+1; k < lenfaces-1; k++){
            // testing if nextfaces[j] is contained in a different nextface
            if(is_subset(nextfaces[j],nextfaces[k], face_length)){
            addfacearray[j] = 0;
            // nextfaces[j] is a subset of nextfaces[k] -> do not add it
            break;
            }
        }
        if (!addfacearray[j]) {
            continue;
        }

        for (k = 0; k < nr_forbidden; k++){
            // we do not want to double count any faces
            if(is_subset(nextfaces[j],forbidden[k], face_length)){
                addfacearray[j] = 0;
                // nextfaces[j] is a subset of forbidden[k]
                // as we have visited all faces of forbidden[k] already, we must
                // have visitied nextfaces[j] before
                break;
            }
        }
    }

    for (j = 0; j < lenfaces -1; j++){
        // let `newfaces2` point to the newfaces we want to consider
        if (!addfacearray[j]) {
            continue;
        }
        nextfaces2[newfacescounter] = nextfaces[j];
        newfacescounter++;
    }
    return newfacescounter;
}


size_t facet_repr_from_bitrep(uint64_t *face, uint64_t **facets, \
                              size_t *output, size_t nr_facets, \
                              size_t face_length){
    // Writes the facet_repr of the current face in output.
    // Returns the length of the representation.
    size_t counter = 0;
    size_t i;
    for (i = 0; i < nr_facets; i++){
        if (is_subset(face, facets[i], face_length)){
            output[counter] = i;
            counter++;
        }
    }
    return counter;
}
