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


// no intrinsics for now,
// yet the following 3 functions are made,
// such that improving them by intrinsics, shall be easy

const size_t chunksize = 64;
inline int is_subset(uint64_t *A, uint64_t *B, size_t face_length){
    // A & ~B == 0
    // returns 1 if A is a subset of B, otherwise returns 0
    // this is done by checking if there is an element in A,
    // which is not in B
    // `face_length` is the length of A and B in terms of uint64_t
    for (size_t i = 0; i < face_length; i++){
        if (A[i] & ~B[i]){
            return 0;
        }
    }
    return 1;
}

inline void intersection(uint64_t *A, uint64_t *B, uint64_t *C, \
                         size_t face_length){
    // C = A & B
    // will set C to be the intersection of A and B
    // `face_length` is the length of A, B and C in terms of uint64_t
    for (size_t i = 0; i < face_length; i++){
        C[i] = A[i] & B[i];
    }
}

inline size_t CountFaceBits(uint64_t* A, size_t face_length) {
    // counts the number of vertices in a face by counting bits set to one
    // `face_length` is the length of A in terms of uint64_t
    unsigned int count = 0;
    for (size_t i=0; i<face_length; i++){
        uint64_t a = A[i];
        while (a){
            count += a & 1;
            a >>= 1;
        }
    }
    return count;
}

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
    int *ommitfacearray = (int *) calloc(lenfaces, sizeof(int));
    // this array has an entry for each entry in nextfaces
    // iff ommitfacearray[i], then we want to ommit the corresponding face
    // newfaces2 will then just contain pointers to all other faces
    for (size_t j = 0; j < lenfaces - 1; j++){
        intersection(faces[j], faces[lenfaces - 1], nextfaces[j], face_length);
    }
    // we have create all possible intersection with the i_th-face,
    // but some of them might not be of exactly one dimension less
    for (size_t j = 0; j < lenfaces-1; j++){
        for(size_t k = 0; k < j; k++){
            // testing if nextfaces[j] is contained in different nextface
            if(is_subset(nextfaces[j], nextfaces[k],face_length)){
                ommitfacearray[j] = 1;
                // nextfaces[j] is a subset of nextfaces[k] -> do not add it
                break;
                }
            }
        if (ommitfacearray[j]) {
            continue;
        }

        for(size_t k = j+1; k < lenfaces-1; k++){
            // testing if nextfaces[j] is contained in a different nextface
            if(is_subset(nextfaces[j],nextfaces[k], face_length)){
            ommitfacearray[j] = 1;
            // nextfaces[j] is a subset of nextfaces[k] -> do not add it
            break;
            }
        }
        if (ommitfacearray[j]) {
            continue;
        }

        for (size_t k = 0; k < nr_forbidden; k++){
            // we do not want to double count any faces
            if(is_subset(nextfaces[j],forbidden[k], face_length)){
                ommitfacearray[j] = 1;
                // nextfaces[j] is a subset of forbidden[k]
                // as we have visited all faces of forbidden[k] already, we must
                // have visitied nextfaces[j] before
                break;
            }
        }
    }

    size_t newfacescounter = 0;  // length of newfaces2
    for (size_t j = 0; j < lenfaces -1; j++){
        // let `newfaces2` point to the newfaces we want to consider
        if (ommitfacearray[j]) {
            continue;
        }
        nextfaces2[newfacescounter] = nextfaces[j];
        newfacescounter++;
    }
    free(ommitfacearray);
    return newfacescounter;
}


size_t facet_repr_from_bitrep(uint64_t *face, uint64_t **facets, \
                              size_t *output, size_t nr_facets, \
                              size_t face_length){
    // Writes the facet_repr of the current face in output.
    // Returns the length of the representation.
    size_t counter = 0;
    for (size_t i = 0; i < nr_facets; i++){
        if (is_subset(face, facets[i], face_length)){
            output[counter] = i;
            counter++;
        }
    }
    return counter;
}
