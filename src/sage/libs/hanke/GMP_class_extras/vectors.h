#ifndef VECTORS_H
#define VECTORS_H

using namespace std;
#include <gmp.h>
#include <gmpxx.h>
#include <valarray>


// Notes:
//    - Sets will mean vectors assumed to be in increasing order, with
//        no duplicates (i.e. strictly increasing).
//    - Routines ending in the word "Ordered" assume that the input is an
//        ordered (i.e. strictly increasing) vector, also referred to as a set.
//        The entries should also be positive integers, but this may not be
//        enforced (I need to look).
//
//    - Types of returns:
//      -----------------
//        Tuples:
//            VectorAppend - append as tuples
//            MakeVector - makes a tuple
//            VectorRemove - makes a tuple
//        Sets:
//            OrderVector - makes a set from a tuple
//            VectorUnion - makes a set
//            VectorIntersection - makes a set
//            VectorComplement - makes a set
//        Boolean (long):
//            IsDisjointOrdered - 0 or 1
//            CheckVectorIntersection - 0 or 1
//            CheckVectorOrdered - 0 or 1




//////////////////////////////////////////////////////////////////////
// Checks if two ordered (strictly increasing) vectors are disjoint //
//////////////////////////////////////////////////////////////////////

bool IsDisjointOrdered(const valarray<int> & v, const valarray<int> & w);


/////////////////////////////////////////////////
// Gives a vector by an arithmetic progression //
/////////////////////////////////////////////////

valarray<int> MakeVector(int length, int start, int increment);


//////////////////////////////////////////////////////////////////////////////////
// Reorders the vector v in strictly increasing order (so no repetitions occur) //
//////////////////////////////////////////////////////////////////////////////////

valarray<int> OrderVector(const valarray<int> & v);


//////////////////////////////////////////////
// Creates the union (v + w) of two vectors //
//////////////////////////////////////////////

valarray<int> VectorAppend(const valarray<int> & v, const valarray<int> & w);


//////////////////////////////////////////////////////////////////
// Finds the ordered (strictly increasing) union of two vectors //
//////////////////////////////////////////////////////////////////
valarray<int> VectorUnion(const valarray<int> & v, const valarray<int> & w);


/////////////////////////////////////////////////////////////////////////////////////
// Finds their ordered union, assuming they are both ordered (strictly increasing) //
/////////////////////////////////////////////////////////////////////////////////////
valarray<int> VectorUnionOrdered(const valarray<int> & v, const valarray<int> & w);


/////////////////////////////////////////////////////////////////////////
// Finds the ordered (strictly increasing) intersection of two vectors //
/////////////////////////////////////////////////////////////////////////
valarray<int> VectorIntersection(const valarray<int> & v, const valarray<int> & w);


/////////////////////////////////////////////////////////////////////////////////
// Finds the ordered (strictly increasing) intersection of two ordered vectors //
/////////////////////////////////////////////////////////////////////////////////
valarray<int> VectorIntersectionOrdered(const valarray<int> & v, const valarray<int> & w);


////////////////////////////////////
// Remove the entries in w from v //
////////////////////////////////////
valarray<int> VectorRemove(const valarray<int> & v, const valarray<int> & w);


//////////////////////////////////////////////////////////////
// Returns an ordered vector with the entries of v not in w //
//////////////////////////////////////////////////////////////
valarray<int> VectorComplement(const valarray<int> & v, const valarray<int> & w);


///////////////////////////////////////////////////////////////
// Returns an ordered vector with the entries of v not in w, //
// assuming both v and w are ordered (strictly increasing)   //
///////////////////////////////////////////////////////////////
valarray<int> VectorComplementOrdered(const valarray<int> & v, const valarray<int> & w);


/////////////////////////////////////////
// Checks if the two vectors intersect //
/////////////////////////////////////////
bool CheckVectorIntersection(const valarray<int> & v, const valarray<int> & w);


///////////////////////////////////////////////////
// Checks if the two (ordered) vectors intersect //
///////////////////////////////////////////////////
bool CheckVectorIntersectionOrdered(const valarray<int> & v, const valarray<int> & w);


/////////////////////////////////////////////////
// Checks if v is a strictly increasing vector //
/////////////////////////////////////////////////
bool CheckVectorOrdered(const valarray<int> & v);



//////////////////////////////////////////
// Finds the GCD of a vector of numbers //
//////////////////////////////////////////
mpz_class GCD(const valarray<mpz_class> & v);




#endif



