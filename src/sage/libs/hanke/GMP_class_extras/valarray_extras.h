#include <iostream>

#if !defined(VALARRAY_EXTRAS_H)
#define VALARRAY_EXTRAS_H



// Easy printing routine for my vectors... =)  -- SHOULD PHASE THIS OUT!
template<class T>
void PrintV(const valarray<T> & v) {
  cout << "[ ";
  for(size_t i = 0; i < v.size(); i++) {
    cout << v[i];
    if (i < v.size() - 1)
      cout << ", ";
  }
  cout << " ]" << endl;
}




//////////////////////////////////////////////////////////////////////////////////



/////////////////////////////////////
// Printing Routines for valarrays //
/////////////////////////////////////


// First part of a valaray<T> printing routine
template<class T>
void PrintVector(ostream & out, valarray<T> v) {

  /*
    cout << " m = " << m << "  n = " << n << endl;
    cout << " length of the valarray = " << M.size() << endl;
  */

  out << " [ ";
  for(size_t j = 0; j < v.size(); j++) {
    out << v[j];
    if (j < (v.size() - 1))
      out << ", ";
  }
  out << " ]" << endl;
}


// Define the << operator for the valarray<T> type
template<class T>
ostream & operator<<(ostream & out, valarray<T> v) {    // Why doesn't "(...., Matrix_mpz & matr)" work?  =|
  PrintVector(out, v);
  return out;
}



////////////////
// VectorTrim //
////////////////


// Trim the non-zero elements from a valarray<T> -- This should be moved to vectors.cc
template <class T>
valarray<T> VectorTrim(valarray<T> & v) {

  size_t n = v.size();
  size_t i = 0;

  while ((i<n) && (v[i] != 0))
    i++;

  valarray<T> newvec(i);
  size_t j;

  for (j=0; j<i; j++)
    newvec[j] = v[j];

  return newvec;
}





// ==========================================================

#include "../max-min.h"

mpq_class Minimum(const valarray<mpq_class> & vec) ;




#endif


