

// --------------------------------------------------------------


// Print the vector V
void PrintV(vector<long> V);


// Print the first n entries of a vector V
void PrintHeadV(vector<long> V, unsigned long n);
void PrintHeadV(vector<double> V, unsigned long n);


// Print the last n entries of a vector V
void PrintTailV(vector<long> V, unsigned long n);
void PrintTailV(vector<double> V, unsigned long n);


// Reads in a vector of longs
vector<long> ReadVector_long(const char* filename);

// Reads in a vector of doubles
vector<double> ReadVector_double(const char* filename);



// ---------------------------------------------------

#if !defined(VEC_UTILS_H)
#define VEC_UTILS_H

#include "vec_utils.cc"

#endif
