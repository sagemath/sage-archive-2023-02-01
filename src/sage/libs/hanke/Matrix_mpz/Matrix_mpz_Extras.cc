

#include "Matrix_mpz.h"



///////////////////////////////////////
// Returns the n x n identity matrix //
///////////////////////////////////////
Matrix_mpz IdentityMatrix(size_t n)
{
  Matrix_mpz M;
  M.IdentityMatrix(n);
  return M;
}





///////////////////////////////////////////////////////////////////
// Swap the ith and jth rows and columns of a symmetric matrix S //
///////////////////////////////////////////////////////////////////
Matrix_mpz SwapSymmetric(Matrix_mpz S, size_t i, size_t j)
{
  // Make a copy of S
  Matrix_mpz T;
  T = S;

  // Do the operation
  T.SwapSymmetric(i,j);

  // Return the result
  return T;
}



//////////////////////////////////////////////////////////////////
// Multiply the ith row and column of a symmetric matrix S by c //
//////////////////////////////////////////////////////////////////

Matrix_mpz MultiplySymmetric(Matrix_mpz S, mpz_class c, size_t i)
{
  // Make a copy of S
  Matrix_mpz T;
  T = S;

  // Do the operation
  T.MultiplySymmetric(c,i);

  // Return the result
  return T;
}



//////////////////////////////////////////////////////////////////
// Multiply the ith row and column of a symmetric matrix S by c //
//////////////////////////////////////////////////////////////////

Matrix_mpz DivideSymmetric(Matrix_mpz S, mpz_class c, size_t i)
{
  // Make a copy of S
  Matrix_mpz T;
  T = S;

  // Do the operation
  T.DivideSymmetric(c,i);

  // Return the result
  return T;
}



//////////////////////////////////////////////////////////////
// Replace the ith rows and columns of a symmetric matrix S //
// with the sum of the ith and c * jth rows and columns     //
//////////////////////////////////////////////////////////////

Matrix_mpz AddSymmetric(Matrix_mpz S, mpz_class c, size_t i, size_t j)
{
  // Make a copy of S
  Matrix_mpz T;
  T = S;

  // Do the operation
  T.AddSymmetric(c,i,j);

  // Return the result
  return T;
}



