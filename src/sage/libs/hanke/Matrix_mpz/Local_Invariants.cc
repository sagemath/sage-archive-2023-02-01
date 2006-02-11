#include "Matrix_mpz_header.h"
#include "../GMP_class_extras/mpz_class_extras.h"
#include "../GMP_class_extras/vectors.h"

// Routines to compute local (p-adic) invariants of a quadratic form Q:
// (Note: Here Q is the matrix so that Q(x) = x^t * Q * x.)
// --------------------------------------------------------------------


//////////////////////////////////////////////////////////////////////
// Finds a diagonal form equivalent to Q over the rational numbers. //
//   (It does this using an integral matrix...)                     //
//////////////////////////////////////////////////////////////////////
Matrix_mpz Matrix_mpz::RationalDiagonal() const {

  // TO DO: Should perform some checks to make sure that Q is a square
  // symmetric matrix.


  // Copy the current matrix into Q
  Matrix_mpz Q;
  Q = (*this);

  size_t n;
  n = Q.NumRows();

  // Construct an integral change of basis matrix T
  // so that T^t * Q * T is diagonal.
  size_t i,j;
  Matrix_mpz temp;
  /*
  Matrix_mpz T;                // Note: T is unnecessary here, but it may be useful to keep track of later.
  T = IdentityMatrix(n);
  */
  mpz_class gcd;

  for(i=1; i<=n; i++) {
    temp.IdentityMatrix(n);
    gcd = GCD(Q.ExtractRow(i));

    for(j=i+1; j<=n; j++) {
      temp(j,j) = Q(i,i);
      temp(i,j) = -Q(i,j);
    }

    /*
    for(j=i+1; j<=n; j++) {
      temp(j,j) = Q(i,i) / gcd;
      temp(i,j) = -Q(i,j) / gcd;
    }
    */

    /*
    cout << endl;
    cout << " Row " << i << ":  Temp matrix is " << endl;
    cout << temp << endl;
    cout << " This changes Q from " << endl;
      cout << Q << endl;
    */
    Q = Q.EvaluateQuadratic(temp);
    /*
    cout << " to " << endl;
    cout << Q << endl;
    */

    /*
    T = T * temp;        // Note: T is unnecessary here, but it may be useful to keep track of later.
    */
  }

  //cout << "\n T is \n" << T << "\n";

  // Thought: Check for square factors in the diagonal elements
  // (This is unnecessary, but it simplifies the answer)

  return Q;
}


////////////////////////////////////////////////////////////////////////
// Finds a diagonal form equivalent to Q over the p-adic numbers Q_p. //
////////////////////////////////////////////////////////////////////////
Matrix_mpz Matrix_mpz::LocalDiagonal(const mpz_class & p) const{
  if (p != 2)
    return (*this).GetLocalNormal(p);  // This is a diagonal matrix for Q over Z_p. =)
  else {
    Matrix_mpz Q2 = (*this).GetLocalNormal(p); // Note: This is upper-triangular

    for(size_t i=1; i<Q2.NumRows(); i++) {
      if (Q2(i, i+1) != 0) {
	// This corresponds to [0 1]  =  [0  1/2]   ===>>>  [1  0]
	//                     [0 0]     [1/2  0]           [0 -1]
	if (Q2(i,i) == 0) {
	  Q2(i, i) = Q2(i, i+1);
	  Q2(i+1, i+1) = -Q2(i, i+1);
	  Q2(i, i+1) = 0;
	}

	// This corresponds to [1 1]  =  [1  1/2]   ===>>>  [1 0]
	//                     [0 1]     [1/2  1]           [0 3]
	else {
	  Q2(i, i+1) = 0;
	  Q2(i+1, i+1) *= 3;
	}
      }
    }

    return Q2;
  }
}




/////////////////////////////////////////////////////
// Computes the Hasse invariant of Q at a prime p. //
/////////////////////////////////////////////////////
long Matrix_mpz::HasseInvariant(const mpz_class & p) const {

  // To Do: Need to deal with the case n=1 separately somewhere!

  Matrix_mpz Diag;
  Diag = (*this).LocalDiagonal(p);

   /*
  cout << "\n Q = " << (*this) << endl;
  cout << "\n Q diagonalized at p = " << p << " gives " << Diag << endl;
   */

  size_t j, k, n;
  long hasse_temp;
  hasse_temp = 1;
  n = Diag.NumRows();

  for (j=1; j<n; j++)
    for (k=j+1; k<=n; k++)
      hasse_temp = hasse_temp * HilbertSymbol(Diag(j,j), Diag(k,k), p);

  return hasse_temp;
}




/////////////////////////////////////////////////////////////
// Checks if Q is anisotropic over the p-adic numbers Q_p. //
/////////////////////////////////////////////////////////////
bool Matrix_mpz::IsAnisotropic(const mpz_class & p) const {

  size_t n = (*this).NumRows();
  mpz_class D = (*this).Determinant();

  // Should check that p is prime and Q is square

  if (n>=5)
    return false;

  if (n==4)
    return (IsPadicSquare(D, p) && ((*this).HasseInvariant(p) == - HilbertSymbol(-1,-1,p)) );

  if (n==3)
    return ((*this).HasseInvariant(p) != HilbertSymbol(-1, -D, p));

  if (n==2)
    return (!IsPadicSquare(-D, p));

  if (n==1)
    return ((*this)(1,1) != 0);

  cerr << "Error in IsAnisotropic " << endl;
  abort();
}



///////////////////////////////////////////////////////////
// Checks if Q is isotropic over the p-adic numbers Q_p. //
///////////////////////////////////////////////////////////
bool Matrix_mpz::IsIsotropic(const mpz_class & p) const {

  return (!(*this).IsAnisotropic(p));

  cout << "\n Error in IsIsotropic:  This case should never occur!\n";
}





//////////////////////////////////////////////////////////////
// Returns a vector with all of the anisotropic primes of Q //
//////////////////////////////////////////////////////////////
valarray<mpz_class> Matrix_mpz::AnisotropicPrimes() const {

  // Look at all prime divisors of 2 * Det(Q) to find the anisotropic primes...

  valarray<mpz_class> AnisoPrimes, possible_primes;

  possible_primes = PrimeDivisors(2 * (*this).Determinant());
  AnisoPrimes.resize(possible_primes.size());

  //cout << " Possible anisotropic primes are: " << possible_primes << endl;

  size_t ptr;
  ptr = 0;

  for(size_t i=0; i < possible_primes.size(); i++)
    if ((*this).IsAnisotropic(possible_primes[i])) {
      AnisoPrimes[ptr] = possible_primes[i];
      ptr++;
    };

  //cout << " leaving AnisotropicPrimes..." << endl;

  return VectorTrim(AnisoPrimes);
}



