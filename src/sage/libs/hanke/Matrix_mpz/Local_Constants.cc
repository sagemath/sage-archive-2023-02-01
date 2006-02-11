#include "Matrix_mpz_header.h"
#include "../GMP_class_extras/mpz_class_extras.h"


///////////////////////////////////////////////////////////////
// Compute the local constants C_p(T) assuming T is p-stable //
///////////////////////////////////////////////////////////////

mpq_class Matrix_mpz::LocalConstantCp(const mpz_class & p, const mpz_class & T) const {

  // Find the dimension of Q assuming it's square
  size_t n;
  n = (*this).NumRows();

  Matrix_mpz Q_normal;
  Q_normal = (*this).GetLocalNormal(p);

  // Note : This only works for n >= 2  ***  This could also be sped up since the Local_Density call is a bit redundant...
  return mpq_class(p^(n-2), (p^(n-2)) - 1) *
        Q_normal.Local_Primitive_Density(p, T*(p*p)) / Q_normal.Local_Density(p,T);

}



//////////////////////////////////////////////////////
// Checks if T is p-stable for the quadratic form Q //
//////////////////////////////////////////////////////

bool Matrix_mpz::IsStable(const mpz_class & p, const mpz_class & T) const {

  bool flag;
  flag = true;

  size_t Stable;
  //  Stable := Ceiling(Valuation(Q.QFLevel(), p) / 2);
  Stable = ((Valuation((*this).QFLevel(), p) + 1) / 2);
                                                  // In fact, this is bigger than we need -- since it would suffice
                                                  // without multipling by T below, but to make it exact is more work
                                                  // with no clear benefit.

  Matrix_mpz Q_normal;
  Q_normal = (*this).GetLocalNormal(p);

  mpq_class Kmax, Kmid;
  Kmax = Q_normal.Local_Primitive_Density(p, T * (p^(2 * Stable)));

  for (size_t i=0; i<=Stable; i++){
    Kmid = Q_normal.Local_Primitive_Density(p, T * (p^(2*i)));
    if (Kmid != Kmax) {
      cout << " Failed at " << T * (p^(2*i)) << ":  Here the stable number is"
	   << T * (p^(2 * Stable)) << " with K = " << Kmax << ", while at "
	   <<  T * (p^(2*i)) << " we have K = " << Kmid;
      flag = false;
    }
  }


  /*
  mpq_class K;
  K = Local_Primitive_Density(Q_normal, p, T);
  */
  //return flag, K, Kmax;   // <--- This was the old return line...
  return flag;

}


