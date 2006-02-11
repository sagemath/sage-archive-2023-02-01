#include "Matrix_mpz_header.h"
#include "../GMP_class_extras/mpz_class_extras.h"

// This is needed in the filter for primitivity...
#include "../max-min.h"



  ////////////////////////////////
  // Private Front-end Routines //
/////////////////////////////////////////////////////////////////////////////////////////////


//////////////////////////////////////////////////////////////////
// Finds the Good-type local density of Q representing m at p.  //
// (Front end routine for its congruence counterpart.)          //
//////////////////////////////////////////////////////////////////

mpq_class Matrix_mpz::Local_Good_Density(const mpz_class & p, const mpz_class & m) const
{
  valarray<int> EmptyVec;
  //  cout << "Doing Good Density with p = " << p << " and m = " << m << endl;
  return (*this).Local_Good_Density_Congruence(p, m, EmptyVec, EmptyVec);
}


//////////////////////////////////////////////////////////////////
// Finds the Zero-type local density of Q representing m at p.  //
// (Front end routine for its congruence counterpart.)          //
//////////////////////////////////////////////////////////////////

mpq_class Matrix_mpz::Local_Zero_Density(const mpz_class & p, const mpz_class & m) const
{
  valarray<int> EmptyVec;
  //  cout << "Doing Zero Density with p = " << p << " and m = " << m << endl;
  return (*this).Local_Zero_Density_Congruence(p, m, EmptyVec, EmptyVec);
}



//////////////////////////////////////////////////////////////////
// Finds the Bad-type local density of Q representing m at p.  //
// (Assuming that p > 2 and Q is given in local diagonal form.) //
//////////////////////////////////////////////////////////////////

mpq_class Matrix_mpz::Local_Bad_Density(const mpz_class & p,const mpz_class & m) const
{
  valarray<int> EmptyVec;
  //  cout << "Doing Bad Density with p = " << p << " and m = " << m << endl;
  return (*this).Local_Bad_Density_Congruence(p, m, EmptyVec, EmptyVec);
}


mpq_class Matrix_mpz::Local_BadI_Density(const mpz_class & p,const mpz_class & m) const
{
  valarray<int> EmptyVec;
  //  cout << "Doing BadI Density with p = " << p << " and m = " << m << endl;
  return (*this).Local_BadI_Density_Congruence(p, m, EmptyVec, EmptyVec);
}


mpq_class Matrix_mpz::Local_BadII_Density(const mpz_class & p,const mpz_class & m) const
{
  valarray<int> EmptyVec;
  //  cout << "Doing BadII Density with p = " << p << " and m = " << m << endl;
  return (*this).Local_BadII_Density_Congruence(p, m, EmptyVec, EmptyVec);
}





// ---------------  These are the important ones, which we'll filter for primitive forms!!!  ------------------


mpq_class Matrix_mpz::Local_Density(const mpz_class & p, const mpz_class & m) const {


  // Set the modulus to check for imprimitive forms at p
  unsigned long p_valuation;      // This assumes that Q has size at least 1! =)
  bool no_val_flag = true;

  // Check for imprimitive forms at p -- ASSUMES THE FORM IS UPPER TRIANGULAR!
  // (NOTE: We could do better if we know it's normalized!)
  for(long i=1; i<=(*this).NumRows(); i++)
    for(long j=i; j<=(*this).NumRows(); j++)
      if ((*this)(i,j) != 0)
	if (no_val_flag == true) {
	  no_val_flag = false;
	  p_valuation = Valuation((*this)(i,j), p);
	}
	else
	  p_valuation = min(p_valuation, Valuation((*this)(i,j), p));


  /*
  cout << " Using the matrix: \n " << (*this) << endl;
  cout << "Valuation(m,p) = " << Valuation(m,p) << endl;
  cout << "p_valuation = " << p_valuation << endl;
  */


  // If m is less p-divisible than the matrix, return zero
  if ((m != 0) && (Valuation(m,p) < p_valuation))   // Note: The (m != 0) condition protects taking the valuation of zero.
    return mpq_class(0);

  // If the form is imprimitive, divide it (and m) by p-powers to get a primitive form
  else {
    if (p_valuation > 0) {

      // Make a new (primitive) matrix
      Matrix_mpz QQ;
      QQ = (*this);
      mpz_class p_mod;

      /*
      // DIAGNOSTIC
      cout << " p = " << p << endl;
      cout << " p_valuation = " << p_valuation << endl;
      cout << " p_mod = " << (p^p_valuation) << endl;
      */

      p_mod = p ^ p_valuation;  // This should give a power...
      for(long i=1; i<=(*this).NumRows(); i++)
	for(long j=i; j<=(*this).NumRows(); j++)
	  QQ(i,j) = (*this)(i,j) / p_mod;

      // Make a new number mm
      mpz_class mm;
      mm = m / p_mod;

      // Then return the densities for the reduced problem
      return QQ.Local_Good_Density(p, mm) + QQ.Local_Zero_Density(p, mm) + QQ.Local_Bad_Density(p, mm);
    }

    // Otherwise, proceed as usual... =)
    else
      return (*this).Local_Good_Density(p, m) + (*this).Local_Zero_Density(p, m) + (*this).Local_Bad_Density(p, m);
  }


  // Diagnostic to trap any missing conditions...
  cout << "Error in Local_Density():  We should never reach this line! =( ..." << endl;
  exit(1);

}




mpq_class Matrix_mpz::Local_Primitive_Density(const mpz_class & p, const mpz_class & m) const {



  // Set the modulus to check for imprimitive forms at p
  unsigned long p_valuation;      // This assumes that Q has size at least 1! =)
  bool no_val_flag = true;

  // Check for imprimitive forms at p -- ASSUMES THE FORM IS UPPER TRIANGULAR!
  // (NOTE: We could do better if we know it's normalized!)
  for(long i=1; i<=(*this).NumRows(); i++)
    for(long j=i; j<=(*this).NumRows(); j++)
      if ((*this)(i,j) != 0)
	if (no_val_flag == true) {
	  no_val_flag = false;
	  p_valuation = Valuation((*this)(i,j), p);
	}
	else
	  p_valuation = min(p_valuation, Valuation((*this)(i,j), p));


  /*
  cout << " Using the matrix: \n " << (*this) << endl;
  cout << "Valuation(m,p) = " << Valuation(m,p) << endl;
  cout << "p_valuation = " << p_valuation << endl;
  */


  // If m is less p-divisible than the matrix, return zero
  if ((m != 0) && (Valuation(m,p) < p_valuation))   // Note: The (m != 0) condition protects taking the valuation of zero.
    return mpq_class(0);

  // If the form is imprimitive, divide it (and m) by p-powers to get a primitive form
  else {
    if (p_valuation > 0) {

      // Make a new (primitive) matrix
      Matrix_mpz QQ;
      QQ = (*this);
      mpz_class p_mod;

      /*
      // DIAGNOSTIC
      cout << " p = " << p << endl;
      cout << " p_valuation = " << p_valuation << endl;
      cout << " p_mod = " << (p^p_valuation) << endl;
      */

      p_mod = p ^ p_valuation;  // This should give a power...
      for(long i=1; i<=(*this).NumRows(); i++)
	for(long j=i; j<=(*this).NumRows(); j++)
	  QQ(i,j) = (*this)(i,j) / p_mod;

      // Make a new number mm
      mpz_class mm;
      mm = m / p_mod;

      // Then return the densities for the reduced problem
      return QQ.Local_Good_Density(p, mm) + QQ.Local_Bad_Density(p, mm);
    }

    // Otherwise, proceed as usual... =)
    else
      return (*this).Local_Good_Density(p, m) + (*this).Local_Bad_Density(p, m);
  }


  // Diagnostic to trap any missing conditions...
  cout << "Error in Local_Primitive_Density():  We should never reach this line! =( ..." << endl;
  exit(1);

}


