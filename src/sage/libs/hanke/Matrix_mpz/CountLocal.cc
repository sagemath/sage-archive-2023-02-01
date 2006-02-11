#include "Matrix_mpz_header.h"
#include "../GMP_class_extras/mpz_class_extras.h"


///////////////////////////////////////
// Increment routine for a ZZ_p type //
///////////////////////////////////////

void Matrix_mpz::Increment(valarray<mpz_class> & v, const mpz_class & R) const
{
    int i;
    i = v.size();

    // Do the carry operations
    while ((i > 0) && (v[i-1] == R-1))   // Assumes that all components satisfy 0 <= v[i] <= R-1
      {
	v[i-1] = 0;
	i--;
      }

    // Only increment if we're not already at the zero vector =)
    if (i > 0)
      v[i-1]++;
}





//////////////////////////////////////////////////////////////
// Naively counts the number of solutions of Q(x) = m mod R //
//////////////////////////////////////////////////////////////

mpz_class Matrix_mpz::CountLocalNaive(const mpz_class & m, const mpz_class & R) const
{

  //  if (Q.IsSymmetric()) {
  if ((*this).IsSymmetric() || ((R % 2) == 0)) {
    unsigned long n;
    n = (*this).NumRows();

    mpz_class count;
    count = 0;

    // Initialize v = (0, ... , 0)
    valarray<mpz_class> v;
    v.resize(n);

    int i;
    /*
      for(i=0; i<n; i++) // This may be redundant... =(
      v[i] = 0;
    */

    // Some declarations to speed up the loop
    mpz_class R_n;
    R_n = R^n;
    mpz_class m1;
    m1 = m % R;

    // Count the local solutions
    for(i=1; i<=R_n; i++) {
      Increment(v, R);
      if ( (*this).EvaluateQuadratic(v,R) == m1 ) {
	count++;
      }
    }

    /*
    cout << "R = " << R << "\n";
    cout << "n = " << n << "\n";
    cout << "R_n = " << R_n << "\n";

    for(i=1; i<=25; i++) {
      cout << "v = " << v << "\n";
      Increment(v,R);
    }
    */

    /*
    cout << "Q = " << Q << "\n";
    cout << "v = " << v << "\n";
    cout << "Q * v = " << Q * v<< "\n";
    */


    return count;
  }
  else
    cout << "Error in CountLocalNaive: Matrix \n" << (*this) << "\n is not symmetric!";

  return -100;
}




//////////////////////////////////////////////////////////////
// Naively counts the number of solutions of Q(x) = m mod R //
//////////////////////////////////////////////////////////////

valarray <mpz_class> Matrix_mpz::CountLocalNaiveValues(const mpz_class & R) const
{

  valarray <mpz_class> value_vector;
  value_vector.resize(R.get_ui());


  //  if (Q.IsSymmetric()) {
  if ((*this).IsSymmetric() || ((R %2) == 0)) {
    unsigned long n;
    n = (*this).NumRows();

    mpz_class count;
    count = 0;

    // Initialize v = (0, ... , 0)
    valarray<mpz_class> v;
    v.resize(n);

    int i;
    /*
      for(i=0; i<n; i++) // This may be redundant... =(
      v[i] = 0;
    */

    // Some declarations to speed up the loop
    mpz_class R_n;
    R_n = R^n;

    // Count the local solutions
    for(i=1; i<=R_n; i++) {
      Increment(v, R);
      value_vector[(*this).EvaluateQuadratic(v,R).get_ui()]++;
    }

    /*
    cout << "R = " << R << "\n";
    cout << "n = " << n << "\n";
    cout << "R_n = " << R_n << "\n";

    for(i=1; i<=25; i++) {
      cout << "v = " << v << "\n";
      Increment(v,R);
    }
    */


    /*
    cout << "Q = " << Q << "\n";
    cout << "v = " << v << "\n";
    cout << "Q * v = " << Q * v<< "\n";
    */

    return value_vector;
  }

  cerr << "Error in CountLocalNaiveValues: Matrix \n" << (*this) << "\n is not symmetric!" << endl;
  abort();

}



///////////////////////////////////////////////////////////////////////////////////
// Private routine to check if a given solution vector w (of Q(w) = m mod p^k)   //
// is of a certain local type and satisfies certain congruence conditions mod p. //
//   (Personal Note: For p=2, we should still use p=2 and not p=8.)              //
///////////////////////////////////////////////////////////////////////////////////

bool Matrix_mpz::IsLocalSolutionType(const mpz_class & p, const valarray<mpz_class> & w, int solntype,
			 const valarray<int> & zero, const valarray<int> & nonzero) const
{

  // Note: Here p is assumed to be a prime >= 2, though the routine still works if not...

  // ToDo?: Add a check that Q is square and has the same size as w.


  bool zero_flag = false;        // Tests the zero mod p congruence conditions
  bool nonzero_flag = false;     // Tests the nonzero congruence conditions

  int i;

  // Check if the solution satisfies the "zero" congruence conditions
  // (either "zero" is empty or its components index the zero vector mod p)
  if (zero.size() == 0)
    zero_flag = true;
  else {
    i=0;
    while( (i < zero.size()) && (w[zero[i]-1] % p == 0) )
      i++;
    if (i == zero.size())
      zero_flag = true;
  }

  /*
  cout << "IsLocalSolutionType: Finished the Zero congruence condition test \n";
  */

  if (zero_flag == false)
    return false;

  /*
    cout << "IsLocalSolutionType: Passed the Zero congruence condition test \n";
  */

  // Check if the solution satisfies the "nonzero" congruence conditions
  // (either "nonzero" is empty or its components index a non-zero vector mod p)
  if (nonzero.size() == 0)
    nonzero_flag = true;
  else {
    i=0;
    while ((nonzero_flag == false) && (i < nonzero.size())) {
      if (w[nonzero[i] - 1] % p != 0)
	nonzero_flag = true;
      i++;
    }
  }

  if (nonzero_flag == false)
    return 0;


  // Check if the solution has the appropriate (local) type

  // 0: Any type
  if (solntype == 0)
    return 1;

  // 1: Good-type
  if (solntype == 1) {
    for (i=1; i <= w.size(); i++)
      if ((w[i-1] % p != 0)  && ((*this)(i,i) % p != 0))
	return true;
    return false;
  }

  // 2: Zero-type
  if (solntype == 2) {
    for (i=1; i <= w.size(); i++)
      if (w[i-1] % p != 0)
	return false;
    return true;
  }

  // 3,4,5: All Bad-types
  if (solntype >= 3) {
    bool is_bad_type = false;
    i=1;
    while ((is_bad_type == false) && (i <= w.size())) {
      if ((w[i-1] % p != 0) && ((*this)(i,i) % p == 0))
	is_bad_type = true;
      i++;
    }

    if (is_bad_type == true)
      return false;

    // 3: Bad-type (I or II)
    if (solntype == 3)
      return true;

    // 4 or 5: Preparations for testing
    if ((solntype == 4) || (solntype == 5)) {
      // Check if w_S1 is zero or not
      bool wS1_nonzero_flag;
      wS1_nonzero_flag = false;
      for (i=1; i<=(*this).NumRows(); i++)
	if ((Valuation((*this)(i,i), p) == 1) && (Valuation(w[i-1], p) > 0))
	  wS1_nonzero_flag = true;

      // 4: Bad-type I
      if (solntype == 4)
	if (wS1_nonzero_flag == true)
	  return true;
	else
	  return false;

      // 5: Bad-type II
      if (solntype == 5)
	if (wS1_nonzero_flag == true)
	  return true;
	else
	  return false;
    }
  }

  cerr << "\n Error in IsLocalSolutionType: Should not execute this line... =( \n";
  cerr << "   Solution Type is " << solntype << endl;
  abort();
}



//////////////////////////////////////////////////////////////////
// Naively counts the number of solutions of Q(x) = m mod p^k   //
// of type solntype, satisfying the mod p congruence conditions //
// at the indices of the vectors "zero" and "nonzero"           //
//////////////////////////////////////////////////////////////////

mpz_class Matrix_mpz::CountLocalTypeNaive(const mpz_class & p, unsigned long k, const mpz_class & m, int solntype,
			      const valarray<int> & zero, const valarray<int> & nonzero) const
{

  if (p==2) {
    cout << " Warning: CountLocalTypeNaive has not been updated to deal with non-diagonal forms when p=2" << endl;
    cout << "   Please use CountAllLocalTypesNaive instead. =) " << endl;
  }


  //  if (Q.IsSymmetric()) {
  if ((*this).IsSymmetric() || (p==2)) {
    unsigned long n;
    n = (*this).NumRows();

    mpz_class R;
    R = p^k;

    mpz_class count;
    count = 0;

    // Initialize v = (0, ... , 0)
    valarray<mpz_class> v;
    v.resize(n);

    int i;
    /*
      for(i=0; i<n; i++) // This may be redundant... =(
      v[i] = 0;
    */

    // Some declarations to speed up the loop
    mpz_class R_n;
    R_n = R^n;
    mpz_class m1;
    m1 = m % R;

    // Count the local solutions
    for(i=1; i<=R_n; i++) {
      Increment(v,R);
      if ( ((*this).EvaluateQuadratic(v,R) == m1 ) && (*this).IsLocalSolutionType(p, v, solntype, zero, nonzero) ) {
	count++;
      }
    }

    /*
    cout << "R = " << R << "\n";
    cout << "n = " << n << "\n";
    cout << "R_n = " << R_n << "\n";

    for(i=1; i<=25; i++) {
      cout << "v = " << v << "\n";
      Increment(v,R);
    }

    cout << "Q1 = " << Q1 << "\n";
    cout << "v = " << v << "\n";
    cout << "Q1 * v = " << Q1 * v<< "\n";
    */

    return count;
  }
  else {
    cout << "Error in CountLocalTypeNaive: Matrix \n" << (*this) << "\n is not symmetric!" << endl;
    cout << "Me too!" << endl;
  }

  return -100;
}


