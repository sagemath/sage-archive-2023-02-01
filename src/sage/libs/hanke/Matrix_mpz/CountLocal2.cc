#include "Matrix_mpz_header.h"
#include "../GMP_class_extras/mpz_class_extras.h"


///////////////////////////////////////////////////////////////////////////////////
// Private routine to check if a given solution vector w (of Q(w) = m mod p^k)   //
// is of a certain local type and satisfies certain congruence conditions mod p. //
//   (Personal Note: For p=2, we should still use p=2 and not p=8.)              //
///////////////////////////////////////////////////////////////////////////////////

int Matrix_mpz::GetLocalSolutionType(const mpz_class & p, const valarray<mpz_class> & w,
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
    while( (i < zero.size()) && ((w[zero[i]-1] % p) == 0) )
      i++;
    if (i == zero.size())
      zero_flag = true;
  }

  /*
  cout << "IsLocalSolutionType: Finished the Zero congruence condition test \n";
  */

  if (zero_flag == false)
    return 0;

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
      if ((w[nonzero[i] - 1] % p) != 0)
	nonzero_flag = true;
      i++;
    }
  }

  if (nonzero_flag == false)
    return 0;


  // Check if the solution has the appropriate (local) type


  // 1: Check Good-type
  for (i=1; i <= w.size(); i++)
    if (((w[i-1] % p) != 0)  && (((*this)(i,i) % p) != 0))
      return 1;
  if (p==2)
    for (i=1; i < w.size(); i++)
      if ((((*this)(i,i+1) % p) != 0) && (((w[i-1] % p) != 0) || ((w[i] % p) != 0)))
	return 1;


  // 2: Check Zero-type
  bool Zero_flag = true;
  for (i=1; i <= w.size(); i++)
    if ((w[i-1] % p) != 0)
      Zero_flag = false;
  if (Zero_flag == true)
    return 2;


  /*  This serves no purpose...
  bool is_bad_type = false;
  i=1;
  while ((is_bad_type == false) && (i <= w.size())) {
    if (((w[i-1] % p) != 0) && (Q(i,i) % p == 0))
      is_bad_type = true;
    i++;
  }

  if (is_bad_type == true) {
  */


  // Check if wS1 is zero or not
  bool wS1_nonzero_flag;
  wS1_nonzero_flag = false;
  for (i=1; i<=(*this).NumRows(); i++) {
    mpz_class val;

    // Compute the valuation of each index, allowing for off-diagonal terms
    if ((*this)(i,i) == 0)
      if (i == 1)
	val = Valuation((*this)(i,i+1), p);  // Look at the term to the right
      else if (i == (*this).NumRows())
	val = Valuation((*this)(i-1,i), p);  // Look at the term above
      else
	val = Valuation((*this)(i,i+1) + (*this)(i-1,i), p);  // Finds the valuation of the off-diagonal term since only one isn't zero
    else
      val = Valuation((*this)(i,i), p);

    // Test each index
    if ((val == 1) && ((w[i-1] % p) != 0))
      wS1_nonzero_flag = true;
  }


  // 4: Check Bad-type I
  if (wS1_nonzero_flag == true) {
    //cout << " Bad I Soln :  " << w;
    return 4;
  }

  //    cout << " Bad II Soln :  " << w << "  wS1_nonzero_flag = " << wS1_nonzero_flag << endl;

  // 5: Check Bad-type II
  if (wS1_nonzero_flag == false) {
    //cout << " Bad II Soln :  " << w;
    return 5;
  }



    /*
  }
    */


  cerr << "\n Error in IsLocalSolutionType: Should not execute this line... =( \n";
  cerr << "   Solution vector is " << w << "\n";
  cerr << "   and Q is \n" << (*this) << "\n\n" << endl;
  abort();
}



//////////////////////////////////////////////////////////////////
// Naively counts the number of solutions of Q(x) = m mod p^k   //
// of type solntype, satisfying the mod p congruence conditions //
// at the indices of the vectors "zero" and "nonzero"           //
//////////////////////////////////////////////////////////////////

valarray <mpz_class> Matrix_mpz::CountAllLocalTypesNaive(const mpz_class & p, unsigned long k, const mpz_class & m,
					     const valarray<int> & zero, const valarray<int> & nonzero) const
{

  /*
  cout << "   --> CountAllLocalTypesNaive is using the form Q " << endl << Q << endl;
  cout << "       p = " << p << "  and   m = " << m << endl;
  */

  //  if (Q.IsSymmetric()) {
  if ((*this).IsSymmetric() || (p==2)) {
    unsigned long n;
    n = (*this).NumRows();

    mpz_class R;
    R = p^k;

    // Initialize the counting vector
    valarray <mpz_class> count_vector(6);
    for (int i=0; i<=5; i++)
      count_vector[i] = 0;


    int solntype;

    // Initialize v = (0, ... , 0)
    valarray<mpz_class> v(n);
    for (int i=0; i<=n-1; i++)
      v[i] = 0;


    /*
      for(int i=0; i<n; i++) // This may be redundant... =(
      v[i] = 0;
    */

    // Some declarations to speed up the loop
    mpz_class R_n;
    R_n = R^n;
    mpz_class m1;
    m1 = m % R;

    // Count the local solutions
    for(int i=1; i<=(R_n).get_ui(); i++) {
      Increment(v,R);
      if ((*this).EvaluateQuadratic(v,R) == m1 ) {
	solntype = (*this).GetLocalSolutionType(p, v, zero, nonzero);
	if (solntype != 0)
	  count_vector[solntype]++;
      }
    }


    // Generate the Bad-type and Total counts
    count_vector[3] = count_vector[4] + count_vector[5];
    count_vector[0] = count_vector[1] + count_vector[2] + count_vector[3];


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

    return count_vector;
  }
  else
    cout << "Error in CountLocalTypeNaive: Matrix \n" << (*this) << "\n is not symmetric!";

  cout << "Ungraceful exit..." << endl;
  exit(0);
}


//////////////////////////////////////////
// Front-ends for our counting routines //
//////////////////////////////////////////

mpz_class Matrix_mpz::CountLocalType(const mpz_class & p, long k, const mpz_class & m, int solntype,
			 const valarray<int> & zero, const valarray<int> & nonzero) const
{
  // Ideally this would use the CountLocalTypeWithSymmetry routine, but this is fine for now. =)
  return (*this).CountAllLocalTypesNaive(p, k, m, zero, nonzero)[0];
}

mpz_class Matrix_mpz::CountLocalGoodType(const mpz_class & p, long k, const mpz_class & m,
			     const valarray<int> & zero, const valarray<int> & nonzero) const
{
  return (*this).CountAllLocalTypesNaive(p, k, m, zero, nonzero)[1];
}

mpz_class Matrix_mpz::CountLocalZeroType(const mpz_class & p, long k, const mpz_class & m,
			     const valarray<int> & zero, const valarray<int> & nonzero) const
{
  return (*this).CountAllLocalTypesNaive(p, k, m, zero, nonzero)[2];
}

mpz_class Matrix_mpz::CountLocalBadType(const mpz_class & p, long k, const mpz_class & m,
			    const valarray<int> & zero, const valarray<int> & nonzero) const
{
  return (*this).CountAllLocalTypesNaive(p, k, m, zero, nonzero)[3];
}

mpz_class Matrix_mpz::CountLocalBadTypeI(const mpz_class & p, long k, const mpz_class & m,
			     const valarray<int> & zero, const valarray<int> & nonzero) const
{
  return (*this).CountAllLocalTypesNaive(p, k, m, zero, nonzero)[4];
}

mpz_class Matrix_mpz::CountLocalBadTypeII(const mpz_class & p, long k, const mpz_class & m,
			      const valarray<int> & zero, const valarray<int> & nonzero) const
{
  return (*this).CountAllLocalTypesNaive(p, k, m, zero, nonzero)[5];
}



