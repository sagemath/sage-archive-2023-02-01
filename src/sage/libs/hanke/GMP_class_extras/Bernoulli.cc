
////////////////////////////////////////
// Computes the n-th Bernoulli vector //
////////////////////////////////////////

void BernoulliNumberVector(valarray<mpq_class> & Bvec, long n) {

  if (n < 0)
    cout << " \n Error in BernoulliVector:  n = " \
	 << n << " < 0. \n";

  if (n == 0) {
    Bvec.resize(1);
    Bvec[0] = mpq_class(1);
  }

  if (n >= 1) {
    Bvec.resize(n+1);
    long i,j;
    mpz_t temp_binom;
    mpz_init (temp_binom);
    Bvec[0] = mpq_class(1);

    for (i=1; i<=n; i++) {
      Bvec[i] = mpq_class(0);
      for (j=0; j<i; j++) {
	mpz_bin_uiui(temp_binom, i+1, j);
	Bvec[i] += Bvec[j] * mpq_class(mpz_class(temp_binom));
      }
      Bvec[i] *= mpq_class(-1, i+1);
    }
  }
}



////////////////////////////////////////
// Computes the n-th Bernoulli number //
////////////////////////////////////////

mpq_class BernoulliNumber(long n) {

  valarray<mpq_class> bern_vector(0);
  BernoulliNumberVector(bern_vector, n);

  return(bern_vector[n]);
}



////////////////////////////////////////////
// Computes the n-th Bernoulli polynomial //
////////////////////////////////////////////

void BernoulliPolynomial(valarray<mpq_class> & bernpoly, long n) {

  valarray<mpq_class> bern_vec;
  BernoulliNumberVector(bern_vec, n);

  // To Do: Add an error message here if n<0... =)

  // Compute B_n(x) for n>=0
  if (n >= 0) {
    bernpoly.resize(n+1);

    long i;
    mpz_t temp_binom;
    mpz_init (temp_binom);

    for (i=0; i<=n; i++) {
      mpz_bin_uiui(temp_binom, n, i);
      bernpoly[n - i] = bern_vec[i] * mpq_class(mpz_class(temp_binom));
    }
  }
}


/////////////////////////////////////////////////////////////////////////////////
// Evaluate a rational "polynomial" (given as a valarray) at a rational number //
/////////////////////////////////////////////////////////////////////////////////

mpq_class EvaluatePolynomial(const valarray<mpq_class> & Coeffs, mpq_class num){

  //    cout << "a" << endl;
  mpq_canonicalize(num.get_mpq_t());
  //    cout << "ab" << endl;

  unsigned int i;
  mpq_class total;
  total = mpq_class(0);

  //  cout << " Polynomial size is " << Coeffs.size() << endl;


  for (i=0; i<Coeffs.size(); i++) {
    total += Coeffs[i] * RationalPower(num, i);
  }

  return total;
}




////////////////////////////////////////////////////////////////////////////
// Computes the k-th generalized Bernoulli number for the character (d/.) //
////////////////////////////////////////////////////////////////////////////

mpq_class QuadraticBernoulliNumber(unsigned long k, long d) {

  mpz_class f;
  f = abs(CoreDiscriminant(mpz_class(d)));

  mpq_class total;
  total = mpq_class(0);

  // Make the k-th Bernoulli polynomial
  valarray<mpq_class> bernpoly;
  BernoulliPolynomial(bernpoly, k);

  unsigned int i;
  for (i=1; i<=f; i++)
    total += KroneckerSymbol(mpz_class(d), mpz_class(i)) * EvaluatePolynomial(bernpoly, mpq_class(i,f));
    // Note: This was ok since i<>0!
  total *= mpq_class(f^(k-1));

  return(total);
}




