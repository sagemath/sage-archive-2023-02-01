





/*
  This computes an upper bound for the precision needed for the
  ternary form to check that our given 4 variable form is represented.
  This bound only depends on the splitting coefficient \f$d\f$ in the
  decomposition \f$Q = T \oplus d w^2\f$.

  The proof of our upper bound for \f$T\f$ is
  \f[
  B > F_4(T) = \prod_p F_4(p)
             = \prod_p \left( \frac{\sqrt{p}}2 \cdot adj(p) \right)
  \hfill \Rightarrow \hfill
   T = \prod_p p \leq \left( \frac{2^{\# p|T} B} {adj(p)} \right)^2.
  \f]

  To bound the maximal difference, we assume that \f$ T = dx^2 \f$
  (which is the worst case) and that we will attempt to find the
  minimal difference \f$t\f$ times.  Then
  \f[
  T - d(x-t)^2 = dx^2 - d(x-t)^2 = d(2tx - t^2) < 2tdx = 2t\sqrt{Td},
  \f]
  since \f$ x = \sqrt{\frac Td} \f$.





*/











/*!  \brief Computes the product of all local densities for comparison with independently computed Eisenstein coefficients.
 *
 *  \todo We fixed the generic factors to compensate for using tha matrix of 2Q, but we need to document this better! =)
 */

////////////////////////////////////////////////////////////////
// Computes the infinite product of local densities of Q at u //
////////////////////////////////////////////////////////////////

mpq_class Matrix_mpz::SiegelProduct(mpz_class u) const {

  size_t n;
  mpz_class d, S;
  n = (*this).NumRows();
  d = (*this).Determinant();      // Warning: This is a factor of 2^n larger than it should be!

  //    cout << "In SiegelProduct:  d = " << d << endl;


  // Product of "bad" places to omit
  S = 2 * d * u;


  size_t m;
  mpz_class d1, f;
  mpq_class genericfactor;

  /*
  cout << "SiegelProduct Break 1. " << endl;
  cout << " u = " << u << endl;
  */

  // Make the odd generic factors
  if ((n % 2) == 1) {
    m = (n-1) / 2;
    d1 = CoreDiscriminant((mpz_class(-1)^m) * 2*d * u);    // Replaced d by 2d here to compensate for the determinant
    f = Abs(d1);                                           // gaining an odd power of 2 by using th matrix of 2Q instead
                                                           // of the matrix of Q.
                                                           //  --> Old d1 = CoreDiscriminant((mpz_class(-1)^m) * d * u);

    // Make the ratio of factorials factor: [(2m)! / m!] * prod_{i=1}^m (2*i-1)
    size_t i;
    mpq_class factor1;
    factor1 = mpq_class(1);
    for (i=1; i<=m; i++)
      factor1 *= 2*i - 1;
    for (i=m+1; i<=2*m; i++)
      factor1 *= i;


    genericfactor = factor1 * RationalPower(mpq_class(u, f), m) *
      RationalSqrt(mpq_class((mpz_class(2)^n) *  f, u * d)) *
      Abs(QuadraticBernoulliNumber(m, d1.get_si()) / BernoulliNumber(2*m));

  }


  //  cout << "SiegelProduct Break 2. " << endl;

  // Make the even generic factor
  if ((n % 2) == 0) {
    m = n / 2;
    d1 = CoreDiscriminant((mpz_class(-1)^m) * d);
    f = Abs(d1);

    /*
        cout << " mpz_class(-1)^m = " << (mpz_class(-1)^m) << " and d = " << d << endl;
        cout << " f = " << f << " and d1 = " << d1 << endl;
    */

    genericfactor = mpq_class(m) / RationalSqrt(f*d)
      * RationalPower(mpq_class(u,2), m-1) * mpq_class(mpz_class(f)^m)
      / mpq_class(Abs(QuadraticBernoulliNumber(m, d1.get_si())))
      * (mpz_class(2) ^ m);                                            // This last factor compensates for using the matrix of 2*Q
  }


  // Omit the generic factors in S and compute them separately
  mpq_class omit, include;
  omit = 1;
  include = 1;

  valarray<mpz_class> S_divisors;
  S_divisors.resize(PrimeDivisors(S).size());
  S_divisors = PrimeDivisors(S);

   /*
  cout << "\n S is " << S << endl;
  cout << " The Prime divisors of S are :";
  PrintV(S_divisors);
    */

  mpz_class p;

  for (size_t i=0; i<S_divisors.size(); i++) {
    Matrix_mpz Q_normal;
    p = S_divisors[i];
    Q_normal = (*this).GetLocalNormal(p);

    //             cout << " p = " << p << " and its Kronecker symbol (d1/p) = (" << d1 << "/" << p << ") is " << KroneckerSymbol(d1, p) << endl;

    omit *= 1 / (1 - mpq_class(KroneckerSymbol(d1, p), p^m));
    omit.canonicalize();

    /*
             cout << " omit = " << omit << endl;
             cout << " Q_normal is \n" << Q_normal << endl;
    */
    /*
    cout << " Q_normal = " << endl << Q_normal;
    cout << " p = " << p << endl;
    cout << " u = " << u << endl;
    cout << " include = " << include << endl;
    */

    include *= Q_normal.Local_Density(p, u);
    include.canonicalize();

    //        cout << " Including the p = " << p << " factor: " << Local_Density(Q_normal, p, u) << endl;

    //        cout << "    ---  Exiting loop " << endl;
  }



  // ****************  Important *******************
  // Additional fix (only included for n=4) to deal
  // with the power of 2 introduced at the real place
  // by working with Q instead of 2*Q.  This needs to
  // be done for all other n as well...
  /*
  if (n==4)
    genericfactor = 4 * genericfactor;
  */


  /*
  cout << endl;
  cout << " generic factor = " << genericfactor << endl;
  cout << " omit = " << omit << endl;
  cout << " include = " << include << endl;
  cout << endl;
  */


  //  cout << "SiegelProduct Break 3. " << endl;

  // Return the final factor (and divide by 2 if n=2)
  if (n==2)
    return (genericfactor * omit * include / 2);
  else
    return (genericfactor * omit * include);

}








    /*
    genericfactor = Factorial(to_ZZ(2 * m)) / Factorial(to_ZZ(m)) \
      * Product(2*i - 1, i=1..m) \  // Need to make this factor before here...
      * (u / f) ^ m * SqrtFraction(f / (u*d)) \  // Need to convert the first part to mpq_class
      * Abs(GeneralizedBernoulliQuadratic(m, d1) / BernoulliNumber(2*m));
    */



    /*
      mpq_class genericfactor
      ZZ CoreDiscriminant(ZZ);

      ZZ Abs(ZZ);
      mpq_class Abs(mpq_class);

      ZZ Factorial(ZZ)

      mpq_class SqrtFraction(mpq_class)

      mpq_class BernoulliNumber(ZZ)
      mpq_class BernoulliNumber(long)

      mpq_class GeneralizedBernoulliQuadratic(ZZ)
      mpq_class GeneralizedBernoulliQuadratic(long)
    */





// Note: This is based on the AnalyzeDensities290() routine in Siegel_Diagnostic/siegel_diagnostic.cc
// To Do: Make this return a bool to say whether it checks out, and be silent otherwise...
void Matrix_mpz::CheckSiegelRange(const long first, const long last, const string & eisfilename) const {

  // Some Sanity Checks
  assert(first <= last);
  assert(last <= 10000);


  cout << " Checking the Siegel Product Formula for coefficients from " << first << " to " << last << ":" << endl;
  cout << " ------------------------------------------------------------------------" << endl;


  // Read in the Eisenstein Series (with the default precision)
  PowerSeries<mpq_class> Eis;
  //  Eis = _GetEisData(eisfilename.c_str(), 10000);         // Use default filename and precision
  Eis = GetMagmaEisSeries(eisfilename.c_str(), 10000);         // Use default filename and precision


  /*
  // Explicitly run through all of the Eisenstein Coefficients
  for (long m = first; m <= last; m++) {
    cout << " m = " << m << "  ==>   Prod = " << (*this).SiegelProduct(m) << "   Eis = " << Eis[(unsigned long) m] << endl;
  }
  */



  // Check the Siegel Product Formula :
  // ---------------------------------
  for (long m = first; m <= last; m++) {
    /*
      // Print the current coefficient (part 1):
      cout << " Checking the coefficient at m = " << m << ": ";
    */
    if ((*this).CheckSiegel(m, Eis) == true) {
      /*
      // Print the current coefficient (part 2):
      cout << "";
      cout << " ok " << endl;
      */
    }
    else {
      cout << " trouble...  at m = " << m << endl;
      assert(0==1);  // Break out of the routine!
      /*
	if (analyze_flag == true)
	break;
      */
    }
  }
  cout << " Finished checking the Siegel Product Formula " << endl;
  cout << endl << endl;


}





// Note: This is based on the AnalyzeDensities290() routine in Siegel_Diagnostic/siegel_diagnostic.cc
// To Do: Make this return a bool to say whether it checks out, and be silent otherwise...
void Matrix_mpz::Check_ComputeTheta_vs_MagmaTheta() const {

  cout << " Checking the theta function coefficients (compute vs. Magma):" << endl;
  cout << " ------------------------------------------------------------------------" << endl;


  // Read in the Magma Theta Series (with the default precision)
  PowerSeries<mpz_class> Theta_Magma;
  Theta_Magma = GetMagmaThetaSeries("", 1000);         // Use default filename and precision


  // Make in the Magma Theta Series (with the default precision)
  PowerSeries<mpz_class> Theta_Computed;
  Theta_Computed = ComputeTheta(1000);         // Use default filename and precisio

  // Check if they agree
  if (Theta_Magma == Theta_Computed)
    cout << " The two series agree! =) " << endl;
  else {
    cout << " The two series don't agree... doing more testing. " << endl;

    // Compare their coefficients
    long i=0;
    while ((i<=1000) && (Theta_Magma[i] == Theta_Computed[i]))
      i++;

    // Write the first mis-match
    if (i<=1000) {
      cout << " Trouble at coefficient " << i << ":" << endl;
      cout << "   Theta_Magma[" << i << "] = " << Theta_Magma[i] << endl;
      cout << "   Theta_Computed[" << i << "] = " << Theta_Computed[i] << endl;
    }

    // Write the first
    cout << endl;
    for(long j=0; j<=min(i+10, 1000); j++)
      cout << "   Theta_Magma[" << j << "] = " << Theta_Magma[j] << "      "
	   << "   Theta_Computed[" << j << "] = " << Theta_Computed[j] << endl;


    // Mis-match ==> Abort
    assert( 0 == 1 );

  }

  cout << " Finished checking the theta functions. " << endl;
  cout << endl << endl;

}
