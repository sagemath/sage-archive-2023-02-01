 ////////////////////////////////////////////////////////////////////////////////
// Checks the Siegel product at m against the pre-computed Eisenstein series. //
////////////////////////////////////////////////////////////////////////////////

bool Matrix_mpz::CheckSiegel(const mpz_class & m, const PowerSeries<mpq_class> & Eis) const {

  mpq_class prod, eis;

  if ((m <= Eis.Precision()) && (m >= 0)) {
    eis = Eis[m.get_ui()];
    prod = (*this).SiegelProduct(m);
    if ( eis == prod )
      return true;
    else {
      cout << " Warning: The product formula and Eisenstein series don't match at m = " << m << "." << endl;
      cout << endl;
      cout << "   Product formula gives:   " << prod << endl;
      cout << "   Eisenstein series gives: " << eis << endl;
      cout << endl;
      return false;
    }
  }
  else {
    cout << " Error in CheckSiegel: The number " << m
	 << " to check is not in the range of the Eisenstein series coeffs: 0 ... " << Eis.Precision() << endl;
    return false;
  }
}


//////////////////////////////////////////////////
// Checks the local density of Q at p and m
// (using the reduction procedure)
// against the counting points approach...
/////////////////////////////////////////////////

void Matrix_mpz::CheckLocalDensity(const mpz_class & p, const mpz_class & m) const {

  /*
  mpq_class reduction_density;
  reduction_density = Local_Density(Q,p,m);
  */

  int n;
  unsigned long pow;
  n = (*this).NumRows();

  mpz_class N;
  N = (*this).QFLevel();

  long Stable;
  Stable = ((Valuation(N, p) + 1) / 2); 	// In fact, this is bigger than we need -- since it would suffice
                                                // without multipling by T below, but to make it exact is more work
                                                // with no clear benefit.
                                                // (Note: This is forced to be an integer by integer division.)

  mpz_class R;
  char ch;
  ch = 'y';

  cout << "\n\n Checking the local density of Q at the prime p = " << p << " and the number m = " << m << endl;
  cout << " The form Q is \n" << endl << (*this) << endl << endl;
  cout << " Level = " << N << endl;
  cout << "\n The stable exponent is = " << Stable << endl << endl;

  mpz_class denom;
  valarray<int> zero, nonzero;

  R = p;
  pow = 1;
  while (ch == 'y') {

    denom = 1;
    for (int i=1; i<=((n-1)*(pow)); i++)
      denom = denom * p;

    // Count the solutions and densities for the original form Q
    valarray<mpz_class> count_vector;
    count_vector.resize(6);
    count_vector = (*this).CountAllLocalTypesNaive(p, pow, m, zero, nonzero);

    valarray <mpq_class> count_densities;
    count_densities.resize(6);

    for(int l=0; l<=5; l++) {
      count_densities[l] = mpq_class(count_vector[l], denom);
      count_densities[l].canonicalize();
    }

    /*
    cout << " Total (original form) Counting density: " << count_densities[0] << endl;
    cout << "   Good Density = " << count_densities[1] << endl;
    cout << "   Zero Density = " << count_densities[2] << endl;
    cout << "   Bad Density  = " << count_densities[3]  << endl;
    cout << "     Bad I Density  = " << count_densities[4]  << endl;
    cout << "     Bad II Density = " << count_densities[5]  << endl;
    */


    // Count the solutions and densities for the normalized form QQ (for our prime p)
    Matrix_mpz QQ;
    QQ = (*this).GetLocalNormal(p);

    cout << endl;

    cout << " It's local normal form at p=" << p << " is \n" << endl << QQ << endl << endl;

    cout << endl;

    cout << " Results for computations mod " << p << "^" << pow << ":" << endl << endl;


    valarray <mpz_class> count_normalized_vector;
    count_normalized_vector.resize(6);
    count_normalized_vector = QQ.CountAllLocalTypesNaive(p, pow, m, zero, nonzero);

    valarray <mpq_class> count_normalized_densities;
    count_normalized_densities.resize(6);

    for(int l=0; l<=5; l++) {
      count_normalized_densities[l] = mpq_class(count_normalized_vector[l], denom);
      count_normalized_densities[l].canonicalize();
    }

    /*
    cout << endl;

    cout << " Total (normalized form) Counting density: " << count_normalized_densities[0] << endl;
    cout << "   Good Density = " << count_normalized_densities[1] << endl;
    cout << "   Zero Density = " << count_normalized_densities[2] << endl;
    cout << "   Bad Density  = " << count_normalized_densities[3]  << endl;
    cout << "     Bad I Density  = " << count_normalized_densities[4]  << endl;
    cout << "     Bad II Density = " << count_normalized_densities[5]  << endl;


    // Print densities from the reduction procedure
    cout << endl;

    cout << " Total Reduction density: " << Local_Density(QQ,p,m) << endl;
    cout << "   Good Density = " << Local_Good_Density(QQ,p,m) << endl;
    cout << "   Zero Density = " << Local_Zero_Density(QQ,p,m) << endl;
    cout << "   Bad Density  = " << Local_Bad_Density(QQ,p,m) << endl;
    cout << "     Bad I Density  = " << Local_BadI_Density(QQ,p,m) << endl;
    cout << "     Bad II Density = " << Local_BadII_Density(QQ,p,m) << endl;

    cout << endl;
    */


    // New combined output...

    cout << endl;

    mpq_class naive_density;
    naive_density = mpq_class((*this).CountLocalNaive(m,R), denom);
    naive_density.canonicalize();

    cout << "                        Reduction     Normalized  Unnormalized    Naive "<< endl;
    cout << "                        Prodecure       Count        Count        Count "<< endl;
    cout << " Total =            ";
    cout.width(10);
    cout << "   " << QQ.Local_Density(p,m);
    cout.width(10);
    cout << "   " << count_normalized_densities[0];
    cout.width(10);
    cout << "   " << count_densities[0];
    cout.width(10);
    cout << "   " << naive_density
	 << endl;
    cout.width(0);

    cout << "   Good Density =   ";
    cout.width(10);
    cout << "   " << QQ.Local_Good_Density(p,m);
    cout.width(10);
    cout << "   " << count_normalized_densities[1];
    cout.width(10);
    cout << "   " << count_densities[1]
	 << endl;
    cout.width(0);

    cout << "   Zero Density =   ";
    cout.width(10);
    cout << "   " << QQ.Local_Zero_Density(p,m);
    cout.width(10);
    cout << "   " << count_normalized_densities[2];
    cout.width(10);
    cout << "   " << count_densities[2]
	 << endl;
    cout.width(0);

    cout << "   Bad Density =    ";
    cout.width(10);
    cout << "   " << QQ.Local_Bad_Density(p,m);
    cout.width(10);
    cout << "   " << count_normalized_densities[3];
    cout.width(10);
    cout << "   " << count_densities[3]
	 << endl;
    cout.width(0);

    cout << "   BadI Density =   ";
    cout.width(10);
    cout << "   " << QQ.Local_BadI_Density(p,m);
    cout.width(10);
    cout << "   " << count_normalized_densities[4];
    cout.width(10);
    cout << "   " << count_densities[4]
	 << endl;
    cout.width(0);

    cout << "   BadII Density =  ";
    cout.width(10);
    cout << "   " << QQ.Local_BadII_Density(p,m);
    cout.width(10);
    cout << "   " << count_normalized_densities[5];
    cout.width(10);
    cout << "   " << count_densities[5]
	 << endl;
    cout.width(0);





    // Extra stuff about counting...
    valarray <mpz_class> countvec;
    countvec = (*this).CountLocalNaiveValues(R);
    cout << " The count vector mod " << R << " is: " << countvec;

    mpz_class sum;
    sum = 0;
    for (int j=0; j<R; j++)
      sum += countvec[j];

    cout << " The sum of these is: " << sum << endl;


    valarray <mpq_class> naive_densities;
    naive_densities.resize(R.get_si());
    for (int j=0; j<R; j++) {
      naive_densities[j] = mpq_class(countvec[j], denom);
      naive_densities[j].canonicalize();
    }

    cout << " The densities vector mod " << R << " is: " << naive_densities;







    /*

    valarray <mpz_class> ct_vec;
    ct_vec.resize(4);
    mpz_class b5;


    cout << "\n\n Testing the increment routine mod 5 for dim=4: " << endl;
    for (int j=0; j<625; j++) {
      b5 = ct_vec[3] + 5*ct_vec[2] + 25*ct_vec[1] + 125*ct_vec[0];
      if (j != b5)
	cout << " error: increment discrepancy in " << j << "  " << b5 << "  " << ct_vec << endl;

      valarray <mpz_class> neg_ct_vec;
      neg_ct_vec.resize(4);
      for (int k=0; k<4; k++)
	neg_ct_vec[k] = mpz_class(5) - ct_vec[k];

      if (Q.EvaluateQuadratic(ct_vec, mpz_class(5)) != Q.EvaluateQuadratic(neg_ct_vec, mpz_class(5)))
	cout << " error: v <-- > -v discrepancy in "
	     << ct_vec << "  " << Q.EvaluateQuadratic(ct_vec, mpz_class(5)) << "  "
	     << neg_ct_vec << "  " << Q.EvaluateQuadratic(neg_ct_vec, mpz_class(5))
	     << endl;

      Increment(ct_vec, mpz_class(5));
    }

    cout << "\n\n Finished testing the increment routine mod 5 for dim=4:  No errors means it worked! =)" << endl;

    */







    // Ask to continue

    cout << " Increment the power and recompute? (y/n) " << endl;

    cin >> ch;

    pow++;
    R *= p;
  }

  cout << endl;

}


