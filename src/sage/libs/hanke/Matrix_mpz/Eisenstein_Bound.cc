
///////////////////////////////////////////////////////////
// Make a vector of primitive square class repn's in Z_p //
///////////////////////////////////////////////////////////

valarray<mpz_class> Matrix_mpz::_PrimSqRepns(const mpz_class & p) const {

  // Should check if p is prime...
  // cout << "Error in UnitSqRepns:  p is not prime!";

  valarray<mpz_class> vec;

  if (p == 2) {

    //    return [1,3,5,7,2,6,10,14];
    vec.resize(8);
    vec[0] = 1;
    vec[1] = 3;
    vec[2] = 5;
    vec[3] = 7;
    vec[4] = 2;
    vec[5] = 6;
    vec[6] = 10;
    vec[7] = 14;
  }

  else {

    //    return [1, NonResidue(p), p, p*NonResidue(p)];
    vec.resize(4);
    vec[0] = 1;
    vec[1] = NonResidue(p);
    vec[2] = p;
    vec[3] = p * NonResidue(p);
  }

  return vec;
}



//////////////////////////////////////////////////////////////////////////////////////
// Takes in a list of representatives for the primitive squareclasses in Z_p and    //
// returns the minimum of beta_p(m) and beta_p(T') C_p(T') within each square class //
//////////////////////////////////////////////////////////////////////////////////////

valarray<mpq_class> Matrix_mpz::_C4_squareclass_constants(const Matrix_mpz & Q, const mpz_class & p, const valarray<mpz_class> & sqclasslist) const {

  mpz_class N;
  N = Q.QFLevel();
  // cout << "Level = " << N << endl;

  size_t list_length;
  list_length = sqclasslist.size();

  valarray<mpq_class> const_list;
  const_list.resize(list_length);

  mpz_class t;
  unsigned long Stable;

  for (size_t i=0; i<list_length; i++) {

    t = sqclasslist[i];

    Stable = ((Valuation(N, p) + 1) / 2); 	// In fact, this is bigger than we need -- since it would suffice
                                                // without multipling by T below, but to make it exact is more work
                                                // with no clear benefit.
                                                // (Note: This is forced to be an integer by integer division.)

    //  density_vector := [ LocalDensity(Q, p, t * p^(2*i)): i in [0..Stable+1] ];  // Should be primitively constant from Stable+1 onwards! =)
    //  prim_density_vector := [ LocalDensityGood(Q, p, t * (p^(2*i))) + LocalDensityBad(Q, p, t * (p^(2*i))) : i in [0..Stable+3] ];
    //  print "For the square-class ", t, "Z_", p, " we have the local density vector ";
    //  print density_vector;
    //  print "   we also have the local primitive density vector ";
    //  print prim_density_vector;
    //  print "";

    /*
      cout << "\n Starting to create the density vector for p = " << p << "\n" << endl;
      cout << " p = " << p << "  and  N = " << N << endl;
      cout << " Using the square-class entry t = " << t << endl;
      cout << " Stable = " << Stable << endl << endl;
    */

    // Assign the product C_p beta_p(T') to the last entry

    valarray<mpq_class> density_vector;
    density_vector.resize(Stable+3);


    Matrix_mpz Q_normal;
    Q_normal = Q.GetLocalNormal(p);


    if (Q.IsIsotropic(p) == true) {

      for (size_t j=0; j<=Stable+2; j++)
	density_vector[j] = Q_normal.Local_Density(p, t * (p^(2*j)));  // Should be primitively constant from Stable+1 onwards! =)


      mpq_class Cp, Cp1;
      Cp = ((density_vector[Stable+2] / density_vector[Stable+1]) * (p*p) - 1) / ((p*p) - 1);
      Cp1 = min(Cp, 1);


      // Check that we don't have any anisotropic primes by mistake!
      assert (Cp1 > 0);


      /*
      // DIAGNOSTIC -- Outputs the square class info (Part 1):
      // -----------------------------------------------------
      cout << " Stable = " << Stable << endl;
      cout << " density_vector = " << density_vector << endl;
      cout << " C_p(T') = " << Cp << endl;
      cout << " C'_p(T') = " << Cp1 << endl << endl;
      */


      density_vector[Stable+2] =  Cp1 * density_vector[Stable+1];
    }
    else {

      // Note: density_vector should be primitively constant from Stable+1 onwards! =)
      for (size_t j=0; j<=Stable+2; j++)
	density_vector[j] = Q_normal.Local_Density(p, t * (p^(2*j))) * (p^(Valuation(t,p) + 2*j));
      //      density_vector[Stable+2] =  density_vector[Stable+1];

    }


    // Find the smallest non-zero entry in density_vector
    mpq_class t_min;
    t_min = 0;
    for (size_t j=0; j < density_vector.size(); j++)
      if (t_min == 0)
	t_min = density_vector[j];
      else
	t_min = min(t_min, density_vector[j]);


    /*
      density_subvector := [density_vector[i] : i in [1..#density_vector] |  density_vector[i] gt 0];
      t_min, t_location := Minimum(density_subvector);
    */


    /*
    // DIAGNOSTIC -- Outputs the square class info (Part 2):
    // -----------------------------------------------------
    cout << " For the square-class " << t << "Z_" << p << " we have the local density vector ";
    cout << density_vector << endl;
    cout << " The minimal value was " << t_min << endl;
    // " at location " << t_location << " of the non-zero subvector ";
    */


    t_min.canonicalize();
    const_list[i] = t_min;


  }


  /*
    cout <<  "Q = " << Q << "  p = " << p << "  list of squareclasses = " << sqclasslist << endl;
  */


  //  return const_list, sqclasslist;
  return const_list;
}



/////////////////  The next 3 routines compute numerical things used in Lambda4test


//////////////////////////////////////////////////////////////////////////////
// Secret function to compute the small adjustment factors for Lambda4_test //
//////////////////////////////////////////////////////////////////////////////

mpq_class Matrix_mpz::_adjustment_secret(const mpz_class & m, const mpz_class & N, const mpz_class & chi_top) const {

  valarray<mpz_class> prime_list;
  prime_list.resize(PrimeDivisors(m).size());
  prime_list = PrimeDivisors(m);

  //  adj_list := [(p-1)/(p+1) : p in PrimeDivisors(m) | (Valuation(N, p) eq 0) and (KroneckerSymbol(chi_top, p) eq -1)];

  mpq_class new_factor;
  new_factor = 1;
  for(size_t i=0; i < prime_list.size(); i++)
    if ((Valuation(N, prime_list[i]) == 0) && (KroneckerSymbol(chi_top, prime_list[i]) == -1))
      new_factor = new_factor * mpq_class(prime_list[i] - 1, prime_list[i] + 1);

  /*
  // DIAGNOSTIC
  if (m == 971)
    cout << " DIAGNOSTIC:  adjustment factor = " << new_factor << endl;
  */

  return new_factor;
}



///////////////////////////////////////////////////////////////////////////////////////////////////
// Auxiliary function to compute the part of an integer m away from a fixed sequesce S of primes //
///////////////////////////////////////////////////////////////////////////////////////////////////

mpz_class Matrix_mpz::_integer_part_away_from_set(const mpz_class & m, const valarray<mpz_class> & S) const {

  mpz_class r;
  r = m;

  for(size_t i=0; i < S.size(); i++)
    r = r / (S[i] ^ Valuation(r, S[i]));

  /*
  // DIAGNOSTIC
  if (m == 971)
    cout << " DIAGNOSTIC:  integer part away from S = " << r << endl;
  */

  return r;
}





///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Compute an upper bound for the constant Lambda_4^(hat) numerically using the first few terms of the Eisenstein series //
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

mpq_class Matrix_mpz::_Lambda4_test(const Matrix_mpz & Q, const PowerSeries<mpq_class> & EE) const {

  mpz_class N;
  N = Q.QFLevel();

  cout << " The level is " << N << endl;

  // Find the anisotropic primes
  valarray <mpz_class> S;
  S.resize(Q.AnisotropicPrimes().size());
  S = Q.AnisotropicPrimes();

  cout << " The anisotropic primes are " << S << endl;

  mpz_class chi_top;
  chi_top = CoreDiscriminant(Q.Determinant());

  cout << " The character is given by " << chi_top << endl;


  /*
  // DIAGNOSTIC:
  cout << " The precision of EE is " << EE.Precision() << endl;
  */


  valarray<mpq_class> eisratios;
  eisratios.resize(EE.Precision() + 1);
  eisratios[0] = 496;  // This is just a dummy entry...
  for(size_t i=1; i<= EE.Precision(); i++) {
    //    cout << " Starting i = " << i << endl;
    eisratios[i] = EE[i] / (mpq_class(_integer_part_away_from_set(i,S)) * _adjustment_secret(i, N, chi_top));
  }


  /*
  // DIAGNOSTIC:
  cout << " The Eisenstein ratios are given by " << eisratios << endl;
  */


  // Return the minimum non-zero entry of Eisratios.
  mpq_class eis_min;
  size_t j_min;

  for(size_t j=1; j < eisratios.size(); j++)
    if (eisratios[j] > 0) {
      eis_min = eisratios[j];
      j_min = j;
    }

  for(size_t j=1; j < eisratios.size(); j++)
    if ((eisratios[j] > 0) && (eisratios[j] < eis_min)) {
      eis_min = eisratios[j];
      j_min = j;
    }

  /*
  // DIAGNOSTIC
  cout << " The minimum was found at the Eisenstein coefficient m = " << j_min << endl;
  cout << "   The original Eisenstein coefficient there is " << EE[j_min] << endl;
  cout << "   The eisenstein ratio there is " << eis_min << endl;
  */

  return eis_min;
}





////////////////////////////////////////////////////////////////////
// Revised Code for finding bounds for quaternary quadratic forms //
////////////////////////////////////////////////////////////////////

mpq_class Matrix_mpz::_new_C4_rational(const Matrix_mpz & Q) const {

  // Sanity Check: Need to make sure we're using a 4x4 matrix.
  assert(Q.NumRows() == 4);
  assert(Q.NumCols() == 4);

  mpz_class N;
  N = Q.QFLevel();
  //print "Level = ", N;

  valarray<mpz_class> N_primes;
  N_primes.resize(PrimeDivisors(N).size());
  N_primes = PrimeDivisors(N);

  //  N_primes_minima := [ -999/2 : i in [1..#N_primes]];


// For each prime and each primitive square class t(Z_p)^2 in Z_p
// we compute the local factors up to the minimal stable number T',
// and the constant C_p(T') (together with the adjustment at
// anisotropic primes).

  mpz_class p;
  valarray <mpq_class> N_primes_minima(N_primes.size());

  for (size_t k=0; k < N_primes.size(); k++){
    p = N_primes[k];
    cout << " p = " << p << endl;

    valarray<mpq_class> temp_C4(_PrimSqRepns(p).size());
    temp_C4 = _C4_squareclass_constants(Q, p, _PrimSqRepns(p));
    cout << "  and  C4 = " <<  temp_C4;

    N_primes_minima[k] = Minimum(temp_C4);
    cout << " which gives a minimum of " << N_primes_minima[k] << endl << endl;
  }


  cout << " Primes dividing the level: " <<  N_primes << endl;
  cout << " Associated minima at p: " << N_primes_minima << endl;

  mpz_class f;
  mpq_class d, big_denom, Bigfactor, Finalfactor;

  d = mpq_class(Q.Determinant(), 16);      // Note: Introduced the factor 16 here to compensate for the determinant being a factor of 16 too large!
  f = CoreDiscriminant(Q.Determinant() * 16);  // Note: CoreDiscriminant only takes an mpz_class, and we need to compensate for using the matrix 2Q

  big_denom = 1;
  for (size_t k=0; k < N_primes.size(); k++){
    p = N_primes[k];
    big_denom *= (1 - mpq_class(KroneckerSymbol(f,p), p*p));
  }

  // Note: Introduced the factor 16 here to compensate for the determinant being a factor of 16 too large!
  Bigfactor = mpq_class(f) * RationalSqrt(mpq_class(f,1) / d) / (QuadraticBernoulliNumber(2, f.get_si()) * big_denom);


  /*
  // DIAGNOSTIC
  cout << " Q = " << endl << Q << endl;
  cout << " Level = " << N << endl;
  cout << " Level primes = " << N_primes << endl;
  cout << " 16 * d = " << Q.Determinant() << endl;
  cout << " d = " << d << endl;
  cout << " f = " << f << endl;
  cout << " big_denom = " << big_denom << endl;
  cout << " Bigfactor = " << Bigfactor << endl;
  */


  Finalfactor = Bigfactor;
  for (size_t k=0; k < N_primes_minima.size(); k++)
    Finalfactor *= N_primes_minima[k];

  Finalfactor.canonicalize();

  return Finalfactor;
}


// ==========================================================================================================


mpq_class Matrix_mpz::GetEisensteinLowerBound(const string & Eis_Bound_Dir) const {


  // Create a place for the bound
  mpq_class Eis_bound;


  // Make the path the the Eisenstein series file
  string EisBoundFilename;
  EisBoundFilename = Eis_Bound_Dir + "Eis_Bound" + "__" + (*this).QF_String() + ".txt";


  // Try to read from the file.
  ifstream eisboundinfile;
  eisboundinfile.open(EisBoundFilename.c_str());
  if (eisboundinfile.is_open() == true) {
    eisboundinfile >> Eis_bound;
    eisboundinfile.close();
  }
  // If this fails, compute it and write the answer!
  else {
    eisboundinfile.close();

    // Compute the Eisenstein lower bound
    Eis_bound = _new_C4_rational(*this);

    // Write the result
    ofstream eisboundoutfile;
    eisboundoutfile.open(EisBoundFilename.c_str());
    eisboundoutfile << Eis_bound;
    eisboundinfile.close();
  }


  // Return the result
  return Eis_bound;

}




mpq_class Matrix_mpz::NumericalEisensteinLowerBound(const string & Eis_Dir, const unsigned long precision) const {

  // Get the Eisenstien series (with the precision specified in _GetEisData())
  PowerSeries<mpq_class> Eis_series;
  Eis_series = GetMagmaEisSeries(Eis_Dir, precision);
  //  Eis_series = _GetEisData(Eis_Dir, precision);

  // Compute the numerical estimate for the Eisenstein coefficient lower bound
  return _Lambda4_test(*this, Eis_series);    // We should also specify the number of terms to use...

}



// ========================= 4 private routines used on reading in the Eisenstein series ==============================

PowerSeries<mpq_class> Matrix_mpz::_GetEisData(const string & Eis_Dir, const unsigned long precision) const {

  // Sanity check that we're only using 4x4 matrices... (because of the filename convention...)
  assert((*this).NumRows() == 4);
  assert((*this).NumCols() == 4);


  // Create a power series with the desired precision
  PowerSeries<mpq_class> Eis_Series(precision);




  // Make the path the the Eisenstein series file
  string EisFilename;
  EisFilename = Eis_Dir + "Eis__" + (*this).QF_String() + "__" + MakeString(precision) + ".txt";


  /*
  // DIAGNOSTIC
  cout << " Using the Eisenstein series filename: " << EisFilename << endl;
  cout << " Here's the QF_String() inside: " << (*this).QF_String() << endl;
  */


  /*
  // Set some constants
  const char EIS_DIR[] = "/home/postdoc/jonhanke/290_Project/EIS_DATA/New_Eis_Data__10-4-04/";


  // Make the path the the Eisenstein series file
  char SeriesFilename[200];
  if (Eis_filename == "") {

    // Make the filename "EIS_DIR[]/Eis__a_b_c_d_e_f_g_h_i_j__precision.txt" for the Eisenstein series file
    sprintf(SeriesFilename, "%sEis__%d_%d_%d_%d_%d_%d_%d_%d_%d_%d__%d.txt",
	    EIS_DIR,
	    (*this)(1,1).get_si() / 2, (*this)(1,2).get_si(), (*this)(1,3).get_si(), (*this)(1,4).get_si(),
	    (*this)(2,2).get_si() / 2, (*this)(2,3).get_si(), (*this)(2,4).get_si(),
	    (*this)(3,3).get_si() / 2, (*this)(3,4).get_si(),
	    (*this)(4,4).get_si() / 2,
	    );
  }
  else {

    // Make the complete filename with the given filename
    sprintf(SeriesFilename, "%s%s", EIS_DIR, Eis_filename);
  }
  */


  // TO DO:
  // Check if any of them exist, and choose the smallest necessary precision


  // If the file doesn't exist, then generate it.
  ifstream seriesfile;
  seriesfile.open(EisFilename.c_str(), ios::in);
  if (seriesfile.is_open() == true) {

    // No need to compute anything, so forget it.
    seriesfile.close();
  }
  else {

    // Compute the Eisenstein Series (by calling Magma), since no file was found...
    _ComputeEisData(EisFilename.c_str(), precision);
  }


  // Read the series (which we know now exists!)
  Eis_Series.ReadSeries(EisFilename.c_str());


  /*
  // DIAGNOSTIC
  cout << " Finished reading the series: " << Eis_Series << endl;
  */


  // Return the result
  return Eis_Series;

}


void Matrix_mpz::_ComputeEisData(const string & EisFilename, const unsigned long Eis_precision) const {

  // Write the file which we send to Magma:
  // --------------------------------------
  char FormFilename[100] = "current_form.magma";

  // Try to open the form file
  ofstream formfile;
  formfile.open(FormFilename, ios::out);

  // Abort if we fail... =(
  if (! formfile.is_open())
    { cout << "_ComputeEisData() Error: Error opening the form file for writing... =("; exit (1); }

  // Write the form
  formfile << "new_Q := [ ";
  for(long i=1; i<=4; i++)
    for(long j=1; j<=4; j++) {
      formfile << (*this)(i,j);
      if ((i != 4) || (j != 4))
	formfile << " ,";
    }
  formfile << " ];" << endl;

  // Write the precision
  formfile << "new_precision :=  " << Eis_precision << ";" << endl;

  // Close the form file
  formfile.close();


  // Send the computation to Magma:
  // ------------------------------

  // Copy the form and precision information
  system("scp current_form.magma jonhanke@math.princeton.edu:__290_Temp_Computations/current_form.magma");

  // Make the Eisenstein series
  system("ssh jonhanke@math.princeton.edu 'magma __290_Make_Eis.magma'");

  // Copy the Eisenstein series back
  char copyback_command[300];
  sprintf(copyback_command, "scp jonhanke@math.princeton.edu:__290_Temp_Computations/Current_Eis.txt %s", EisFilename.c_str());
  system(copyback_command);

  // Erase the temporary files
  system("rm -f current_form.magma");
  system("ssh jonhanke@math.princeton.edu 'rm -f __290_Temp_Computations/current_form.magma'");
  system("ssh jonhanke@math.princeton.edu 'rm -f __290_Temp_Computations/Current_Eis.txt'");


  // ----------------------------------------------------------------
  // Note: The instructions for setting up a remote key pair is at:
  // http://www.ale.org/archive/ale/ale-2002-07/msg01029.html
  // ----------------------------------------------------------------
}


































// Tries to find the Eisenstein series in a file, and if none is found, asks Magma to compute one!

PowerSeries<mpq_class> Matrix_mpz::GetMagmaEisSeries(const string & Eis_Dir, const unsigned long precision) const {

  // Sanity check that we're only using 4x4 matrices... (because of the filename convention...)
  assert((*this).NumRows() == 4);
  assert((*this).NumCols() == 4);


  // Create a power series with the desired precision
  PowerSeries<mpq_class> Eis_Series(precision);


  // Make the path to the Eisenstein series file
  string EisFilename;
  EisFilename = Eis_Dir + "Eis__" + (*this).QF_String() + "__" + MakeString(precision) + ".txt";


  // If the file doesn't exist, then generate it.
  ifstream seriesfile;
  seriesfile.open(EisFilename.c_str());
  if ((Eis_Dir != "") && (seriesfile.is_open() == true)) {

    // No need to compute anything, so forget it.
    seriesfile.close();

    // Read the series (which we know now exists!)
    Eis_Series.ReadSeries(EisFilename.c_str());

  }
  else {

    // Make the variable input string for the Magma script
    ostringstream var_stringstream;

    // Input: Write the size of the form
    var_stringstream << "new_n := " << NumRows() << ";" << endl;

    // Input: Write the form as a vector
    var_stringstream << "new_Q := [ ";
    for(long i=1; i<=NumRows(); i++)
      for(long j=1; j<=NumRows(); j++) {
	var_stringstream << (*this)(i,j);
	if ((i != NumRows()) || (j != NumRows()))
	  var_stringstream << " ,";
      }
    var_stringstream << " ];" << endl;

    // Input: Write the precision
    var_stringstream << "new_precision :=  " << precision << ";" << endl;


    // Compute the Eisenstein Series (by calling Magma), since no file was found...
    string tmp_filename;
    tmp_filename = _GetMagmaComputation("Eis","__Make_Magma_Eisenstein_Series.m", var_stringstream.str());

    // Read the answer
    Eis_Series.ReadSeries(tmp_filename.c_str());


    // Deal with the temporary file
    if (Eis_Dir == "") {

      // Erase the temporary file
      string command_string;
      command_string = "rm -f " + tmp_filename;
      system(command_string.c_str());
    }
    else {

      // Move the temporary file to its new home =)
      string command_string;
      command_string = "mv " + tmp_filename + " " + EisFilename;
      system(command_string.c_str());
    }


  }


  /*
  // DIAGNOSTIC
  cout << " Finished reading the series: " << Eis_Series << endl;
  */


  // Return the result
  return Eis_Series;

}




// Tries to find the theta series in a file, and if none is found, asks Magma to compute one!

PowerSeries<mpz_class> Matrix_mpz::GetMagmaThetaSeries(const string & Theta_Dir, const unsigned long precision) const {

  // This sanity check seems unnecessary now that we wrote QF_String()
  /*
  // Sanity check that we're only using 4x4 matrices... (because of the filename convention...)
  assert((*this).NumRows() == 4);
  assert((*this).NumCols() == 4);
  */

  // Create a power series with the desired precision
  PowerSeries<mpz_class> Theta_Series(precision);


  // Make the path to the theta series file
  string ThetaFilename;
  ThetaFilename = Theta_Dir + "Theta__" + (*this).QF_String() + "__" + MakeString(precision) + ".txt";


  // If the file doesn't exist, then generate it.
  ifstream seriesfile;
  seriesfile.open(ThetaFilename.c_str());
  if ((Theta_Dir != "") && (seriesfile.is_open() == true)) {

    // No need to compute anything, so forget it.
    seriesfile.close();

    // Read the series (which we know now exists!)
    Theta_Series.ReadSeries(ThetaFilename.c_str());

  }
  else {

    // Make the variable input string for the Magma script
    ostringstream var_stringstream;

    // Input: Write the size of the form
    var_stringstream << "new_n := " << NumRows() << ";" << endl;

    // Input: Write the form as a vector
    var_stringstream << "new_Q := [ ";
    for(long i=1; i<=NumRows(); i++)
      for(long j=1; j<=NumRows(); j++) {
	var_stringstream << (*this)(i,j);
	if ((i != NumRows()) || (j != NumRows()))
	  var_stringstream << " ,";
      }
    var_stringstream << " ];" << endl;

    // Input: Write the precision
    var_stringstream << "new_precision :=  " << precision << ";" << endl;


    // Compute the Eisenstein Series (by calling Magma), since no file was found...
    string tmp_filename;
    tmp_filename = _GetMagmaComputation("Theta","__Make_Magma_Theta_Series.m", var_stringstream.str());

    // Read the answer
    Theta_Series.ReadSeries(tmp_filename.c_str());


    // Deal with the temporary file
    if (Theta_Dir == "") {

      // Erase the temporary file
      string command_string;
      command_string = "rm -f " + tmp_filename;
      system(command_string.c_str());
    }
    else {

      // Move the temporary file to its new home =)
      string command_string;
      command_string = "mv " + tmp_filename + " " + ThetaFilename;
      system(command_string.c_str());
    }

  }


  /*
  // DIAGNOSTIC
  cout << " Finished reading the series: " << Theta_Series << endl;
  */


  // Return the result
  return Theta_Series;

}








// Calls Magma to do a quadratic form computation, returning a filename to it.

string Matrix_mpz::_GetMagmaComputation(const string & Filename_fragment, const string & Client_Scriptname, const string & var_string) const {

  // TO DO: The remote directory is assumed to exist, but we should really check/generate it!


  // An often used variable...
  string command_string;


  // Set the temporary directory filename
  string temp_dir = GetAbsolutePath("~/QF_Project_Data/__Temp_MAGMA_Computations/");


  // Check for the existence of a temporary directory, and make it if necessary
  if (DirectoryExists(temp_dir) == false) {
    command_string = "mkdir " + temp_dir;
    system(command_string.c_str());
  }


  // Make the Magma client filename(s)
  string magma_base_filename, local_magma_input_filename, local_magma_output_filename;
  magma_base_filename = "Magma_" + Filename_fragment + "_client__" + QF_String() + "__" + MakeUniqueFilename() + ".magma";
  local_magma_input_filename = temp_dir + magma_base_filename;
  local_magma_output_filename = temp_dir + "Output__" + magma_base_filename;

  string remote_magma_input_filename, remote_magma_output_filename;
  string remote_dir = "__290_Temp_Computations/";                               // The directory where it all happens on the remote host
  remote_magma_input_filename = remote_dir + magma_base_filename;            // The client filename on the remote host
  remote_magma_output_filename = remote_dir + "Output__" + magma_base_filename;    // The output filename on the remote host


  // Write the file which we send to Magma:
  // --------------------------------------

  // Try to open the form file
  ofstream formfile;
  formfile.open(local_magma_input_filename.c_str());

  // Abort if we fail... =(
  if (! formfile.is_open())
    { cout << "_ComputeEisData() Error: Error opening the form file for writing... =("; exit (1); }

  // Write the variable string passed into the routine
  formfile << var_string << endl;

  // Write the output filename
  formfile << "output_filename := \"" << remote_magma_output_filename << "\";" << endl;

  // Close the form file
  formfile.close();


  // Append the rest of the Magma client script:
  // -------------------------------------------
  string client_command;
  client_command = "cat " + temp_dir + "Magma_Client_Scripts/" + Client_Scriptname + " >> " + local_magma_input_filename;
  system(client_command.c_str());


  // Send/Receive the computation to/from Magma:
  // -------------------------------------------

  // Copy the form and precision information
  command_string = "scp " + local_magma_input_filename + " jonhanke@math.princeton.edu:" + remote_magma_input_filename;
  system(command_string.c_str());
  //  system("scp current_form.magma jonhanke@math.princeton.edu:__290_Temp_Computations/current_form.magma");

  // Make the Eisenstein series
  command_string = "ssh jonhanke@math.princeton.edu 'magma " + remote_magma_input_filename + "'";
  system(command_string.c_str());
  //  system("ssh jonhanke@math.princeton.edu 'magma __290_Make_Eis.magma'");

  // Copy the Eisenstein series back
  command_string = "scp jonhanke@math.princeton.edu:" + remote_magma_output_filename + " " + local_magma_output_filename;
  system(command_string.c_str());
  //  sprintf(copyback_command, "scp jonhanke@math.princeton.edu:__290_Temp_Computations/Current_Eis.txt %s", EisFilename.c_str());



  // Erase the temporary files:
  // --------------------------
  command_string = "rm -f " + local_magma_input_filename;
  system(command_string.c_str());
  //  system("rm -f current_form.magma");

  command_string = "ssh jonhanke@math.princeton.edu 'rm -f " + remote_magma_input_filename + "'";
  system(command_string.c_str());
  //  system("ssh jonhanke@math.princeton.edu 'rm -f __290_Temp_Computations/current_form.magma'");

  command_string = "ssh jonhanke@math.princeton.edu 'rm -f " + remote_magma_output_filename + "'";
  system(command_string.c_str());
  //  system("ssh jonhanke@math.princeton.edu 'rm -f __290_Temp_Computations/Current_Eis.txt'");



  // ----------------------------------------------------------------
  // Note: The instructions for setting up a remote key pair is at:
  // http://www.ale.org/archive/ale/ale-2002-07/msg01029.html
  // ----------------------------------------------------------------



  // Return the result
  return local_magma_output_filename;

}


