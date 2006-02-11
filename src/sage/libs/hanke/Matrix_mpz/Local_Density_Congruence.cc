#include "Matrix_mpz_header.h"
#include "../GMP_class_extras/mpz_class_extras.h"
#include "../GMP_class_extras/vectors.h"


///////////////////////////////////////////////////////////////////////////
// Creates a new index vector which points to the extracted rows/columns
// in the extracted matrix where the ...?
//
// (Note: This is probably not very efficient,
// and could be improved by using vectors,
// but we're using valarrays for now...)
////////////////////////////////////////////////////////////////////////////

valarray<int> Matrix_mpz::ReindexVectorFromExtraction(const valarray<int> & Original, const valarray<int> & Extracted) const
{
  valarray<int> Reindex;
  Reindex.resize(Original.size());
  Reindex = 0;

  // Replace all Original indices entries with the position of the matching Extraced index
  unsigned long i, j, ind;
  ind = 0;
  for (i=0; i<Original.size(); i++)
    for (j=0; j<Extracted.size(); j++)
      if (Original[i] == Extracted[j]) {
	Reindex[ind] = j+1;
	ind++;
      }

  // Copy these to a new vecarray of the appropriate length -- Since Reindex may be too big
  valarray<int> Final(ind);
  for (i=0; i<ind; i++)
      Final[i] = Reindex[i];

	 return Final;
}









///////////////////////////////////////////////////////////////////////////////
// Computes the Gauss sum for the number of local solutions Q(x) = m mod p>2 //
///////////////////////////////////////////////////////////////////////////////

mpz_class Matrix_mpz::GaussLocal(int n1, const mpz_class & p, const mpz_class & m, const mpz_class & Qdet) const
{
  mpz_class count;
  mpz_class neg1 = -1;

  // Sanity Check
  assert (n1 >= 0);
  unsigned long n = n1;

  /*
  cout << " n = " << n << endl;
  cout << " neg1 = " << neg1 << endl;
  cout << " p^(n-1) = " << (p^(n-1)) << endl;
  */

  if ((p != 2) && (n >= 1))  // Check that p is an odd prime and n >= 1  --  To Do: Need to check if p is a prime...
    {
      if (m % p == 0)
	if (n % 2 != 0)
	  count = (p^(n-1));
	else
	  count = (p^(n-1)) + (p-1) * (p^((n-2)/2)) * KroneckerSymbol(((neg1^(n/2)) * Qdet) % p, p);
      else
	if (n % 2 != 0)
	  count = (p^(n-1)) + (p^((n-1)/2)) * KroneckerSymbol(((neg1^((n-1)/2)) * Qdet * m) % p, p);
	else
	  count = (p^(n-1)) - (p^((n-2)/2)) * KroneckerSymbol(((neg1^(n/2)) * Qdet) % p, p);
    }
  else
    cout << "\n Error in GaussLocal: Either p is not prime or n<1 \n";

  return(count);
}




//////////////////////////////////////////////////////////////////
// Finds the Good-type local density of Q representing m at p.  //
// (Assuming that p > 2 and Q is given in local diagonal form.) //
//////////////////////////////////////////////////////////////////

mpq_class Matrix_mpz::Local_Good_Density_Congruence_Odd(const mpz_class & p, const mpz_class & m,
					    const valarray<int> & Zvec, const valarray<int> & NZvec) const
{

     // To Do: Check to see if Q_int is diagonal (only ok for p>2)

  int n;
  n = (*this).NumRows();

  Matrix_mpz Qtrim;
  int i;

  // Assuming Q is diagonal, trim it to ensure it's non-degenerate mod p
  valarray<int> Qtrimvec_temp(n);
  for (i=0; i<n; i++)
    Qtrimvec_temp[i] = 0;

  int ptr = 0;

  for (i=1; i<=n; i++)
    if (((*this)(i,i) % p) != 0) {
      Qtrimvec_temp[ptr] = i;
      ptr++;
    }

  /*
  cout << " Q = " << endl << Q << endl;
  cout << " n = " << n << "  and  Qtrimvec_temp  = " << Qtrimvec_temp << endl;
  */

  valarray<int> Qtrimvec(VectorTrim(Qtrimvec_temp));


  // DEBUGGING:  Check the matrix is primitive --
  //   (WARNING: We may be passed imprimitive form from the BI-reduction...)  <==  This should be fixed now. =)
  assert(Qtrimvec.size() > 0);


  /*
  cout << " Stage 1" << endl;
  cout << " Qtrimvec_temp = " << Qtrimvec_temp << endl;
  cout << " Qtrimvec = " << Qtrimvec << endl;
  cout << " Q = \n" << Q << endl;
  */

  Qtrim = (*this).ExtractSquareSubmatrixOrdered(Qtrimvec);


  // cout << "We're here now..." << endl;


  // Construct the big and small matrices:
  // -------------------------------------
  unsigned int smallfreelength, bigfreelength;
  // int j, ind;


  // Construct new congruence condition indices for the trimmed matrix
  valarray<int> trimZvec(ReindexVectorFromExtraction(Zvec, Qtrimvec));

  /*
  cout << "\n Zvec = " << Zvec << "\n";
  cout << " NZvec = " << NZvec << "\n";
  cout << " Qtrim = " << Qtrim << "\n";
  cout << " Qtrimvec = " << Qtrimvec << "\n";
  cout << " trimZvec = " << trimZvec << "\n";
  cout << " Vector indexing Qtrim = " << MakeVector(Qtrimvec.size(), 1, 1) << "\n";
  cout << " and its complement by trimZvec  = " << \
    VectorComplement(MakeVector(Qtrimvec.size(), 1, 1), trimZvec) << "\n";
  */


  // cout << "Getting closer..." << endl;


  // Make the big trim vector (all trim indices not somewhere in Zvec)
  //   and the big free vector (all non-trim indices not somewhere in Zvec)
  //   [To Do: This could probably be faster if we assumed Zvec was ordered.]
  valarray<int> bigtrimvec(VectorComplement(MakeVector(Qtrimvec.size(), 1, 1), trimZvec));
  bigfreelength = (n - Qtrimvec.size()) - (Zvec.size() - trimZvec.size());




  // Make the small vector from the big one (all indices not somewhere in Zvec or NZvec)
  valarray<int> new_vec(VectorUnion(Zvec, NZvec)),
    trim_new_vec(ReindexVectorFromExtraction(new_vec, Qtrimvec));

  valarray<int> smalltrimvec(VectorComplement(MakeVector(Qtrimvec.size(), 1, 1), new_vec));
  smallfreelength = (n - Qtrimvec.size()) - (new_vec.size() - trim_new_vec.size());




  // Make the big and small matrices
  Matrix_mpz bigmatrix, smallmatrix;
  bigmatrix = Qtrim.ExtractSquareSubmatrixOrdered(bigtrimvec);
  smallmatrix = Qtrim.ExtractSquareSubmatrixOrdered(smalltrimvec);


  /*
  cout << "\n Q is : " << Q << "\n";
  cout << " m is : " << m << "\n";
  cout << " p is : " << p << "\n";
  cout << " Qtrim is : " << Qtrim << "\n";
  cout << " bigtrimvec is : " << bigtrimvec << "\n";
  cout << " bigfreelength is : " << bigfreelength << "\n";
  cout << " smalltrimvec is : " << smalltrimvec << "\n";
  cout << " smallfreelength is : " << smallfreelength << "\n";
  cout << " bigmatrix is : \n" << bigmatrix << "\n";
  cout << " smallmatrix is : \n" << smallmatrix << "\n";
  */


  // cout << "And closer still..." << endl;


  // Number of representations
  mpz_class big_factor, small_factor;
  mpz_class total;
  mpq_class good_density;

  return mpq_class(0);

  // Adjust the local solutions to count only the good-type ones
  if (bigtrimvec.size() == 0)    // Check if we have the zero-matrix...

    big_factor = 0;

  else

    if (m%p != 0)
      big_factor = (p^bigfreelength)
	* GaussLocal(bigmatrix.NumRows(), p, m, bigmatrix.Determinant());
    else
      big_factor = (p^bigfreelength)
	* (GaussLocal(bigmatrix.NumRows(), p, m, bigmatrix.Determinant()) - 1);


  // Similarly for the smallmatrix if it exists
  if (smalltrimvec.size() == 0)    // Check if we have the zero-matrix...

    small_factor = 0;

  else

    if (NZvec.size() > 0)
      if (m%p != 0)
	small_factor = (p^smallfreelength)
	  * GaussLocal(smallmatrix.NumRows(), p, m, smallmatrix.Determinant());
      else
	small_factor = (p^smallfreelength)
	  * (GaussLocal(smallmatrix.NumRows(), p, m, smallmatrix.Determinant()) - 1);

    else
      small_factor = 0;

  total = big_factor - small_factor;
  good_density = mpq_class(total, mpz_power(p, n-1));
  //  good_density = mpq_class(to_double(total), to_double(power(p,(n-1)));
  good_density.canonicalize();


  return(good_density);
}



//////////////////////////////////////////////////////////////////
// Finds the Good-type local density of Q representing m at p.  //
// (Assuming that p > 2 and Q is given in local diagonal form.) //
//////////////////////////////////////////////////////////////////

mpq_class Matrix_mpz::Local_Good_Density_Congruence_Even(const mpz_class & p, const mpz_class & m,
					     const valarray<int> & Zvec, const valarray<int> & NZvec) const
{

  //  cout << " Break 0" << endl;

  /*
  cout << "\n Q is : " << Q << "\n";
  */

  int n = (*this).NumRows();

  Matrix_mpz Qtrim;
  int i;

  // Assuming Q is diagonal, trim it to ensure it's non-degenerate mod 8
  valarray<int> Qtrimvec(n);
  for (i=0; i<n; i++)
    Qtrimvec[i] = 0;

  int Qtrimvec_ptr = 0;

  //    cout << " Break 0.1" << endl;

  // Find the indices of the non-zero blocks mod 8
    for (i=1; i<=n; i++) {

      /*
      cout << " i = " << i << endl;
      cout << " n = " << n << endl;
      cout << " Qtrimvec_ptr = " << Qtrimvec_ptr << endl;
      cout << " Qtrimvec = " << Qtrimvec << endl;
      */

      bool nz_flag = false;
      if  (((*this)(i,i) % 8) != 0)
	nz_flag = true;
      else {
	//	cout << " here 1" << endl;
	if ((i==1) && (((*this)(i,i+1) % 8) != 0))
	  nz_flag = true;
	else {
	  //	  cout << " here 2" << endl;
	  if ((i==n) && (((*this)(i-1,i) % 8) != 0))
	    nz_flag = true;
	  else {
	    //	    cout << " here 3" << endl;
	    if ( (i > 1)  &&  (i < n)  &&  ((((*this)(i,i+1) % 8) != 0) || (((*this)(i-1,i) % 8) != 0)) )
	      nz_flag = true;
	  }
	}
      }


      if (nz_flag == true) {
	Qtrimvec[Qtrimvec_ptr] = i;
	Qtrimvec_ptr++;
      }

      // DEBUGGING: Check the pointer is always in range...
      assert(Qtrimvec_ptr <= Qtrimvec.size());
    }






    //  cout << " Break 1" << endl;



  valarray<int> QtrimvecNew(VectorTrim(Qtrimvec));


  // DEBUGGING: Quick tests to make sure the form isn't zero mod 8
  assert(Qtrimvec.size() > 0);
  assert(QtrimvecNew.size() > 0);



  Qtrim = (*this).ExtractSquareSubmatrixOrdered(QtrimvecNew);


  // DEBUGGING: Quick tests to make sure the extracted form isn't zero mod 8
  assert(Qtrim.NumRows() > 0);


  //  cout << " Break 2" << endl;



  // Construct the big and small matrices:
  // -------------------------------------

  // Construct new congruence condition indices for the trimmed matrix
  valarray<int> trimZvec(ReindexVectorFromExtraction(Zvec, QtrimvecNew)),
    trimNZvec(ReindexVectorFromExtraction(NZvec, QtrimvecNew));


  // Make the trimmed congruence vector
  valarray<int> new_vec(VectorUnion(Zvec, NZvec)),
    trim_new_vec(ReindexVectorFromExtraction(new_vec, QtrimvecNew));




  /*
  cout << "\n Q is : " << Q << "\n";
  cout << " m is : " << m << "\n";
  cout << " Qtrim is : " << Qtrim << "\n";
  cout << " Qtrimvec is : " << Qtrimvec << "\n";
  cout << " trimZvec is : " << trimZvec << "\n";
  cout << " trimNZvec is : " << trimNZvec << "\n";
  */

  // DEBUGGING: Check that partlyfreenum is in range...
  assert(new_vec.size() >= trim_new_vec.size());



  // Compute the number of different free components
  int partlyfreenum, veryfreenum;

  partlyfreenum = new_vec.size() - trim_new_vec.size();
  veryfreenum = (n - QtrimvecNew.size()) - partlyfreenum;


  // In the free part, each component with a congruence condition contrubite a factor of 4,
  // while components with no congruence conditions contribute a factor of 8.
  mpz_class total;
  mpq_class good_density;

  total = mpz_power(4, partlyfreenum) * mpz_power(8, veryfreenum)
    * Qtrim.CountLocalGoodType(2, 3, m, trimZvec, trimNZvec);

  good_density = mpq_class(total, mpz_power(8, n-1));
  good_density.canonicalize();


  /*
  cout << " partlyfreenum = " << partlyfreenum << "\n";
  cout << " veryfreenum = " << veryfreenum << "\n";
  cout << " CountLocalGoodType(Qtrim, 2, 3, m, trimZvec, trimNZvec) = " << CountLocalGoodType(Qtrim, 2, 3, m, trimZvec, trimNZvec)  << "\n";
  cout << " total = " << total << "\n";
  cout << " denominator = " << (mpz_class(8)^(n-1)) << "\n";
  cout << " Good Density = " << good_density << "\n";
  */

  return(good_density);
}



     // Note: Assumes all forms are primitive


//////////////////////////////////////////////////////////////////
// Finds the Good-type local density of Q representing m at p.  //
// (Front end routine for parity specific routines for p.)      //
//////////////////////////////////////////////////////////////////

mpq_class Matrix_mpz::Local_Good_Density_Congruence(const mpz_class & p, const mpz_class & m,
					const valarray<int> & Zvec, const valarray<int> & NZvec) const
{

  /*
  // For debugging purposes:
  cout << " In Local_Good_Density_Congruence with " << endl;
  cout << " Q is: " << endl << Q << endl;
  cout << " p = " << p << endl;
  cout << " m = " << m << endl;
  cout << " Zvec = " << Zvec << endl;
  cout << " NZvec = " << NZvec << endl;
  */


  // Check that Q is in local normal form -- should replace this with a diagonalization check?
  //   (it often may not be since the reduction procedure
  //   often mixes up the order of the valuations...)
  /*
  if  (Q != LocalNormal(Q, p))
    cout << "Warning in Local_Good_Density_Congruence: Q is not in local normal form! \n";
  */


  // Check that the congruence conditions don't overlap
  if (CheckVectorIntersection(Zvec, NZvec) == true)
    return (mpq_class(0));

  //  cout << "We're here!" << endl;


  // Decide which routine to use to compute the Good-type density
  if (p>2)
    return (*this).Local_Good_Density_Congruence_Odd(p, m, Zvec, NZvec);

  if (p==2) {
    // cout << "\n Using the (p=2) Local_Good_Density_Even routine! \n";
    return (*this).Local_Good_Density_Congruence_Even(p, m, Zvec, NZvec);
  }

  cout <<  "\n Error in Local_Good_Density: The 'prime' p = " << p << " is < 2. \n";
  return -1;
}






///////////////////////////////////////////////////////////
// Finds the Zero-type local density of Q representing   //
// m at p, allowing certain congruence conditions mod p. //
///////////////////////////////////////////////////////////

mpq_class Matrix_mpz::Local_Zero_Density_Congruence(const mpz_class & p, const mpz_class & m,
					const valarray<int> & Zvec, const valarray<int> & NZvec) const
{

  /*
  // For debugging purposes:
  cout << " In Local_Zero_Density_Congruence with " << endl;
  cout << " Q is: " << endl << Q << endl;
  cout << " p = " << p << endl;
  cout << " m = " << m << endl;
  cout << " Zvec = " << Zvec << endl;
  cout << " NZvec = " << NZvec << endl;
  */


  mpz_class p2;
  p2=p*p;

  if ((m % (p2) != 0) || (NZvec.size() > 0))
    return 0;
  else {
    valarray<int> EmptyVec;
    // Need 2 steps since we can't take negative powers... =|
    if ((*this).NumRows() > 2)
      return mpq_class(1, mpz_power(p, (*this).NumRows() - 2)) \
	* (*this).Local_Density_Congruence(p, m / (p2), EmptyVec, EmptyVec);
    else
      return mpq_class(mpz_power(p, 2 - (*this).NumRows())) \
	* (*this).Local_Density_Congruence(p, m / (p2), EmptyVec, EmptyVec);
  }
}




///////////////////////////////////////////////////////////////////
// Finds the Bad-type I local density of Q representing m at p.  //
// (Assuming that p > 2 and Q is given in local diagonal form.)  //
///////////////////////////////////////////////////////////////////

mpq_class Matrix_mpz::Local_BadI_Density_Congruence(const mpz_class & p, const mpz_class & m,
					const valarray<int> & Zvec, const valarray<int> & NZvec) const
{

  /*
  // For debugging purposes:
  cout << " In Local_BadI_Density_Congruence with " << endl;
  cout << " Q is: " << endl << Q << endl;
  cout << " p = " << p << endl;
  cout << " m = " << m << endl;
  cout << " Zvec = " << Zvec << endl;
  cout << " NZvec = " << NZvec << endl;
  */


  int n = (*this).NumRows();
  assert(n > 0);
  valarray<int> S0(n);
  bool S1_empty_flag = true;    // This is used to check if we should be computing BI solutions at all!
                                // (We should really to this earlier, but S1 must be non-zero to proceed.)

  // Define the indexing sets S_i
  int S0_ptr = 0;

  int i, j;

  mpz_class val;

  for(i=1; i<=n; i++) {

    // Compute the valuation of each index, allowing for off-diagonal terms
    if ((*this)(i,i) == 0)
      {
	if (i==1)
	  {
	    val = Valuation((*this)(i,i+1), p);  // Look at the term to the right
	  }
	else
	  {
	    if (i==n)
	      {
		val = Valuation((*this)(i-1,i), p);  // Look at the term above
	      }
	    else
	      {
		val = Valuation((*this)(i,i+1) + (*this)(i-1,i), p);  // Finds the valuation of the off-diagonal term since only one isn't zero
	      }
	  }
      }
    else
      {
	val = Valuation((*this)(i,i), p);
      }

    if (val == 0) {
      assert(S0_ptr < S0.size());  // Check that the pointer is in range...
      S0[S0_ptr] = i;
      S0_ptr++;
    }

    if (val == 1)
      S1_empty_flag = false;   // Need to have a non-empty S1 set to proceed with BI reduction...

  }

  // Check that S1 is non-empty to proceed, otherwise return no solutions.
  if (S1_empty_flag == true)
    return 0;


  // Check that the form is primitive...
  if (S0.size() <= 0)
    {
      cerr << " Error in Local_BadI_Density_Congruence: The form is not primitive!" << endl;
      cerr << " Using Q = " << (*this) << endl;
      cerr << " and p = " << p << endl;
      exit(0);
    }


  valarray<int> S0New(VectorTrim(S0));
  //  cerr << S0New.size() << endl;
  assert(S0New.size() > 0);     // Check that the S0 vector is is non-empty...





  // Note: The following lines assume that NZvec is an ordered vector.  Should check this here...


  /*
  cout << " m = " << m << "   p = " << p << "\n";
  cout << " S0 = " << S0 << "   NZvec = " << NZvec << "   IsDisjoint = " << IsDisjointOrdered(S0, NZvec)  << "\n";
  cout << " S0.size = " << S0.size() << "\n";
  */

  // Make the form Qnew for the reduction procedure
  if ((m % p == 0) && IsDisjointOrdered(S0New, NZvec) && (S0New.size() != n)) {
    //cout << "Test 1.1 \n";
    Matrix_mpz Qnew;
    Qnew = (*this);
    j=1;
    for(i=1; i<=n; i++)
      {
	// assert(j-1 >= 0);
	// assert(j-1 < S0New.size());
	// Short circuit && necessary here:
	if ((j - 1 < S0New.size()) && (i == S0New[j-1])) {                 // i is in S0
	  j++;
	  Qnew(i,i) = p * Qnew(i,i);
	  if ((p == 2) && (i < n)) {
	    Qnew(i+1,i) = p * Qnew(i+1,i);
	    Qnew(i,i+1) = p * Qnew(i,i+1);
	  }
	}
	else {                            // i not in S0;  Since S0 is increasing, we know i < S0(j)
	  Qnew(i,i) = Qnew(i,i) / p;

	  //	cout << "  dividing in row " << i << endl;

	  if ((p == 2) && (i < n)) {
	    Qnew(i+1,i) = Qnew(i+1,i) / p;
	    Qnew(i,i+1) = Qnew(i,i+1) / p;
	  }
	}
      }

    /*
    cout << "\n\n Check of Bad-type I reduction: \n";
    cout << " Q is " << Q << "\n";
    cout << " Qnew is " << Qnew << "\n";
    cout << " p = " << p << endl;
    cout << " m / p = " << (m/p) << endl;
    cout << " VectorComplement(Zvec, S0) is " << VectorComplement(Zvec, S0) << endl;
    cout << " NZvec " << NZvec << endl;
    */


    // Do the reduction
    // (Need 2 steps since we can't take negative powers... =| )
    if (S0New.size() > 1)
      return mpq_class(1, p^(S0New.size() - 1)) \
	* Qnew.Local_Good_Density_Congruence(p, m / p, VectorComplement(Zvec, S0New), NZvec);
    else
      return mpq_class(p^(1 - S0New.size())) \
      * Qnew.Local_Good_Density_Congruence(p, m / p, VectorComplement(Zvec, S0New), NZvec);
  }
  else
    return 0;

}



////////////////////////////////////////////////////////////////////
// Finds the Bad-type II local density of Q representing m at p.  //
// (Assuming that p > 2 and Q is given in local diagonal form.)   //
////////////////////////////////////////////////////////////////////

mpq_class Matrix_mpz::Local_BadII_Density_Congruence(const mpz_class & p, const mpz_class & m,
					 const valarray<int> & Zvec, const valarray<int> & NZvec) const
{


  /*
  // For debugging purposes:
  cout << " In Local_BadII_Density_Congruence with " << endl;
  cout << " Q is: " << endl << Q << endl;
  cout << " p = " << p << endl;
  cout << " m = " << m << endl;
  cout << " Zvec = " << Zvec << endl;
  cout << " NZvec = " << NZvec << endl;
  */


  int n;
  n = (*this).NumRows();

  int i, j;

  // Define the indexing sets S_i
  mpz_class val;
  valarray<int> S0, S1, S2;
  S0.resize(n);
  S1.resize(n);
  S2.resize(n);


  int S0_ptr, S1_ptr, S2_ptr;
  S0_ptr=0;
  S1_ptr=0;
  S2_ptr=0;

  for(i=1; i<=n; i++) {

    // Compute the valuation of each index, allowing for off-diagonal terms
    if ((*this)(i,i) == 0)
      if (i==1)
	val = Valuation((*this)(i,i+1), p);  // Look at the term to the right
      else if (i==n)
	val = Valuation((*this)(i-1,i), p);  // Look at the term above
      else
	val = Valuation((*this)(i,i+1) + (*this)(i-1,i), p);  // Finds the valuation of the off-diagonal term since only one isn't zero
    else
      val = Valuation((*this)(i,i), p);

    if (val == 0) {
      S0[S0_ptr] = i;
      S0_ptr++;
    }
    if (val == 1) {
      S1[S1_ptr] = i;
      S1_ptr++;
    }
    if (val >= 2) {
      S2[S2_ptr] = i;
      S2_ptr++;
    }
  }


  // Check that the form is primitive...
  if (S0.size() <= 0)
    {
      cerr << " Error in Local_BadI_Density_Congruence: The form is not primitive!" << endl;
      cerr << " Using Q = " << (*this) << endl;
      cerr << " and p = " << p << endl;
      exit(0);
    }


  S0 = VectorTrim(S0);
  S1 = VectorTrim(S1);
  S2 = VectorTrim(S2);

  assert(S0.size() > 0);     // Check that the S0 vector is is non-empty...


  /*
  cout << "\n Entering BII routine " << endl;
  cout << " S0 is " << S0 << endl;
  cout << " S1 is " << S1 << endl;
  cout << " S2 is " << S2 << endl;
  */

  // Note: The following lines assume that NZvec is an ordered vector.  Should check this here...


    /*
	cout << " m = " << m << "   p = " << p << "\n";
	cout << " S0 = " << S0 << "   NZvec = " << NZvec << "   IsDisjoint = " << IsDisjointOrdered(S0, NZvec)  << "\n";
	cout << " S0.size = " << S0.size() << "\n";
    */

  mpz_class p2;
  p2 = p*p;

  // Make the form Qnew for the reduction procedure
  if ((m % (p2) == 0) && (S2.size() != 0) && IsDisjointOrdered(VectorUnion(S0,S1), NZvec)) {  //  <=====  CHECK THIS!!! ****  Should this be (S0 U S1) instead of S0 ???
    //cout << "Test 1.1 \n";
    Matrix_mpz Qnew;
    Qnew = (*this);
    j=1;
    for(i=1; i<=n; i++)
      // Short circuit && required here:
      if ((j-1 < S2.size()) && (i == S2[j-1])) {                 // i is in S2
	j++;
	Qnew(i,i) = Qnew(i,i) / p2;
	if ((p == 2) && (i < n)) {
	  Qnew(i+1,i) = Qnew(i+1,i) / p2;
	  Qnew(i,i+1) = Qnew(i,i+1) / p2;
	}
      }

    /*
    cout << "\n\n Check of Bad-type II reduction: \n";
    cout << " Q is " << Q << "\n";
    cout << " Qnew is " << Qnew << "\n";
    */

    // Do the reduction
    // (Need 2 steps for each case since we can't take negative powers... =| )
    // ------------------------------------------------------------------------
    valarray<int> new_Zvec;
    new_Zvec = VectorComplement(Zvec, VectorUnion(S0, S1));


      /*
      cout << "  m = " << m << ",  m/p2 = " << (m / p2) << endl;
      cout << "  new_Zvec = " << new_Zvec << endl;
      cout << "  NZvec = " << NZvec << endl;
      cout << "  Local_Density_Congruence(Qnew, p, m / p2, new_Zvec, NZvec) = "
	   << Local_Density_Congruence(Qnew, p, m / p2, new_Zvec, NZvec) << endl;
      cout << "  Local_Density_Congruence(Qnew, p, m / p2, VectorUnion(S2, new_Zvec), NZvec)) = "
	   << Local_Density_Congruence(Qnew, p, m / p2, VectorUnion(S2, new_Zvec), NZvec) << endl;
      cout << "  CountLocalType(Qnew, 2, 2, m=1, 0, trimZvec, trimNZvec) = "
	   << CountLocalType(Qnew, 2, 2, 1, 0, Zvec, NZvec) << endl;
      cout << "  CountLocalType(Qnew, 2, 3, m=1, 0, trimZvec, trimNZvec) = "
	   << CountLocalType(Qnew, 2, 3, 1, 0, Zvec, NZvec) << endl;
      cout << "  CountLocalType(Qnew, 2, 4, m=1, 0, trimZvec, trimNZvec) = "
	   << CountLocalType(Qnew, 2, 4, 1, 0, Zvec, NZvec) << endl;
      cout << "  CountLocalType(Qnew, 2, 5, m=1, 0, trimZvec, trimNZvec) = "
	   << CountLocalType(Qnew, 2, 5, 1, 0, Zvec, NZvec) << endl;
      cout << "  CountLocalType(Q, 2, 2, m=2, 0, trimZvec, trimNZvec) = "
	   << CountLocalType(Q, 2, 2, 2, 0, Zvec, NZvec) << endl;
      cout << "  CountLocalType(Q, 2, 3, m=2, 0, trimZvec, trimNZvec) = "
	   << CountLocalType(Q, 2, 3, 2, 0, Zvec, NZvec) << endl;
      cout << "  CountLocalType(Q, 2, 4, m=2, 0, trimZvec, trimNZvec) = "
	   << CountLocalType(Q, 2, 4, 2, 0, Zvec, NZvec) << endl;
      cout << "  CountLocalType(Q, 2, 5, m=2, 0, trimZvec, trimNZvec) = "
	   << CountLocalType(Q, 2, 5, 2, 0, Zvec, NZvec) << endl;
      */

      if (n > S2.size() + 2)
	return mpq_class(1, p^(n - S2.size() - 2))
	  * (Qnew.Local_Density_Congruence(p, m / p2, new_Zvec, NZvec)
	     - Qnew.Local_Density_Congruence(p, m / p2, VectorUnion(S2, new_Zvec), NZvec));
      else
	return mpq_class(p^(S2.size() + 2 - n))
	  * (Qnew.Local_Density_Congruence(p, m / p2, new_Zvec, NZvec)
	     - Qnew.Local_Density_Congruence(p, m / p2, VectorUnion(S2, new_Zvec), NZvec));

  }
  else
    return 0;
}




///////////////////////////////////////////////////////////
// Finds the Bad-type local density of Q representing   //
// m at p, allowing certain congruence conditions mod p. //
///////////////////////////////////////////////////////////

mpq_class Matrix_mpz::Local_Bad_Density_Congruence(const mpz_class & p, const mpz_class & m,
				       const valarray<int> & Zvec, const valarray<int> & NZvec) const
{
  return (*this).Local_BadI_Density_Congruence(p, m, Zvec, NZvec) + (*this).Local_BadII_Density_Congruence(p, m, Zvec, NZvec);
}










/////////////////////////////////////////////////////////
// Local_Density and Local_Density_Congruence routines //
/////////////////////////////////////////////////////////


mpq_class Matrix_mpz::Local_Density_Congruence(const mpz_class & p, const mpz_class & m,
				   const valarray<int> & Zvec, const valarray<int> & NZvec) const
{
  return (*this).Local_Good_Density_Congruence(p, m, Zvec, NZvec) + (*this).Local_Zero_Density_Congruence(p, m, Zvec, NZvec) + (*this).Local_Bad_Density_Congruence(p, m, Zvec, NZvec);
}


// Note: The following routine is not used internally, but is included for consistency.
mpq_class Matrix_mpz::Local_Primitive_Density_Congruence(const mpz_class & p, const mpz_class & m,
					     const valarray<int> & Zvec, const valarray<int> & NZvec) const
{
  return (*this).Local_Good_Density_Congruence(p, m, Zvec, NZvec) + (*this).Local_Bad_Density_Congruence(p, m, Zvec, NZvec);
}








////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////

