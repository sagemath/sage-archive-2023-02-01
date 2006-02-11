#include "Matrix_mpz_header.h"
#include "../GMP_class_extras/mpz_class_extras.h"


///////////////////////////////////////////////////////////////////
// Swap the ith and jth rows and columns of a symmetric matrix S //
///////////////////////////////////////////////////////////////////

void Matrix_mpz::SwapSymmetric(int i, int j) {

  if ((*this).IsSymmetric()) {
    int k;
    mpz_class temp;

    // Swap Rows
    for(k=1; k<=(*this).NumRows(); k++) {
      temp = (*this)(i,k);
      (*this)(i,k) = (*this)(j,k);
      (*this)(j,k) = temp;
    }

    // Swap Columns
    for(k=1; k<=(*this).NumRows(); k++) {
      temp = (*this)(k,i);
      (*this)(k,i) = (*this)(k,j);
      (*this)(k,j) = temp;
    }
  }
  else {
    cerr << "Error in SwapSymmetric: Matrix \n" << (*this) << "\n is not symmetric!" << endl;
    abort();
  }

}



//////////////////////////////////////////////////////////////////
// Multiply the ith row and column of a symmetric matrix S by c //
//////////////////////////////////////////////////////////////////

void Matrix_mpz::MultiplySymmetric(mpz_class c, int i) {

  if ((*this).IsSymmetric()) {
      int k;

      // Row operation
      for(k=1; k<=(*this).NumRows(); k++)
	(*this)(i,k) = c * (*this)(i,k);


      // Column operation
      for(k=1; k<=(*this).NumRows(); k++)
	(*this)(k,i) = c * (*this)(k,i);

    }
  else {
    cerr << "Error in MultiplySymmetric: Matrix \n" << (*this) << "\n is not symmetric!" << endl;
    abort();
  }

}



////////////////////////////////////////////////////////////////
// Divide the ith row and column of a symmetric matrix by c //
////////////////////////////////////////////////////////////////

void Matrix_mpz::DivideSymmetric(mpz_class c, int i) {

  if ((*this).IsSymmetric()) {
    int k;

    // Rows operation
    for(k=1; k<=(*this).NumRows(); k++)
      (*this)(i,k) = (*this)(i,k) / c;

    // Column operation
    for(k=1; k<=(*this).NumRows(); k++)
      (*this)(k,i) = (*this)(k,i) / c;

  }

  else {
    cerr << "Error in MultiplySymmetric: Matrix \n" << (*this) << "\n is not symmetric!" << endl;
    abort();
  }

}



//////////////////////////////////////////////////////////////
// Replace the ith rows and columns of a symmetric matrix S //
// with the sum of the ith and c * jth rows and columns     //
//////////////////////////////////////////////////////////////

void Matrix_mpz::AddSymmetric(mpz_class c, int i, int j) {

  if ((*this).IsSymmetric()) {
    int k;

    // Row operation
    for(k=1; k<=(*this).NumRows(); k++)
      (*this)(i,k) = (*this)(i,k) + c * (*this)(j,k);

    // Column operation
    for(k=1; k<=(*this).NumRows(); k++)
      (*this)(k,i) = (*this)(k,i) + c * (*this)(k,j);
  }
  else {
    cerr << "Error in AddSymmetric: Matrix \n" << (*this) << "\n is not symmetric!" << endl;
    abort();
  }

}



//////////////////////////////////////////////////////
// Returns a normalized version of Q at the prime p //
//////////////////////////////////////////////////////

Matrix_mpz Matrix_mpz::GetLocalNormal(const mpz_class & p) const{

  // Copy the current matrix to Q
  Matrix_mpz Q;
  Q = (*this);


  if (Q.IsSymmetric()) {
    int n;
    n = Q.NumRows();

    int upper_left;
    upper_left = 1;

    int block_size;
    while(upper_left <= n){
      block_size = 0;

      // Step 1: Find the minimally p-divisible matrix entry, preferring diagonals
      // -------------------------------------------------------------------------
      unsigned long min_val = 0, tmp_val;
      bool min_defined = false;   // Flag to say whether we have found a non-zero matrix entry... (initally false)
      int i,j;
      int min_i = 0, min_j = 0;

      if (Q(upper_left, upper_left) != 0) {
	min_defined = true;
	min_val = Valuation(Q(upper_left, upper_left),p);
	min_i = min_j = upper_left;
      }

      /*
      cout << "\n Q(1,1) = " << Q(upper_left, upper_left) << "\n";
      cout << "\n min_val = " << min_val << "\n";
      */


      // Step 1a: Check the diagonal
      for(i=upper_left; i<=n; i++) {
	if (Q(i,i) != 0) {
	  tmp_val = Valuation(Q(i,i),p);
	  if ((min_defined == 0) || (tmp_val < min_val)) {
	    min_defined = true;
	    min_val = tmp_val;
	    min_i = min_j = i;
	  }
	}
      }

      /*
	cout << "\n Finished Step 1a \n";
	cout << "\n  Q = " << Q << "\n";
	cout << "upper_left = " << upper_left << "\n";
	cout << "block_size = " << block_size << "\n";
	cout << "min_i = " << min_i << "\n";
	cout << "min_j = " << min_j << "\n";
	cout << "min_val = " << min_val << "\n";
      */


      // Step 1b: Check below the diagonal
      for(i=upper_left; i<=n; i++)
	for(j=i+1; j<=n; j++) {
	  if (Q(i,j) != 0) {
	    tmp_val = Valuation(Q(i,j),p);
	    if ((min_defined == 0) || (tmp_val < min_val)) {
	      min_val = tmp_val;
	      min_i = i;
	      min_j = j;
	    }
	  }
	}


      /*
      cout << "\n Finished Step 1 \n";
      cout << "\n  Q = " << Q << "\n";
      cout << "upper_left = " << upper_left << "\n";
      cout << "block_size = " << block_size << "\n";
      cout << "min_i = " << min_i << "\n";
      cout << "min_j = " << min_j << "\n";
      cout << "min_val = " << min_val << "\n";
      */


      // Error if we still haven't seen non-zero coefficients!
      if  (min_defined == 0) {
	cout << "Error in LocalNormal: The original matrix is degeneate! \n";
      }


      // Step 2: Arrange for the upper leftmost entries to have minimal valuation
      // ------------------------------------------------------------------------

      mpz_class min_scale;
      min_scale = p ^ min_val;


      if (min_i == min_j) {
	block_size = 1;
	Q.SwapSymmetric(upper_left, min_i);
      }

      else {
	// Work in the upper-left 2x2 block, and replace it by its Z_2-equivalent form
	Q.SwapSymmetric(upper_left, min_i);
	Q.SwapSymmetric(upper_left+1, min_j);


	// 1x1 => make upper left the smallest
	if(p != 2) {
	  block_size = 1;
	  Q.AddSymmetric(1, upper_left, upper_left+1);
	}
	// 2x2 => replace it with the appropriate 2x2 matrix
	else
	  block_size = 2;
      }

      /*
      cout << "\n Finished Step 2 \n";
      cout << "\n Q is: \n" << Q << "\n\n";
      cout << "  p is: " << p << "\n";
      cout << "  min_val is: " << min_val << "\n";
      cout << "  min_scale is: " << min_scale << "\n";
      cout << "  block_size is: " << block_size << "\n";
      cout << "\n Starting Step 3 \n";
      */



      // Step 3: Clear out the remaining entries
      // ---------------------------------------

      mpz_class a,b,g;

      // Perform cancellation over Z by ensuring divisibility
      if(block_size == 1) {
	a = Q(upper_left, upper_left);
	for(j = upper_left + block_size; j<=n; j++) {
	  b = Q(upper_left, j);
	  g = GCD(a, b);

	  Q.MultiplySymmetric(a/g, j);  // Ensures that the new b entry is divisible by a
	  Q.AddSymmetric(-b/g, j, upper_left);  // Performs the cancellation
	}
      }

      mpz_class a1, a2, b1, b2, big_det, small_det;

      if(block_size == 2) {
	a1 = Q(upper_left, upper_left);
	a2 = Q(upper_left, upper_left+1);
	b1 = Q(upper_left+1, upper_left);
	b2 = Q(upper_left+1, upper_left+1);
	/*
	cout << "\n I'm here 1! \n";
	cout << "Big Det = " << (a1*b2 - a2*b1) << "\n";
	cout << "Small Det = " << ((a1*b2 - a2*b1) / (min_scale * min_scale)) << "\n";
	cout << "min_scale = " << min_scale << "\n";
	*/

	big_det = (a1*b2 - a2*b1);
	small_det = big_det / (min_scale * min_scale);

	//cout << "\n I'm here 2! \n";


	// Cancels out the rows/columns of the 2x2 block
	for(j = upper_left + block_size; j<=n; j++) {
	  a = Q(upper_left, j);
	  b = Q(upper_left+1, j);

	  // Ensures an integral result
	  Q.MultiplySymmetric(big_det, j);

	  // Performs the cancellation
	  Q.AddSymmetric(-(a*b2 - b*a2), j, upper_left);
	  Q.AddSymmetric(-(-a*b1 + b*a1), j, upper_left+1);

	  // Now remove the extra factor of big_det we introduced above
	  //	  Q = DivideSymmetric(Q, big_det, j);
	  Q.DivideSymmetric(min_scale * min_scale, j);


	}




	/*
	cout << "\n small_det = " << small_det << "\n";
	cout << "\n (1 + small_det) % 8 = " << ((1 + small_det) % 8) << "\n";
	*/


	// Uses Cassels's proof to replace the remaining 2 x 2 block
	if (((1 + small_det) % 8) == 0) {
	  Q(upper_left, upper_left) = 0;
	  Q(upper_left+1 ,upper_left+1) = 0;
	  Q(upper_left ,upper_left+1) = min_scale;
	  Q(upper_left+1 ,upper_left) = min_scale;
	}
	  else
	    if (((5 + small_det) % 8) == 0) {
	      Q(upper_left, upper_left) = 2 * min_scale;
	      Q(upper_left+1 ,upper_left+1) = 2 * min_scale;
	      Q(upper_left ,upper_left+1) = min_scale;
	      Q(upper_left+1 ,upper_left) = min_scale;
	    }
	    else {
	      cout << "Error in LocalNormal: Impossible behavior for a 2x2 block! \n";
	      block_size = 1000;
	    }



      }


      upper_left = upper_left + block_size;

      /*
      cout << "\n Finished Step 3 \n";
      cout << "\n Q is: \n" << Q << "\n\n";
      int dummy;
      cin >> dummy;
      */

    }


    // Finally, if p = 2 then make it upper triangular and take half of the matrix
    // (since we assume the given matrix was for 2*Q, but the local output will be Q)
    if (p==2) {
      for(int i=1; i<=n-1; i++)
	for(int j=i+1; j<=n; j++) {
	  //	  Q(i,j) += Q(j,i);  //  <<=====  By removing this, since the original matrix was symmetric we take halve the off-diagonal entries
	  Q(j,i) = 0;
	}
    }


    // Now take half of the diagonal
    // (this should be done whether p=2 or p>2!)
    for(int i=1; i<=n; i++)
      Q(i,i) = Q(i,i) / 2;


    return Q;
  }

  cerr << "Error in LocalNormal: Matrix \n" << Q << "\n is not symmetric!" << endl;
  abort();
}

