//////////////////////////////////////////////////////////////////////////
//
// To Do: (9/24/04)
//   - Replace the determinant routine with something working in any ring.
//   - Make this a templated class
//   - Change EvaluateQuadratic(v,m) to EvaluateQuadraticMod(v,m)
//   - Allow vector and valarray inputs in our routines
//
//////////////////////////////////////////////////////////////////////////


using namespace std;

#include <gmp.h>
#include <gmpxx.h>
#include <valarray>
#include <iostream>
#include <sstream>

#include "Matrix_mpz_header.h"
#include "../GMP_class_extras/mpz_class_extras.h"

Matrix_mpz::Matrix_mpz(int r, int s) {
    M.resize(r*s);

    // Clear all entries
    for(int i=0; i<(r*s); i++)
      M[i] = 0;

    m = r;
    n = s;
}

Matrix_mpz::Matrix_mpz() {
  M.resize(0);
    m = 0;
    n = 0;
}


///////////////////////////////////////////////////////////////////
// Copy Constructor -- Needed so that the valarray is resized! =)
///////////////////////////////////////////////////////////////////
void Matrix_mpz::operator=(const Matrix_mpz & source) {
  // Protect against self-assignment
  if (this != &source) {
    m = source.m;
    n = source.n;

    /*
    // Redid this since I heard there was a valarray resize bug!
    valarray<mpz_class> N(m*n,0);
    for(int i=0; i < m*n; i++)
    N[i] = source[i];
    M = N;
    */

    //      /*
    // This was the original way...
    //    M.resize(m*n, 0);
    M = source.M;
    //      */
  }
}



// Comparison operator
bool Matrix_mpz::operator==(const Matrix_mpz & source) const {

  // Check the sizes are the same
  if ((m != source.m) || (n != source.n))
    return false;

  // Check the entries are the same
  for(int i=1; i<=m; i++)
    for(int j=1; j<=n; j++)
      if (source(i,j) != (*this)(i,j))
	return false;

  // If so, then return true
  return true;
}



/* UNABLE TO CORRECTLY DEFINE THE ADDITION OF TWO MATRICES!!

// Define the addition of two matrices
//  void operator+(const Matrix_mpz B, const Matrix_mpz C) {
//    Matrix_mpz operator+(const Matrix_mpz B, const Matrix_mpz C) {
//    Matrix_mpz operator+(const Matrix_mpz & B, const Matrix_mpz & C) {
//    Matrix_mpz & operator+ (Matrix_mpz B, Matrix_mpz C) {
//  void operator+(Matrix_mpz & A, const Matrix_mpz & B, const Matrix_mpz & C) {
//    void operator+(Matrix_mpz & B, const Matrix_mpz & C) {
//  void operator+(Matrix_mpz B, Matrix_mpz C) {

// Check sizes are the same before adding
if ((B.m == C.m) && (B.n == C.n)) {
Matrix_mpz A;
A.m = B.m;
A.n = B.n;

A.M.resize(m*n);
A.M = B.M + C.M;
return A;
}
else
cout << "Error in matrix addition: They aren't the same size!" << endl;

}

*/


int Matrix_mpz::NumRows() const {
  return m;
}

int Matrix_mpz::NumCols() const {
  return n;
}

int Matrix_mpz::Length() const {
  return m*n;
}



//////////////////////////////
// Allow the notation M % R
///////////////////////////////
Matrix_mpz Matrix_mpz::operator%(const mpz_class & R) const {

  // Check the modulus is > 0
  if (R < 0) {
    cout << "Error in % operator: The modulus must be positive!" << endl;
    exit(1);
  }

  Matrix_mpz M_mod(m,n);

  for(int i=1; i<=m; i++)
    for(int j=1; j<=n; j++)
      M_mod(i,j) = (((*this)(i,j) % R) + R) % R;

  return M_mod;
}




////////////////////////////////
// Allow the notation M[n*i+j]
////////////////////////////////
mpz_class & Matrix_mpz::operator[](int ind) {
  return M[ind];
}

const mpz_class & Matrix_mpz::operator[](int ind) const {
  return M[ind];
}


/////////////////////////////////
// Allow the notation M(n*i+j)
////////////////////////////////
mpz_class & Matrix_mpz::operator()(int ind) {
  return M[ind - 1];
}

const mpz_class & Matrix_mpz::operator()(int ind) const {
  return M[ind - 1];
}

void Matrix_mpz::set(int row, int col, const mpz_class& x) {
  M[row*n + col] = x;
}

mpz_class Matrix_mpz::get(int row, int col) {
  return M[row*n + col];
}


///////////////////////////////
// Allow the notation M(i,j)
///////////////////////////////
mpz_class & Matrix_mpz::operator()(int row, int col) {
  if ((1 <= row) && (1 <= col) && (row <= m) && (col <= n))   // This does some basic error checking
    return M[n*(row-1) + (col - 1)];


  cerr << "Error in matrix read: row " << row << " and col " << col  << " are out of range... =|" << endl;
  abort();
}
const mpz_class & Matrix_mpz::operator()(int row, int col) const {
  if ((1 <= row) && (1 <= col) && (row <= m) && (col <= n))   // This does some basic error checking
    return M[n*(row-1) + (col - 1)];


  cerr << "Error in matrix read: row " << row << " and col " << col  << " are out of range... =|" << endl;
  abort();
}



/*
// WARNING: We can't do this since it takes only one argument!!!

// Allow the notation M[i,j]
  mpz_class & Matrix_mpz::operator[](int row, int col) {
    if ((0 <= row < m) && (0 <= col < n))   // This does some basic error checking
      return M[n * row + col];

    cerr << "Error in matrix read: row " << row << " and col " << col  << " are out of range... =|" << endl;
    abort();
  }

  // Allow the notation M[i,j]
  const mpz_class & Matrix_mpz::operator[](int row, int col) const {
    if ((0 <= row < m) && (0 <= col < n))   // This does some basic error checking
      return M[n * row + col];

    cerr << "Error in matrix read: row " << row << " and col " << col  << " are out of range... =|" << endl;
    abort();
  }
*/



// Define Matrix multiplication
Matrix_mpz Matrix_mpz::operator*(const Matrix_mpz & B) {

  // Abort if we can't multiply them
  if (NumCols() != B.NumRows()) {
    cout << "Error in Matrix Multiplication: Trying to multiply a "
	 << NumRows() << " x " << NumCols() << " matrix by a "
	 << B.NumRows() << " x " << B.NumCols() << " matrix." ;
    exit(1);
  }

  // Do the multiplication
  Matrix_mpz C((*this).NumRows(), B.NumCols());
  for (long i=1; i<=C.NumRows(); i++)
    for (long j=1; j<=C.NumCols(); j++) {
      C(i,j) = 0;
      //      cout << "The (" << i << "," << j << ") entry is " << C(i,j) << endl;
      for (long k=1; k<=(*this).NumCols(); k++)
	C(i,j) += (*this)(i,k) * B(k,j);
      //      cout << "The (" << i << "," << j << ") entry is " << C(i,j) << endl;
    }

  // Return the product
  return C;
}



//////////////////////
// Prints the Matrix
//////////////////////
void Matrix_mpz::Print(ostream & out) const {

  /*
    cout << " m = " << m << "  n = " << n << endl;
    cout << " length of the valarray = " << M.size() << endl;
  */

  for(int i = 1; i <= m; i++) {
    out << " [ ";
    for(int j = 1; j <= n; j++) {
      out << (*this)(i,j);
      if (j <= n - 1)
	out << ", ";
    }
    out << " ]" << endl;
  }
}

void Matrix_mpz::Print() const {
  Print(cout);
}

string Matrix_mpz::repr() const {
  ostringstream tmp;
  Print(tmp);
  return tmp.str();
}

void Matrix_mpz::PrintM(ostream & out) const {
  Print(out);
}







////////////////////////////////////
// Writes a matrix to the ostream
///////////////////////////////////
void Matrix_mpz::_FileOutput(ostream & out) const {

  // Run through the rows
  out << " [ ";
  for(int i = 1; i <= m; i++) {

    // Run through the columns
    for(int j = 1; j <= n; j++) {
      out << (*this)(i,j);
      if (j <= n-1)
	out << ", ";
    }
    if (i <= m-1)
      out << " ; ";
  }
  out << " ] ";

}


////////////////////////////////////
// Reads a matrix from the istream
////////////////////////////////////
void Matrix_mpz::_FileInput(istream & in) {

  // Read the matrix entries into a temporary vector first
  vector<mpz_class> tmp_vec;
  long rows, cols;
  char ch;
  mpz_class num;

  // Read the opening " [ "
  in >> ch;

  // Quick check that the matrix isn't empty
  // (if so, we'll never enter the loop!)
  in >> ch;
  in.putback(ch);

  // Search for the closing bracket
  while (ch != ']') {

    // Read the next matrix entry and append it to the vector
    in >> num;
    tmp_vec.push_back(num);

    // Read the separator (',' or ';' or ']')
    in >> ch;
    assert((ch == ',') || (ch == ';') || (ch == ']'));

    // Check for ';' to define the number of columns
    if (ch == ';') {
      if (n == 0)
	n = tmp_vec.size();                   // Set the number of columns if it's the first ';'
      else
	assert( tmp_vec.size() % n == 0 );    // Check that the # of entries is a multiple of the # of cols.
    }

  }

  // Set the matrix size and entries
  M = tmp_vec;
  assert( tmp_vec.size() % n == 0 );    // Check that the # of entries is a multiple of the # of cols.
  m = tmp_vec.size() / n;

}







////////////////////////////////////////
// Makes the QF filename for this form
////////////////////////////////////////
string Matrix_mpz::QF_String() const {

  // Make the string identifying this form
  string form_string;
  for(long i=1; i<=NumRows(); i++)
    for(long j=1; j<=i; j++)
      form_string.append(MakeString((*this)(i,j).get_si()) + "_");

  // Erase the extra "_" at the end.
  form_string.erase(form_string.size() - 1);

  return form_string;
}




// Check this: would like to zero it out (now)                // TO DO: Get rid of this routine!!! =)
//   or to preserve the existing matrix... (later?)

void Matrix_mpz::SetDims(int r, int s){
  Matrix_mpz QQ(r,s);
  (*this) = QQ;
}



//////////////////////////////////////////////////////
// Resize the matrix M to size m x n, preserving    //
// existing entries and extending by zero elsewhere //
//////////////////////////////////////////////////////

void Matrix_mpz::SafeResize(int m, int n)
{
  Matrix_mpz New(m,n);

  // Copy matrix entries
  for(int i=1; i <= min((*this).NumRows(), m); i++)
    for(int j=1; j <= min((*this).NumCols() ,n); j++)
      New(i,j) = (*this)(i,j);

  (*this) = New;
}








bool Matrix_mpz::IsSquare() const {
  return (m == n);
}


bool Matrix_mpz::IsSymmetric() const {
  bool flag = false;

  if (m == n) {
    flag = true;
    for(int i=1; i<=n; i++)
      for(int j=1; j<=m; j++)
	if ((*this)(i,j) != (*this)(j,i))
	  flag = false;
  }

  return flag;
  }



bool Matrix_mpz::IsQuadraticForm() const {

  // Check if it's symmetric
  if ((*this).IsSymmetric() == false)
      return false;

  // Check that the diagonal is even
  for(int i=1; i<=n; i++)
    if ((*this)(i,i) % 2 != 0)
      return false;

  // Ok, it's a quadratic form
  return true;
}




////////////////////////////////////////////////////////
// Makes a Diagonal matrix from a valarray<mpz_class> //
////////////////////////////////////////////////////////

void Matrix_mpz::DiagonalMatrix(const valarray<mpz_class> & v) {

  Matrix_mpz QQ(v.size(), v.size());

  for(int i=1; i<=v.size(); i++)
    QQ(i,i) = v[i-1];

  (*this) = QQ;

}



////////////////////////////////////////////////////////////////
// Converts the current matrix into the n x n identity matrix //
////////////////////////////////////////////////////////////////

void Matrix_mpz::IdentityMatrix(int n)
{
  Matrix_mpz QQ(n,n);

  for(int i=1; i<=n; i++)
    QQ(i,i) = 1;

  (*this) = QQ;
}















void Matrix_mpz::Transpose() {

  vector<mpz_class> N(n*m);  // WANT TO AVOID USING AN EXPLICIT VECTOR TYPE HERE!
  // cout << " m = " << m << "  n = " << n << endl;

  for(int i=1; i<=n; i++)
    for(int j=1; j<=m; j++)
      N[m*(i-1) + (j-1)] = (*this)(j,i);

  M = N;
  swap(m,n);


  /*
  // This looks nicer, but it doesn't work...
  Matrix_mpz N((*this).NumCols(),(*this).NumRows());

  for(int i=1; i<=(*this).NumRows(); i++)
    for(int j=1; j<=(*this).NumCols(); j++)
      N(j,i) = (*this(i,j));

  *this = N;
  */
}


Matrix_mpz Matrix_mpz::GetTranspose() const {

  Matrix_mpz Trans;
  Trans = *this;
  Trans.Transpose();

  return Trans;
}




//////////////////////////////////////////////////////////////////////////////
// Quickly evaluate the expression v^t * M * v (mod R)   <-- Should change the name to EvaluateQuadraticMod()...
//////////////////////////////////////////////////////////////////////////////
mpz_class Matrix_mpz::EvaluateQuadratic(const valarray<mpz_class> & v, const mpz_class & R) const {
  mpz_class total;
  total = 0;

  // Check the modulus is > 0
  if (R < 0) {
    cout << "Error in EvaluateQuadratic(): The modulus must be positive!" << endl;
    exit(1);
  }

  // Check the matrix is square
  if ((*this).NumRows() != (*this).NumCols()) {
    cout << "Error in EvaluateQuadratic(): The matrix is not square!" << endl;
    exit(1);
  }

  // Check the valarray has the appropriate size
  if (v.size() != (*this).NumRows()) {
    cout << "Error in EvaluateQuadratic(): The valarray is the wrong size!" << endl;
    exit(1);
  }

  // Evaluate v^T * Q * v (mod R)
  for(int i=1; i<=m; i++)
    for(int j=1; j<=n; j++)
      total = (total + (v[i-1] * (*this)(i,j) * v[j-1])) % R;

  /*
    if (total >= R)
    cout << " Error in EvaluateQuadratic: R exceeded for vector: " << v << endl;

    if (total < 0)
    cout << " Error in EvaluateQuadratic: Negative value for vector: " << v << endl;
  */

  return (total + R) % R;  // This is necessary for some strange GMP reason...  (see 2/18/04 Notes...)
}



/////////////////////////////////////////
// Evaluates T^t * Q * T for a matrix T
////////////////////////////////////////
Matrix_mpz Matrix_mpz::EvaluateQuadratic(const Matrix_mpz & T) const {

  // Check the matrix is square
  if ((*this).NumRows() != (*this).NumCols()) {
    cout << "Error in EvaluateQuadratic(): The matrix is not square!" << endl;
    exit(1);
  }

  // Check that the new matrix T has the correct number of rows
  if (T.NumRows() != (*this).NumRows()) {
    cout << "Error in EvaluateQuadratic(): The new matrix has the wrong number of rows!" << endl;
    exit(1);
  }

  Matrix_mpz ans(T.NumRows(), T.NumRows());

  for(int i=1; i<=T.NumRows(); i++)
    for(int j=1; j<=T.NumRows(); j++)
      for(int k=1; k<=m; k++)
	for(int l=1; l<=n; l++)
	  ans(i,j) += T(k,i) * (*this)(k,l) * T(l,j);

  return ans;
}


/////////////////////////////////////////////
// Evaluates v^t * Q * v for a valarray v
/////////////////////////////////////////////
mpz_class Matrix_mpz::EvaluateQuadratic(const valarray<mpz_class> & v) const{

  // Check the matrix is square
  if ((*this).NumRows() != (*this).NumCols()) {
    cout << "Error in EvaluateQuadratic(): The matrix is not square!" << endl;
    exit(1);
  }

  // Check the valarray has the appropriate size
  if (v.size() != (*this).NumRows()) {
    cout << "Error in EvaluateQuadratic(): The valarray is the wrong size!" << endl;
    exit(1);
  }
  // Should check that *this is a square matrix, and T has the same # of rows.

  mpz_class ans;
  ans = 0;

  for(int k=1; k<=m; k++)
    for(int l=1; l<=n; l++)
      ans += v[k-1] * (*this)(k,l) * v[l-1];

  return ans;
}


///////////////////////////////////////////////////////
// Evaluates v^t * Q * v for a vector v (of mpz_class)
///////////////////////////////////////////////////////
mpz_class Matrix_mpz::EvaluateQuadratic(const vector<mpz_class> & v) const{

  // Check the matrix is square
  if ((*this).NumRows() != (*this).NumCols()) {
    cout << "Error in EvaluateQuadratic(): The matrix is not square!" << endl;
    exit(1);
  }

  // Check the vector has the appropriate size
  if (v.size() != (*this).NumRows()) {
    cout << "Error in EvaluateQuadratic(): The valarray is the wrong size!" << endl;
    exit(1);
  }

  mpz_class ans;
  ans = 0;

  for(int k=1; k<=(*this).NumRows(); k++)
    for(int l=1; l<=(*this).NumCols(); l++)
      ans += v[k-1] * (*this)(k,l) * v[l-1];

  return ans;
}



////////////////////////////////////////////////////
// Evaluates v^t * Q * v for a vector v (of longs)
/////////////////////////////////////////////////////
mpz_class Matrix_mpz::EvaluateQuadratic(const vector<long> & v) const{

  // Check the matrix is square
  if ((*this).NumRows() != (*this).NumCols()) {
    cout << "Error in EvaluateQuadratic(): The matrix is not square!" << endl;
    exit(1);
  }

  // Check the vector has the appropriate size
  if (v.size() != (*this).NumRows()) {
    cout << "Error in EvaluateQuadratic(): The valarray is the wrong size!" << endl;
    exit(1);
  }

  mpz_class ans;
  ans = 0;

  for(int k=1; k<=(*this).NumRows(); k++)
    for(int l=1; l<=(*this).NumCols(); l++)
      ans += v[k-1] * (*this)(k,l) * v[l-1];

  return ans;
}





// Extracts a square matrix according to the entries in Index
// (Note: The indexing in Index starts at 1, not at zero.)
Matrix_mpz Matrix_mpz::ExtractSquareSubmatrix(const valarray<int> & Index) const {
  int len;
  len = Index.size();

  int i, j, max_ind;

  // Check that the biggest entry of Index is in range of the matrix M
  max_ind = Index.max();
  if ( (max_ind > NumRows()) || (max_ind > NumCols()) ) {
    cout << "\n Error in ExtractSquareSubmatrix: The index vector is out of range! \n\n";
    return Matrix_mpz();
  }

  // Extract the appropriate entries into Mnew
  Matrix_mpz Mnew(len, len);

  /*
    cout << " Just created Mnew: " << endl;
    Mnew.PrintM(cout);
    cout << endl;
  */

  for (i=1; i<=len; i++)
    for (j=1; j<=len; j++)
      Mnew(i,j) = (*this)(Index[i-1],Index[j-1]);

  /*
    cout << " Now we have Mnew: " << endl;
    Mnew.PrintM(cout);
    cout << endl;
  */

  return Mnew;
}



// Extracts a column vector from a matrix

valarray<mpz_class> Matrix_mpz::ExtractColumn(int j) const {

  valarray<mpz_class> col;
  int i;

  if ( (j > n) || (j <= 0) ) {
    cout << "\n Error in ExtractColumn: The column index " \
	 << j << " exceeds the number of columns " << m << ".\n";
    return col;
  }

  col.resize(n);
  for(i=1; i<=m; i++)
    col[i-1] = (*this)(i,j);

  return col;
}



// Extracts a row vector from a matrix

valarray<mpz_class> Matrix_mpz::ExtractRow(int i) const {

  valarray<mpz_class> row;
  int j;

  if ( (i > m) || (i <= 0) ) {
    cout << "\n Error in ExtractColumn: The row index " \
	 << i << " exceeds the number of rows " << n << ".\n";
    return row;
  }

  row.resize(m);
  for(j=1; j<=n; j++)
    row[j-1] = (*this)(i,j);

  return row;
}



// Extracts a square submatrix with indices associated to the given vector

Matrix_mpz Matrix_mpz::ExtractSquareSubmatrixOrdered(const valarray<int> & Index) const
{
  int len;
  len = Index.size();

  int i, j, max;
  // Check that the biggest entry of Index is in range of the matrix M
  max = 0;
  for (i=1; i<=len; i++)
    if (Index[i-1] > max)
      max = Index[i-i];

  /*
    cout << "\n ExtractSquareSubmatrixOrdered Status: \n";
    cout << " Matrix M = " << M << "\n";
    cout << " Index vector = " << Index << "\n";
    cout << " maximum index vector entry= " << max << "\n";
  */

  if ( (max > m) || (max > n) ) {
    cout << "\n Error in ExtractSquareSubmatrixOrdered: The index vector is out of range! \n\n";
    return *this;
  }

  // Extract the appropriate entries
  Matrix_mpz Mnew(len,len);
  for (i=1; i<=len; i++)
    for (j=1; j<=len; j++)
      Mnew(i,j) = (*this)(Index[i-1],Index[j-1]);

  return Mnew;
}






// Swaps the ith and jth rows of our matrix

void Matrix_mpz::SwapRows(int i, int j){

  mpz_class temp;
  int k;

  for (k=1; k<=n; k++) {
    temp = (*this)(i,k);
    (*this)(i,k) = (*this)(j,k);
    (*this)(j,k) = temp;
  }

}





/*!  \brief Computes the determinant of a Matrix_mpz using explicit row reduction over mpq_class.
 *
 *
 * \todo This should be rewritten to only work over mpz_class, so we can template the Matrix class! =)
 *
 */

// Computes the Determinant
mpz_class Matrix_mpz::Determinant() const {

  // Check that the matrix is square
  if (IsSquare() && (m>0)) {
    int i, j, k;
    Matrix_mpz Temp;
    Temp = *this;
    mpz_class extra = 1;

    // Check for the 1 x 1 case first
    if (m == 1)
      return Temp(1,1);

    mpz_class A_ij, A_jj;
    // Go through the columns in order
    for (j=1; j<=n-1; j++) {
      for (i=j+1; i<=m; i++) {

	/*
	  cout << "Before: " << endl;
	  cout << " i = " << i << "  j = " << j << endl;
	  cout << " Extra factor = " << extra << endl;
	  Temp.Print(cout);
	  cout << endl;
	*/

	// Check to see if our pivot entry is zero
	if (Temp(j,j) == 0) {
	  for(k=j+1; k<=n; k++)
	    if (Temp(k,j) != 0) {
	      Temp.SwapRows(j,k);
	      extra *= -1;
	      break;
	    }
	}

	A_ij = Temp(i,j);
	A_jj = Temp(j,j);


	// If a non-zero pivot is found, use it to do an elementary row operation
	if (Temp(j,j) != 0) {

	  // Perform the linear combimation [Row i -> - Temp(i,j) * Row i + Temp(j,j) * Row j]
	  for(k=1; k<=n; k++) {
	    //    cout << " - " << Temp(i,j) << " * " <<  Temp(j,k) << " + " << Temp(j,j) << " * " << Temp(i,k) << endl;

	    Temp(i,k) = - A_ij * Temp(j,k) + A_jj * Temp(i,k);
	  }

	  extra *= Temp(j,j);
	}

	/*
	  cout << "After: " << endl;
	  cout << " i = " << i << "  j = " << j << endl;
	  cout << " Extra factor = " << extra << endl;
	  Temp.Print(cout);
	  cout << endl;
	*/

      }
    }

    // Now multiply by the determinant of the upper triangular matrix Temp

    mpz_class diag_det;
    diag_det = 1;
    for(k=1; k<=m; k++)
      diag_det *= Temp(k,k);

    //      cout << "the det is " << (diag_det / extra) << endl;

    return (diag_det / extra);
  }


  cerr << " Error in Determinant method: The matrix is not square or it has size 0! =(" << endl;
  cerr << " Using Q = ";
  PrintM(cout);
  cerr << endl;
  cerr << " where m = " << m << "  and  n = " << n << endl;
  abort();
}






/*!  \brief Computes the adjoint of a Matrix_mpz.
 *
 *
 */


// Finds the adjoint of a Matrix_mpz
Matrix_mpz Matrix_mpz::Adjoint() const {

  // Note: The Adjoint matrix should be square...
  Matrix_mpz Adj(m,n);

  for(int i=1; i<=m; i++)
    for(int j=1; j<=n; j++){

      int k_ptr, l_ptr;
      Matrix_mpz Minor_mat(m-1, n-1);

      //	cout << "\n i = " << i << "  and  j = " << j << endl << endl;

      // Loop to construct the (m-1) x (n-1) submatrix
      k_ptr = 1;
      for(int k=1; k<=m-1; k++) {

	if (k == i)   // Skip the i-th row
	  k_ptr++;

	l_ptr = 1;
	for(int l=1; l<=n-1; l++){

	  if (l == j)   // Skip the j-th column
	    l_ptr++;

	  /*
	      cout << " k = " << k << "  and  l = " << l << endl;
	      cout << " k_ptr = " << k_ptr << "  and  l_ptr = " << l_ptr << endl;
	  */

	  Minor_mat(k, l) = (*this)(k_ptr, l_ptr);

	  l_ptr++;
	}

	k_ptr++;
      }

      /*
	cout << " The matrix minor is \n";
	Minor_mat.PrintM(cout);
	cout << endl;
      */

      mpz_class sign;
      if ((i+j) % 2 == 0)
	sign = mpz_class(1);
      else
	sign = mpz_class(-1);


      Adj(j,i) = sign * (Minor_mat).Determinant();

    }

  return Adj;

}



//  Matrix_mpz Inverse()  //  <---- Need an mpq_class type to do this...













/*!  \brief Computes the level of the (global) form $Q$ given by the matrix $2*Q$.
 *
 *
 * \todo This should take an input number \f$t\f$ for the number of tries.
 *
 */

// Computes the level of the form 2*Q
  mpz_class Matrix_mpz::QFLevel() const {  // WARNING: FORGOT THE FACTOR OF 2

    // Sanity Check
    assert( m == n );


    mpz_class temp_lvl, det;
    det = (*this).Determinant();

    // Check the determinant isn't zero
    if (det == 0){
      cout << " Error in QFLevel(): The determinant of the matrix is zero." << endl;
      cout << "   The (singular) matrix is: " << endl << (*this) << endl;
      exit(1);
    }

    /*
    // DIAGNSOTIC
    cout << " Test1 " << endl;
    cout << "   Matrix = " << endl << (*this) << endl;
    cout << "   m = " << m << endl;
    cout << "   n = " << n << endl;
    cout << "   det = " << det << endl;
    */

    temp_lvl = 1;
    Matrix_mpz Adj;
    Adj = (*this).Adjoint();

    /*
    // DIAGNSOTIC
    cout << " Test2 " << endl;
    cout << "   Adjoint = " << endl << Adj << endl;
    */

    // Find the LCM of the denominators of the upper triangular entries
    for(int i=1; i<=m-1; i++)
      for(int j=i+1; j<=n; j++) {
	/*
	// DIAGNOSTIC
	cout << "A" << endl;
	cout << "i = " << i << endl;
	cout << "j = " << j << endl;
	cout << "Adj(i,j) = " << Adj(i,j) << endl;
	cout << "det = " << det << endl;
	cout << "mpq_class(Adj(i,j), det) = " << mpq_class(Adj(i,j), det) << endl;
	*/
	mpq_class temp_inv_entry = mpq_class(Adj(i,j), det);
	//cout << "B" << endl;
	temp_inv_entry.canonicalize();
	//cout << "C" << endl;
	temp_inv_entry = temp_inv_entry.get_den();
	//cout << "D" << endl;
	temp_lvl = LCM(temp_lvl, temp_inv_entry);
	//cout << "E" << endl;
      }

    /*
    // DIAGNSOTIC
    cout << " Test3 " << endl;
    */

    // Find the LCM of the denominators of the diagonal entries
    for(int i=1; i<=m; i++) {
      mpq_class temp_inv_entry = mpq_class(Adj(i,i), det);
      temp_inv_entry.canonicalize();
      temp_inv_entry = temp_inv_entry.get_den();
      temp_lvl = LCM(temp_lvl, 4 * temp_inv_entry);
    }

    /*
    // DIAGNSOTIC
    cout << " Test4 " << endl;
    */


    // Remove a factor of 2 here since we are using the global matrix for 2*Q
    // (Warning: This means that we are assuming this routine is only used globally!)
    temp_lvl = temp_lvl / 2;


    return temp_lvl;

  }





/////////////////////////////////////////////////////////////
// Finds the Cholesky decomposition of a quadratic form -- as an upper-triangular matrix!
// (It's assumed to be global, hence twice the form it refers to.)
/////////////////////////////////////////////////////////////

vector< vector<double> > Matrix_mpz::CholeskyDecomposition() const {

  // Declare some initial variables
  Matrix_mpz QQ = (*this);
  long n = QQ.NumRows();

  // Make Q a vector which is (n+1) by (n+1)       <--- Only the positive indices (1 --> n) are used... =)
  vector< vector<double> > Q;
  Q.resize(n+1);
  for(long i=0; i<=n; i++)           // Note: We set Q[0] just for consistency...
    Q[i].resize(n+1);


  // Check the matrix is square, and determine its dimension.
  if ((QQ.NumRows() != QQ.NumCols()) || (QQ.NumRows() <= 0)) {
    cout << "Error in CholeskyDecomposition():  The matrix is empty or not square!" << endl;
    exit(1);
  }


  // 1. Initialize (from a symmetric matrix QQ) -- using 2 * QQ right now...
  long counter = 0;
  for (long i=1; i<=n; i++)
    for (long j=i; j<=n; j++) {
      Q[i][j] = 0.5 * QQ(i,j).get_d();
      counter++;
    }
  long i = 0;


  // 1a. Print the resulting form
  /*
  cout << endl;
  for (long r=1; r<=n; r++) {
    cout << "[ ";
    for (long s=1; s<=n; s++)
      cout << Q[r][s] << " ";
    cout << "]" << endl;
  }
  cout << endl;
  */


  // 2. Loop on i
  i = i + 1;
  while (i != n) {
    for (long j=i+1; j<=n; j++) {
      Q[j][i] = Q[i][j];
      Q[i][j] = Q[i][j] / Q[i][i];
    }

    // 2a. Print the resulting form
    /*
    cout << endl;
    for (long r=1; r<=n; r++) {
      cout << "[ ";
      for (long s=1; s<=n; s++)
        cout << Q[r][s] << " ";
      cout << "]" << endl;
    }
    cout << endl;
    */

    // 3. Main Loop
    for (long k=i+1; k<=n; k++)
      for (long l=k; l<=n; l++)
        Q[k][l] = Q[k][l] - Q[k][i] * Q[i][l];


    // 3a. Print the resulting form
    /*
    cout << endl;
    for (long r=1; r<=n; r++) {
      cout << "[ ";
      for (long s=1; s<=n; s++)
      cout << Q[r][s] << " ";
      cout << "]" << endl;
    }
    cout << endl;
    */

    i = i + 1;
  }


  // 4. Zero out the strictly lower-triangular entries
  for (long i=1; i<=n; i++)
    for (long j=1; j<i; j++)
      Q[i][j] = 0;



  // 4a. Print the resulting form
  /*
  cout << endl;
  cout << " We obtain the Cholesky Decomposition of:" << endl;
  for (long r=1; r<=n; r++) {
    cout << "[ ";
    for (long s=1; s<=n; s++)
      cout << Q[r][s] << " ";
    cout << "]" << endl;
  }
  cout << endl;
  */


  // Return the Cholesky decomposition
  return Q;

}





//////////////////////////////////////////////////////////////////////
// Writes the theta function of a ternary form of desired precision //
//////////////////////////////////////////////////////////////////////

PowerSeries<mpz_class> Matrix_mpz::ComputeTheta(const unsigned long & precision) const {
  /*
  // Print the ternary form and its level
  cout << endl;
  cout << " The form is:  [ "
       << QQ[0] << ", "
       << QQ[1] << ", "
       << QQ[2] << ", "
       << QQ[3] << ", "
       << QQ[4] << ", "
       << QQ[5] << " ] " << endl;
  cout << " The level of the form is " << QF_Ternary_Level(QQ) << endl;
  */

  // Make the power series for the theta function
  PowerSeries<mpz_class> theta(precision);


  // Find the (lower-triangular) Cholesky Decomposition (uses indices 1 --> n)
  vector<vector<double> > Cholesky;
  Cholesky = CholeskyDecomposition();
  //  cout << " The Cholesky decomposition is: " << endl << Cholesky << endl;
  //  cout << " Returned with the decomposition..." << endl;


  /*
  cout << " Computing the theta function " << endl;
  PrintTime();

  cout << endl << "Entering FastBinaryTheta" << endl;
  cout << " Using precision = " << precision() << endl;
  */

  // ERROR CHECKING: Check that C+1 fits in an unsigned long.



  const long n = NumRows();

  // Make the (constant) matrix Q from QQ  -- NOTE: We only use indices 1 --> n.
  double Q[n+1][n+1];
  for (long i=1; i<=n; i++)
    for (long j=1; j<=n; j++)
      Q[i][j] = 0;   // Clear the matrix
  long counter = 0;
  for (long i=1; i<=n; i++)
    for (long j=i; j<=n; j++) {
      Q[i][j] = Cholesky[i][j];  // Put Cholesky in the upper triangular part
      counter++;
    }


  /*
  // Print Q
  cout << "Using Q = " << endl;
  for (long i=1; i<=n; i++) {
    for (long j=1; j<=n; j++)
      cout << Q[i][j] << " ";
    cout << endl;
  }
  */


  // 1. Initialize
  long i = n;
  vector<double> T(n+1, 0);  // Note: We use n+1 so we can index the entries as 1 --> n
  vector<double> U(n+1, 0);
  T[i] = (double) precision;
  U[i] = 0;
  double Z;
  vector<long> L(n+1, 0);
  vector<long> x(n+1, 0);


  // 2. Compute bounds
  Z = sqrt(T[i] / Q[i][i]);
  L[i] = long(floor(Z - U[i]));  // Check this is ok...
  /*
  cout << " L[i] float gives :     " << floor(Z - U[i]) << endl;
  cout << " L[i] mpz_class gives : " << L[i] << endl;
  */

  x[i] = long(ceil(-Z - U[i]) - 1);  // Check this is ok...
  /*
  cout << " x[i] float gives :     " << ceil(-Z - U[i] -1) << endl;
  cout << " x[i] mpz_class gives : " << x[i] << endl;
  cout << endl;
  */


  bool done_flag = false;
  double Q_val_double;
  unsigned long Q_val;                 // WARNING: Still need a good way of checking overflow for this value...


  // Big loop which runs through all vectors
  while (done_flag == false) {

    // Loop through until we get to i=1 (so we defined a vector x)
    do {

      // 3a. Main loop
      x[i] = x[i] + 1;
      while (x[i] > L[i]) {
	i = i + 1;
	x[i] = x[i] + 1;
      }

      // 3b. Main loop
      if (i>1) {
	/*
	cout << " i = " << i << endl;
	cout << " T[i] = " << T[i] << endl;
	cout << " Q[i][i] = " << Q[i][i] << endl;
	cout << " x[i] = " << x[i] << endl;
	cout << " U[i] = " << U[i] << endl;
	cout << " x[i] + U[i] = " << (x[i] + U[i]) << endl;
	cout << " T[i-1] = " << T[i-1] << endl;
	*/
	T[i-1] = T[i] - Q[i][i] * (x[i] + U[i]) * (x[i] + U[i]);
	/*
	cout << " T[i-1] = " << T[i-1] << endl;
	cout << endl;
	*/
	i = i - 1;
	U[i] = 0;
	for(long j=i+1; j<=n; j++)
	  U[i] = U[i] + Q[i][j] * x[j];

	// Now go back and compute the bounds...
	// 2. Compute bounds
	Z = sqrt(T[i] / Q[i][i]);
	L[i] = long(floor(Z - U[i]));
	x[i] = long(ceil(-Z - U[i]) - 1);
      }

    } while (i > 1);


    // 4. Solution found (This happens when i=1)
    /*
    cout << " x = [ " << x[1] << ", " << x[2] << ", " << x[3] << " ]" << endl;
    cout << " Q_val = Q(x) = " << Q_val << endl;
    */
    Q_val_double = precision - T[1] + Q[1][1] * (x[1] + U[1]) * (x[1] + U[1]);
    Q_val = (unsigned long) round(Q_val_double);

    //    cout << " Float = " << Q_val_double << "   Long = " << Q_val << "  XX " << endl;
    /*
    cout << " The float value is " << Q_val_double << endl;
    cout << " The associated long value is " << Q_val << endl;
    cout << endl;
    */

    if (Q_val <= precision) {
      theta[Q_val] = theta[Q_val] + 2;
    }



    // 5. Check if x = 0, for exit condition. =)
    long j=1;
    done_flag = true;
    while (j<=n) {
      if (x[j] != 0)
	done_flag = false;
      j++;
    }
  }


  // Set the value: theta[0] = 1
  theta[0] = 1;

  cout << " The precision of theta is: " << theta.Precision() << endl;

  /*
  cout << "Leaving ComputeTheta" << endl << endl;
  */

  /*
  // DIAGNOSTIC
  for (long rr=0; rr<=10; rr++)
    cout << " theta[" << rr << "] = " << theta[rr] << endl;
  */


  // Return the series
  return theta;

  /*
  // DIAGNOSTIC:
   cout << " The last two longs are: " << endl;
  cout << "  i = " << ((precision() >> 5) - 1) << " theta[i] = " << _theta[(precision() >> 5) - 1] << endl;
  cout << "  i = " << ((precision() >> 5) + 0) << " theta[i] = " << _theta[(precision() >> 5) + 0] << endl;
  cout << endl;
  */

}







// Define the << operator for the Matrix_mpz type
ostream & operator<<(ostream & out, const Matrix_mpz & matr) {    // Why doesn't "(...., Matrix_mpz & matr)" work?  =|
  matr.Print(out);
  return out;
}





// ================================ Strict I/O Front-end Routines ==========================================


///////////////////////////////////////////////
// Writes the local conditions to the ostream
///////////////////////////////////////////////
void FileOutput(const Matrix_mpz & elt, ostream & out) {

  elt._FileOutput(out);

}


////////////////////////////////////////////////
// Reads the local conditions from the istream
////////////////////////////////////////////////
void FileInput(Matrix_mpz & elt, istream & in) {

  elt._FileInput(in);

}

