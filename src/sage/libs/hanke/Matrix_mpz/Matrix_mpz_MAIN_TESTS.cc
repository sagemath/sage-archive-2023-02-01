

#include "Matrix_mpz.h"
//#include "Matrix_mpz.cc"

#include <gmp.h>
#include <gmpxx.h>

#include <vector>
#include <iostream>

using namespace std;



// Define the << operator for a valarray<T> type
template<class T>
ostream & operator<<(ostream & out, const valarray<T> & v) {
  out << "[ ";
  for(long i=0; i < v.size() - 1; i++)
    out << v[i] << ", ";
  out << v[v.size() - 1] << " ]";

  return out;
}


// Define the << operator for a vector<T> type
template<class T>
ostream & operator<<(ostream & out, const vector<T> & v) {
  out << "[ ";
  for(long i=0; i < v.size() - 1; i++)
    out << v[i] << ", ";
  out << v[v.size() - 1] << " ]";

  return out;
}




int main() {


  // --------------------------------------------------------------------------------

  // Testing:
  //    Matrix_mpz(long a, long b) Constructor
  //    Transpose()
  //    GetTranspose();


  // Make a 2 x 3 matrix with the "sized" constructor
  Matrix_mpz MM(2,3);

  MM(1,1) = 1;
  MM(1,2) = 2;
  MM(1,3) = 3;
  MM(2,1) = 4;
  MM(2,2) = 5;
  MM(2,3) = 6;

  /*
  // This gives linking errors...
  cout << "Made matrix " << endl << MM << endl;
  */

  /*
  // So does this...
  cout << "Made matrix " << endl;
  MM.Print(cout);
  cout << endl;
  */

  mpz_class x = 0;
  cout << " This is an mpz_class " << x << " for the linker... " << endl;

  // /*
  // Now it's ok! =)
  cout << "Made matrix MM " << endl << MM;
  cout << " It has " << MM.NumRows() << " rows and " << MM.NumCols() << " columns. " << endl << endl;
  //   */


  // First do Transpose() (which modifies the matrix)
  MM.Transpose();
  cout << " This is the matrix MM transposed: " << endl << MM;
  cout << " It has " << MM.NumRows() << " rows and " << MM.NumCols() << " columns. " << endl << endl;


  // Now do GetTranspose() (which copies the matrix)
  Matrix_mpz NN;
  NN = MM.GetTranspose();
  cout << " This is the matrix MM re-transposed (now called NN): " << endl << NN;
  cout << " It has " << NN.NumRows() << " rows and " << NN.NumCols() << " columns. " << endl << endl;

  // Trying it again...
  cout << " Here it is again: " << endl << MM.GetTranspose() << endl << endl;



  // ---------------------------------------------------------------
  cout << " ------------------------------------------------- " << endl;

  // Testing:
  //   Matrix_mpz() Constructor
  //   operator==


  // Make a 2 x 3 matrix with the "unsized" constructor
  Matrix_mpz MM1;
  MM1.SetDims(2,3);

  MM1(1,1) = 1;
  MM1(1,2) = 2;
  MM1(1,3) = 3;
  MM1(2,1) = 4;
  MM1(2,2) = 5;
  MM1(2,3) = 6;

  cout << "Made matrix " << endl << MM1 << endl << endl;

  // Is this the same as the "sized" construction?
  cout << " Are the two starting matrices the same? " << endl;
  cout << " Ans: " << (NN == MM1) << endl;


  // ----------------------------------------------------------------
  cout << " ------------------------------------------------- " << endl;

  // Testing
  //   NumRows()
  //   NumCols()
  //   operator(a,b)


  cout << " This matrix NN has " << NN.NumRows() << " rows and " << NN.NumCols() << " columns. " << endl;


  cout << " Here is the matrix using NN(i,j):" << endl;
  cout << NN(1,1) << " "
       << NN(1,2) << " "
       << NN(1,3) << " " << endl;
  cout << NN(2,1) << " "
       << NN(2,2) << " "
       << NN(2,3) << " " << endl << endl;



  // ---------------------------------------------------------------
  cout << " ------------------------------------------------- " << endl;

  // Testing
  //   operator %

  Matrix_mpz M1(2,3);
  M1(1,1) = 17;
  M1(1,2) = -17;
  M1(1,3) = 100;
  M1(2,1) = 0;
  M1(2,2) = -3;
  M1(2,3) = -5;

  mpz_class m;

  cout << " Started with the matrix M1:" << endl;
  cout << M1 << endl;

  m = 5;
  cout << " Modulo " << m << " this gives:" << endl;
  cout << (M1 % m) << endl;

  /*
  // This gives an error like it's supposed to! =)
  m = -5;
  cout << " Modulo " << m << " this gives:" << endl;
  cout << (M1 % m) << endl;
  */



  // ---------------------------------------------------------------
  cout << " ------------------------------------------------- " << endl;

  // Testing
  //   EvaluateQuadratic()   <-- Should change the name to EvaluateQuadraticMod()...

  Matrix_mpz M2(3,3);
  M2(1,1) = 17;
  M2(1,2) = -17;
  M2(1,3) = 100;
  M2(2,1) = 0;
  M2(2,2) = -3;
  M2(2,3) = -5;
  M2(3,1) = 1;
  M2(3,2) = 2;
  M2(3,3) = -3;

  mpz_class m2;
  m2 = 5;
  valarray<mpz_class> v2;
  v2.resize(3);
  v2[0] = 1;
  v2[1] = 0;
  v2[2] = 0;


  cout << " Started with the matrix M2:" << endl;
  cout << M2 << endl;

  cout << " Evaluated at the valarray " << v2 << " Modulo " << m2
       << ", which gives: " << M2.EvaluateQuadratic(v2, m2) << endl;


  /*
  // Correctly aborts with a modulus error
  m2 = -5;
  cout << " Evaluated at the valarray " << v2 << " Modulo " << m2
       << ", which gives: " << M2.EvaluateQuadratic(v2, m2) << endl;
  */


  /*
  // Correctly aborts with a valarray size error
  v2.resize(2);
  cout << " Evaluated at the valarray " << v2 << " Modulo " << m2
       << ", which gives: " << M2.EvaluateQuadratic(v2, m2) << endl;
  */


  /*
  // Correctly aborts with a non-square matrix error
  M2.SetDims(3,2);
  cout << " Evaluated at the valarray " << v2 << " Modulo " << m2
       << ", which gives: " << M2.EvaluateQuadratic(v2, m2) << endl;
  */



  // ---------------------------------------------------------------
  cout << " ------------------------------------------------- " << endl;

  // Testing
  //   EvaluateQuadratic()   (for vector<mpz_class, long> and valarray<mpz_class>

  Matrix_mpz M3(3,3);
  M3(1,1) = 17;
  M3(1,2) = -17;
  M3(1,3) = 100;
  M3(2,1) = 0;
  M3(2,2) = -3;
  M3(2,3) = -5;
  M3(3,1) = 1;
  M3(3,2) = 2;
  M3(3,3) = -3;

  valarray<mpz_class> v3;
  v3.resize(3);
  v3[0] = 1;
  v3[1] = 0;
  v3[2] = 0;

  vector<mpz_class> v4;
  v4.resize(3);
  v4[0] = 1;
  v4[1] = 0;
  v4[2] = 0;


  vector<long> v5;
  v5.resize(3);
  v5[0] = 1;
  v5[1] = 0;
  v5[2] = 0;



  cout << " Started with the matrix M3:" << endl;
  cout << M3 << endl;

  cout << " Evaluated at the valarray<mpz_class> " << v3
       << " gives: " << M3.EvaluateQuadratic(v3) << endl;

  cout << " Evaluated at the vector<mpz_class> " << v4
       << " gives: " << M3.EvaluateQuadratic(v4) << endl;

  cout << " Evaluated at the vector<long> " << v5
       << " gives: " << M3.EvaluateQuadratic(v5) << endl;



  // ---------------------------------------------------------------
  cout << " ------------------------------------------------- " << endl;

  // Testing
  //   Matrix Multiplication

  Matrix_mpz aa(2,3);
  aa(1,1) = 1;
  aa(1,2) = 2;
  aa(1,3) = 3;
  aa(2,1) = -1;
  aa(2,2) = -2;
  aa(2,3) = -3;

  Matrix_mpz bb(3,2);
  bb(1,1) = -1;
  bb(1,2) = 0;
  bb(2,1) = -2;
  bb(2,2) = 2;
  bb(3,1) = 6;
  bb(3,2) = 1;

  cout << " Starting with matrices aa: " << endl << aa
       << " and bb: " << endl << bb << endl;

  cout << " Then aa*bb is " << endl << (aa * bb) << endl;
  cout << " and bb*aa is " << endl << (bb * aa) << endl;


  /*
  // Correctly gives a "wrong size" error
  cout << " Also aa*a is " << endl << (aa * aa) << endl;
  */



  // -----------------------------------------------------------------

  // Remaining Tests To Do:
  //   - ExtractRow/Column
  //   - ExtractSquareSubmatrix(Ordered)
  //   - Det, Adj, QFLevel





  cout << " Finishing up now... =)" << endl;

  return 1;

}
