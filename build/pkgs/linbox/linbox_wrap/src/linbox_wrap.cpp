#include <iostream>
#include <gmp.h>

#include <cstdlib>
#include <vector>

#include "linbox_wrap.h"

#include <linbox/util/commentator.h>

#include <linbox/blackbox/sparse.h>

#include "linbox/element/givaro-polynomial.h"

#include <linbox/matrix/blas-matrix.h>
#include <linbox/matrix/sparse.h>
#include <linbox/vector/sparse.h>

#include <linbox/algorithms/blas-domain.h>
#include <linbox/algorithms/echelon-form.h>
#include "linbox/algorithms/gauss.h"
#include "linbox/algorithms/smith-form-adaptive.h"

#include <linbox/solutions/rank.h>
#include <linbox/solutions/det.h>
#include <linbox/solutions/solve.h>
#include "linbox/solutions/methods.h"
#include <linbox/solutions/minpoly.h>
#include <linbox/solutions/charpoly.h>

#include <linbox/integer.h>
#include <linbox/field/gmp-integers.h>
#include <linbox/field/gmp-rational.h>
#include <linbox/ring/givaro-polynomial.h>
#include <linbox/field/modular.h>

using namespace LinBox;
using namespace std;

/*************************************************************************
   dense modulo Z/nZ
*************************************************************************/

//we are using Modular<double> here as it seems to be best supported
typedef Modular<double> ModInt;
typedef GivPolynomialRing<ModInt::Element,Dense> ModIntPolRing;

static DenseMatrix<ModInt> linbox_new_modn_matrix(mod_int modulus, mod_int** matrix, size_t nrows, size_t ncols);
static void linbox_set_modn_matrix(mod_int** matrix, DenseMatrix<ModInt>& A, size_t nrows, size_t ncols);


/* NOTE: There are many echelon form functions, possible base rings, etc.  Strangely,
   most don't build.  This combination below does though.
*/
int linbox_modn_dense_echelonize(mod_int modulus,
				 mod_int** matrix, size_t nrows, size_t ncols) {


    //this is the way Clement suggested, need to figure out if it is way faster for very
    //large inputs (approx. 2 seconds slower for 5000x5000 on my system) (malb)
    //     typedef Modular<double> Field;
    //     Field F((double)modulus);
    //     BlasMatrix<Field::Element> A(nrows, ncols);
    //     BlasMatrix<Field::Element> E(A.rowdim(),A.coldim());
    //     EchelonFormDomain<Modular<double> > EF (F);

    ModInt F((double)modulus);
    EchelonFormDomain< ModInt > EF(F);
    BlasMatrix<ModInt::Element> A( nrows, ncols);
    BlasMatrix<ModInt::Element> E( nrows, ncols);

    mod_int* row;
    for (size_t i=0; i < nrows; i++) {
	row = matrix[i];
	for (size_t j=0; j < ncols; j++)
	    A.setEntry(i, j, (double)row[j]);
	}
    int rank = EF.rowReducedEchelon(E, A);
    for (size_t i=0; i < nrows; i++) {
	row = matrix[i];
	for (size_t j=0; j < ncols; j++)
	    row[j] = (mod_int)E.getEntry(i, j);
	}
    return rank;
}

int linbox_modn_dense_rank(mod_int modulus,
			   mod_int** matrix, size_t nrows, size_t ncols) {

    ModInt F((double)modulus);
    EchelonFormDomain< ModInt > EF(F);
    DenseMatrix<ModInt> A(F, nrows, ncols);

    mod_int* row;
    for (size_t i=0; i < nrows; i++) {
	row = matrix[i];
	for (size_t j=0; j < ncols; j++)
	    A.setEntry(i, j, (double)row[j]);
	}

    unsigned long r;
    rank(r, A);
    return r;
}


void linbox_modn_dense_minpoly(mod_int modulus, mod_int **mp, size_t* degree, size_t n, mod_int **matrix, int do_minpoly) {

    ModInt F((double)modulus);

    size_t m = n;

    DenseMatrix<ModInt> A(linbox_new_modn_matrix( modulus, matrix, m, m));

    GivPolynomial<ModInt::Element> m_A;

    if (do_minpoly)
	minpoly(m_A, A);
    else
        charpoly(m_A, A);

    (*mp) = new mod_int[m_A.size()];
    *degree = m_A.size() - 1;
    for (size_t i=0; i <= *degree; i++) {
	(*mp)[i] = (mod_int)m_A[i];
    }

}

void linbox_modn_dense_delete_array(mod_int *f) {
    delete[] f;
}

static DenseMatrix<ModInt> linbox_new_modn_matrix(mod_int modulus, mod_int** matrix, size_t nrows, size_t ncols) {

    ModInt F((double)modulus);

    DenseMatrix<ModInt> A ( F, nrows, ncols);

    size_t i, j, k;
    for (i=0; i < nrows; i++) {
	for (j=0; j < ncols; j++) {
	    A.setEntry(i, j, (double)matrix[i][j]);
	}
    }
    return A;
};

static void linbox_set_modn_matrix(mod_int** matrix, BlasMatrix<ModInt::Element>& A, size_t nrows, size_t ncols) {
    size_t i, j, k;
    for (i=0; i < nrows; i++) {
	for (j=0; j < ncols; j++) {
	    matrix[i][j] = (mod_int)A.getEntry(i,j);
	}
    }
};

int linbox_modn_dense_matrix_matrix_multiply(mod_int modulus, mod_int **ans, mod_int **A, mod_int **B,
					     size_t A_nr, size_t A_nc, size_t B_nr, size_t B_nc)
{

    ModInt F((double)modulus);

    BlasMatrix<ModInt::Element> AA(linbox_new_modn_matrix(modulus, A, A_nr, A_nc));
    BlasMatrix<ModInt::Element> BB(linbox_new_modn_matrix(modulus, B, B_nr, B_nc));
    if (A_nc != B_nr)
	return -1;   // error
    BlasMatrix<ModInt::Element> CC( A_nr, B_nc);

    BlasMatrixDomain<ModInt> MD(F);

    MD.mul(CC, AA, BB);

    linbox_set_modn_matrix(ans, CC, A_nr, B_nc);

    return 0;
}

/*************************************************************************
    dense over ZZ
*************************************************************************/

typedef PID_integer Integers;

template <class Field, class Polynomial>
void printPolynomial (const Field &F, const Polynomial &v)
{
        for (int i = v.size () - 1; i >= 0; i--) {
                F.write (cout, v[i]);
                if (i > 0)
                        cout << " x^" << i << " + ";
        }
        cout << endl;
}

GMP_Integers ZZ;
SpyInteger spy;
typedef GivPolynomialRing<GMP_Integers,Dense> IntPolRing;

DenseMatrix<GMP_Integers> new_matrix(mpz_t** matrix, size_t nrows, size_t ncols) {
  DenseMatrix<GMP_Integers> A ( ZZ, nrows, ncols);

     size_t i, j, k;
     for (i=0; i < nrows; i++) {
	 for (j=0; j < ncols; j++) {
	     GMP_Integers::Element t;
	     mpz_set(spy.get_mpz(t), matrix[i][j]);
	     A.setEntry(i, j, t);
	 }
     }
     return A;
}

DenseMatrix<Integers> new_matrix_integers(mpz_t** matrix, size_t nrows, size_t ncols) {
     Integers ZZ;
     DenseMatrix<Integers> A ( ZZ,nrows, ncols);

     size_t i, j, k;
     for (i=0; i < nrows; i++) {
	 for (j=0; j < ncols; j++) {
	     Integers::Element t;
	     mpz_set(spy.get_mpz(t), matrix[i][j]);
	     A.setEntry(i, j, t);
	 }
     }
     return A;
}

template<class Element>
void set_matrix(mpz_t** matrix, BlasMatrix<Element>& A, size_t nrows, size_t ncols) {
     size_t i, j, k;
     for (i=0; i < nrows; i++) {
	 for (j=0; j < ncols; j++) {
	     mpz_set(matrix[i][j], spy.get_mpz(A.getEntry(i,j)));
	 }
     }
}

void linbox_integer_dense_minpoly_hacked(mpz_t* *mp, size_t* degree, size_t n, mpz_t** matrix, int do_minpoly) {
     /* We program around a bizarre bug in linbox, where minpoly doesn't work
	on matrices that are n x n with n divisible by 4!
     */
     size_t m;
     if (n % 4 == 0 || !do_minpoly) {
	 m = n + 1;
     } else {
	 m = n;
     }

     DenseMatrix<GMP_Integers> A( ZZ, m, m);

     size_t i, j;
     GMP_Integers::Element t;
     for (i=0; i < n; i++) {
	 for (j=0; j < n; j++) {
	     mpz_set(spy.get_mpz(t), matrix[i][j]);
	     A.setEntry(i, j, t);
	 }
     }

 //    vector<GMP_Integers::Element> m_A;
     IntPolRing::Element m_A;

     if (do_minpoly)
	 minpoly(m_A, A);
     else
	 charpoly(m_A, A);

     if (n%4 == 0 || !do_minpoly) {
	 /* Program around the bug.
	    It is OK that this code is crappy and redundant, since it will get replaced
	    when linbox gets fixed. */
	 int divide_by_x;

	 if (!do_minpoly)
	     divide_by_x = 1;
	 else {
	     long unsigned int r;
	     rank(r, A);
	     divide_by_x = (r==n);
	 }
	 if (divide_by_x) {
	     /* x was not a factor of the charpoly after all. */
	     (*mp) = new mpz_t[m_A.size()-1];
	     *degree = m_A.size() - 2;
	     for (size_t i=0; i <= *degree; i++) {
		 mpz_init((*mp)[i]);
		 mpz_set((*mp)[i], spy.get_mpz(m_A[i+1]));
	     }
	     return;
	 }
     }

     (*mp) = new mpz_t[m_A.size()];
     *degree = m_A.size() - 1;
     for (size_t i=0; i <= *degree; i++) {
	 mpz_init((*mp)[i]);
	 mpz_set((*mp)[i], spy.get_mpz(m_A[i]));
     }

}

 void linbox_integer_dense_charpoly(mpz_t* *mp, size_t* degree, size_t n, mpz_t** matrix) {
     /* THIS IS Broken when n % 4 == 0!!!!  Use above function instead. */
   /*    linbox_integer_dense_minpoly(mp, degree, n, matrix, 0); */

     DenseMatrix<GMP_Integers> A(new_matrix(matrix, n, n));
     IntPolRing::Element m_A;
     charpoly(m_A, A);

     (*mp) = new mpz_t[m_A.size()];
     *degree = m_A.size() - 1;
     for (size_t i=0; i <= *degree; i++) {
	 mpz_init((*mp)[i]);
	 mpz_set((*mp)[i], spy.get_mpz(m_A[i]));
     }

 }

 void linbox_integer_dense_minpoly(mpz_t* *mp, size_t* degree, size_t n, mpz_t** matrix) {
     /* THIS IS Broken when n % 4 == 0!!!!  Use above function instead. */
   /*    linbox_integer_dense_minpoly(mp, degree, n, matrix, 0); */

   DenseMatrix<GMP_Integers> A(new_matrix(matrix, n, n));
     IntPolRing::Element m_A;
     minpoly(m_A, A);

     (*mp) = new mpz_t[m_A.size()];
     *degree = m_A.size() - 1;
     for (size_t i=0; i <= *degree; i++) {
	 mpz_init((*mp)[i]);
	 mpz_set((*mp)[i], spy.get_mpz(m_A[i]));
     }

 }

 void linbox_integer_dense_delete_array(mpz_t* f) {
     delete[] f;
 }

 int linbox_integer_dense_matrix_matrix_multiply(mpz_t** ans, mpz_t **A, mpz_t **B,
				   size_t A_nr, size_t A_nc, size_t B_nr, size_t B_nc)
 {
   typedef PID_integer Integers;
   Integers ZZ;

     BlasMatrix<Integers::Element> AA(new_matrix_integers(A, A_nr, A_nc));
     BlasMatrix<Integers::Element> BB(new_matrix_integers(B, B_nr, B_nc));
     if (A_nc != B_nr)
	 return -1;   // error
     BlasMatrix<Integers::Element> CC( A_nr, B_nc);

     MatrixDomain<Integers> MD(ZZ);

     MD.mul(CC, AA, BB);

     set_matrix(ans, CC, A_nr, B_nc);

     return 0;
 }

 unsigned long linbox_integer_dense_rank(mpz_t** matrix, size_t nrows,
					 size_t ncols) {
     DenseMatrix<GMP_Integers> A(new_matrix(matrix, nrows, ncols));
     unsigned long r;
     rank(r, A);
     return r;
 }

 void linbox_integer_dense_det(mpz_t ans, mpz_t** matrix, size_t nrows,
			       size_t ncols) {
   commentator.setMaxDetailLevel(0);
   commentator.setMaxDepth (0);

   DenseMatrix<Integers> A(new_matrix_integers(matrix, nrows, ncols));
   GMP_Integers::Element d;
   det(d, A);
   mpz_set(ans, spy.get_mpz(d));
}


DenseMatrix<NTL_ZZ> new_matrix_integer_dense_ntl(mpz_t** matrix, size_t nrows, size_t ncols) {
     NTL_ZZ ZZ;
     DenseMatrix<NTL_ZZ> A (ZZ,nrows, ncols);
     size_t i, j, k;
     for (i=0; i < nrows; i++) {
	 for (j=0; j < ncols; j++) {
	   NTL_ZZ::Element t;
	   PID_integer::Element s;
	   mpz_set(spy.get_mpz(s), matrix[i][j]);
	   ZZ.init(t, s);
	   A.setEntry(i, j, t);
	 }
     }
     return A;
}

/*
This won't build on OS X PPC.

void linbox_integer_dense_smithform(mpz_t **v,
                                    mpz_t **matrix, size_t nrows, size_t ncols) {
  typedef NTL_ZZ Ints;
  Ints Z;
  DenseMatrix<Ints> M(new_matrix_integer_dense_ntl(matrix, nrows, ncols));
  vector<integer> w(ncols);
  SmithFormAdaptive::smithForm(w, M);

  (*v) = new mpz_t[ncols];
  for (size_t i=0; i < ncols; i++) {
    mpz_init((*v)[i]);
    mpz_set((*v)[i], spy.get_mpz(w[i]));
  }
}
*/

/*************************************************************************
   sparse modulo Z/nZ
*************************************************************************/

struct c_vector_modint_linbox {
  // copy of the declaration in vector_modn_sparse.pxi
  int *entries;
  int p;
  size_t *positions;
  size_t degree;
  size_t num_nonzero;
};

typedef Modular<unsigned int> GFp;
typedef GFp::Element  Element;
typedef std::vector <pair <size_t, Element> > SparseSeqVectorGFp;
typedef SparseMatrix<GFp, SparseSeqVectorGFp> SparseMatrixGFp;

static SparseMatrixGFp linbox_new_modn_sparse_matrix(mod_int modulus, size_t numrows, size_t numcols, void *rows) {
  GFp F(modulus);
  SparseMatrixGFp M(F, numrows, numcols);

  struct c_vector_modint_linbox *A = static_cast<struct c_vector_modint_linbox *>(rows);

  for(int i = 0; i < numrows; i++) {
    for(int j = 0; j < A[i].num_nonzero; j++) {
      M.setEntry(i, A[i].positions[j], A[i].entries[j]);
    }
  }
  return M;
}

static vector<Element> linbox_new_modn_sparse_vector(mod_int modulus, size_t len, void *_vec) {
  GFp F(modulus);

  vector<GFp::Element> A(len);

  if (_vec==NULL) {
    return A;
  }

  struct c_vector_modint_linbox *vec = static_cast<struct c_vector_modint_linbox*>(_vec);
  for(int i = 0; i < vec->num_nonzero; i++) {
    A[vec->positions[i]] = vec->entries[i];
  }
  return A;
}


unsigned long linbox_modn_sparse_matrix_rank(mod_int modulus, size_t numrows, size_t numcols,  void *rows, int gauss) {
  GFp F(modulus);
  unsigned long M_rank;
  Element M_det;
  GaussDomain<GFp> dom(F);

  SparseMatrixGFp M( linbox_new_modn_sparse_matrix(modulus, numrows, numcols, rows) );

  if(!gauss) {
    dom.InPlaceLinearPivoting(M_rank, M_det, M, numrows, numcols);
  } else {
    dom.NoReordering(M_rank, M_det, M, numrows, numcols);
  }

  //*pivots = (int*)calloc(sizeof(int), dom.pivots.size());

//   int j=0;
//   for(vector<int>::const_iterator i= dom.pivots.begin(); i!= dom.pivots.end(); i++, j++){
//     (*pivots)[j] = *i;
//   }

  return M_rank;
}


vector<unsigned int> linbox_modn_sparse_matrix_solve(mod_int p, size_t numrows, size_t numcols, void *_a, void *b, int method) {
  // solve ax = b, for x, a matrix, b vector, x vector
  GFp F(p);

  vector<Element> X( numrows);
  vector<Element> B( linbox_new_modn_sparse_vector(p, numcols, b));

  SparseMatrixGFp A(linbox_new_modn_sparse_matrix(p, numrows, numcols, _a));

  switch(method) {
  case 1:
    solve(X, A, B, Method::BlasElimination());
    break;

  case 2:
    solve(X, A, B, Method::Blackbox());
    break;

  case 3:
    solve(X, A, B, Method::Wiedemann());
    break;

  default:
    solve(X, A, B);
  }
  return X;
}
