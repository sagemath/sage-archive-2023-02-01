#include <iostream>
#include <gmp.h>
#include "linbox_wrap.h"




/*************************************************************************
   dense modulo Z/nZ
*************************************************************************/

#include <linbox/integer.h>
#include <linbox/matrix/blas-matrix.h>
#include <linbox/matrix/matrix-domain.h>
#include <linbox/field/gmp-rational.h>
#include <linbox/blackbox/dense.h>
#include <linbox/algorithms/echelon-form.h>
#include <linbox/solutions/minpoly.h>
#include <linbox/solutions/charpoly.h>
#include <linbox/ring/givaro-polynomial.h>
#include <linbox/field/modular.h>

using namespace LinBox;
using namespace std;

#include "linbox/util/commentator.h"


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
    DenseMatrix<ModInt> A(F, nrows, ncols);
    DenseMatrix<ModInt> E(F, nrows, ncols);

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
    /* We program around a bizarre bug in linbox, where minpoly doesn't work
       on matrices that are n x n with n divisible by 4!
    */

    ModInt F((double)modulus);

    size_t m = n;
//     if (n % 4 == 0 || !do_minpoly) {
//         m = n + 1;
//     } else {
// 	m = n;
//     }

    DenseMatrix<ModInt> A(linbox_new_modn_matrix( modulus, matrix, m, m));

    GivPolynomial<ModInt::Element> m_A;

    if (do_minpoly)
	minpoly(m_A, A);
    else
        charpoly(m_A, A);


//     if (n%4 == 0 || !do_minpoly) {
//         /* Program around the bug.
// 	   It is OK that this code is crappy and redundant, since it will get replaced
//            when linbox gets fixed. */
// 	int divide_by_x;

// 	if (!do_minpoly)
// 	    divide_by_x = 1;
// 	else {
// 	    long unsigned int r;
// 	    rank(r, A);
// 	    divide_by_x = (r==n);
// 	}
// 	if (divide_by_x) {
// 	    /* x was not a factor of the charpoly after all. */
// 	    (*mp) = new mod_int[m_A.size()-1];
// 	    *degree = m_A.size() - 2;
// 	    for (size_t i=0; i <= *degree; i++) {
// 		(*mp)[i] = m_A[i+1];
// 	    }
// 	    return;
// 	}
//     }

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

    DenseMatrix<ModInt> A (F, nrows, ncols);

    size_t i, j, k;
    for (i=0; i < nrows; i++) {
	for (j=0; j < ncols; j++) {
	    A.setEntry(i, j, (double)matrix[i][j]);
	}
    }
    return A;
};

static void linbox_set_modn_matrix(mod_int** matrix, DenseMatrix<ModInt>& A, size_t nrows, size_t ncols) {
    size_t i, j, k;
    for (i=0; i < nrows; i++) {
	for (j=0; j < ncols; j++) {
	    matrix[i][j] = (mod_int)A.getEntry(i,j);
	}
    }
};

/*************************************************************************
    Modular<int> versions of everything: much faster for this.
**********************************************************************/

typedef Modular<int> Mod_int;

static DenseMatrix<Mod_int> linbox_new_modn_matrix2(mod_int modulus, mod_int** matrix, size_t nrows, size_t ncols) {

    Mod_int F((double)modulus);

    DenseMatrix<Mod_int> A (F, nrows, ncols);

    size_t i, j, k;
    for (i=0; i < nrows; i++) {
	for (j=0; j < ncols; j++) {
	    A.setEntry(i, j, matrix[i][j]);
	}
    }
    return A;
};

static void linbox_set_modn_matrix2(mod_int** matrix, DenseMatrix<Mod_int>& A, size_t nrows, size_t ncols) {
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

    Mod_int F(modulus);

    DenseMatrix<Mod_int> AA(linbox_new_modn_matrix2(modulus, A, A_nr, A_nc));
    DenseMatrix<Mod_int> BB(linbox_new_modn_matrix2(modulus, B, B_nr, B_nc));
    if (A_nc != B_nr)
	return -1;   // error
    DenseMatrix<Mod_int> CC(F, A_nr, B_nc);

    MatrixDomain<Mod_int> MD(F);

    MD.mul(CC, AA, BB);

    linbox_set_modn_matrix2(ans, CC, A_nr, B_nc);

    return 0;
}




/*************************************************************************
   sparse modulo Z/nZ
*************************************************************************/


int linbox_modn_sparse_rank(mod_int modulus,
			   mod_int** matrix, size_t nrows, size_t ncols) {

  /*    typedef Modular<double> Field;
    Field F(modulus);
    SparseMatrix<Field, Vector<Field>::SparseSeq > B (ms);

    mod_int* row;
    for (size_t i=0; i < nrows; i++) {
	row = matrix[i];
	for (size_t j=0; j < ncols; j++)
	    A.setEntry(i, j, (double)row[j]);
	}

    unsigned long r;
    rank(r, A);
    return r;
  */
}



/*************************************************************************
    dense over ZZ
*************************************************************************/

#include "gmp.h"

#include<iostream>

#include <linbox/integer.h>
#include <linbox/matrix/matrix-domain.h>
#include "linbox/field/gmp-integers.h"
#include "linbox/solutions/rank.h"
#include "linbox/solutions/minpoly.h"
#include "linbox/solutions/charpoly.h"
#include "linbox/ring/givaro-polynomial.h"
#include "linbox/element/givaro-polynomial.h"

using namespace LinBox;
using namespace std;

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
     DenseMatrix<GMP_Integers> A (ZZ, nrows, ncols);

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
     DenseMatrix<Integers> A (ZZ, nrows, ncols);

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

template<class Field>
void set_matrix(mpz_t** matrix, DenseMatrix<Field>& A, size_t nrows, size_t ncols) {
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

     DenseMatrix<GMP_Integers> A(ZZ, m, m);

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

     DenseMatrix<Integers> AA(new_matrix_integers(A, A_nr, A_nc));
     DenseMatrix<Integers> BB(new_matrix_integers(B, B_nr, B_nc));
     if (A_nc != B_nr)
	 return -1;   // error
     DenseMatrix<Integers> CC(ZZ, A_nr, B_nc);

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


#include "linbox/algorithms/smith-form-adaptive.h"

DenseMatrix<NTL_ZZ> new_matrix_integer_dense_ntl(mpz_t** matrix, size_t nrows, size_t ncols) {
     NTL_ZZ ZZ;
     DenseMatrix<NTL_ZZ> A (ZZ, nrows, ncols);
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
