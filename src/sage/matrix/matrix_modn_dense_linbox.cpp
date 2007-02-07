#include "matrix_modn_dense_linbox.h"


#include <iostream>
#include <gmp.h>
#include <linbox/integer.h>
#include <linbox/matrix/blas-matrix.h>
#include <linbox/matrix/matrix-domain.h>
#include <linbox/field/gmp-rational.h>
#include <linbox/blackbox/dense.h>
#include <linbox/algorithms/echelon-form.h>
#include <linbox/solutions/minpoly.h>
#include <linbox/solutions/charpoly.h>
#include <linbox/ring/givaro-polynomial.h>


using namespace LinBox;
using namespace std;

#include <linbox/field/modular.h>

/** local header **/

//we are using Modular<double> here as it seems to be best supported
typedef Modular<double> ModInt;
typedef GivPolynomialRing<ModInt::Element,Dense> ModIntPolRing;


static DenseMatrix<ModInt> linbox_new_modn_matrix(unsigned long moddulus, mod_int** matrix, size_t nrows, size_t ncols);
static void linbox_set_modn_matrix(mod_int** matrix, DenseMatrix<ModInt>& A, size_t nrows, size_t ncols);


/* NOTE: There are many echelon form functions, possible base rings, etc.  Strangely,
   most don't build.  This combination below does though.
*/
int linbox_modn_dense_echelonize(unsigned long modulus,
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

    unsigned long* row;
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

int linbox_modn_dense_rank(unsigned long modulus,
			   mod_int** matrix, size_t nrows, size_t ncols) {

    ModInt F((double)modulus);
    EchelonFormDomain< ModInt > EF(F);
    DenseMatrix<ModInt> A(F, nrows, ncols);

    unsigned long* row;
    for (size_t i=0; i < nrows; i++) {
	row = matrix[i];
	for (size_t j=0; j < ncols; j++)
	    A.setEntry(i, j, (double)row[j]);
	}

    unsigned long r;
    rank(r, A);
    return r;
}


void linbox_modn_dense_minpoly(unsigned long modulus, mod_int **mp, size_t* degree, size_t n, mod_int **matrix, int do_minpoly) {
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

static DenseMatrix<ModInt> linbox_new_modn_matrix(unsigned long modulus, mod_int** matrix, size_t nrows, size_t ncols) {

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


int linbox_modn_dense_matrix_matrix_multiply(unsigned long modulus, mod_int **ans, mod_int **A, mod_int **B,
					     size_t A_nr, size_t A_nc, size_t B_nr, size_t B_nc)
{

    ModInt F((double)modulus);

    DenseMatrix<ModInt> AA(linbox_new_modn_matrix(modulus, A, A_nr, A_nc));
    DenseMatrix<ModInt> BB(linbox_new_modn_matrix(modulus, B, B_nr, B_nc));
    if (A_nc != B_nr)
	return -1;   // error
    DenseMatrix<ModInt> CC(F, A_nr, B_nc);

    MatrixDomain<ModInt> MD(F);

    MD.mul(CC, AA, BB);

    linbox_set_modn_matrix(ans, CC, A_nr, B_nc);

    return 0;
}



