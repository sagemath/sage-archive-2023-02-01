#include "gmp.h"

#include<iostream>

#include "matrix_integer_dense_linbox.h"

#include <linbox/integer.h>
#include "linbox/matrix/blas-matrix.h"
#include <linbox/matrix/matrix-domain.h>

#include "linbox/field/modular-double.h"
#include "linbox/field/gmp-integers.h"

#include "linbox/blackbox/dense.h"
#include "linbox/solutions/charpoly.h"
#include "linbox/solutions/minpoly.h"
#include "linbox/solutions/rank.h"

using namespace LinBox;
using namespace std;


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

void set_matrix(mpz_t** matrix, DenseMatrix<GMP_Integers>& A, size_t nrows, size_t ncols) {
    size_t i, j, k;
    for (i=0; i < nrows; i++) {
	for (j=0; j < ncols; j++) {
	    mpz_set(matrix[i][j], spy.get_mpz(A.getEntry(i,j)));
	}
    }
}

void linbox_integer_dense_minpoly(mpz_t* *mp, size_t* degree, size_t n, mpz_t** matrix, int do_minpoly) {
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

    vector<GMP_Integers::Element> m_A;

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

/* broken when n % 4 == 0 */
void linbox_integer_dense_charpoly(mpz_t* *mp, size_t* degree, size_t n, mpz_t** matrix) {
    DenseMatrix<GMP_Integers> A(new_matrix(matrix, n, n));
    vector<GMP_Integers::Element> m_A;
    charpoly(m_A, A);

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
    DenseMatrix<GMP_Integers> AA(new_matrix(A, A_nr, A_nc));
    DenseMatrix<GMP_Integers> BB(new_matrix(B, B_nr, B_nc));
    if (A_nc != B_nr)
	return -1;   // error
    DenseMatrix<GMP_Integers> CC(ZZ, A_nr, B_nc);

    MatrixDomain<GMP_Integers> MD(ZZ);

    MD.mul(CC, AA, BB);

    set_matrix(ans, CC, A_nr, B_nc);

    return 0;
}
