#include<iostream>

#include "gmp.h"
#include "matrix_rational_dense_linbox.h"
#include <linbox/integer.h>
#include "linbox/matrix/blas-matrix.h"
#include <linbox/matrix/matrix-domain.h>
#include "linbox/field/gmp-rational.h"
#include "linbox/blackbox/dense.h"
#include "linbox/algorithms/echelon-form.h"

using namespace LinBox;
using namespace std;

GMPRationalField QQ;
SpyInteger spy;

void new_gmp_matrix(DenseMatrix<GMPRationalField>& A, mpq_t** matrix, size_t nrows, size_t ncols) {
    size_t i, j;
    for (i=0; i < nrows; i++) {
	for (j=0; j < ncols; j++) {
	    GMPRationalField::Element t;
	    t = A.getEntry(i, j);
	    integer x;
	    mpz_set(spy.get_mpz(QQ.get_num(x, t)), mpq_numref(matrix[i][j]));
	    mpz_set(spy.get_mpz(QQ.get_den(x, t)), mpq_denref(matrix[i][j]));
	    A.setEntry(i, j, t);
	}
    }
}

void linbox_rational_dense_echelon_form_2(mpq_t** matrix, size_t nr, size_t nc)
{
/*    EchelonFormDomain<GMPRationalField> EF(QQ);
    DenseMatrix<GMPRationalField> A(QQ,nr, nc);
    DenseMatrix<GMPRationalField> E(QQ,nr, nc);
    new_gmp_matrix(A, matrix, nr, nc);
    cout << "made matrix\n";
    EF.rowEchelon(E, A); */
}

#include "linbox/field/givaro-rational.h"
#include "linbox/field/modular.h"

typedef Modular<int> ModInt;
ModInt F(32771);

void linbox_rational_dense_echelon_form(mpq_t** matrix, size_t nr, size_t nc)
{
    EchelonFormDomain< ModInt > EF(F);
    DenseMatrix<ModInt> A(F,nr, nc);
    DenseMatrix<ModInt> E(F,nr, nc);
    size_t k;
    k = 19;
    for (size_t i=0; i < nr; i++)
	for (size_t j=0; j < nc; j++) {
	    A.setEntry(i,j,5+i*i+j-j*j+i*i*i+k*k);
	    k += i*i + j*j + 17;
	}
    // new_gmp_matrix(A, matrix, nr, nc);
    cout << "made matrix\n";
    EF.rowEchelon(E, A);
    cout << "did echelon\n";
}
