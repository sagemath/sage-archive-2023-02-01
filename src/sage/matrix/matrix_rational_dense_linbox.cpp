#include "gmp.h"
#include <linbox/integer.h>
#include <linbox/matrix/matrix-domain.h>
#include "linbox/field/gmp-integers.h"
#include "linbox/field/gmp-rational.h"

#include "linbox/blackbox/dense.h"

#include "linbox/algorithms/echelon-form.h"

using namespace LinBox;
using namespace std;

GMPRationalField QQ;
SpyInteger spy;

DenseMatrix<GMPRationalField> new_gmp_matrix(mpq_t** matrix, size_t nrows, size_t ncols) {
    DenseMatrix<GMPRationalField> A (QQ, nrows, ncols);

    size_t i, j, k;
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
    return A;
}

void linbox_rational_dense_echelon_form(mpq_t** matrix, size_t nr, size_t nc)
{
    EchelonFormDomain<GMPRationalField> EF(QQ);
    DenseMatrix<GMPRationalField> AA(new_matrix(matrix, nr, nc));
}
