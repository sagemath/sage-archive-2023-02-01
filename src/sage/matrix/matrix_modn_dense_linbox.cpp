#include "matrix_modn_dense_linbox.h"


#include<iostream>
#include "gmp.h"
#include <linbox/integer.h>
#include "linbox/matrix/blas-matrix.h"
#include <linbox/matrix/matrix-domain.h>
#include "linbox/field/gmp-rational.h"
#include "linbox/blackbox/dense.h"
#include "linbox/algorithms/echelon-form.h"

using namespace LinBox;
using namespace std;

#include "linbox/field/modular.h"

/* NOTE: There are many echelon form functions, possible base rings, etc.  Strangely,
   most don't build.  This combination below does though.
*/

typedef Modular<int> ModInt;

int linbox_matrix_modn_dense_echelonize(unsigned long modulus,
					unsigned long** matrix, size_t nrows, size_t ncols) {
    ModInt F(modulus);
    EchelonFormDomain< ModInt > EF(F);
    DenseMatrix<ModInt> A(F, nrows, ncols);
    DenseMatrix<ModInt> E(F, nrows, ncols);

    unsigned long* row;
    for (size_t i=0; i < nrows; i++) {
	row = matrix[i];
	for (size_t j=0; j < ncols; j++)
	    A.setEntry(i, j, row[j]);
	}
    int rank = EF.rowReducedEchelon(E, A);
    for (size_t i=0; i < nrows; i++) {
	row = matrix[i];
	for (size_t j=0; j < ncols; j++)
	    row[j] = E.getEntry(i, j);
	}
    return rank;
}


