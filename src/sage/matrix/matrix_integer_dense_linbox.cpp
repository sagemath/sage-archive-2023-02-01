#include "gmp.h"

#include<iostream>

#include "matrix_integer_dense_linbox.h"

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

void linbox_minpoly(mpz_t* *mp, size_t* degree, size_t n, mpz_t** matrix) {
    GMP_Integers ZZ;

    /* We program around a bizarre bug in linbox, where minpoly doesn't work
       on matrices that are n x n with n divisible by 4!
    */
    size_t m;
    if (n % 4 == 0) {
        m = n + 1;
    } else {
	m = n;
    }

    DenseMatrix<GMP_Integers> A (ZZ, m, m);
    size_t i, j, k;
    SpyInteger spy;
    for (i=0; i < n; i++) {
	for (j=0; j < n; j++) {
	    GMP_Integers::Element t;
	    mpz_set(spy.get_mpz(t), matrix[i][j]);
	    A.setEntry(i, j, t);
	}
    }

    vector<GMP_Integers::Element> m_A;

    minpoly(m_A, A);
    //charpoly(m_A, A);

    if (n%4 == 0) {
        /* Program around the bug.
	   It is OK that this code is crappy and redundant, since it will get replaced
           when linbox gets fixed. */
	long unsigned int r;
	rank(r, A);
	if (r == n) {
	    /* x was not a factor of the charpoly after all. */
	    (*mp) = new mpz_t[m_A.size()-1];
	    *degree = m_A.size() - 2;
	    for (i=0; i <= *degree; i++) {
		mpz_init((*mp)[i]);
		mpz_set((*mp)[i], spy.get_mpz(m_A[i+1]));
	    }
	    return;
	}
    }

    (*mp) = new mpz_t[m_A.size()];
    *degree = m_A.size() - 1;
    for (i=0; i <= *degree; i++) {
	mpz_init((*mp)[i]);
	mpz_set((*mp)[i], spy.get_mpz(m_A[i]));
    }

}

void linbox_delete_array(mpz_t* f) {
    delete[] f;
}

