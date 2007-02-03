#include "gmp.h"

#include<iostream>

#include "linbox_wrap.h"

#include "linbox/field/modular-double.h"
#include "linbox/field/gmp-integers.h"
#include "linbox/blackbox/dense.h"
#include "linbox/solutions/charpoly.h"
#include "linbox/solutions/minpoly.h"

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

mpz_t* linbox_minpoly(size_t n, mpz_t** matrix) {
    std::cout << "Starting linbox minpoly\n";
    GMP_Integers ZZ;
    DenseMatrix<GMP_Integers> A (ZZ, n, n);
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

    cout << "computed minpoly" << endl;

    cout << "now converting to list" << endl;
    printPolynomial (ZZ, m_A);

    mpz_t* f = new mpz_t[m_A.size()];
    cout << m_A.size() << endl;
    for (i=0; i < m_A.size(); i++) {
	mpz_init(f[i]);
	mpz_set(f[i], spy.get_mpz(m_A[i]));
    }

    return f;
}

