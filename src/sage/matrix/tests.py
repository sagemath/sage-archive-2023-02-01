"""
TESTS:

We test various degenerate cases of kernel computation:

    sage: matrix(ZZ,1,0).kernel()
    Ambient free module of rank 1 over the principal ideal domain Integer Ring
    sage: matrix(QQ,1,0).kernel()
    Vector space of dimension 1 over Rational Field
    sage: matrix(GF(7),1,0).kernel()
    Vector space of dimension 1 over Finite Field of size 7
    sage: matrix(Frac(QQ['x']),1,0).kernel()
    Vector space of dimension 1 over Fraction Field of Univariate Polynomial Ring in x over Rational Field

    sage: matrix(ZZ,0,1).kernel()
    Free module of degree 0 and rank 0 over Integer Ring
    Echelon basis matrix:
    []
    sage: matrix(QQ,0,1).kernel()
    Vector space of degree 0 and dimension 0 over Rational Field
    Basis matrix:
    []
    sage: matrix(GF(7),0,1).kernel()
    Vector space of degree 0 and dimension 0 over Finite Field of size 7
    Basis matrix:
    []
    sage: matrix(Frac(QQ['x']),0,1).kernel()
    Vector space of degree 0 and dimension 0 over Fraction Field of Univariate Polynomial Ring in x over Rational Field
    Basis matrix:
    []
"""
