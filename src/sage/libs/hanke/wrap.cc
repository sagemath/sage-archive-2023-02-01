#include <iostream>
#include <sstream>
using namespace std;
#include "wrap.h"

char* stringstream_to_char(ostringstream& instore) {
  int n = strlen(instore.str().data());
  char* buf = new char[n+1];
  strcpy(buf, instore.str().data());
  return buf;
}

char* mpz_to_char(const mpz_class& x) {
  ostringstream instore;
  instore << x;
  return stringstream_to_char(instore);
}

char* mpq_to_char(const mpq_class& x) {
  ostringstream instore;
  instore << x;
  return stringstream_to_char(instore);
}

struct Matrix_mpz* Matrix_mpz_new(int r, int s)
{
  return new Matrix_mpz(r, s);
}

void Matrix_mpz_del(struct Matrix_mpz* m)
{
  delete m;
}

char* Matrix_mpz_repr(struct Matrix_mpz* m)
{
  ostringstream instore;
  instore << (*m);
  return stringstream_to_char(instore);
}

int Matrix_mpz_nrows(struct Matrix_mpz* m)
{
  return m->NumRows();
}

int Matrix_mpz_ncols(struct Matrix_mpz* m)
{
  return m->NumCols();
}


void Matrix_mpz_setitem(struct Matrix_mpz* x, int i, int j, char* z)
{
  mpz_class y(z);
  x->set(i, j, y);
}

char* Matrix_mpz_getitem(struct Matrix_mpz* x, int i, int j)
{
  mpz_class y = x->get(i, j);
  ostringstream instore;
  instore << y;
  return stringstream_to_char(instore);
}

char* Matrix_mpz_determinant(const struct Matrix_mpz* x)
{
  mpz_class d = x->Determinant();
  ostringstream instore;
  instore << d;
  return stringstream_to_char(instore);
}

Matrix_mpz* Matrix_mpz_adjoint(const struct Matrix_mpz* x)
{
  Matrix_mpz* t = new Matrix_mpz();
  *t = x->Adjoint();
  return t;
}

char* Matrix_mpz_Local_Density(struct Matrix_mpz* x, char* p, char* m)
{
  mpz_class _p(*p), _m(*m);
  return mpq_to_char(x->Local_Density(_p, _m));
}


char* Matrix_mpz_Local_Primitive_Density(struct Matrix_mpz* x, char* p, char* m)
{
  mpz_class _p(*p), _m(*m);
  return mpq_to_char(x->Local_Primitive_Density(_p, _m));
}

char* Matrix_mpz_level(Matrix_mpz* x)
{
  return mpz_to_char(x->QFLevel());
}


void Matrix_mpz_symmetric_swap(Matrix_mpz* x, int i, int j)
{
  x->SwapSymmetric(i, j);
}


void Matrix_mpz_symmetric_multiply(Matrix_mpz* x, int i, char* y)
{
  mpz_class c(*y);
  x->MultiplySymmetric(c, i);
}


void Matrix_mpz_symmetric_divide(Matrix_mpz* x, int i, char* y)
{
  mpz_class c(*y);
  x->DivideSymmetric(c, i);
}


void Matrix_mpz_symmetric_add(Matrix_mpz* x, int i, int j, char* y)
{
  mpz_class c(*y);
  x->AddSymmetric(c, i, j);
}


Matrix_mpz* Matrix_mpz_local_normal_form(Matrix_mpz* x, char* p)
{
  mpz_class _p(*p);
  Matrix_mpz* t = new Matrix_mpz();
  *t = x->GetLocalNormal(_p);
  return t;
}


Matrix_mpz* Matrix_mpz_local_diagonal_form(Matrix_mpz* x, char* p)
{
  mpz_class _p(*p);
  Matrix_mpz* t = new Matrix_mpz();
  *t = x->LocalDiagonal(_p);
  return t;
}


long Matrix_mpz_hasse_invariant(Matrix_mpz* x, char* p)
{
  mpz_class _p(*p);
  return x->HasseInvariant(_p);
}


int Matrix_mpz_is_anisotropic(Matrix_mpz* x, char* p)
{
  mpz_class _p(*p);
  return x->IsAnisotropic(_p);
}


int Matrix_mpz_is_isotropic(Matrix_mpz* x, char* p)
{
  mpz_class _p(*p);
  return x->IsIsotropic(_p);
}

int Matrix_mpz_is_quadratic_form(Matrix_mpz* x)
{
  return x->IsQuadraticForm();
}

int Matrix_mpz_is_symmetric(Matrix_mpz* x)
{
  return x->IsSymmetric();
}

char* Matrix_mpz_anisotropic_primes(Matrix_mpz* x)
{
  valarray<mpz_class> v = x->AnisotropicPrimes();
  ostringstream instore;
  for (unsigned int i=0; i<v.size(); i++) {
    instore << v[i];
    if (i+1 < v.size())
      instore << ", ";
  }
  return stringstream_to_char(instore);
}


char* Matrix_mpz_local_constants(Matrix_mpz* x, char* p, char* T)
{
  mpz_class _p(*p), _T(*T);
  ostringstream instore;
  instore << x->LocalConstantCp(_p, _T);
  return stringstream_to_char(instore);
}


int Matrix_mpz_is_stable(Matrix_mpz* x, char* p, char* T)
{
  mpz_class _p(*p), _T(*T);
  return x->IsStable(_p, _T);
}


int Matrix_mpz_cmp(struct Matrix_mpz* x, struct Matrix_mpz* y)
{
  return *x == *y;
}


