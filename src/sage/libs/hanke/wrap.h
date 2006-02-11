#ifdef __cplusplus
#include "Matrix_mpz/Matrix_mpz_header.h"
#endif

#ifdef __cplusplus
#define EXTERN extern "C"
#else
#define EXTERN
#endif

EXTERN struct Matrix_mpz;

EXTERN struct Matrix_mpz* Matrix_mpz_new(int r, int s);

EXTERN void Matrix_mpz_del(struct Matrix_mpz* m);

EXTERN char* Matrix_mpz_repr(struct Matrix_mpz* m);

EXTERN int Matrix_mpz_nrows(struct Matrix_mpz* m);
EXTERN int Matrix_mpz_ncols(struct Matrix_mpz* m);

EXTERN void Matrix_mpz_setitem(struct Matrix_mpz* x, int i, int j, char* z);

EXTERN char* Matrix_mpz_getitem(struct Matrix_mpz* x, int i, int j);

EXTERN char* Matrix_mpz_determinant(const struct Matrix_mpz* x);

EXTERN struct Matrix_mpz* Matrix_mpz_adjoint(const struct Matrix_mpz* x);

EXTERN char* Matrix_mpz_Local_Density(struct Matrix_mpz* x, char* p, char* m);

EXTERN char* Matrix_mpz_Local_Primitive_Density(struct Matrix_mpz* x, char* p, char* m);

EXTERN char* Matrix_mpz_level(struct Matrix_mpz* x);

EXTERN void Matrix_mpz_symmetric_swap(struct Matrix_mpz* x, int i, int j);

EXTERN void Matrix_mpz_symmetric_multiply(struct Matrix_mpz* x, int i, char* y);

EXTERN void Matrix_mpz_symmetric_divide(struct Matrix_mpz* x, int i, char* y);

EXTERN void Matrix_mpz_symmetric_add(struct Matrix_mpz* x, int i, int j, char* y);

EXTERN struct Matrix_mpz* Matrix_mpz_local_normal_form(struct Matrix_mpz* x, char* p);

EXTERN struct Matrix_mpz* Matrix_mpz_local_diagonal_form(struct Matrix_mpz* x, char* p);

EXTERN long Matrix_mpz_hasse_invariant(struct Matrix_mpz* x, char* p);

EXTERN int Matrix_mpz_is_anisotropic(struct Matrix_mpz* x, char* p);

EXTERN int Matrix_mpz_is_isotropic(struct Matrix_mpz* x, char* p);

EXTERN int Matrix_mpz_is_quadratic_form(struct Matrix_mpz* x);

EXTERN int Matrix_mpz_is_symmetric(struct Matrix_mpz* x);

EXTERN char* Matrix_mpz_anisotropic_primes(struct Matrix_mpz* x);

EXTERN char* Matrix_mpz_local_constants(struct Matrix_mpz* x, char* p, char* T);

EXTERN int Matrix_mpz_is_stable(struct Matrix_mpz* x, char* p, char* T);

EXTERN int Matrix_mpz_cmp(struct Matrix_mpz* x, struct Matrix_mpz* y);



