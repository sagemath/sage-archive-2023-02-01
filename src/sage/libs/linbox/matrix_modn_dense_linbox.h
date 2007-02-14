#ifdef __cplusplus
#define EXTERN extern "C"
#else
#define EXTERN
#endif


#include<stddef.h>
typedef size_t mod_int;

EXTERN int linbox_modn_dense_echelonize(mod_int modulus,
					       mod_int** matrix, size_t nrows, size_t ncols);


EXTERN void linbox_modn_dense_minpoly(mod_int modulus, mod_int **mp, size_t* degree,
				      size_t n, mod_int **matrix, int do_minpoly);

EXTERN void linbox_modn_dense_delete_array(mod_int *f);

EXTERN int linbox_modn_dense_matrix_matrix_multiply(mod_int modulus, mod_int **ans, mod_int **A, mod_int **B,
					     size_t A_nr, size_t A_nc, size_t B_nr, size_t B_nc);

EXTERN int linbox_modn_dense_rank(mod_int modulus,
				  mod_int** matrix, size_t nrows, size_t ncols);
