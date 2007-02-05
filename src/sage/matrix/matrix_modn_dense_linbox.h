#ifdef __cplusplus
#define EXTERN extern "C"
#else
#define EXTERN
#endif

#include <stdlib.h>

EXTERN int linbox_matrix_modn_dense_echelonize(unsigned long modulus,
					       unsigned long** matrix, size_t nrows, size_t ncols);


