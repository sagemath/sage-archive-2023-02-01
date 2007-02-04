#ifdef __cplusplus

#endif

#ifdef __cplusplus
#define EXTERN extern "C"
#else
#define EXTERN
#endif

/* linbox_minpoly allocates space for minpoly, so you have to call linbox_delete_array
   to free it up afterwards. */
EXTERN void linbox_minpoly(mpz_t** minpoly, size_t* degree, size_t n, mpz_t** matrix);
EXTERN void linbox_delete_array(mpz_t* f);

/* ans must be a pre-allocated and pre-initialized array of GMP ints. */
EXTERN int linbox_matrix_matrix_multiply(mpz_t** ans, mpz_t **A, mpz_t **B,
					  size_t A_nr, size_t A_nc, size_t B_nr, size_t B_nc);
