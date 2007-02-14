#ifdef __cplusplus
#define EXTERN extern "C"
#else
#define EXTERN
#endif

/* linbox_minpoly allocates space for minpoly, so you have to call linbox_delete_array
   to free it up afterwards. */
EXTERN void linbox_integer_dense_minpoly_hacked(mpz_t** minpoly, size_t* degree,
                  size_t n, mpz_t** matrix, int do_minpoly);
EXTERN void linbox_integer_dense_minpoly(mpz_t** minpoly, size_t* degree,
                  size_t n, mpz_t** matrix);
EXTERN void linbox_integer_dense_charpoly(mpz_t** charpoly, size_t* degree,
                  size_t n, mpz_t** matrix);
EXTERN void linbox_integer_dense_delete_array(mpz_t* f);

/* ans must be a pre-allocated and pre-initialized array of GMP ints. */
EXTERN int linbox_integer_dense_matrix_matrix_multiply(mpz_t** ans, mpz_t **A, mpz_t **B,
			      size_t A_nr, size_t A_nc, size_t B_nr, size_t B_nc);

EXTERN unsigned long linbox_integer_dense_rank(mpz_t** matrix, size_t nrows,
					       size_t ncols);
