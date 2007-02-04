#ifdef __cplusplus

#endif

#ifdef __cplusplus
#define EXTERN extern "C"
#else
#define EXTERN
#endif

EXTERN void linbox_minpoly(mpz_t** minpoly, size_t* degree, size_t n, mpz_t** matrix);
EXTERN void linbox_delete_array(mpz_t* f);
