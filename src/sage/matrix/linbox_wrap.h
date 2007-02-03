#ifdef __cplusplus

#endif

#ifdef __cplusplus
#define EXTERN extern "C"
#else
#define EXTERN
#endif

EXTERN mpz_t* linbox_minpoly(size_t n, mpz_t** matrix);
