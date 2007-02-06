#ifdef __cplusplus
#define EXTERN extern "C"
#else
#define EXTERN
#endif


EXTERN void linbox_rational_dense_echelon_form(mpq_t** matrix, size_t nr, size_t nc);

