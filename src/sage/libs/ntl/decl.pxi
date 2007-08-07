cdef extern from "Python.h":
    object PyString_FromString(char *v)
    char* PyString_AsString(object string)

cdef extern from "stdlib.h":
    void free(void *ptr)

cdef extern from "gmp.h":
    ctypedef struct mpz_t:
        pass

cdef extern from "ntl_wrap.h":

    # really, this is from NTL/ZZ.h
    ctypedef struct ZZ_c "struct ZZ":
        pass

    # Some boiler-plate
    ZZ_c* ZZ_new "New<ZZ>"()
    ZZ_c* ZZ_construct "Construct<ZZ>"(void *mem)
    void ZZ_destruct "Destruct<ZZ>"(ZZ_c *mem)
    void ZZ_delete "Delete<ZZ>"(ZZ_c *mem)
    void ZZ_from_str "_from_str<ZZ>"(ZZ_c* dest, char* s)

    void conv_ZZ_int "conv"(ZZ_c x, int i)
    void add_ZZ "add"( ZZ_c x, ZZ_c a, ZZ_c b)
    void sub_ZZ "sub"( ZZ_c x, ZZ_c a, ZZ_c b)
    void mul_ZZ "mul"( ZZ_c x, ZZ_c a, ZZ_c b)
    void mul_ZZ_long "mul"( ZZ_c x, ZZ_c a, long b)
    void div_ZZ_ZZ "div"( ZZ_c x, ZZ_c a, ZZ_c b)
    void GCD_ZZ "GCD"(ZZ_c d, ZZ_c a, ZZ_c b)
    void negate(ZZ_c x, ZZ_c a)
    long sign(ZZ_c d)
#    ZZ_c *ZZ_factory "new NTL::ZZ" ()

    #### ZZ_c
    ZZ_c* new_ZZ()
    ZZ_c* str_to_ZZ(char* s)
    void del_ZZ(ZZ_c* n)
    char* ZZ_to_str(ZZ_c* x)
    int ZZ_equal(ZZ_c* x, ZZ_c* y)
    ZZ_c* ZZ_add(ZZ_c* x, ZZ_c* y)
    ZZ_c* ZZ_sub(ZZ_c* x, ZZ_c* y)
    ZZ_c* ZZ_mul(ZZ_c* x, ZZ_c* y)
    ZZ_c* ZZ_pow(ZZ_c* x, long e)
    int ZZ_is_one(ZZ_c* x)
    int ZZ_is_zero(ZZ_c* x)
    ZZ_c* ZZ_neg(ZZ_c* x)
    ZZ_c* ZZ_copy(ZZ_c* x)
    # Random-number generation
    void setSeed(ZZ_c* x)
    ZZ_c* ZZ_randomBnd(ZZ_c* x)
    ZZ_c* ZZ_randomBits(long n)
    void ZZ_to_mpz(mpz_t* output, ZZ_c* x)
    void mpz_to_ZZ(ZZ_c *output, mpz_t* x)
    cdef int ZZ_to_int(ZZ_c* x)
    cdef ZZ_c* int_to_ZZ(int value)
    cdef void ZZ_set_from_int(ZZ_c* x, int value)

    #### ZZ_p_c
    ctypedef struct ZZ_p_c "struct ZZ_p":
        pass

    # Some boiler-plate
    ZZ_p_c* ZZ_p_new "New<ZZ_p_c>"()
    ZZ_p_c* ZZ_p_construct "Construct<ZZ_p_c>"(void *mem)
    void ZZ_p_destruct "Destruct<ZZ_p_c>"(ZZ_p_c *mem)
    void ZZ_p_delete "Delete<ZZ_p_c>"(ZZ_p_c *mem)
    void ZZ_p_from_str "_from_str<ZZ_p_c>"(ZZ_p_c* dest, char* s)

    ZZ_p_c* new_ZZ_p()
    ZZ_p_c* str_to_ZZ_p(char* s)
    void del_ZZ_p(ZZ_p_c* n)
    char* ZZ_p_to_str(ZZ_p_c* x)
    ZZ_p_c* ZZ_p_add(ZZ_p_c* x, ZZ_p_c* y)
    ZZ_p_c* ZZ_p_sub(ZZ_p_c* x, ZZ_p_c* y)
    ZZ_p_c* ZZ_p_mul(ZZ_p_c* x, ZZ_p_c* y)
    ZZ_p_c* ZZ_p_pow(ZZ_p_c* x, long e)
    int ZZ_p_is_one(ZZ_p_c* x)
    int ZZ_p_is_zero(ZZ_p_c* x)
    ZZ_p_c* ZZ_p_neg(ZZ_p_c* x)
    void ntl_ZZ_set_modulus(ZZ_c* x)
    int ZZ_p_eq( ZZ_p_c* x,  ZZ_p_c* y)
    ZZ_p_c* ZZ_p_inv(ZZ_p_c* x)
    ZZ_p_c* ZZ_p_random()

    #### ZZX_c

    # really, this is from NTL/ZZX.h
    ctypedef struct ZZX_c "struct ZZX":
        pass

    # Some boiler-plate
    ZZX_c* ZZX_new "New<ZZX>"()
    ZZX_c* ZZX_construct "Construct<ZZX>"(void *mem)
    void ZZX_destruct "Destruct<ZZX>"(ZZX_c *mem)
    void ZZX_delete "Delete<ZZX>"(ZZX_c *mem)
    void ZZX_from_str "_from_str<ZZX>"(ZZX_c* dest, char* s)

    void GetCoeff(ZZ_c x, ZZX_c a, long i)
    void SetCoeff(ZZX_c x, long i, ZZ_c a)
    long IsZero_ZZX "IsZero"(ZZX_c a)
    # f must be monic!
    void MulMod_ZZX "MulMod"(ZZX_c x, ZZX_c a, ZZX_c b, ZZX_c f)
    void mul_ZZX_long "mul"( ZZX_c x, ZZX_c a, long b)
    void mul_ZZX_ZZ "mul"( ZZX_c x, ZZX_c a, ZZ_c b)
    void mul_ZZX "mul"( ZZX_c x, ZZX_c a, ZZX_c b)
    void add_ZZX "add"( ZZX_c x, ZZX_c a, ZZX_c b)
    void sub_ZZX "sub"( ZZX_c x, ZZX_c a, ZZX_c b)
    void div_ZZX_ZZ "div"( ZZX_c x, ZZX_c a, ZZ_c b)
    long deg( ZZX_c x )
    void rem_ZZX "rem"(ZZX_c r, ZZX_c a, ZZX_c b)
    void XGCD_ZZX "XGCD"(ZZ_c r, ZZX_c s, ZZX_c t, ZZX_c a, ZZX_c b, long deterministic)
    void content(ZZ_c d, ZZX_c f)
#    ZZX_c *ZZX_factory "new NTL::ZZX" ()

    void ZZX_square_free_decomposition(ZZX_c*** v, long** e, long* n, ZZX_c* x)

    ZZX_c* ZZX_init()
    ZZX_c* str_to_ZZX(char* s)
    char* ZZX_repr(ZZX_c* x)
    void ZZX_dealloc(ZZX_c* x)
    ZZX_c* ZZX_copy(ZZX_c* x)
    void ZZX_setitem(ZZX_c* x, long i, char* a)
    char* ZZX_getitem(ZZX_c* x, long i)
    ZZX_c* ZZX_add(ZZX_c* x, ZZX_c* y)
    ZZX_c* ZZX_sub(ZZX_c* x, ZZX_c* y)
    ZZX_c* ZZX_mul(ZZX_c* x, ZZX_c* y)
    ZZX_c* ZZX_div(ZZX_c* x, ZZX_c* y, int* divisible)
    ZZX_c* ZZX_mod(ZZX_c* x, ZZX_c* y)
    void ZZX_quo_rem(ZZX_c* x, ZZX_c* other, ZZX_c** r, ZZX_c** q)
    ZZX_c* ZZX_square(ZZX_c* x)
    int ZZX_equal(ZZX_c* x, ZZX_c* y)
    int ZZX_is_zero(ZZX_c* x)
    int ZZX_is_one(ZZX_c* x)
    int ZZX_is_monic(ZZX_c* x)
    ZZX_c* ZZX_neg(ZZX_c* x)
    ZZX_c* ZZX_left_shift(ZZX_c* x, long n)
    ZZX_c* ZZX_right_shift(ZZX_c* x, long n)
    char* ZZX_content(ZZX_c* x)
    ZZX_c* ZZX_primitive_part(ZZX_c* x)
    void ZZX_pseudo_quo_rem(ZZX_c* x, ZZX_c* y, ZZX_c** r, ZZX_c** q)
    ZZX_c* ZZX_gcd(ZZX_c* x, ZZX_c* y)
    ZZX_c* ZZX_xgcd(ZZX_c* x, ZZX_c* y, ZZ_c** r, ZZX_c** s, ZZX_c** t, int proof)
    long ZZX_degree(ZZX_c* x)
    ZZ_c* ZZX_leading_coefficient(ZZX_c* x)
    char* ZZX_constant_term(ZZX_c* x)
    void ZZX_set_x(ZZX_c* x)
    int ZZX_is_x(ZZX_c* x)
    ZZX_c* ZZX_derivative(ZZX_c* x)
    ZZX_c* ZZX_reverse(ZZX_c* x)
    ZZX_c* ZZX_reverse_hi(ZZX_c* x, long hi)
    ZZX_c* ZZX_truncate(ZZX_c* x, long m)
    ZZX_c* ZZX_multiply_and_truncate(ZZX_c* x, ZZX_c* y, long m)
    ZZX_c* ZZX_square_and_truncate(ZZX_c* x, long m)
    ZZX_c* ZZX_invert_and_truncate(ZZX_c* x, long m)
    ZZX_c* ZZX_multiply_mod(ZZX_c* x, ZZX_c* y,  ZZX_c* modulus)
    ZZ_c* ZZX_trace_mod(ZZX_c* x, ZZX_c* y)
    char* ZZX_trace_list(ZZX_c* x)
    ZZ_c* ZZX_resultant(ZZX_c* x, ZZX_c* y, int proof)
    ZZ_c* ZZX_norm_mod(ZZX_c* x, ZZX_c* y, int proof)
    ZZ_c* ZZX_discriminant(ZZX_c* x, int proof)
    ZZ_c* ZZX_polyeval(ZZX_c* x, ZZ_c* a)
    ZZX_c* ZZX_charpoly_mod(ZZX_c* x, ZZX_c* y, int proof)
    ZZX_c* ZZX_minpoly_mod(ZZX_c* x, ZZX_c* y)
    void ZZX_clear(ZZX_c* x)
    void ZZX_preallocate_space(ZZX_c* x, long n)
    void ZZX_getitem_as_mpz(mpz_t* output, ZZX_c* x, long i)
    cdef void ZZX_setitem_from_int(ZZX_c* x, long i, int value)
    cdef int ZZX_getitem_as_int(ZZX_c* x, long i)

    #### ZZ_pX_c
    ctypedef struct ZZ_pX_c "struct ZZ_pX":
        pass

    ZZ_pX_c* ZZ_pX_init()
    ZZ_pX_c* str_to_ZZ_pX(char* s)
    char* ZZ_pX_repr(ZZ_pX_c* x)
    void ZZ_pX_dealloc(ZZ_pX_c* x)
    ZZ_pX_c* ZZ_pX_copy(ZZ_pX_c* x)
    void ZZ_pX_setitem(ZZ_pX_c* x, long i, char* a)
    char* ZZ_pX_getitem(ZZ_pX_c* x, long i)
    ZZ_pX_c* ZZ_pX_add(ZZ_pX_c* x, ZZ_pX_c* y)
    ZZ_pX_c* ZZ_pX_sub(ZZ_pX_c* x, ZZ_pX_c* y)
    ZZ_pX_c* ZZ_pX_mul(ZZ_pX_c* x, ZZ_pX_c* y)
    ZZ_pX_c* ZZ_pX_div(ZZ_pX_c* x, ZZ_pX_c* y, int* divisible)
    ZZ_pX_c* ZZ_pX_mod(ZZ_pX_c* x, ZZ_pX_c* y)
    void ZZ_pX_quo_rem(ZZ_pX_c* x, ZZ_pX_c* other, ZZ_pX_c** r, ZZ_pX_c** q)
    ZZ_pX_c* ZZ_pX_square(ZZ_pX_c* x)
    int ZZ_pX_equal(ZZ_pX_c* x, ZZ_pX_c* y)
    int ZZ_pX_is_zero(ZZ_pX_c* x)
    int ZZ_pX_is_one(ZZ_pX_c* x)
    int ZZ_pX_is_monic(ZZ_pX_c* x)
    ZZ_pX_c* ZZ_pX_neg(ZZ_pX_c* x)
    ZZ_pX_c* ZZ_pX_left_shift(ZZ_pX_c* x, long n)
    ZZ_pX_c* ZZ_pX_right_shift(ZZ_pX_c* x, long n)
    ZZ_pX_c* ZZ_pX_gcd(ZZ_pX_c* x, ZZ_pX_c* y)
    ZZ_pX_c* ZZ_pX_xgcd(ZZ_pX_c** d, ZZ_pX_c** s, ZZ_pX_c** t, ZZ_pX_c* a, ZZ_pX_c* b)
    ZZ_pX_c* ZZ_pX_plain_xgcd(ZZ_pX_c** d, ZZ_pX_c** s, ZZ_pX_c** t, ZZ_pX_c* a, ZZ_pX_c* b)
    long ZZ_pX_degree(ZZ_pX_c* x)
    ZZ_p_c* ZZ_pX_leading_coefficient(ZZ_pX_c* x)
    char* ZZ_pX_constant_term(ZZ_pX_c* x)
    void ZZ_pX_set_x(ZZ_pX_c* x)
    int ZZ_pX_is_x(ZZ_pX_c* x)
    ZZ_pX_c* ZZ_pX_derivative(ZZ_pX_c* x)
    ZZ_pX_c* ZZ_pX_reverse(ZZ_pX_c* x)
    ZZ_pX_c* ZZ_pX_reverse_hi(ZZ_pX_c* x, long hi)
    ZZ_pX_c* ZZ_pX_truncate(ZZ_pX_c* x, long m)
    ZZ_pX_c* ZZ_pX_multiply_and_truncate(ZZ_pX_c* x, ZZ_pX_c* y, long m)
    ZZ_pX_c* ZZ_pX_square_and_truncate(ZZ_pX_c* x, long m)
    ZZ_pX_c* ZZ_pX_invert_and_truncate(ZZ_pX_c* x, long m)
    ZZ_pX_c* ZZ_pX_multiply_mod(ZZ_pX_c* x, ZZ_pX_c* y,  ZZ_pX_c* modulus)
    ZZ_p_c* ZZ_pX_trace_mod(ZZ_pX_c* x, ZZ_pX_c* y)
    char* ZZ_pX_trace_list(ZZ_pX_c* x)
    ZZ_p_c* ZZ_pX_resultant(ZZ_pX_c* x, ZZ_pX_c* y)
    ZZ_p_c* ZZ_pX_norm_mod(ZZ_pX_c* x, ZZ_pX_c* y)
    ZZ_p_c* ZZ_pX_discriminant(ZZ_pX_c* x)
    ZZ_pX_c* ZZ_pX_charpoly_mod(ZZ_pX_c* x, ZZ_pX_c* y)
    ZZ_pX_c* ZZ_pX_minpoly_mod(ZZ_pX_c* x, ZZ_pX_c* y)
    void ZZ_pX_clear(ZZ_pX_c* x)
    void ZZ_pX_preallocate_space(ZZ_pX_c* x, long n)

    void ZZ_pX_factor(ZZ_pX_c*** v, long** e, long* n, ZZ_pX_c* x, long verbose)
    void ZZ_pX_linear_roots(ZZ_p_c*** v, long* n, ZZ_pX_c* x)

    #### mat_ZZ_c
    ctypedef struct mat_ZZ_c "struct mat_ZZ":
        pass

    # Some boiler-plate
    mat_ZZ_c* mat_ZZ_construct "Construct<mat_ZZ_c>"(void *mem)
    void mat_ZZ_destruct "Destruct<mat_ZZ_c>"(mat_ZZ_c *mem)
    void mat_ZZ_delete "Delete<mat_ZZ_c>"(mat_ZZ_c *mem)

    mat_ZZ_c* new_mat_ZZ(long nrows, long ncols)
    void del_mat_ZZ(mat_ZZ_c* n)
    char* mat_ZZ_to_str(mat_ZZ_c* x)
    mat_ZZ_c* mat_ZZ_add(mat_ZZ_c* x, mat_ZZ_c* y)
    mat_ZZ_c* mat_ZZ_sub(mat_ZZ_c* x, mat_ZZ_c* y)
    mat_ZZ_c* mat_ZZ_mul(mat_ZZ_c* x, mat_ZZ_c* y)
    mat_ZZ_c* mat_ZZ_pow(mat_ZZ_c* x, long e)
    long mat_ZZ_nrows(mat_ZZ_c* x)
    long mat_ZZ_ncols(mat_ZZ_c* x)
    void mat_ZZ_setitem(mat_ZZ_c* x, int i, int j, ZZ_c* z)
    ZZ_c* mat_ZZ_getitem(mat_ZZ_c* x, int i, int j)
    ZZ_c* mat_ZZ_determinant(mat_ZZ_c* x, long deterministic)
    mat_ZZ_c* mat_ZZ_HNF(mat_ZZ_c* A, ZZ_c* D)
    ZZX_c* mat_ZZ_charpoly(mat_ZZ_c* A)
    long mat_ZZ_LLL(ZZ_c **det, mat_ZZ_c *x, long a, long b, long verbose)
    long mat_ZZ_LLL_U(ZZ_c **det, mat_ZZ_c *x, mat_ZZ_c *U, long a, long b, long verbose)

    #### GF2X_c
    ctypedef struct GF2X_c "struct GF2X":
        pass

    GF2X_c* GF2X_new "New<GF2X_c>"()
    GF2X_c* GF2X_construct "Construct<GF2X_c>"(void *mem)
    void GF2X_destruct "Destruct<GF2X_c>"(GF2X_c *mem)
    void GF2X_delete "Delete<GF2X_c>"(GF2X_c *mem)
    void GF2X_from_str "_from_str<GF2X_c>"(GF2X_c* dest, char* s)

    GF2X_c* new_GF2X()
    GF2X_c* str_to_GF2X(char* s)
    void del_GF2X(GF2X_c* n)
    char* GF2X_to_str(GF2X_c* x)
    GF2X_c* GF2X_add(GF2X_c* x, GF2X_c* y)
    GF2X_c* GF2X_sub(GF2X_c* x, GF2X_c* y)
    GF2X_c* GF2X_mul(GF2X_c* x, GF2X_c* y)
    GF2X_c* GF2X_pow(GF2X_c* x, long e)
    int GF2X_eq( GF2X_c* x,  GF2X_c* y)
    int GF2X_is_one(GF2X_c* x)
    int GF2X_is_zero(GF2X_c* x)
    GF2X_c* GF2X_neg(GF2X_c* x)
    GF2X_c* GF2X_copy(GF2X_c* x)
    long GF2X_deg(GF2X_c* x)
    void GF2X_hex(long h)
    char *GF2X_to_bin(GF2X_c* x)
    char *GF2X_to_hex(GF2X_c* x)


    #### GF2E_c
    ctypedef struct GF2E_c "struct GF2E":
        pass

    GF2E_c* GF2E_new "New<GF2E_c>"()
    GF2E_c* GF2E_construct "Construct<GF2E_c>"(void *mem)
    void GF2E_destruct "Destruct<GF2E_c>"(GF2E_c *mem)
    void GF2E_delete "Delete<GF2E_c>"(GF2E_c *mem)
    void GF2E_from_str "_from_str<GF2E_c>"(GF2E_c* dest, char* s)

    void ntl_GF2E_set_modulus(GF2X_c *x)
    GF2E_c* new_GF2E()
    GF2E_c* str_to_GF2E(char *s)
    void del_GF2E(GF2E_c *n)
    char *GF2E_to_str(GF2E_c* x)
    GF2E_c *GF2E_add(GF2E_c *x, GF2E_c *y)
    GF2E_c *GF2E_sub(GF2E_c *x, GF2E_c *y)
    GF2E_c *GF2E_mul(GF2E_c *x, GF2E_c *y)
    GF2E_c *GF2E_pow(GF2E_c *x, long e)
    int GF2E_eq( GF2E_c* x,  GF2E_c* y)
    int GF2E_is_one(GF2E_c *x)
    int GF2E_is_zero(GF2E_c *x)
    GF2E_c *GF2E_neg(GF2E_c *x)
    GF2E_c *GF2E_copy(GF2E_c *x)
    long GF2E_degree()
    GF2X_c *GF2E_modulus()
    GF2E_c *GF2E_random()
    long  GF2E_trace(GF2E_c *x)
    GF2X_c *GF2E_ntl_GF2X(GF2E_c *x)

    #### GF2EX_c
    ctypedef struct GF2EX_c "struct GF2EX":
        pass

    GF2EX_c* GF2EX_new "New<GF2EX_c>"()
    GF2EX_c* GF2EX_construct "Construct<GF2EX_c>"(void *mem)
    void GF2EX_destruct "Destruct<GF2EX_c>"(GF2EX_c *mem)
    void GF2EX_delete "Delete<GF2EX_c>"(GF2EX_c *mem)
    void GF2EX_from_str "_from_str<GF2EX_c>"(GF2EX_c* dest, char* s)

    GF2EX_c* new_GF2EX()
    GF2EX_c* str_to_GF2EX(char* s)
    void del_GF2EX(GF2EX_c* n)
    char* GF2EX_to_str(GF2EX_c* x)
    GF2EX_c* GF2EX_add(GF2EX_c* x, GF2EX_c* y)
    GF2EX_c* GF2EX_sub(GF2EX_c* x, GF2EX_c* y)
    GF2EX_c* GF2EX_mul(GF2EX_c* x, GF2EX_c* y)
    GF2EX_c* GF2EX_pow(GF2EX_c* x, long e)
    int GF2EX_is_one(GF2EX_c* x)
    int GF2EX_is_zero(GF2EX_c* x)
    GF2EX_c* GF2EX_neg(GF2EX_c* x)
    GF2EX_c* GF2EX_copy(GF2EX_c* x)


    #### mat_GF2E_c
    ctypedef struct mat_GF2E_c "struct mat_GF2E":
        pass

    mat_GF2E_c* mat_GF2E_new "New<mat_GF2E_c>"()
    mat_GF2E_c* mat_GF2E_construct "Construct<mat_GF2E_c>"(void *mem)
    void mat_GF2E_destruct "Destruct<mat_GF2E_c>"(mat_GF2E_c *mem)
    void mat_GF2E_delete "Delete<mat_GF2E_c>"(mat_GF2E_c *mem)
    void mat_GF2E_from_str "_from_str<mat_GF2E_c>"(mat_GF2E_c* dest, char* s)

    mat_GF2E_c* new_mat_GF2E(long nrows, long ncols)
    void del_mat_GF2E(mat_GF2E_c* n)
    char* mat_GF2E_to_str(mat_GF2E_c* x)
    mat_GF2E_c* mat_GF2E_add(mat_GF2E_c* x, mat_GF2E_c* y)
    mat_GF2E_c* mat_GF2E_sub(mat_GF2E_c* x, mat_GF2E_c* y)
    mat_GF2E_c* mat_GF2E_mul(mat_GF2E_c* x, mat_GF2E_c* y)
    mat_GF2E_c* mat_GF2E_pow(mat_GF2E_c* x, long e)
    long mat_GF2E_nrows(mat_GF2E_c* x)
    long mat_GF2E_ncols(mat_GF2E_c* x)
    void mat_GF2E_setitem(mat_GF2E_c* x, int i, int j, GF2E_c* z)
    GF2E_c* mat_GF2E_getitem(mat_GF2E_c* x, int i, int j)
    GF2E_c* mat_GF2E_determinant(mat_GF2E_c* x)
    long mat_GF2E_gauss(mat_GF2E_c *x, long w)
    long mat_GF2E_is_zero(mat_GF2E_c *x)
    mat_GF2E_c* mat_GF2E_transpose(mat_GF2E_c *x)

