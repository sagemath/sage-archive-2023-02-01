include "../../ext/cdefs.pxi"

cdef extern from "stdlib.h":
    void free(void *ptr)

cdef extern from "ntl_wrap.h":

    # really, this is from NTL/ZZ.h
    #### ZZ_c
    ctypedef struct ZZ_c "struct ZZ":
        pass

    ctypedef struct vec_ZZ_c "vec_ZZ":
        ZZ_c RawGet(long i)
        ZZ_c *elts()
        long length()

    void del_charstar(char*)

    # Some boiler-plate
    ZZ_c* ZZ_new "New<ZZ>"()
    ZZ_c* ZZ_construct "Construct<ZZ>"(void *mem)
    void ZZ_destruct "Destruct<ZZ>"(ZZ_c *mem)
    void ZZ_delete "Delete<ZZ>"(ZZ_c *mem)
    void ZZ_from_str "_from_str<ZZ>"(ZZ_c* dest, char* s)
    object ZZ_to_PyString "_to_PyString<ZZ>"(ZZ_c *x)
    int ZZ_equal "_equal<ZZ>"(ZZ_c x, ZZ_c y)

    void ZZ_conv_from_int "conv"(ZZ_c x, int i)
    void ZZ_conv_to_int "conv"(int i, ZZ_c x)
    void ZZ_conv_from_long "conv"(ZZ_c x, long l)
    void ZZ_conv_to_long "conv"(long l, ZZ_c x)
    void ZZ_to_mpz(mpz_t* output, ZZ_c* x)
    void mpz_to_ZZ(ZZ_c *output, mpz_t* x)
    cdef int ZZ_to_int(ZZ_c* x)
    cdef ZZ_c* int_to_ZZ(int value)
    cdef void ZZ_set_from_int(ZZ_c* x, int value)

    long ZZ_sign "sign"(ZZ_c a)
    long ZZ_IsZero "IsZero"(ZZ_c a)
    long ZZ_IsOne "IsOne"(ZZ_c a)
    long ZZ_compare "compare"(ZZ_c a, ZZ_c b)
    void ZZ_add "add"( ZZ_c x, ZZ_c a, ZZ_c b)
    void ZZ_add_long "add"(ZZ_c x, ZZ_c a, long b)
    void ZZ_sub "sub"( ZZ_c x, ZZ_c a, ZZ_c b)
    void ZZ_sub_long "sub"(ZZ_c x, long a, ZZ_c b)
    void ZZ_negate "negate"(ZZ_c x, ZZ_c a)
    void ZZ_abs "abs"(ZZ_c x, ZZ_c a)

    void ZZ_mul "mul"( ZZ_c x, ZZ_c a, ZZ_c b)
    void ZZ_mul_long "mul"( ZZ_c x, ZZ_c a, long b)
    void ZZ_sqr "sqr"(ZZ_c x, ZZ_c a)

    void ZZ_DivRem "DivRem"(ZZ_c q, ZZ_c r, ZZ_c a, ZZ_c b)
    void ZZ_div "div"( ZZ_c x, ZZ_c a, ZZ_c b)
    void ZZ_rem "rem"( ZZ_c r, ZZ_c a, ZZ_c b)
    long ZZ_DivRem_long "DivRem"(ZZ_c q, ZZ_c a, long b)
    long ZZ_rem_long "rem"(ZZ_c a, long b)
    long ZZ_divide "divide"(ZZ_c q, ZZ_c a, ZZ_c b)
    long ZZ_divide_long "divide"(ZZ_c q, ZZ_c a, long b)
    long ZZ_divide_test "divide"(ZZ_c a, ZZ_c b)
    long ZZ_divide_test_long "divide"(ZZ_c a, long b)

    void ZZ_GCD "GCD"(ZZ_c d, ZZ_c a, ZZ_c b)
    void ZZ_XGCD "XGCD"(ZZ_c d, ZZ_c s, ZZ_c t, ZZ_c a, ZZ_c b)

    void ZZ_AddMod "AddMod"(ZZ_c x, ZZ_c a, ZZ_c b, ZZ_c n)
    void ZZ_SubMod "SubMod"(ZZ_c x, ZZ_c a, ZZ_c b, ZZ_c n)
    void ZZ_NegateMod "NegateMod"(ZZ_c x, ZZ_c a, ZZ_c n)
    void ZZ_MulMod "MulMod"(ZZ_c x, ZZ_c a, ZZ_c b, ZZ_c n)
    void ZZ_SqrMod "SqrMod"(ZZ_c x, ZZ_c a, ZZ_c n)
    void ZZ_InvMod "InvMod"(ZZ_c x, ZZ_c a, ZZ_c n)
    long ZZ_InvModStatus "InvModStatus"(ZZ_c x, ZZ_c a, ZZ_c n)
    void ZZ_PowerMod "PowerMod"(ZZ_c x, ZZ_c a, ZZ_c e, ZZ_c n)
    void ZZ_PowerMod_long "PowerMod"(ZZ_c x, ZZ_c a, long e, ZZ_c n)

    void ZZ_LeftShift "LeftShift"(ZZ_c x, ZZ_c a, long n)
    void ZZ_RightShift "RightShift"(ZZ_c x, ZZ_c a, long n)

    long ZZ_MakeOdd "MakeOdd"(ZZ_c x)
    long ZZ_NumTwos "NumTwos"(ZZ_c x)
    long ZZ_IsOdd "IsOdd"(ZZ_c a)
    long ZZ_NumBits "NumBits"(ZZ_c a)
    long ZZ_bit "bit"(ZZ_c a, long k)
    void ZZ_trunc "trunc"(ZZ_c x, ZZ_c a, long k)
    void ZZ_SetBit "SetBit"(ZZ_c x, long p)
    void ZZ_SwitchBit "SwitchBit"(ZZ_c x, long p)
    void ZZ_weight "weight"(ZZ_c a)
    void ZZ_bit_and "bit_and"(ZZ_c x, ZZ_c a, ZZ_c b)
    void ZZ_bit_or "bit_or"(ZZ_c x, ZZ_c a, ZZ_c b)
    void ZZ_bit_xor "bit_xor"(ZZ_c x, ZZ_c a, ZZ_c b)

    void ZZ_SetSeed "SetSeed"(ZZ_c s)
    void ZZ_RandomBnd "RandomBnd"(ZZ_c x, ZZ_c n)
    void ZZ_RandomBits "RandomBits"(ZZ_c x, long l)

    long ZZ_CRT "CRT"(ZZ_c a, ZZ_c p, ZZ_c A, ZZ_c P)

    long ZZ_ReconstructRational "ReconstructRational"(ZZ_c a, ZZ_c b, ZZ_c x, ZZ_c m, ZZ_c a_bound, ZZ_c b_bound)

    void ZZ_GenPrime "GenPrime"(ZZ_c n, long l, long err)
    void ZZ_GenGermainPrime "GenGermainPrime"(ZZ_c n, long l, long err)
    long ZZ_ProbPrime "ProbPrime"(ZZ_c n, long NumTrials)
    void ZZ_RandomPrime "RandomPrime"(ZZ_c n, long l, long NumTrials)
    void ZZ_NextPrime "NextPrime"(ZZ_c n, ZZ_c m, long NumTrials)
    long ZZ_MillerWitness "MillerWitness"(ZZ_c n, ZZ_c w)

    void ZZ_power "power"(ZZ_c t, ZZ_c x, long e)
    void ZZ_power2 "power2"(ZZ_c x, long e)

    void ZZ_SqrRoot "SqrRoot"(ZZ_c x, ZZ_c a)

    long ZZ_Jacobi "Jacobi"(ZZ_c a, ZZ_c n)
    void ZZ_SqrRootMod "SqrRootMod"(ZZ_c x, ZZ_c a, ZZ_c n)

    long ZZ_remove(ZZ_c x, ZZ_c a, ZZ_c p) # a la mpz_remove.  Written in ntl_wrap.cpp.

    # Random-number generation
    #void setSeed(ZZ_c* x)
    #ZZ_c* ZZ_randomBnd(ZZ_c* x)
    #ZZ_c* ZZ_randomBits(long n)


    #### ZZ_pContext_c
    ctypedef struct ZZ_pContext_c "struct ZZ_pContext":
        void (*restore)()

    ZZ_pContext_c* ZZ_pContext_new "New<ZZ_pContext>"()
    ZZ_pContext_c* ZZ_pContext_construct "Construct<ZZ_pContext>"(void *mem)
    ZZ_pContext_c* ZZ_pContext_new_ZZ "ZZ_pContext_new"(ZZ_c* p)
    ZZ_pContext_c* ZZ_pContext_construct_ZZ "ZZ_pContext_construct"(void *mem, ZZ_c* p)
    void ZZ_pContext_destruct "Destruct<ZZ_pContext>"(ZZ_pContext_c *mem)
    void ZZ_pContext_delete "Delete<ZZ_pContext>"(ZZ_pContext_c *mem)

    #### ZZ_p_c
    ctypedef struct ZZ_p_c "struct ZZ_p":
        pass

    # Some boiler-plate
    ZZ_p_c* ZZ_p_new "New<ZZ_p>"()
    ZZ_p_c* ZZ_p_construct "Construct<ZZ_p>"(void *mem)
    void ZZ_p_destruct "Destruct<ZZ_p>"(ZZ_p_c *mem)
    void ZZ_p_delete "Delete<ZZ_p>"(ZZ_p_c *mem)
    void ZZ_p_from_str "_from_str<ZZ_p>"(ZZ_p_c* dest, char* s)
    object ZZ_p_to_PyString "_to_PyString<ZZ_p>"(ZZ_p_c *x)
    int ZZ_p_equal "_equal<ZZ_p>"(ZZ_p_c x, ZZ_p_c y)

    char* ZZ_p_to_str(ZZ_p_c* x)
    void ZZ_p_add "add"( ZZ_p_c x, ZZ_p_c a, ZZ_p_c b)
    void ZZ_p_sub "sub"( ZZ_p_c x, ZZ_p_c a, ZZ_p_c b)
    void ZZ_p_mul "mul"( ZZ_p_c x, ZZ_p_c a, ZZ_p_c b)
    void ZZ_p_mul_long "mul"( ZZ_p_c x, ZZ_p_c a, long b)
    void ZZ_p_div "div"( ZZ_p_c x, ZZ_p_c a, ZZ_p_c b)
    void ZZ_p_negate "negate"(ZZ_p_c x, ZZ_p_c a)
    void ZZ_p_power "power"(ZZ_p_c t, ZZ_p_c x, long e)
    int ZZ_p_IsOne "IsOne"(ZZ_p_c x)
    int ZZ_p_IsZero "IsZero"(ZZ_p_c x)
    ZZ_c ZZ_p_rep "rep"(ZZ_p_c z)
    ZZ_p_c* ZZ_p_neg(ZZ_p_c* x)
    void ntl_ZZ_set_modulus(ZZ_c* x)
    void ZZ_p_inv "inv"(ZZ_p_c r, ZZ_p_c x)
    void ZZ_p_random "random"(ZZ_p_c r)
    ZZ_p_c long_to_ZZ_p "to_ZZ_p"(long i)
    ZZ_p_c ZZ_to_ZZ_p "to_ZZ_p"(ZZ_c i)
    int ZZ_p_to_int(ZZ_p_c x)
    ZZ_p_c int_to_ZZ_p(int i)
    void ZZ_p_modulus(ZZ_c* mod, ZZ_p_c* x)

    #### vec_ZZ_p_c
    ctypedef struct vec_ZZ_p_c "struct vec_ZZ_p":
        pass

    #### ZZX_c

    # really, this is from NTL/ZZX.h
    ctypedef struct ZZX_c "struct ZZX":
        vec_ZZ_c rep

    ctypedef struct pair_ZZX_long_c "pair_ZZX_long":
        ZZX_c a
        long b

    ctypedef struct vec_pair_ZZX_long_c "vec_pair_ZZX_long":
        pair_ZZX_long_c RawGet(long i)
        long length()

    # Some boiler-plate
    ZZX_c* ZZX_new "New<ZZX>"()
    ZZX_c* ZZX_construct "Construct<ZZX>"(void *mem)
    void ZZX_destruct "Destruct<ZZX>"(ZZX_c *mem)
    void ZZX_swap "swap"(ZZX_c x, ZZX_c y)
    void ZZX_delete "Delete<ZZX>"(ZZX_c *mem)
    void ZZX_from_str "_from_str<ZZX>"(ZZX_c* dest, char* s)
    object ZZX_to_PyString "_to_PyString<ZZX>"(ZZX_c *x)
    int ZZX_equal "_equal<ZZX>"(ZZX_c x, ZZX_c y)

    void ZZX_PseudoRem "PseudoRem"(ZZX_c x, ZZX_c a, ZZX_c b)
    ZZ_c ZZX_LeadCoeff "LeadCoeff"(ZZX_c x)
    ZZ_c ZZX_ConstTerm "ConstTerm"(ZZX_c x)
    ZZ_c ZZX_coeff "coeff"(ZZX_c a, long i)
    void ZZX_SetCoeff "SetCoeff"(ZZX_c x, long i, ZZ_c a)
    long IsZero_ZZX "IsZero"(ZZX_c a)
    # f must be monic!
    void ZZX_MulMod "MulMod"(ZZX_c x, ZZX_c a, ZZX_c b, ZZX_c f)
    void ZZX_mul_long "mul"( ZZX_c x, ZZX_c a, long b)
    void ZZX_mul_ZZ "mul"( ZZX_c x, ZZX_c a, ZZ_c b)
    void ZZX_mul "mul"( ZZX_c x, ZZX_c a, ZZX_c b)
    void ZZX_add "add"( ZZX_c x, ZZX_c a, ZZX_c b)
    void ZZX_sub "sub"( ZZX_c x, ZZX_c a, ZZX_c b)
    void ZZX_negate "negate"( ZZX_c x, ZZX_c a)
    void ZZX_div_ZZ "div"( ZZX_c x, ZZX_c a, ZZ_c b)
    long ZZX_deg "deg"( ZZX_c x )
    void ZZX_rem "rem"(ZZX_c r, ZZX_c a, ZZX_c b)
    void ZZX_XGCD "XGCD"(ZZ_c r, ZZX_c s, ZZX_c t, ZZX_c a, ZZX_c b, long deterministic)
    void ZZX_content "content"(ZZ_c d, ZZX_c f)
    void ZZX_factor "factor"(ZZ_c c, vec_pair_ZZX_long_c factors, ZZX_c f, long verbose, long bnd)

    void ZZX_squarefree_decomposition(ZZX_c*** v, long** e, long* n, ZZX_c* x)

    char* ZZX_repr(ZZX_c* x)
    ## for cleaning up after ZZX_repr:
    void cpp_delete_array "delete []"(char *str)
    ZZX_c* ZZX_copy(ZZX_c* x)
    ZZX_c* ZZX_div(ZZX_c* x, ZZX_c* y, int* divisible)
    void ZZX_quo_rem(ZZX_c* x, ZZX_c* other, ZZX_c** r, ZZX_c** q)
    ZZX_c* ZZX_square(ZZX_c* x)
    int ZZX_IsZero "IsZero"(ZZX_c x)
    int ZZX_IsOne "IsOne"(ZZX_c x)
    int ZZX_is_monic(ZZX_c* x)
    ZZX_c* ZZX_neg(ZZX_c* x)
    ZZX_c* ZZX_left_shift(ZZX_c* x, long n)
    ZZX_c* ZZX_right_shift(ZZX_c* x, long n)
    ZZX_c* ZZX_primitive_part(ZZX_c* x)
    void ZZX_pseudo_quo_rem(ZZX_c* x, ZZX_c* y, ZZX_c** r, ZZX_c** q)
    ZZX_c* ZZX_gcd(ZZX_c* x, ZZX_c* y)
    ZZX_c* ZZX_xgcd(ZZX_c* x, ZZX_c* y, ZZ_c** r, ZZX_c** s, ZZX_c** t, int proof)
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
        void *rep
        void (* SetMaxLength)(long n)

    ZZ_pX_c* ZZ_pX_new "New<ZZ_pX>"()
    ZZ_pX_c* ZZ_pX_construct "Construct<ZZ_pX>"(void *mem)
    void ZZ_pX_destruct "Destruct<ZZ_pX>"(ZZ_pX_c *mem)
    void ZZ_pX_delete "Delete<ZZ_pX>"(ZZ_pX_c *mem)
    void ZZ_pX_from_str "_from_str<ZZ_pX>"(ZZ_pX_c* dest, char* s)
    object ZZ_pX_to_PyString "_to_PyString<ZZ_pX>"(ZZ_pX_c *x)
    #int ZZ_pX_equal "_equal<ZZ_pX>"(ZZ_pX_c x, ZZ_pX_c y)

    long ZZ_pX_equal "operator=="(ZZ_pX_c a, ZZ_pX_c b)
    long ZZ_pX_IsZero "IsZero"(ZZ_pX_c a)
    long ZZ_pX_IsOne "IsOne"(ZZ_pX_c a)

    void ZZ_pX_add "add"( ZZ_pX_c x, ZZ_pX_c a, ZZ_pX_c b)
    void ZZ_pX_add_long "add"(ZZ_pX_c x, ZZ_pX_c a, long b)
    void ZZ_pX_sub "sub"( ZZ_pX_c x, ZZ_pX_c a, ZZ_pX_c b)
    void ZZ_pX_sub_long "sub"(ZZ_pX_c x, long a, ZZ_pX_c b)
    void ZZ_pX_negate "negate"(ZZ_pX_c x, ZZ_pX_c a)

    void ZZ_pX_mul "mul"( ZZ_pX_c x, ZZ_pX_c a, ZZ_pX_c b)
    void ZZ_pX_mul_long "mul"( ZZ_pX_c x, ZZ_pX_c a, long b)
    void ZZ_pX_mul_ZZ_p "mul"( ZZ_pX_c x, ZZ_pX_c a, ZZ_p_c b)
    void ZZ_pX_rmul "mul"( ZZ_pX_c x, ZZ_pX_c a, ZZ_p_c b)
    void ZZ_pX_sqr "sqr"( ZZ_pX_c x, ZZ_pX_c a)
    long ZZ_pX_power "power"( ZZ_pX_c x, ZZ_pX_c a, long e)

    void ZZ_pX_LeftShift "LeftShift"(ZZ_pX_c x, ZZ_pX_c a, long n)
    void ZZ_pX_RightShift "RightShift"(ZZ_pX_c x, ZZ_pX_c a, long n)

    void ZZ_pX_DivRem "DivRem"(ZZ_pX_c q, ZZ_pX_c r, ZZ_pX_c a, ZZ_pX_c b)
    void ZZ_pX_div_ZZ "div"( ZZ_pX_c q, ZZ_pX_c a, ZZ_pX_c b)
    void ZZ_pX_div "div"( ZZ_pX_c q, ZZ_pX_c a, ZZ_pX_c b)
    long ZZ_pX_divide "divide"( ZZ_pX_c q, ZZ_pX_c a, ZZ_pX_c b)
    void ZZ_pX_rem "rem"(ZZ_pX_c r, ZZ_pX_c a, ZZ_pX_c b)

    void ZZ_pX_GCD "GCD"(ZZ_pX_c x, ZZ_pX_c a, ZZ_pX_c b)
    void ZZ_pX_PlainXGCD "PlainXGCD"(ZZ_pX_c d, ZZ_pX_c s, ZZ_pX_c t, ZZ_pX_c a, ZZ_pX_c b)
    void ZZ_pX_XGCD "XGCD"(ZZ_pX_c d, ZZ_pX_c s, ZZ_pX_c t, ZZ_pX_c a, ZZ_pX_c b)

    long ZZ_pX_deg "deg"( ZZ_pX_c x )
    ZZ_p_c ZZ_pX_coeff "coeff"(ZZ_pX_c a, long i)
    ZZ_p_c ZZ_pX_LeadCoeff "LeadCoeff"(ZZ_pX_c x)
    ZZ_p_c ZZ_pX_ConstTerm "ConstTerm"(ZZ_pX_c x)
    void ZZ_pX_SetCoeff "SetCoeff"(ZZ_pX_c x, long i, ZZ_p_c a)
    void ZZ_pX_SetCoeff_long "SetCoeff"(ZZ_pX_c x, long i, long a)
    void ZZ_pX_SetX "SetX"(ZZ_pX_c x)
    long ZZ_pX_IsX "IsX"(ZZ_pX_c a)
    void ZZ_pX_diff "diff"(ZZ_pX_c x, ZZ_pX_c a)
    void ZZ_pX_MakeMonic "MakeMonic"(ZZ_pX_c x)
    void ZZ_pX_reverse_hi "reverse"(ZZ_pX_c x, ZZ_pX_c a, long hi)
    void ZZ_pX_reverse "reverse"(ZZ_pX_c x, ZZ_pX_c a)
    void ZZ_pX_VectorCopy "VectorCopy"(vec_ZZ_p_c x, ZZ_pX_c a, long n)

    void ZZ_pX_random "random"(ZZ_pX_c x, long n)

    void ZZ_pX_BuildFromRoots "BuildFromRoots"(ZZ_pX_c x, vec_ZZ_p_c a)
    void ZZ_pX_eval "eval"(ZZ_p_c fa, ZZ_pX_c f, ZZ_p_c a)
    void ZZ_pX_eval_vec "eval"(vec_ZZ_p_c fa, ZZ_pX_c f, vec_ZZ_p_c a)
    void ZZ_pX_interpolate "interpolate"(ZZ_pX_c f, vec_ZZ_p_c a, vec_ZZ_p_c b)

    void ZZ_pX_trunc "trunc"(ZZ_pX_c x, ZZ_pX_c a, long n)
    void ZZ_pX_MulTrunc "MulTrunc"(ZZ_pX_c x, ZZ_pX_c a, ZZ_pX_c b, long n)
    void ZZ_pX_SqrTrunc "SqrTrunc"(ZZ_pX_c x, ZZ_pX_c a, long n)
    void ZZ_pX_InvTrunc "InvTrunc"(ZZ_pX_c x, ZZ_pX_c a, long n)

    void ZZ_pX_MulMod "MulMod"(ZZ_pX_c x, ZZ_pX_c a, ZZ_pX_c b, ZZ_pX_c f)
    void ZZ_pX_SqrMod "SqrMod"(ZZ_pX_c x, ZZ_pX_c a, ZZ_pX_c f)
    void ZZ_pX_MulByXMod "MulByXMod"(ZZ_pX_c x, ZZ_pX_c a, ZZ_pX_c f)
    void ZZ_pX_InvMod "InvMod"(ZZ_pX_c x, ZZ_pX_c a, ZZ_pX_c f)
    long ZZ_pX_InvModStatus "InvModStatus"(ZZ_pX_c x, ZZ_pX_c a, ZZ_pX_c f)

    ctypedef struct ZZ_pX_Modulus_c "struct ZZ_pXModulus":
        ZZ_pX_c (* val) ( )
    ZZ_pX_Modulus_c* ZZ_pX_Modulus_new "New<ZZ_pXModulus>"()
    ZZ_pX_Modulus_c* ZZ_pX_Modulus_construct "Construct<ZZ_pXModulus>"(void *mem)
    void ZZ_pX_Modulus_destruct "Destruct<ZZ_pXModulus>"(ZZ_pX_Modulus_c *mem)
    void ZZ_pX_Modulus_delete "Delete<ZZ_pXModulus>"(ZZ_pX_Modulus_c *mem)
    void ZZ_pX_Modulus_from_str "_from_str<ZZ_pXModulus>"(ZZ_pX_Modulus_c* dest, char* s)
    void ZZ_pX_Modulus_build "build"(ZZ_pX_Modulus_c F, ZZ_pX_c f) # MUST be called before using the modulus
    long ZZ_pX_Modulus_deg "deg"(ZZ_pX_Modulus_c F)

    void ZZ_pX_MulMod_pre "MulMod"(ZZ_pX_c x, ZZ_pX_c a, ZZ_pX_c b, ZZ_pX_Modulus_c F)
    void ZZ_pX_SqrMod_pre "SqrMod"(ZZ_pX_c x, ZZ_pX_c a, ZZ_pX_Modulus_c F)
    void ZZ_pX_PowerMod_pre "PowerMod"(ZZ_pX_c x, ZZ_pX_c a, ZZ_c e, ZZ_pX_Modulus_c F)
    void ZZ_pX_PowerMod_long_pre "PowerMod"(ZZ_pX_c x, ZZ_pX_c a, long e, ZZ_pX_Modulus_c F)
    void ZZ_pX_PowerXMod_pre "PowerXMod"(ZZ_pX_c x, ZZ_c e, ZZ_pX_Modulus_c F)
    void ZZ_pX_PowerXMod_long_pre "PowerXMod"(ZZ_pX_c x, long e, ZZ_pX_Modulus_c F)
    void ZZ_pX_PowerXPlusAMod_pre "PowerXPlusAMod"(ZZ_pX_c x, ZZ_p_c a, ZZ_c e, ZZ_pX_Modulus_c F)
    void ZZ_pX_PowerXPlusAMod_long_pre "PowerXPlusAMod"(ZZ_pX_c x, ZZ_p_c a, long e, ZZ_pX_Modulus_c F)
    void ZZ_pX_rem_pre "rem"(ZZ_pX_c x, ZZ_pX_c a, ZZ_pX_Modulus_c F)
    void ZZ_pX_DivRem_pre "DivRem"(ZZ_pX_c q, ZZ_pX_c r, ZZ_pX_c a, ZZ_pX_Modulus_c F)
    void ZZ_pX_div_pre "div"(ZZ_pX_c q, ZZ_pX_c a, ZZ_pX_Modulus_c F)
    void ZZ_pX_InvMod_pre "InvMod"(ZZ_pX_c x, ZZ_pX_c a, ZZ_pX_Modulus_c F)

    ctypedef struct ZZ_pX_Multiplier_c "struct ZZ_pXMultiplier":
        ZZ_pX_c (* val) ( )
    ZZ_pX_Multiplier_c* ZZ_pX_Multiplier_new "New<ZZ_pXMultiplier>"()
    ZZ_pX_Multiplier_c* ZZ_pX_Multiplier_construct "Construct<ZZ_pXMultiplier>"(void *mem)
    void ZZ_pX_Multiplier_destruct "Destruct<ZZ_pXMultiplier>"(ZZ_pX_Multiplier_c *mem)
    void ZZ_pX_Multiplier_delete "Delete<ZZ_pXMultiplier>"(ZZ_pX_Multiplier_c *mem)
    void ZZ_pX_Multiplier_from_str "_from_str<ZZ_pXMultiplier>"(ZZ_pX_Multiplier_c* dest, char* s)
    void ZZ_pX_Multiplier_build "build"(ZZ_pX_Multiplier_c F, ZZ_pX_c b, ZZ_pX_Modulus_c F) # MUST be called before using the multiplier
    void ZZ_pX_MulMod_premul "MulMod"(ZZ_pX_c x, ZZ_pX_c a, ZZ_pX_Multiplier_c B, ZZ_pX_Modulus_c F)

    void ZZ_pX_CompMod "CompMod"(ZZ_pX_c x, ZZ_pX_c g, ZZ_pX_c h, ZZ_pX_Modulus_c F)
    # Add other composition functions here
    # Add power projection routines
    # Add minimum polynomials of recurrence sequences
    void ZZ_pX_MinPolyMod "MinPolyMod"(ZZ_pX_c h, ZZ_pX_c g, ZZ_pX_c f)
    void ZZ_pX_MinPolyMod_pre "MinPolyMod"(ZZ_pX_c h, ZZ_pX_c g, ZZ_pX_Modulus_c F)

    void ZZ_pX_TraceMod "TraceMod"(ZZ_p_c x, ZZ_pX_c a, ZZ_pX_c f)
    void ZZ_pX_TraceMod_pre "TraceMod"(ZZ_p_c x, ZZ_pX_c a, ZZ_pX_Modulus_c F)
    void ZZ_pX_TraceVec "TraceVec"(vec_ZZ_p_c S, ZZ_pX_c f)
    void ZZ_pX_NormMod "NormMod"(ZZ_p_c x, ZZ_pX_c a, ZZ_pX_c f)
    void ZZ_pX_resultant "resultant"(ZZ_p_c x, ZZ_pX_c a, ZZ_pX_c b)
    void ZZ_pX_CharPolyMod "CharPolyMod"(ZZ_pX_c g, ZZ_pX_c a, ZZ_pX_c f)

    void ZZ_pX_clear "clear"(ZZ_pX_c x)
    void ZZ_pX_set "set"(ZZ_pX_c x)

    void ZZ_pX_to_ZZX "conv"(ZZX_c x, ZZ_pX_c a)
    void ZZX_to_ZZ_pX "conv"(ZZ_pX_c x, ZZX_c a)

    #char* ZZ_pX_repr(ZZ_pX_c* x)
    #ZZ_pX_c* ZZ_pX_copy(ZZ_pX_c* x)
    #ZZ_pX_c* ZZ_pX_mod(ZZ_pX_c* x, ZZ_pX_c* y)
    #void ZZ_pX_quo_rem(ZZ_pX_c* x, ZZ_pX_c* other, ZZ_pX_c** r, ZZ_pX_c** q)
    #ZZ_pX_c* ZZ_pX_square(ZZ_pX_c* x)
    #int ZZ_pX_is_monic(ZZ_pX_c* x)
    #ZZ_pX_c* ZZ_pX_neg(ZZ_pX_c* x)
    #ZZ_pX_c* ZZ_pX_left_shift(ZZ_pX_c* x, long n)
    #ZZ_pX_c* ZZ_pX_right_shift(ZZ_pX_c* x, long n)
    #ZZ_pX_c* ZZ_pX_gcd(ZZ_pX_c* x, ZZ_pX_c* y)
    #ZZ_pX_c* ZZ_pX_xgcd(ZZ_pX_c** d, ZZ_pX_c** s, ZZ_pX_c** t, ZZ_pX_c* a, ZZ_pX_c* b)
    #ZZ_pX_c* ZZ_pX_plain_xgcd(ZZ_pX_c** d, ZZ_pX_c** s, ZZ_pX_c** t, ZZ_pX_c* a, ZZ_pX_c* b)
    #void ZZ_pX_set_x(ZZ_pX_c* x)
    #int ZZ_pX_is_x(ZZ_pX_c* x)
    #ZZ_pX_c* ZZ_pX_derivative(ZZ_pX_c* x)
    #ZZ_pX_c* ZZ_pX_reverse(ZZ_pX_c* x)
    #ZZ_pX_c* ZZ_pX_reverse_hi(ZZ_pX_c* x, long hi)
    #ZZ_pX_c* ZZ_pX_truncate(ZZ_pX_c* x, long m)
    #ZZ_pX_c* ZZ_pX_multiply_and_truncate(ZZ_pX_c* x, ZZ_pX_c* y, long m)
    #ZZ_pX_c* ZZ_pX_square_and_truncate(ZZ_pX_c* x, long m)
    #ZZ_pX_c* ZZ_pX_invert_and_truncate(ZZ_pX_c* x, long m)
    #ZZ_pX_c* ZZ_pX_multiply_mod(ZZ_pX_c* x, ZZ_pX_c* y,  ZZ_pX_c* modulus)
    #ZZ_p_c* ZZ_pX_trace_mod(ZZ_pX_c* x, ZZ_pX_c* y)
    char* ZZ_pX_trace_list(ZZ_pX_c* x)
    #ZZ_p_c* ZZ_pX_resultant(ZZ_pX_c* x, ZZ_pX_c* y)
    #ZZ_p_c* ZZ_pX_norm_mod(ZZ_pX_c* x, ZZ_pX_c* y)
    #ZZ_p_c* ZZ_pX_discriminant(ZZ_pX_c* x)
    #ZZ_pX_c* ZZ_pX_charpoly_mod(ZZ_pX_c* x, ZZ_pX_c* y)
    #ZZ_pX_c* ZZ_pX_minpoly_mod(ZZ_pX_c* x, ZZ_pX_c* y)
    #void ZZ_pX_preallocate_space(ZZ_pX_c* x, long n)

    void ZZ_pX_factor(ZZ_pX_c*** v, long** e, long* n, ZZ_pX_c* x, long verbose)
    void ZZ_pX_linear_roots(ZZ_p_c*** v, long* n, ZZ_pX_c* x)

    # The following are ZZ_pX functions written in ntl_wrap, used for padics.

    void ZZ_pX_conv_modulus(ZZ_pX_c fout, ZZ_pX_c fin, ZZ_pContext_c c)
    void ZZ_pX_min_val_coeff(long valuation, long index, ZZ_pX_c f, ZZ_c p)
    long ZZ_pX_get_val_coeff(ZZ_pX_c f, ZZ_c p, long i)
    void ZZ_pX_left_pshift(ZZ_pX_c x, ZZ_pX_c a, ZZ_c pn, ZZ_pContext_c c)
    void ZZ_pX_right_pshift(ZZ_pX_c x, ZZ_pX_c a, ZZ_c pn, ZZ_pContext_c c)
    void ZZ_pX_InvMod_newton_unram(ZZ_pX_c x, ZZ_pX_c a, ZZ_pX_Modulus_c F, ZZ_pContext_c cpn, ZZ_pContext_c cp)
    void ZZ_pX_InvMod_newton_ram(ZZ_pX_c x, ZZ_pX_c a, ZZ_pX_Modulus_c F, ZZ_pContext_c cpn)

    # The following are ZZ_pX functions written in ntl_wrap, used for padics.

    void ZZ_pX_conv_modulus(ZZ_pX_c fout, ZZ_pX_c fin, ZZ_pContext_c c)
    void ZZ_pX_min_val_coeff(long valuation, long index, ZZ_pX_c f, ZZ_c p)
    long ZZ_pX_get_val_coeff(ZZ_pX_c f, ZZ_c p, long i)
    void ZZ_pX_left_pshift(ZZ_pX_c x, ZZ_pX_c a, ZZ_c pn, ZZ_pContext_c c)
    void ZZ_pX_right_pshift(ZZ_pX_c x, ZZ_pX_c a, ZZ_c pn, ZZ_pContext_c c)
    void ZZ_pX_InvMod_newton(ZZ_pX_c x, ZZ_pX_c a, ZZ_pX_Modulus_c F, ZZ_pContext_c cpn, ZZ_pContext_c cp)
    void ZZ_pX_eis_shift(ZZ_pX_c x, ZZ_pX_c a, long n, ZZ_pX_Multiplier_c* low_shifter, ZZ_pX_Multiplier_c* high_shifter, ZZ_pX_Modulus_c modulus, ZZ_c p, ZZ_pContext_c cupper, ZZ_pContext_c clower)

    #### zz_p_c
    ctypedef struct zz_p_c "struct zz_p":
        void *rep
    void zz_p_construct "Construct<zz_p>"(void *mem)
    void zz_p_destruct "Destruct<zz_p>"(zz_p_c *mem)
    long zz_p_rep "rep"(zz_p_c x)
    long zz_p_isZero "IsZero"(zz_p_c x)
    void zz_p_set_from_long(zz_p_c x, long a)
    void zz_p_add "add"(zz_p_c x, zz_p_c a, zz_p_c b)
    void zz_p_sub "sub"(zz_p_c x, zz_p_c a, zz_p_c b)
    void zz_p_negate "negate"(zz_p_c x, zz_p_c a)
    void zz_p_mul "mul"(zz_p_c x, zz_p_c a, zz_p_c b)
    void zz_p_div "div"(zz_p_c x, zz_p_c a, zz_p_c b)
    void zz_p_inv "inv"(zz_p_c x, zz_p_c a)
    void zz_p_power "power"(zz_p_c x, zz_p_c a, long e)
    void zz_p_sqr "sqr"(zz_p_c x, zz_p_c a)
    void zz_p_clear "clear"(zz_p_c x)
    void zz_p_set_zero "clear"(zz_p_c x)
    void zz_p_set_one "set"(zz_p_c x)
    void zz_p_swap "swap"(zz_p_c x, zz_p_c y)

    bint NTL_zz_p_DOUBLE_EQUALS(zz_p_c x, zz_p_c y)

    #### zz_pContext_c
    ctypedef struct zz_pContext_c "struct zz_pContext":
        pass

    zz_pContext_c* zz_pContext_new "New<zz_pContext>"()
    zz_pContext_c* zz_pContext_construct "Construct<zz_pContext>"(void *mem)
    zz_pContext_c* zz_pContext_new_long "zz_pContext_new"(long p)
    zz_pContext_c* zz_pContext_construct_long "zz_pContext_construct"(void *mem, long p)
    void zz_pContext_destruct "Destruct<zz_pContext>"(zz_pContext_c *mem)
    void zz_pContext_delete "Delete<zz_pContext>"(zz_pContext_c *mem)

    void zz_pContext_restore(zz_pContext_c *c)


    #### zz_pX_c
    ctypedef struct zz_pX_c "struct zz_pX":
        void *rep
        void (* SetMaxLength)(long n)

    void zz_pX_construct "Construct<zz_pX>"(void *mem)
    void zz_pX_destruct "Destruct<zz_pX>"(zz_pX_c *mem)
    char* zz_pX_repr(zz_pX_c* x)
    void zz_pX_SetCoeff_long "SetCoeff"(zz_pX_c x, long i, long a)
    zz_p_c zz_pX_GetCoeff "coeff"(zz_pX_c x, long i)
    void zz_pX_add "add"(zz_pX_c x, zz_pX_c a, zz_pX_c b)
    void zz_pX_sub "sub"(zz_pX_c x, zz_pX_c a, zz_pX_c b)
    void zz_pX_mul "mul"(zz_pX_c x, zz_pX_c a, zz_pX_c b)
    void zz_pX_rmul "mul"(zz_pX_c x, zz_pX_c a, long b)
    void zz_pX_lmul "mul"(zz_pX_c x, long a, zz_pX_c b)
    long zz_pX_div "div"(zz_pX_c x, zz_pX_c a, zz_pX_c b)
    long zz_pX_divide "divide"(zz_pX_c x, zz_pX_c a, zz_pX_c b)
    void zz_pX_mod "rem"(zz_pX_c x, zz_pX_c a, zz_pX_c b)
    void zz_pX_divrem "DivRem"(zz_pX_c q, zz_pX_c r, zz_pX_c a, zz_pX_c b)
    void zz_pX_LeftShift "LeftShift"(zz_pX_c x, zz_pX_c a, long b)
    void zz_pX_RightShift "RightShift"(zz_pX_c x, zz_pX_c a, long b)
    void zz_pX_negate "negate"(zz_pX_c x, zz_pX_c a)
    zz_p_c zz_pX_LeadCoeff "LeadCoeff"(zz_pX_c x)
    zz_p_c zz_pX_ConstTerm "ConstTerm" (zz_pX_c x)
    void zz_pX_negate "negate"(zz_pX_c x, zz_pX_c a)
    void zz_pX_trunc "trunc"(zz_pX_c x, zz_pX_c a, long n) ## x = a % X^n
    void zz_pX_MulTrunc "MulTrunc"(zz_pX_c x, zz_pX_c a, zz_pX_c b, long n)
    void zz_pX_SqrTrunc "SqrTrunc"(zz_pX_c x, zz_pX_c a, long n)
    void zz_pX_InvTrunc "InvTrunc"(zz_pX_c x, zz_pX_c a, long n)
    void zz_pX_sqr "sqr"(zz_pX_c x, zz_pX_c a)
    void zz_pX_power "power"(zz_pX_c x, zz_pX_c a, long e)
    void zz_pX_clear "clear"(zz_pX_c x)
    void zz_pX_SetX "SetX"(zz_pX_c x)
    bint zz_pX_IsX "IsX"(zz_pX_c x)
    bint zz_pX_IsZero "IsZero"(zz_pX_c x)
    bint zz_pX_IsOne "IsOne"(zz_pX_c x)
    long zz_pX_deg "deg"(zz_pX_c x)
    zz_pX_c zz_pX_zero "zz_pX::zero"()
    void zz_pX_diff "diff"(zz_pX_c x, zz_pX_c a)
    void zz_pX_reverse "reverse"(zz_pX_c x, zz_pX_c a)
    void zz_pX_eval "eval" (zz_p_c fa, zz_pX_c f, zz_p_c a)
    void zz_pX_MakeMonic "MakeMonic"(zz_pX_c x)

    long NTL_SP_BOUND
    bint NTL_zz_pX_DOUBLE_EQUALS(zz_pX_c x, zz_pX_c y)


    #### ZZ_pEContext_c
    ctypedef struct ZZ_pEContext_c "struct ZZ_pEContext":
        void (*restore)()

    ZZ_pEContext_c* ZZ_pEContext_new "New<ZZ_pEContext>"()
    ZZ_pEContext_c* ZZ_pEContext_construct "Construct<ZZ_pEContext>"(void *mem)
    ZZ_pEContext_c* ZZ_pEContext_new_ZZ_pX "ZZ_pEContext_new"(ZZ_pX_c* p)
    ZZ_pEContext_c* ZZ_pEContext_construct_ZZ_pX "ZZ_pEContext_construct"(void *mem, ZZ_pX_c* p)
    void ZZ_pEContext_destruct "Destruct<ZZ_pEContext>"(ZZ_pEContext_c *mem)
    void ZZ_pEContext_delete "Delete<ZZ_pEContext>"(ZZ_pEContext_c *mem)

    void ZZ_pEContext_restore(ZZ_pEContext_c *c)

    #### ZZ_pE_c
    ctypedef struct ZZ_pE_c "struct ZZ_pE":
        pass

    # Some boiler-plate
    ZZ_pE_c* ZZ_pE_new "New<ZZ_pE>"()
    ZZ_pE_c* ZZ_pE_construct "Construct<ZZ_pE>"(void *mem)
    void ZZ_pE_destruct "Destruct<ZZ_pE>"(ZZ_pE_c *mem)
    void ZZ_pE_delete "Delete<ZZ_pE>"(ZZ_pE_c *mem)
    void ZZ_pE_from_str "_from_str<ZZ_pE>"(ZZ_pE_c* dest, char* s)
    object ZZ_pE_to_PyString "_to_PyString<ZZ_pE>"(ZZ_pE_c *x)
    int ZZ_pE_equal "_equal<ZZ_pE>"(ZZ_pE_c x, ZZ_pE_c y)

    #ZZ_pE_c* str_to_ZZ_pE(char* s)
    #void del_ZZ_pE(ZZ_pE_c* n)
    #void ZZ_pE_to_str(char** s, ZZ_pE_c* x)
    void ZZ_pE_add "add"( ZZ_pE_c x, ZZ_pE_c a, ZZ_pE_c b)
    void ZZ_pE_add_long "add"( ZZ_pE_c x, ZZ_pE_c a, long b)
    void ZZ_pE_add_ZZ_p "add"( ZZ_pE_c x, ZZ_pE_c a, ZZ_p_c b)
    void ZZ_pE_sub "sub"( ZZ_pE_c x, ZZ_pE_c a, ZZ_pE_c b)
    void ZZ_pE_sub_long "add"( ZZ_pE_c x, ZZ_pE_c a, long b)
    void ZZ_pE_sub_ZZ_p "sub"( ZZ_pE_c x, ZZ_pE_c a, ZZ_p_c b)
    void ZZ_pE_mul "mul"( ZZ_pE_c x, ZZ_pE_c a, ZZ_pE_c b)
    void ZZ_pE_mul_long "mul"( ZZ_pE_c x, ZZ_pE_c a, long b)
    void ZZ_pE_mul_ZZ_p "mul"( ZZ_pE_c x, ZZ_pE_c a, ZZ_p_c b)
    void ZZ_pE_negate "negate"(ZZ_pE_c x, ZZ_pE_c a)
    void ZZ_pE_power "power"(ZZ_pE_c t, ZZ_pE_c x, long e)
    int ZZ_pE_IsOne "IsOne"(ZZ_pE_c x)
    int ZZ_pE_IsZero "IsZero"(ZZ_pE_c x)
    ZZ_pX_c ZZ_pE_rep "rep"(ZZ_pE_c z)
    #void ntl_ZZ_pE_set_modulus(ZZ_pX_c* x)
    #int ZZ_pE_eq( ZZ_pE_c* x,  ZZ_pE_c* y)
    #int ZZ_pE_eq_ZZ_p( ZZ_pE_c* x, ZZ_p_c* y)
    void ZZ_pE_inv "inv"( ZZ_pE_c x, ZZ_pE_c a) # raises an error if a not invertible
    void ZZ_pE_div "div"( ZZ_pE_c x, ZZ_pE_c a, ZZ_pE_c b) # raises an error if b not invertible
    void ZZ_pE_div_ZZ_p "div"( ZZ_pE_c x, ZZ_pE_c a, ZZ_p_c b) # raises an error if b not invertible
    void ZZ_pE_ZZ_p_div "div"( ZZ_pE_c x, ZZ_p_c a, ZZ_pE_c b) # raises an error if b not invertible
    void ZZ_pE_random "random"( ZZ_pE_c x)
    void ZZ_pE_trace "trace"( ZZ_p_c x, ZZ_pE_c a)
    void ZZ_pE_norm "norm"( ZZ_p_c x, ZZ_pE_c a)

    ZZ_pE_c long_to_ZZ_pE "to_ZZ_pE"(long i)
    ZZ_pE_c ZZ_to_ZZ_pE "to_ZZ_pE"(ZZ_c i)
    ZZ_pE_c ZZ_p_to_ZZ_pE "to_ZZ_pE"(ZZ_p_c i)
    ZZ_pE_c ZZ_pX_to_ZZ_pE "to_ZZ_pE"(ZZ_pX_c i)
    ZZ_pX_c ZZ_pE_to_ZZ_pX(ZZ_pE_c x)

    #ZZ_pX_c rep(ZZ_pE_c x)

    #### vec_ZZ_pE_c
    ctypedef struct vec_ZZ_pE_c "struct vec_ZZ_pE":
        pass

    #### ZZ_pEX_c
    ctypedef struct ZZ_pEX_c "struct ZZ_pEX":
        void *rep
        void (* SetMaxLength)(long n)

    ZZ_pEX_c* ZZ_pEX_new "New<ZZ_pEX>"()
    ZZ_pEX_c* ZZ_pEX_construct "Construct<ZZ_pEX>"(void *mem)
    void ZZ_pEX_destruct "Destruct<ZZ_pEX>"(ZZ_pEX_c *mem)
    void ZZ_pEX_delete "Delete<ZZ_pEX>"(ZZ_pEX_c *mem)
    void ZZ_pEX_from_str "_from_str<ZZ_pEX>"(ZZ_pEX_c* dest, char* s)
    object ZZ_pEX_to_PyString "_to_PyString<ZZ_pEX>"(ZZ_pEX_c *x)
    #int ZZ_pEX_equal "_equal<ZZX>"(ZZ_pEX_c x, ZZ_pEX_c y)

    #ZZ_pEX_c* str_to_ZZ_pEX(char* s)
    #char* ZZ_pEX_to_str(ZZ_pEX_c* x)

    long ZZ_pEX_equal "operator=="(ZZ_pEX_c a, ZZ_pEX_c b)
    long ZZ_pEX_IsZero "IsZero"(ZZ_pEX_c a)
    long ZZ_pEX_IsOne "IsOne"(ZZ_pEX_c a)

    void ZZ_pEX_add "add"( ZZ_pEX_c x, ZZ_pEX_c a, ZZ_pEX_c b)
    void ZZ_pEX_add_ZZ_p "add"( ZZ_pEX_c x, ZZ_pEX_c a, ZZ_p_c b)
    void ZZ_pEX_add_ZZ_pE "add"( ZZ_pEX_c x, ZZ_pEX_c a, ZZ_pE_c b)
    void ZZ_pEX_sub "sub"( ZZ_pEX_c x, ZZ_pEX_c a, ZZ_pEX_c b)
    void ZZ_pEX_sub_ZZ_p "sub"( ZZ_pEX_c x, ZZ_pEX_c a, ZZ_p_c b)
    void ZZ_pEX_sub_ZZ_pE "sub"( ZZ_pEX_c x, ZZ_pEX_c a, ZZ_pE_c b)
    void ZZ_pEX_negate "negate"(ZZ_pEX_c x, ZZ_pEX_c a)

    void ZZ_pEX_mul "mul"( ZZ_pEX_c x, ZZ_pEX_c a, ZZ_pEX_c b)
    void ZZ_pEX_mul_long "mul"( ZZ_pEX_c x, ZZ_pEX_c a, long b)
    void ZZ_pEX_mul_ZZ_p "mul"( ZZ_pEX_c x, ZZ_pEX_c a, ZZ_p_c b)
    void ZZ_pEX_mul_ZZ_pE "mul"( ZZ_pEX_c x, ZZ_pEX_c a, ZZ_pE_c b)
    void ZZ_pEX_sqr "sqr"( ZZ_pEX_c x, ZZ_pEX_c a)
    void ZZ_pEX_power "power"( ZZ_pEX_c x, ZZ_pEX_c a, long e)

    void ZZ_pEX_LeftShift "LeftShift"(ZZ_pEX_c x, ZZ_pEX_c a, long n)
    void ZZ_pEX_RightShift "RightShift"(ZZ_pEX_c x, ZZ_pEX_c a, long n)

    void ZZ_pEX_DivRem "DivRem"(ZZ_pEX_c q, ZZ_pEX_c r, ZZ_pEX_c a, ZZ_pEX_c b)
    void ZZ_pEX_div_ZZ_pEX "div"(ZZ_pEX_c q, ZZ_pEX_c a, ZZ_pEX_c b)
    void ZZ_pEX_div_ZZ_pE "div"(ZZ_pEX_c q, ZZ_pEX_c a, ZZ_pE_c b)
    void ZZ_pEX_div_ZZ_p "div"(ZZ_pEX_c q, ZZ_pEX_c a, ZZ_p_c b)
    void ZZ_pEX_div_long "div"( ZZ_pEX_c x, ZZ_pEX_c a, long b)
    void ZZ_pEX_rem "rem"(ZZ_pEX_c r, ZZ_pEX_c a, ZZ_pEX_c b)
    long ZZ_pEX_divide "divide"( ZZ_pEX_c x, ZZ_pEX_c a, ZZ_pEX_c b)

    void ZZ_pEX_GCD "GCD"(ZZ_pEX_c r, ZZ_pEX_c a, ZZ_pEX_c b)
    void ZZ_pEX_XGCD "XGCD"(ZZ_pEX_c d, ZZ_pEX_c s, ZZ_pEX_c t, ZZ_pEX_c a, ZZ_pEX_c b)

    long ZZ_pEX_deg "deg"( ZZ_pEX_c x )
    ZZ_pE_c ZZ_pEX_coeff "coeff"(ZZ_pEX_c a, long i)
    ZZ_pE_c ZZ_pEX_LeadCoeff "LeadCoeff"(ZZ_pEX_c a)
    void ZZ_pEX_SetCoeff "SetCoeff"(ZZ_pEX_c x, long i, ZZ_pE_c a)
    void ZZ_pEX_SetCoeff_ZZ_p "SetCoeff"(ZZ_pEX_c x, long i, ZZ_p_c a)
    void ZZ_pEX_SetCoeff_long "SetCoeff"(ZZ_pEX_c x, long i, long a)
    void ZZ_pEX_SetCoeff_one "SetCoeff"(ZZ_pEX_c x, long i)
    void ZZ_pEX_SetX "SetX"(ZZ_pEX_c x)
    long ZZ_pEX_IsX "IsX"(ZZ_pEX_c x)
    void ZZ_pEX_diff "diff"(ZZ_pEX_c x, ZZ_pEX_c a)
    void ZZ_pEX_MakeMonic "MakeMonic"(ZZ_pEX_c x)
    void ZZ_pEX_reverse_hi "reverse"(ZZ_pEX_c x, ZZ_pEX_c a, long hi)
    void ZZ_pEX_reverse "reverse"(ZZ_pEX_c x, ZZ_pEX_c a)

    void ZZ_pEX_random "random"(ZZ_pEX_c x, long n)

    void ZZ_pEX_BuildFromRoots "BuildFromRoots"(ZZ_pEX_c x, vec_ZZ_pE_c a)
    void ZZ_pEX_eval "eval"(ZZ_pE_c b, ZZ_pEX_c f, ZZ_pE_c a)
    void ZZ_pEX_eval_vec "eval"(vec_ZZ_pE_c b, ZZ_pEX_c f, vec_ZZ_pE_c a)
    void ZZ_pEX_interpolate "interpolate"(ZZ_pEX_c f, vec_ZZ_pE_c a, vec_ZZ_pE_c b)

    void ZZ_pEX_trunc "trunc"(ZZ_pEX_c x, ZZ_pEX_c a, long n)
    void ZZ_pEX_MulTrunc "MulTrunc"(ZZ_pEX_c x, ZZ_pEX_c a, ZZ_pEX_c b, long n)
    void ZZ_pEX_SqrTrunc "SqrTrunc"(ZZ_pEX_c x, ZZ_pEX_c a, long n)
    void ZZ_pEX_InvTrunc "InvTrunc"(ZZ_pEX_c x, ZZ_pEX_c a, long n)

    void ZZ_pEX_MulMod "MulMod"(ZZ_pEX_c x, ZZ_pEX_c a, ZZ_pEX_c b, ZZ_pEX_c f)
    void ZZ_pEX_SqrMod "SqrMod"(ZZ_pEX_c x, ZZ_pEX_c a, ZZ_pEX_c f)
    void ZZ_pEX_MulByXMod "MulByXMod"(ZZ_pEX_c x, ZZ_pEX_c a, ZZ_pEX_c f)
    void ZZ_pEX_InvMod "InvMod"(ZZ_pEX_c x, ZZ_pEX_c a, ZZ_pEX_c f)
    long ZZ_pEX_InvModStatus "InvModStatus"(ZZ_pEX_c x, ZZ_pEX_c a, ZZ_pEX_c f)

    ctypedef struct ZZ_pEX_Modulus_c "struct ZZ_pEXModulus":
        ZZ_pEX_c val()
    ZZ_pEX_Modulus_c* ZZ_pEX_Modulus_new "New<ZZ_pEXModulus>"()
    ZZ_pEX_Modulus_c* ZZ_pEX_Modulus_construct "Construct<ZZ_pEXModulus>"(void *mem)
    void ZZ_pEX_Modulus_destruct "Destruct<ZZ_pEXModulus>"(ZZ_pEX_Modulus_c *mem)
    void ZZ_pEX_Modulus_delete "Delete<ZZ_pEXModulus>"(ZZ_pEX_Modulus_c *mem)
    void ZZ_pEX_Modulus_from_str "_from_str<ZZ_pEXModulus>"(ZZ_pEX_Modulus_c* dest, char* s)
    void ZZ_pEX_Modulus_build "build"(ZZ_pEX_Modulus_c F, ZZ_pEX_c f)
    long ZZ_pEX_Modulus_deg "deg"(ZZ_pEX_Modulus_c F)

    void ZZ_pEX_MulMod_pre "MulMod"(ZZ_pEX_c x, ZZ_pEX_c a, ZZ_pEX_c b, ZZ_pEX_Modulus_c F)
    void ZZ_pEX_SqrMod_pre "SqrMod"(ZZ_pEX_c x, ZZ_pEX_c a, ZZ_pEX_Modulus_c F)
    void ZZ_pEX_PowerMod_pre "PowerMod"(ZZ_pEX_c x, ZZ_pEX_c a, long e, ZZ_pEX_Modulus_c F)
    void ZZ_pEX_PowerMod_ZZ_pre "PowerMod"(ZZ_pEX_c x, ZZ_pEX_c a, ZZ_c e, ZZ_pEX_Modulus_c F)
    void ZZ_pEX_PowerXMod_pre "PowerXMod"(ZZ_pEX_c x, long e, ZZ_pEX_Modulus_c F)
    void ZZ_pEX_PowerXMod_ZZ_pre "PowerXMod"(ZZ_pEX_c x, ZZ_c e, ZZ_pEX_Modulus_c F)
    void ZZ_pEX_rem_pre "rem"(ZZ_pEX_c x, ZZ_pEX_c a, ZZ_pEX_Modulus_c F)
    void ZZ_pEX_DivRem_pre "DivRem"(ZZ_pEX_c q, ZZ_pEX_c r, ZZ_pEX_c a, ZZ_pEX_Modulus_c F)
    void ZZ_pEX_div_pre "div"(ZZ_pEX_c q, ZZ_pEX_c a, ZZ_pEX_Modulus_c F)

    void ZZ_pEX_MinPolyMod "MinPolyMod"(ZZ_pEX_c h, ZZ_pEX_c g, ZZ_pEX_c f)
    void ZZ_pEX_MinPolyMod_pre "MinPolyMod"(ZZ_pEX_c h, ZZ_pEX_c g, ZZ_pEX_Modulus_c F)

    void ZZ_pEX_TraceMod "TraceMod"(ZZ_pE_c x, ZZ_pEX_c a, ZZ_pEX_c f)
    void ZZ_pEX_TraceMod_pre "TraceMod"(ZZ_pE_c x, ZZ_pEX_c a, ZZ_pEX_Modulus_c F)
    void ZZ_pEX_TraceVec "TraceVec"(vec_ZZ_pE_c x, ZZ_pEX_c f)
    void ZZ_pEX_NormMod "NormMod"(ZZ_pE_c x, ZZ_pEX_c a, ZZ_pEX_c f)
    void ZZ_pEX_resultant "resultant"(ZZ_pE_c x, ZZ_pEX_c a, ZZ_pEX_c b)

    void ZZ_pEX_clear "clear"(ZZ_pEX_c x)
    void ZZ_pEX_set "set"(ZZ_pEX_c x)

    void ZZ_pEX_conv_modulus(ZZ_pEX_c fout, ZZ_pEX_c fin, ZZ_pContext_c c)


    #### mat_ZZ_c
    ctypedef struct mat_ZZ_c "struct mat_ZZ":
        pass

    # Some boiler-plate
    mat_ZZ_c* mat_ZZ_construct "Construct<mat_ZZ>"(void *mem)
    void mat_ZZ_destruct "Destruct<mat_ZZ>"(mat_ZZ_c *mem)
    void mat_ZZ_delete "Delete<mat_ZZ>"(mat_ZZ_c *mem)
    object mat_ZZ_to_PyString "_to_PyString<mat_ZZ>"(mat_ZZ_c *x)
    int mat_ZZ_equal "_equal<mat_ZZ>"(mat_ZZ_c x, mat_ZZ_c y)

    void mat_ZZ_mul "mul"( mat_ZZ_c x, mat_ZZ_c a, mat_ZZ_c b)
    void mat_ZZ_add "add"( mat_ZZ_c x, mat_ZZ_c a, mat_ZZ_c b)
    void mat_ZZ_sub "sub"( mat_ZZ_c x, mat_ZZ_c a, mat_ZZ_c b)
    void mat_ZZ_power "power"( mat_ZZ_c x, mat_ZZ_c a, long e)
    void mat_ZZ_CharPoly "CharPoly"(ZZX_c r, mat_ZZ_c m)
    void mat_ZZ_SetDims(mat_ZZ_c* mZZ, long nrows, long ncols)

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

    #### GF2_c
    ctypedef struct GF2_c "struct GF2":
        pass

    GF2_c* GF2_new "New<GF2>"()
    GF2_c* GF2_construct "Construct<GF2>"(void *mem)
    void GF2_destruct "Destruct<GF2>"(GF2_c *mem)
    void GF2_delete "Delete<GF2>"(GF2_c *mem)
    void GF2_from_str "_from_str<GF2>"(GF2_c* dest, char* s)
    object GF2_to_PyString "_to_PyString<GF2>"(GF2_c *x)
    int GF2_equal "_equal<GF2>"(GF2_c x, GF2_c y)
    int GF2_IsOne "IsOne"(GF2_c x)
    int GF2_IsZero "IsZero"(GF2_c x)

    void GF2_add "add"( GF2_c x, GF2_c a, GF2_c b)
    void GF2_sub "sub"( GF2_c x, GF2_c a, GF2_c b)
    void GF2_mul "mul"( GF2_c x, GF2_c a, GF2_c b)
    void GF2_div "div"( GF2_c x, GF2_c a, GF2_c b)
    void GF2_negate "negate"(GF2_c x, GF2_c a)
    void GF2_power "power"(GF2_c t, GF2_c x, long e)
    long GF2_deg "deg"(GF2_c x)

    void GF2_conv_long "conv" (GF2_c x, long i)
    long GF2_conv_to_long "rep" (GF2_c x)

    #### GF2X_c

    ctypedef struct GF2X_c "struct GF2X":
        pass

    long *GF2XHexOutput_c "(&GF2X::HexOutput)" # work-around for Cython bug

    GF2X_c* GF2X_new "New<GF2X>"()
    GF2X_c* GF2X_construct "Construct<GF2X>"(void *mem)
    void GF2X_destruct "Destruct<GF2X>"(GF2X_c *mem)
    void GF2X_delete "Delete<GF2X>"(GF2X_c *mem)
    void GF2X_from_str "_from_str<GF2X>"(GF2X_c* dest, char* s)
    object GF2X_to_PyString "_to_PyString<GF2X>"(GF2X_c *x)
    int GF2X_equal "_equal<GF2X>"(GF2X_c x, GF2X_c y)
    int GF2X_IsOne "IsOne"(GF2X_c x)
    int GF2X_IsZero "IsZero"(GF2X_c x)

    void GF2X_add "add"( GF2X_c x, GF2X_c a, GF2X_c b)
    void GF2X_sub "sub"( GF2X_c x, GF2X_c a, GF2X_c b)
    void GF2X_mul "mul"( GF2X_c x, GF2X_c a, GF2X_c b)
    void GF2X_negate "negate"(GF2X_c x, GF2X_c a)
    void GF2X_power "power"(GF2X_c t, GF2X_c x, long e)
    long GF2X_deg "deg"(GF2X_c x)

    void GF2X_conv_long "conv" (GF2X_c x, long a)
    void GF2X_conv_GF2 "conv" (GF2X_c x, GF2_c a)

    void GF2X_LeftShift "LeftShift"( GF2X_c r, GF2X_c a, long offset)
    void GF2X_RightShift "RightShift"( GF2X_c r, GF2X_c a, long offset)

    void GF2X_DivRem "DivRem"(GF2X_c q, GF2X_c r, GF2X_c a, GF2X_c b)
    void GF2X_div "div" (GF2X_c q, GF2X_c a, GF2X_c b)
    void GF2X_rem "rem" (GF2X_c r, GF2X_c a, GF2X_c b)
    long GF2X_divide "divide"(GF2X_c q, GF2X_c a, GF2X_c b)

    GF2X_c GF2X_GCD "GCD" (GF2X_c a, GF2X_c b)
    void GF2X_XGCD "XGCD" (GF2X_c r, GF2X_c s, GF2X_c t, GF2X_c a, GF2X_c b)

    void GF2XFromBytes(GF2X_c a, unsigned char *p, long n)
    void BytesFromGF2X "BytesFromGF2X" (unsigned char *p, GF2X_c a, long n)

    GF2_c GF2X_coeff "coeff"(GF2X_c a, long i)
    GF2_c GF2X_LeadCoeff "LeadCoeff"(GF2X_c a)
    GF2_c GF2X_ConstTerm "ConstTerm"(GF2X_c a)
    void GF2X_SetCoeff "SetCoeff"(GF2X_c x, long i, GF2_c a)

    GF2X_c GF2X_diff "diff"(GF2X_c a)
    GF2X_c GF2X_reverse "reverse"(GF2X_c a, long hi)

    long GF2X_weight "weight"(GF2X_c a)
    long GF2X_NumBits "NumBits" (GF2X_c a)
    long GF2X_NumBytes "NumBytes"(GF2X_c a)

    #### GF2XModulus_c
    ctypedef struct GF2XModulus_c "struct GF2XModulus":
        pass

    GF2X_c GF2XModulus_GF2X "GF2X" (GF2XModulus_c m)

    #### GF2EContext_c

    ctypedef struct GF2EContext_c "struct GF2EContext":
        void (*restore)()

    GF2EContext_c* GF2EContext_new "New<GF2EContext>"()
    GF2EContext_c* GF2EContext_construct "Construct<GF2EContext>"(void *mem)
    GF2EContext_c* GF2EContext_new_GF2X "GF2EContext_new"(GF2X_c* p)
    GF2EContext_c* GF2EContext_construct_GF2X "GF2EContext_construct"(void *mem, GF2X_c* p)
    void GF2EContext_destruct "Destruct<GF2EContext>"(GF2EContext_c *mem)
    void GF2EContext_delete "Delete<GF2EContext>"(GF2EContext_c *mem)

    #### GF2E_c
    ctypedef struct GF2E_c "struct GF2E":
        pass

    void GF2E_init "GF2E::init"(GF2X_c x)
    long GF2E_degree "GF2E::degree"()
    GF2XModulus_c GF2E_modulus "GF2E::modulus"()

    GF2E_c* GF2E_new "New<GF2E>"()
    GF2E_c* GF2E_construct "Construct<GF2E>"(void *mem)
    void GF2E_destruct "Destruct<GF2E>"(GF2E_c *mem)
    void GF2E_delete "Delete<GF2E>"(GF2E_c *mem)
    void GF2E_from_str "_from_str<GF2E>"(GF2E_c* dest, char* s)
    object GF2E_to_PyString "_to_PyString<GF2E>"(GF2E_c *x)
    int GF2E_equal "_equal<GF2E>"(GF2E_c x, GF2E_c y)
    int GF2E_IsOne "IsOne"(GF2E_c x)
    int GF2E_IsZero "IsZero"(GF2E_c x)

    void GF2E_add "add"( GF2E_c x, GF2E_c a, GF2E_c b)
    void GF2E_sub "sub"( GF2E_c x, GF2E_c a, GF2E_c b)
    void GF2E_mul "mul"( GF2E_c x, GF2E_c a, GF2E_c b)
    void GF2E_div "div"( GF2E_c x, GF2E_c a, GF2E_c b)
    void GF2E_power "power"(GF2E_c t, GF2E_c x, long e)
    long GF2E_deg "deg"(GF2E_c x)

    void GF2E_conv_GF2X "conv" (GF2E_c out, GF2X_c inp)
    void GF2E_conv_long "conv" (GF2E_c out, long inp)
    void GF2E_conv_ZZ "conv" (GF2E_c out, ZZ_c inp)
    void GF2E_conv_GF2 "conv" (GF2E_c out, GF2_c inp)
    GF2X_c GF2E_rep "rep"(GF2E_c x)

    GF2E_c GF2E_random "random_GF2E"()

    GF2_c GF2E_trace "trace"(GF2E_c x)

    #### GF2EX_c
    ctypedef struct GF2EX_c "struct GF2EX":
        pass

    GF2EX_c* GF2EX_new "New<GF2EX>"()
    GF2EX_c* GF2EX_construct "Construct<GF2EX>"(void *mem)
    void GF2EX_destruct "Destruct<GF2EX>"(GF2EX_c *mem)
    void GF2EX_delete "Delete<GF2EX>"(GF2EX_c *mem)
    void GF2EX_from_str "_from_str<GF2EX>"(GF2EX_c* dest, char* s)
    object GF2EX_to_PyString "_to_PyString<GF2EX>"(GF2EX_c *x)
    int GF2EX_equal "_equal<GF2EX>"(GF2EX_c x, GF2EX_c y)
    void GF2EX_add "add"( GF2EX_c x, GF2EX_c a, GF2EX_c b)
    void GF2EX_sub "sub"( GF2EX_c x, GF2EX_c a, GF2EX_c b)
    void GF2EX_mul "mul"( GF2EX_c x, GF2EX_c a, GF2EX_c b)
    void GF2EX_negate "negate"(GF2EX_c x, GF2EX_c a)
    void GF2EX_power "power"(GF2EX_c t, GF2EX_c x, long e)
    int GF2EX_IsOne "IsOne"(GF2EX_c x)
    int GF2EX_IsZero "IsZero"(GF2EX_c x)

    #### vec_GF2E_c
    ctypedef struct vec_GF2E_c "struct vec_GF2E":
        pass

    vec_GF2E_c* vec_GF2E_new "New<vec_GF2E>"()
    vec_GF2E_c* vec_GF2E_construct "Construct<vec_GF2E>"(void *mem)
    void vec_GF2E_destruct "Destruct<vec_GF2E>"(vec_GF2E_c *mem)
    void vec_GF2E_delete "Delete<vec_GF2E>"(vec_GF2E_c *mem)
    void vec_GF2E_from_str "_from_str<vec_GF2E>"(vec_GF2E_c* dest, char* s)
    object vec_GF2E_to_PyString "_to_PyString<vec_GF2E>"(vec_GF2E_c *x)

    #### mat_GF2E_c
    ctypedef struct mat_GF2E_c "struct mat_GF2E":
        void (*SetDims)(long nrows, long ncols)
        long (*NumRows)()
        long (*NumCols)()
        GF2E_c (*get "operator()") (long i, long j)

    mat_GF2E_c* mat_GF2E_new "New<mat_GF2E>"()
    mat_GF2E_c* mat_GF2E_construct "Construct<mat_GF2E>"(void *mem)
    void mat_GF2E_destruct "Destruct<mat_GF2E>"(mat_GF2E_c *mem)
    void mat_GF2E_delete "Delete<mat_GF2E>"(mat_GF2E_c *mem)
    void mat_GF2E_from_str "_from_str<mat_GF2E>"(mat_GF2E_c* dest, char* s)
    object mat_GF2E_to_PyString "_to_PyString<mat_GF2E>"(mat_GF2E_c *x)
    int mat_GF2E_equal "_equal<mat_GF2E>"(mat_GF2E_c x, mat_GF2E_c y)
    void mat_GF2E_add "add"( mat_GF2E_c x, mat_GF2E_c a, mat_GF2E_c b)
    void mat_GF2E_sub "sub"( mat_GF2E_c x, mat_GF2E_c a, mat_GF2E_c b)
    void mat_GF2E_mul "mul"( mat_GF2E_c x, mat_GF2E_c a, mat_GF2E_c b)
    void mat_GF2E_negate "negate"(mat_GF2E_c x, mat_GF2E_c a)
    void mat_GF2E_power "power"(mat_GF2E_c t, mat_GF2E_c x, long e)
    GF2E_c mat_GF2E_determinant "determinant"(mat_GF2E_c m)
    void mat_GF2E_transpose "transpose"(mat_GF2E_c r, mat_GF2E_c m)
    long mat_GF2E_IsZero "IsZero"(mat_GF2E_c x)
    void mat_GF2E_setitem(mat_GF2E_c* x, int i, int j, GF2E_c* z)

    long mat_GF2E_gauss "gauss"(mat_GF2E_c A, long w)
    void mat_GF2E_solve "solve"(GF2E_c d, vec_GF2E_c X, mat_GF2E_c A, vec_GF2E_c b)
    void mat_GF2E_inv "inv" (mat_GF2E_c X, mat_GF2E_c A)


    long mat_GF2E_IsIdent "IsIdent"(mat_GF2E_c A, long n)
    long mat_GF2E_IsDiag "IsDiag"(mat_GF2E_c A, long n, GF2E_c d)


    void mat_GF2E_image "image"(mat_GF2E_c X, mat_GF2E_c A)
    void mat_GF2E_kernel "kernel" (mat_GF2E_c X, mat_GF2E_c A)

    void vec_GF2E_conv_mat_GF2E "conv" (vec_GF2E_c out, mat_GF2E_c inp)
    void mat_GF2E_conv_vec_GF2E(mat_GF2E_c out, vec_GF2E_c inp)

cdef extern from "ZZ_pylong.h":
    object ZZ_get_pylong(ZZ_c z)
    int ZZ_set_pylong(ZZ_c z, object ll)
