cdef extern from "gmp.h":

    ### Type Declarations ###

    # Underlying typedefs
    ctypedef unsigned long mp_limb_t
    ctypedef unsigned long mp_bitcnt_t
    ctypedef long mp_size_t
    ctypedef long mp_exp_t

    ctypedef long mp_limb_t
    ctypedef mp_limb_t* mp_ptr

    # User-level types
    ctypedef void* mpz_t
    ctypedef void* mpq_t
    ctypedef void* mpf_t
    ctypedef void* gmp_randstate_t

    # This internal structure is not guaranteed to stay the same with future releases of gmp.
    # We declare mpz_t as a void* because one is supposed to treat it as a black box.
    ctypedef struct __mpz_struct:
        int _mp_alloc
        int _mp_size
        mp_ptr _mp_d
