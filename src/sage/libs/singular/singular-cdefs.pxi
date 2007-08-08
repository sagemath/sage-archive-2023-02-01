"""
C level declarations of symbols in the SINGULAR libary.

AUTHOR: Martin Albrecht <malb@informatik.uni-bremen.de>

NOTE: our ring, poly etc. types are not the Singular ring, poly,
etc. types. They are deferences. So a Singular ring is a ring* here.

"""
include "../../ext/cdefs.pxi"

cdef extern from "stdsage.h":
    ctypedef void PyObject

    # Global tuple -- useful optimization
    void init_global_empty_tuple()
    object PY_NEW(object t)
    void* PY_TYPE(object o)
    int PY_TYPE_CHECK(object o, object t)
    object IS_INSTANCE(object o, object t)


# Shared Library Loading,
cdef extern from "dlfcn.h":
    void *dlopen(char *, long)
    char *dlerror()

cdef extern from "stdlib.h":
    void *calloc(size_t nmemb, size_t size)
    void free(void *ptr)
    void delete "delete" (void *ptr)

cdef extern from "libsingular.h":

    # all possible ring orders per block
    ctypedef struct number "snumber":
        mpz_t z
        mpz_t n
        int s

    ctypedef struct napoly "polyrec"

    cdef enum tHomog:
        isNotHomog
        isHomog
        testHomog

    cdef enum rRingOrder_t:
        ringorder_no
        ringorder_a
        ringorder_a64 # for int64 weights
        ringorder_c
        ringorder_C
        ringorder_M
        ringorder_S
        ringorder_s
        ringorder_lp
        ringorder_dp
        ringorder_rp
        ringorder_Dp
        ringorder_wp
        ringorder_Wp
        ringorder_ls
        ringorder_ds
        ringorder_Ds
        ringorder_ws
        ringorder_Ws
        ringorder_L

    # This is the basis ring datatype of SINGULAR
    # Please note: This is not the SINGULAR ring, it has one pointer
    # layer less.

    ctypedef struct ring "ip_sring":
        int  *order  # array of orderings
        int  *block0 # starting pos
        int  *block1 # ending pos
        int  **wvhdl
        int  OrdSgn
        int  ShortOut
        int  CanShortOut
        number *minpoly
        char **names
        char **parameter
        ring *algring
        short N
        short P
        int ch

    # This is the basic polynomial datatype of SINGULAR

    # Please note: This is not the SINGULAR poly, it has one pointer
    # layer less.

    ctypedef struct poly "polyrec":
        poly *next

    # This is the basic ideal/matrix datatype of SINGULAR

    # Please note: This is not the SINGULAR ideal, it has one pointer
    # layer less.

    ctypedef struct ideal "ip_sideal":
        poly **m
        long rank
        int nrows
        int ncols

    # comment from the SINGULAR code:

    #  'keiner (ausser obachman) darf das folgenden benutzen !!!'
    #  in English: 'nobody (except obachman) may use the following!!!'
    #
    # I feel so obachman today.

    ctypedef struct intvec:
        int *(*ivGetVec)()
        int (*rows)()
        int (*cols)()

    # oMalloc Bins
    ctypedef struct omBin "omBin_s"


    # SINGULAR Init
    # ------------------------
    void feInitResources(char *name)

    void rChangeCurrRing(ring *r)
    cdef ring *currRing
    cdef omBin *rnumber_bin
    cdef omBin *sip_sring_bin
    cdef int (*pLDeg)(poly *p, int *l, ring *r)

    int siInit(char *)

    # OMalloc
    void *omAlloc0(size_t size)
    void *omAllocBin(omBin *bin)
    void *omAlloc0Bin(omBin *bin)
    char *omStrDup(char *)
    void omFree(void *)


    # rDefault accepts the characteristic, the number of variables,
    # and an array of variable names.
    ring *rDefault(int char, int nvars, char **names)

    void rDelete(ring *r) # Destructor
    int rChar(ring *r) # Returns the characteristic of the ring
    char* rRingVar(short i, ring *r) # Returns the name of the i-th variable of the ring r.

    # if you change anything in the ring, like the term ordering
    # this needs to be called
    void rUnComplete(ring *r)

    # this enables the changes done after rUnComplete
    void rComplete(ring *r, int force)

    ring *rCopy0(ring *)


    ###

    ### Polynomials
    ####

    ### The rule of thumb is: p_XXX accepts a ring as parameter, while
    ### pXXX doesn't. It relies on currRing.

    ## Constructions / Destructors
    ## -------------------------------

    # new empty monomial
    poly *p_Init(ring *r)

    # const polynomial from int
    poly *p_ISet(int i, ring *r)

    # const polynomial from number (coefficient)
    poly *p_NSet(number *n,ring *r)

    # destructor
    void p_Delete(poly **p, ring *r)

    # set the coeff n for the current list element p in r
    int p_SetCoeff(poly *p, number *n, ring *r)

    # get the coeff of the current list element p in r
    number *p_GetCoeff(poly *p, ring *r)

    # sets the exponent e at index v for the list element (monomial) p in r
    # v starts counting at 1
    int p_SetExp(poly *p, int v, int e, ring *r)

    # get the exponent at index v of the monomial p in r
    int p_GetExp(poly *p, int v, ring *r)

    # if SetExp is called on p, p_Setm needs to be called afterwards
    # to finalize the change.
    void p_Setm(poly *p, ring *r)

    # gets a component out of a polynomial vector
    poly *pTakeOutComp1(poly **, int)

    # copies p
    poly *p_Copy(poly *p, ring *r)

    # homogenizes p by multiplying certain powers of the varnum-th variable
    poly *pHomogen (poly *p, int varnum)

    # returns whether a polynomial is homogenous.
    int pIsHomogeneous(poly *p)

    char *p_String(poly *p, ring *r, ring *r)

    ## Arithmetic
    ## -------------------------------

    # returns -p, p is destroyed
    poly *p_Neg(poly *p, ring *r)

    # returns p*n, p is const (i.e. copied)
    poly *pp_Mult_nn(poly *p, number *n, ring *r)

    # returns p*m, does neither destroy p nor m
    poly *pp_Mult_mm(poly *p, poly *m, ring *r)

    # returns p+q, destroys p and q
    poly *p_Add_q(poly *p, poly *q, ring *r)

    # return p - m*q, destroys p; const: q,m
    poly *p_Minus_mm_Mult_qq(poly *p, poly *m, poly *q, ring *r)

    # returns p + m*q destroys p, const: q, m
    poly *p_Plus_mm_Mult_qq(poly *p, poly *m, poly *q, ring *r)

    # returns p*q, does neither destroy p nor q
    poly *pp_Mult_qq(poly *p, poly *q, ring *r)

    # returns p*q, destroys p and q
    poly *p_Mult_q(poly *p, poly *q, ring *r)

    poly *pDivide(poly *,poly *)

    # returns the i-th power of p; p will be destroyed, requires global ring
    poly *pPower(poly *p, int exp)

    # returns newly allocated copy of Lm(p), coef is copied,
    # next=NULL, p might be NULL
    poly *p_Head(poly *p, ring *r)

    # returns TRUE, if leading monom of a divides leading monom of b
    # i.e., if there exists a expvector c > 0, s.t. b = a + c;
    int p_DivisibleBy(poly *a, poly *b, ring *r)

    # like pDivisibleBy, except that it is assumed that a!=NULL, b!=NULL
    int p_LmDivisibleBy(poly *a, poly *b, ring *r)

    # least common multiplies for MONOMIALS only, result is written to m
    # p_Setm must be called on m afterwards.
    void pLcm(poly *a, poly *b, poly *m)

    # total degree of p
    long pTotaldegree(poly *p, ring *r)

    # iterates through the monomials of p
    poly *pNext(poly *p)

    int p_Cmp(poly *l, poly *r, ring *r)
    int p_ExpVectorEqual(poly *p, poly *m, ring *r)

    int p_IsConstant(poly *, ring *)
    int p_LmIsConstant(poly *p, ring *)
    int p_IsUnit(poly *, ring *)

    poly *pSubst(poly *, int varidx, poly *value)

    poly *pInvers(int n, poly *, intvec *)

    # gcd of f and g
    poly *singclap_gcd ( poly *f, poly *g )

    # resultants of f and g in x
    poly *singclap_resultant ( poly *f, poly *g , poly *x)

    # extended GVD
    int singclap_extgcd( poly *f, poly *g, poly *res, poly *pa, poly *pb )

    # polynomial division (as opposed to monomial division)
    poly *singclap_pdivide ( poly *f, poly *g )

    # factorization
    ideal *singclap_factorize ( poly *f, intvec ** v , int with_exps)

    # is square free
    int singclap_isSqrFree(poly *f)

    # normal form calculation of p with respect to F, Q is quotient
    # ring.
    poly *kNF(ideal *F, ideal *Q, poly *p)

    # General Numbers
    poly *pDiff(poly *, int)

    number *n_Init(int n, ring *r)
    void n_Delete(number **n, ring *r)
    int nInt(number *n)
    number *n_Div(number *a, number *b, ring *r)
    int n_GreaterZero(number *a, ring *r)
    int n_IsZero(number *a, ring *r)

    number *n_Sub(number *a, number *b, ring *r)
    number *nInvers(number *n)

    # Rational Numbers
    number *nlInit(int)
    number *nlInit2gmp(mpz_t i, mpz_t j)
    number *nlInit2(int i, int j)
    number *nlGetNom(number *n, ring *r)
    number *nlGetDenom(number *n, ring *r)
    number *nlRInit(int)


    # Algebraic Numbers
    number *naPar(int)
    void naPower(number *, int, number **)
    number *naMult(number *, number *)
    number *naAdd(number *, number *)
    number *naCopy(number *)
    number *naInit(int)
    void naDelete(number **, ring*)
    int naIsZero(number *)
    char * naRead(char *s, number *)
    int naIsOne(number *)
    int naIsZero(number *)

    number *napGetCoeff(napoly *z)
    int napGetExp(napoly *, int)
    napoly *napIter(napoly *)

    # Integer Numbers
    cdef long SR_INT
    long SR_TO_INT(number *)
    long SR_HDL(mpz_t )


    ctypedef struct napoly "polyrec"

    ctypedef struct lnumber "slnumber":
        napoly *z
        napoly *n
        int s


    # Ideals
    ideal *idInit(int size, int rank)
    void id_Delete(ideal **, ring *)
    ideal *fast_map(ideal *, ring *, ideal *, ring *)
    ideal *idLift(ideal *mod, ideal *submod, ideal **rest, int goodShape, int isSB, int divide)
    void idShow(ideal *i)
    int IDELEMS(ideal *i)


    void idSkipZeroes (ideal *ide)
    long idRankFreeModule(ideal *m, ring *r)
    ideal *kStd(ideal *i, ideal *q, tHomog h, intvec *w)
    ideal *t_rep_gb(ring *r,ideal *arg_I, int syz_comp, int F4_mode)
