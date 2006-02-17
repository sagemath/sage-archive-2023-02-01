include "interrupt.pxi"

ctypedef unsigned long ulong

cdef extern from "stdlib.h":
    ctypedef unsigned long size_t
    void free(void *ptr)

cdef extern from 'pari/pari.h':
    ctypedef unsigned long pari_sp
    void pari_init(size_t parisize, ulong maxprime)
    ctypedef long* GEN
    pari_sp avma
    extern unsigned long precdl
    GEN gen_1, gen_0
    void output(GEN)
    char*   GENtostr(GEN x)

cdef extern from "ellcurv.h":
    GEN sage_docurve(char *INP,
                      int dorank, int dox0isog, int domoddeg,
                      int doanalrank, int doarderiv)
    ctypedef long long int int64

    int ELLACC
    int ROOTNO
    GEN TAMA    # product of the local Tamagawa number
    GEN COND    # conductor
    GEN SSCOND  # symmetric square conductor
    GEN CURVE   # the curve
    GEN CURVES  # list of curves in the isogeny class
    GEN MODDEG  # the modular degree
    GEN X1VOL   # the volume of period lattice of minimal faltings height curve
                # the "X1" suggests "X1-optimal"
    GEN TWPROD  # the product of the primes you twist by to get to the
                # input curve, except if you have CM

    int PLACE   # position of input curve in the isogeny class isomorphic to input curve, with 1-based indexing
    int ISPRIME # is the conductor prime?  maybe minimal twist conductor ?
    int ISSQFREE
    int ISSETZER # the curve
    int CM       # does curve have CM?  (either 0 or -something, but maybe wrong); nothing uses this
    int ISOG # the maximal degree of an cyclic isogeny between any two curves in the isogeny class, up to sign;  sign for 4 or 9 tells you whether or not some kind of spookie thing is going on (should be in Cremona-Watkins; tells whether 2^5 or 3^3 divides the conductor)
    int AISOG # absolute value
    int NI  # number of isogenies
    int X0NUM # where in the CURVES list the X_0-optimal curve probably is
    int64 PMAX
    int PSIZE
    int TRACE
    int PCOUNT
    int PRINT

cdef object GEN_to_str(GEN g):
    cdef char* c
    c = GENtostr(g)
    s = str(c)
    free(c)
    return s

def ec(ainvs):
    r"""
    Compute information about a given elliptic curve using Mark Watkins's ec C library.

    \note{This only works on 32-bit platforms.  It raises a RuntimeError exception
    on 64-bit platforms.}

    EXAMPLES:
        sage: import sage.libs.ec.all as ec
        sage: ec.ec([1,2,3,4,5])
        {'conductor': '10351', 'Modular degree': '464', 'Curve1': '[1, -1, 0, 4, 3, -3, 8, 12, -25, -183, -3429, -10351, 6128487/10351]', 'minimal equation': '[1, -1, 0, 4, 3, -3, 8, 12, -25, -183, -3429, -10351, 6128487/10351]', 'X0NUM': '1', 'lratio': '1.2653937917645377009', 'analytic rank': '1', 'X0CURVE': '[1, -1, 0, 4, 3, -3, 8, 12, -25, -183, -3429, -10351, 6128487/10351]', 'X0CURVE moddeg': '464'} # 32-bit
        Traceback (most recent call last):                       # 64-bit
        ...                                                      # 64-bit
        RuntimeError: ec does not work on 64-bit platforms       # 64-bit
    """
    if not isinstance(ainvs, list) or len(ainvs) != 5:
        raise TypeError, "ainvs must be a list of length 5"

    import sage.misc.all
    if sage.misc.all.is_64_bit:
        raise RuntimeError, "ec does not work on 64-bit platforms"

    if (not precdl):
        pari_init(10000000, 1000000)

    cdef pari_sp save_avma
    global avma
    save_avma = avma

    s = str(ainvs)
    cdef GEN g
    _sig_on
    g = sage_docurve(s, True, True, True, True, False)
    _sig_off
    st = GEN_to_str(g)
    avma = save_avma
    v = st.split("|")[:-1]
    X = {}
    i = 0
    while i < len(v)-1:
        X[v[i]] = v[i+1]
        i = i + 2
    return X


