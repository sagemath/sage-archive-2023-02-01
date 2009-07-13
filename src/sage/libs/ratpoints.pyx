r"""
Hyperelliptic Curve Point Finding, via ratpoints.

"""

include '../ext/cdefs.pxi'
include '../ext/stdsage.pxi'
from sage.rings.integer cimport Integer
cdef extern from "ratpoints.h":
    long RATPOINTS_MAX_DEGREE
    long RATPOINTS_ARRAY_SIZE
    long RATPOINTS_DEFAULT_SP1
    long RATPOINTS_DEFAULT_SP2
    long RATPOINTS_DEFAULT_NUM_PRIMES
    long RATPOINTS_DEFAULT_MAX_FORBIDDEN
    long RATPOINTS_DEFAULT_STURM
    long RATPOINTS_NON_SQUAREFREE
    long RATPOINTS_BAD_ARGS

    # for args flags:
    long RATPOINTS_NO_CHECK # when set, do not check whether the surviving
                            # x-coordinates give rise to rational points
    long RATPOINTS_NO_Y # when set, only list x coordinates instead of actual points
    long RATPOINTS_NO_REVERSE # when set, do not modify the mpz_t array
    long RATPOINTS_NO_JACOBI # when set, prevent use of Jacobi symbol test
    long RATPOINTS_VERBOSE # when set, print some output on what ratpoints is doing
    # define RATPOINTS_FLAGS_INPUT_MASK \
    # (RATPOINTS_NO_CHECK | RATPOINTS_NO_Y | RATPOINTS_NO_REVERSE | \
    #  RATPOINTS_NO_JACOBI | RATPOINTS_VERBOSE)


    ctypedef struct ratpoints_interval:
        double low
        double up
    ctypedef struct ratpoints_args:
        mpz_t *cof
        long degree
        long height
        ratpoints_interval *domain
        long num_inter
        long b_low
        long b_high
        long sp1
        long sp2
        long array_size
        long sturm
        long num_primes
        long max_forbidden
        unsigned int flags
        # from here: private data
        # mpz_t *work
        # void *se_buffer
        # void *se_next
        # void *ba_buffer
        # void *ba_next
        # int *int_buffer
        # int *int_next
        # void *sieve_list
    long find_points(ratpoints_args*, int proc(long, long, mpz_t, void*, int*), void*)
    void find_points_init(ratpoints_args*)
    long find_points_work(ratpoints_args*, int proc(long, long, mpz_t, void*, int*), void*)
    void find_points_clear(ratpoints_args*)

ctypedef struct point_list:
    long *xes
    mpz_t *ys
    long *zs
    long array_size
    long num_points
    long max_num_points

cdef int process(long x, long z, mpz_t y, void *info0, int *quit):
    # ratpoints calls this function when it finds a point [x : y : z]
    # info0 is the pointer passed to ratpoints originally
    # if quit[0] is set to a nonzero value, ratpoints will abort immediately
    cdef point_list *plist = <point_list *> info0
    cdef long i
    if plist.array_size == plist.num_points:
        i = plist.array_size
        plist.array_size *= 2
        plist.xes = <long *> sage_realloc(plist.xes, plist.array_size * sizeof(long))
        plist.ys = <mpz_t *> sage_realloc(plist.ys, plist.array_size * sizeof(mpz_t))
        plist.zs = <long *> sage_realloc(plist.zs, plist.array_size * sizeof(long))
        while i < plist.array_size:
            mpz_init(plist.ys[i])
            i += 1
    plist.xes[plist.num_points] = x
    mpz_set(plist.ys[plist.num_points], y)
    plist.zs[plist.num_points] = z
    plist.num_points += 1
    if plist.max_num_points > 0:
        if plist.max_num_points == plist.num_points:
            quit[0] = -1
    return 1 # weight for counting the points

def ratpoints(list coeffs, long H, verbose=False, long max=0):
    """
    Access the ratpoints library to find points on the hyperelliptic curve:

    `y^2 = a_n x^n + \cdots + a_1 x + a_0.`

    INPUT::

        coeffs  -- list of integer coefficients a_0, a_1, ..., a_n

        H       -- the bound for the denominator and the absolute value of the
                   numerator of the x-coordinate

        verbose -- if True, ratpoints will print comments about its progress

        max     -- maximum number of points to find (if 0, find all of them)

    OUTPUT::

    The points output by this program are points in (1, ceil(n/2), 1)-weighted
    projective space. If n is even, then the associated homogeneous equation is
    `y^2 = a_n x^n + \cdots + a_1 x z^{n-1} + a_0 z^n` while if n is odd, it is
    `y^2 = a_n x^n z + \cdots + a_1 x z^n + a_0 z^{n+1}`.

    EXAMPLE::

        sage: from sage.libs.ratpoints import ratpoints
        sage: for x,y,z in ratpoints([1..6], 200):
        ...    print -1*y^2 + 1*z^6 + 2*x*z^5 + 3*x^2*z^4 + 4*x^3*z^3 + 5*x^4*z^2 + 6*x^5*z
        0
        0
        0
        0
        0
        0
        0
        sage: for x,y,z in ratpoints([1..5], 200):
        ...    print -1*y^2 + 1*z^4 + 2*x*z^3 + 3*x^2*z^2 + 4*x^3*z + 5*x^4
        0
        0
        0
        0
        0
        0
        0
        0

    """
    cdef ratpoints_args args
    cdef long i, total, verby
    cdef Integer sage_int, s_x, s_y, s_z
    cdef point_list *plist

    coeffs = [Integer(a) for a in coeffs]
    assert len(coeffs)-1 <= RATPOINTS_MAX_DEGREE

    verby = ~0 if verbose else 0

    args.degree = len(coeffs)-1

    args.cof = <mpz_t *> sage_malloc((args.degree+1) * sizeof(mpz_t))
                                          # example uses RATPOINTS_MAX_DEGREE -- necessary?
    args.domain = <ratpoints_interval *> sage_malloc(2*args.degree * sizeof(ratpoints_interval))
    plist = <point_list *> sage_malloc(sizeof(point_list))
    if max == 0:
        plist.array_size = 64
    else:
        plist.array_size = max
    plist.xes = <long *> sage_malloc(plist.array_size * sizeof(long))
    plist.ys = <mpz_t *> sage_malloc(plist.array_size * sizeof(mpz_t))
    for i from 0 <= i < plist.array_size:
        mpz_init(plist.ys[i])
    plist.zs = <long *> sage_malloc(plist.array_size * sizeof(long))
    plist.num_points = 0
    plist.max_num_points = max

    args.height = H
    args.num_inter = 0
    args.b_low = 1
    args.b_high = H
    args.sp1 = RATPOINTS_DEFAULT_SP1
    args.sp2 = RATPOINTS_DEFAULT_SP2
    args.array_size = RATPOINTS_ARRAY_SIZE
    args.sturm = RATPOINTS_DEFAULT_STURM
    args.num_primes = RATPOINTS_DEFAULT_NUM_PRIMES
    args.max_forbidden = RATPOINTS_DEFAULT_MAX_FORBIDDEN
    args.flags = (RATPOINTS_VERBOSE & verby)

    for i from 0 <= i <= args.degree:
        mpz_init(args.cof[i])
        sage_int = <Integer> coeffs[i]
        mpz_set(args.cof[i], sage_int.value)

    total = find_points(&args, process, <void *>plist)
    if total == RATPOINTS_NON_SQUAREFREE:
        raise RuntimeError('Polynomial must be square-free')
    if total == RATPOINTS_BAD_ARGS:
        raise RuntimeError('Bad arguments to ratpoints')

    for i from 0 <= i <= args.degree:
        mpz_clear(args.cof[i])

    sage_free(args.cof)
    sage_free(args.domain)

    cdef list L = []
    for i from 0 <= i < plist.num_points:
        s_x = Integer(0)
        s_y = Integer(0)
        s_z = Integer(0)
        mpz_set_si(s_x.value, plist.xes[i])
        mpz_set(s_y.value, plist.ys[i])
        mpz_set_si(s_z.value, plist.zs[i])
        L.append((s_x,s_y,s_z))

    for i from 0 <= i < plist.array_size:
        mpz_clear(plist.ys[i])
    sage_free(plist.xes)
    sage_free(plist.ys)
    sage_free(plist.zs)
    sage_free(plist)

    return L


