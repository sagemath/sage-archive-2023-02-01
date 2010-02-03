r"""
Hyperelliptic Curve Point Finding, via ratpoints.

"""

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

        sage: for x,y,z in ratpoints([1..200], 1000):
        ...    print x,y,z
        1 0 0
        0 1 1
        0 -1 1
        201 25353012004564588029934064107520000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000 200
        201 -25353012004564588029934064107520000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000 200


    """
    cdef ratpoints_args args
    cdef long i, total, verby
    cdef Integer sage_int, s_x, s_y, s_z
    cdef point_list *plist

    coeffs = [Integer(a) for a in coeffs]

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

    _sig_on
    total = find_points(&args, process, <void *>plist)
    _sig_off
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

cdef int process_exists_only(long x, long z, mpz_t y, void *info0, int *quit):
    cdef info_struct_exists_only *info_s = <info_struct_exists_only *>info0
    cdef Integer YY
    if info_s.verbose:
        YY = Integer(0); mpz_set(YY.value, y)
        print 'Found point [ %d : %d : %d ], quitting'%(x,YY,z)
    quit[0] = -1
    return 1

cdef int ratpoints_mpz_exists_only(mpz_t *coeffs, long H, int degree, bint verbose) except -1:
    """
    Direct call to ratpoints to search for existence only.

    WARNING - The coefficient array will be modified by ratpoints.
    """
    cdef ratpoints_args args
    cdef info_struct_exists_only info_s
    cdef long total, verby = ~0 if verbose else 0
    info_s.verbose = verbose
    assert degree <= RATPOINTS_MAX_DEGREE
    args.degree = degree
    args.cof = coeffs
    args.domain = <ratpoints_interval *> sage_malloc(2*args.degree * sizeof(ratpoints_interval))
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
    _sig_on
    total = find_points(&args, process_exists_only, <void *>(&info_s))
    _sig_off
    sage_free(args.domain)
    if total == RATPOINTS_NON_SQUAREFREE:
        raise RuntimeError('Polynomial must be square-free')
    if total == RATPOINTS_BAD_ARGS:
        raise RuntimeError('Bad arguments to ratpoints')
    return 1 if (total > 0) else 0




