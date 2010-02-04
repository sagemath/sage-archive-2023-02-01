"""
EXAMPLES:
    sage: from sage.libs.mwrank.mwrank import _Curvedata, _mw
    sage: c = _Curvedata(1,2,3,4,5)

    sage: print c
    [1,2,3,4,5]
    b2 = 9       b4 = 11         b6 = 29         b8 = 35
    c4 = -183           c6 = -3429
    disc = -10351       (# real components = 1)
    #torsion not yet computed

    sage: t= _mw(c)
    sage: t.search(10)
    sage: t
    [[1:2:1]]
"""

import os

include '../../ext/interrupt.pxi'

cdef extern from "stdlib.h":
    void free(void *ptr)

cdef extern from "wrap.h":
    ### misc functions ###
    void mwrank_set_precision(long n)
    void mwrank_initprimes(char* pfilename, int verb)


    ### bigint ###
    struct bigint
    bigint* new_bigint()
    void del_bigint(bigint* x)
    bigint* str_to_bigint(char* s)
    char* bigint_to_str(bigint* x)

    ### Curvedata ###
    struct Curvedata
    Curvedata* Curvedata_new(bigint* a1, bigint* a2,
                             bigint* a3, bigint* a4,
                             bigint* a6, int min_on_init)
    void Curvedata_del(Curvedata* curve)
    char* Curvedata_repr(Curvedata* curve)
    double Curvedata_silverman_bound(Curvedata* curve)
    double Curvedata_cps_bound(Curvedata* curve)
    double Curvedata_height_constant(Curvedata* curve)
    char* Curvedata_getdiscr(Curvedata* curve)
    char* Curvedata_conductor(Curvedata* m)
    char* Curvedata_isogeny_class(Curvedata* E, int verbose)

    ## mw ##
    struct mw
    mw* mw_new(Curvedata* curve, int verb, int pp, int maxr)
    void mw_del(mw* m)
    int mw_process(Curvedata* curve, mw* m,
                   bigint* x, bigint* y,
                   bigint* z, int sat)
    char* mw_getbasis(mw* m)
    char* mw_regulator(mw* m)
    int mw_rank(mw* m)
    int mw_saturate(mw* m, bigint* index, char** unsat,
                    long sat_bd, int odd_primes_only)
    void mw_search(mw* m, char* h_lim, int moduli_option, int verb)

    ### two_descent ###
    struct two_descent
    two_descent* two_descent_new(Curvedata* curve,
                                 int verb, int sel,
                                 long firstlim, long secondlim,
                                 long n_aux, int second_descent)

    void two_descent_del(two_descent* t)
    int two_descent_ok(two_descent* t)
    long two_descent_get_certain(two_descent* t)
    char* two_descent_get_basis(two_descent* t)
    char* two_descent_regulator(two_descent* t)
    long two_descent_get_rank(two_descent* t)
    long two_descent_get_rank_bound(two_descent* t)
    long two_descent_get_selmer_rank(two_descent* t)
    void two_descent_saturate(two_descent* t, long sat_bd)



cdef object string_sigoff(char* s):
    _sig_off
    # Makes a python string and deletes what is pointed to by s.
    t = str(s)
    free(s)
    return t

class __init:
    pass

_INIT = __init()
# set the default
mwrank_set_precision(50)

def set_precision(n):
    mwrank_set_precision(n)

def initprimes(filename, verb=False):
    """
    mwrank_initprimes(filename, verb=False):

    INPUT:
        filename -- (string) the name of a file of primes
        verb -- (bool: default False) verbose or not?

    EXAMPLES:
        sage: file= SAGE_TMP + '/PRIMES'
        sage: open(file,'w').write(' '.join([str(p) for p in prime_range(10^6)]))
        sage: mwrank_initprimes(file, verb=False)
        sage: mwrank_initprimes("x" + file, True)
        Traceback (most recent call last):
        ...
        IOError: No such file or directory: ...
    """
    if not os.path.exists(filename):
        raise IOError, 'No such file or directory: %s'%filename
    mwrank_initprimes(filename, verb)

############# bigint #################

cdef class _bigint:
    cdef bigint* x

    def __init__(self, x="0"):
        if not (x is _INIT):
            s = str(x)
            if s.isdigit() or s[0] == "-" and s[1:].isdigit():
                self.x = str_to_bigint(s)
            else:
                raise ValueError, "invalid _bigint: %s"%x

    def __dealloc__(self):
        del_bigint(self.x)

    def __repr__(self):
        _sig_on
        return string_sigoff(bigint_to_str(self.x))


cdef make_bigint(bigint* x):
    cdef _bigint y
    _sig_off
    y = _bigint(_INIT)
    y.x = x
    return y

############# Curvedata #################

cdef class _Curvedata:
    cdef Curvedata* x

    def __init__(self, a1, a2, a3,
                 a4, a6, min_on_init=0):
        cdef _bigint _a1, _a2, _a3, _a4, _a6
        _a1 = _bigint(a1)
        _a2 = _bigint(a2)
        _a3 = _bigint(a3)
        _a4 = _bigint(a4)
        _a6 = _bigint(a6)
        self.x = Curvedata_new(_a1.x, _a2.x, _a3.x, _a4.x, _a6.x, min_on_init)
        if self.discriminant() == 0:
            raise ArithmeticError, "Invariants (= %s) do not describe an elliptic curve."%([a1,a2,a3,a4,a6])

    def __dealloc__(self):
        Curvedata_del(self.x)

    def __repr__(self):
        _sig_on
        return string_sigoff(Curvedata_repr(self.x))[:-1]

    def silverman_bound(self):
        return Curvedata_silverman_bound(self.x)

    def cps_bound(self):
        cdef double x
        # We declare x so there are *no* Python library
        # calls within the _sig_on/_sig_off.
        _sig_on
        x = Curvedata_cps_bound(self.x)
        _sig_off
        return x

    def height_constant(self):
        return Curvedata_height_constant(self.x)

    def discriminant(self):
        _sig_on
        return int(string_sigoff(Curvedata_getdiscr(self.x)))

    def conductor(self):
        _sig_on
        return int(string_sigoff(Curvedata_conductor(self.x)))

    def isogeny_class(self, verbose=False):
        _sig_on
        s = string_sigoff(Curvedata_isogeny_class(self.x, verbose))
        _sig_off
        return eval(s)


############# _mw #################

cdef class _mw:
    cdef mw* x
    cdef Curvedata* curve

    def __init__(self, _Curvedata curve, verb=False, pp=1, maxr=999):
        self.curve = curve.x
        self.x = mw_new(curve.x, verb, pp, maxr)

    def __dealloc__(self):
        mw_del(self.x)

    def __repr__(self):
        _sig_on
        return string_sigoff(mw_getbasis(self.x))


    def process(self, point, sat=True):
        if not isinstance(point, (tuple, list)) and len(point) == 3:
            raise TypeError, "point must be a list or tuple of length 3."
        cdef _bigint x,y,z
        _sig_on
        x,y,z = _bigint(point[0]), _bigint(point[1]), _bigint(point[2])
        r = mw_process(self.curve, self.x, x.x, y.x, z.x, sat)
        _sig_off
        if r != 0:
            raise ArithmeticError, "point (=%s) not on curve."%point

    def getbasis(self):
        _sig_on
        s = string_sigoff(mw_getbasis(self.x))
        _sig_off
        return s

    def regulator(self):
        cdef float f
        _sig_on
        f = float(string_sigoff(mw_regulator(self.x)))
        _sig_off
        return f

    def rank(self):
        _sig_on
        r = mw_rank(self.x)
        _sig_off
        return r

    def saturate(self, int sat_bd=-1, int odd_primes_only=0):
        cdef _bigint index
        cdef char* s
        cdef int ok
        _sig_on
        index = _bigint()
        ok = mw_saturate(self.x, index.x, &s, sat_bd, odd_primes_only)
        unsat = string_sigoff(s)   # includes _sig_off
        return ok, index, unsat

    def search(self, h_lim, int moduli_option=0, int verb=0):
        cdef char* _h_lim

        h_lim = str(h_lim)
        _h_lim = h_lim

        _sig_on
        mw_search(self.x, _h_lim, moduli_option, verb)
        _sig_off


############# two_descent #################
cdef class _two_descent:
    cdef two_descent* x

    def __init__(self):
        self.x = <two_descent*> 0

    def __dealloc__(self):
        if self.x:
            two_descent_del(self.x)

    def do_descent(self, _Curvedata curve,
                 int verb = 1,
                 int sel = 0,
                 int firstlim = 20,
                 int secondlim = 8,
                 int n_aux = -1,
                 int second_descent = 1):
        _sig_on
        self.x = two_descent_new(curve.x, verb, sel, firstlim, secondlim, n_aux, second_descent)
        _sig_off

    def getrank(self):
        cdef int r
        _sig_on
        r = two_descent_get_rank(self.x)
        _sig_off
        return r

    def getrankbound(self):
        cdef int r
        _sig_on
        r = two_descent_get_rank_bound(self.x)
        _sig_off
        return r

    def getselmer(self):
        _sig_on
        r = two_descent_get_selmer_rank(self.x)
        _sig_off
        return r

    def ok(self):
        return two_descent_ok(self.x)

    def getcertain(self):
        return two_descent_get_certain(self.x)

    def saturate(self, saturation_bound=0):
        _sig_on
        two_descent_saturate(self.x, saturation_bound)
        _sig_off

    def getbasis(self):
        _sig_on
        return string_sigoff(two_descent_get_basis(self.x))

    def regulator(self):
        _sig_on
        return float(string_sigoff(two_descent_regulator(self.x)))
