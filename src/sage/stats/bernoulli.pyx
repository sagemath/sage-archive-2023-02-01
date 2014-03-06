from sage.rings.real_mpfr cimport RealNumber, RealField
from sage.rings.integer cimport Integer
from sage.misc.randstate cimport randstate, current_randstate

from dgs cimport dgs_bern_exp_mp_init, dgs_bern_exp_mp_call, dgs_bern_exp_mp_clear
from dgs cimport dgs_bern_mp_init, dgs_bern_mp_call, dgs_bern_mp_clear
from dgs cimport DGS_BERN_UNIFORM_MAX_LENGTH, dgs_bern_uniform_init, dgs_bern_uniform_call, dgs_bern_uniform_clear

cdef class BernoulliBase(SageObject):
    def __init__(self):
        self.b = tuple([Integer(0), Integer(1)])

    def integer(self, *args, **kwds):
        return self.b[self.raw_call(*args, **kwds)]

    def boolean(self, *args, **kwds):
        return bool(self.raw_call(*args, **kwds))

    def __call__(self, *args, **kwds):
        return self.b[self.raw_call(*args, **kwds)]

cdef class BernoulliUniformSampler(BernoulliBase):
    def __init__(self, sample_length=0):
        if sample_length < 0 or sample_length > DGS_BERN_UNIFORM_MAX_LENGTH:
            raise ValueError("sample_length must be larger than zero and smaller than or equal to %d"%DGS_BERN_UNIFORM_MAX_LENGTH)

        BernoulliBase.__init__(self)
        self._gen = dgs_bern_uniform_init(sample_length)

    def __clear__(self):
        if self._gen:
            dgs_bern_uniform_clear(self._gen)

    cpdef int raw_call(self):
        cdef randstate rstate = current_randstate()
        return dgs_bern_uniform_call(self._gen, rstate.gmp_state)

cdef class BernoulliSampler(BernoulliBase):
    def __init__(self, c):
        if c <= 0 or c >= 1:
            raise ValueError("Parameter c must be between 0 and 1 but got %f"%c)

        if not isinstance(c, RealNumber):
            RR = RealField()
            c = RR(c)

        BernoulliBase.__init__(self)
        self._gen = dgs_bern_mp_init((<RealNumber>c).value)

    cpdef int raw_call(self):
        cdef randstate rstate = current_randstate()
        return dgs_bern_mp_call(self._gen, rstate.gmp_state)
        
    def __clear__(self):
        if self._gen:
            dgs_bern_mp_clear(self._gen)

cdef class BernoulliExpSampler(BernoulliBase):
    def __init__(self, f, l=None):
        if l and l < 0:
            raise ValueError("Parameter l must be >=0 but got %d"%l)

        if f <= 0.0:
            raise ValueError("Parameter f must be > 0 but got %f"%f)
            
        if not isinstance(f, RealNumber):
            RR = RealField()
            c = RR(f)

        BernoulliBase.__init__(self)
        if l is None:
            l = 0
        self._gen = dgs_bern_exp_mp_init((<RealNumber>f).value, l)

    cpdef int raw_call(self, x):
        cdef randstate rstate = current_randstate()
        if not isinstance(x, Integer):
            x = Integer(x)
        return dgs_bern_exp_mp_call(self._gen, (<Integer>x).value, rstate.gmp_state)
        
    def __clear__(self):
        if self._gen:
            dgs_bern_exp_mp_clear(self._gen)

