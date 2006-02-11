"""nodoctest
"""
include 'interrupt.pxi'
include 'decl.pxi'


cdef GEN toGEN(object x) except NULL:
    cdef gen.gen _x
    if isinstance(x, gen.gen):
        _x = x
        return _x.g
    s = str(x)
    cdef GEN g
    _sig_on
    g = flisseq(s)
    _sig_off
    return g

# temp variables
cdef GEN t0,t1,t2,t3,t4

cdef t0GEN(x):
    global t0
    t0 = toGEN(x)
cdef t1GEN(x):
    global t1
    t1 = toGEN(x)
cdef t2GEN(x):
    global t2
    t2 = toGEN(x)
cdef t3GEN(x):
    global t3
    t3 = toGEN(x)
cdef t4GEN(x):
    global t4
    t4 = toGEN(x)


import sage.libs._pari.gen_py

cdef class testclass(gen.gen):
    cdef __init__(self):
        gen.gen.__init__(self)

    def _init(self, g):
        t0GEN(g)
        cdef pari_sp address
        cdef GEN f
        f = self._deepcopy_to_python_heap(t0, &address)
        self.init(f, address)

    def python(self, unsigned long precision=0):
        """
        Return Python eval of self.
        """
        global prec
        if precision:
            s = str(prec)
            p = str(precision)
            _sig_on
            sd_realprecision(p, 2)
            _sig_off
            x = sage.libs._pari.gen_py.python(self)
            _sig_on
            sd_realprecision(s, 2)
            _sig_off
            return x
        else:
            return sage.libs._pari.gen_py.python(self)




    def polcyclo(self, long n, v=-1):
        cdef long var
        var = self.get_var(v)
        _sig_on
        return self.new_gen(cyclo(n, var))


    def elllseries2(self, s, A=1, long prec=0):
        t0GEN(s); t1GEN(A)
        _sig_on
        return self.new_gen(lseriesell(self.g, t0, t1, prec))
