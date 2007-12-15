## NOT DONE

include 'vector_rational_sparse_c.pxi'


from sage.rings.rational cimport Rational
cimport free_module_element


cdef class Vector_mpq

cdef void Vector_mpq_rescale(Vector_mpq w, mpq_t x):
    mpq_vector_scale(&w.v, x)

cdef class Vector_mpq:
    """
    Vector_mpq -- a sparse vector of GMP rationals.
    """
    cdef mpq_vector v

    # We are doing coersion all the time. Should that be a flag?
    def __init__(self, int degree, int num_nonzero=0, entries=[], sort=True):
        cdef int i
        cdef Rational z
        mpq_vector_init(&self.v, degree, num_nonzero)
        if entries != []:
            if len(entries) != num_nonzero:
                raise ValueError, "length of entries (=%s) must equal num_nonzero (=%s)"%(len(entries), num_nonzero)
            if sort:
                entries = list(entries) # copy so as not to modify passed in list
                entries.sort()
            for i from 0 <= i < num_nonzero:
                z = Rational(entries[i])
                mpq_set(self.v.entries[i], z.value)
                self.v.positions[i] = entries[i][0]

    def __dealloc__(self):
        clear_mpq_vector(&self.v)

    def __getitem__(self, int n):
        cdef mpq_t x
        cdef Rational a
        mpq_init(x)
        mpq_vector_get_entry(&x, &self.v, n)
        a = PY_NEW(Rational)
        mpq_set(a.value, x)
        mpq_clear(x)
        return a

    def cmp(self, Vector_mpq other):
        return mpq_vector_cmp(&self.v, &other.v)

    def __richcmp__(Vector_mpq self, x, int op):
        if not isinstance(x, Vector_mpq):
            return -1
        cdef int n
        n = self.cmp(x)
        if op == 0:
            return (n < 0)
        elif op == 1:
            return (n <= 0)
        elif op == 2:
            return (n == 0)
        elif op == 3:
            return (n != 0)
        elif op == 4:
            return (n > 0)
        else:
            return (n >= 0)

    def __setitem__(self, Py_ssize_t n, x):
        if not self._is_mutable:
            raise ValueError, "vector is immutable; please change a copy instead (use self.copy())"
        mpq_vector_set_entry(&self.v, n, (<Rational> x).value)

    def __repr__(self):
        return str(list(self))

    def degree(self):
        return self.v.degree

    def num_nonzero(self):
        return self.v.num_nonzero

    def list(self):
        return mpq_vector_to_list(&self.v)

    cdef void rescale(self, mpq_t x):
        scale_mpq_vector(&self.v, x)

    def __add__(Vector_mpq self, Vector_mpq other):
        cdef mpq_vector z1, *z2
        cdef Vector_mpq w
        cdef mpq_t ONE

        mpq_init(ONE)
        mpq_set_si(ONE,1,1)

        add_mpq_vector_init(&z1, &self.v, &other.v, ONE)
        mpq_clear(ONE)

        w = Vector_mpq(self.v.degree)
        z2 = &(w.v)
        clear_mpq_vector(z2)   # free memory wasted on allocated w
        z2.entries = z1.entries
        z2.positions = z1.positions
        z2.num_nonzero = z1.num_nonzero
        # at this point we do *not* free z1, since it is referenced by w.
        return w

    def __sub__(Vector_mpq self, Vector_mpq other):
        return self + other*(-1)

    def copy(self):
        cdef int i
        cdef Vector_mpq w
        w = Vector_mpq(self.v.degree, self.v.num_nonzero)
        for i from 0 <= i < self.v.num_nonzero:
            mpq_set(w.v.entries[i], self.v.entries[i])
            w.v.positions[i] = self.v.positions[i]
        return w

    def __mul__(x, y):
        if isinstance(x, Vector_mpq):
            self = x
            other = y
        elif isinstance(y, Vector_mpq):
            self = y
            other = x
        else:
            raise TypeError, "Invalid types."
        cdef object s, z
        cdef mpq_t t
        cdef Rational r

        z = self.copy()
        r = Rational(y)
        Vector_mpq_rescale(z, r.value)
        return z

    def randomize(self, int sparcity, bound=3):
        """
        randomize(self, int sparcity, exact=False):

        The sparcity is a bound on the number of nonzeros per row.
        """
        cdef int i
        for i from 0 <= i < sparcity:
            self[random.randrange(self.v.degree)] = random.randrange(1,bound)

