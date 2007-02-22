include 'vector_rational_sparse_c.pxi'


cdef class Vector_mpq

cdef void Vector_mpq_rescale(Vector_mpq w, mpq_t x):
    scale_mpq_vector(&w.v, x)

cdef class Vector_mpq:
    """
    Vector_mpq -- a sparse vector of GMP rationals.  This is a Python
    extension type that wraps the C implementation of sparse vectors
    modulo a small prime.
    """
    cdef mpq_vector v

    def __init__(self, int degree, int num_nonzero=0, entries=[], sort=True):
        cdef int i
        init_mpq_vector(&self.v, degree, num_nonzero)
        if entries != []:
            if len(entries) != num_nonzero:
                raise ValueError, "length of entries (=%s) must equal num_nonzero (=%s)"%(len(entries), num_nonzero)
            if sort:
                entries = list(entries) # copy so as not to modify passed in list
                entries.sort()
            for i from 0 <= i < num_nonzero:
                s = str(entries[i][1])
                mpq_set_str(self.v.entries[i], s, 0)
                self.v.positions[i] = entries[i][0]

    def __dealloc__(self):
        clear_mpq_vector(&self.v)

    def __getitem__(self, int n):
        cdef mpq_t x
        cdef sage.rings.rational.Rational a
        mpq_init(x)
        mpq_vector_get_entry(&x, &self.v, n)
        a = sage.rings.rational.Rational()
        a.set_from_mpq(x)
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
            return bool(n < 0)
        elif op == 1:
            return bool(n <= 0)
        elif op == 2:
            return bool(n == 0)
        elif op == 3:
            return bool(n != 0)
        elif op == 4:
            return bool(n > 0)
        elif op == 5:
            return bool(n >= 0)

    def __setitem__(self, int n, x):
        cdef object s
        s = str(x)
        mpq_vector_set_entry_str(&self.v, n, s)

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
        z = self.copy()
        mpq_init(t)
        s = str(other)
        mpq_set_str(t, s, 0)
        Vector_mpq_rescale(z, t)
        mpq_clear(t)
        return z

    def randomize(self, int sparcity, bound=3):
        """
        randomize(self, int sparcity, exact=False):

        The sparcity is a bound on the number of nonzeros per row.
        """
        cdef int i
        for i from 0 <= i < sparcity:
            self[random.randrange(self.v.degree)] = random.randrange(1,bound)

