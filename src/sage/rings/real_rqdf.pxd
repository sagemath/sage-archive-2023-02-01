from sage.structure.element cimport RingElement, ModuleElement, Element, FieldElement
from sage.rings.ring cimport Field
from sage.structure.parent  cimport Parent

cdef extern from *:
    bint finite(double)
    bint isinf(double)
    double NAN
    double INFINITY

cdef extern from "stdsage.h":
    object PY_NEW(object t)
    void* PY_TYPE(object o)
    int PY_TYPE_CHECK(object o, object t)
    void PY_SET_TP_NEW(object t1, object t2)

cdef extern from "Python.h":
    ctypedef struct PyTypeObject
    ctypedef struct PyObject
    bint PyObject_TypeCheck(object o, PyTypeObject *t)

cdef extern from "qd/qd_real.h":

    ctypedef struct qd "qd_real":
        # Members
        double *x
        qd _pi
        qd _log2
        qd _nan
        qd _e
        # Methods
        void (* write)(char *s, int prec, int show, int upper)
        void (* to_digits)(char *s, int expn, int prec)
        int (* is_negative)()
    # Init
    qd *qd_from_str "new qd_real"(char *a)
    qd *qd_from_double "new qd_real"(double a)
    qd *new_qd_real "new qd_real"()
    qd *qd_from_qd "new qd_real"(double x0, double x1, double x2, double x3)
    qd *qd_from_int "new qd_real"(int i)

    # Conversions
    int qd_to_int "to_int"(qd a)
    double qd_to_double "to_double"(qd a)

    # Nan, inf
    bint qd_is_nan "isnan"(qd a)
    bint qd_is_inf "isinf"(qd a)

# C interface to qd
cdef extern from "qd/c_qd.h":
    # Arithmetic
    void c_qd_add(double *a,double *b, double *c)
    void c_qd_mul(double *a,double *b, double *c)
    void c_qd_sub(double *a,double *b, double *c)
    void c_qd_div(double *a,double *b, double *c)

    void c_qd_abs(double *a,double *b)
    void c_qd_neg(double *a,double *b)

    # Powers, etc
    void c_qd_sqr(double *a,double *b)
    void c_qd_sqrt(double *a,double *b)
    void c_qd_npwr(double *a, int b, double *c)
    void c_qd_nroot(double *a, int b, double *c)

    # Rounding
    void c_qd_floor(double *a,double *b)
    void c_qd_ceil(double *a,double *b)

    # exp, log
    void c_qd_exp(double *a,double *b)
    void c_qd_log(double *a,double *b)
    void c_qd_log10(double *a,double *b)

    # trigonometric
    void c_qd_sin(double *a,double *b)
    void c_qd_cos(double *a,double *b)
    void c_qd_tan(double *a,double *b)

    void c_qd_asin(double *a,double *b)
    void c_qd_acos(double *a,double *b)
    void c_qd_atan(double *a,double *b)
    void c_qd_atan2(double *a,double *b)

    # hyperbolic
    void c_qd_sinh(double *a,double *b)
    void c_qd_cosh(double *a,double *b)
    void c_qd_tanh(double *a,double *b)

    void c_qd_asinh(double *a,double *b)
    void c_qd_acosh(double *a,double *b)
    void c_qd_atanh(double *a,double *b)

    # comparison
    void c_qd_comp(double *a,double *b, int *result)

    # random number generation
    void c_qd_rand(double *a)

    # SAGE-specific
    void delete "delete "(void *o)
    qd qd_deref "*"(qd *q)

# For controlling the round-to-double bit on x86 machines
# These functions have no effects on other platforms
cdef extern from "qd/fpu.h":
    void fpu_fix_start(unsigned int *old_cw)
    void fpu_fix_end(unsigned int *old_cw)

cdef class RealQuadDoubleField_class(Field):
    # so it is possible to make weakrefs to this finite field
    cdef object __weakref__
    # round-to-double bit
    cdef unsigned int *cwf

cdef class QuadDoubleElement(FieldElement):
    cdef qd *initptr #
    # round-to-double bit
    cdef unsigned int *cw

    cdef _set(self,x)
    cdef _new(self)
    cdef _new_c(self,qd a)
    #cdef abs(QuadDoubleElement s)
