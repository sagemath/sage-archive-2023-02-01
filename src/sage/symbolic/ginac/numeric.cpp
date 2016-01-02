/** @file numeric.cpp
 *
 *  This file contains the interface to the underlying numerics packages.
 *  Its most important design principle is to completely hide the inner
 *  working of that other package from the user of GiNaC.  It must either 
 *  provide implementation of arithmetic operators and numerical evaluation
 *  of special functions or implement the interface to the numericsd package. */

/*
 *  This is a modified version of code included with Ginac.  
 *  The modifications and modified version is:
 * 
 *      GiNaC-SAGE Copyright (C) 2008 William Stein
 *      Copyright (C) 2009 Burcin Erocal <burcin@erocal.org>
 *                (C) 2015 Ralf Stephan <ralf@ark.in-berlin.de>
 *
 *  This program is free software; you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation; either version 2 of the License, or
 *  (at your option) any later version.
 *
 *  This program is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with this program; if not, write to the Free Software
 *  Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301  USA
 */


/*  The original copyright:
 * 
 *  GiNaC Copyright (C) 1999-2008 Johannes Gutenberg University Mainz, Germany
 *
 *  This program is free software; you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation; either version 2 of the License, or
 *  (at your option) any later version.
 *
 *  This program is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with this program; if not, write to the Free Software
 *  Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301  USA
 */

#include "numeric.h"
#include "operators.h"
#include "power.h"
#include "archive.h"
#include "tostring.h"
#include "utils.h"

//#define Logging_refctr
#if defined(Logging_refctr)
#undef Py_INCREF
#undef Py_DECREF

#define Py_INCREF(op) (                         \
    _Py_INC_REFTOTAL  _Py_REF_DEBUG_COMMA       \
    ((PyObject*)(op))->ob_refcnt++) ; \
std::cerr << "+ " << long(op) << ", " << Py_REFCNT(op) << ", " << Py_TYPE(op)->tp_name << std::endl; std::cerr.flush();

#define Py_DECREF(op)                                   \
    do {                                                \
std::cerr << "- " << long(op) << ", " << Py_REFCNT(op) << ", " << Py_TYPE(op)->tp_name << std::endl; std::cerr.flush(); \
        if (_Py_DEC_REFTOTAL  _Py_REF_DEBUG_COMMA       \
        --((PyObject*)(op))->ob_refcnt != 0)            \
            _Py_CHECK_REFCNT(op)                        \
        else                                            \
        _Py_Dealloc((PyObject *)(op));                  \
    } while (0); std::cerr << "Done" << std::endl; std::cerr.flush();
#endif

const GiNaC::numeric& to_numeric(GiNaC::ex& e)
{
        return GiNaC::ex_to<GiNaC::numeric>(e);
}

// Call the Python function f on *this as input and return the result
// as a PyObject*.
#define PY_RETURN(f)  PyObject *a = to_pyobject();		 \
  PyObject *ans = f(a);						 \
  Py_DECREF(a);                                                  \
  if (!ans) py_error("error calling function");			 \
  return ans; 

// Call the Python function f on *this and return the result
// as a PyObject*.
#define PY_RETURN2(f, b)  PyObject *aa = to_pyobject();	 \
  PyObject* bb = b.to_pyobject();				 \
  PyObject *ans = f(aa, bb);					 \
  if (!ans) py_error("error calling function");			 \
  Py_DECREF(aa); Py_DECREF(bb); return ans; 

// Call the Python function f on *this and b, and get back
// a 2-tuple (z,w).  Set c = w, where c should be a 
// reference in the caller, and return z.  This is used
// to return two inputs from a Python function call.  See
// its usage in code below. 
#define PY_RETURN3(f, b, c)						\
  PyObject *aa = to_pyobject();					\
  PyObject* bb = b.to_pyobject();					\
  PyObject *ans = f(aa, bb);						\
  if (!ans) py_error("error calling function");				\
  if (!PyTuple_CheckExact(ans) || PyTuple_GET_SIZE(ans) != 2) 		\
    py_error("error calling function -- return not a 2-tuple.");	\
  PyObject* z =  PyTuple_GET_ITEM(ans, 0); Py_INCREF(z);                \
  PyObject* w =  PyTuple_GET_ITEM(ans, 1); Py_INCREF(w);          \
  c = w;                                                         \
  Py_DECREF(aa); Py_DECREF(bb); Py_DECREF(ans);                   \
  return z; 

//#define DEBUG
//#define VERBOSE

#ifdef DEBUG
#define todo(s) std::cerr << "TODO: " << s << std::endl;
#define stub(s) { std::cerr << "Hit STUB: " << s << std::endl; throw std::runtime_error("stub"); }
#define fake(s) std::cerr << "fake: " << s << std::endl;
#define ASSERT(s, msg) if (!s) { std::cerr << "Failed assertion: " << msg << std::endl; }
#else
#define todo(s)
#define stub(s) { std::cerr << "** Hit STUB**: " << s << std::endl; throw std::runtime_error("stub"); }
#define fake(s)
#endif

#ifdef VERBOSE
#define verbose(s) std::cerr << s << std::endl;
#define verbose2(s,t) std::cerr << s << " " << t << std::endl;
#define verbose3(s,t,u) std::cerr << s << " " << t << ", " << u << std::endl;
#else
#define verbose(s)
#define verbose2(s,t)
#define verbose3(s,t,u)
#endif


//////////////////////////////////////////////////////////////
// Python Interface
//////////////////////////////////////////////////////////////

void py_error(const char*) {
        if (PyErr_Occurred()) {
                throw std::runtime_error("");
        }
}

#if PY_MAJOR_VERSION < 3
#define PyNumber_TrueDivide PyNumber_Divide

inline int Pynac_PyObj_Cmp(PyObject *optr1, PyObject *optr2, const char *errmsg) {
        int result;
        if (PyObject_Cmp(optr1, optr2, &result) == -1)
                py_error(errmsg);
        return result;
}

inline bool Pynac_PyObj_RichCmp(PyObject *optr1, PyObject *optr2, int opid, const char *errmsg) {
        int result;
        if (PyObject_Cmp(optr1, optr2, &result) == -1)
                py_error(errmsg);
        switch (result) {
                case -1: return opid == Py_LT || opid == Py_LE || opid == Py_NE;
                case 0: return opid == Py_LE || opid == Py_EQ || opid == Py_GE;
                case 1: return opid == Py_GT || opid == Py_GE || opid == Py_NE;
                default: return false;
        }
}
#else
#define PyInt_Check PyLong_Check
#define PyInt_AsLong PyLong_AsLong
#define PyInt_FromLong PyLong_FromLong

inline int Pynac_PyObj_Cmp(PyObject *optr1, PyObject *optr2, const char *errmsg) {
        int result = PyObject_RichCompareBool(optr1, optr2, Py_LT);
        if (result == 1)
                return -1;
        else if (result == -1)
                py_error(errmsg);
        else { // result == 0
                result = PyObject_RichCompareBool(optr1, optr2, Py_GT);
                if (result == -1)
                        py_error(errmsg);
        }
        return result;
}

inline bool Pynac_PyObj_RichCmp(PyObject *optr1, PyObject *optr2, int opid, const char *errmsg) {
        int result = PyObject_RichCompareBool(optr1, optr2, opid);
        if (result == -1)
                py_error(errmsg);
        return result;
}
#endif

// The following variable gets changed to true once
// this library has been imported by the Python
// interpreter.  This is because the interpreter calls
// ginac_pyinit_I below, which sets this to true.
// Once this is done we can call all the py_* functions
// defined in Cython modules, which is often much faster
// than what we can do without those. 
static bool initialized = false;

static PyObject* pyfunc_Float = nullptr;

void ginac_pyinit_Float(PyObject* f) {
        Py_INCREF(f);
        pyfunc_Float = f;
}

void ginac_pyinit_I(PyObject* z) {
        initialized = true;
        Py_INCREF(z);
        GiNaC::I = z; // I is a global constant defined below.

        PyObject* m = PyImport_ImportModule("sage.rings.real_mpfr");
        if (!m)
                py_error("Error importing sage.rings.real_mpfr");
        PyObject* obj = PyObject_GetAttrString(m, "RR");
        if (!obj)
                py_error("Error getting RR attribute");
        Py_INCREF(obj);
        GiNaC::RR = obj;
        m = PyImport_ImportModule("sage.rings.complex_field");
        if (!m)
                py_error("Error importing sage.complex_field");
        obj = PyObject_GetAttrString(m, "ComplexField");
        if (!obj)
                py_error("Error getting ComplexField attribute");
        obj = PyObject_CallObject(obj, NULL);
        if (!obj)
                py_error("Error getting CC attribute");
        Py_INCREF(obj);
        GiNaC::CC = obj;
}

static PyObject* pyfunc_Integer = nullptr;

void ginac_pyinit_Integer(PyObject* f) {
        Py_INCREF(f);
        pyfunc_Integer = f;
}

PyObject* Integer(const long int& x) {
        if (initialized)
                return GiNaC::py_funcs.py_integer_from_long(x);

        // Slow version since we can't call Cython-exported code yet.
        PyObject* m = PyImport_ImportModule("sage.rings.integer");
        if (!m)
                py_error("Error importing sage.rings.integer");
        PyObject* Integer = PyObject_GetAttrString(m, "Integer");
        if (!Integer)
                py_error("Error getting Integer attribute");
        PyObject* ans = PyObject_CallFunction(Integer, const_cast<char*> ("l"), x);
        Py_DECREF(m);
        Py_DECREF(Integer);
        return ans;
}

namespace GiNaC {

numeric I;
PyObject *RR, *CC;

///////////////////////////////////////////////////////////////////////////////
// class numeric
///////////////////////////////////////////////////////////////////////////////

PyObject* ZERO = PyInt_FromLong(0); // todo: never freed
PyObject* ONE = PyInt_FromLong(1); // todo: never freed
PyObject* TWO = PyInt_FromLong(2); // todo: never freed

std::ostream& operator<<(std::ostream& os, const numeric& s) {
        PyObject* o;
        switch (s.t) {
                case MPZ:
                {
                        std::vector<char> cp(2 + mpz_sizeinbase(s.v._bigint, 10));
                        mpz_get_str(&cp[0], 10, s.v._bigint);
                        return os << &cp[0];
                }
                case MPQ:
                {
                        size_t size = mpz_sizeinbase(mpq_numref(s.v._bigrat), 10)
                             + mpz_sizeinbase(mpq_denref(s.v._bigrat), 10) + 5;
                        std::vector<char> cp(size);
                        mpq_get_str(&cp[0], 10, s.v._bigrat);
                        return os << &cp[0];
                }
                case DOUBLE:
                        return os << s.v._double;
                case PYOBJECT:
                        // TODO: maybe program around Python's braindead L suffix?
                        // PyLong_Check(s.v._pyobject) 
                        o = PyObject_Repr(s.v._pyobject);
                        if (!o) {
                                PyErr_Clear();
                                throw (std::runtime_error(
                                        "operator<<(ostream, numeric): exception printing python object"));
                        } else {
#if PY_MAJOR_VERSION < 3
                                os << PyString_AsString(o);
#else
                                os << PyObject_REPR(o);
#endif
                                Py_DECREF(o);
                        }
                        return os;
                default:
                        stub("operator <<: type not yet handled");
        }
}

const numeric& numeric::operator=(const numeric& x) {
        switch (t) {
                case MPZ:
                        mpz_clear(v._bigint);
                        break;
                case MPQ:
                        mpq_clear(v._bigrat);
                        break;
                case PYOBJECT:
                        Py_DECREF(v._pyobject);
                        break;
                case DOUBLE:
                        break;
        }
        t = x.t;
        hash = x.hash;
        switch (x.t) {
                case MPZ:
                        mpz_init(v._bigint);
                        mpz_set(v._bigint, x.v._bigint);
                        break;
                case MPQ:
                        mpq_init(v._bigrat);
                        mpq_set(v._bigrat, x.v._bigrat);
                        break;
                case PYOBJECT:
                        v = x.v;
                        Py_INCREF(v._pyobject);
                        break;
                case DOUBLE:
                        break;
        }
        return *this;
}

//numeric numeric::operator()(const int& x) {
//        verbose("operator()(const int)");
//        return numeric(x);
//}

int numeric::compare_same_type(const numeric& right) const {
        verbose("compare_same_type");
        if (t != right.t) {
                numeric a, b;
                coerce(a, b, *this, right);
                return a.compare_same_type(b);
        }
        int ret;
        switch (t) {
                case DOUBLE:
                        return (v._double < right.v._double) ? -1 : (v._double > right.v._double);
                case MPZ:
                        ret = mpz_cmp(v._bigint, right.v._bigint);
                        if (ret > 0)
                                ret = 1;
                        else if (ret < 0)
                                ret = -1;
                        return ret;
                case MPQ:
                        ret = mpq_cmp(v._bigrat, right.v._bigrat);
                        if (ret > 0)
                                ret = 1;
                        else if (ret < 0)
                                ret = -1;
                        return ret;
                case PYOBJECT:
                        return Pynac_PyObj_Cmp(v._pyobject, right.v._pyobject, "compare_same_type");
                default:
                        stub("invalid type: compare_same_type type not handled");
        }

}

/* By convention hashes of PyObjects must be identical
   with their Python hashes, this applies to our MPZ
   and MPQ objects too. */
static long _mpz_pythonhash(mpz_t the_int)
{
    mp_limb_t h1=0, h0;
    size_t n = mpz_size(the_int);
    for (unsigned i=0; i<n; ++i) {
        h0 = h1;
        h1 += mpz_getlimbn(the_int, i);
        if (h1 < h0)
            ++h1;
        }
    long h = h1;
    if (mpz_sgn(the_int) < 0)
        h = -h;
    if (h == -1)
        return -2;
    return h;
}

static long _mpq_pythonhash(mpq_t the_rat)
{
    mpq_t rat;
    mpq_init(rat);
    mpq_set(rat, the_rat);
    long n = _mpz_pythonhash(mpq_numref(rat));
    long d = _mpz_pythonhash(mpq_denref(rat));
    if (d != 1L)
        n = n + (d-1) * 7461864723258187525;
    mpq_clear(rat);
    return n;
}

///////////////////////////////////////////////////////////////////////////////
// class numeric
///////////////////////////////////////////////////////////////////////////////

GINAC_IMPLEMENT_REGISTERED_CLASS_OPT(numeric, basic,
        print_func<print_context>(&numeric::do_print).
        print_func<print_latex>(&numeric::do_print_latex).
        print_func<print_csrc>(&numeric::do_print_csrc).
        print_func<print_tree>(&numeric::do_print_tree).
        print_func<print_python_repr>(&numeric::do_print_python_repr))


//////////
// default constructor
//////////

/** default constructor. Numerically it initializes to an integer zero. */
numeric::numeric() : basic(&numeric::tinfo_static) {
        t = MPZ;
        mpz_init(v._bigint);
        mpz_set_ui(v._bigint, 0);
        hash = _mpz_pythonhash(v._bigint);
        setflag(status_flags::evaluated | status_flags::expanded);
}

//////////
// other constructors
//////////

// public

numeric::numeric(const numeric& other) : basic(&numeric::tinfo_static) {
        t = other.t;
        hash = other.hash;
        switch (t) {
                case PYOBJECT:
                        v = other.v;
                        Py_INCREF(v._pyobject);
                        return;
                case DOUBLE:
                        v = other.v;
                        return;
                case MPZ:
                        mpz_init(v._bigint);
                        mpz_set(v._bigint, other.v._bigint);
                        return;
                case MPQ:
                        mpq_init(v._bigrat);
                        mpq_set(v._bigrat, other.v._bigrat);
                        return;
        }
}

numeric::numeric(PyObject* o, bool force_py) : basic(&numeric::tinfo_static) {
        if (!o) py_error("Error");
        if (not force_py) {
                if (PyInt_Check(o)) {
                        t = MPZ;
                        mpz_init(v._bigint);
                        mpz_set_si(v._bigint, PyInt_AS_LONG(o));
                        hash = _mpz_pythonhash(v._bigint);
                        setflag(status_flags::evaluated | status_flags::expanded);
                        Py_DECREF(o);
                        return;
                }
                if (initialized) {
                        if (py_funcs.py_is_Integer(o)) {
                                t = MPZ;
                                mpz_init_set(v._bigint, py_funcs.py_mpz_from_integer(o));
                                hash = _mpz_pythonhash(v._bigint);
                                setflag(status_flags::evaluated | status_flags::expanded);
                                Py_DECREF(o);
                                return;
                        }
                        else if (py_funcs.py_is_Rational(o)) {
                                t = MPQ;
                                mpq_init(v._bigrat);
                                mpq_set(v._bigrat, py_funcs.py_mpq_from_rational(o));
                                hash = _mpq_pythonhash(v._bigrat);
                                setflag(status_flags::evaluated | status_flags::expanded);
                                Py_DECREF(o);
                                return;
                        }
                }
        }

        t = PYOBJECT;
        hash = (long)PyObject_Hash(o);
        if (hash == -1 && PyErr_Occurred()) {
            // error is thrown on first hash request
            PyErr_Clear();
            is_hashable = false;
            }
        v._pyobject = o; // STEAL a reference
        setflag(status_flags::evaluated | status_flags::expanded);

}

numeric::numeric(int i) : basic(&numeric::tinfo_static) {
        t = MPZ;
        mpz_init(v._bigint);
        mpz_set_si(v._bigint, i);
        hash = _mpz_pythonhash(v._bigint);
        setflag(status_flags::evaluated | status_flags::expanded);
}

numeric::numeric(unsigned int i) : basic(&numeric::tinfo_static) {
        t = MPZ;
        mpz_init(v._bigint);
        mpz_set_ui(v._bigint, i);
        hash = _mpz_pythonhash(v._bigint);
        setflag(status_flags::evaluated | status_flags::expanded);
}

numeric::numeric(long i) : basic(&numeric::tinfo_static) {
        t = MPZ;
        mpz_init(v._bigint);
        mpz_set_si(v._bigint, i);
        hash = _mpz_pythonhash(v._bigint);
        setflag(status_flags::evaluated | status_flags::expanded);
}

numeric::numeric(unsigned long i) : basic(&numeric::tinfo_static) {
        t = MPZ;
        mpz_init(v._bigint);
        mpz_set_ui(v._bigint, i);
        hash = _mpz_pythonhash(v._bigint);
        setflag(status_flags::evaluated | status_flags::expanded);
}

numeric::numeric(mpz_t bigint) : basic(&numeric::tinfo_static) {
        t = MPZ;
        mpz_init_set(v._bigint, bigint);
        mpz_clear(bigint);
        hash = _mpz_pythonhash(v._bigint);
        setflag(status_flags::evaluated | status_flags::expanded);
}

numeric::numeric(mpq_t bigrat) : basic(&numeric::tinfo_static) {
        t = MPQ;
        mpq_init(v._bigrat);
        mpq_set(v._bigrat, bigrat);
        hash = _mpq_pythonhash(v._bigrat);
        mpq_clear(bigrat);
        setflag(status_flags::evaluated | status_flags::expanded);
}

/** Constructor for rational numerics a/b.
 *
 *  @exception overflow_error (division by zero) */
numeric::numeric(long num, long den) : basic(&numeric::tinfo_static) {
        if (!den)
                throw std::overflow_error("numeric::div(): division by zero");
        if ((num%den) == 0) {
                t = MPZ;
                mpz_init(v._bigint);
                mpz_set_si(v._bigint, num/den);                
                hash = _mpz_pythonhash(v._bigint);
        }
        else
        {
                t = MPQ;
                mpq_init(v._bigrat);
                mpq_set_si(v._bigrat, num, den);
                mpq_canonicalize(v._bigrat);
                hash = _mpq_pythonhash(v._bigrat);
        }
        setflag(status_flags::evaluated | status_flags::expanded);
}

numeric::numeric(double d) : basic(&numeric::tinfo_static) {
        t = PYOBJECT;
        if (!(v._pyobject = PyFloat_FromDouble(d)))
                py_error("Error creating double");
        setflag(status_flags::evaluated | status_flags::expanded);
}

numeric::~numeric() {
        switch (t) {
                case PYOBJECT:
                        Py_DECREF(v._pyobject);
                        return;
                case DOUBLE:
                        return;
                case MPZ:
                        mpz_clear(v._bigint);
                        return;
                case MPQ:
                        mpq_clear(v._bigrat);
                        return;
        }       
}

//////////
// archiving
//////////

numeric::numeric(const archive_node &n, lst &sym_lst) :
inherited(n, sym_lst) {
        // get type information
        unsigned int t_tmp;
        if (!n.find_unsigned(std::string("T"), t_tmp))
                throw std::runtime_error("archive error: cannot read type info");
        t = Type(t_tmp);
        std::string str;
        // read object to a string
        if (!n.find_string("S", str))
                throw (std::runtime_error("archive error: cannot read object data"));
        PyObject *arg;
        switch (t) {
                case MPZ:
                        mpz_init(v._bigint);
                        mpz_set_str(v._bigint, str.c_str(), 10);
                        hash = _mpz_pythonhash(v._bigint);
                        return;
                case MPQ:
                        mpq_init(v._bigrat);
                        mpq_set_str(v._bigrat, str.c_str(), 10);
                        hash = _mpq_pythonhash(v._bigrat);
                        return;
                case PYOBJECT:
                        // read pickled python object to a string
                        if (!n.find_string("S", str))
                                throw (std::runtime_error("archive error: cannot read pyobject data"));
                        arg = Py_BuildValue("s#", str.c_str(), str.size());
                        // unpickle
                        v._pyobject = py_funcs.py_loads(arg);
                        Py_DECREF(arg);
                        if (PyErr_Occurred()) {
                                throw (std::runtime_error("archive error: caught exception in py_loads"));
                        }
                        hash = (long)PyObject_Hash(v._pyobject);
                        if (hash == -1 && PyErr_Occurred()) {
                            PyErr_Clear();
                            is_hashable = false;
                            }
                        Py_INCREF(v._pyobject);
                        return;
                default:
                        stub("unarchiving numeric");
                        return;
        }
        setflag(status_flags::evaluated | status_flags::expanded);
}

void numeric::archive(archive_node &n) const {
        // store type information
        n.add_unsigned("T", t);

        // create a string representation of this object
        std::string *tstr;
        switch (t) {
                case MPZ:
                {
                        std::vector<char> cp(2 + mpz_sizeinbase(v._bigint, 10));
                        mpz_get_str(&cp[0], 10, v._bigint);
                        tstr = new std::string(&cp[0]);
                        break;
                }
                case MPQ:
                {
                        size_t size = mpz_sizeinbase(mpq_numref(v._bigrat), 10)
                             + mpz_sizeinbase(mpq_denref(v._bigrat), 10) + 5;
                        std::vector<char> cp(size);
                        mpq_get_str(&cp[0], 10, v._bigrat);
                        tstr = new std::string(&cp[0]);
                        break;
                }
                case PYOBJECT:
                        tstr = py_funcs.py_dumps(v._pyobject);
                        if (PyErr_Occurred()) {
                                throw (std::runtime_error("archive error: exception in py_dumps"));
                        }
                        break;
                default:
                        stub("archive numeric");
        }

        n.add_string("S", *tstr);
        delete tstr;
        inherited::archive(n);
}

DEFAULT_UNARCHIVE(numeric)

//////////
// functions overriding virtual functions from base classes
//////////

template<typename T1, typename T2>
static inline bool coerce(T1& dst, const T2& arg);

void numeric::print_numeric(const print_context & c, const char*,
                            const char*, const char*, const char*,
                            unsigned level, bool latex = false) const {
        std::string* out;
        if (latex) {
                out = py_funcs.py_latex(to_pyobject(), level);
        } else {
                out = py_funcs.py_repr(to_pyobject(), level);
        }
        c.s << *out;
        delete out;
        return;
}

void numeric::do_print(const print_context & c, unsigned level) const {
        print_numeric(c, "(", ")", "I", "*", level, false);
}

void numeric::do_print_latex(const print_latex & c, unsigned level) const {
        print_numeric(c, "{(", ")}", "i", " ", level, true);
}

void numeric::do_print_csrc(const print_csrc & c, unsigned) const {
        // TODO: not really needed?
        stub("print_csrc");
}

void numeric::do_print_tree(const print_tree & c, unsigned level) const {
        c.s << std::string(level, ' ') << *this
                << " (" << class_name() << ")" << " @" << this
                << std::hex << ", hash=0x" << hashvalue << ", flags=0x" << flags << std::dec
                << std::endl;
}

void numeric::do_print_python_repr(const print_python_repr & c, unsigned level) const {
        c.s << class_name() << "('";
        print_numeric(c, "(", ")", "I", "*", level);
        c.s << "')";
}

bool numeric::info(unsigned inf) const {
        switch (inf) {
                case info_flags::numeric:
                case info_flags::polynomial:
                case info_flags::rational_function:
                case info_flags::expanded:
                        return true;
                case info_flags::real:
                        return is_real();
                case info_flags::rational:
                case info_flags::rational_polynomial:
                        return is_rational();
                case info_flags::crational:
                case info_flags::crational_polynomial:
                        return is_crational();
                case info_flags::integer:
                case info_flags::integer_polynomial:
                        return is_integer();
                case info_flags::cinteger:
                case info_flags::cinteger_polynomial:
                        return is_cinteger();
                case info_flags::positive:
                        return is_positive();
                case info_flags::negative:
                        return is_negative();
                case info_flags::nonnegative:
                        return is_zero() || is_positive();
                case info_flags::nonzero:
                        return not is_zero();
                case info_flags::posint:
                        return is_pos_integer();
                case info_flags::negint:
                        return is_integer() && is_negative();
                case info_flags::nonnegint:
                        return is_nonneg_integer();
                case info_flags::even:
                        return is_even();
                case info_flags::odd:
                        return is_odd();
                case info_flags::prime:
                        return is_prime();
                case info_flags::algebraic:
                        return !is_real();
                case info_flags::infinity:
                case info_flags::has_indices:
                        return false;
                case info_flags::inexact:
                        return not is_exact();
                default:
                        throw (std::runtime_error("numeric::info()"));
        }
        return false;
}

bool numeric::is_polynomial(const ex & var) const {
        return true;
}

int numeric::degree(const ex & s) const {
        // In sage deg (0) != 0 !!!
        return 0;
}

int numeric::ldegree(const ex & s) const {
        return 0;
}

ex numeric::coeff(const ex & s, int n) const {
        return n == 0 ? * this : _ex0;
}

/** Disassemble real part and imaginary part to scan for the occurrence of a
 *  single number.  Also handles the imaginary unit.  It ignores the sign on
 *  both this and the argument, which may lead to what might appear as funny
 *  results:  (2+I).has(-2) -> true.  But this is consistent, since we also
 *  would like to have (-2+I).has(2) -> true and we want to think about the
 *  sign as a multiplicative factor. */
bool numeric::has(const ex &other, unsigned) const {
        if (!is_exactly_a<numeric>(other))
                return false;
        const numeric &o = ex_to<numeric>(other);
        if (this->is_equal(o) || this->is_equal(-o))
                return true;
        if (o.imag().is_zero()) { // e.g. scan for 3 in -3*I
                if (!this->real().is_equal(*_num0_p))
                        if (this->real().is_equal(o) || this->real().is_equal(-o))
                                return true;
                if (!this->imag().is_equal(*_num0_p))
                        if (this->imag().is_equal(o) || this->imag().is_equal(-o))
                                return true;
                return false;
        } else {
                if (o.is_equal(I)) // e.g scan for I in 42*I
                        return !this->is_real();
                if (o.real().is_zero()) // e.g. scan for 2*I in 2*I+1
                        if (!this->imag().is_equal(*_num0_p))
                                if (this->imag().is_equal(o * I) || this->imag().is_equal(-o * I))
                                        return true;
        }
        return false;
}

/** Evaluation of numbers doesn't do anything at all. */
ex numeric::eval(int) const {
        // Warning: if this is ever gonna do something, the ex constructors from all kinds
        // of numbers should be checking for status_flags::evaluated.
        return this->hold();
}

/** Cast numeric into a floating-point object.  For example exact numeric(1) is
 *  returned as a 1.0000000000000000000000 and so on according to how Digits is
 *  currently set.  In case the object already was a floating point number the
 *  precision is trimmed to match the currently set default.
 *
 *  @param level  ignored, only needed for overriding basic::evalf.
 *  @return  an ex-handle to a numeric. */
ex numeric::evalf(int, PyObject* parent) const {
        PyObject *a = to_pyobject();
        PyObject *ans = py_funcs.py_float(a, parent);
        Py_DECREF(a);
        if (!ans)
                throw (std::runtime_error("numeric::evalf(): error calling py_float()"));

        return ans;
}

ex numeric::conjugate() const {
        PY_RETURN(py_funcs.py_conjugate);
}

ex numeric::real_part() const {
        return real();
}

ex numeric::imag_part() const {
        return imag();
}

// protected

/** This method establishes a canonical order on all numbers.  For complex
 *  numbers this is not possible in a mathematically consistent way but we need
 *  to establish some order and it ought to be fast.  So we simply define it
 *  to be compatible with our method csgn.
 *
 *  @return csgn(*this-other) */
int numeric::compare_same_type(const basic &other) const {
        GINAC_ASSERT(is_exactly_a<numeric>(other));
        const numeric &o = static_cast<const numeric &> (other);

        return compare_same_type(o);
}

bool numeric::is_equal_same_type(const basic &other) const {
        GINAC_ASSERT(is_exactly_a<numeric>(other));
        const numeric &o = static_cast<const numeric &> (other);

        return this->is_equal(o);
}

long numeric::calchash() const {
        switch (t) {
                case DOUBLE:
                        return (long) v._double;
                case MPZ:
                case MPQ:
                case PYOBJECT:
                        if (is_hashable)
                            return hash;
                        throw (std::runtime_error("Python object not hashable"));
                default:
                        stub("invalid type: ::hash() type not handled");
        }
}


//////////
// new virtual functions which can be overridden by derived classes
//////////

// none

//////////
// non-virtual functions in this class
//////////

// public

/** Numerical addition method.  Adds argument to *this and returns result as
 *  a numeric object. */
const numeric numeric::add(const numeric &other) const {
        verbose("operator+");
        // We will have to write a coercion model :-(
        // Or just a nested switch of switches that covers all
        // combinations of long,double,mpfr,mpz,mpq,mpfc,mpqc.  Yikes.
        if (t != other.t) {
                numeric a, b;
                coerce(a, b, *this, other);
                return a + b;
        }
        switch (t) {
                case DOUBLE:
                        return v._double + other.v._double;
                case MPZ:
                        mpz_t bigint;
                        mpz_init(bigint);
                        mpz_add(bigint, v._bigint, other.v._bigint);
                        return bigint;
                case MPQ:
                        mpq_t bigrat;
                        mpq_init(bigrat);
                        mpq_add(bigrat, v._bigrat, other.v._bigrat);
                        return bigrat;
                case PYOBJECT:
                        return PyNumber_Add(v._pyobject, other.v._pyobject);
                default:
                        stub("invalid type: operator+() type not handled");
        }
}

/** Numerical subtraction method.  Subtracts argument from *this and returns
 *  result as a numeric object. */
const numeric numeric::sub(const numeric &other) const {
        verbose("operator-");
        if (t != other.t) {
                numeric a, b;
                coerce(a, b, *this, other);
                return a - b;
        }
        switch (t) {
                case DOUBLE:
                        return v._double - other.v._double;
                case MPZ:
                        mpz_t bigint;
                        mpz_init(bigint);
                        mpz_sub(bigint, v._bigint, other.v._bigint);
                        return bigint;
                case MPQ:
                        mpq_t bigrat;
                        mpq_init(bigrat);
                        mpq_sub(bigrat, v._bigrat, other.v._bigrat);
                        return bigrat;
                case PYOBJECT:
                        return PyNumber_Subtract(v._pyobject, other.v._pyobject);
                default:
                        stub("invalid type: operator-() type not handled");
        }
}

/** Numerical multiplication method.  Multiplies *this and argument and returns
 *  result as a numeric object. */
const numeric numeric::mul(const numeric &other) const {
        verbose("operator*");
        if (t != other.t) {
                numeric a, b;
                coerce(a, b, *this, other);
                return a * b;
        }
        switch (t) {
                case DOUBLE:
                        return v._double * other.v._double;
                case MPZ:
                        mpz_t bigint;
                        mpz_init(bigint);
                        mpz_mul(bigint, v._bigint, other.v._bigint);
                        return bigint;
                case MPQ:
                        mpq_t bigrat;
                        mpq_init(bigrat);
                        mpq_mul(bigrat, v._bigrat, other.v._bigrat);
                        return bigrat;
                case PYOBJECT:
                        return PyNumber_Multiply(v._pyobject, other.v._pyobject);
                default:
                        stub("invalid type: operator*() type not handled");
        }
}

/** Numerical division method.  Divides *this by argument and returns result as
 *  a numeric object.
 *
 *  @exception overflow_error (division by zero) */
const numeric numeric::div(const numeric &other) const {
        verbose("operator/");
        if (other.is_zero())
                throw std::overflow_error("numeric::div(): division by zero");
        if (t != other.t) {
                numeric a, b;
                coerce(a, b, *this, other);
                return a / b;
        }
        switch (t) {
                case DOUBLE:
                        return v._double / other.v._double;
                case MPZ:
                        if (mpz_divisible_p(v._bigint, other.v._bigint)) {
                                mpz_t bigint;
                                mpz_init(bigint);
                                mpz_divexact(bigint, v._bigint, other.v._bigint);
                                return bigint;
                        }
                        else {
                                mpq_t bigrat, obigrat;
                                mpq_init(bigrat);
                                mpq_init(obigrat);
                                mpq_set_z(bigrat, v._bigint);
                                mpq_set_z(obigrat, other.v._bigint);
                                mpq_div(bigrat, bigrat, obigrat);
                                mpq_clear(obigrat);
                                return bigrat;
                        }
                case MPQ:
                                mpq_t bigrat;
                                mpq_init(bigrat);
                                mpq_div(bigrat, v._bigrat, other.v._bigrat);
                                return bigrat;

                case PYOBJECT:
#if PY_MAJOR_VERSION < 3
                        if (PyObject_Compare(other.v._pyobject, ONE) == 0) {
                                return *this;
                        }
                        if (PyInt_Check(v._pyobject)) {
                                if (PyInt_Check(other.v._pyobject)) {
                                        // This branch happens at startup.
                                        PyObject* o = PyNumber_TrueDivide(Integer(PyInt_AS_LONG(v._pyobject)),
                                                Integer(PyInt_AS_LONG(other.v._pyobject)));
                                        // I don't 100% understand why I have to incref this, 
                                        // but if I don't, Sage crashes on exit.
                                        Py_INCREF(o);
                                        return o;
                                } else if (PyLong_Check(other.v._pyobject)) {
                                        PyObject* d = py_funcs.py_integer_from_python_obj(other.v._pyobject);
                                        PyObject* ans = PyNumber_TrueDivide(v._pyobject, d);
                                        Py_DECREF(d);
                                        return ans;
                                }
                        } else if (PyLong_Check(v._pyobject)) {
#else
                        if (PyLong_Check(v._pyobject)) {
#endif
                                PyObject* n = py_funcs.py_integer_from_python_obj(v._pyobject);
                                PyObject* ans = PyNumber_TrueDivide(n, other.v._pyobject);
                                Py_DECREF(n);
                                return ans;
                        }
                        return PyNumber_TrueDivide(v._pyobject, other.v._pyobject);

                default:
                        stub("invalid type: operator/() type not handled");
        }
}

/** Numerical exponentiation.  Raises *this to the integer power given as argument and
 *  returns result as numeric. */
const numeric numeric::power(const numeric &exponent) const {
        verbose("pow");
        numeric ex(exponent);
        if (exponent.t == PYOBJECT and PyInt_Check(exponent.v._pyobject)) {
                ex.t = MPZ;
                long si = PyInt_AsLong(exponent.v._pyobject);
                mpz_set_si(ex.v._bigint, si);
        }
        if (ex.t == MPZ) {
                if (not mpz_fits_sint_p(ex.v._bigint)) {
                        throw std::runtime_error("numeric::power(): exponent doesn't fit in signed long");
                }
                signed long int exp_si = mpz_get_si(ex.v._bigint);
                PyObject *o, *r;
                switch (t) {
                        case DOUBLE:
                                return ::pow(v._double, double(exp_si));
                        case MPZ:
                                if (exp_si >= 0) {
                                        mpz_t bigint;
                                        mpz_init(bigint);
                                        mpz_pow_ui(bigint, v._bigint, exp_si);
                                        return bigint;
                                }
                                else {
                                        mpz_t bigint;
                                        mpz_init_set(bigint, v._bigint);
                                        mpz_pow_ui(bigint, bigint, -exp_si);
                                        mpq_t bigrat;
                                        mpq_init(bigrat);
                                        mpq_set_z(bigrat, bigint);
                                        mpq_inv(bigrat, bigrat);
                                        mpz_clear(bigint);
                                        return bigrat;
                                }
                        case MPQ:
                                mpz_t bigint;
                                mpq_t bigrat, obigrat;
                                mpz_init(bigint);
                                mpq_init(bigrat);
                                mpq_init(obigrat);
                                if (exp_si >= 0) {
                                        mpz_pow_ui(bigint, mpq_numref(v._bigrat), exp_si);
                                        mpq_set_z(bigrat, bigint);
                                        mpz_pow_ui(bigint, mpq_denref(v._bigrat), exp_si);
                                        mpq_set_z(obigrat, bigint);
                                        mpq_div(bigrat, bigrat, obigrat);
                                }
                                else {
                                        mpz_pow_ui(bigint, mpq_denref(v._bigrat), -exp_si);
                                        mpq_set_z(bigrat, bigint);
                                        mpz_pow_ui(bigint, mpq_numref(v._bigrat), -exp_si);
                                        mpq_set_z(obigrat, bigint);
                                        mpq_div(bigrat, bigrat, obigrat);
                                }
                                mpz_clear(bigint);
                                mpq_clear(obigrat);
                                return bigrat;
                        case PYOBJECT:
                                o = Integer(exp_si);
                                r = PyNumber_Power(v._pyobject, o, Py_None);
                                Py_DECREF(o);
                                return r;
                        default:
                                stub("invalid type: pow numeric");
                }
        }
        if (t != exponent.t) {
                numeric a, b;
                coerce(a, b, *this, exponent);
                return pow(a, b);
        }
        if (t == MPQ) {
                mpq_t basis;
                mpq_init(basis);
                mpq_set(basis, v._bigrat);
                PyObject *r1 = py_funcs.py_rational_from_mpq(basis);
                PyObject *r2 = py_funcs.py_rational_from_mpq(ex.v._bigrat);
                PyObject *r = PyNumber_Power(r1, r2, Py_None);
                Py_DECREF(r1);
                Py_DECREF(r2);
                mpq_clear(basis);
                numeric p(r, true);
                return p;
        }
        return PyNumber_Power(v._pyobject, exponent.v._pyobject, Py_None);
}


/** Numerical addition method.  Adds argument to *this and returns result as
 *  a numeric object on the heap.  Use internally only for direct wrapping into
 *  an ex object, where the result would end up on the heap anyways. */
const numeric &numeric::add_dyn(const numeric &other) const {
        // Efficiency shortcut: trap the neutral element by pointer.  This hack
        // is supposed to keep the number of distinct numeric objects low.
        if (this == _num0_p)
                return other;
        else if (&other == _num0_p)
                return *this;

        return static_cast<const numeric &> ((new numeric(*this + other))->
                setflag(status_flags::dynallocated));
}

/** Numerical subtraction method.  Subtracts argument from *this and returns
 *  result as a numeric object on the heap.  Use internally only for direct
 *  wrapping into an ex object, where the result would end up on the heap
 *  anyways. */
const numeric &numeric::sub_dyn(const numeric &other) const {
        // Efficiency shortcut: trap the neutral exponent (first by pointer).  This
        // hack is supposed to keep the number of distinct numeric objects low.
        if (&other == _num0_p || (other.is_zero()))
                return *this;

        return static_cast<const numeric &> ((new numeric(*this - other))->
                setflag(status_flags::dynallocated));
}

/** Numerical multiplication method.  Multiplies *this and argument and returns
 *  result as a numeric object on the heap.  Use internally only for direct
 *  wrapping into an ex object, where the result would end up on the heap
 *  anyways. */
const numeric &numeric::mul_dyn(const numeric &other) const {
        // Efficiency shortcut: trap the neutral element by pointer.  This hack
        // is supposed to keep the number of distinct numeric objects low.
        if (this == _num1_p)
                return other;
        else if (&other == _num1_p)
                return *this;

        return static_cast<const numeric &> ((new numeric(*this * other))->
                setflag(status_flags::dynallocated));
}

/** Numerical division method.  Divides *this by argument and returns result as
 *  a numeric object on the heap.  Use internally only for direct wrapping
 *  into an ex object, where the result would end up on the heap
 *  anyways.
 *
 *  @exception overflow_error (division by zero) */
const numeric &numeric::div_dyn(const numeric &other) const {
        // Efficiency shortcut: trap the neutral element by pointer.  This hack
        // is supposed to keep the number of distinct numeric objects low.
        if (&other == _num1_p)
                return *this;
        if (other.is_zero())
                throw std::overflow_error("division by zero");
        return static_cast<const numeric &> ((new numeric(*this / other))->
                setflag(status_flags::dynallocated));
}

/** Numerical exponentiation.  Raises *this to the power given as argument and
 *  returns result as a numeric object on the heap.  Use internally only for
 *  direct wrapping into an ex object, where the result would end up on the
 *  heap anyways. */
const numeric &numeric::power_dyn(const numeric &other) const {
        // Efficiency shortcut: trap the neutral exponent (first try by pointer, then
        // try harder, since calls to cln::expt() below may return amazing results for
        // floating point exponent 1.0).
        if (&other == _num1_p || (other == *_num1_p))
                return *this;

        return static_cast<const numeric &> ((new numeric(pow(*this, other)))->
                setflag(status_flags::dynallocated));
}

const numeric &numeric::operator=(int i) {
        return operator=(numeric(i));
}

const numeric &numeric::operator=(unsigned int i) {
        return operator=(numeric(i));
}

const numeric &numeric::operator=(long i) {
        return operator=(numeric(i));
}

const numeric &numeric::operator=(unsigned long i) {
        return operator=(numeric(i));
}

const numeric &numeric::operator=(double d) {
        return operator=(numeric(d));
}

const numeric numeric::negative() const {
        verbose("operator-");
        switch (t) {
                case DOUBLE:
                        return -v._double;
                case MPZ:
                        mpz_t bigint;
                        mpz_init_set(bigint, v._bigint);
                        mpz_neg(bigint, bigint);
                        return bigint;
                case MPQ:
                        mpq_t bigrat;
                        mpq_init(bigrat);
                        mpq_set(bigrat, v._bigrat);
                        mpq_neg(bigrat, bigrat);
                        return bigrat;
                case PYOBJECT:
                        return PyNumber_Negative(v._pyobject);
                default:
                        stub("invalid type: operator-() type not handled");
        }
}

/** Inverse of a number. */
const numeric numeric::inverse() const {
        if (is_zero())
                throw std::overflow_error("numeric::inverse(): division by zero");
        return numeric(1) / *this;
}

/** Return the step function of a numeric. The imaginary part of it is
 *  ignored because the step function is generally considered real but
 *  a numeric may develop a small imaginary part due to rounding errors.
 */
const numeric numeric::step() const {
        PY_RETURN(py_funcs.py_step);
}

/** Return the complex half-plane (left or right) in which the number lies.
 *  csgn(x)==0 for x==0, csgn(x)==1 for Re(x)>0 or Re(x)=0 and Im(x)>0,
 *  csgn(x)==-1 for Re(x)<0 or Re(x)=0 and Im(x)<0.
 *
 *  @see numeric::compare(const numeric &other) */
int numeric::csgn() const {
        verbose("csgn");
        switch (t) {
                case DOUBLE:
                        if (v._double < 0)
                                return -1;
                        if (v._double == 0)
                                return 0;
                        return 1;
                case MPZ:
                        return mpz_sgn(v._bigint);
                case MPQ:
                        return mpq_sgn(v._bigrat);
                case PYOBJECT:
                        int result;
                        if (is_real()) {
                                result = Pynac_PyObj_Cmp(v._pyobject, ZERO, "csgn");
                        } else {
                                PyObject *tmp = py_funcs.py_real(v._pyobject);
                                result = Pynac_PyObj_Cmp(tmp, ZERO, "csgn");
                                Py_DECREF(tmp);
                                if (result == 0) {
                                        tmp = py_funcs.py_imag(v._pyobject);
                                        result = Pynac_PyObj_Cmp(tmp, ZERO, "csgn");
                                        Py_DECREF(tmp);
                                }
                        }
                        return result;
                default:
                        stub("invalid type: csgn() type not handled");
        }
}

int numeric::compare(const numeric &other) const {
        return (*this-other).csgn();
}

bool numeric::is_equal(const numeric &other) const {
        return *this == other;
}

/** True if object is zero. */
bool numeric::is_zero() const {
        verbose("is_zero");
        int a;
        switch (t) {
                case DOUBLE:
                        return v._double == 0;
                case MPZ:
                        return mpz_cmp_si(v._bigint, 0) == 0;
                case MPQ:
                        return mpq_cmp_si(v._bigrat, 0, 1) == 0;
                case PYOBJECT:
                        a = PyObject_Not(v._pyobject);
                        if (a == -1)
                                py_error("is_zero");
                        return a == 1;
                default:
                        std::cerr << "type = " << t << "\n";
                        stub("invalid type: is_zero() type not handled");
        }
}

/** True if object is not complex and greater than zero. */
bool numeric::is_positive() const {
        verbose("is_positive");
        switch (t) {
                case DOUBLE:
                        return v._double > 0;
                case MPZ:
                        return mpz_cmp_si(v._bigint, 0) > 0;
                case MPQ:
                        return mpq_cmp_si(v._bigrat, 0, 1) > 0;
                case PYOBJECT:
                        return is_real() and Pynac_PyObj_RichCmp(v._pyobject, ZERO, Py_GT, "is_positive");
                default:
                        stub("invalid type: is_positive() type not handled");
        }
}

/** True if object is not complex and less than zero. */
bool numeric::is_negative() const {
        verbose("is_negative");
        switch (t) {
                case DOUBLE:
                        return v._double < 0;
                case MPZ:
                        return mpz_cmp_si(v._bigint, 0) < 0;
                case MPQ:
                        return mpq_cmp_si(v._bigrat, 0, 1) < 0;
                case PYOBJECT:
                        return is_real() and Pynac_PyObj_RichCmp(v._pyobject, ZERO, Py_LT, "is_negative");
                default:
                        stub("invalid type: is_negative() type not handled");
        }
}

/** True if object is a non-complex integer. */
bool numeric::is_integer() const {
        verbose2("is_integer", *this);

        bool ret;
        switch (t) {
                case DOUBLE:
                        return false;
                case MPZ:
                        return true;
                case MPQ:
                        mpq_t bigrat;
                        mpq_init(bigrat);
                        mpq_set(bigrat, v._bigrat);
                        mpq_canonicalize(bigrat);
                        ret = mpz_cmp_ui(mpq_denref(bigrat), 1) ==0;
                        mpq_clear(bigrat);
                        return ret;
                case PYOBJECT:
                        return py_funcs.py_is_integer(v._pyobject) != 0;
                default:
                        stub("invalid type: is_integer() type not handled");
        }
}

/** True if object is an exact integer greater than zero. */
bool numeric::is_pos_integer() const {
        verbose("is_pos_integer");
        switch (t) {
                case DOUBLE:
                        return false;
                case MPZ:
                        return is_positive();
                case MPQ:
                        return (is_integer() && is_positive());
                case PYOBJECT:
                        return (is_integer() && is_positive());
                default:
                        stub("invalid type: is_pos_integer() type not handled");
        }
}

/** True if object is an exact integer greater or equal zero. */
bool numeric::is_nonneg_integer() const {
        verbose("is_nonneg_integer");
        switch (t) {
                case DOUBLE:
                        return false;
                case MPZ:
                        return is_positive() or is_zero();
                case MPQ:
                        return (is_integer() and (is_positive() or is_zero()));
                case PYOBJECT:
                        return is_integer() and Pynac_PyObj_RichCmp(v._pyobject, ZERO, Py_GE, "is_nonneg_integer");
                default:
                        stub("invalid type: is_nonneg_integer() type not handled");
        }
}

/** True if object is an exact even integer. */
bool numeric::is_even() const {
        verbose("is_even");

        if (!is_integer())
                return false;

        switch (t) {
                case DOUBLE:
                        return false;
                case MPZ:
                        return mpz_tstbit(v._bigint, 0) == 0;
                case MPQ:
                        return is_integer()
                                and mpz_tstbit(mpq_numref(v._bigrat), 0) == 0;
                case PYOBJECT:
                        return py_funcs.py_is_even(v._pyobject) != 0;
                default:
                        stub("invalid type: is_even() type not handled");
        }
}

/** True if object is an exact odd integer. */
bool numeric::is_odd() const {
        switch (t) {
                case DOUBLE:
                        return false;
                case MPZ:
                        return mpz_tstbit(v._bigint, 0) == 1;
                case MPQ:
                        return is_integer()
                                and mpz_tstbit(mpq_numref(v._bigrat), 0) == 1;
                case PYOBJECT:
                        return !is_even();
                default:
                        stub("invalid type: is_odd() type not handled");
        }
}

/** Probabilistic primality test.
 *
 *  @return  true if object is exact integer and prime. */
bool numeric::is_prime() const {
        verbose("is_prime");
        switch (t) {
                case DOUBLE:
                        return false;
                case MPZ:
                        return mpz_probab_prime_p(v._bigint, 25) > 0;
                case MPQ:
                        return is_integer() and
                                mpz_probab_prime_p(mpq_numref(v._bigrat), 25) > 0;
                case PYOBJECT:
                        return py_funcs.py_is_prime(v._pyobject) != 0;
                default:
                        stub("invalid type: is_prime() type not handled");
        }
}

/** True if object is an exact rational number, may even be complex
 *  (denominator may be unity). */
bool numeric::is_rational() const {
        verbose("is_rational");
        switch (t) {
                case DOUBLE:
                        return false;
                case MPZ:
                        return true;
                case MPQ:
                        return true;
                case PYOBJECT:
                        return py_funcs.py_is_rational(v._pyobject) != 0;
                default:
                        stub("invalid type -- is_rational() type not handled");
        }
}

/** True if object is a real integer, rational or float (but not complex). */
bool numeric::is_real() const {
        verbose("is_real");
        switch (t) {
                case DOUBLE:
                case MPZ:
                        return true;
                case MPQ:
                        return true;
                case PYOBJECT:
                        return py_funcs.py_is_real(v._pyobject) != 0;
                default:
                        stub("invalid type -- is_real() type not handled");
        }
}

/** True if the parent of the object has positive characteristic. */
bool numeric::is_parent_pos_char() const {
        return get_parent_char() > 0;
}

/** Returns the characteristic of the parent of this object. */
int numeric::get_parent_char() const {
        verbose("get_parent_char");
        switch (t) {
                case DOUBLE:
                        return 0;
                case MPZ:
                        return 0;
                case MPQ:
                        return 0;
                case PYOBJECT:
                {
                        int c = py_funcs.py_get_parent_char(v._pyobject);
                        if (c == -1) py_error("error in py_get_parent_char");
                        return c;
                }
                default:
                        stub("invalid type -- is_parent_pos_char() type not handled");
        }
}

/** True if object is an exact rational number, may even be complex
 *  (denominator may be unity). */
bool numeric::is_exact() const {
        verbose("is_exact");
        switch (t) {
                case DOUBLE:
                        return false;
                case MPZ:
                        return true;
                case MPQ:
                        return true;
                case PYOBJECT:
                        return py_funcs.py_is_exact(v._pyobject) != 0;
                default:
                        stub("invalid type -- is_exact() type not handled");
        }
}

bool numeric::operator==(const numeric &right) const {
        verbose3("operator==", *this, right);
        if (t != right.t) {
                numeric a, b;
                coerce(a, b, *this, right);
                return a == b;
        }
        switch (t) {
                case DOUBLE:
                        return v._double == right.v._double;
                case MPZ:
                        return mpz_cmp(v._bigint, right.v._bigint) ==0;
                case MPQ:
                        return mpq_equal(v._bigrat, right.v._bigrat) !=0;
                case PYOBJECT:
                        return py_funcs.py_is_equal(v._pyobject, right.v._pyobject) != 0;
                default:
                        stub("invalid type: operator== type not handled");
        }
}

bool numeric::operator!=(const numeric &right) const {
        verbose("operator!=");
        if (t != right.t) {
                numeric a, b;
                coerce(a, b, *this, right);
                return a != b;
        }
        switch (t) {
                case DOUBLE:
                        return v._double != right.v._double;
                case MPZ:
                        return mpz_cmp(v._bigint, right.v._bigint) !=0;
                case MPQ:
                        return mpq_equal(v._bigrat, right.v._bigrat) ==0;
                case PYOBJECT:
                        return (!py_funcs.py_is_equal(v._pyobject, right.v._pyobject));
                default:
                        stub("invalid type: operator!= type not handled");
        }
}

/** True if object is element of the domain of integers extended by I, i.e. is
 *  of the form a+b*I, where a and b are integers. */
bool numeric::is_cinteger() const {
        verbose("is_crational");
        switch (t) {
                case DOUBLE:
                        return false;
                case MPZ:
                        return true;
                case MPQ:
                        return is_integer();
                case PYOBJECT:
                        return py_funcs.py_is_cinteger(v._pyobject) != 0;
                default:
                        stub("invalid type -- is_cinteger() type not handled");
        }
}

/** True if object is an exact rational number, may even be complex
 *  (denominator may be unity). */
bool numeric::is_crational() const {
        verbose("is_crational");
        switch (t) {
                case DOUBLE:
                        return false;
                case MPZ:
                        return true;
                case MPQ:
                        return true;
                case PYOBJECT:
                        return py_funcs.py_is_crational(v._pyobject) != 0;
                default:
                        stub("invalid type -- is_crational() type not handled");
        }
}

/** Numerical comparison: less.
 *
 *  @exception invalid_argument (complex inequality) */
bool numeric::operator<(const numeric &right) const {
        verbose("operator<");
        if (t != right.t) {
                numeric a, b;
                coerce(a, b, *this, right);
                return a < b;
        }

        switch (t) {
                case DOUBLE:
                        return v._double < right.v._double;
                case MPZ:
                        return mpz_cmp(v._bigint, right.v._bigint) < 0;
                case MPQ:
                        return mpq_cmp(v._bigrat, right.v._bigrat) < 0;
                case PYOBJECT:
                        return Pynac_PyObj_RichCmp(v._pyobject, right.v._pyobject, Py_LT, "<");
                default:
                        stub("invalid type: operator< type not handled");
        }
}

/** Numerical comparison: less or equal.
 *
 *  @exception invalid_argument (complex inequality) */
bool numeric::operator<=(const numeric &right) const {
        verbose("operator<=");
        if (t != right.t) {
                numeric a, b;
                coerce(a, b, *this, right);
                return a <= b;
        }
        switch (t) {
                case DOUBLE:
                        return v._double <= right.v._double;
                case MPZ:
                        return mpz_cmp(v._bigint, right.v._bigint) <= 0;
                case MPQ:
                        return mpq_cmp(v._bigrat, right.v._bigrat) <= 0;
                case PYOBJECT:
                        return Pynac_PyObj_RichCmp(v._pyobject, right.v._pyobject, Py_LE, "<=");
                default:
                        stub("invalid type: operator<= type not handled");
        }
}

/** Numerical comparison: greater.
 *
 *  @exception invalid_argument (complex inequality) */
bool numeric::operator>(const numeric &right) const {
        verbose("operator>");
        if (t != right.t) {
                numeric a, b;
                coerce(a, b, *this, right);
                return a > b;
        }
        switch (t) {
                case DOUBLE:
                        return v._double > right.v._double;
                case MPZ:
                        return mpz_cmp(v._bigint, right.v._bigint) > 0;
                case MPQ:
                        return mpq_cmp(v._bigrat, right.v._bigrat) > 0;
                case PYOBJECT:
                        return Pynac_PyObj_RichCmp(v._pyobject, right.v._pyobject, Py_GT, ">");
                default:
                        stub("invalid type: operator> type not handled");
        }
}

/** Numerical comparison: greater or equal.
 *
 *  @exception invalid_argument (complex inequality) */
bool numeric::operator>=(const numeric &right) const {
        verbose("operator>=");
        if (t != right.t) {
                numeric a, b;
                coerce(a, b, *this, right);
                return a >= b;
        }
        switch (t) {
                case DOUBLE:
                        return v._double >= right.v._double;
                case MPZ:
                        return mpz_cmp(v._bigint, right.v._bigint) >= 0;
                case MPQ:
                        return mpq_cmp(v._bigrat, right.v._bigrat) >= 0;
                case PYOBJECT:
                        return Pynac_PyObj_RichCmp(v._pyobject, right.v._pyobject, Py_GE, ">=");
                default:
                        stub("invalid type: operator!= type not handled");
        }
}

/** Converts numeric types to machine's long.  You should check with
 *  is_integer() if the number is really an integer before calling this method.
 *  You may also consider checking the range first. */
long numeric::to_long() const {
        GINAC_ASSERT(this->is_integer());
        verbose("operator long int");
        signed long int n;
        switch (t) {
                case DOUBLE:
                        return (long int) v._double;
                case MPZ:
                        return (long int) mpz_get_si(v._bigint);
                case MPQ:
                        mpz_t bigint;
                        mpz_init(bigint);
                        mpz_fdiv_q(bigint, mpq_numref(v._bigrat), mpq_denref(v._bigrat));
                        n = mpz_get_si(bigint);
                        mpz_clear(bigint);
                        return n;
                case PYOBJECT:
                        n = PyInt_AsLong(v._pyobject);
                        if (n == -1 && PyErr_Occurred()) {
                                PyErr_Print();
                                py_error("Overfloat converting to long int");
                        }
                        return n;
                default:
                        stub("invalid type: operator long int() type not handled");
        }
}

/* Return the underlying Python object corresponding to this
   numeric.  If this numeric isn't implemented using a Python
   object, the corresponding Python object is constructed on
   the fly.

   Returns a NEW REFERENCE.
 */
PyObject* numeric::to_pyobject() const {
        // Returns a New Reference
        PyObject* o;
        switch (t) {
                case MPZ:
                        mpz_t bigint;
                        mpz_init_set(bigint, v._bigint);
                        o = py_funcs.py_integer_from_mpz(bigint);
                        mpz_clear(bigint);
                        return o;
                case MPQ:
                        mpq_t bigrat;
                        mpq_init(bigrat);
                        mpq_set(bigrat, v._bigrat);
                        mpq_canonicalize(bigrat);
                        o = py_funcs.py_rational_from_mpq(bigrat);
                        mpq_clear(bigrat);
                        return o;
                case DOUBLE:
                        if (!(o = PyFloat_FromDouble(v._double)))
                                py_error("Error creating double");
                        return o;
                        //if (!(o = PyObject_CallFunction(pyfunc_Float, "d", x.v._double))) {
                        //  py_error("Error coercing a long to an Integer");

                case PYOBJECT:
                        Py_INCREF(v._pyobject);
                        return v._pyobject;

                default:
                        std::cout << t << std::endl;
                        stub("numeric::to_pyobject -- not able to do conversion to pyobject; everything else will be nonsense");
                        return nullptr;
        }
}


/** Converts numeric types to machine's double. You should check with is_real()
 *  if the number is really not complex before calling this method. */
double numeric::to_double() const {
        GINAC_ASSERT(this->is_real());
        verbose("operator double");
        double d;
        switch (t) {
                case DOUBLE:
                        return v._double;
                case MPZ:
                        return mpz_get_d(v._bigint);
                case MPQ:
                        return mpq_get_d(v._bigrat);
                case PYOBJECT:
                        d = PyFloat_AsDouble(v._pyobject);
                        if (d == -1 && PyErr_Occurred())
                                py_error("Error converting to a double.");
                        return d;
                default:
                        std::cerr << "type = " << t << std::endl;
                        stub("invalid type: operator double() type not handled");
        }
}

/** Real part of a number. */
const numeric numeric::real() const {
        verbose("real_part(a)");
        PyObject *ans;
        switch (t) {
                case DOUBLE:
                        return *this;
                case MPZ:
                        return *this;
                case MPQ:
                        return *this;
                case PYOBJECT:
                        ans = py_funcs.py_real(v._pyobject);
                        if (!ans) py_error("real_part");
                        return ans;
                default:
                        std::cerr << "type = " << t << std::endl;
                        stub("invalid type: operator double() type not handled");
        }
}

/** Imaginary part of a number. */
const numeric numeric::imag() const {
        if (is_real())
                return 0;
        verbose("imag_part(a)");
        PyObject *a = to_pyobject();
        PyObject *ans = py_funcs.py_imag(a);
        if (!ans) py_error("imag_part");
        Py_DECREF(a);
        return ans;
}

/** Numerator.  Computes the numerator of rational numbers, rationalized
 *  numerator of complex if real and imaginary part are both rational numbers
 *  (i.e numer(4/3+5/6*I) == 8+5*I), the number carrying the sign in all other
 *  cases. */
const numeric numeric::numer() const {
        verbose2("numer -- in:", *this);
        numeric ans;
        PyObject* a;

        switch (t) {

                case DOUBLE:
                        return *this;
                case MPZ:
                        return *this;
                case MPQ:
                        mpz_t bigint;
                        mpz_init_set(bigint, mpq_numref(v._bigrat));
                        return bigint;
                case PYOBJECT:
                        a = py_funcs.py_numer(v._pyobject);
                        if (!a) py_error("numer");
                        ans = a;
                        break;
                default:
                        stub("invalid type -- numer() type not handled");
                        ans = *this;
        }
        verbose2("numer -- out:", ans);
        return ans;
}

/** Denominator.  Computes the denominator of rational numbers, common integer
 *  denominator of complex if real and imaginary part are both rational numbers
 *  (i.e denom(4/3+5/6*I) == 6), one in all other cases. */
const numeric numeric::denom() const {
        verbose2("denom -- in:", *this);
        numeric ans;
        PyObject* a;

        switch (t) {
                case DOUBLE:
                case MPZ:
                        return 1;
                case MPQ:
                        mpz_t bigint;
                        mpz_init_set(bigint, mpq_denref(v._bigrat));
                        return bigint;
                case PYOBJECT:
                        a = py_funcs.py_denom(v._pyobject);
                        if (!a) py_error("denom");
                        ans = a;
                        break;

                default:
                        stub("invalid type -- denom() type not handled");
                        ans = ONE;
        }
        verbose2("denom -- out:", ans);
        return ans;
}

const numeric numeric::fibonacci() const {
    PY_RETURN(py_funcs.py_fibonacci);
}

const numeric numeric::sin() const {
        PY_RETURN(py_funcs.py_sin);
}

const numeric numeric::cos() const {
        PY_RETURN(py_funcs.py_cos);
}

const numeric numeric::zeta() const {
        PY_RETURN(py_funcs.py_zeta);
}

const numeric numeric::exp() const {
        PY_RETURN(py_funcs.py_exp);
}

const numeric numeric::log() const {
        PY_RETURN(py_funcs.py_log);
}

const numeric numeric::tan() const {
        PY_RETURN(py_funcs.py_tan);
}

const numeric numeric::asin() const {
        PY_RETURN(py_funcs.py_asin);
}

const numeric numeric::acos() const {
        PY_RETURN(py_funcs.py_acos);
}

const numeric numeric::atan() const {
        PY_RETURN(py_funcs.py_atan);
}

const numeric numeric::atan(const numeric& y) const {
        PY_RETURN2(py_funcs.py_atan2, y);
}

const numeric numeric::sinh() const {
        PY_RETURN(py_funcs.py_sinh);
}

const numeric numeric::cosh() const {
        PY_RETURN(py_funcs.py_cosh);
}

const numeric numeric::tanh() const {
        PY_RETURN(py_funcs.py_tanh);
}

const numeric numeric::asinh() const {
        PY_RETURN(py_funcs.py_asinh);
}

const numeric numeric::acosh() const {
        PY_RETURN(py_funcs.py_acosh);
}

const numeric numeric::atanh() const {
        PY_RETURN(py_funcs.py_atanh);
}

const numeric numeric::Li2(const numeric &n, PyObject* parent) const {
        PyObject *aa = to_pyobject();
        PyObject* nn = n.to_pyobject();
        PyObject *ans = py_funcs.py_li(aa, nn, parent);
        if (!ans) py_error("error calling function");
        Py_DECREF(aa);
        Py_DECREF(nn);
        return ans;
}

const numeric numeric::Li2() const {
        PY_RETURN(py_funcs.py_li2);
}

const numeric numeric::lgamma() const {
        PY_RETURN(py_funcs.py_lgamma);
}

const numeric numeric::tgamma() const {
        PY_RETURN(py_funcs.py_tgamma);
}

const numeric numeric::psi() const {
        PY_RETURN(py_funcs.py_psi);
}

const numeric numeric::psi(const numeric& y) const {
        PY_RETURN2(py_funcs.py_psi2, y);
}

const numeric numeric::factorial() const {
        PY_RETURN(py_funcs.py_factorial);
}

const numeric numeric::doublefactorial() const {
        PY_RETURN(py_funcs.py_doublefactorial);
}

const numeric numeric::binomial(const numeric &k) const {
        PY_RETURN2(py_funcs.py_binomial, k);
}

const numeric numeric::bernoulli() const {
        PY_RETURN(py_funcs.py_bernoulli);
}

const numeric numeric::isqrt() const {
        PY_RETURN(py_funcs.py_isqrt);
}

const numeric numeric::sqrt() const {
        PY_RETURN(py_funcs.py_sqrt);
}

const numeric numeric::abs() const {
        if (t == MPZ) {
                mpz_t bigint;
                mpz_init(bigint);
                mpz_abs(bigint, v._bigint);
                return bigint;
        }
        else if (t == MPQ) {
                mpq_t bigrat;
                mpq_init(bigrat);
                mpq_abs(bigrat, v._bigrat);
                return bigrat;
        }
        
        PY_RETURN(py_funcs.py_abs);
}

const numeric numeric::mod(const numeric &b) const {
        if (t == MPZ) {
                mpz_t bigint;
                mpz_init(bigint);
                mpz_mod(bigint, v._bigint, b.v._bigint);
                return bigint;
        }
        
        PY_RETURN2(py_funcs.py_mod, b);
}

const numeric numeric::_smod(const numeric &b) const {
        PY_RETURN2(py_funcs.py_smod, b);
}

const numeric numeric::irem(const numeric &b) const {
        if (t == MPZ) {
                mpz_t bigint;
                mpz_init(bigint);
                mpz_fdiv_r(bigint, v._bigint, b.v._bigint);
                return bigint;
        }
        
        PY_RETURN2(py_funcs.py_irem, b);
}

const numeric numeric::iquo(const numeric &b) const {
        if (t == MPZ) {
                mpz_t bigint;
                mpz_init(bigint);
                mpz_fdiv_q(bigint, v._bigint, b.v._bigint);
                return bigint;
        }
        
        PY_RETURN2(py_funcs.py_iquo, b);
}

const numeric numeric::iquo(const numeric &b, numeric& r) const {
        PY_RETURN3(py_funcs.py_iquo2, b, r);
}

const numeric numeric::gcd(const numeric &b) const {
        if (t == MPZ and b.t == MPZ) {
                mpz_t bigint;
                mpz_init(bigint);
                mpz_gcd(bigint, v._bigint, b.v._bigint);
                return bigint;
        }
        
        PY_RETURN2(py_funcs.py_gcd, b);
}

const numeric numeric::lcm(const numeric &b) const {
        if (t == MPZ and b.t == MPZ) {
                mpz_t bigint;
                mpz_init(bigint);
                mpz_lcm(bigint, v._bigint, b.v._bigint);
                return bigint;
        }
        
        PY_RETURN2(py_funcs.py_lcm, b);
}

/** Size in binary notation.  For integers, this is the smallest n >= 0 such
 *  that -2^n <= x < 2^n. If x > 0, this is the unique n > 0 such that
 *  2^(n-1) <= x < 2^n.
 *
 *  @return  number of bits (excluding sign) needed to represent that number
 *  in two's complement if it is an integer, 0 otherwise. */
int numeric::int_length() const {
        PyObject* a = to_pyobject();
        int n = py_funcs.py_int_length(a);
        Py_DECREF(a);
        if (n == -1)
                py_error("int_length");
        return n;
}

//////////
// global 
//////////

void coerce(numeric& new_left, numeric& new_right, const numeric& left, const numeric& right) {
        verbose("coerce");
        // Return a numeric 
        if (left.t == right.t) {
                new_left = left;
                new_right = right;
                return;
        }
        mpq_t bigrat;
        PyObject *o;
        switch (left.t) {
                case MPZ:
                        switch (right.t) {
                                case DOUBLE:
                                        new_left = left.to_double();
                                        new_right = right;
                                        return;
                                case MPQ:
                                        mpq_init(bigrat);
                                        mpq_set_z(bigrat, left.v._bigint);
                                        new_left = numeric(bigrat);
                                        new_right = right;
                                        return;
                                case PYOBJECT:
                                        mpz_t bigint;
                                        mpz_init_set(bigint, left.v._bigint);
                                        o = py_funcs.py_integer_from_mpz(bigint);
                                        new_left = numeric(o, true);
                                        new_right = right;
                                        mpz_clear(bigint);
                                        return;
                                default:
                                        std::cerr << "type = " << right.t << "\n";
                                        stub("** invalid coercion -- left MPZ**");
                        }
                case MPQ:
                        switch (right.t) {
                                case DOUBLE:
                                        new_left = left.to_double();
                                        new_right = right;
                                        return;
                                case MPZ:
                                        mpq_init(bigrat);
                                        mpq_set_z(bigrat, right.v._bigint);
                                        new_left = left;
                                        new_right = numeric(bigrat);
                                        return;
                                case PYOBJECT:
                                        mpq_init(bigrat);
                                        mpq_set(bigrat, left.v._bigrat);
                                        o = py_funcs.py_rational_from_mpq(bigrat);
                                        mpq_clear(bigrat);
                                        new_left = numeric(o, true);
                                        new_right = right;
                                        return;
                                default:
                                        std::cerr << "type = " << right.t << "\n";
                                        stub("** invalid coercion -- left MPQ**");
                        }
                case DOUBLE:
                        switch (right.t) {
                                case MPZ:
                                        new_left = left;
                                        new_right = right.to_double();
                                        return;
                                case MPQ:
                                        new_left = left;
                                        new_right = right.to_double();
                                        return;
                                case PYOBJECT:
                                        new_left = PyFloat_FromDouble(left.to_double());
                                        new_right = right;
                                        return;
                                default:
                                        std::cerr << "type = " << right.t << "\n";
                                        stub("** invalid coercion -- left DOUBLE ** ");
                        }
                case PYOBJECT:
                        new_left = left;
                        switch (right.t) {
                                case MPZ:
                                        mpz_t bigint;
                                        mpz_init_set(bigint, right.v._bigint);
                                        o = py_funcs.py_integer_from_mpz(bigint);
                                        mpz_clear(bigint);
                                        new_right = numeric(o, true);
                                        return;
                                case MPQ:
                                        mpq_init(bigrat);
                                        mpq_set(bigrat, right.v._bigrat);
                                        o = py_funcs.py_rational_from_mpq(bigrat);
                                        mpq_clear(bigrat);
                                        new_right = numeric(o, true);
                                        return;
                                case DOUBLE:
                                        new_left = PyFloat_FromDouble(right.to_double());
                                        return;
                                default:
                                        std::cerr << "type = " << right.t << "\n";
                                        stub("** invalid coercion -- left PYOBJECT**");
                        }
        }
        std::cerr << "type = " << left.t << "\n";
        stub("** invalid coercion **");
}

/** Imaginary unit.  This is not a constant but a numeric since we are
 *  natively handing complex numbers anyways, so in each expression containing
 *  an I it is automatically eval'ed away anyhow. */

/** Exponential function.
 *
 *  @return  arbitrary precision numerical exp(x). */
const numeric exp(const numeric &x) {
        return x.exp();
}

/** Natural logarithm.
 *
 *  @param x complex number
 *  @return  arbitrary precision numerical log(x).
 *  @exception pole_error("log(): logarithmic pole",0) */
const numeric log(const numeric &x) {
        /*
  if (x.is_zero())
    throw pole_error("log(): logarithmic pole",0);
         */
        return x.log();
}

/** Numeric sine (trigonometric function).
 *
 *  @return  arbitrary precision numerical sin(x). */
const numeric sin(const numeric &x) {
        return x.sin();
}

/** Numeric cosine (trigonometric function).
 *
 *  @return  arbitrary precision numerical cos(x). */
const numeric cos(const numeric &x) {
        return x.cos();
}

/** Numeric tangent (trigonometric function).
 *
 *  @return  arbitrary precision numerical tan(x). */
const numeric tan(const numeric &x) {
        return x.tan();
}

/** Numeric inverse sine (trigonometric function).
 *
 *  @return  arbitrary precision numerical asin(x). */
const numeric asin(const numeric &x) {
        return x.asin();
}

/** Numeric inverse cosine (trigonometric function).
 *
 *  @return  arbitrary precision numerical acos(x). */
const numeric acos(const numeric &x) {
        return x.acos();
}

/** Numeric arcustangent.
 *
 *  @param x complex number
 *  @return atan(x)
 *  @exception pole_error("atan(): logarithmic pole",0) if x==I or x==-I. */
const numeric atan(const numeric &x) {
        if (!x.is_real() &&
                x.real().is_zero() &&
                abs(x.imag()).is_equal(*_num1_p))
                throw pole_error("atan(): logarithmic pole", 0);
        return x.atan();
}

/** Numeric arcustangent of two arguments, analytically continued in a suitable way.
 *
 *  @param y complex number
 *  @param x complex number
 *  @return -I*log((x+I*y)/sqrt(x^2+y^2)), which is equal to atan(y/x) if y and
 *    x are both real.
 *  @exception pole_error("atan(): logarithmic pole",0) if y/x==+I or y/x==-I. */
const numeric atan(const numeric &y, const numeric &x) {
        return x.atan(y);
}

/** Numeric hyperbolic sine (trigonometric function).
 *
 *  @return  arbitrary precision numerical sinh(x). */
const numeric sinh(const numeric &x) {
        return x.sinh();
}

/** Numeric hyperbolic cosine (trigonometric function).
 *
 *  @return  arbitrary precision numerical cosh(x). */
const numeric cosh(const numeric &x) {
        return x.cosh();
}

/** Numeric hyperbolic tangent (trigonometric function).
 *
 *  @return  arbitrary precision numerical tanh(x). */
const numeric tanh(const numeric &x) {
        return x.tanh();
}

/** Numeric inverse hyperbolic sine (trigonometric function).
 *
 *  @return  arbitrary precision numerical asinh(x). */
const numeric asinh(const numeric &x) {
        return x.asinh();
}

/** Numeric inverse hyperbolic cosine (trigonometric function).
 *
 *  @return  arbitrary precision numerical acosh(x). */
const numeric acosh(const numeric &x) {
        return x.acosh();
}

/** Numeric inverse hyperbolic tangent (trigonometric function).
 *
 *  @return  arbitrary precision numerical atanh(x). */
const numeric atanh(const numeric &x) {
        return x.atanh();
}


/** Numeric evaluation of Dilogarithm within circle of convergence (unit
 *  circle) using a power series. */

/** Numeric evaluation of Dilogarithm.  The domain is the entire complex plane,
 *  the branch cut lies along the positive real axis, starting at 1 and
 *  continuous with quadrant IV.
 *
 *  @return  arbitrary precision numerical Li2(x). */
const numeric Li(const numeric &x) {
        return x.Li2();
}

const numeric Li2(const numeric &x) {
        return x.Li2();
}

const numeric Li2(const numeric &x, const numeric &n, PyObject* parent) {
        return x.Li2(n, parent);
}

/** Evaluation of Riemann's Zeta function.  */
const numeric zeta(const numeric &x) {
        return x.zeta();
}

/** The Gamma function.
 *  Use the Lanczos approximation. If the coefficients used here are not
 *  sufficiently many or sufficiently accurate, more can be calculated
 *  using the program doc/examples/lanczos.cpp. In that case, be sure to
 *  read the comments in that file. */
const numeric lgamma(const numeric &x) {
        return x.lgamma();
}

const numeric tgamma(const numeric &x) {
        return x.tgamma();
}

/** The psi function (aka polygamma function). */
const numeric psi(const numeric &x) {
        return x.psi();
}

/** The psi functions (aka polygamma functions). */
const numeric psi(const numeric &n, const numeric &x) {
        return n.psi(x);
}

/** Factorial combinatorial function.
 *
 *  @param n  integer argument >= 0
 *  @exception range_error (argument must be integer >= 0) */
const numeric factorial(const numeric &n) {
        return n.factorial();
}

/** The double factorial combinatorial function.  (Scarcely used, but still
 *  useful in cases, like for exact results of tgamma(n+1/2) for instance.)
 *
 *  @param n  integer argument >= -1
 *  @return n!! == n * (n-2) * (n-4) * ... * ({1|2}) with 0!! == (-1)!! == 1
 *  @exception range_error (argument must be integer >= -1) */
const numeric doublefactorial(const numeric &n) {
        return n.doublefactorial();
}

/** The Binomial coefficients.  It computes the binomial coefficients.  For
 *  integer n and k and positive n this is the number of ways of choosing k
 *  objects from n distinct objects.  If n is negative, the formula
 *  binomial(n,k) == (-1)^k*binomial(k-n-1,k) is used to compute the result. */
const numeric binomial(const numeric &n, const numeric &k) {
        return n.binomial(k);
}

/** Bernoulli number.  The nth Bernoulli number is the coefficient of x^n/n!
 *  in the expansion of the function x/(e^x-1).
 *
 *  @return the nth Bernoulli number (a rational number).
 *  @exception range_error (argument must be integer >= 0) */
const numeric bernoulli(const numeric &n) {
        return n.bernoulli();
}

/** Fibonacci number.  The nth Fibonacci number F(n) is defined by the
 *  recurrence formula F(n)==F(n-1)+F(n-2) with F(0)==0 and F(1)==1.
 *
 *  @param n an integer
 *  @return the nth Fibonacci number F(n) (an integer number)
 *  @exception range_error (argument must be an integer) */
const numeric fibonacci(const numeric &n) {
        return n.fibonacci();
}

/** Absolute value. */
const numeric abs(const numeric& x) {
        return x.abs();
}

/** Modulus (in positive representation).
 *  In general, mod(a,b) has the sign of b or is zero, and rem(a,b) has the
 *  sign of a or is zero. This is different from Maple's modp, where the sign
 *  of b is ignored. It is in agreement with Mathematica's Mod.
 *
 *  @return a mod b in the range [0,abs(b)-1] with sign of b if both are
 *  integer, 0 otherwise. */
const numeric mod(const numeric &a, const numeric &b) {
        return a.mod(b);
}

///** Modulus (in symmetric representation).
// *  Equivalent to Maple's mods.
// *
// *  @return a mod b in the range [-iquo(abs(b)-1,2), iquo(abs(b),2)]. */
const numeric smod(const numeric &a, const numeric &b) {
        return a._smod(b);
}

/** Numeric integer remainder.
 *  Equivalent to Maple's irem(a,b) as far as sign conventions are concerned.
 *  In general, mod(a,b) has the sign of b or is zero, and irem(a,b) has the
 *  sign of a or is zero.
 *
 *  @return remainder of a/b if both are integer, 0 otherwise.
 *  @exception overflow_error (division by zero) if b is zero. */
const numeric irem(const numeric &a, const numeric &b) {
        return a.irem(b);
}

/** Numeric integer quotient.
 *  Equivalent to Maple's iquo as far as sign conventions are concerned.
 *  
 *  @return truncated quotient of a/b if both are integer, 0 otherwise.
 *  @exception overflow_error (division by zero) if b is zero. */
const numeric iquo(const numeric &a, const numeric &b) {
        return a.iquo(b);
}

/** Numeric integer quotient.
 *  Equivalent to Maple's iquo(a,b,'r') it obeyes the relation
 *  r == a - iquo(a,b,r)*b.
 *
 *  @return truncated quotient of a/b and remainder stored in r if both are
 *  integer, 0 otherwise.
 *  @exception overflow_error (division by zero) if b is zero. */
const numeric iquo(const numeric &a, const numeric &b, numeric &r) {
        return a.iquo(b, r);
}

/** Greatest Common Divisor.
 *   
 *  @return  The GCD of two numbers if both are integer, a numerical 1
 *  if they are not. */
const numeric gcd(const numeric &a, const numeric &b) {
        return a.gcd(b);
}

/** Least Common Multiple.
 *   
 *  @return  The LCM of two numbers if both are integer, the product of those
 *  two numbers if they are not. */
const numeric lcm(const numeric &a, const numeric &b) {
        return a.lcm(b);
}

/** Numeric square root.
 *  If possible, sqrt(x) should respect squares of exact numbers, i.e. sqrt(4)
 *  should return integer 2.
 *
 *  @param x numeric argument
 *  @return square root of x. Branch cut along negative real axis, the negative
 *  real axis itself where imag(x)==0 and real(x)<0 belongs to the upper part
 *  where imag(x)>0. */
const numeric sqrt(const numeric &x) {
        return x.sqrt();
}

/** Integer numeric square root. */
const numeric isqrt(const numeric &x) {
        return x.isqrt();
}

/** Floating point evaluation of Sage's constants. */
ex ConstantEvalf(unsigned serial, PyObject* dict) {
        if (dict == nullptr) {
                dict = PyDict_New();
                PyDict_SetItemString(dict, "parent", CC);
        }
        PyObject* x = py_funcs.py_eval_constant(serial, dict);
        if (!x) py_error("error getting digits of constant");
        return x;
}

ex UnsignedInfinityEvalf(unsigned serial, PyObject* parent) {
        PyObject* x = py_funcs.py_eval_unsigned_infinity();
        return x;
}

ex InfinityEvalf(unsigned serial, PyObject* parent) {
        PyObject* x = py_funcs.py_eval_infinity();
        return x;
}

ex NegInfinityEvalf(unsigned serial, PyObject* parent) {
        PyObject* x = py_funcs.py_eval_neg_infinity();
        return x;
}


} // namespace GiNaC
