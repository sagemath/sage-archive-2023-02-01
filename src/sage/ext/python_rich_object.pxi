from cpython.ref cimport PyObject, PyTypeObject

cdef extern from "Python.h":
    ctypedef void (*freefunc)(void *)
    ctypedef void (*destructor)(PyObject *)
    ctypedef PyObject *(*getattrfunc)(PyObject *, char *)
    ctypedef PyObject *(*getattrofunc)(PyObject *, PyObject *)
    ctypedef int (*setattrfunc)(PyObject *, char *, PyObject *)
    ctypedef int (*setattrofunc)(PyObject *, PyObject *, PyObject *)
    ctypedef int (*cmpfunc)(PyObject *, PyObject *)
    ctypedef PyObject *(*reprfunc)(PyObject *)
    ctypedef long (*hashfunc)(PyObject *)
    ctypedef PyObject *(*richcmpfunc) (PyObject *, PyObject *, int)
    ctypedef PyObject *(*getiterfunc) (PyObject *)
    ctypedef PyObject *(*iternextfunc) (PyObject *)
    ctypedef PyObject *(*descrgetfunc) (PyObject *, PyObject *, PyObject *)
    ctypedef int (*descrsetfunc) (PyObject *, PyObject *, PyObject *)
    ctypedef int (*initproc)(PyObject *, PyObject *, PyObject *)
    ctypedef object (*newfunc)(PyTypeObject *, PyObject *, PyObject *)
    ctypedef PyObject *(*allocfunc)(PyTypeObject *, Py_ssize_t)

    # We need a PyTypeObject with elements so we can
    # get and set tp_new, tp_dealloc, tp_flags, and tp_basicsize
    ctypedef struct RichPyTypeObject "PyTypeObject":
        long tp_dictoffset

        allocfunc tp_alloc
        newfunc tp_new
        freefunc tp_free
        destructor tp_dealloc
        hashfunc tp_hash
        richcmpfunc tp_richcompare

        #sizeof(Object)
        Py_ssize_t tp_basicsize

        long tp_flags


    cdef long Py_TPFLAGS_HAVE_GC

    # Allocation
    PyObject* PyObject_MALLOC(int)

    # Free
    void PyObject_FREE(PyObject*)
