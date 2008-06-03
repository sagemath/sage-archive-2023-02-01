cdef extern from "Python.h":
    ctypedef void PyObject
    ctypedef void PyTypeObject
    ctypedef struct FILE
    ctypedef void _typeobject
    ctypedef void (*freefunc)(void *)
    ctypedef void (*destructor)(PyObject *)
    #ctypedef int (*printfunc)(PyObject *, FILE *, int)
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
    ctypedef PyObject *(*newfunc)(PyTypeObject *, PyObject *, PyObject *)
    ctypedef PyObject *(*allocfunc)(PyTypeObject *, Py_ssize_t)

    # We need a PyTypeObject with elements so we can
    # get and set tp_new, tp_dealloc, tp_flags, and tp_basicsize
    ctypedef struct RichPyTypeObject "PyTypeObject":

        #We replace this one
        allocfunc tp_alloc
        newfunc tp_new
        freefunc tp_free
        destructor tp_dealloc

        #sizeof(Object)
        Py_ssize_t tp_basicsize

        long tp_flags

    # We need a PyObject where we can get/set the refcnt directly
    # and access the type.
    ctypedef struct RichPyObject "PyObject":
        int ob_refcnt
        RichPyTypeObject* ob_type

    cdef long Py_TPFLAGS_HAVE_GC

    # Allocation
    RichPyObject* PyObject_MALLOC(int)

    # Useful for debugging, see below
    void PyObject_INIT(RichPyObject *, RichPyTypeObject *)

    ## Free
    void PyObject_FREE(PyObject*)
