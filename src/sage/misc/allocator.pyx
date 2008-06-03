include "../ext/interrupt.pxi"  # ctrl-c interrupt block support
include "../ext/stdsage.pxi"
include "../ext/python.pxi"
include "../ext/python_list.pxi"
include "../ext/python_number.pxi"
include "../ext/python_int.pxi"
include "../ext/python_rich_object.pxi"

#def time_alloc_list(n):
    #"""
    #Allocate n a list of n SAGE integers using PY_NEW.
    #(Used for timing purposes.)

    #EXAMPLES:
       #sage: from sage.rings.integer import time_alloc_list
       #sage: v = time_alloc_list(100)
    #"""
    #cdef int i
    #l = []
    #for i from 0 <= i < n:
        #l.append(PY_NEW(Integer))

    #return l

#def time_alloc(n):
    #"""
    #Time allocating n integers using PY_NEW.
    #Used for timing purposes.

    #EXAMPLES:
       #sage: from sage.rings.integer import time_alloc
       #sage: time_alloc(100)
    #"""
    #cdef int i
    #for i from 0 <= i < n:
        #z = PY_NEW(Integer)

#def pool_stats():
    #"""
    #Returns information about the Integer object pool.

    #EXAMPLES:
        #sage: from sage.rings.integer import pool_stats
        #sage: pool_stats()            # random-ish output
        #Used pool 0 / 0 times
        #Pool contains 3 / 100 items
    #"""
    #return ["Used pool %s / %s times" % (use_pool, total_alloc),
            #"Pool contains %s / %s items" % (integer_pool_count, integer_pool_size)]

cdef public hook_tp_functions(object global_dummy, allocfunc tp_alloc, newfunc tp_new, freefunc tp_free, destructor tp_dealloc, bint useGC):
    """
    Initialize the fast integer creation functions.
    """

    cdef long flag

    cdef RichPyObject* o
    o = <RichPyObject*>global_dummy

    # Make sure this never, ever gets collected
    Py_INCREF(global_dummy)

    # By default every object created in Pyrex is garbage
    # collected. This means it may have references to other objects
    # the Garbage collector has to look out for. We remove this flag
    # as the only reference an Integer has is to the global Integer
    # ring. As this object is unique we don't need to garbage collect
    # it as we always have a module level reference to it. If another
    # attribute is added to the Integer class this flag removal so as
    # the alloc and free functions may not be used anymore.
    # This object will still be reference counted.
    if not useGC:
        flag = Py_TPFLAGS_HAVE_GC
        o.ob_type.tp_flags = <long>(o.ob_type.tp_flags & (~flag))

    # Finally replace the functions called when an Integer needs
    # to be constructed/destructed.
    o.ob_type.tp_new = tp_new
    o.ob_type.tp_dealloc = tp_dealloc
