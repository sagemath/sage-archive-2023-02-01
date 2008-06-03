include "../ext/python_rich_object.pxi"
cdef public hook_tp_functions(object global_dummy, allocfunc tp_alloc, newfunc tp_new, freefunc tp_free, destructor tp_dealloc, bint useGC)
