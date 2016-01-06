cpdef dict_key(o)
cpdef cache_key(o)

cdef class CachedFunction(object):
    cdef public str __name__
    cdef public str __module__
    cdef _argument_fixer
    cdef public _fix_to_pos
    cdef public f
    cdef public cache  # not always of type <dict>
    cdef _default_key
    cdef bint is_classmethod
    cdef argfix_init(self)
    cdef key

cdef class CachedMethod(object):
    cdef str _cache_name
    cdef public str __name__
    cdef public str __module__
    cdef CachedFunction _cachedfunc
    cdef int nargs
    cpdef dict _get_instance_cache(self, inst)

cdef class CachedMethodCaller(CachedFunction):
    cdef public _instance
    cdef public bint _inst_in_key
    cdef public CachedMethod _cachedmethod

cdef class CachedMethodCallerNoArgs(CachedFunction):
    cdef public _instance
