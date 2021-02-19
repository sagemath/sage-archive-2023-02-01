cdef class MapReduceEngine():
    cdef public list worker_results
    cdef public object input_iter
    cdef public int mp_thresh

    cpdef map_caller(self,mp_params,mapper,extra_args=*)

