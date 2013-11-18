cdef str pari_error_string

cdef void _pari_init_error_handling()
cdef void _pari_check_warning "_pari_check_warning"()
cdef int _pari_handle_exception(long err) except 0
cdef void _pari_err_recover(long err)
