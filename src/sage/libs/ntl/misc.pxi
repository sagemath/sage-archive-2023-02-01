# Unset the signal handler and create a string from the buffer,
# then free the memory in the buffer.
cdef object string(char* s):
    _sig_off
    # Makes a python string and deletes what is pointed to by s.
    t = str(s)
    free(s)
    return t

_INIT = None

