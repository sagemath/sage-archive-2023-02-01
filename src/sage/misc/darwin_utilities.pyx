
IF UNAME_SYSNAME == "Darwin":
    cdef extern from "darwin_memory_usage.h":
        cdef unsigned long long darwin_virtual_size()

def darwin_memory_usage():
    r"""
    On Darwin, returns the virtual size of the process in bytes. This will match what "top" reports
    as the VSIZE. Raises on all other platforms

    EXAMPLES
        sage: from sage.misc.darwin_utilities import darwin_memory_usage
        sage: if os.uname()[0] == 'Darwin':
        ...       memory_usage = darwin_memory_usage()


        sage: from sage.misc.darwin_utilities import darwin_memory_usage
        sage: try:
        ...     if os.uname()[0] != 'Darwin':
        ...         memory_usage = darwin_memory_usage()
        ...     else:
        ...         raise NotImplementedError
        ... except NotImplementedError:
        ...     print "NotImplementedError"
        NotImplementedError

    """
    IF UNAME_SYSNAME == "Darwin":
        return darwin_virtual_size()
    ELSE:
        import os
        uname = os.uname()
        platform_name = uname[0]
        raise NotImplementedError, "darwin_memory_usage() not implemented on platform %s"%platform_name
