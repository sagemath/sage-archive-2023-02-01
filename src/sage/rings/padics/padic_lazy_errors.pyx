from sage.rings.padics.precision_error import PrecisionError


cdef inline int ERROR_ABANDON      = 1
cdef inline int ERROR_NOTDEFINED   = 1 << 1
cdef inline int ERROR_PRECISION    = 1 << 2
cdef inline int ERROR_OVERFLOW     = 1 << 3
cdef inline int ERROR_NOTSQUARE    = 1 << 4  # maybe we should have something more generic here
cdef inline int ERROR_DIVISION     = 1 << 5
cdef inline int ERROR_CIRCULAR     = 1 << 6

cdef inline int ERROR_UNKNOWN      = 1 << 30


def raise_error(error):
    if error & ERROR_UNKNOWN:  # should never occur
        raise RuntimeError
    if error & ERROR_CIRCULAR:
        raise RecursionError("definition looks circular")
    if error & ERROR_DIVISION:
        raise ZeroDivisionError("cannot divide by something indistinguishable from zero")
    if error & ERROR_NOTSQUARE:
        raise ValueError("not a square")
    if error & ERROR_OVERFLOW:
        raise OverflowError
    if error & (ERROR_PRECISION | ERROR_NOTDEFINED):
        raise PrecisionError("not enough precision")
    if error & ERROR_ABANDON:
        raise PrecisionError("computation has been abandonned; try to increase precision by hand")


def error_to_str(error):
    if error & ERROR_UNKNOWN:  # should never occur
        raise RuntimeError
    if error & ERROR_CIRCULAR:
        return "Error: definition looks circular"
    if error & ERROR_DIVISION:
        return "Error: cannot divide by something indistinguishable from zero"
    if error & ERROR_NOTSQUARE:
        return "Error: not a square"
    if error & ERROR_OVERFLOW:
        return "Error: overflow"
    if error & (ERROR_PRECISION | ERROR_NOTDEFINED):
        return "Error: not enough precision"
    if error & ERROR_ABANDON:
        return "Abandon; try to increase precision by hand"
