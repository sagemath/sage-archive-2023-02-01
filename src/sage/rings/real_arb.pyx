from sage.libs.arb.arb cimport *

def _test_arb_():
    """
    EXAMPLES::

        sage: from sage.rings.real_arb import _test_arb_
        sage: _test_arb_()
        (948239929664858513 * 2^-59) +/- (692955552 * 2^-58)
    """

    cdef arb_t x
    cdef arb_t y
    arb_init(x)
    arb_init(y)
    arb_set_ui(x, 2)
    arb_zeta(y, x, 53)
    arb_print(y)
    arb_clear(x)
    arb_clear(y)
