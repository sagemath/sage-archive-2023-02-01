"""
Utility functions to handle the PARI stack and copy objects from it.
"""

#*****************************************************************************
#       Copyright (C) 2016 Luca De Feo <luca.defeo@polytechnique.edu>
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 2 of the License, or
# (at your option) any later version.
#                  http://www.gnu.org/licenses/
#*****************************************************************************

from __future__ import absolute_import, division, print_function

include "cysignals/signals.pxi"
include "cysignals/memory.pxi"

cdef extern from *:
    int sig_on_count "cysigs.sig_on_count"

from .paridecl cimport pari_mainstack, avma, paristack_setsize, gsizebyte, gcopy_avma, gnil


cdef inline void clear_stack():
    """
    Call ``sig_off()``. If we are leaving the outermost
    ``sig_on() ... sig_off()`` block, then clear the PARI stack.
    """
    global avma
    if sig_on_count <= 1:
        avma = pari_mainstack.top
    sig_off()

cdef inline GEN deepcopy_to_python_heap(GEN x, void** address) except NULL:
    cdef size_t s = gsizebyte(x)
    cdef void* tmp_bot = check_malloc(s)
    address[0] = tmp_bot
    cdef pari_sp tmp_top = <pari_sp>tmp_bot + s
    return gcopy_avma(x, &tmp_top)

cdef inline Gen new_gen(GEN x):
    """
    Create a new Gen wrapping `x`, then call ``clear_stack()``.
    Except if `x` is ``gnil``, then we return ``None`` instead.
    """
    cdef Gen g
    if x is gnil:
        g = None
    else:
        g = new_gen_noclear(x)
    clear_stack()
    return g

cdef inline Gen new_gen_noclear(GEN x):
    """
    Create a new gen, but don't free any memory on the stack and don't
    call sig_off().
    """
    cdef Gen y = Gen.__new__(Gen)
    y.g = deepcopy_to_python_heap(x, &y.chunk)
    return y
