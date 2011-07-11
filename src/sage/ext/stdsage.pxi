"""
Standard SAGE Pyrex Helper Code
"""

################################################################################
# stdsage.pxi
#   Standard useful stuff for SAGE modules to include:
#   See stdsage.h for macros and stdsage.c for C functions.
#
#   Each module currently gets its own copy of this, which is why
#   we call the initialization code below.
#
################################################################################

###############################################################################
#   SAGE: System for Algebra and Geometry Experimentation
#       Copyright (C) 2005, 2006 William Stein <wstein@gmail.com>
#  Distributed under the terms of the GNU General Public License (GPL)
#  The full text of the GPL is available at:
#                  http://www.gnu.org/licenses/
###############################################################################

cdef extern from "stdsage.h":
    ctypedef void PyObject

    # Global tuple -- useful optimization
    void init_global_empty_tuple()
    object PY_NEW(object t)
    object PY_NEW_SAME_TYPE(object t)

    void* PY_TYPE(object o)
    bint PY_TYPE_CHECK(object o, object t)
    bint PY_TYPE_CHECK_EXACT(object o, object t)

    object IS_INSTANCE(object o, object t)
    void PY_SET_TP_NEW(object t1, object t2)
    bint HAS_DICTIONARY(object o)
    bint PY_IS_NUMERIC(object o)


# Memory management
cdef extern from "stdsage.h":
    void  sage_free(void *p)
    void* sage_realloc(void *p, size_t n)
    void* sage_malloc(size_t)
    void* sage_calloc(size_t nmemb, size_t size)
    void  init_csage()
    void  init_csage_module()
    void  init_memory_functions()


# Do this for every single module that links in stdsage.
init_csage_module()
