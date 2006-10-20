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
    void init_global_empty_tuple()

init_global_empty_tuple()
