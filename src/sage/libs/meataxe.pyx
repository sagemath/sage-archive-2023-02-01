#*****************************************************************************
#       Copyright (C) 2017 Simon King <simon.king@uni-jena.de>
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 2 of the License, or
# (at your option) any later version.
#                  http://www.gnu.org/licenses/
#*****************************************************************************

from cpython.exc cimport PyErr_SetObject
from cysignals.signals cimport sig_block, sig_unblock
cdef dict ErrMsg = {
    "Not enough memory": MemoryError,
    "Time limit exceeded": RuntimeError,
    "Division by zero": ZeroDivisionError,
    "Bad file format": OSError,
    "No such file or directory": OSError,
    "Bad argument": ValueError,
    "Argument out of range": IndexError,

    "Matrix not in echelon form": ValueError,
    "Matrix not square": ArithmeticError,
    "Incompatible objects": TypeError,

    "Bad syntax, try `-help'": SyntaxError,
    "Bad usage of option, try `-help'": ValueError,
    "Bad number of arguments, try `-help'": ValueError,

    "Not a matrix": TypeError,
    "Not a permutation": TypeError
}

###############################################################
## It is needed to do some initialisation. Since the meataxe
## version used in Sage (SharedMeatAxe) is a dynamic (shared)
## library, it sufficed to do this initialisation only once.
## For convenience, the meataxe_init() function is called in
## this module. Hence, `import sage.libs.meataxe` is enough
## to make sure that MeatAxe is initialised.

from cpython.bytes cimport PyBytes_AsString

cdef void sage_meataxe_error_handler(const MtxErrorRecord_t *err):
    sig_block()
    cdef bytes ErrText = err.Text
    PyErr_SetObject(ErrMsg.get(ErrText.split(': ')[-1], RuntimeError), "{} in file {} (line {})".format(ErrText, err.FileInfo.BaseName, err.LineNo))
    sig_unblock()

cdef inline meataxe_init():
    ## Assign to a variable that enables MeatAxe to find
    ## its multiplication tables.
    import os
    from sage.env import DOT_SAGE
    global MtxLibDir
    MtxLibDir = PyBytes_AsString(os.path.join(DOT_SAGE,'meataxe'))
    ## Error handling for MeatAxe, to prevent immediate exit of the program
    MtxSetErrorHandler(sage_meataxe_error_handler)

meataxe_init()
