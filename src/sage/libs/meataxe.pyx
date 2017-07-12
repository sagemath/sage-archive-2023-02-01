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
    "Bad file format": IOError,
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

cdef void sage_meataxe_error_handler(const MtxErrorRecord_t *err):
    sig_block()
    PyErr_SetObject(ErrMsg.get(err.Text, SystemError), "{} in file {} (line {})".format(err.Text, err.FileInfo.BaseName, err.LineNo))
    sig_unblock()
