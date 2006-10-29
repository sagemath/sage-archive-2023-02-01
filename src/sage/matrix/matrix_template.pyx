###################################################################################
#  Copyright (C) 2006 William Stein                                               #
#  Distributed under the terms of the GNU General Public License (GPL) Version 2. #
#  The full text of the GPL is available at: http://www.gnu.org/licenses/         #
###################################################################################

cimport matrix

cdef class Matrix_generic_dense(matrix.Matrix):
    def __init__(self, parent, entries, coerce, copy):
        matrix.Matrix.__init__(self, parent)
        cdef size_t i, n

        if entries is None:
            entries = 0

        if not isinstance(entries, list):
            try:
                x = parent.base_ring()(entries)
                is_list = False
            except TypeError:
                try:
                    entries = list(entries)
                    is_list = True
                except TypeError:
                    raise TypeError, "entries must be coercible to a list or the basering"

        else:
            is_list = True

        if is_list:

            if len(entries) != self._nrows * self._ncols:
                raise TypeError, "entries has the wrong length"

            if not (coerce or copy):
                self.__entries = entries
            else:
                self.__entries = [None]*(self._nrows*self._ncols)
                n = len(entries)
                if coerce:
                    R = parent.base_ring()
                    for i from 0 <= i < n:
                        self.__entries[i] = R(entries[i])
                else:
                    for i from 0 <= i < n:
                        self.__entries[i] = entries[i]

        else:

            self.__entries = [None]*(self._nrows*self._ncols)
            zero = parent.base_ring()(0)
            for i from 0 <= i < self._nrows * self._ncols:
                self.__entries[i] = zero

            if x != zero:
                if self._nrows != self._ncols:
                    raise TypeError, "nonzero scalar matrix must be square"
                for i from 0 <= i < self._nrows:
                    self.__entries[i*i+i] = x

