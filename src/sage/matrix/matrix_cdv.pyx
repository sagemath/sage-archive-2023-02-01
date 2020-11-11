r"""
Special methods for matrices over discrete valuation rings/fields.
"""

# ****************************************************************************
#       Copyright (C) 2020 Xavier Caruso <xavier.caruso@normalesup.org>
#                          Raphaël Pagès <raphael.pages@u-bordeaux.fr>
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 2 of the License, or
# (at your option) any later version.
#                  https://www.gnu.org/licenses/
# ****************************************************************************


from sage.misc.misc import walltime

from sage.matrix.matrix_generic_dense cimport Matrix_generic_dense
from sage.structure.element cimport RingElement


# We assume that M is square
cpdef charpoly_cdv(Matrix_generic_dense M):
    r"""
    Return the characteristic polynomial of M.

    EXAMPLES::

        sage: ###
    """
    # We compute the Hessenberg form
    # by selecting pivots with minimal valuation on each column
    cdef Py_ssize_t n, m, i, j, k
    R = M.base_ring()
    K = R.fraction_field()
    cdef Matrix_generic_dense H = M.change_ring(K)
    cdef RingElement pivot, inv, scalar

    timeval = 0
    timeop = 0
    timeswap = 0
    timehes = walltime()

    n = H.nrows()
    for j in range(n-1):
        tme = walltime()
        k = j + 1
        mini = H.get_unsafe(k, j).valuation()  # _valuation_c()
        i = j + 2
        while mini != 0 and i < n:
            entry = H.get_unsafe(i, j)
            if entry:
                m = entry.valuation()
                if m < mini:
                    mini = m
                    k = i
            i += 1
        timeval += walltime(tme)

        tme = walltime()
        if k != j + 1:
            H.swap_rows_c(j+1, k)
            H.swap_columns_c(j+1, k)
        timeswap = walltime(tme)
        tme = walltime()
        pivot = H.get_unsafe(j+1, j)
        if pivot:
            inv = ~pivot
            for i in range(j+2, n):
                scalar = inv * H.get_unsafe(i, j)
                H.add_multiple_of_row_c(i, j+1, -scalar, 0)
                H.add_multiple_of_column_c(j+1, i, scalar, 0)
        timeop += walltime(tme)

    timehes = walltime(timehes)

    print("timeval = %s" % timeval)
    print("timeswap = %s" % timeswap)
    print("timeop = %s" % timeop)
    print("timehes = %s" % timehes)

    cdef Matrix_generic_dense c
    c = H.new_matrix(nrows=n+1, ncols=n+1)    # the 0 matrix
    one = H._coerce_element(1)
    c.set_unsafe(0, 0, one)

    for m from 1 <= m <= n:
        # Set the m-th row of c to (x - H[m-1,m-1])*c[m-1] = x*c[m-1] - H[m-1,m-1]*c[m-1]
        # We do this by hand by setting the m-th row to c[m-1]
        # shifted to the right by one.  We then add
        # -H[m-1,m-1]*c[m-1] to the resulting m-th row.
        for i from 1 <= i <= n:
            c.set_unsafe(m, i, c.get_unsafe(m-1,i-1))
        c.add_multiple_of_row_c(m, m-1, -H.get_unsafe(m-1, m-1), 0)
        t = one
        for i from 1 <= i < m:
            t = t * H.get_unsafe(m-i,m-i-1)
            # Set the m-th row of c to c[m] - t*H[m-i-1,m-1]*c[m-i-1]
            c.add_multiple_of_row_c(m, m-i-1, - t*H.get_unsafe(m-i-1,m-1), 0)

    # The answer is now the n-th row of c.
    v = [ c.get_unsafe(n, i) for i in range(n+1) ]

    return R['x'](v)
