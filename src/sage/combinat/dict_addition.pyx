r"""
Linear arithmetic on dictionaries

Provides low-level functions for linear arithmetic of dictionaries
with values in a common ring. Specifically this is used by
:class:`CombinatorialFreeModule`.
"""

#*****************************************************************************
#       Copyright (C) 2010 Christian Stump christian.stump@univie.ac.at
#
#  Distributed under the terms of the GNU General Public License (GPL)
#                  http://www.gnu.org/licenses/
#*****************************************************************************

from cpython cimport PyDict_Copy

cpdef dict dict_addition(dict_iter):
    r"""
    Return the pointwise addition of dictionaries with coefficients.

    INPUT:

    - ``dict_iter`` -- iterator of dictionaries whose values are in
      a common ring and all values are non-zero

    OUTPUT:

    - a dictionary containing all keys of the dictionaries in ``dict_list``
      with values being the sum of the values in the different dictionaries
      (keys with zero value are omitted)

    EXAMPLES::

        sage: from sage.combinat.dict_addition import dict_addition
        sage: D = {0: 1, 1: 1}; D
        {0: 1, 1: 1}
        sage: dict_addition(D for x in range(5))
        {0: 5, 1: 5}

        sage: D1 = {0: 1, 1: 1}; D2 = {0: -1, 1: 1}
        sage: dict_addition([D1, D2])
        {1: 2}
    """
    cdef dict D, D_tmp
    cdef list for_removal

    D = {}
    for D_tmp in dict_iter:
        if D:
            dict_iadd(D, D_tmp, remove_zeros=False)
        else:
            D = PyDict_Copy(D_tmp)

    for_removal = [key for key in D if not D[key]]
    for key in for_removal:
        del D[key]
    return D

cpdef dict_iadd(dict D, dict D2, bint remove_zeros=True, bint negative=False):
    r"""
    Return the inplace pointwise addition of dictionaries ``D`` and ``D2``.

    INPUT:

    - ``D`` -- dictionary that gets mutated whose values are all non-zero
    - ``D2`` -- dictionary whose values are in the same ring as ``D`` and
      are all non-zero
    - ``remove_zeros`` -- boolean; remove the zeros after the addition
      has been performed
    - ``negative`` -- boolean (default: ``False``); add the negative of ``D2``

    OUTPUT:

    None; ``D`` has been mutated.

    EXAMPLES::

        sage: from sage.combinat.dict_addition import dict_iadd
        sage: D1 = {0: 1, 1: 1}
        sage: D2 = {0: -1, 1: 1}
        sage: dict_iadd(D1, D2)
        sage: D1
        {1: 2}

    ::

        sage: D1 = {0: 1, 1: 1}
        sage: D2 = {0: -1, 1: 1}
        sage: dict_iadd(D1, D2, negative=True)
        sage: D1 = {1: 1}

    ::

        sage: D1 = {0: 1, 1: 1}
        sage: D2 = {0: -1, 1: 1}
        sage: dict_iadd(D1, D2, remove_zeros=False, negative=True)
        sage: D1 = {0: 0, 1: 1}
    """
    for key in D2:
        value = D2[key]
        if key in D:
            if negative:
                D[key] -= value
            else:
                D[key] += value
            if remove_zeros and not D[key]:
                del D[key]
        elif value:
            if negative:
                D[key]  = -value
            else:
                D[key]  = value

cpdef dict dict_add(dict D, dict D2, bint negative=False):
    r"""
    Return the pointwise addition of dictionaries ``D`` and ``D2``.

    INPUT:

    - ``D``, ``D2`` -- dictionaries whose values are in a common ring
      and all values are non-zero
    - ``negative`` -- boolean (default: ``False``); add the negative of ``D2``

    EXAMPLES::

        sage: from sage.combinat.dict_addition import dict_add
        sage: D1 = {0: 1, 1: 1}
        sage: D2 = {0: -1, 1: 1}
        sage: dict_add(D1, D2)
        {1: 2}
        sage: D1
        {0: 1, 1: 1}
    """
    # Copy the larger dict and check over the smaller one if adding
    if not negative and len(D) < len(D2):
        D, D2 = D2, D
    cdef dict ret = PyDict_Copy(D)
    dict_iadd(ret, D2, negative=negative)
    return ret

cpdef dict dict_negate(dict D):
    r"""
    Return the negation of the dictionary ``D``.

    EXAMPLES::

        sage: from sage.combinat.dict_addition import dict_negate
        sage: D1 = {0: 1, 1: 1}
        sage: dict_negate(D1)
        {0: -1, 1: -1}
    """
    return {key: -D[key] for key in D}

cpdef dict dict_linear_combination(dict_factor_iter, bint factor_on_left=True):
    r"""
    Return the pointwise addition of dictionaries with coefficients.

    INPUT:

    - ``dict_factor_iter`` -- iterator of pairs ``D``, ``coeff``, where

      * the ``D``'s are dictionaries with values in a common ring
      * the ``coeff``'s are coefficients in this ring

    - ``factor_on_left`` -- boolean (default: ``True``); if ``True``,
      the coefficients are multiplied on the left, otherwise they are
      multiplied on the right

    OUTPUT:

    - a dictionary containing all keys of dictionaries in ``dict_list``
      with values being the sum of the values in the different
      dictionaries and each one first multiplied by the given factor
      (keys with zero value are omitted)

    EXAMPLES::

        sage: from sage.combinat.dict_addition import dict_linear_combination
        sage: D = { 0:1, 1:1 }; D
        {0: 1, 1: 1}
        sage: dict_linear_combination( (D,i) for i in range(5) )
        {0: 10, 1: 10}
        sage: dict_linear_combination( [(D,1), (D,-1)] )
        {}
    """
    cdef dict D = {}
    cdef dict D_tmp
    cdef list for_removal

    for D_tmp, fac_tmp in dict_factor_iter:
        if not fac_tmp: # We multiply by 0, so nothing to do
            continue
        if not D and fac_tmp == 1:
            D = PyDict_Copy(D_tmp)
        elif fac_tmp == 1:
            dict_iadd(D, D_tmp, remove_zeros=False)
        elif fac_tmp == -1:
            dict_iadd(D, D_tmp, remove_zeros=False, negative=True)
        else:
            if factor_on_left:
                for key in D_tmp:
                    value = D_tmp[key]
                    if key in D:
                        D[key] += fac_tmp * value
                    else:
                        D[key]  = fac_tmp * value
            else:
                for key in D_tmp:
                    value = D_tmp[key]
                    if key in D:
                        D[key] += value * fac_tmp
                    else:
                        D[key]  = value * fac_tmp

    for_removal = [key for key in D if not D[key]]
    for key in for_removal:
        del D[key]

    return D

