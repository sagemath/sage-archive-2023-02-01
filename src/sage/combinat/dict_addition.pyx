r"""
Pointwise addition of dictionaries

Provides function to add dictionaries pointwise with values in a common ring and to compute linear combinations

EXAMPLES::

    sage: from sage.combinat.dict_addition import dict_addition
    sage: D1 = { 0:1, 1:1 }; D2 = { 0:-1, 1:1 }
    sage: dict_addition( [D1,D2] )
    {1: 2}
"""
#*****************************************************************************
#       Copyright (C) 2010 Christian Stump christian.stump@univie.ac.at
#
#  Distributed under the terms of the GNU General Public License (GPL)
#                  http://www.gnu.org/licenses/
#*****************************************************************************
from cpython cimport PyDict_Copy

cpdef dict_addition(dict_iter):
    r"""
    Returns the pointwise addition of dictionaries with coefficients.

    :param dict_iter: iterator of dictionaries with values in a common ring.

    OUTPUT:

    - a dictionary containing all keys of dictionaries in ``dict_list``, with values being the sum of the values in the different dictionaries (keys with zero value are omitted)

    EXAMPLES::

        sage: from sage.combinat.dict_addition import dict_addition
        sage: D = { 0:1, 1:1 }; D
        {0: 1, 1: 1}
        sage: dict_addition( D for _ in range(5) )
        {0: 5, 1: 5}

        sage: D1 = { 0:1, 1:1 }; D2 = { 0:-1, 1:1 }
        sage: dict_addition( [D1,D2] )
        {1: 2}
    """
    cdef dict D, D_tmp

    D = {}
    for D_tmp in dict_iter:
        if D == {}:
            D = PyDict_Copy( D_tmp )
        else:
            for key in D_tmp:
                value = D_tmp[key]
                if key in D:
                    D[key] += value
                else:
                    D[key]  = value
    for_removal = [key for key in D if not D[key]]
    for key in for_removal:
        del D[key]
    return D

cpdef dict_linear_combination( dict_factor_iter, factor_on_left=True ):
    r"""
    Returns the pointwise addition of dictionaries with coefficients.

    :param dict_factor_iter: iterator of pairs D, coeff, where
        - the D's are dictionaries with values in a common ring
        - the coeff's are coefficients in this ring

    :param factor_on_left: if True, the coefficients are multiplied on the left, otherwise they are multiplied on the right

    :type factor_on_left: boolean; optional, default ``True``

    OUTPUT:

    - a dictionary containing all keys of dictionaries in ``dict_list``, with values being the sum of the values in the different dictionaries, each one first multiplied by the given factor (keys with zero value are omitted)

    EXAMPLES::

        sage: from sage.combinat.dict_addition import dict_linear_combination
        sage: D = { 0:1, 1:1 }; D
        {0: 1, 1: 1}
        sage: dict_linear_combination( (D,i) for i in range(5) )
        {0: 10, 1: 10}
        sage: dict_linear_combination( [(D,1),(D,-1)] )
        {}
    """
    D = {}
    for D_tmp, fac_tmp in dict_factor_iter:
        if D == {} and fac_tmp == 1:
            D = PyDict_Copy(D_tmp)
        elif fac_tmp == 1:
            for key in D_tmp:
                value = D_tmp[key]
                if key in D:
                    D[ key ] += value
                else:
                    D[ key ]  = value
        elif fac_tmp == -1:
            for key in D_tmp:
                value = D_tmp[key]
                if key in D:
                    D[ key ] -= value
                else:
                    D[ key ]  = -value
        else:
            if factor_on_left:
                for key in D_tmp:
                    value = D_tmp[key]
                    if key in D:
                        D[ key ] += fac_tmp * value
                    else:
                        D[ key ]  = fac_tmp * value
            else:
                for key in D_tmp:
                    value = D_tmp[key]
                    if key in D:
                        D[ key ] += value * fac_tmp
                    else:
                        D[ key ]  = value * fac_tmp

    for_removal = [key for key in D if not D[key]]
    for key in for_removal:
        del D[key]

    return D
