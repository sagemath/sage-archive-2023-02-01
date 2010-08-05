"""
Provides function to add dictionaries pointwise with values in a common ring and to compute linear combinations
"""
#*****************************************************************************
#       Copyright (C) 2010 Christian Stump christian.stump@univie.ac.at
#
#  Distributed under the terms of the GNU General Public License (GPL)
#                  http://www.gnu.org/licenses/
#*****************************************************************************

def dict_addition( dict_iter ):
    """
    returns the pointwise addition of dictionaries with coefficients

    INPUT:
        dict_iter         -- iterator of dictionaries with values in a common ring R

    OUTPUT:
        dictionary containing all keys of dictionaries in dict_list (and non-zero values) being the the sum of the values in the different
        dictionaries

    EXAMPLES:
        sage: from sage.combinat.dict_addition import dict_addition
        sage: D = { 0:1, 1:1 }; D
        {0: 1, 1: 1}
        sage: dict_addition( D for _ in range(5) )
        {0: 5, 1: 5}
    """
    cdef dict D, D_tmp

    D = {}
    for D_tmp in dict_iter:
        if D == {}:
            D = D_tmp.copy()
        else:
            for key, value in D_tmp.iteritems():
                if key in D:
                    D[ key ] += value
                else:
                    D[ key ]  = value

    for_removal = [ key for key, value in D.iteritems() if value == 0 ]
    for key in for_removal:
        del D[key]

    return D

def dict_linear_combination( dict_factor_iter, factor_on_left=True ):
    """
    returns the pointwise addition of dictionaries with coefficients

    INPUT:
        dict_factor_iter            -- iterator of pairs D, coeff, where
                                        the D's are dictionaries with values in a common ring R
                                        the coeff's are coefficients in R
        factor_on_left(optional)    -- if True, the coefficients are multiplied on the left, otherwise they are multiplied on the right

    OUTPUT:
        dictionary containing all keys of dictionaries in dict_list (and non-zero values) being the the sum of the values in the different
        dictionaries each one first multiplied by the given factor

    EXAMPLES:
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
            D = D_tmp.copy()
        elif fac_tmp == 1:
            for key, value in D_tmp.iteritems():
                if key in D:
                    D[ key ] += value
                else:
                    D[ key ]  = value
        elif fac_tmp == -1:
            for key, value in D_tmp.iteritems():
                if key in D:
                    D[ key ]  -= value
                else:
                    D[ key ]   = -value
        else:
            if factor_on_left:
                for key, value in D_tmp.iteritems():
                    if key in D:
                        D[ key ] += fac_tmp * value
                    else:
                        D[ key ]  = fac_tmp * value
            else:
                for key, value in D_tmp.iteritems():
                    if key in D:
                        D[ key ] += value * fac_tmp
                    else:
                        D[ key ]  = value * fac_tmp

    for_removal = [ key for key, value in D.iteritems() if value == 0 ]
    for key in for_removal:
        del D[key]

    return D
