r"""
The cythonized low-level methods used in the implementation of the universal cyclotomic field with the Zumbroich basis.

.. SEEALSO::

    :class:`UniversalCyclotomicField`

AUTHORS:

- Christian Stump
"""
#*****************************************************************************
#       Copyright (C) 2012 Christian Stump <christian.stump@univie.ac.at>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#                  http://www.gnu.org/licenses/
#*****************************************************************************
include "sage/ext/stdsage.pxi"
include "sage/ext/cdefs.pxi"
include "sage/ext/gmp.pxi"

import sys
import operator

from cpython cimport bool, PyDict_Copy
from itertools import product
from sage.combinat.dict_addition import dict_linear_combination, dict_addition
from sage.categories.map cimport Map
from sage.rings.arith import prod, factor, lcm
from sage.rings.finite_rings.integer_mod cimport mod_inverse_int
from sage.rings.integer import LCM_list
from sage.rings.all import QQ
from sage.rings.rational cimport Rational

import numpy as np
cimport numpy as np

import cython

DTYPE = np.int
ctypedef np.int_t DTYPE_t

cdef dict factor_cache = dict()
cpdef cached_factor(int n):
    r"""
    Returns a cached version of the dict of a factorization.

    :param n: a positive integer

    EXAMPLES::

        sage: from sage.rings.universal_cyclotomic_field.universal_cyclotomic_field_c import cached_factor
        sage: cached_factor(12)
        {2: 2, 3: 1}

    .. WARNING::

        Internal function, not to be used directly!
    """
    if n in factor_cache:
        return factor_cache[n]
    else:
        X = dict(factor(n))
        factor_cache[n] = X
        return X

cpdef ZumbroichIndexSet(int p, int k):
    """
    Returns a list of integers depending on p and k, which is needed in the :meth:`ZumbroichBasisCython`.

    :param p: a positive integer
    :param k: a positive integer

    OUTPUT:

    - list of indices as defined for J_{k,p} [B97] p. 283, Remark 1

    EXAMPLES::

        sage: from sage.rings.universal_cyclotomic_field.universal_cyclotomic_field_c import ZumbroichIndexSet
        sage: ZumbroichIndexSet(6,0)
        [1, 2, 3, 4, 5]
        sage: ZumbroichIndexSet(6,1)
        [-2, -1, 0, 1, 2]
    """
    if k == 0:
        if p == 2:
            return [ 0 ]
        else:
            return range(1,p)
    else:
        if p == 2:
            return [0,1]
        else:
            return range(-(p-1)/2, (p-1)/2+1)

cdef dict zumbroich_basis_cache = dict()
cpdef ZumbroichBasisCython(int n, int m = 1):
    """
    Returns the Zumbroich basis of the cyclotomics of order `n` over the cyclotomics of order `m`.

    :param n: positive integer
    :param m: positive integer dividing ``n``
    :type m: optional, default:``1``

    OUTPUT:

    - '\{i_1,\ldots,i_k\}` such that `\{ E(n)^{i_j} \}` is the Zumbroich basis of the cyclotomic field
      of order `n` over the cyclotomic field of order `m`.

    .. NOTE::

        The computation is based on the description in [Bre97]_.

    EXAMPLES:

    We compute the Zumbroich basis of the cyclotomics of order `12` over the cyclotomics of order `3`, given by `\{E(12)^0,E(12)^3\}`::

        sage: from sage.rings.universal_cyclotomic_field.universal_cyclotomic_field_c import ZumbroichBasisCython
        sage: ZumbroichBasisCython(12,3)
        array([0, 3])
    """
    if (n,m) in zumbroich_basis_cache:
        return zumbroich_basis_cache[(n,m)]
    cdef int mu, n_fac, m_fac, k, p, i
    cdef list list_of_index_sets
    cdef tuple index_tuple
    cdef dict n_fac_dict, quo_fac_dict

    n_fac_dict = cached_factor(n)
    quo_fac_dict = cached_factor(n/m)

    list_of_index_sets = []

    for p in quo_fac_dict:
        mu = quo_fac_dict[p]
        n_fac = n_fac_dict[ p ]
        m_fac = n_fac - mu
        for k in xrange(m_fac, n_fac):
            list_of_index_sets.append([ n*i / p**(k+1) for i in ZumbroichIndexSet(p, k) ])
    index_tuples = product(*list_of_index_sets)
    cdef int size = prod([ len(index_set) for index_set in list_of_index_sets ])

    cdef np.ndarray[DTYPE_t,ndim=1] L = np.ndarray([size],dtype=DTYPE)
    cdef tuple X
    for i from 0 <= i < size:
        X = index_tuples.next()
        L[i] = sum(X)%n
    L.sort()
    zumbroich_basis_cache[(n,m)] = L
    return L

cdef inline int find_array(int a, int * ar, int size):
    r"""
    Returns the position of a in ar.

    .. WARNING::

        Internal function, not to be used directly!
    """
    cdef unsigned int i
    for i from 0 <= i < size:
        if a == ar[i]:
            return i
    return size

cpdef ZumbroichDecomposition(int n, int i):
    """
    Returns the decomposition of `E(n,i)` in the Zumbroich basis of the cyclotomics of order `n`.

    :param n: positive integer
    :param i: positive integer `<` ``n``

    OUTPUT:

    - The decomposition of `E(n,i)` in the Zumbroich basis, as a dictionary

    EXAMPLES::

        sage: from sage.rings.universal_cyclotomic_field.universal_cyclotomic_field_c import ZumbroichDecomposition
        sage: ZumbroichDecomposition(6, 1)
        {(6, 4): -1}
    """
    cdef int mu, p, k, j
    cdef dict n_fac, i_mod_n_dict, return_dict
    cdef list mod_prime, list_of_i_mod_n_dict, new_list_of_i_mod_n_dict

    n_fac = cached_factor(n)
    mod_prime = mod_prime_list(n,i,n_fac)

    cdef int c = -1 if len(mod_prime) %2 else 1

    prime_index_list = [ range(1,p) for p in mod_prime ]
    prime_index_tuples = product(*prime_index_list)

    return_dict = {}

    for tup in prime_index_tuples:
        k = (i + sum([ n * tup[j] / mod_prime[j] for j in
                       xrange(len(mod_prime))] )) % n
        return_dict[ (n,k) ] = c

    return return_dict

cpdef list ZumbroichDecomposition_list(int n, int i):
    """
    Returns the decomposition of `E(n,i)` in the Zumbroich basis of the cyclotomics of order `n`.

    :param n: positive integer
    :param i: positive integer `<` ``n``

    OUTPUT:

    - The decomposition of `E(n,i)` in the Zumbroich basis, as an array of indices of the Zumbroich basis, together with the length of mod_prime_list at the end.
        It is used as follows: `\zeta_n^i = (-1)^l \sum \zeta_n^ZB(k)`, where l is the last entry of the last of the array, ZB(k) is the `k`-th entry in the Zumbroich
        basis and where the sum ranges over all but the last entry in the array.


    EXAMPLES::

        sage: from sage.rings.universal_cyclotomic_field.universal_cyclotomic_field_c import ZumbroichDecomposition_list
        sage: from sage.rings.universal_cyclotomic_field.universal_cyclotomic_field_c import ZumbroichBasisCython
        sage: ZumbroichDecomposition_list(6, 1)
        [1, 1]

    This reads as `\zeta_6 = (-1)^1 \zeta_6^ZB(1) = -\zeta_6^4`
    """
    cdef dict n_fac = cached_factor(n)
    cdef list mod_prime = mod_prime_list(n,i,n_fac)
    cdef int list_size = len(mod_prime)

    cdef np.ndarray[DTYPE_t] ZB_list
    cdef int ZB_list_size
    ZB_list = ZumbroichBasisCython(n)
    ZB_size = len(ZB_list)
    cdef int *ZB_array = <int*>sage_malloc(sizeof(int) * ZB_size)
    cdef int j
    for j from 0 <= j < ZB_size:
        ZB_array[j] = ZB_list[j]

    cdef int p
    cdef list prime_index_list = [ range(1,p) for p in mod_prime ]
    prime_index_tuples = product(*prime_index_list)

    cdef int prime_prod = 1
    for j from 0 <= j < list_size:
        p = mod_prime[j]
        prime_prod *= p - 1
        mod_prime[j] = n / p
    cdef list return_list = []
    cdef int tup_sum, fac, k
    cdef tuple tup
    for tup in prime_index_tuples:
        tup_sum = i
        for j from 0 <= j < list_size:
            p = mod_prime[j]
            fac = tup[j]
            tup_sum += p * fac
        k = tup_sum % n
        if k < 0:
            k += n
        return_list.append(find_array(k, ZB_array, ZB_size))
    return_list.append(list_size)

    return return_list

cpdef mod_prime_list(int n, int i, dict n_fac):
    r"""
    Returns the list of prime factors of ``n`` which need to be modified in order to express `\zeta_n^i` in the Zumbroich basis.

    For details see [Bre97]_.

    :param n: positive integer
    :param i: positive integer `<` ``n``

    EXAMPLES::

        sage: from sage.rings.universal_cyclotomic_field.universal_cyclotomic_field_c import mod_prime_list
        sage: n = 12
        sage: n_fac = dict(n.factor())
        sage: for k in [1..12]: print k, mod_prime_list(n,k,n_fac)
        1 [2]
        2 [2]
        3 [3]
        4 []
        5 [2]
        6 [2, 3]
        7 []
        8 []
        9 [2, 3]
        10 [2]
        11 []
        12 [3]
    """
    cdef list mod_primes_list

    cdef Py_ssize_t factor, count,mu,k,p,m

    mod_primes_list = []
    cdef  Py_ssize_t max_mu = 0

    for p in n_fac:
        mu = n_fac[p]
        if mu > max_mu:
            max_mu = mu
    cdef int *tmp_list = <int *>sage_malloc(max_mu * sizeof(int))
    if not tmp_list:
            raise MemoryError()

    for p in n_fac:
        mu = n_fac[p]
        factor = p**mu

        m = n / factor

        count = (i * mod_inverse_int(m,factor)) % factor
        if count < 0:
            count += factor
        i -= count * m

        if p == 2:
            if 1 == count/(factor/p):
                mod_primes_list.append(p)
            factor /= p**mu
        else:
            for k from 0 <= k < mu:
                factor /= p
                tmp_list[k] = count / factor
                count = count % factor

            for k from mu-1 >= k > 0:
                if tmp_list[k] > (p-1) / 2:
                    tmp_list[k] -= p
                    if k == 1 and tmp_list[k-1] == p-1:
                        tmp_list[k-1] = 0
                    else:
                        tmp_list[k-1] += 1
            if tmp_list[0] == 0:
                mod_primes_list.append(p)
    sage_free(tmp_list)
    return mod_primes_list

cdef res_classes(list L, int k):
    """
    Groups the elements of ``L`` by residue class modulo ``k``.

    :param L: a list of integers
    :param k: positive integer

    OUTPUT:

    - a dictionary mapping residue classes to lists of elements of ``L``.

    .. WARNING::

        Internal function, not to be used directly!
    """
    cdef int x, tmp
    cdef dict res

    res = {}
    for i in range(len(L)):
        x = L[i]
        tmp = x % k
        if tmp in res:
            res[ tmp ].append(x)
        else:
            res[ tmp ] = [ x ]
    return res

cpdef push_to_higher_field(dict D, int n, int l):
    """
    Pushes an element of the universal cyclotomic field of order `n`, represented by its dictionary,
    to the cyclotomic field of order l. This method is used e.g. for adding/multiplying two elements in different cyclotomics.

    :param D: a dictionary representing an element in the universal cyclotomic field
    :param n: the order of the cyclotomics containing ``D``
    :param l: a multiple of ``n``, the order of the field ``D`` is pushed

    EXAMPLES::

        sage: from sage.rings.universal_cyclotomic_field.universal_cyclotomic_field_c import push_to_higher_field
        sage: UCF.<E> = UniversalCyclotomicField()

        sage: D = E(6)._dict_(); D
        {(3, 2): -1}

        sage: push_to_higher_field(D,3,6)
        {(6, 4): -1}
    """
    cdef int k
    if n == l:
        return D
    else:
        return dict_linear_combination([ (ZumbroichDecomposition(l, k*l/n), D[ (n,k) ]) for _,k in D ])

cpdef zip_key_dict(list keys, dict key_dict, int n, int m, int p, bool positive=True):
    """
    TESTS::

        sage: from sage.rings.universal_cyclotomic_field.universal_cyclotomic_field_c import push_down_cython
        sage: D = { (6,4):-1 }
        sage: push_down_cython(6,D) # indirect doctest
        {(3, 2): -1}

    .. WARNING::

        Internal function, not to be used directly!
    """
    cdef int fac
    cdef int k
    cdef dict G
    cdef tuple key

    G = {}
    if positive:
        fac = 1
    else:
        fac = -1

    for key in keys:
        G[ key ] = fac * key_dict[ (n, (m+key[1]*p)%n) ]
    return G

cpdef push_down_cython(int n, dict dict_in_basis):
    """
    Returns `x` in the smallest cyclotomic field containing it, as a dictionary.
    Implementation of the push down algorithm described in [B97].

    :param n: a positive integer
    :param dict_in_basis: a dictionary representing an element of
       the cyclotomic field of order `n` expressed in the Zumbroich basis

    EXAMPLES::

        sage: from sage.rings.universal_cyclotomic_field.universal_cyclotomic_field_c import push_down_cython
        sage: D = { (6,4):-1 }
        sage: push_down_cython(6,D)
        {(3, 2): -1}
    """
    cdef int mu
    cdef int p, m, k, _
    cdef dict res
    cdef bool has_changed, should_change

    n_fac = cached_factor(n)

    for p in n_fac:
        mu = n_fac[p]
        has_changed = True
        while has_changed and mu > 0:
            has_changed = False
            m = n/p

            should_change = True
            if mu > 1:
                for _,k in dict_in_basis:
                    if k % p != 0:
                        should_change = False
                if should_change:
                    dict_in_basis = dict([ ((m,k/p), dict_in_basis[ (n,k) ]) for n,k in dict_in_basis ])
                    has_changed = True
            else:
                res = res_classes([ L[1] for L in dict_in_basis ], m)
                for k in res:
                    if len(res[ k ]) != p - 1:
                        should_change = False
                if should_change:
                    for k in res:
                        res[ k ] = k
                        while res[ k ] % p != 0:
                            res[ k ] += m
                        res[ k ] /= p
                    keys = [ (m, res[k]) for k in res ]

                    if p == 2:
                        dict_in_basis = zip_key_dict(keys, dict_in_basis, n, 0, p, positive=True)
                        has_changed = True
                    else:
                        for _,k in keys:
                            if should_change:
                                for i in range(2, p):
                                    if dict_in_basis[ (n, (m*i+k*p)%n) ] != dict_in_basis[ (n,(m+k*p)%n) ]:
                                        should_change = False
                        if should_change:
                            dict_in_basis = zip_key_dict(keys, dict_in_basis, n, m, p, positive=False)
                            has_changed = True
            if has_changed:
                n = m
                mu -= 1
    return dict_in_basis

cpdef galois_conjugates_cython(dict D, int n, int m, list coprimes):
    """
    Returns all Galois conjugates of ``elem`` which lives in the cyclotomics of order ``n``, and with the ``coprimes`` already given.
    Observe that the product of all Galois conjugates is 1, and thus, this method is used to compute the inverse for elements in
    the universal cyclotomic field.

    :param D: a dictionary representing an element in the universal cyclotomic field
    :param n: order of the element in ``D``
    :param m: a multiple of ``n`` in which field the Galois conjugates should be considered
    :param coprimes: list of coprimes of ``n``

    EXAMPLES::

        sage: from sage.rings.universal_cyclotomic_field.universal_cyclotomic_field_c import galois_conjugates_cython
        sage: UCF.<E> = UniversalCyclotomicField()

        sage: galois_conjugates_cython(E(12).value._monomial_coefficients, 12, 12, [1,5,7,11])
        [{(12, 7): -1}, {(12, 11): -1}, {(12, 7): 1}, {(12, 11): 1}]

        sage: galois_conjugates_cython(E(12).value._monomial_coefficients, 12, 24, [1, 5, 7, 11, 13, 17, 19, 23])
        [{(12, 7): -1}, {(12, 11): -1}, {(12, 7): 1}, {(12, 11): 1}, {(12, 7): -1}, {(12, 11): -1}, {(12, 7): 1}, {(12, 11): 1}]

        sage: galois_conjugates_cython(E(9).value._monomial_coefficients, 9, 9, [1, 2, 4, 5, 7, 8])
        [{(9, 4): -1, (9, 7): -1}, {(9, 2): 1}, {(9, 4): 1}, {(9, 5): 1}, {(9, 7): 1}, {(9, 2): -1, (9, 5): -1}]
    """
    cdef int i, k
    cdef list conjugates

    if n != m:
        D = push_to_higher_field(D, n, m)
    if len(coprimes) == m-1 and m > 2:
        conjugates = []
        for i in coprimes:
            conjugates.append(dict([ ((m,(k*i) % m), D[(m,k)]) for m,k in D]))
    else:
        conjugates = []
        for i in coprimes:
            conjugates.append(dict_linear_combination([(ZumbroichDecomposition(m, (k*i) % m), D[(m,k)]) for m,k in D]))
    if n != m:
        conjugates = [ push_down_cython(m,D) for D in conjugates ]
    return conjugates

cpdef dict_vector_multiplication(int n, tuple l1, tuple l2):
    r"""
    Returns the vector multiplication of ``l1`` and ``l2``.

    :param l1, l2: tuples of the same length containing dictionaries representing elements in the universal cyclotomics

    EXAMPLES::

        sage: from sage.rings.universal_cyclotomic_field.universal_cyclotomic_field_c import dict_vector_multiplication
        sage: UCF.<E> = UniversalCyclotomicField()

        sage: D1 = E(3).value._monomial_coefficients; D1
        {(3, 1): 1}

        sage: D2 = E(4).value._monomial_coefficients; D2
        {(4, 1): 1}

        sage: l = (D1,D2)
        sage: dict_vector_multiplication(12,l,l)
        {(12, 4): 1, (12, 8): 2}
    """
    cdef int a0,a1,b0,b1,i,k
    cdef list out = []
    cdef dict empty_dict = {}
    cdef dict D1, D2

    for i in xrange(len(l1)):
        D1 = l1[i]
        D2 = l2[i]
        if D1 == empty_dict or D2 == empty_dict:
            out.append(empty_dict)
        else:
            n1 = iter(D1).next()[0]
            n2 = iter(D2).next()[0]
            out.append(dict_multiplication(D1,D2,n1,n2,n))
    return dict_addition(out)

cdef int *linear_combination_pointer
cdef add_to_basis_linear_combination(list decomposition, int coeff, int size):
    r"""
    Adds the entries in ``decomposition`` `\times` ``coeff`` to ``linear_combination_pointer``.

    .. WARNING::

        Internal function, not to be used directly!
    """
    global linear_combination_pointer
    cdef int value, add_value,key

    cdef int counter = 0
    for key in decomposition:
        if counter < size:
            add_value = linear_combination_pointer[key]
            add_value += coeff
            linear_combination_pointer[key] = add_value
            counter += 1

cdef inline as_rat(int a):
    r"""
    Returns the integer ``a`` as a rational.

    .. WARNING::

        Internal function, not to be used directly!
    """
    cdef Rational rat = <Rational> PY_NEW(Rational)
    mpq_set_si(rat.value,a,1)
    return rat

cpdef dict dict_multiplication(dict D1, dict D2, int n1, int n2, int n):
    r"""
    Returns the multiplication of two elements in the universal cyclotomics, both represented as a dictionary, again as a dict.

    EXAMPLES::

        sage: from sage.rings.universal_cyclotomic_field.universal_cyclotomic_field_c import dict_multiplication
        sage: UCF.<E> = UniversalCyclotomicField()

        sage: D1 = E(3).value._monomial_coefficients; D1
        {(3, 1): 1}

        sage: D2 = E(4).value._monomial_coefficients; D2
        {(4, 1): 1}

        sage: dict_multiplication(D1,D2,3,4,12)
        {(12, 7): 1}
    """
    # the overhead of the method seems to be big, but most of the part is spent in the double for loop, and I don't know how
    # to improve this piece for now

    # setting L1, L2, and Lnew to the indices of the Zumbroich basis for n1, n2, and n (which are supposed to be the order of the two dicts and their lcm)
    # also setting L1_size, L2_size, and Lnew_size to their sizes
    cdef np.ndarray[DTYPE_t] L1, L2, Lnew
    cdef int L1_size,L2_size,Lnew_size
    L1 = ZumbroichBasisCython(n1)
    L1_size = len(L1)
    L2 = ZumbroichBasisCython(n2)
    L2_size = len(L2)
    Lnew = ZumbroichBasisCython(n)
    Lnew_size = len(Lnew)

    # writing L1 and L2 into c arrays
    cdef int *L1_array = <int*>sage_malloc(sizeof(int) * L1_size)
    cdef int *L2_array = <int*>sage_malloc(sizeof(int) * L2_size)
    cdef int j
    for j from 0 <= j < L1_size:
        L1_array[j] = L1[j]
    for j from 0 <= j < L2_size:
        L2_array[j] = L2[j]

    # setting up pointers for the two input dicts
    cdef Py_ssize_t key

    cdef int *D1_pointer = <int *>sage_malloc(L1_size * sizeof(int))
    for key from 0 <= key < L1_size:
        D1_pointer[key] = 0
    cdef int *D2_pointer = <int *>sage_malloc(L2_size * sizeof(int))
    for key from 0 <= key < L2_size:
        D2_pointer[key] = 0

    # setting up pointer for the output dict
    global linear_combination_pointer
    linear_combination_pointer = <int *>sage_malloc(Lnew_size * sizeof(int))
    for key from 0 <= key < Lnew_size:
        linear_combination_pointer[key] = 0

    #computing the lcm of the denominators in order to work with ints rather than with rationals
    #FIXME: there might be a faster method doing this
    cdef list LCM = [None]*len(D1)
    key = 0
    for X in D1:
        X = D1[X]
        if type(X) == Rational:
            LCM[key] = X.denom()
        else:
            LCM[key] = X
        key += 1
    cdef int denom1 = LCM_list(LCM)
    cdef Rational denom_rat1 = as_rat(denom1)
    LCM = [None]*len(D2)
    key = 0
    for X in D2:
        X = D2[X]
        if type(X) == Rational:
            LCM[key] = X.denom()
        else:
            LCM[key] = X
        key += 1
    cdef int denom2 = LCM_list(LCM)
    cdef Rational denom_rat2 = as_rat(denom2)

    # filling the pointers for the dicts with ints from the dicts multiplied with the lcm
    # we use positions in the Zumbroich basis with find_array
    cdef Rational old_val, new_val1, new_val2
    new_val1 = <Rational> PY_NEW(Rational)
    new_val2 = <Rational> PY_NEW(Rational)

    cdef int k

    for X in D1:
        k = find_array(X[1],L1_array,L1_size)
        old_val = D1[X]
        mpq_mul(new_val1.value,old_val.value,denom_rat1.value)
        D1_pointer[k] = new_val1
    for X in D2:
        k = find_array(X[1],L2_array,L2_size)
        old_val = D2[X]
        mpq_mul(new_val2.value,old_val.value,denom_rat2.value)
        D2_pointer[k] = new_val2

    # finally, everything is set up such that we can do the actual multiplication
    cdef unsigned int size
    cdef int fac
    cdef int fac1, fac2

    # these pieces are used to cache ZumbroichDecomposition_list properly,
    # in the c array, the k's are stored, and in the python list the actual lists
    cdef int *zumbroich_decomposition_chache_index = <int *>sage_malloc(n * sizeof(int))
    cdef list zumbroich_decomposition_chache_values = []
    cdef int zumbroich_decomposition_chache_size = 0
    cdef int index

    cdef list L

    cdef Py_ssize_t counter1, couter2
    n1 = n / n1
    n2 = n / n2
    for counter1 from 0 <= counter1 < L1_size:
        fac1 = D1_pointer[counter1]
        if fac1 != 0:
            for counter2 from 0 <= counter2 < L2_size:
                fac2 = D2_pointer[counter2]
                if fac2 != 0:
                    k = ((n1*L1[counter1])+(n2*L2[counter2])) % n
                    if k < 0:
                        k += n
                    # this is the actual caching part
                    index = find_array(k,zumbroich_decomposition_chache_index,zumbroich_decomposition_chache_size)
                    if index < zumbroich_decomposition_chache_size:
                        L = zumbroich_decomposition_chache_values[index]
                    else:
                        L = ZumbroichDecomposition_list(n,k)
                        zumbroich_decomposition_chache_index[zumbroich_decomposition_chache_size] = k
                        zumbroich_decomposition_chache_values.append(L)
                        zumbroich_decomposition_chache_size += 1
                    fac = fac1 * fac2
                    size = len(L)
                    k = L[size-1]
                    if k%2 == 1:
                        fac = -fac
                    # using a different method to add the entries in L to the c array linear_combination_pointer
                    add_to_basis_linear_combination(L,fac,size-1)

    # in the end, the data is stored in a dictionary
    cdef dict D = {}
    cdef Rational new_value
    for key from 0 <= key < Lnew_size:
        if linear_combination_pointer[key] != 0:
            old_val = as_rat(linear_combination_pointer[key])
            new_val = <Rational> PY_NEW(Rational)
            mpq_div(new_val.value,old_val.value,denom_rat1.value)
            mpq_div(new_val.value,new_val.value,denom_rat2.value)
            D[(n,Lnew[key])] = new_val

    sage_free(linear_combination_pointer)
    sage_free(D1_pointer)
    sage_free(D2_pointer)
    return D
