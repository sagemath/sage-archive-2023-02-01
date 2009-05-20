"""
Partition backtrack functions for lists -- a simple example of using partn_ref.

DOCTEST:
    sage: import sage.groups.perm_gps.partn_ref.refinement_lists

"""

#*****************************************************************************
#      Copyright (C) 2006 - 2008 Robert L. Miller <rlmillster@gmail.com>
#      Copyright (C) 2009 Nicolas Borie <nicolas.borie@math.u-psud.fr>
#
# Distributed  under  the  terms  of  the  GNU  General  Public  License (GPL)
#                         http://www.gnu.org/licenses/
#*****************************************************************************

include 'data_structures_pyx.pxi' # includes bitsets

def is_isomorphic(self, other):
    r"""
    Return the bijection as a permutation if two lists are isomorphic, return
    False otherwise.

    EXAMPLE::

        sage: from sage.groups.perm_gps.partn_ref.refinement_lists import is_isomorphic
        sage: is_isomorphic([0,0,1],[1,0,0])
        [1, 2, 0]

    """
    cdef int **part, i, j
    cdef int *output, *ordering
    partition = [range(len(self))]
    part = <int **> sage_malloc((len(partition)+1) * sizeof(int *))
    ordering = <int *> sage_malloc((len(self)) * sizeof(int))
    if part is NULL or ordering is NULL:
        if part is not NULL: sage_free(part)
        if ordering is not NULL: sage_free(ordering)
        raise MemoryError
    for i from 0 <= i < len(partition):
        part[i] = <int *> sage_malloc((len(partition[i])+1) * sizeof(int))
        if part[i] is NULL:
            for j from 0 <= j < i:
                sage_free(part[j])
            sage_free(part)
            raise MemoryError
        for j from 0 <= j < len(partition[i]):
            part[i][j] = partition[i][j]
        part[i][len(partition[i])] = -1
    part[len(partition)] = NULL
    for i from 0 <= i < (len(self)):
        ordering[i] = i

    output = double_coset(self, other, part, ordering, (len(self)), &all_list_children_are_equivalent, &refine_list, &compare_lists)

    for i from 0 <= i < len(partition):
        sage_free(part[i])
    sage_free(part)
    sage_free(ordering)

    if output is NULL:
        return False
    else:
        output_py = [output[i] for i from 0 <= i < (len(self))]
        sage_free(output)
        return output_py

cdef bint all_list_children_are_equivalent(PartitionStack *PS, object S):
    return 0

cdef int refine_list(PartitionStack *PS, object S, int *cells_to_refine_by, int ctrb_len):
    return 0

cdef int compare_lists(int *gamma_1, int *gamma_2, object S1, object S2):
    r"""
    Compare two lists according to the lexicographic order.
    """
    cdef list MS1 = <list> S1
    cdef list MS2 = <list> S2
    cdef int i, j
    for i from 0 <= i < len(MS1):
        j = cmp(MS1[gamma_1[i]], MS2[gamma_2[i]])
        if j != 0: return j
    return 0
