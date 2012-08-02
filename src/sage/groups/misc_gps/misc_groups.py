r"""
Miscellaneous Groups

This is a collection of groups that may not fit into some of the other infinite families described elsewhere.
"""

# Implementation note:
#     When adding a group here make an entry
#     in the  misc_groups_catalog.py module
#     or a more closely-related catalog module

def QuaternionMatrixGroupGF3():
    r"""
    The quaternion group as a set of `2\times 2` matrices over `GF(3)`.

    OUTPUT:

    A matrix group consisting of `2\times 2` matrices with
    elements from the finite field of order 3.  The group is
    the quaternion group, the nonabelian group of order 8 that
    is not isomorphic to the group of symmetries of a square
    (the dihedral group `D_4`).

    .. note::
        This group is most easily available via ``groups.matrix.QuaternionGF3()``.

    EXAMPLES:

    The generators are the matrix representations of the
    elements commonly called `I` and `J`, while `K`
    is the product of `I` and `J`. ::

        sage: from sage.groups.misc_gps.misc_groups import QuaternionMatrixGroupGF3
        sage: Q = QuaternionMatrixGroupGF3()
        sage: Q.order()
        8
        sage: aye = Q.gens()[0]; aye
        [1 1]
        [1 2]
        sage: jay = Q.gens()[1]; jay
        [2 1]
        [1 1]
        sage: kay = aye*jay; kay
        [0 2]
        [1 0]

    TESTS::

        sage: groups.matrix.QuaternionGF3()
        Matrix group over Finite Field of size 3 with 2 generators:
            [[[1, 1], [1, 2]], [[2, 1], [1, 1]]]

        sage: Q = QuaternionMatrixGroupGF3()
        sage: QP = Q.as_permutation_group()
        sage: QP.is_isomorphic(QuaternionGroup())
        True
        sage: H = DihedralGroup(4)
        sage: H.order()
        8
        sage: QP.is_abelian(), H.is_abelian()
        (False, False)
        sage: QP.is_isomorphic(H)
        False
    """
    from sage.rings.finite_rings.constructor import FiniteField
    from sage.matrix.matrix_space import MatrixSpace
    from sage.groups.matrix_gps.matrix_group import MatrixGroup
    MS = MatrixSpace(FiniteField(3), 2)
    aye = MS([1,1,1,2])
    jay = MS([2,1,1,1])
    return MatrixGroup([aye, jay])

