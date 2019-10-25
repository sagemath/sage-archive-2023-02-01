from __future__ import absolute_import

from . import permgroup_element

def PermutationGroupElement(g, parent=None, check=True):
    r"""
    Build a permutation on {1, 2, ..., n}.

    EXAMPLES:

    Initialization as a list of images::

        sage: p = PermutationGroupElement([1,4,2,3])
        sage: p
        (2,4,3)
        sage: p.parent()
        Symmetric group of order 4! as a permutation group

    Initialization as a list of cycles::

        sage: p = PermutationGroupElement([(3,5),(4,6,9)])
        sage: p
        (3,5)(4,6,9)
        sage: p.parent()
        Symmetric group of order 9! as a permutation group

    Initialization as a string representing a cycle decomposition::

        sage: p = PermutationGroupElement('(2,4)(3,5)')
        sage: p
        (2,4)(3,5)
        sage: p.parent()
        Symmetric group of order 5! as a permutation group
    """
    if isinstance(g, permgroup_element.PermutationGroupElement):
        if parent is None or g.parent() is parent:
            return g

    if parent is None:
        from sage.groups.perm_gps.permgroup_named import SymmetricGroup
        from sage.combinat.permutation import from_cycles

        convert_dict = parent._domain_to_gap if parent is not None else None
        try:
            v = permgroup_element.standardize_generator(g, convert_dict)
        except KeyError:
            raise ValueError("Invalid permutation vector: %s" % g)
        degree = max([1] + [max(cycle+(1,)) for cycle in v])
        v = from_cycles(degree, v)

        if parent is None:
            parent = SymmetricGroup(len(v))

    return parent.element_class(g, parent, check)


