"""
Fast functions for the category framework


AUTHOR:

- Simon King (initial version)

"""

#*****************************************************************************
#  Copyright (C) 2014      Simon King <simon.king@uni-jena.de>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#                  http://www.gnu.org/licenses/
#*****************************************************************************
include 'sage/ext/python.pxi'

cpdef inline tuple category_sort_key(object category):
    """
    Return ``category._cmp_key``.

    This helper function is used for sorting lists of categories.

    It is semantically equivalent to
    :func:`operator.attrgetter` ``("_cmp_key")``, but currently faster.

    EXAMPLES::

        sage: from sage.categories.category import category_sort_key
        sage: category_sort_key(Rings()) is Rings()._cmp_key
        True
    """
    return category._cmp_key

cpdef tuple _sort_uniq(categories):
    """
    Return the categories after sorting them and removing redundant categories.

    Redundant categories include duplicates and categories which
    are super categories of other categories in the input.

    INPUT:

    - ``categories`` -- a list (or iterable) of categories

    OUTPUT: a sorted tuple of mutually incomparable categories

    EXAMPLES::

        sage: Category._sort_uniq([Rings(), Monoids(), Coalgebras(QQ)])
        (Category of rings, Category of coalgebras over Rational Field)

    Note that, in the above example, ``Monoids()`` does not appear
    in the result because it is a super category of ``Rings()``.
    """
    cdef tuple cats = tuple(sorted(categories, key=category_sort_key, reverse=True))
    cdef list result = []
    cdef bint append
    for category in cats:
        append = True
        for cat in result:
            if cat.is_subcategory(category):
                append = False
                break
        if append:
            result.append(category)
    return tuple(result)


#############################################
## Axiom related functions

cdef class AxiomContainer(dict):
    def add(self, axiom):
        self[axiom] = len(self)
    def __iadd__(self, L):
        for axiom in L:
            self.add(axiom)
        return self

cpdef get_axiom_index(AxiomContainer all_axioms, str axiom):
    return <object>PyDict_GetItemString(all_axioms, PyString_AsString(axiom))

cpdef tuple canonicalize_axioms(AxiomContainer all_axioms, axioms):
    r"""
    Canonicalize a set of axioms.

    INPUT:

     - ``axioms`` -- a set (or iterable) of axioms

    OUTPUT:

    A set of axioms as a tuple sorted according to the order of the
    tuple ``all_axioms`` in :mod:`sage.categories.category_with_axiom`.

    EXAMPLES::

        sage: from sage.categories.category_with_axiom import canonicalize_axioms
        sage: canonicalize_axioms(["Commutative", "Connected", "WithBasis", "Finite"])
        ('Finite', 'Connected', 'WithBasis', 'Commutative')
        sage: canonicalize_axioms(["Commutative", "Connected", "Commutative", "WithBasis", "Finite"])
        ('Finite', 'Connected', 'WithBasis', 'Commutative')
    """
    cdef list L = list(set(axioms))
    L.sort(key = (all_axioms).__getitem__)
    return tuple(L)
