"""
The C3 algorithm

The C3 algorithm is used as method resolution order for new style
classes in Python. The implementation here is used to order the list
of super categories of a category.

AUTHOR:

- Simon King (2011-11): initial version.
"""
#*****************************************************************************
#  Copyright (C) 2011    Simon King <simon.king@uni-jena.de>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#                  http://www.gnu.org/licenses/
#******************************************************************************

include "../ext/python.pxi"

cpdef list C3_algorithm(object start, str bases, str attribute, bint proper):
    """
    An implementation of the C3 algorithm.

    C3 is the algorithm used by Python to construct the method
    resolution order for new style classes involving multiple
    inheritance.

    This implementation is used to order the list of super categories
    of a category; see
    :meth:`~sage.categories.category.Category.all_super_categories`.
    The purpose is to ensure that list of super categories matches
    with the method resolution order of the parent or element classes
    of a category.

    INPUT:

    - ``start`` -- an object; the returned list is built upon data
      provided by certain attributes of ``start``.
    - ``bases`` -- a string; the name of an attribute of ``start``
      providing a list of objects.
    - ``attribute`` -- a string; the name of an attribute of the
      objects provided in ``getattr(start,bases)``. That attribute is
      supposed to provide a list.

    ASSUMPTIONS:

    Our implementation of the algorithm only works on lists of
    objects that compare equal if and only if they are identical.

    OUTPUT:

    A list, the result of the C3 algorithm applied to the list
    ``[getattr(X,attribute) for X in getattr(start,bases)]``.

    EXAMPLES:

    We start with some categories having an inconsistent inheritance
    order::

        sage: class X(Category):
        ...    def super_categories(self):
        ...        return [Objects()]
        sage: class Y(Category):
        ...    def super_categories(self):
        ...        return [Objects()]
        sage: class A(Category):
        ...    def super_categories(self):
        ...        return [X(), Y()]
        sage: class B(Category):
        ...    def super_categories(self):
        ...        return [Y(), X()]
        sage: class Foo(Category):
        ...    def super_categories(self):
        ...       return [A(), B()]
        sage: F = Foo()

    Python is not able to create a consistent method resolution order
    for the parent class::

        sage: F.parent_class
        Traceback (most recent call last):
        ...
        TypeError: Cannot create a consistent method resolution
        order (MRO) for bases ....parent_class, ....parent_class

    Since the C3 algorithm is used for determining the list of
    all super categories (by trac ticket #11943), a similar error
    arises here::

        sage: F.all_super_categories()
        Traceback (most recent call last):
        ...
        ValueError: Can not merge the items Category of x, Category of y.

    Next, we demonstrate how our implementation of the C3 algorithm
    is used to compute the list of all super categories::

        sage: C = Category.join([HopfAlgebrasWithBasis(QQ), FiniteEnumeratedSets()])
        sage: from sage.misc.c3 import C3_algorithm
        sage: C3_algorithm(C,'_super_categories','_all_super_categories',True) == C._all_super_categories_proper
        True
        sage: C3_algorithm(C,'_super_categories','_all_super_categories',False) == C._all_super_categories
        True

    By trac ticket #11943, the following consistency tests are part
    of the test suites of categories (except for hom categories)::

        sage: C.parent_class.mro() == [x.parent_class for x in C.all_super_categories()]+[object]
        True
        sage: C.element_class.mro() == [x.element_class for x in C.all_super_categories()]+[object]
        True

    """
    # The lists in the arguments are reverted,
    # so that we can do pop() in lieue of pop(0).
    # In addition, containedness in the tail is tested using lists.
    cdef list out
    if proper:
        out = []
    else:
        out = [start]
    cdef list args = getattr(start,bases)
    if not args:
        return out
    cdef list curr_tail, tmp_tail
    cdef set curr_set, tmp_set
    cdef object curr_obj
    cdef bint next_item_found
    cdef list heads = []
    cdef list tails = []
    cdef list tailsets = []
    for curr_obj in args:
        curr_tail = getattr(curr_obj, attribute)
        heads.append(curr_tail[0])
        tmp_tail = PyList_GetSlice(curr_tail,1,PyList_GET_SIZE(curr_tail))
        PyList_Reverse(tmp_tail)
        tails.append(tmp_tail)
        tailsets.append(set(tmp_tail))
    cdef int i,j, lenargs
    lenargs = len(heads)
    for i from 0<=i<lenargs:
        curr_tail = <list>PyList_GET_ITEM(tails,i)
        if curr_tail is None:
            continue
        curr_set = <set>PyList_GET_ITEM(tailsets,i)
        O = <object>PyList_GET_ITEM(heads,i)
        next_item_found=True
        for j from 0<=j<i:
            tmp_tail = <list>PyList_GET_ITEM(tails,j)
            if tmp_tail is None:
                continue
            tmp_set = <set>PyList_GET_ITEM(tailsets,j)
            X = <object>PyList_GET_ITEM(heads,j)
            if X is O:
                try:
                    X = tmp_tail.pop()
                    heads[j] = X
                    tmp_set.remove(X)
                except IndexError:
                    tails[j] = None
            elif O in tmp_set:
                next_item_found=False
                break
        if next_item_found:
            for j from i<j<lenargs:
                tmp_tail = <list>PyList_GET_ITEM(tails,j)
                if tmp_tail is None:
                    continue
                tmp_set = <set>PyList_GET_ITEM(tailsets,j)
                X = <object>PyList_GET_ITEM(heads,j)
                if X is O:
                    try:
                        X = tmp_tail.pop()
                        heads[j] = X
                        tmp_set.remove(X)
                    except IndexError:
                        tails[j] = None
                elif O in tmp_set:
                    next_item_found=False
                    break
        if next_item_found:
            out.append(O)
            try:
                O = curr_tail.pop()
                heads[i] = O
                curr_set.remove(O)
            except IndexError:
                tails[i] = None

            i = -1
    # Either we need to raise an error, or the list is done.
    if tails.count(None)<lenargs:
        raise ValueError, "Can not merge the items %s."%', '.join([repr(heads[i]) for i,t in enumerate(tails) if t is not None])
    return out
