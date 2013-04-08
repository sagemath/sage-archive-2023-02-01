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

    TESTS:

    Regression test for bug #1 of :trac:`13501`::

        sage: class C(object): pass
        sage: class F(object): pass
        sage: class G(object): pass
        sage: class B(C,F):    pass
        sage: class D(F,G):    pass
        sage: class E(F):      pass
        sage: class A(B,D,E):  pass
        sage: [cls.__name__ for cls in A.mro()]
        ['A', 'B', 'C', 'D', 'E', 'F', 'G', 'object']

        sage: class Cs(Category):
        ...     def super_categories(self): return []
        sage: class Fs(Category):
        ...       def super_categories(self): return []
        sage: class Gs(Category):
        ...       def super_categories(self): return []
        sage: class Bs(Category):
        ...       def super_categories(self): return [Cs(), Fs()]
        sage: class Ds(Category):
        ...       def super_categories(self): return [Fs(), Gs()]
        sage: class Es(Category):
        ...       def super_categories(self): return [Fs()]
        sage: class As(Category):
        ...       def super_categories(self): return [Bs(), Ds(), Es()]
        sage: As().all_super_categories()
        [Category of as, Category of bs, Category of cs, Category of ds, Category of es, Category of fs, Category of gs]
        sage: TestSuite(As()).run(skip=["_test_pickling"])

    Regression test for bug #2 of :trac:`13501`. The following should
    fail since ``A`` asks for ``B`` to come before ``C``, where as
    ``B`` is a super class of ``C``::

        sage: class B(object): pass
        sage: class C(B): pass
        sage: class A(B, C): pass
        Traceback (most recent call last):
        ...
        TypeError: Error when calling the metaclass bases
            Cannot create a consistent method resolution
        order (MRO) for bases ...

        sage: class Bs(Category):
        ...       def super_categories(self): return []
        sage: class Cs(Category):
        ...       def super_categories(self): return [Bs()]
        sage: class As(Category):
        ...       def super_categories(self): return [Bs(), Cs()]
        ...       class subcategory_class(object): # Quick hack to skip the failure when computing the mro for subcategory_class
        ...            pass
        sage: As()
        Category of as
        sage: As().all_super_categories()
        Traceback (most recent call last):
        ...
        ValueError: Can not merge the items Category of bs, Category of cs, Category of bs.
    """
    cdef list out
    if proper:
        out = []
    else:
        out = [start]
    cdef list args = getattr(start,bases)
    if not args:
        return out
    # Data structure / invariants:
    # We will be working with the MRO's of the super objects
    # together with the list of bases of ``self``.
    # Each list is split between its head (in ``heads``) and tail (in
    # ``tails'') . Each tail is stored reversed, so that we can use a
    # cheap pop() in lieue of pop(0). A duplicate of the tail is
    # stored as a set in ``tailsets`` for cheap membership testing.
    # Since we actually want comparison by identity, not equality,
    # what we store is the set of memory locations of objects
    cdef object O, X
    cdef list tails = [getattr(obj, attribute) for obj in args]
    tails.append(args)
    tails              = [list(reversed(tail))                   for tail in tails]
    cdef list heads    = [tail.pop()                             for tail in tails]
    cdef list tailsets = [set([<size_t><void *>O for O in tail]) for tail in tails]

    cdef int i, j, nbheads
    nbheads = len(heads)
    cdef bint next_item_found
    cdef list tail_list

    while nbheads:
        for i from 0 <= i < nbheads:
            O = heads[i]
            # Does O appear in none of the tails?  ``all(O not in tail for tail in tailsets)``
            next_item_found = True
            for j from 0 <= j < nbheads:
                if j == i:
                    continue
                if <size_t><void *>O in <set>tailsets[j]:
                    next_item_found = False
                    break
            if next_item_found:
                out.append(O)
                # Clear O from other heads, removing the line altogether
                # if the tail is already empty.
                # j goes down so that ``del heads[j]`` does not screw up the numbering
                for j from nbheads > j >= 0:
                    if heads[j] is O:
                        tail_list = tails[j]
                        if tail_list:
                            X = tail_list.pop()
                            heads[j] = X
                            <set>tailsets[j].remove(<size_t><void *>X)
                        else:
                            del heads[j]
                            del tails[j]
                            del tailsets[j]
                            nbheads -= 1
                break
        if not next_item_found:
            # No head is available
            raise ValueError, "Can not merge the items %s."%', '.join([repr(head) for head in heads])
    return out
