r"""
Crystals of letters
"""

#*****************************************************************************
#       Copyright (C) 2007 Anne Schilling <anne at math.ucdavis.edu>
#                          Nicolas M. Thiery <nthiery at users.sf.net>
#                          Daniel Bump    <bump at match.stanford.edu>
#                          Brant Jones    <brant at math.ucdavis.edu>
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 2 of the License, or
# (at your option) any later version.
#                  http://www.gnu.org/licenses/
#*****************************************************************************

from cpython.object cimport Py_EQ, Py_NE, Py_LE, Py_GE, Py_LT, Py_GT
from sage.misc.cachefunc import cached_method
from sage.misc.lazy_attribute import lazy_attribute
from sage.structure.unique_representation import UniqueRepresentation
from sage.structure.parent import Parent
from sage.structure.element cimport Element
from sage.categories.enumerated_sets import EnumeratedSets
from sage.categories.classical_crystals import ClassicalCrystals
from sage.combinat.root_system.cartan_type import CartanType

def CrystalOfLetters(cartan_type, element_print_style=None, dual=None):
    r"""
    Return the crystal of letters of the given type.

    For classical types, this is a combinatorial model for the crystal
    with highest weight `\Lambda_1` (the first fundamental weight).

    Any irreducible classical crystal appears as the irreducible
    component of the tensor product of several copies of this crystal
    (plus possibly one copy of the spin crystal, see
    :class:`~sage.combinat.crystals.spins.CrystalOfSpins`).
    See [KN94]_. Elements of this irreducible component have a fixed shape,
    and can be fit inside a tableau shape. Otherwise said, any irreducible
    classical crystal is isomorphic to a crystal of tableaux with cells
    filled by elements of the crystal of letters (possibly tensored with
    the crystal of spins).

    INPUT:

    - ``T`` -- A Cartan type

    REFERENCES:

    .. [KN94] \M. Kashiwara and T. Nakashima.
       Crystal graphs for representations of the `q`-analogue of classical Lie
       algebras.
       J. Algebra **165**, no. 2, pp. 295--345, 1994.

    EXAMPLES::

        sage: C = crystals.Letters(['A',5])
        sage: C.list()
        [1, 2, 3, 4, 5, 6]
        sage: C.cartan_type()
        ['A', 5]

    For type `E_6`, one can also specify how elements are printed.
    This option is usually set to None and the default representation is used.
    If one chooses the option 'compact', the elements are printed in the more
    compact convention with 27 letters ``+abcdefghijklmnopqrstuvwxyz`` and
    the 27 letters ``-ABCDEFGHIJKLMNOPQRSTUVWXYZ`` for the dual crystal.

    EXAMPLES::

        sage: C = crystals.Letters(['E',6], element_print_style = 'compact')
        sage: C
        The crystal of letters for type ['E', 6]
        sage: C.list()
        [+, a, b, c, d, e, f, g, h, i, j, k, l, m, n, o, p, q, r, s, t, u, v, w, x, y, z]
        sage: C = crystals.Letters(['E',6], element_print_style = 'compact', dual = True)
        sage: C
        The crystal of letters for type ['E', 6] (dual)
        sage: C.list()
        [-, A, B, C, D, E, F, G, H, I, J, K, L, M, N, O, P, Q, R, S, T, U, V, W, X, Y, Z]
    """
    ct = CartanType(cartan_type)
    if ct.letter == 'A':
        return ClassicalCrystalOfLetters(ct, Crystal_of_letters_type_A_element)
    elif ct.letter == 'B':
        return ClassicalCrystalOfLetters(ct, Crystal_of_letters_type_B_element)
    elif ct.letter == 'C':
        return ClassicalCrystalOfLetters(ct, Crystal_of_letters_type_C_element)
    elif ct.letter == 'D':
        return ClassicalCrystalOfLetters(ct, Crystal_of_letters_type_D_element)
    elif ct.letter == 'E' and ct.rank() == 6:
        if dual is None:
            return ClassicalCrystalOfLetters(ct,
                                             Crystal_of_letters_type_E6_element,
                                             element_print_style)
        else:
            return ClassicalCrystalOfLetters(ct,
                                             Crystal_of_letters_type_E6_element_dual,
                                             element_print_style, dual = True)
    elif ct.letter == 'E' and ct.rank() == 7:
        return ClassicalCrystalOfLetters(ct, Crystal_of_letters_type_E7_element)
    elif ct.letter == 'G':
        return ClassicalCrystalOfLetters(ct, Crystal_of_letters_type_G_element)
    else:
        raise NotImplementedError

class ClassicalCrystalOfLetters(UniqueRepresentation, Parent):
    r"""
    A generic class for classical crystals of letters.

    All classical crystals of letters should be instances of this class
    or of subclasses. To define a new crystal of letters, one only
    needs to implement a class for the elements (which subclasses
    :class:`~sage.combinat.crystals.Letter`), with appropriate
    `e_i` and `f_i` operations. If the module generator is not `1`, one also
    needs to define the subclass
    :class:`~sage.combinat.crystals.letters.ClassicalCrystalOfLetters` for the
    crystal itself.

    The basic assumption is that crystals of letters are small, but
    used intensively as building blocks. Therefore, we explicitly build
    in memory the list of all elements, the crystal graph and its
    transitive closure, so as to make the following operations constant
    time: ``list``, ``cmp``, (todo: ``phi``, ``epsilon``, ``e``, and
    ``f`` with caching)
    """
    def __init__(self, cartan_type, element_class, element_print_style = None, dual = None):
        """
        EXAMPLES::

            sage: C = crystals.Letters(['A',5])
            sage: C.category()
            Category of classical crystals
            sage: TestSuite(C).run()
        """
        self.Element = element_class
        Parent.__init__(self, category = ClassicalCrystals())
        self._cartan_type = CartanType(cartan_type)
        self.rename("The crystal of letters for type %s"%self._cartan_type)
        if cartan_type.type() == 'E':
            if cartan_type.rank() == 6:
                if dual:
                    self.module_generators = (self._element_constructor_((6,)),)
                    self._ambient = CrystalOfLetters(CartanType(['E',6]))
                    self.rename("%s (dual)"%self)
                else:
                    self.module_generators = (self._element_constructor_((1,)),)
            elif cartan_type.rank() == 7:
                self.module_generators = (self._element_constructor_((7,)),)
            self._list = [x for x in super(ClassicalCrystalOfLetters, self).__iter__()]
        else:
            self.module_generators = (self._element_constructor_(1),)
            if cartan_type.type() == 'G':
                self._list = [self._element_constructor_(1),
                              self._element_constructor_(2),
                              self._element_constructor_(3),
                              self._element_constructor_(0),
                              self._element_constructor_(-3),
                              self._element_constructor_(-2),
                              self._element_constructor_(-1)]
            else:
                self._list = [self._element_constructor_(i)
                                for i in xrange(1, cartan_type.rank()+1)]
                if cartan_type.type() == 'B':
                    self._list.append(self._element_constructor_(0))
                if cartan_type.type() != 'A':
                    self._list += [self._element_constructor_(-i)
                                     for i in xrange(cartan_type.rank(), 0, -1)]
                else:
                    self._list.append(self._element_constructor_(cartan_type.rank()+1))
        self._element_print_style = element_print_style

    def __call__(self, value):
        """
        Parse input to valid values to give to ``_element_constructor_()``.

        EXAMPLES::

            sage: C = crystals.Letters(['E',6])
            sage: c = C((1,))
            sage: C([1]) == c
            True
        """
        if value.__class__ == self.element_class and value.parent() is self:
            return value
        if isinstance(value, list):
            return self._element_constructor_(tuple(value))
        return self._element_constructor_(value)

    @cached_method
    def _element_constructor_(self, value):
        """
        Convert ``value`` into an element of ``self``.

        EXAMPLES::

            sage: C = crystals.Letters(['A',5])
            sage: c = C(1); c
            1
            sage: c.parent()
            The crystal of letters for type ['A', 5]
            sage: c is C(c)
            True
        """
        if value == 'E':
            return EmptyLetter(self)
        else: # Should do sanity checks!
            return self.element_class(self, value)

    def __iter__(self):
        """
        Iterate through ``self``.

        EXAMPLES::

            sage: C = crystals.Letters(['A',5])
            sage: [x for x in C]
            [1, 2, 3, 4, 5, 6]
        """
        return iter(self._list)

    def list(self):
        """
        Return a list of the elements of ``self``.

        EXAMPLES::

            sage: C = crystals.Letters(['A',5])
            sage: C.list()
            [1, 2, 3, 4, 5, 6]
        """
        return self._list

    @lazy_attribute
    def _digraph_closure(self):
        """
        The transitive closure of the directed graph associated to ``self``.

        EXAMPLES::

            sage: crystals.Letters(['A',5])._digraph_closure
            Transitive closure of : Digraph on 6 vertices
        """
        return self.digraph().transitive_closure()

    def __contains__(self, x):
        """
        EXAMPLES::

            sage: C = crystals.Letters(['A',5])
            sage: 1 in C
            False
            sage: C(1) in C
            True
        """
        return x in self._list

    def lt_elements(self, x, y):
        r"""
        Return ``True`` if and only if there is a path from ``x`` to ``y`` in
        the crystal graph, when ``x`` is not equal to ``y``.

        Because the crystal graph is classical, it is a directed acyclic
        graph which can be interpreted as a poset. This function implements
        the comparison function of this poset.

        EXAMPLES::

            sage: C = crystals.Letters(['A', 5])
            sage: x = C(1)
            sage: y = C(2)
            sage: C.lt_elements(x,y)
            True
            sage: C.lt_elements(y,x)
            False
            sage: C.lt_elements(x,x)
            False
            sage: C = crystals.Letters(['D', 4])
            sage: C.lt_elements(C(4),C(-4))
            False
            sage: C.lt_elements(C(-4),C(4))
            False
        """
        if x.parent() is not self or y.parent() is not self:
            raise ValueError("Cannot compare elements of different parents")
        if self._digraph_closure.has_edge(x,y):
            return True
        return False

    # temporary woraround while an_element is overriden by Parent
    _an_element_ = EnumeratedSets.ParentMethods._an_element_

# Utility. Note: much of this class should be factored out at some point!
cdef class Letter(Element):
    r"""
    A class for letters.

    Like :class:`ElementWrapper`, plus delegates ``__lt__`` (comparison)
    to the parent.

    EXAMPLES::

        sage: from sage.combinat.crystals.letters import Letter
        sage: a = Letter(ZZ, 1)
        sage: Letter(ZZ, 1).parent()
        Integer Ring

        sage: Letter(ZZ, 1)._repr_()
        '1'

        sage: parent1 = ZZ  # Any fake value ...
        sage: parent2 = QQ  # Any fake value ...
        sage: l11 = Letter(parent1, 1)
        sage: l12 = Letter(parent1, 2)
        sage: l21 = Letter(parent2, 1)
        sage: l22 = Letter(parent2, 2)
        sage: l11 == l11
        True
        sage: l11 == l12
        False
        sage: l11 == l21 # not tested
        False

        sage: C = crystals.Letters(['B', 3])
        sage: C(0) != C(0)
        False
        sage: C(1) != C(-1)
        True
    """
    cdef readonly int value

    def __init__(self, parent, int value):
        """
        EXAMPLES::

            sage: C = crystals.Letters(['B',4])
            sage: a = C(3)
            sage: TestSuite(a).run()
        """
        self.value = value
        Element.__init__(self, parent)

    def __setstate__(self, state):
        r"""
        Used in unpickling old pickles.

        EXAMPLES::

            sage: C = crystals.Letters(['B',4])
            sage: a = C(3)
            sage: loads(dumps(a)) == a
            True
        """
        self._set_parent(state[0])
        self.value = state[1]['value']

    def __reduce__(self):
        r"""
        Used in pickling crystal of letters elements.

        EXAMPLES::

            sage: C = crystals.Letters(['A',3])
            sage: a = C(1)
            sage: a.__reduce__()
            (The crystal of letters for type ['A', 3], (1,))
        """
        return (self._parent, (self.value,))

    def _repr_(self):
        """
        Return a string representation of ``self``.

        EXAMPLES::

            sage: C = crystals.Letters(['B', 3])
            sage: C(0)
            0
            sage: C(1)
            1
            sage: C(-1)
            -1
        """
        return repr(self.value)

    def _latex_(self):
        r"""
        A latex representation of ``self``.

        EXAMPLES::

            sage: C = crystals.Letters(['D', 4])
            sage: latex(C(2))
            2
            sage: latex(C(-3))
            \overline{3}
        """
        if self.value < 0:
            return "\\overline{" + repr(-self.value) + "}"
        return repr(self.value)

    def __hash__(self):
        """
        Return the hash value of ``self``.

        EXAMPLES::

            sage: C = crystals.Letters(['D', 4])
            sage: hash(C(4)) == hash(4)
            True
        """
        return self.value

    cpdef _richcmp_(left, Element right, int op):
        """
        Return ``True`` if ``left`` compares with ``right`` based on ``op``.

        EXAMPLES::

            sage: C = crystals.Letters(['D', 4])
            sage: C(4) > C(-4) # indirect doctest
            False
            sage: C(4) < C(-3)
            True
            sage: C(4) == C(4)
            True

        TESTS::

            sage: C = crystals.Letters(['C', 3])
            sage: C('E') == C(2)
            False
            sage: C(2) == C('E')
            False
            sage: C('E') == C('E')
            True
        """
        # Special case for the empty letter
        if isinstance(left, EmptyLetter):
            return isinstance(right, EmptyLetter) \
                   and (op == Py_EQ or op == Py_LE or op == Py_GE)
        if isinstance(right, EmptyLetter):
            return op == Py_NE

        cdef Letter self, x
        self = left
        x = right
        if op == Py_EQ:
            return self.value == x.value
        if op == Py_NE:
            return self.value != x.value
        if op == Py_LT:
            return self._parent.lt_elements(self, x)
        if op == Py_GT:
            return x.parent().lt_elements(x, self)
        if op == Py_LE:
            return self.value == x.value or self._parent.lt_elements(self, x)
        if op == Py_GE:
            return self.value == x.value or x.parent().lt_elements(x, self)
        return False

cdef class EmptyLetter(Element):
    """
    The affine letter `\emptyset` thought of as a classical crystal letter
    in classical type `B_n` and `C_n`.

    .. WARNING::

        This is not a classical letter.

    Used in the rigged configuration bijections.
    """
    cdef readonly str value

    def __init__(self, parent):
        """
        Initialize ``self``.

        EXAMPLES::

            sage: C = crystals.Letters(['C', 3])
            sage: TestSuite(C('E')).run()
        """
        self.value = 'E'
        Element.__init__(self, parent)

    def __reduce__(self):
        r"""
        Used in pickling crystal of letters elements.

        EXAMPLES::

            sage: C = crystals.Letters(['C',3])
            sage: a = C('E')
            sage: a.__reduce__()
            (The crystal of letters for type ['C', 3], ('E',))
        """
        return (self._parent, ('E',))

    def _repr_(self):
        """
        Return a string representation of ``self``.

        EXAMPLES::

            sage: C = crystals.Letters(['C', 3])
            sage: C('E')
            E
        """
        return 'E'

    def _latex_(self):
        """
        Return a latex representation of ``self``.

        EXAMPLES::

            sage: C = crystals.Letters(['C', 3])
            sage: latex(C('E'))
            \emptyset
        """
        return "\\emptyset"

    def __hash__(self):
        """
        Return the hash value of ``self``.

        EXAMPLES::

            sage: C = crystals.Letters(['D', 4])
            sage: hash(C('E')) == hash('E')
            True
        """
        return hash(self.value)

    cpdef _richcmp_(left, Element right, int op):
        """
        Return ``True`` if ``left`` compares with ``right`` based on ``op``.

        EXAMPLES::

            sage: C = crystals.Letters(['C', 3])
            sage: C('E') == C(2)
            False
            sage: C('E') < C(2)
            False
            sage: C('E') <= C(2)
            False
            sage: C('E') != C(2)
            True
            sage: C('E') == C('E')
            True
            sage: C('E') != C('E')
            False
            sage: C('E') >= C('E')
            True
            sage: C('E') < C('E')
            False
        """
        if isinstance(left, EmptyLetter) and isinstance(right, EmptyLetter):
           return op == Py_EQ or op == Py_LE or op == Py_GE
        return op == Py_NE

    def weight(self):
        """
        Return the weight of ``self``.

        EXAMPLES::

            sage: C = crystals.Letters(['C', 3])
            sage: C('E').weight()
            (0, 0, 0)
        """
        return self.parent().weight_lattice_realization().zero()

    cpdef e(self, int i):
        """
        Return `e_i` of ``self`` which is ``None``.

        EXAMPLES::

            sage: C = crystals.Letters(['C', 3])
            sage: C('E').e(1)
        """
        return None

    cpdef f(self, int i):
        """
        Return `f_i` of ``self`` which is ``None``.

        EXAMPLES::

            sage: C = crystals.Letters(['C', 3])
            sage: C('E').f(1)
        """
        return None

    cpdef int epsilon(self, int i):
        r"""
        Return `\varepsilon_i` of ``self``.

        EXAMPLES::

            sage: C = crystals.Letters(['C', 3])
            sage: C('E').epsilon(1)
            0
        """
        return 0

    cpdef int phi(self, int i):
        r"""
        Return `\varphi_i` of ``self``.

        EXAMPLES::

            sage: C = crystals.Letters(['C', 3])
            sage: C('E').phi(1)
            0
        """
        return 0

#########################
# Type A
#########################

cdef class Crystal_of_letters_type_A_element(Letter):
    r"""
    Type `A` crystal of letters elements.

    TESTS::

        sage: C = crystals.Letters(['A',3])
        sage: C.list()
        [1, 2, 3, 4]
        sage: [ [x < y for y in C] for x in C ]
        [[False, True, True, True],
         [False, False, True, True],
         [False, False, False, True],
         [False, False, False, False]]

    ::

        sage: C = crystals.Letters(['A',5])
        sage: C(1) < C(1), C(1) < C(2), C(1) < C(3), C(2) < C(1)
        (False, True, True, False)

    ::

        sage: TestSuite(C).run()
    """
    def weight(self):
        """
        Return the weight of ``self``.

        EXAMPLES::

            sage: [v.weight() for v in crystals.Letters(['A',3])]
            [(1, 0, 0, 0), (0, 1, 0, 0), (0, 0, 1, 0), (0, 0, 0, 1)]
        """
        return self._parent.weight_lattice_realization().monomial(self.value-1)

    cpdef Letter e(self, int i):
        r"""
        Return the action of `e_i` on ``self``.

        EXAMPLES::

            sage: C = crystals.Letters(['A',4])
            sage: [(c,i,c.e(i)) for i in C.index_set() for c in C if c.e(i) is not None]
            [(2, 1, 1), (3, 2, 2), (4, 3, 3), (5, 4, 4)]
        """
        if self.value == i+1:
            return self._parent._element_constructor_(self.value-1)
        else:
            return None

    cpdef Letter f(self, int i):
        r"""
        Return the action of `f_i` on ``self``.

        EXAMPLES::

            sage: C = crystals.Letters(['A',4])
            sage: [(c,i,c.f(i)) for i in C.index_set() for c in C if c.f(i) is not None]
            [(1, 1, 2), (2, 2, 3), (3, 3, 4), (4, 4, 5)]
        """
        if self.value == i:
            return self._parent._element_constructor_(self.value+1)
        else:
            return None

    cpdef int epsilon(self, int i):
        r"""
        Return `\varepsilon_i` of ``self``.

        EXAMPLES::

            sage: C = crystals.Letters(['A',4])
            sage: [(c,i) for i in C.index_set() for c in C if c.epsilon(i) != 0]
            [(2, 1), (3, 2), (4, 3), (5, 4)]
        """
        if self.value == i+1:
            return 1
        return 0

    cpdef int phi(self, int i):
        r"""
        Return `\varphi_i` of ``self``.

        EXAMPLES::

            sage: C = crystals.Letters(['A',4])
            sage: [(c,i) for i in C.index_set() for c in C if c.phi(i) != 0]
            [(1, 1), (2, 2), (3, 3), (4, 4)]
        """
        if self.value == i:
            return 1
        return 0

#########################
# Type B
#########################

cdef class Crystal_of_letters_type_B_element(Letter):
    r"""
    Type `B` crystal of letters elements.

    TESTS::

        sage: C = crystals.Letters(['B',3])
        sage: TestSuite(C).run()
    """
    def weight(self):
        """
        Return the weight of ``self``.

        EXAMPLES::

            sage: [v.weight() for v in crystals.Letters(['B',3])]
            [(1, 0, 0),
             (0, 1, 0),
             (0, 0, 1),
             (0, 0, 0),
             (0, 0, -1),
             (0, -1, 0),
             (-1, 0, 0)]
        """
        if self.value > 0:
            return self._parent.weight_lattice_realization().monomial(self.value-1)
        elif self.value < 0:
            return -self._parent.weight_lattice_realization().monomial(-self.value-1)
        else:
            return self._parent.weight_lattice_realization()(0)

    cpdef Letter e(self, int i):
        r"""
        Return the action of `e_i` on ``self``.

        EXAMPLES::

            sage: C = crystals.Letters(['B',4])
            sage: [(c,i,c.e(i)) for i in C.index_set() for c in C if c.e(i) is not None]
            [(2, 1, 1),
             (-1, 1, -2),
             (3, 2, 2),
             (-2, 2, -3),
             (4, 3, 3),
             (-3, 3, -4),
             (0, 4, 4),
             (-4, 4, 0)]
        """
        if self.value == i+1:
            return self._parent._element_constructor_(i)
        elif self.value == 0 and i == self._parent._cartan_type.n:
            return self._parent._element_constructor_(self._parent._cartan_type.n)
        elif self.value == -i:
            if i == self._parent._cartan_type.n:
                return self._parent._element_constructor_(0)
            else:
                return self._parent._element_constructor_(-i-1)
        else:
            return None

    cpdef Letter f(self, int i):
        r"""
        Return the actions of `f_i` on ``self``.

        EXAMPLES::

            sage: C = crystals.Letters(['B',4])
            sage: [(c,i,c.f(i)) for i in C.index_set() for c in C if c.f(i) is not None]
            [(1, 1, 2),
             (-2, 1, -1),
             (2, 2, 3),
             (-3, 2, -2),
             (3, 3, 4),
             (-4, 3, -3),
             (4, 4, 0),
             (0, 4, -4)]
        """
        if self.value == i:
            if i < self._parent._cartan_type.n:
                return self._parent._element_constructor_(i+1)
            else:
                return self._parent._element_constructor_(0)
        elif self.value == 0 and i == self._parent._cartan_type.n:
            return self._parent._element_constructor_(-self._parent._cartan_type.n)
        elif self.value == -i-1:
            return self._parent._element_constructor_(-i)
        else:
            return None

    cpdef int epsilon(self, int i):
        r"""
        Return `\varepsilon_i` of ``self``.

        EXAMPLES::

            sage: C = crystals.Letters(['B',3])
            sage: [(c,i) for i in C.index_set() for c in C if c.epsilon(i) != 0]
            [(2, 1), (-1, 1), (3, 2), (-2, 2), (0, 3), (-3, 3)]
        """
        cdef int n = self._parent._cartan_type.n
        if self.value == 0:
            if i == n:
                return 1
            return 0
        if i == n and self.value == -n:
            return 2
        if self.value == i+1 or self.value == -i:
            return 1
        return 0

    cpdef int phi(self, int i):
        r"""
        Return `\varphi_i` of ``self``.

        EXAMPLES::

            sage: C = crystals.Letters(['B',3])
            sage: [(c,i) for i in C.index_set() for c in C if c.phi(i) != 0]
            [(1, 1), (-2, 1), (2, 2), (-3, 2), (3, 3), (0, 3)]
        """
        cdef int n = self._parent._cartan_type.n
        if self.value == 0:
            if i == n:
                return 1
            return 0
        if i == n and self.value == n:
            return 2
        if self.value == i or self.value == -i-1:
            return 1
        return 0

#########################
# Type C
#########################

cdef class Crystal_of_letters_type_C_element(Letter):
    r"""
    Type `C` crystal of letters elements.

    TESTS::

        sage: C = crystals.Letters (['C',3])
        sage: C.list()
        [1, 2, 3, -3, -2, -1]
        sage: [ [x < y for y in C] for x in C ]
        [[False, True, True, True, True, True],
         [False, False, True, True, True, True],
         [False, False, False, True, True, True],
         [False, False, False, False, True, True],
         [False, False, False, False, False, True],
         [False, False, False, False, False, False]]
        sage: TestSuite(C).run()
    """
    def weight(self):
        """
        Return the weight of ``self``.

        EXAMPLES::

            sage: [v.weight() for v in crystals.Letters(['C',3])]
            [(1, 0, 0), (0, 1, 0), (0, 0, 1), (0, 0, -1), (0, -1, 0), (-1, 0, 0)]
        """
        if self.value > 0:
            return self._parent.weight_lattice_realization().monomial(self.value-1)
        elif self.value < 0:
            return -self._parent.weight_lattice_realization().monomial(-self.value-1)
        else:
            return self._parent.weight_lattice_realization()(0)

    cpdef Letter e(self, int i):
        r"""
        Return the action of `e_i` on ``self``.

        EXAMPLES::

            sage: C = crystals.Letters(['C',4])
            sage: [(c,i,c.e(i)) for i in C.index_set() for c in C if c.e(i) is not None]
            [(2, 1, 1),
             (-1, 1, -2),
             (3, 2, 2),
             (-2, 2, -3),
             (4, 3, 3),
             (-3, 3, -4),
             (-4, 4, 4)]
        """
        if self.value == -self._parent._cartan_type.n and self.value == -i:
            return self._parent._element_constructor_(-self.value)
        elif self.value == i+1 or self.value == -i:
            return self._parent._element_constructor_(self.value-1)
        else:
            return None

    cpdef Letter f(self, int i):
        r"""
        Return the action of `f_i` on ``self``.

        EXAMPLES::

            sage: C = crystals.Letters(['C',4])
            sage: [(c,i,c.f(i)) for i in C.index_set() for c in C if c.f(i) is not None]
            [(1, 1, 2), (-2, 1, -1), (2, 2, 3),
             (-3, 2, -2), (3, 3, 4), (-4, 3, -3), (4, 4, -4)]
        """
        if self.value == self._parent._cartan_type.n and self.value == i:
            return self._parent._element_constructor_(-self.value)
        elif self.value == i or self.value == -i-1:
            return self._parent._element_constructor_(self.value+1)
        else:
            return None

    cpdef int epsilon(self, int i):
        r"""
        Return `\varepsilon_i` of ``self``.

        EXAMPLES::

            sage: C = crystals.Letters(['C',3])
            sage: [(c,i) for i in C.index_set() for c in C if c.epsilon(i) != 0]
            [(2, 1), (-1, 1), (3, 2), (-2, 2), (-3, 3)]
        """
        if self.value == i+1 or self.value == -i:
            return 1
        return 0

    cpdef int phi(self, int i):
        r"""
        Return `\varphi_i` of ``self``.

        EXAMPLES::

            sage: C = crystals.Letters(['C',3])
            sage: [(c,i) for i in C.index_set() for c in C if c.phi(i) != 0]
            [(1, 1), (-2, 1), (2, 2), (-3, 2), (3, 3)]
        """
        if self.value == i or self.value == -i-1:
            return 1
        return 0

#########################
# Type D
#########################

cdef class Crystal_of_letters_type_D_element(Letter):
    r"""
    Type `D` crystal of letters elements.

    TESTS::

        sage: C = crystals.Letters(['D',4])
        sage: C.list()
        [1, 2, 3, 4, -4, -3, -2, -1]
        sage: TestSuite(C).run()
    """
    def weight(self):
        """
        Return the weight of ``self``.

        EXAMPLES::

            sage: [v.weight() for v in crystals.Letters(['D',4])]
            [(1, 0, 0, 0),
             (0, 1, 0, 0),
             (0, 0, 1, 0),
             (0, 0, 0, 1),
             (0, 0, 0, -1),
             (0, 0, -1, 0),
             (0, -1, 0, 0),
             (-1, 0, 0, 0)]
        """
        if self.value > 0:
            return self._parent.weight_lattice_realization().monomial(self.value-1)
        elif self.value < 0:
            return -self._parent.weight_lattice_realization().monomial(-self.value-1)
        else:
            return self._parent.weight_lattice_realization()(0)

    cpdef Letter e(self, int i):
        r"""
        Return the action of `e_i` on ``self``.

        EXAMPLES::

            sage: C = crystals.Letters(['D',5])
            sage: [(c,i,c.e(i)) for i in C.index_set() for c in C if c.e(i) is not None]
            [(2, 1, 1),
             (-1, 1, -2),
             (3, 2, 2),
             (-2, 2, -3),
             (4, 3, 3),
             (-3, 3, -4),
             (5, 4, 4),
             (-4, 4, -5),
             (-5, 5, 4),
             (-4, 5, 5)]
        """
        if i == self._parent._cartan_type.n:
            if self.value == -i:
                return self._parent._element_constructor_(i-1)
            elif self.value == -(i-1):
                return self._parent._element_constructor_(i)
            else:
                return None
        elif self.value == i+1:
            return self._parent._element_constructor_(i)
        elif self.value == -i:
            return self._parent._element_constructor_(-(i+1))
        else:
            return None

    cpdef Letter f(self, int i):
        r"""
        Return the action of `f_i` on ``self``.

        EXAMPLES::

            sage: C = crystals.Letters(['D',5])
            sage: [(c,i,c.f(i)) for i in C.index_set() for c in C if c.f(i) is not None]
            [(1, 1, 2),
             (-2, 1, -1),
             (2, 2, 3),
             (-3, 2, -2),
             (3, 3, 4),
             (-4, 3, -3),
             (4, 4, 5),
             (-5, 4, -4),
             (4, 5, -5),
             (5, 5, -4)]
        """
        if i == self.value:
            if i == self._parent._cartan_type.n:
                return self._parent._element_constructor_(-(i-1))
            else:
                return self._parent._element_constructor_(i+1)
        elif self.value == -(i+1):
            return self._parent._element_constructor_(-i)
        elif self.value == self._parent._cartan_type.n-1 and i == self.value+1:
            return self._parent._element_constructor_(-i)
        else:
            return None

    cpdef int epsilon(self, int i):
        r"""
        Return `\varepsilon_i` of ``self``.

        EXAMPLES::

            sage: C = crystals.Letters(['D',4])
            sage: [(c,i) for i in C.index_set() for c in C if c.epsilon(i) != 0]
            [(2, 1), (-1, 1), (3, 2), (-2, 2), (4, 3), (-3, 3), (-4, 4), (-3, 4)]
        """
        if self.value == i+1 or self.value == -i:
            return 1
        cdef int n = self._parent._cartan_type.n
        if i == n and self.value == -n+1:
            return 1
        return 0

    cpdef int phi(self, int i):
        r"""
        Return `\varphi_i` of ``self``.

        EXAMPLES::

            sage: C = crystals.Letters(['D',4])
            sage: [(c,i) for i in C.index_set() for c in C if c.phi(i) != 0]
            [(1, 1), (-2, 1), (2, 2), (-3, 2), (3, 3), (-4, 3), (3, 4), (4, 4)]
        """
        if self.value == i or self.value == -i-1:
            return 1
        cdef int n = self._parent._cartan_type.n
        if i == n and self.value == n-1:
            return 1
        return 0

#########################
# Type G2
#########################

cdef class Crystal_of_letters_type_G_element(Letter):
    r"""
    Type `G_2` crystal of letters elements.

    TESTS::

        sage: C = crystals.Letters(['G',2])
        sage: C.list()
        [1, 2, 3, 0, -3, -2, -1]
        sage: TestSuite(C).run()
    """
    def weight(self):
        """
        Return the weight of ``self``.

        EXAMPLES::

            sage: [v.weight() for v in crystals.Letters(['G',2])]
            [(1, 0, -1), (1, -1, 0), (0, 1, -1), (0, 0, 0), (0, -1, 1), (-1, 1, 0), (-1, 0, 1)]
        """
        if self.value == 1:
            return self._parent.weight_lattice_realization()((1, 0, -1))
        elif self.value == 2:
            return self._parent.weight_lattice_realization()((1, -1, 0))
        elif self.value == 3:
            return self._parent.weight_lattice_realization()((0, 1, -1))
        elif self.value == 0:
            return self._parent.weight_lattice_realization()((0, 0, 0))
        elif self.value == -3:
            return self._parent.weight_lattice_realization()((0, -1, 1))
        elif self.value == -2:
            return self._parent.weight_lattice_realization()((-1, 1, 0))
        elif self.value == -1:
            return self._parent.weight_lattice_realization()((-1, 0, 1))
        else:
            raise RuntimeError("G2 crystal of letters element %d not valid"%self.value)

    cpdef Letter e(self, int i):
        r"""
        Return the action of `e_i` on ``self``.

        EXAMPLES::

            sage: C = crystals.Letters(['G',2])
            sage: [(c,i,c.e(i)) for i in C.index_set() for c in C if c.e(i) is not None]
            [(2, 1, 1),
             (0, 1, 3),
             (-3, 1, 0),
             (-1, 1, -2),
             (3, 2, 2),
             (-2, 2, -3)]
        """
        if i == 1:
            if self.value == 2:
                return self._parent._element_constructor_(1)
            elif self.value == 0:
                return self._parent._element_constructor_(3)
            elif self.value == -3:
                return self._parent._element_constructor_(0)
            elif self.value == -1:
                return self._parent._element_constructor_(-2)
            else:
                return None
        else:
            if self.value == 3:
                return self._parent._element_constructor_(2)
            elif self.value == -2:
                return self._parent._element_constructor_(-3)
            else:
                return None

    cpdef Letter f(self, int i):
        r"""
        Return the action of `f_i` on ``self``.

        EXAMPLES::

            sage: C = crystals.Letters(['G',2])
            sage: [(c,i,c.f(i)) for i in C.index_set() for c in C if c.f(i) is not None]
            [(1, 1, 2),
             (3, 1, 0),
             (0, 1, -3),
             (-2, 1, -1),
             (2, 2, 3),
             (-3, 2, -2)]
        """
        if i == 1:
            if self.value == 1:
                return self._parent._element_constructor_(2)
            elif self.value == 3:
                return self._parent._element_constructor_(0)
            elif self.value == 0:
                return self._parent._element_constructor_(-3)
            elif self.value == -2:
                return self._parent._element_constructor_(-1)
            else:
                return None
        else:
            if self.value == 2:
                return self._parent._element_constructor_(3)
            elif self.value == -3:
                return self._parent._element_constructor_(-2)
            else:
                return None

    cpdef int epsilon(self, int i):
        r"""
        Return `\varepsilon_i` of ``self``.

        EXAMPLES::

            sage: C = crystals.Letters(['G',2])
            sage: [(c,i,c.epsilon(i)) for i in C.index_set() for c in C if c.epsilon(i) != 0]
            [(2, 1, 1), (0, 1, 1), (-3, 1, 2), (-1, 1, 1), (3, 2, 1), (-2, 2, 1)]
        """
        if i == 1:
            if self.value in (2,0,-1):
                return 1
            if self.value == -3:
                return 2
            return 0
        if self.value == 3 or self.value == -2: # i must be 2
            return 1
        return 0

    cpdef int phi(self, int i):
        r"""
        Return `\varphi_i` of ``self``.

        EXAMPLES::

            sage: C = crystals.Letters(['G',2])
            sage: [(c,i,c.phi(i)) for i in C.index_set() for c in C if c.phi(i) != 0]
            [(1, 1, 1), (3, 1, 2), (0, 1, 1), (-2, 1, 1), (2, 2, 1), (-3, 2, 1)]
        """
        if i == 1:
            if self.value in (1,0,-2):
                return 1
            if self.value == 3:
                return 2
            return 0
        if self.value == -3 or self.value == 2: # i must be 2
            return 1
        return 0

#########################
# Type E Letter
#########################

cdef class LetterTuple(Element):
    """
    Abstract class for type `E` letters.
    """
    cdef readonly tuple value

    def __init__(self, parent, tuple value):
        """
        Initialize ``self``.

        EXAMPLES::

            sage: C = crystals.Letters(['E',6])
            sage: a = C((1,-3))
            sage: TestSuite(a).run()
        """
        self.value = value
        Element.__init__(self, parent)

    def __setstate__(self, state):
        r"""
        Used in unpickling old pickles.

        EXAMPLES::

            sage: C = crystals.Letters(['E',6])
            sage: a = C((1,-3))
            sage: loads(dumps(a)) == a
            True
        """
        self._set_parent(state[0])
        self.value = tuple(state[1]['value'])

    def __reduce__(self):
        """
        Used in pickling of letters.

        EXAMPLES::

            sage: C = crystals.Letters(['E',6])
            sage: a = C((1,-3))
            sage: a.__reduce__()
            (The crystal of letters for type ['E', 6], ((1, -3),))
        """
        return (self._parent, (self.value,))

    def __hash__(self):
        """
        Return the hash value of ``self``.

        EXAMPLES::

            sage: C = crystals.Letters(['E', 6])
            sage: hash(C((1, -3))) == hash((1, -3))
            True
        """
        return hash(self.value)

    cpdef _richcmp_(left, Element right, int op):
        """
        Check comparison between ``left`` and ``right`` based on ``op``

        EXAMPLES::

            sage: C = crystals.Letters(['E', 6])
            sage: C((1,)) < C((-1, 3)) # indirect doctest
            True
            sage: C((6,)) < C((1,))
            False
            sage: C((-1, 3)) == C((-1, 3))
            True
        """
        cdef LetterTuple self, x
        self = left
        x = right
        if op == Py_EQ:
            return self.value == x.value
        if op == Py_NE:
            return self.value != x.value
        if op == Py_LT:
            return self._parent.lt_elements(self, x)
        if op == Py_GT:
            return x._parent.lt_elements(x, self)
        if op == Py_LE:
            return self.value == x.value or self._parent.lt_elements(self, x)
        if op == Py_GE:
            return self.value == x.value or x._parent.lt_elements(x, self)
        return False

    def _repr_(self):
        """
        Return a string representation of ``self``.

        EXAMPLES::

            sage: C = crystals.Letters(['E', 6])
            sage: C((-1, 3))
            (-1, 3)
        """
        return repr(self.value)

    def _latex_(self):
        r"""
        A latex representation of ``self``.

        EXAMPLES::

            sage: C = crystals.Letters(['E', 6])
            sage: latex(C((-1, 3)))
            \left(\overline{1}, 3\right)
        """
        ret = "\\left("
        first = True
        for v in self.value:
            if not first:
                ret += ", "
            else:
                first = False
            if v < 0:
                ret += "\\overline{" + repr(-v) + "}"
            else:
                ret+= repr(v)
        return ret + "\\right)"

    cpdef int epsilon(self, int i):
        r"""
        Return `\varepsilon_i` of ``self``.

        EXAMPLES::

            sage: C = crystals.Letters(['E',6])
            sage: C((-6,)).epsilon(1)
            0
            sage: C((-6,)).epsilon(6)
            1
        """
        if -i in self.value:
            return 1
        return 0

    cpdef int phi(self, int i):
        r"""
        Return `\varphi_i` of ``self``.

        EXAMPLES::

            sage: C = crystals.Letters(['E',6])
            sage: C((1,)).phi(1)
            1
            sage: C((1,)).phi(6)
            0
        """
        if i in self.value:
            return 1
        return 0

#########################
# Type E6
#########################

cdef class Crystal_of_letters_type_E6_element(LetterTuple):
    r"""
    Type `E_6` crystal of letters elements. This crystal corresponds to the highest weight
    crystal `B(\Lambda_1)`.

    TESTS::

        sage: C = crystals.Letters(['E',6])
        sage: C.module_generators
        ((1,),)
        sage: C.list()
        [(1,), (-1, 3), (-3, 4), (-4, 2, 5), (-2, 5), (-5, 2, 6), (-2, -5, 4, 6),
        (-4, 3, 6), (-3, 1, 6), (-1, 6), (-6, 2), (-2, -6, 4), (-4, -6, 3, 5),
        (-3, -6, 1, 5), (-1, -6, 5), (-5, 3), (-3, -5, 1, 4), (-1, -5, 4), (-4, 1, 2),
        (-1, -4, 2, 3), (-3, 2), (-2, -3, 4), (-4, 5), (-5, 6), (-6,), (-2, 1), (-1, -2, 3)]
        sage: TestSuite(C).run()
        sage: all(b.f(i).e(i) == b for i in C.index_set() for b in C if b.f(i) is not None)
        True
        sage: all(b.e(i).f(i) == b for i in C.index_set() for b in C if b.e(i) is not None)
        True
        sage: G = C.digraph()
        sage: G.show(edge_labels=true, figsize=12, vertex_size=1)
    """

    def _repr_(self):
        """
        In their full representation, the vertices of this crystal are labeled
        by their weight. For example vertex (-5,2,6) indicates that a 5-arrow
        is coming into this vertex, and a 2-arrow and 6-arrow is leaving the vertex.
        Specifying element_print_style = 'compact' for a given crystal C, labels the
        vertices of this crystal by the 27 letters +abcdefghijklmnopqrstuvwxyz.

        EXAMPLES::

            sage: C = crystals.Letters(['E',6], element_print_style = 'compact')
            sage: C.list()
            [+, a, b, c, d, e, f, g, h, i, j, k, l, m, n, o, p, q, r, s, t, u, v, w, x, y, z]
        """
        if self._parent._element_print_style == 'compact':
            l=['+','a','b','c','d','e','f','g','h','i','j','k','l','m','n','o','p','q','r','s','t','u','v','w','x','y','z']
            return l[self._parent.list().index(self)]
        return repr(self.value)

    def weight(self):
        """
        Return the weight of ``self``.

        EXAMPLES::

            sage: [v.weight() for v in crystals.Letters(['E',6])]
            [(0, 0, 0, 0, 0, -2/3, -2/3, 2/3),
             (-1/2, 1/2, 1/2, 1/2, 1/2, -1/6, -1/6, 1/6),
             (1/2, -1/2, 1/2, 1/2, 1/2, -1/6, -1/6, 1/6),
             (1/2, 1/2, -1/2, 1/2, 1/2, -1/6, -1/6, 1/6),
             (-1/2, -1/2, -1/2, 1/2, 1/2, -1/6, -1/6, 1/6),
             (1/2, 1/2, 1/2, -1/2, 1/2, -1/6, -1/6, 1/6),
             (-1/2, -1/2, 1/2, -1/2, 1/2, -1/6, -1/6, 1/6),
             (-1/2, 1/2, -1/2, -1/2, 1/2, -1/6, -1/6, 1/6),
             (1/2, -1/2, -1/2, -1/2, 1/2, -1/6, -1/6, 1/6),
             (0, 0, 0, 0, 1, 1/3, 1/3, -1/3),
             (1/2, 1/2, 1/2, 1/2, -1/2, -1/6, -1/6, 1/6),
             (-1/2, -1/2, 1/2, 1/2, -1/2, -1/6, -1/6, 1/6),
             (-1/2, 1/2, -1/2, 1/2, -1/2, -1/6, -1/6, 1/6),
             (1/2, -1/2, -1/2, 1/2, -1/2, -1/6, -1/6, 1/6),
             (0, 0, 0, 1, 0, 1/3, 1/3, -1/3),
             (-1/2, 1/2, 1/2, -1/2, -1/2, -1/6, -1/6, 1/6),
             (1/2, -1/2, 1/2, -1/2, -1/2, -1/6, -1/6, 1/6),
             (0, 0, 1, 0, 0, 1/3, 1/3, -1/3),
             (1/2, 1/2, -1/2, -1/2, -1/2, -1/6, -1/6, 1/6),
             (0, 1, 0, 0, 0, 1/3, 1/3, -1/3),
             (1, 0, 0, 0, 0, 1/3, 1/3, -1/3),
             (0, -1, 0, 0, 0, 1/3, 1/3, -1/3),
             (0, 0, -1, 0, 0, 1/3, 1/3, -1/3),
             (0, 0, 0, -1, 0, 1/3, 1/3, -1/3),
             (0, 0, 0, 0, -1, 1/3, 1/3, -1/3),
             (-1/2, -1/2, -1/2, -1/2, -1/2, -1/6, -1/6, 1/6),
             (-1, 0, 0, 0, 0, 1/3, 1/3, -1/3)]
        """
        R=self._parent.weight_lattice_realization().fundamental_weights()
        return sum(cmp(i,0)*R[abs(i)] for i in self.value)

    cpdef LetterTuple e(self, int i):
        r"""
        Return the action of `e_i` on ``self``.

        EXAMPLES::

            sage: C = crystals.Letters(['E',6])
            sage: C((-1,3)).e(1)
            (1,)
            sage: C((-2,-3,4)).e(2)
            (-3, 2)
            sage: C((1,)).e(1)
        """
        if self.value == (-1, 3) and i == 1:
            return self._parent._element_constructor_((1,))
        if self.value == (-3, 4) and i == 3:
            return self._parent._element_constructor_((-1, 3))
        if self.value == (-4, 2, 5) and i == 4:
            return self._parent._element_constructor_((-3, 4))
        if self.value == (-5, 2, 6) and i == 5:
            return self._parent._element_constructor_((-4, 2, 5))
        if self.value == (-2, 5) and i == 2:
            return self._parent._element_constructor_((-4, 2, 5))
        if self.value == (-6, 2) and i == 6:
            return self._parent._element_constructor_((-5, 2, 6))
        if self.value == (-2, -5, 4, 6) and i == 2:
            return self._parent._element_constructor_((-5, 2, 6))
        if self.value == (-2, -6, 4) and i == 2:
            return self._parent._element_constructor_((-6, 2))
        if self.value == (-2, -5, 4, 6) and i == 5:
            return self._parent._element_constructor_((-2, 5))
        if self.value == (-2, -6, 4) and i == 6:
            return self._parent._element_constructor_((-2, -5, 4, 6))
        if self.value == (-4, 3, 6) and i == 4:
            return self._parent._element_constructor_((-2, -5, 4, 6))
        if self.value == (-4, -6, 3, 5) and i == 4:
            return self._parent._element_constructor_((-2, -6, 4))
        if self.value == (-4, -6, 3, 5) and i == 6:
            return self._parent._element_constructor_((-4, 3, 6))
        if self.value == (-3, 1, 6) and i == 3:
            return self._parent._element_constructor_((-4, 3, 6))
        if self.value == (-5, 3) and i == 5:
            return self._parent._element_constructor_((-4, -6, 3, 5))
        if self.value == (-3, -6, 1, 5) and i == 3:
            return self._parent._element_constructor_((-4, -6, 3, 5))
        if self.value == (-3, -5, 1, 4) and i == 3:
            return self._parent._element_constructor_((-5, 3))
        if self.value == (-3, -6, 1, 5) and i == 6:
            return self._parent._element_constructor_((-3, 1, 6))
        if self.value == (-1, 6) and i == 1:
            return self._parent._element_constructor_((-3, 1, 6))
        if self.value == (-3, -5, 1, 4) and i == 5:
            return self._parent._element_constructor_((-3, -6, 1, 5))
        if self.value == (-1, -6, 5) and i == 1:
            return self._parent._element_constructor_((-3, -6, 1, 5))
        if self.value == (-4, 1, 2) and i == 4:
            return self._parent._element_constructor_((-3, -5, 1, 4))
        if self.value == (-1, -5, 4) and i == 1:
            return self._parent._element_constructor_((-3, -5, 1, 4))
        if self.value == (-2, 1) and i == 2:
            return self._parent._element_constructor_((-4, 1, 2))
        if self.value == (-1, -4, 2, 3) and i == 1:
            return self._parent._element_constructor_((-4, 1, 2))
        if self.value == (-1, -2, 3) and i == 1:
            return self._parent._element_constructor_((-2, 1))
        if self.value == (-1, -6, 5) and i == 6:
            return self._parent._element_constructor_((-1, 6))
        if self.value == (-1, -5, 4) and i == 5:
            return self._parent._element_constructor_((-1, -6, 5))
        if self.value == (-1, -4, 2, 3) and i == 4:
            return self._parent._element_constructor_((-1, -5, 4))
        if self.value == (-1, -2, 3) and i == 2:
            return self._parent._element_constructor_((-1, -4, 2, 3))
        if self.value == (-3, 2) and i == 3:
            return self._parent._element_constructor_((-1, -4, 2, 3))
        if self.value == (-2, -3, 4) and i == 3:
            return self._parent._element_constructor_((-1, -2, 3))
        if self.value == (-2, -3, 4) and i == 2:
            return self._parent._element_constructor_((-3, 2))
        if self.value == (-4, 5) and i == 4:
            return self._parent._element_constructor_((-2, -3, 4))
        if self.value == (-5, 6) and i == 5:
            return self._parent._element_constructor_((-4, 5))
        if self.value == (-6,) and i == 6:
            return self._parent._element_constructor_((-5, 6))
        else:
            return None

    cpdef LetterTuple f(self, int i):
        r"""
        Return the action of `f_i` on ``self``.

        EXAMPLES::

            sage: C = crystals.Letters(['E',6])
            sage: C((1,)).f(1)
            (-1, 3)
            sage: C((-6,)).f(1)
        """
        if self.value == (1,) and i == 1:
            return self._parent._element_constructor_((-1, 3))
        if self.value == (-1, 3) and i == 3:
            return self._parent._element_constructor_((-3, 4))
        if self.value == (-3, 4) and i == 4:
            return self._parent._element_constructor_((-4, 2, 5))
        if self.value == (-4, 2, 5) and i == 5:
            return self._parent._element_constructor_((-5, 2, 6))
        if self.value == (-4, 2, 5) and i == 2:
            return self._parent._element_constructor_((-2, 5))
        if self.value == (-5, 2, 6) and i == 6:
            return self._parent._element_constructor_((-6, 2))
        if self.value == (-5, 2, 6) and i == 2:
            return self._parent._element_constructor_((-2, -5, 4, 6))
        if self.value == (-6, 2) and i == 2:
            return self._parent._element_constructor_((-2, -6, 4))
        if self.value == (-2, 5) and i == 5:
            return self._parent._element_constructor_((-2, -5, 4, 6))
        if self.value == (-2, -5, 4, 6) and i == 6:
            return self._parent._element_constructor_((-2, -6, 4))
        if self.value == (-2, -5, 4, 6) and i == 4:
            return self._parent._element_constructor_((-4, 3, 6))
        if self.value == (-2, -6, 4) and i == 4:
            return self._parent._element_constructor_((-4, -6, 3, 5))
        if self.value == (-4, 3, 6) and i == 6:
            return self._parent._element_constructor_((-4, -6, 3, 5))
        if self.value == (-4, 3, 6) and i == 3:
            return self._parent._element_constructor_((-3, 1, 6))
        if self.value == (-4, -6, 3, 5) and i == 5:
            return self._parent._element_constructor_((-5, 3))
        if self.value == (-4, -6, 3, 5) and i == 3:
            return self._parent._element_constructor_((-3, -6, 1, 5))
        if self.value == (-5, 3) and i == 3:
            return self._parent._element_constructor_((-3, -5, 1, 4))
        if self.value == (-3, 1, 6) and i == 6:
            return self._parent._element_constructor_((-3, -6, 1, 5))
        if self.value == (-3, 1, 6) and i == 1:
            return self._parent._element_constructor_((-1, 6))
        if self.value == (-3, -6, 1, 5) and i == 5:
            return self._parent._element_constructor_((-3, -5, 1, 4))
        if self.value == (-3, -6, 1, 5) and i == 1:
            return self._parent._element_constructor_((-1, -6, 5))
        if self.value == (-3, -5, 1, 4) and i == 4:
            return self._parent._element_constructor_((-4, 1, 2))
        if self.value == (-3, -5, 1, 4) and i == 1:
            return self._parent._element_constructor_((-1, -5, 4))
        if self.value == (-4, 1, 2) and i == 2:
            return self._parent._element_constructor_((-2, 1))
        if self.value == (-4, 1, 2) and i == 1:
            return self._parent._element_constructor_((-1, -4, 2, 3))
        if self.value == (-2, 1) and i == 1:
            return self._parent._element_constructor_((-1, -2, 3))
        if self.value == (-1, 6) and i == 6:
            return self._parent._element_constructor_((-1, -6, 5))
        if self.value == (-1, -6, 5) and i == 5:
            return self._parent._element_constructor_((-1, -5, 4))
        if self.value == (-1, -5, 4) and i == 4:
            return self._parent._element_constructor_((-1, -4, 2, 3))
        if self.value == (-1, -4, 2, 3) and i == 2:
            return self._parent._element_constructor_((-1, -2, 3))
        if self.value == (-1, -4, 2, 3) and i == 3:
            return self._parent._element_constructor_((-3, 2))
        if self.value == (-1, -2, 3) and i == 3:
            return self._parent._element_constructor_((-2, -3, 4))
        if self.value == (-3, 2) and i == 2:
            return self._parent._element_constructor_((-2, -3, 4))
        if self.value == (-2, -3, 4) and i == 4:
            return self._parent._element_constructor_((-4, 5))
        if self.value == (-4, 5) and i == 5:
            return self._parent._element_constructor_((-5, 6))
        if self.value == (-5, 6) and i == 6:
            return self._parent._element_constructor_((-6,))
        else:
            return None

cdef class Crystal_of_letters_type_E6_element_dual(LetterTuple):
    r"""
    Type `E_6` crystal of letters elements. This crystal corresponds to the highest weight
    crystal `B(\Lambda_6)`. This crystal is dual to `B(\Lambda_1)` of type `E_6`.

    TESTS::

        sage: C = crystals.Letters(['E',6], dual = True)
        sage: C.module_generators
        ((6,),)
        sage: all(b==b.retract(b.lift()) for b in C)
        True
        sage: C.list()
        [(6,), (5, -6), (4, -5), (2, 3, -4), (3, -2), (1, 2, -3), (2, -1), (1, 4, -2, -3),
         (4, -1, -2), (1, 5, -4), (3, 5, -1, -4), (5, -3), (1, 6, -5), (3, 6, -1, -5), (4, 6, -3, -5),
         (2, 6, -4), (6, -2), (1, -6), (3, -1, -6), (4, -3, -6), (2, 5, -4, -6), (5, -2, -6), (2, -5),
         (4, -2, -5), (3, -4), (1, -3), (-1,)]
        sage: TestSuite(C).run()
        sage: all(b.f(i).e(i) == b for i in C.index_set() for b in C if b.f(i) is not None)
        True
        sage: all(b.e(i).f(i) == b for i in C.index_set() for b in C if b.e(i) is not None)
        True
        sage: G = C.digraph()
        sage: G.show(edge_labels=true, figsize=12, vertex_size=1)
    """

    def _repr_(self):
        """
        In their full representation, the vertices of this crystal are labeled
        by their weight. For example vertex (-2,1) indicates that a 2-arrow
        is coming into this vertex, and a 1-arrow is leaving the vertex.
        Specifying the option element_print_style = 'compact' for a given crystal C,
        labels the vertices of this crystal by the 27 letters -ABCDEFGHIJKLMNOPQRSTUVWXYZ

        EXAMPLES::

            sage: C = crystals.Letters(['E',6], element_print_style = 'compact', dual = True)
            sage: C.list()
            [-, A, B, C, D, E, F, G, H, I, J, K, L, M, N, O, P, Q, R, S, T, U, V, W, X, Y, Z]
            """
        if self._parent._element_print_style == 'compact':
            l=['-','A','B','C','D','E','F','G','H','I','J','K','L','M','N','O','P','Q','R','S','T','U','V','W','X','Y','Z']
            return l[self._parent.list().index(self)]
        return repr(self.value)

    cpdef LetterTuple lift(self):
        """
        Lift an element of ``self`` to the crystal of letters
        ``crystals.Letters(['E',6])`` by taking its inverse weight.

        EXAMPLES::

            sage: C = crystals.Letters(['E',6], dual = True)
            sage: b = C.module_generators[0]
            sage: b.lift()
            (-6,)
        """
        # Because a generators are not supported and the element constuctor
        #  being a cached method can't take lists as input, we have to make a
        #  tuple from a list
        return self._parent._ambient(tuple([-i for i in self.value]))

    cpdef LetterTuple retract(self, LetterTuple p):
        """
        Retract element ``p``, which is an element in
        ``crystals.Letters(['E',6])`` to an element in
        ``crystals.Letters(['E',6], dual=True)`` by taking its inverse weight.

        EXAMPLES::

            sage: C = crystals.Letters(['E',6])
            sage: Cd = crystals.Letters(['E',6], dual = True)
            sage: b = Cd.module_generators[0]
            sage: p = C((-1,3))
            sage: b.retract(p)
            (1, -3)
            sage: b.retract(None)
        """
        if p is None:
            return None
        # Because a generators are not supported and the element constuctor
        #  being a cached method can't take lists as input, we have to make a
        #  tuple from a list
        return self._parent._element_constructor_(tuple([-i for i in p.value]))

    cpdef LetterTuple e(self, int i):
        r"""
        Return the action of `e_i` on ``self``.

        EXAMPLES::

            sage: C = crystals.Letters(['E',6], dual = True)
            sage: C((-1,)).e(1)
            (1, -3)
        """
        return self.retract(self.lift().f(i))

    cpdef LetterTuple f(self, int i):
        r"""
        Return the action of `f_i` on ``self``.

        EXAMPLES::

            sage: C = crystals.Letters(['E',6], dual = True)
            sage: C((6,)).f(6)
            (5, -6)
            sage: C((6,)).f(1)
        """
        return self.retract(self.lift().e(i))

    def weight(self):
        """
        Return the weight of ``self``.

        EXAMPLES::

            sage: C = crystals.Letters(['E',6], dual = True)
            sage: b=C.module_generators[0]
            sage: b.weight()
            (0, 0, 0, 0, 1, -1/3, -1/3, 1/3)
            sage: [v.weight() for v in C]
            [(0, 0, 0, 0, 1, -1/3, -1/3, 1/3),
            (0, 0, 0, 1, 0, -1/3, -1/3, 1/3),
            (0, 0, 1, 0, 0, -1/3, -1/3, 1/3),
            (0, 1, 0, 0, 0, -1/3, -1/3, 1/3),
            (-1, 0, 0, 0, 0, -1/3, -1/3, 1/3),
            (1, 0, 0, 0, 0, -1/3, -1/3, 1/3),
            (1/2, 1/2, 1/2, 1/2, 1/2, 1/6, 1/6, -1/6),
            (0, -1, 0, 0, 0, -1/3, -1/3, 1/3),
            (-1/2, -1/2, 1/2, 1/2, 1/2, 1/6, 1/6, -1/6),
            (0, 0, -1, 0, 0, -1/3, -1/3, 1/3),
            (-1/2, 1/2, -1/2, 1/2, 1/2, 1/6, 1/6, -1/6),
            (1/2, -1/2, -1/2, 1/2, 1/2, 1/6, 1/6, -1/6),
            (0, 0, 0, -1, 0, -1/3, -1/3, 1/3),
            (-1/2, 1/2, 1/2, -1/2, 1/2, 1/6, 1/6, -1/6),
            (1/2, -1/2, 1/2, -1/2, 1/2, 1/6, 1/6, -1/6),
            (1/2, 1/2, -1/2, -1/2, 1/2, 1/6, 1/6, -1/6),
            (-1/2, -1/2, -1/2, -1/2, 1/2, 1/6, 1/6, -1/6),
            (0, 0, 0, 0, -1, -1/3, -1/3, 1/3),
            (-1/2, 1/2, 1/2, 1/2, -1/2, 1/6, 1/6, -1/6),
            (1/2, -1/2, 1/2, 1/2, -1/2, 1/6, 1/6, -1/6),
            (1/2, 1/2, -1/2, 1/2, -1/2, 1/6, 1/6, -1/6),
            (-1/2, -1/2, -1/2, 1/2, -1/2, 1/6, 1/6, -1/6),
            (1/2, 1/2, 1/2, -1/2, -1/2, 1/6, 1/6, -1/6),
            (-1/2, -1/2, 1/2, -1/2, -1/2, 1/6, 1/6, -1/6),
            (-1/2, 1/2, -1/2, -1/2, -1/2, 1/6, 1/6, -1/6),
            (1/2, -1/2, -1/2, -1/2, -1/2, 1/6, 1/6, -1/6),
            (0, 0, 0, 0, 0, 2/3, 2/3, -2/3)]
        """
        return -self.lift().weight()


#########################
# Type E7
#########################

cdef class Crystal_of_letters_type_E7_element(LetterTuple):
    r"""
    Type `E_7` crystal of letters elements. This crystal corresponds to the highest weight
    crystal `B(\Lambda_7)`.

    TESTS::

        sage: C = crystals.Letters(['E',7])
        sage: C.module_generators
        ((7,),)
        sage: C.list()
        [(7,), (-7, 6), (-6, 5), (-5, 4), (-4, 2, 3), (-2, 3), (-3, 1, 2), (-1,
        2), (-3, -2, 1, 4), (-1, -2, 4), (-4, 1, 5), (-4, -1, 3, 5), (-3, 5),
        (-5, 6, 1), (-5, -1, 3, 6), (-5, -3, 4, 6), (-4, 2, 6), (-2, 6), (-6, 7,
        1), (-1, -6, 3, 7), (-6, -3, 7, 4), (-6, -4, 2, 7, 5), (-6, -2, 7, 5),
        (-5, 7, 2), (-5, -2, 4, 7), (-4, 7, 3), (-3, 1, 7), (-1, 7), (-7, 1),
        (-1, -7, 3), (-7, -3, 4), (-4, -7, 2, 5), (-7, -2, 5), (-5, -7, 6, 2),
        (-5, -2, -7, 4, 6), (-7, -4, 6, 3), (-3, -7, 1, 6), (-7, -1, 6), (-6,
        2), (-2, -6, 4), (-6, -4, 5, 3), (-3, -6, 1, 5), (-6, -1, 5), (-5, 3),
        (-3, -5, 4, 1), (-5, -1, 4), (-4, 1, 2), (-1, -4, 3, 2), (-3, 2), (-2,
        -3, 4), (-4, 5), (-5, 6), (-6, 7), (-7,), (-2, 1), (-2, -1, 3)]
        sage: TestSuite(C).run()
        sage: all(b.f(i).e(i) == b for i in C.index_set() for b in C if b.f(i) is not None)
        True
        sage: all(b.e(i).f(i) == b for i in C.index_set() for b in C if b.e(i) is not None)
        True
        sage: G = C.digraph()
        sage: G.show(edge_labels=true, figsize=12, vertex_size=1)
    """

    def weight(self):
        """
        Return the weight of ``self``.

        EXAMPLES::

            sage: [v.weight() for v in crystals.Letters(['E',7])]
            [(0, 0, 0, 0, 0, 1, -1/2, 1/2), (0, 0, 0, 0, 1, 0, -1/2, 1/2), (0, 0, 0,
            1, 0, 0, -1/2, 1/2), (0, 0, 1, 0, 0, 0, -1/2, 1/2), (0, 1, 0, 0, 0, 0,
            -1/2, 1/2), (-1, 0, 0, 0, 0, 0, -1/2, 1/2), (1, 0, 0, 0, 0, 0, -1/2,
            1/2), (1/2, 1/2, 1/2, 1/2, 1/2, 1/2, 0, 0), (0, -1, 0, 0, 0, 0, -1/2,
            1/2), (-1/2, -1/2, 1/2, 1/2, 1/2, 1/2, 0, 0), (0, 0, -1, 0, 0, 0, -1/2,
            1/2), (-1/2, 1/2, -1/2, 1/2, 1/2, 1/2, 0, 0), (1/2, -1/2, -1/2, 1/2,
            1/2, 1/2, 0, 0), (0, 0, 0, -1, 0, 0, -1/2, 1/2), (-1/2, 1/2, 1/2, -1/2,
            1/2, 1/2, 0, 0), (1/2, -1/2, 1/2, -1/2, 1/2, 1/2, 0, 0), (1/2, 1/2,
            -1/2, -1/2, 1/2, 1/2, 0, 0), (-1/2, -1/2, -1/2, -1/2, 1/2, 1/2, 0, 0),
            (0, 0, 0, 0, -1, 0, -1/2, 1/2), (-1/2, 1/2, 1/2, 1/2, -1/2, 1/2, 0, 0),
            (1/2, -1/2, 1/2, 1/2, -1/2, 1/2, 0, 0), (1/2, 1/2, -1/2, 1/2, -1/2, 1/2,
            0, 0), (-1/2, -1/2, -1/2, 1/2, -1/2, 1/2, 0, 0), (1/2, 1/2, 1/2, -1/2,
            -1/2, 1/2, 0, 0), (-1/2, -1/2, 1/2, -1/2, -1/2, 1/2, 0, 0), (-1/2, 1/2,
            -1/2, -1/2, -1/2, 1/2, 0, 0), (1/2, -1/2, -1/2, -1/2, -1/2, 1/2, 0, 0),
            (0, 0, 0, 0, 0, 1, 1/2, -1/2), (0, 0, 0, 0, 0, -1, -1/2, 1/2), (-1/2,
            1/2, 1/2, 1/2, 1/2, -1/2, 0, 0), (1/2, -1/2, 1/2, 1/2, 1/2, -1/2, 0, 0),
            (1/2, 1/2, -1/2, 1/2, 1/2, -1/2, 0, 0), (-1/2, -1/2, -1/2, 1/2, 1/2,
            -1/2, 0, 0), (1/2, 1/2, 1/2, -1/2, 1/2, -1/2, 0, 0), (-1/2, -1/2, 1/2,
            -1/2, 1/2, -1/2, 0, 0), (-1/2, 1/2, -1/2, -1/2, 1/2, -1/2, 0, 0), (1/2,
            -1/2, -1/2, -1/2, 1/2, -1/2, 0, 0), (0, 0, 0, 0, 1, 0, 1/2, -1/2), (1/2,
            1/2, 1/2, 1/2, -1/2, -1/2, 0, 0), (-1/2, -1/2, 1/2, 1/2, -1/2, -1/2, 0,
            0), (-1/2, 1/2, -1/2, 1/2, -1/2, -1/2, 0, 0), (1/2, -1/2, -1/2, 1/2,
            -1/2, -1/2, 0, 0), (0, 0, 0, 1, 0, 0, 1/2, -1/2), (-1/2, 1/2, 1/2, -1/2,
            -1/2, -1/2, 0, 0), (1/2, -1/2, 1/2, -1/2, -1/2, -1/2, 0, 0), (0, 0, 1,
            0, 0, 0, 1/2, -1/2), (1/2, 1/2, -1/2, -1/2, -1/2, -1/2, 0, 0), (0, 1, 0,
            0, 0, 0, 1/2, -1/2), (1, 0, 0, 0, 0, 0, 1/2, -1/2), (0, -1, 0, 0, 0, 0,
            1/2, -1/2), (0, 0, -1, 0, 0, 0, 1/2, -1/2), (0, 0, 0, -1, 0, 0, 1/2,
            -1/2), (0, 0, 0, 0, -1, 0, 1/2, -1/2), (0, 0, 0, 0, 0, -1, 1/2, -1/2),
            (-1/2, -1/2, -1/2, -1/2, -1/2, -1/2, 0, 0), (-1, 0, 0, 0, 0, 0, 1/2,
            -1/2)]
        """
        R=self._parent.weight_lattice_realization().fundamental_weights()
        return sum(cmp(i,0)*R[abs(i)] for i in self.value)

    cpdef LetterTuple e(self, int i):
        r"""
        Return the action of `e_i` on ``self``.

        EXAMPLES::

            sage: C = crystals.Letters(['E',7])
            sage: C((7,)).e(7)
            sage: C((-7,6)).e(7)
            (7,)
        """
        if self.value ==  (-7, 6)  and i ==  7 :
            return self._parent._element_constructor_( (7,) )
        if self.value ==  (-6, 5)  and i ==  6 :
            return self._parent._element_constructor_( (-7, 6) )
        if self.value ==  (-5, 4)  and i ==  5 :
            return self._parent._element_constructor_( (-6, 5) )
        if self.value ==  (-4, 2, 3)  and i ==  4 :
            return self._parent._element_constructor_( (-5, 4) )
        if self.value ==  (-2, 3)  and i ==  2 :
            return self._parent._element_constructor_( (-4, 2, 3) )
        if self.value ==  (-3, 1, 2)  and i ==  3 :
            return self._parent._element_constructor_( (-4, 2, 3) )
        if self.value ==  (-3, -2, 1, 4)  and i ==  3 :
            return self._parent._element_constructor_( (-2, 3) )
        if self.value ==  (-1, 2)  and i ==  1 :
            return self._parent._element_constructor_( (-3, 1, 2) )
        if self.value ==  (-3, -2, 1, 4)  and i ==  2 :
            return self._parent._element_constructor_( (-3, 1, 2) )
        if self.value ==  (-1, -2, 4)  and i ==  1 :
            return self._parent._element_constructor_( (-3, -2, 1, 4) )
        if self.value ==  (-4, 1, 5)  and i ==  4 :
            return self._parent._element_constructor_( (-3, -2, 1, 4) )
        if self.value ==  (-7, 1)  and i ==  7 :
            return self._parent._element_constructor_( (-6, 7, 1) )
        if self.value ==  (-1, -6, 3, 7)  and i ==  1 :
            return self._parent._element_constructor_( (-6, 7, 1) )
        if self.value ==  (-1, -2, 4)  and i ==  2 :
            return self._parent._element_constructor_( (-1, 2) )
        if self.value ==  (-4, -1, 3, 5)  and i ==  4 :
            return self._parent._element_constructor_( (-1, -2, 4) )
        if self.value ==  (-4, -1, 3, 5)  and i ==  1 :
            return self._parent._element_constructor_( (-4, 1, 5) )
        if self.value ==  (-5, 6, 1)  and i ==  5 :
            return self._parent._element_constructor_( (-4, 1, 5) )
        if self.value ==  (-3, 5)  and i ==  3 :
            return self._parent._element_constructor_( (-4, -1, 3, 5) )
        if self.value ==  (-5, -1, 3, 6)  and i ==  5 :
            return self._parent._element_constructor_( (-4, -1, 3, 5) )
        if self.value ==  (-5, -3, 4, 6)  and i ==  5 :
            return self._parent._element_constructor_( (-3, 5) )
        if self.value ==  (-6, 7, 1)  and i ==  6 :
            return self._parent._element_constructor_( (-5, 6, 1) )
        if self.value ==  (-5, -1, 3, 6)  and i ==  1 :
            return self._parent._element_constructor_( (-5, 6, 1) )
        if self.value ==  (-5, -3, 4, 6)  and i ==  3 :
            return self._parent._element_constructor_( (-5, -1, 3, 6) )
        if self.value ==  (-1, -6, 3, 7)  and i ==  6 :
            return self._parent._element_constructor_( (-5, -1, 3, 6) )
        if self.value ==  (-4, 2, 6)  and i ==  4 :
            return self._parent._element_constructor_( (-5, -3, 4, 6) )
        if self.value ==  (-6, -3, 7, 4)  and i ==  6 :
            return self._parent._element_constructor_( (-5, -3, 4, 6) )
        if self.value ==  (-6, -2, 7, 5)  and i ==  6 :
            return self._parent._element_constructor_( (-2, 6) )
        if self.value ==  (-6, -3, 7, 4)  and i ==  3 :
            return self._parent._element_constructor_( (-1, -6, 3, 7) )
        if self.value ==  (-1, -7, 3)  and i ==  7 :
            return self._parent._element_constructor_( (-1, -6, 3, 7) )
        if self.value ==  (-7, -3, 4)  and i ==  7 :
            return self._parent._element_constructor_( (-6, -3, 7, 4) )
        if self.value ==  (-6, -4, 2, 7, 5)  and i ==  4 :
            return self._parent._element_constructor_( (-6, -3, 7, 4) )
        if self.value ==  (-2, 6)  and i ==  2 :
            return self._parent._element_constructor_( (-4, 2, 6) )
        if self.value ==  (-6, -4, 2, 7, 5)  and i ==  6 :
            return self._parent._element_constructor_( (-4, 2, 6) )
        if self.value ==  (-6, -2, 7, 5)  and i ==  2 :
            return self._parent._element_constructor_( (-6, -4, 2, 7, 5) )
        if self.value ==  (-4, -7, 2, 5)  and i ==  7 :
            return self._parent._element_constructor_( (-6, -4, 2, 7, 5) )
        if self.value ==  (-7, -4, 6, 3)  and i ==  7 :
            return self._parent._element_constructor_( (-4, 7, 3) )
        if self.value ==  (-3, 1, 7)  and i ==  3 :
            return self._parent._element_constructor_( (-4, 7, 3) )
        if self.value ==  (-1, 7)  and i ==  1 :
            return self._parent._element_constructor_( (-3, 1, 7) )
        if self.value ==  (-3, -7, 1, 6)  and i ==  7 :
            return self._parent._element_constructor_( (-3, 1, 7) )
        if self.value ==  (-1, -7, 3)  and i ==  1 :
            return self._parent._element_constructor_( (-7, 1) )
        if self.value ==  (-7, -2, 5)  and i ==  2 :
            return self._parent._element_constructor_( (-4, -7, 2, 5) )
        if self.value ==  (-5, -7, 6, 2)  and i ==  5 :
            return self._parent._element_constructor_( (-4, -7, 2, 5) )
        if self.value ==  (-5, -2, -7, 4, 6)  and i ==  5 :
            return self._parent._element_constructor_( (-7, -2, 5) )
        if self.value ==  (-5, -7, 6, 2)  and i ==  7 :
            return self._parent._element_constructor_( (-5, 7, 2) )
        if self.value ==  (-5, -2, 4, 7)  and i ==  2 :
            return self._parent._element_constructor_( (-5, 7, 2) )
        if self.value ==  (-7, -3, 4)  and i ==  3 :
            return self._parent._element_constructor_( (-1, -7, 3) )
        if self.value ==  (-5, 7, 2)  and i ==  5 :
            return self._parent._element_constructor_( (-6, -4, 2, 7, 5) )
        if self.value ==  (-6, 2)  and i ==  6 :
            return self._parent._element_constructor_( (-5, -7, 6, 2) )
        if self.value ==  (-5, -2, -7, 4, 6)  and i ==  2 :
            return self._parent._element_constructor_( (-5, -7, 6, 2) )
        if self.value ==  (-7, -2, 5)  and i ==  7 :
            return self._parent._element_constructor_( (-6, -2, 7, 5) )
        if self.value ==  (-5, -2, 4, 7)  and i ==  5 :
            return self._parent._element_constructor_( (-6, -2, 7, 5) )
        if self.value ==  (-4, 7, 3)  and i ==  4 :
            return self._parent._element_constructor_( (-5, -2, 4, 7) )
        if self.value ==  (-5, -2, -7, 4, 6)  and i ==  7 :
            return self._parent._element_constructor_( (-5, -2, 4, 7) )
        if self.value ==  (-4, -7, 2, 5)  and i ==  4 :
            return self._parent._element_constructor_( (-7, -3, 4) )
        if self.value ==  (-7, -4, 6, 3)  and i ==  4 :
            return self._parent._element_constructor_( (-5, -2, -7, 4, 6) )
        if self.value ==  (-2, -6, 4)  and i ==  6 :
            return self._parent._element_constructor_( (-5, -2, -7, 4, 6) )
        if self.value ==  (-6, -4, 5, 3)  and i ==  6 :
            return self._parent._element_constructor_( (-7, -4, 6, 3) )
        if self.value ==  (-3, -7, 1, 6)  and i ==  3 :
            return self._parent._element_constructor_( (-7, -4, 6, 3) )
        if self.value ==  (-3, -6, 1, 5)  and i ==  6 :
            return self._parent._element_constructor_( (-3, -7, 1, 6) )
        if self.value ==  (-6, -1, 5)  and i ==  6 :
            return self._parent._element_constructor_( (-7, -1, 6) )
        if self.value ==  (-2, -6, 4)  and i ==  2 :
            return self._parent._element_constructor_( (-6, 2) )
        if self.value ==  (-6, -4, 5, 3)  and i ==  4 :
            return self._parent._element_constructor_( (-2, -6, 4) )
        if self.value ==  (-7, -1, 6)  and i ==  1 :
            return self._parent._element_constructor_( (-3, -7, 1, 6) )
        if self.value ==  (-5, 3)  and i ==  5 :
            return self._parent._element_constructor_( (-6, -4, 5, 3) )
        if self.value ==  (-3, -6, 1, 5)  and i ==  3 :
            return self._parent._element_constructor_( (-6, -4, 5, 3) )
        if self.value ==  (-6, -1, 5)  and i ==  1 :
            return self._parent._element_constructor_( (-3, -6, 1, 5) )
        if self.value ==  (-3, -5, 4, 1)  and i ==  5 :
            return self._parent._element_constructor_( (-3, -6, 1, 5) )
        if self.value ==  (-5, -1, 4)  and i ==  5 :
            return self._parent._element_constructor_( (-6, -1, 5) )
        if self.value ==  (-3, -5, 4, 1)  and i ==  3 :
            return self._parent._element_constructor_( (-5, 3) )
        if self.value ==  (-4, 1, 2)  and i ==  4 :
            return self._parent._element_constructor_( (-3, -5, 4, 1) )
        if self.value ==  (-5, -1, 4)  and i ==  1 :
            return self._parent._element_constructor_( (-3, -5, 4, 1) )
        if self.value ==  (-1, -4, 3, 2)  and i ==  4 :
            return self._parent._element_constructor_( (-5, -1, 4) )
        if self.value ==  (-1, -4, 3, 2)  and i ==  1 :
            return self._parent._element_constructor_( (-4, 1, 2) )
        if self.value ==  (-2, 1)  and i ==  2 :
            return self._parent._element_constructor_( (-4, 1, 2) )
        if self.value ==  (-3, 2)  and i ==  3 :
            return self._parent._element_constructor_( (-1, -4, 3, 2) )
        if self.value ==  (-2, -1, 3)  and i ==  2 :
            return self._parent._element_constructor_( (-1, -4, 3, 2) )
        if self.value ==  (-2, -1, 3)  and i ==  1 :
            return self._parent._element_constructor_( (-2, 1) )
        if self.value ==  (-7, -1, 6)  and i ==  7 :
            return self._parent._element_constructor_( (-1, 7) )
        if self.value ==  (-2, -3, 4)  and i ==  3 :
            return self._parent._element_constructor_( (-2, -1, 3) )
        if self.value ==  (-2, -3, 4)  and i ==  2 :
            return self._parent._element_constructor_( (-3, 2) )
        if self.value ==  (-4, 5)  and i ==  4 :
            return self._parent._element_constructor_( (-2, -3, 4) )
        if self.value ==  (-5, 6)  and i ==  5 :
            return self._parent._element_constructor_( (-4, 5) )
        if self.value ==  (-6, 7)  and i ==  6 :
            return self._parent._element_constructor_( (-5, 6) )
        if self.value ==  (-7,)  and i ==  7 :
            return self._parent._element_constructor_( (-6, 7) )
        else:
            return None

    cpdef LetterTuple f(self, int i):
        r"""
        Return the action of `f_i` on ``self``.

        EXAMPLES::

            sage: C = crystals.Letters(['E',7])
            sage: C((-7,)).f(7)
            sage: C((7,)).f(7)
            (-7, 6)
        """
        if self.value ==  (7,)  and i ==  7 :
            return self._parent._element_constructor_( (-7, 6) )
        if self.value ==  (-7, 6)  and i ==  6 :
            return self._parent._element_constructor_( (-6, 5) )
        if self.value ==  (-6, 5)  and i ==  5 :
            return self._parent._element_constructor_( (-5, 4) )
        if self.value ==  (-5, 4)  and i ==  4 :
            return self._parent._element_constructor_( (-4, 2, 3) )
        if self.value ==  (-4, 2, 3)  and i ==  2 :
            return self._parent._element_constructor_( (-2, 3) )
        if self.value ==  (-4, 2, 3)  and i ==  3 :
            return self._parent._element_constructor_( (-3, 1, 2) )
        if self.value ==  (-2, 3)  and i ==  3 :
            return self._parent._element_constructor_( (-3, -2, 1, 4) )
        if self.value ==  (-3, 1, 2)  and i ==  1 :
            return self._parent._element_constructor_( (-1, 2) )
        if self.value ==  (-3, 1, 2)  and i ==  2 :
            return self._parent._element_constructor_( (-3, -2, 1, 4) )
        if self.value ==  (-3, -2, 1, 4)  and i ==  1 :
            return self._parent._element_constructor_( (-1, -2, 4) )
        if self.value ==  (-3, -2, 1, 4)  and i ==  4 :
            return self._parent._element_constructor_( (-4, 1, 5) )
        if self.value ==  (-6, 7, 1)  and i ==  7 :
            return self._parent._element_constructor_( (-7, 1) )
        if self.value ==  (-6, 7, 1)  and i ==  1 :
            return self._parent._element_constructor_( (-1, -6, 3, 7) )
        if self.value ==  (-1, 2)  and i ==  2 :
            return self._parent._element_constructor_( (-1, -2, 4) )
        if self.value ==  (-1, -2, 4)  and i ==  4 :
            return self._parent._element_constructor_( (-4, -1, 3, 5) )
        if self.value ==  (-4, 1, 5)  and i ==  1 :
            return self._parent._element_constructor_( (-4, -1, 3, 5) )
        if self.value ==  (-4, 1, 5)  and i ==  5 :
            return self._parent._element_constructor_( (-5, 6, 1) )
        if self.value ==  (-4, -1, 3, 5)  and i ==  3 :
            return self._parent._element_constructor_( (-3, 5) )
        if self.value ==  (-4, -1, 3, 5)  and i ==  5 :
            return self._parent._element_constructor_( (-5, -1, 3, 6) )
        if self.value ==  (-3, 5)  and i ==  5 :
            return self._parent._element_constructor_( (-5, -3, 4, 6) )
        if self.value ==  (-5, 6, 1)  and i ==  6 :
            return self._parent._element_constructor_( (-6, 7, 1) )
        if self.value ==  (-5, 6, 1)  and i ==  1 :
            return self._parent._element_constructor_( (-5, -1, 3, 6) )
        if self.value ==  (-5, -1, 3, 6)  and i ==  3 :
            return self._parent._element_constructor_( (-5, -3, 4, 6) )
        if self.value ==  (-5, -1, 3, 6)  and i ==  6 :
            return self._parent._element_constructor_( (-1, -6, 3, 7) )
        if self.value ==  (-5, -3, 4, 6)  and i ==  4 :
            return self._parent._element_constructor_( (-4, 2, 6) )
        if self.value ==  (-5, -3, 4, 6)  and i ==  6 :
            return self._parent._element_constructor_( (-6, -3, 7, 4) )
        if self.value ==  (-2, 6)  and i ==  6 :
            return self._parent._element_constructor_( (-6, -2, 7, 5) )
        if self.value ==  (-1, -6, 3, 7)  and i ==  3 :
            return self._parent._element_constructor_( (-6, -3, 7, 4) )
        if self.value ==  (-1, -6, 3, 7)  and i ==  7 :
            return self._parent._element_constructor_( (-1, -7, 3) )
        if self.value ==  (-6, -3, 7, 4)  and i ==  7 :
            return self._parent._element_constructor_( (-7, -3, 4) )
        if self.value ==  (-6, -3, 7, 4)  and i ==  4 :
            return self._parent._element_constructor_( (-6, -4, 2, 7, 5) )
        if self.value ==  (-4, 2, 6)  and i ==  2 :
            return self._parent._element_constructor_( (-2, 6) )
        if self.value ==  (-4, 2, 6)  and i ==  6 :
            return self._parent._element_constructor_( (-6, -4, 2, 7, 5) )
        if self.value ==  (-6, -4, 2, 7, 5)  and i ==  2 :
            return self._parent._element_constructor_( (-6, -2, 7, 5) )
        if self.value ==  (-6, -4, 2, 7, 5)  and i ==  7 :
            return self._parent._element_constructor_( (-4, -7, 2, 5) )
        if self.value ==  (-4, 7, 3)  and i ==  7 :
            return self._parent._element_constructor_( (-7, -4, 6, 3) )
        if self.value ==  (-4, 7, 3)  and i ==  3 :
            return self._parent._element_constructor_( (-3, 1, 7) )
        if self.value ==  (-3, 1, 7)  and i ==  1 :
            return self._parent._element_constructor_( (-1, 7) )
        if self.value ==  (-3, 1, 7)  and i ==  7 :
            return self._parent._element_constructor_( (-3, -7, 1, 6) )
        if self.value ==  (-7, 1)  and i ==  1 :
            return self._parent._element_constructor_( (-1, -7, 3) )
        if self.value ==  (-4, -7, 2, 5)  and i ==  2 :
            return self._parent._element_constructor_( (-7, -2, 5) )
        if self.value ==  (-4, -7, 2, 5)  and i ==  5 :
            return self._parent._element_constructor_( (-5, -7, 6, 2) )
        if self.value ==  (-7, -2, 5)  and i ==  5 :
            return self._parent._element_constructor_( (-5, -2, -7, 4, 6) )
        if self.value ==  (-5, 7, 2)  and i ==  7 :
            return self._parent._element_constructor_( (-5, -7, 6, 2) )
        if self.value ==  (-5, 7, 2)  and i ==  2 :
            return self._parent._element_constructor_( (-5, -2, 4, 7) )
        if self.value ==  (-1, -7, 3)  and i ==  3 :
            return self._parent._element_constructor_( (-7, -3, 4) )
        if self.value ==  (-6, -4, 2, 7, 5)  and i ==  5 :
            return self._parent._element_constructor_( (-5, 7, 2) )
        if self.value ==  (-5, -7, 6, 2)  and i ==  6 :
            return self._parent._element_constructor_( (-6, 2) )
        if self.value ==  (-5, -7, 6, 2)  and i ==  2 :
            return self._parent._element_constructor_( (-5, -2, -7, 4, 6) )
        if self.value ==  (-6, -2, 7, 5)  and i ==  7 :
            return self._parent._element_constructor_( (-7, -2, 5) )
        if self.value ==  (-6, -2, 7, 5)  and i ==  5 :
            return self._parent._element_constructor_( (-5, -2, 4, 7) )
        if self.value ==  (-5, -2, 4, 7)  and i ==  4 :
            return self._parent._element_constructor_( (-4, 7, 3) )
        if self.value ==  (-5, -2, 4, 7)  and i ==  7 :
            return self._parent._element_constructor_( (-5, -2, -7, 4, 6) )
        if self.value ==  (-7, -3, 4)  and i ==  4 :
            return self._parent._element_constructor_( (-4, -7, 2, 5) )
        if self.value ==  (-5, -2, -7, 4, 6)  and i ==  4 :
            return self._parent._element_constructor_( (-7, -4, 6, 3) )
        if self.value ==  (-5, -2, -7, 4, 6)  and i ==  6 :
            return self._parent._element_constructor_( (-2, -6, 4) )
        if self.value ==  (-7, -4, 6, 3)  and i ==  6 :
            return self._parent._element_constructor_( (-6, -4, 5, 3) )
        if self.value ==  (-7, -4, 6, 3)  and i ==  3 :
            return self._parent._element_constructor_( (-3, -7, 1, 6) )
        if self.value ==  (-3, -7, 1, 6)  and i ==  6 :
            return self._parent._element_constructor_( (-3, -6, 1, 5) )
        if self.value ==  (-7, -1, 6)  and i ==  6 :
            return self._parent._element_constructor_( (-6, -1, 5) )
        if self.value ==  (-6, 2)  and i ==  2 :
            return self._parent._element_constructor_( (-2, -6, 4) )
        if self.value ==  (-2, -6, 4)  and i ==  4 :
            return self._parent._element_constructor_( (-6, -4, 5, 3) )
        if self.value ==  (-3, -7, 1, 6)  and i ==  1 :
            return self._parent._element_constructor_( (-7, -1, 6) )
        if self.value ==  (-6, -4, 5, 3)  and i ==  5 :
            return self._parent._element_constructor_( (-5, 3) )
        if self.value ==  (-6, -4, 5, 3)  and i ==  3 :
            return self._parent._element_constructor_( (-3, -6, 1, 5) )
        if self.value ==  (-3, -6, 1, 5)  and i ==  1 :
            return self._parent._element_constructor_( (-6, -1, 5) )
        if self.value ==  (-3, -6, 1, 5)  and i ==  5 :
            return self._parent._element_constructor_( (-3, -5, 4, 1) )
        if self.value ==  (-6, -1, 5)  and i ==  5 :
            return self._parent._element_constructor_( (-5, -1, 4) )
        if self.value ==  (-5, 3)  and i ==  3 :
            return self._parent._element_constructor_( (-3, -5, 4, 1) )
        if self.value ==  (-3, -5, 4, 1)  and i ==  4 :
            return self._parent._element_constructor_( (-4, 1, 2) )
        if self.value ==  (-3, -5, 4, 1)  and i ==  1 :
            return self._parent._element_constructor_( (-5, -1, 4) )
        if self.value ==  (-5, -1, 4)  and i ==  4 :
            return self._parent._element_constructor_( (-1, -4, 3, 2) )
        if self.value ==  (-4, 1, 2)  and i ==  1 :
            return self._parent._element_constructor_( (-1, -4, 3, 2) )
        if self.value ==  (-4, 1, 2)  and i ==  2 :
            return self._parent._element_constructor_( (-2, 1) )
        if self.value ==  (-1, -4, 3, 2)  and i ==  3 :
            return self._parent._element_constructor_( (-3, 2) )
        if self.value ==  (-1, -4, 3, 2)  and i ==  2 :
            return self._parent._element_constructor_( (-2, -1, 3) )
        if self.value ==  (-2, 1)  and i ==  1 :
            return self._parent._element_constructor_( (-2, -1, 3) )
        if self.value ==  (-1, 7)  and i ==  7 :
            return self._parent._element_constructor_( (-7, -1, 6) )
        if self.value ==  (-2, -1, 3)  and i ==  3 :
            return self._parent._element_constructor_( (-2, -3, 4) )
        if self.value ==  (-3, 2)  and i ==  2 :
            return self._parent._element_constructor_( (-2, -3, 4) )
        if self.value ==  (-2, -3, 4)  and i ==  4 :
            return self._parent._element_constructor_( (-4, 5) )
        if self.value ==  (-4, 5)  and i ==  5 :
            return self._parent._element_constructor_( (-5, 6) )
        if self.value ==  (-5, 6)  and i ==  6 :
            return self._parent._element_constructor_( (-6, 7) )
        if self.value ==  (-6, 7)  and i ==  7 :
            return self._parent._element_constructor_( (-7,) )
        else:
            return None

