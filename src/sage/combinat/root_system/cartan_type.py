r"""
Cartan types

Loosely speaking, Dynkin diagrams (or equivalently Cartan matrices)
are graphs which are used to classify root systems, Coxeter and Weyl
groups, Lie algebras, Lie groups, crystals, etc. up to an
isomorphism. *Cartan types* are a standard set of names for those
Dynkin diagrams.

        http://en.wikipedia.org/wiki/Dynkin_diagram

Let us consider for example, the Cartan type `A_4`::

    sage: T = CartanType(['A', 4])
    sage: T
    ['A', 4]

It is the name of the following Dynkin diagram::

    sage: DynkinDiagram(T)
    O---O---O---O
    1   2   3   4
    A4

Note: for convenience, the following shortcuts are available::

    sage: DynkinDiagram(['A',4])
    O---O---O---O
    1   2   3   4
    A4
    sage: DynkinDiagram('A4')
    O---O---O---O
    1   2   3   4
    A4
    sage: T.dynkin_diagram()
    O---O---O---O
    1   2   3   4
    A4

See :class:`~sage.combinat.root_system.dynkin_diagram.DynkinDiagram`
for how to further manipulate Dynkin diagrams.

From this data (the *Cartan datum*), one can construct the associated
root system::

    sage: RootSystem(T)
    Root system of type ['A', 4]

The associated Weyl group is the symmetric group `S_{n+1}`::

    sage: W = WeylGroup(T)
    sage: W
    Weyl Group of type ['A', 4] (as a matrix group acting on the ambient space)
    sage: W.cardinality()
    120

while the Lie algebra is `sl_{n+1}`, and the Lie group `SL_{n+1}`
(TODO: illustrate this once this is implemented)

One may also construct crystals associated to various Dynkin diagrams. For example::

    sage: C = CrystalOfLetters(T)
    sage: C
    The crystal of letters for type ['A', 4]
    sage: C.list()
    [1, 2, 3, 4, 5]

    sage: C = CrystalOfTableaux(T, shape=[2])
    sage: C
    The crystal of tableaux of type ['A', 4] and shape(s) [[2]]
    sage: C.cardinality()
    15

Here is a sample of all the finite irreducible crystalographic Cartan
types::

    sage: CartanType.samples(finite = True, crystalographic = True)
    [['A', 1], ['A', 5], ['B', 1], ['B', 5], ['C', 1], ['C', 5], ['D', 2], ['D', 3], ['D', 5],
     ['E', 6], ['E', 7], ['E', 8], ['F', 4], ['G', 2]]

Non-crystallographic Cartan types are also partially supported::

    sage: CartanType.samples(finite = True, crystalographic = False)
    [['I', 5], ['H', 3], ['H', 4]]

In Sage, a Cartan type is used as a database of type-specific
information and algorithms (see e.g. :mod:`sage.combinat.root_system.type_A`).
This database includes how to construct the Dynkin diagram, the ambient space for the root
system (see http://en.wikipedia.org/wiki/Root_system), and further
mathematical properties::

    sage: T.is_finite(), T.is_simply_laced(), T.is_affine(), T.is_crystalographic()
    (True, True, False, True)

It will eventually include Coxeter numbers, etc.

In particular, a Sage Cartan type is endowed with a fixed choice of
labels for the nodes of the Dynkin diagram. This choice follows the
conventions of Nicolas Bourbaki, Lie Groups and Lie Algebras: Chapter 4-6, Elements of Mathematics,
Springer (2002). ISBN 978-3540426509. For example::

    sage: T = CartanType(['D', 4])
    sage: DynkinDiagram(T)
        O 4
        |
        |
    O---O---O
    1   2   3
    D4

    sage: E6 = CartanType(['E',6])
    sage: DynkinDiagram(E6)
            O 2
            |
            |
    O---O---O---O---O
    1   3   4   5   6
    E6


If desired, other node labelling conventions can be achieved. For
example the Kac labelling for type 'E_6` can be obtained via::

    sage: E6.relabel({1:1,2:6,3:2,4:3,5:4,6:5}).dynkin_diagram()
            O 6
            |
            |
    O---O---O---O---O
    1   2   3   4   5
    E6 relabelled by {1: 1, 2: 6, 3: 2, 4: 3, 5: 4, 6: 5}

Contributions implementing other conventions are very welcome.

Another option is to build from scratch a new Dynkin diagram.  The
architecture has been designed to make it fairly easy to add other
labelling conventions. In particular, we strived at choosing type free
algorithms whenever possible, so in principle most feature should
remain available even with custom Cartan types. This has not been used
much yet, so some rough corners certainly remain.

Here, we construct the hyperbolic example of Exercise 4.9 p. 57 of
Kac, Infinite Dimensional Lie Algebras. We start with an empty Dynkin
diagram, and add a couple nodes::

    sage: g = DynkinDiagram()
    sage: g.add_vertices([1,2,3])

Note that the diagonal of the Cartan matrix is already initialized::

    sage: g.cartan_matrix()
    [2 0 0]
    [0 2 0]
    [0 0 2]

Then we add a couple edges::

    sage: g.add_edge(1,2,2)
    sage: g.add_edge(1,3)
    sage: g.add_edge(2,3)

and we get the desired Cartan matrix::

    sage: g.cartan_matrix()
    [ 2 -1 -1]
    [-2  2 -1]
    [-1 -1  2]

Note that backward edges have been automatically added::

    sage: g.edges()
    [(1, 2, 2), (1, 3, 1), (2, 1, 1), (2, 3, 1), (3, 1, 1), (3, 2, 1)]

Caveat: the Dynkin diagram should not be modified after having been
used; this is not checked currently (Todo: add a method
:meth:`.set_mutable`, as for matrices).

.. rubric:: Reducible Cartan types

Reducible Cartan types can be specified by passing a sequence
or list of irreducible Cartan types::

    sage: CartanType(['A',2],['B',2])
    A2xB2
    sage: CartanType([['A',2],['B',2]])
    A2xB2
    sage: CartanType(['A',2],['B',2]).is_reducible()
    True

or using the following short hand notation::

    sage: CartanType("A2xB2")
    A2xB2
    sage: CartanType("A2","B2") == CartanType("A2xB2")
    True

.. rubric:: Degenerate cases

When possible, type `I_n` is automatically converted to the isomorphic
crystalographic Cartan types (any reason not to do so?)::

    sage: CartanType(["I",1])
    A1xA1
    sage: CartanType(["I",3])
    ['A', 2]
    sage: CartanType(["I",4])
    ['C', 2]
    sage: CartanType(["I",6])
    ['G', 2]

The Dynkin diagrams for types `B_1`, `C_1`, `D_2`, and `D_3` are
isomorphic to that for `A_1`, `A_1`, `A_1 \times A_1`, and `A_3`,
respectively. However their natural ambient space realizations (stemming
from the corresponding infinite families of Lie groups) are different.
Therefore, the Cartan types are considered as distinct::

    sage: CartanType(['B',1])
    ['B', 1]
    sage: CartanType(['C',1])
    ['C', 1]
    sage: CartanType(['D',2])
    ['D', 2]
    sage: CartanType(['D',3])
    ['D', 3]

.. rubric:: Affine Cartan types

For affine types, we use the usual conventions for affine Coxeter groups: each affine type
is either untwisted (that is arise from the natural affinisation
of a finite cartan type)::

    sage: CartanType(["A", 4, 1]).dynkin_diagram()
    0
    O-----------+
    |           |
    |           |
    O---O---O---O
    1   2   3   4
    A4~
    sage: CartanType(["B", 4, 1]).dynkin_diagram()
        O 0
        |
        |
    O---O---O=>=O
    1   2   3   4
    B4~

or dual thereof::

    sage: CartanType(["B", 4, 1]).dual().dynkin_diagram()
        O 0
        |
        |
    O---O---O=<=O
    1   2   3   4
    B4~*

or is of type `\widetilde{BC}_n` (which yields an irreducible, but nonreduced root system)::

    sage: CartanType(["BC", 4, 2]).dynkin_diagram()
    O=<=O---O---O=<=O
    0   1   2   3   4
    BC4~

This includes the two degenerate cases::

    sage: CartanType(["A", 1, 1]).dynkin_diagram()
    O<=>O
    0   1
    A1~
    sage: CartanType(["BC", 1, 2]).dynkin_diagram()
      4
    O=<=O
    0   1
    BC1~

For the user convenience, Kac's notations for twisted affine types are
automatically translated into the previous ones::

    sage: CartanType(["A", 9, 2])
    ['B', 5, 1]^*
    sage: CartanType(["A", 9, 2]).dynkin_diagram()
        O 0
        |
        |
    O---O---O---O=<=O
    1   2   3   4   5
    B5~*
    sage: CartanType(["A", 10, 2]).dynkin_diagram()
    O=<=O---O---O---O=<=O
    0   1   2   3   4   5
    BC5~
    sage: CartanType(["D", 5, 2]).dynkin_diagram()
    O=<=O---O---O=>=O
    0   1   2   3   4
    C4~*
    sage: CartanType(["D", 4, 3]).dynkin_diagram()
      3
    O=>=O---O
    2   1   0
    G2~* relabelled by {0: 0, 1: 2, 2: 1}
    sage: CartanType(["E", 6, 2]).dynkin_diagram()
    O---O---O=<=O---O
    0   1   2   3   4
    F4~*
"""
#*****************************************************************************
#       Copyright (C) 2007 Mike Hansen <mhansen@gmail.com>,
#       Copyright (C) 2008-2009 Nicolas M. Thiery <nthiery at users.sf.net>,
#
#  Distributed under the terms of the GNU General Public License (GPL)
#                  http://www.gnu.org/licenses/
#*****************************************************************************
from sage.misc.cachefunc import cached_method
from sage.misc.abstract_method import abstract_method
from sage.rings.all import ZZ
from sage.structure.sage_object import SageObject
from sage.structure.unique_representation import UniqueRepresentation
from sage.sets.family import Family

# TODO:
# Implement type_dual.AmbientSpace
# Implement dual ambient space
# Implement affinisation of ambient space
# Implement the Kac conventions by relabeling/dual/... of the above
# Implement coxeter diagrams for non crystalographic


# Intention: we want simultaneously CartanType to be a factory for
# the various subtypes of CartanType_abstract, as in:
#     CartanType(["A",4,1])
# and to behaves as a "module" for some extra utilities:
#     CartanType.samples()
#
# Implementation: CartanType is the unique instance of this class
# CartanTypeFactory. Is there a better/more standard way to do it?

class CartanTypeFactory(SageObject):
    def __call__(self, *args):
        """
        Constructs a Cartan type object.

        INPUT:
         - ``[letter, rank]``: letter is one of 'A','B','C','D','E','F','G' and rank is an integer
         - ``[letter, rank, twist]``: letter is one of 'A','B','C','D','E','F','G', 'BC', and rank and twist are integers
         - ``str``: a string
         - ``object``: a cartan type, or an object with a cartan type method

        EXAMPLES:

        We construct the Cartan type `D_4`::

            sage: d4 = CartanType(['D',4])
            sage: d4
            ['D', 4]

        or, for short::

            sage: CartanType("D4")
            ['D', 4]

        SEE ALSO: func:`~sage.combinat.root_system.cartan_type.CartanType`

        """
        if len(args) == 1:
            t = args[0]
        else:
            t = args
        if isinstance(t, CartanType_abstract):
            return t
        if hasattr(t, "cartan_type"):
            return  t.cartan_type()

        if type(t)==str:
            if "x" in t:
                import type_reducible
                return type_reducible.CartanType([CartanType(u) for u in t.split("x")])
            else:
                return CartanType([t[0], eval(t[1:])])

        t = list(t)

        if type(t[0]) == str and t[1] in ZZ and t[1] >= 0:
            letter, n = t[0], t[1]
            if len(t) == 2:
                if letter == "A":
                    if n >= 0:
                        import type_A
                        return type_A.CartanType(n)
                if letter == "B":
                    if n >= 1:
                        import type_B
                        return type_B.CartanType(n)
                if letter == "C":
                    if n >= 1:
                        import type_C
                        return type_C.CartanType(n)
                if letter == "D":
                    import type_D
                    if n >= 2:
                        return type_D.CartanType(n)
                if letter == "E":
                    if n >= 6 and n <= 8:
                        import type_E
                        return type_E.CartanType(n)
                if letter == "F":
                    if n == 4:
                        import type_F
                        return type_F.CartanType()
                if letter == "G":
                    if n == 2:
                        import type_G
                        return type_G.CartanType()
                if letter == "H":
                    if n in [3, 4]:
                        import type_H
                        return type_H.CartanType(n)
                if letter == "I":
                    if n == 1:
                        return CartanType([["A", 1], ["A", 1]])
                    if n == 3:
                        return CartanType(["A", 2])
                    if n == 4:
                        return CartanType(["C", 2])
                    if n == 6:
                        return CartanType(["G", 2])
                    if n >= 1:
                        import type_I
                        return type_I.CartanType(n)
            if len(t) == 3:
                if t[2] == 1: # Untwisted affind
                    if letter == "A":
                        if n >= 1:
                            import type_A_affine
                            return type_A_affine.CartanType(n)
                    if letter == "B":
                        if n >= 1:
                            import type_B_affine
                            return type_B_affine.CartanType(n)
                    if letter == "C":
                        if n >= 1:
                            import type_C_affine
                            return type_C_affine.CartanType(n)
                    if letter == "D":
                        import type_D_affine
                        if n >= 3:
                            return type_D_affine.CartanType(n)
                    if letter == "E":
                        if n >= 6 and n <= 8:
                            import type_E_affine
                            return type_E_affine.CartanType(n)
                    if letter == "F":
                        if n == 4:
                            import type_F_affine
                            return type_F_affine.CartanType()
                    if letter == "G":
                        if n == 2:
                            import type_G_affine
                            return type_G_affine.CartanType()
                if t[2] in [2,3]:
                    if letter == "BC" and t[2] == 2:
                        if n >= 1:
                            import type_BC_affine
                            return type_BC_affine.CartanType(n)
                    if letter == "A" and t[2] == 2:
                        if n%2 == 0: # Kac' A_2n^(2)
                            return CartanType(["BC", ZZ(n/2), 2])
                        else:        # Kac' A_2n-1^(2)
                            return CartanType(["B", ZZ((n+1)/2), 1]).dual()
                    if letter == "D" and t[2] == 2:
                        return CartanType(["C", n-1, 1]).dual()
                    if letter == "D" and t[2] == 3 and n == 4:
                        return CartanType(["G", 2, 1]).dual().relabel([0,2,1])
                    if letter == "E" and t[2] == 2 and n == 6:
                        return CartanType(["F", 4, 1]).dual()
            raise ValueError, "%s is not a valid cartan type"%t
        import type_reducible
        return type_reducible.CartanType([ CartanType(subtype) for subtype in t ])

    def _repr_(self):
        """
            sage: CartanType
            CartanType
        """
        return "CartanType"

    def samples(self, finite=None, affine=None, crystalographic=None):
        """
        Returns a sample of the available Cartan types.

        INPUT:
         - ``finite``: a boolean or None; defaults to None.
         - ``affine``: a boolean or None; defaults to None.
         - ``crystalographic``: a boolean or None; defaults to None.

        The sample contains all the exceptional finite and affine
        Cartan types, as well as typical representatives of the
        infinite families.

        EXAMPLES::

            sage: CartanType.samples()
            [['A', 1], ['A', 5], ['B', 1], ['B', 5], ['C', 1], ['C', 5], ['D', 2], ['D', 3], ['D', 5],
             ['E', 6], ['E', 7], ['E', 8], ['F', 4], ['G', 2], ['I', 5], ['H', 3], ['H', 4],
             ['A', 1, 1], ['A', 5, 1], ['B', 1, 1], ['B', 5, 1],
             ['C', 1, 1], ['C', 5, 1], ['D', 3, 1], ['D', 5, 1],
             ['E', 6, 1], ['E', 7, 1], ['E', 8, 1], ['F', 4, 1], ['G', 2, 1],
             ['B', 5, 1]^*, ['C', 4, 1]^*, ['F', 4, 1]^*, ['G', 2, 1]^*,
             ['BC', 1, 2], ['BC', 5, 2]]

        The finite, affine and crystalographic options allow
        respectively for restricting to (non) finite, (non) affine,
        and (non) crystalographic Cartan types::

            sage: CartanType.samples(finite=True)
            [['A', 1], ['A', 5], ['B', 1], ['B', 5], ['C', 1], ['C', 5], ['D', 2], ['D', 3], ['D', 5],
             ['E', 6], ['E', 7], ['E', 8], ['F', 4], ['G', 2], ['I', 5], ['H', 3], ['H', 4]]

            sage: CartanType.samples(affine=True)
            [['A', 1, 1], ['A', 5, 1], ['B', 1, 1], ['B', 5, 1],
             ['C', 1, 1], ['C', 5, 1], ['D', 3, 1], ['D', 5, 1],
             ['E', 6, 1], ['E', 7, 1], ['E', 8, 1], ['F', 4, 1], ['G', 2, 1],
             ['B', 5, 1]^*, ['C', 4, 1]^*, ['F', 4, 1]^*, ['G', 2, 1]^*,
             ['BC', 1, 2], ['BC', 5, 2]]

            sage: CartanType.samples(crystalographic=True)
            [['A', 1], ['A', 5], ['B', 1], ['B', 5], ['C', 1], ['C', 5], ['D', 2], ['D', 3], ['D', 5],
             ['E', 6], ['E', 7], ['E', 8], ['F', 4], ['G', 2],
             ['A', 1, 1], ['A', 5, 1], ['B', 1, 1], ['B', 5, 1],
             ['C', 1, 1], ['C', 5, 1], ['D', 3, 1], ['D', 5, 1],
             ['E', 6, 1], ['E', 7, 1], ['E', 8, 1], ['F', 4, 1], ['G', 2, 1],
             ['B', 5, 1]^*, ['C', 4, 1]^*, ['F', 4, 1]^*, ['G', 2, 1]^*,
             ['BC', 1, 2], ['BC', 5, 2]]

            sage: CartanType.samples(crystalographic=False)
            [['I', 5], ['H', 3], ['H', 4]]

        Todo: add some reducible Cartan types (suggestions?)

        """
        result = self._samples()
        if crystalographic is not None:
            result = [t for t in result if t.is_crystalographic() == crystalographic ]
        if finite is not None:
            result = [t for t in result if t.is_finite() == finite]
        if affine is not None:
            result = [t for t in result if t.is_affine() == affine]
        return result

    @cached_method
    def _samples(self):
        finite_crystalographic = \
            [CartanType (t)       for t in [['A', 1], ['A', 5], ['B', 1], ['B', 5],
                                            ['C', 1], ['C', 5], ['D', 2], ['D', 3], ['D', 5],
                                            ["E", 6], ["E", 7], ["E", 8],
                                            ["F", 4],
                                            ["G", 2]]]

        # Support for hand constructed Dynkin diagrams as Cartan types is not yet ready enough for including an example here.
        # from sage.combinat.root_system.dynkin_diagram import DynkinDiagram_class
        # g = DynkinDiagram_class.an_instance()
        return finite_crystalographic + \
            [CartanType(t)        for t in [["I", 5], ["H", 3], ["H", 4]]] + \
            [t.affine()           for t in finite_crystalographic if t.is_irreducible()] + \
            [CartanType(t).dual() for t in [["B", 5, 1], ["C", 4, 1], ["F", 4, 1], ["G", 2, 1]]] + \
            [CartanType(t)        for t in [["BC", 1, 2], ["BC", 5, 2], ]] #+ \
            # [ g ]


CartanType = CartanTypeFactory()

CartanType.__doc__ = __doc__

class CartanType_abstract(object):
    r"""
    Abstract class for Cartan types

    Subclasses should implement::

    - dynkin_diagram()

    - cartan_matrix()

    - is_finite()

    - is_affine()

    - is_irreducible()
    """

    def type(self):
        r"""
        Returns the type of self, or None if unknown. This method should be
        overridden in any subclass.

        EXAMPLES::

            sage: from sage.combinat.root_system.cartan_type import CartanType_abstract
            sage: C = CartanType_abstract()
            sage: C.type() is None
            True
        """
        return None

    def _add_abstract_superclass(self, cls):
        from sage.structure.dynamic_class import dynamic_class
        self.__class__ = dynamic_class(self.__class__.__name__+"_with_superclass", (self.__class__, cls))

    @abstract_method
    def rank(self):
        """
        Returns the rank of self, that is the number of nodes of the Dynkin diagram

        EXAMPLES::

            sage: CartanType(['A', 4]).rank()
            4
            sage: CartanType(['A', 7, 2]).rank()
            5
            sage: CartanType(['I', 8]).rank()
            2
        """
        #return len(self.index_set())

    @abstract_method
    def index_set(self):
        """
        Returns the index set for self. This is the list of the nodes
        of the associated Coxeter / Dynkin diagram.

        EXAMPLES::

            sage: CartanType(['A', 3, 1]).index_set()
            [0, 1, 2, 3]
            sage: CartanType(['D', 4]).index_set()
            [1, 2, 3, 4]
            sage: CartanType(['A', 7, 2]).index_set()
            [0, 1, 2, 3, 4]
            sage: CartanType(['A', 7, 2]).index_set()
            [0, 1, 2, 3, 4]
            sage: CartanType(['A', 6, 2]).index_set()
            [0, 1, 2, 3]
            sage: CartanType(['D', 6, 2]).index_set()
            [0, 1, 2, 3, 4, 5]
            sage: CartanType(['E', 6, 1]).index_set()
            [0, 1, 2, 3, 4, 5, 6]
            sage: CartanType(['E', 6, 2]).index_set()
            [0, 1, 2, 3, 4]
            sage: CartanType(['A', 2, 2]).index_set()
            [0, 1]
            sage: CartanType(['G', 2, 1]).index_set()
            [0, 1, 2]
            sage: CartanType(['F', 4, 1]).index_set()
            [0, 1, 2, 3, 4]
        """

    @abstract_method(optional = True)
    def coxeter_diagram(self):
        """
        Returns the Coxeter diagram for self.
        """

    def dual(self):
        """
        Returns the dual cartan type, possibly just as a formal dual.

        EXAMPLES::

            sage: CartanType(['A',3]).dual()
            ['A', 3]
            sage: CartanType(["B", 3]).dual()
            ['C', 3]
            sage: CartanType(['C',2]).dual()
            ['B', 2]
            sage: CartanType(['D',4]).dual()
            ['D', 4]
            sage: CartanType(['E',8]).dual()
            ['E', 8]
            sage: CartanType(['F',4]).dual()
            ['F', 4]^*
        """
        import type_dual
        return type_dual.CartanType(self)

    def relabel(self, relabelling):
        """
        INPUT:

         - ``type`` -- a Cartan type
         - ``relabelling`` -- a function (or a list, or a dictionary)

        Returns an isomorphic Cartan type obtained by relabelling the
        nodes of the dynkin diagram. Namely the node with label ``i``
        is relabelled ``f(i)`` (or, by ``f[i]`` if f is a list or
        dictionary).

        EXAMPLES::

           sage: CartanType(['F',4]).relabel({ 1:4, 2:3, 3:2, 4:1 }).dynkin_diagram()
           O---O=>=O---O
           4   3   2   1
           F4 relabelled by {1: 4, 2: 3, 3: 2, 4: 1}

        """
        import type_relabel
        return type_relabel.CartanType(self, relabelling)

    def is_reducible(self):
        """
        Report whether the root system is reducible (i.e. not simple), that
        is whether it can be factored as a product of root systems.

        EXAMPLES::

            sage: CartanType("A2xB3").is_reducible()
            True
            sage: CartanType(['A',2]).is_reducible()
            False
        """
        return not self.is_irreducible()

    def is_irreducible(self):
        """
        Report whether this Cartan type is irreducible (i.e. simple). This
        should be overridden in any subclass.

        This returns False by default. Derived class should override this appropriately.

        EXAMPLES::

            sage: from sage.combinat.root_system.cartan_type import CartanType_abstract
            sage: C = CartanType_abstract()
            sage: C.is_irreducible()
            False

        """
        return False

    def is_atomic(self):
        """
        This method is usually equivalent to :meth:'.is_reducible',
        except for the Cartan type `D_2`.

        `D_2` is not a standard Cartan type. It is equivalent to type
        `A_1\times A_1` which is reducible; however the isomorphism
        from its ambient space (for the orthogonal group of degree 4)
        to that of `A_1\times A_1` is non trivial, and it is useful to
        have it.

        From a programming point of view its implementation is more
        similar to the irreducible types, and so the method
        is_atomic() is supplied.

        EXAMPLES::

            sage: CartanType("D2").is_atomic()
            True
            sage: CartanType("D2").is_irreducible()
            False

        TESTS::

            sage: all( T.is_irreducible() == T.is_atomic() for T in CartanType.samples() if T != CartanType("D2"))
            True
        """
        return self.is_irreducible()

    def is_compound(self):
        """
        A short hand for not :meth:`.is_atomic`.

        TESTS::

            sage: all( T.is_compound() == (not T.is_atomic()) for T in CartanType.samples())
            True
        """
        return not self.is_atomic()

    @abstract_method
    def is_finite(self):
        """
        Returns whether this Cartan type is finite.

        EXAMPLES::

            sage: from sage.combinat.root_system.cartan_type import CartanType_abstract
            sage: C = CartanType_abstract()
            sage: C.is_finite()
            Traceback (most recent call last):
            ...
            NotImplementedError: <abstract method is_finite at ...>

        ::

            sage: CartanType(['A',4]).is_finite()
            True
            sage: CartanType(['A',4, 1]).is_finite()
            False
        """

    @abstract_method
    def is_affine(self):
        """
        Returns whether self is affine.

        EXAMPLES::

            sage: CartanType(['A', 3]).is_affine()
            False
            sage: CartanType(['A', 3, 1]).is_affine()
            True
        """

    def is_crystalographic(self):
        """
        Returns whether this Cartan type is crystalographic.

        This returns False by default. Derived class should override this appropriately.

        EXAMPLES::

            sage: [ [t, t.is_crystalographic() ] for t in CartanType.samples(finite=True) ]
            [[['A', 1], True], [['A', 5], True],
             [['B', 1], True], [['B', 5], True],
             [['C', 1], True], [['C', 5], True],
             [['D', 2], True], [['D', 3], True], [['D', 5], True],
             [['E', 6], True], [['E', 7], True], [['E', 8], True],
             [['F', 4], True], [['G', 2], True],
             [['I', 5], False], [['H', 3], False], [['H', 4], False]]

        TESTS::

            sage: all(t.is_crystalographic() for t in CartanType.samples(affine=True))
            True
        """
        return False

    def is_simply_laced(self):
        """
        Returns whether this Cartan type is simply laced

        This returns False by default. Derived class should override this appropriately.

        EXAMPLES::

            sage: [ [t, t.is_simply_laced() ] for t in CartanType.samples() ]
            [[['A', 1], True], [['A', 5], True],
             [['B', 1], True], [['B', 5], False],
             [['C', 1], True], [['C', 5], False],
             [['D', 2], True], [['D', 3], True], [['D', 5], True],
             [['E', 6], True], [['E', 7], True], [['E', 8], True],
             [['F', 4], False], [['G', 2], False], [['I', 5], False], [['H', 3], False], [['H', 4], False],
             [['A', 1, 1], False], [['A', 5, 1], True],
             [['B', 1, 1], False], [['B', 5, 1], False],
             [['C', 1, 1], False], [['C', 5, 1], False],
             [['D', 3, 1], True], [['D', 5, 1], True],
             [['E', 6, 1], True], [['E', 7, 1], True], [['E', 8, 1], True],
             [['F', 4, 1], False], [['G', 2, 1], False],
             [['B', 5, 1]^*, False], [['C', 4, 1]^*, False], [['F', 4, 1]^*, False], [['G', 2, 1]^*, False],
             [['BC', 1, 2], False], [['BC', 5, 2], False]]
        """
        return False

    def is_implemented(self):
        try:
            self.dynkin_diagram()
            return True
        except StandardError:
            return False

    def root_system(self):
        """
        Returns the root system associated to self.

        EXAMPLES::

            sage: CartanType(['A',4]).root_system()
            Root system of type ['A', 4]
        """
        from sage.combinat.root_system.root_system import RootSystem
        return RootSystem(self)

class CartanType_crystalographic(CartanType_abstract):
    """
    An abstract class for crystalographic cartan types
    """
    @abstract_method
    def dynkin_diagram(self):
        """
        Returns the Dynkin diagram associated with self.

        EXAMPLES::

            sage: CartanType(['A',4]).dynkin_diagram()
            O---O---O---O
            1   2   3   4
            A4

        Note: derived subclasses should typically implement this as a cached method
        """

    def cartan_matrix(self):
        """
        Returns the Cartan matrix associated with self.

        EXAMPLES::

            sage: CartanType(['A',4]).cartan_matrix()
            [ 2 -1  0  0]
            [-1  2 -1  0]
            [ 0 -1  2 -1]
            [ 0  0 -1  2]
        """
        from sage.combinat.root_system.cartan_matrix import cartan_matrix
        return cartan_matrix(self)

    def is_crystalographic(self):
        """
        Implements :meth:`CartanType_abstract.is_crystalographic`
        by returning `True`.

        EXAMPLES::

            sage: CartanType(['A', 3, 1]).is_crystalographic()
            True
        """
        return True

    def symmetrizer(self):
        """
        A Cartan matrix `M` is symmetrizable if there exists a non
        trivial diagonal matrix `D` such that `DM` is a symmetric
        matrix, that is `DM = M^tD`. In that case, `D` is unique, up
        to a scalar factor for each connected component of the Dynkin
        diagram.

        This method currently assumes that the Cartan type is irreducible.

        This method computes the unique minimal such `D` with non
        negative integral coefficients. If `D` exists, it is returned
        as a family. Otherwise None is returned.

        EXAMPLES::

            sage: CartanType(["B",5]).symmetrizer()
            Finite family {1: 2, 2: 2, 3: 2, 4: 2, 5: 1}

        Here is a neat trick to visualize it better::

            sage: T = CartanType(["B",5])
            sage: print T.ascii_art(T.symmetrizer().__getitem__)
            O---O---O---O=>=O
            2   2   2   2   1

            sage: T = CartanType(["BC",5, 2])
            sage: print T.ascii_art(T.symmetrizer().__getitem__)
            O=<=O---O---O---O=<=O
            1   2   2   2   2   4

        Property: up to an overall scalar factor, this gives the norm
        of the simple roots in the ambient space::

            sage: T = CartanType(["C",5])
            sage: print T.ascii_art(T.symmetrizer().__getitem__)
            O---O---O---O=<=O
            1   1   1   1   2

            sage: alpha = RootSystem(T).ambient_space().simple_roots()
            sage: print T.ascii_art(lambda i: alpha[i].scalar(alpha[i]))
            O---O---O---O=<=O
            2   2   2   2   4

        """
        assert self.is_irreducible()
        from sage.matrix.constructor import matrix, diagonal_matrix
        m = self.cartan_matrix()
        n = m.nrows()
        M = matrix(ZZ, n, n*n, sparse = True)
        for (i,j) in m.nonzero_positions():
            M[i, n * i + j]  = m[i,j]
            M[j, n * i + j] -= m[j,i]
        kern = M.integer_kernel()
        assert kern.dimension() <= 1
        if kern.dimension() == 0:
            return None
        D = kern.basis()[0]
        assert diagonal_matrix(D) * m == m.transpose() * diagonal_matrix(D)
        I = self.index_set()
        return Family( dict( (I[i], D[i]) for i in range(n) ) )


class CartanType_simply_laced(CartanType_crystalographic):
    """
    An abstract class for simply laced cartan types
    """
    def is_simply_laced(self):
        """
        EXAMPLES::

            sage: CartanType(['A',3,1]).is_simply_laced()
            True
            sage: CartanType(['A',2]).is_simply_laced()
            True
        """
        return True

    def dual(self):
        """
        Simply laced cartan types are self dual

        EXAMPLES::

            sage: CartanType(["A", 3]).dual()
            ['A', 3]
            sage: CartanType(["A", 3, 1]).dual()
            ['A', 3, 1]
            sage: CartanType(["D", 3]).dual()
            ['D', 3]
            sage: CartanType(["D", 4, 1]).dual()
            ['D', 4, 1]
            sage: CartanType(["E", 6]).dual()
            ['E', 6]
            sage: CartanType(["E", 6, 1]).dual()
            ['E', 6, 1]
        """
        return self

class CartanType_simple(CartanType_abstract):
    """
    An abstract class for simple Cartan types
    """
    def is_irreducible(self):
        """
        EXAMPLES::

            sage: CartanType(['A', 3]).is_irreducible()
            True
        """
        return True

class CartanType_finite(CartanType_abstract):
    """
    An abstract class for simple affine Cartan types
    """
    def is_finite(self):
        """
        EXAMPLES::

            sage: CartanType(["A", 3]).is_finite()
            True
        """
        return True

    def is_affine(self):
        """
        EXAMPLES::

            sage: CartanType(["A", 3]).is_affine()
            False
        """
        return False

class CartanType_affine(CartanType_simple, CartanType_crystalographic):
    """
    An abstract class for simple affine Cartan types
    """
    def is_finite(self):
        """
        EXAMPLES::

            sage: CartanType(['A', 3, 1]).is_finite()
            False
        """
        return False

    def is_affine(self):
        """
        EXAMPLES::

            sage: CartanType(['A', 3, 1]).is_affine()
            True
        """
        return True

    def is_untwisted_affine(self):
        """
        Returns whether self is untwisted affine, that is if it is the
        affine extension of some finite type.

        Every affine type is either untwisted affine, or dual thereof
        or of type ``BC``.

        EXAMPLES::

            sage: CartanType(['A', 3, 1]).is_untwisted_affine()
            True
            sage: CartanType(['A', 3, 1]).dual().is_untwisted_affine() # this one is self dual!
            True
            sage: CartanType(['B', 3, 1]).dual().is_untwisted_affine()
            False
            sage: CartanType(['BC', 3, 2]).is_untwisted_affine()
            False
        """
        return False

    @abstract_method
    def special_node(self):
        r"""

        A special node is a node of the Dynkin diagram such that
        pruning it yields a Dynkin diagram for the associated
        classical type (see :meth:`.classical`).

        This method returns the label of some special node. This is
        usually `0` in the standard conventions.

        EXAMPLES::

            sage: CartanType(['A', 3, 1]).special_node()
            0

        The choice is guaranteed to be consistent with the indexing of
        the nodes of the classical Dynkin diagram::

            sage: CartanType(['A', 3, 1]).index_set()
            [0, 1, 2, 3]
            sage: CartanType(['A', 3, 1]).classical().index_set()
            [1, 2, 3]

        """

    @abstract_method
    def classical(self):
        r"""
        Returns the classical Cartan type associated with this affine Cartan type.

        EXAMPLES::

            sage: CartanType(['A', 1, 1]).classical()
            ['A', 1]
            sage: CartanType(['A', 3, 1]).classical()
            ['A', 3]
            sage: CartanType(['B', 3, 1]).classical()
            ['B', 3]

            sage: CartanType(['A', 2, 2]).classical()
            ['C', 1]
            sage: CartanType(['BC', 1, 2]).classical()
            ['C', 1]
            sage: CartanType(['A', 4, 2]).classical()
            ['C', 2]
            sage: CartanType(['BC', 2, 2]).classical()
            ['C', 2]
            sage: CartanType(['A', 10, 2]).classical()
            ['C', 5]
            sage: CartanType(['BC', 5, 2]).classical()
            ['C', 5]

            sage: CartanType(['D', 5, 2]).classical()
            ['B', 4]
            sage: CartanType(['E', 6, 1]).classical()
            ['E', 6]
            sage: CartanType(['G', 2, 1]).classical()
            ['G', 2]
            sage: CartanType(['E', 6, 2]).classical() # todo: double check
            ['F', 4]^*
            sage: CartanType(['D', 4, 3]).classical() # todo: double check
            ['G', 2]^* relabelled by {1: 2, 2: 1}

        We check that :meth:`classical`,
        :meth:`sage.combinat.root_system.cartan_type.CartanType_crystalographic.dynkin_diagram`,
        and :meth:`.special_node` are consistent::

            sage: for ct in CartanType.samples(affine = True):
            ...       g1 = ct.classical().dynkin_diagram()
            ...       g2 = ct.dynkin_diagram()
            ...       g2.delete_vertex(ct.special_node())
            ...       assert sorted(g1.vertices()) == sorted(g2.vertices())
            ...       assert sorted(g1.edges()) == sorted(g2.edges())

        """

    def row_annihilator(self, m = None):
        r"""
        Return the unique minimal non trivial annihilating linear
        combination of `\alpha_0, \alpha1, \ldots, \alpha_n` with
        nonnegative coefficients (or alternatively, the unique minimal
        non trivial annihilating linear combination of the rows of the
        Cartan matrix with non-negative coefficients).

        Throw an error if the existence of uniqueness does not hold

        The optional argument ``m`` is for internal use only.

        EXAMPLES::

            sage: RootSystem(['C',2,1]).cartan_type().acheck()
            Finite family {0: 1, 1: 1, 2: 1}
            sage: RootSystem(['D',4,1]).cartan_type().acheck()
            Finite family {0: 1, 1: 1, 2: 2, 3: 1, 4: 1}
            sage: RootSystem(['F',4,1]).cartan_type().acheck()
            Finite family {0: 1, 1: 2, 2: 3, 3: 2, 4: 1}
            sage: RootSystem(['BC',4,2]).cartan_type().acheck()
            Finite family {0: 1, 1: 2, 2: 2, 3: 2, 4: 2}

        FIXME:
         - The current implementation assumes that the Cartan matrix
           is indexed by `[0,1,...]`, in the same order as the index set.
         - This really should be a method of CartanMatrix
        """
        if m is None:
            m = self.cartan_matrix()
        assert self.index_set() == list(range(m.ncols()))
        annihilator_basis = m.integer_kernel().gens()
        assert(len(annihilator_basis) == 1)
        assert(all(coef > 0 for coef in annihilator_basis[0]))

        return Family(dict((i,annihilator_basis[0][i])for i in self.index_set()))

    acheck = row_annihilator

    def col_annihilator(self):
        r"""
        Return the unique minimal non trivial annihilating linear
        combination of `\alpha^\vee_0, \alpha^\vee, \ldots, \alpha^\vee` with
        nonnegative coefficients (or alternatively, the unique minimal
        non trivial annihilating linear combination of the columns of the
        Cartan matrix with non-negative coefficients).

        Throw an error if the existence or uniqueness does not hold

        FIXME: the current implementation assumes that the Cartan
        matrix is indexed by `[0,1,...]`, in the same order as the
        index set.

        EXAMPLES::

            sage: RootSystem(['C',2,1]).cartan_type().a()
            Finite family {0: 1, 1: 2, 2: 1}
            sage: RootSystem(['D',4,1]).cartan_type().a()
            Finite family {0: 1, 1: 1, 2: 2, 3: 1, 4: 1}
            sage: RootSystem(['F',4,1]).cartan_type().a()
            Finite family {0: 1, 1: 2, 2: 3, 3: 4, 4: 2}
            sage: RootSystem(['BC',4,2]).cartan_type().a()
            Finite family {0: 2, 1: 2, 2: 2, 3: 2, 4: 1}
        """
        return self.row_annihilator(self.cartan_matrix().transpose())

    a = col_annihilator

    def c(self):
        r"""
        Returns the family (c_i)_i of integer coefficients defined by
        `c_i=max(1, a_i/a^vee_i)` (see e.g. [FSS07]_ p. 3)

        FIXME: the current implementation assumes that the Cartan
        matrix is indexed by `[0,1,...]`, in the same order as the
        index set.

        EXAMPLES::

            sage: RootSystem(['C',2,1]).cartan_type().c()
            Finite family {0: 1, 1: 2, 2: 1}
            sage: RootSystem(['D',4,1]).cartan_type().c()
            Finite family {0: 1, 1: 1, 2: 1, 3: 1, 4: 1}
            sage: RootSystem(['F',4,1]).cartan_type().c()
            Finite family {0: 1, 1: 1, 2: 1, 3: 2, 4: 2}
            sage: RootSystem(['BC',4,2]).cartan_type().c()
            Finite family {0: 2, 1: 1, 2: 1, 3: 1, 4: 1}

        TESTS::

            sage: CartanType(["B", 3, 1]).c().map(parent)
            Finite family {0: Integer Ring, 1: Integer Ring, 2: Integer Ring, 3: Integer Ring}

        REFERENCES:

        .. [FSS07] G. Fourier, A. Schilling, and M. Shimozono,
           Demazure structure inside Kirillov-Reshetikhin crystals,
           J. Algebra, Vol. 309, (2007), p. 386-404
           http://arxiv.org/abs/math/0605451
        """
        a = self.a()
        acheck = self.acheck()
        return Family(dict((i, max(ZZ(1), a[i] // acheck[i]))
                           for i in self.index_set()))

    def translation_factors(self):
        r"""
        Returns the translation factors for ``self``. Those are the
        smallest factors `t_i` such that the translation by `t_i
        \alpha_i` maps the fundamental polygon to another polygon in
        the alcove picture.

        OUTPUT: a dictionary from ``self.index_set()`` to `\ZZ`
        (or `\QQ` for affine type `BC`)

        Those coefficients are all `1` for dual untwisted, and in
        particular for simply laced. They coincide with the usual
        `c_i` coefficients (see :meth:`c`) for untwisted and dual
        thereof. See the discussion below for affine type `BC`.

        Note: one usually realizes the alcove picture in the coweight
        lattice, with translations by coroots; in that case, one will
        use the translation factors for the dual Cartan type.

        FIXME: the current implementation assumes that the Cartan
        matrix is indexed by `[0,1,...]`, in the same order as the
        index set.

        EXAMPLES::

            sage: CartanType(['C',2,1]).translation_factors()
            Finite family {0: 1, 1: 2, 2: 1}
            sage: CartanType(['C',2,1]).dual().translation_factors()
            Finite family {0: 1, 1: 1, 2: 1}
            sage: CartanType(['D',4,1]).translation_factors()
            Finite family {0: 1, 1: 1, 2: 1, 3: 1, 4: 1}
            sage: CartanType(['F',4,1]).translation_factors()
            Finite family {0: 1, 1: 1, 2: 1, 3: 2, 4: 2}
            sage: CartanType(['BC',4,2]).translation_factors()
            Finite family {0: 1, 1: 1, 2: 1, 3: 1, 4: 1/2}

        We proceed with systematic tests taken from MuPAD-Combinat's
        testsuite::

            sage: list(CartanType(["A", 1, 1]).translation_factors())
            [1, 1]
            sage: list(CartanType(["A", 5, 1]).translation_factors())
            [1, 1, 1, 1, 1, 1]
            sage: list(CartanType(["B", 5, 1]).translation_factors())
            [1, 1, 1, 1, 1, 2]
            sage: list(CartanType(["C", 5, 1]).translation_factors())
            [1, 2, 2, 2, 2, 1]
            sage: list(CartanType(["D", 5, 1]).translation_factors())
            [1, 1, 1, 1, 1, 1]
            sage: list(CartanType(["E", 6, 1]).translation_factors())
            [1, 1, 1, 1, 1, 1, 1]
            sage: list(CartanType(["E", 7, 1]).translation_factors())
            [1, 1, 1, 1, 1, 1, 1, 1]
            sage: list(CartanType(["E", 8, 1]).translation_factors())
            [1, 1, 1, 1, 1, 1, 1, 1, 1]
            sage: list(CartanType(["F", 4, 1]).translation_factors())
            [1, 1, 1, 2, 2]
            sage: list(CartanType(["G", 2, 1]).translation_factors())
            [1, 3, 1]
            sage: list(CartanType(["A", 2, 2]).translation_factors())
            [1, 1/2]
            sage: list(CartanType(["A", 2, 2]).dual().translation_factors())
            [1/2, 1]
            sage: list(CartanType(["A", 10, 2]).translation_factors())
            [1, 1, 1, 1, 1, 1/2]
            sage: list(CartanType(["A", 10, 2]).dual().translation_factors())
            [1/2, 1, 1, 1, 1, 1]
            sage: list(CartanType(["A", 9, 2]).translation_factors())
            [1, 1, 1, 1, 1, 1]
            sage: list(CartanType(["D", 5, 2]).translation_factors())
            [1, 1, 1, 1, 1]
            sage: list(CartanType(["D", 4, 3]).translation_factors())
            [1, 1, 1]
            sage: list(CartanType(["E", 6, 2]).translation_factors())
            [1, 1, 1, 1, 1]

        We conclude with a discussion of the appropriate value for
        affine type `BC`. Let us consider the alcove picture realized
        in the weight lattice. It is obtained by taking the level-`1`
        affine hyperplane in the weight lattice, and projecting it
        along `\Lambda_0`::

            sage: R = RootSystem(["BC",2,2])
            sage: alpha = R.weight_space().simple_roots()
            sage: alphacheck = R.coroot_space().simple_roots()
            sage: Lambda = R.weight_space().fundamental_weights()

        Here are the levels of the fundamental weights::

            sage: Lambda[0].level(), Lambda[1].level(), Lambda[2].level()
            (1, 2, 2)

        So the "center" of the fundamental polygon at level `1` is::

            sage: O = Lambda[0]
            sage: O.level()
            1

        We take the projection `\omega_1` at level `0` of `\Lambda_1`
        as unit vector on the `x`-axis, and the projection `\omega_2`
        at level 0 of `\Lambda_2` as unit vector of the `y`-axis::

            sage: omega1 = Lambda[1]-2*Lambda[0]
            sage: omega2 = Lambda[2]-2*Lambda[0]
            sage: omega1.level(), omega2.level()
            (0, 0)

        The projections of the simple roots can be read off::

            sage: alpha[0]
            2*Lambda[0] - Lambda[1]
            sage: alpha[1]
            -2*Lambda[0] + 2*Lambda[1] - Lambda[2]
            sage: alpha[2]
            -2*Lambda[1] + 2*Lambda[2]

        Namely `\alpha_0 = -\omega_1`, `\alpha_1 = 2\omega_1 -
        \omega_2` and `\alpha_2 = -2 \omega_1 + 2 \omega_2`.

        The reflection hyperplane defined by `\alpha_0^\vee` goes
        through the points `O+1/2 \omega_1` and `O+1/2 \omega_2`::

            sage: (O+(1/2)*omega1).scalar(alphacheck[0])
            0
            sage: (O+(1/2)*omega2).scalar(alphacheck[0])
            0

        Hence, the fundamental alcove is the triangle `(O, O+1/2
        \omega_1, O+1/2 \omega_2)`. By successive reflections, one can
        tile the full plane. This induces a tiling of the full plane
        by translates of the fundamental polygon.

        TODO: add the picture here, once root system plots in the
        weight lattice will be implemented. In the mean time, the
        reader may look up the dual picture on Figure 2 of [HST09]_
        which was produced by MuPAD-Combinat.

        From this picture, one can read that translations by
        `\alpha_0`, `\alpha_1`, and `1/2\alpha_2` map the fundamental
        polygon to translates of it in the alcove picture, and are
        smallest with this property. Hence, the translation factors
        for affine type `BC` are `t_0=1, t_1=1, t_2=1/2`::

            sage: CartanType(['BC',2,2]).translation_factors()
            Finite family {0: 1, 1: 1, 2: 1/2}

        TESTS::

            sage: CartanType(["B", 3, 1]).translation_factors().map(parent)
            Finite family {0: Integer Ring, 1: Integer Ring, 2: Integer Ring, 3: Integer Ring}
            sage: CartanType(["BC", 3, 2]).translation_factors().map(parent)
            Finite family {0: Integer Ring, 1: Integer Ring, 2: Integer Ring, 3: Rational Field}

        REFERENCES:

        .. [HST09] F. Hivert, A. Schilling, and N. M. Thiery,
           Hecke group algebras as quotients of affine Hecke
           algebras at level 0, JCT A, Vol. 116, (2009) p. 844-863
           http://arxiv.org/abs/0804.3781
        """
        a = self.a()
        acheck = self.acheck()
        if set([1/ZZ(2), 2]).issubset( set(a[i]/acheck[i] for i in self.index_set()) ):
            # The test above and the formula below are rather meaningless
            # But they detect properly type BC or dual and return the correct value
            return Family(dict((i, min(ZZ(1), a[i] / acheck[i]))
                               for i in self.index_set()))

        else:
            return self.c()

##############################################################################
# Concrete base classes

class CartanType_standard_finite(UniqueRepresentation, SageObject, CartanType_finite):
    """
    A concrete base class for the finite standard Cartan types `A_3`, `D_4`, `E_8`
    """
    def __init__(self, letter, n):
        """
        EXAMPLES::

            sage: ct = CartanType(['A',4])

        TESTS::

            sage: TestSuite(ct).run(verbose = True)
            running ._test_category() . . . pass
            running ._test_not_implemented_methods() . . . pass
            running ._test_pickling() . . . pass
        """
#         assert(t[0] in ['A', 'B', 'C', 'D', 'E', 'F', 'G', 'H', 'I'])
#         assert(t[1] in ZZ and t[1] >= 0)
#         if t[0] in ['B', 'C']:
#             assert(t[1] >= 2)
#         if t[0] == 'D':
#             assert(t[1] >= 3)
#         if t[0] == 'E':
#             assert(t[1] <= 8)
#         if t[0] == 'F':
#             assert(t[1] <= 4)
#         if t[0] == 'G':
#             assert(t[1] <= 2)
#         if t[0] == 'H':
#             assert(t[1] <= 4)
        self.letter = letter
        self.n = n

    # Technical methods
    def _repr_(self, compact = False):
        """
        TESTS::

            sage: ct = CartanType(['A',3])
            sage: repr(ct)
            "['A', 3]"
        """
        format = '%s%s' if compact else "['%s', %s]"
        return format%(self.letter, self.n)

    def __reduce__(self):
        """
        TESTS::

            sage: T = CartanType(['D', 4])
            sage: T.__reduce__()
            (CartanType, ('D', 4))
            sage: T == loads(dumps(T))
            True

        """
        from cartan_type import CartanType
        return (CartanType, (self.letter, self.n))

    def __cmp__(self, other):
         """
         TESTS::

             sage: ct1 = CartanType(['A',4])
             sage: ct2 = CartanType(['A',4])
             sage: ct3 = CartanType(['A',5])
             sage: ct1 == ct2
             True
             sage: ct1 != ct3
             True
         """
         if other.__class__ != self.__class__:
             return cmp(self.__class__, other.__class__)
         if other.letter != self.letter:
             return cmp(self.letter, other.letter)
         return cmp(self.n, other.n)

    def __hash__(self):
        """
        EXAMPLES::

            sage: ct = CartanType(['A',2])
            sage: hash(ct) #random
            -5684143898951441983
        """
        return hash((self.n, self.letter))

    def __getitem__(self, i):
        """
        EXAMPLES::

            sage: t = CartanType(['A', 3, 1])
            sage: t[0]
            'A'
            sage: t[1]
            3
            sage: t[2]
            1
            sage: t[3]
            Traceback (most recent call last):
            ...
            IndexError: index out of range
        """
        if i == 0:
            return self.letter
        elif i==1:
            return self.n
        else:
            raise IndexError("index out of range")

    def __len__(self):
        """
        EXAMPLES::

            sage: len(CartanType(['A',4]))
            2
        """
        return 2

    # mathematical methods

    def index_set(self):
        """
        Implements :meth:`CartanType_abstract.index_set`.

        The index set for all standard finite Cartan types is of the form `{1,\dots,n}`.
        (but see type_I for a slight abuse of this).

        EXAMPLES::
            sage: CartanType(['A', 5]).index_set()
            [1, 2, 3, 4, 5]

        """
        return range(1,self.n+1)

    def rank(self):
        """
        EXAMPLES::

            sage: CartanType(['A', 3]).rank()
            3
            sage: CartanType(['B', 3]).rank()
            3
            sage: CartanType(['C', 3]).rank()
            3
            sage: CartanType(['D', 4]).rank()
            4
            sage: CartanType(['E', 6]).rank()
            6
        """
        return self.n

    def affine(self):
        """
        Returns the corresponding untwisted affine Cartan type

        EXAMPLES::

            sage: CartanType(['A',3]).affine()
            ['A', 3, 1]
        """
        return CartanType([self.letter, self.n, 1])

    def type(self):
        """
        Returns the type of self.

        EXAMPLES::

            sage: CartanType(['A', 4]).type()
            'A'
            sage: CartanType(['A', 4, 1]).type()
            'A'
        """
        return self.letter


##########################################################################
class CartanType_standard_affine(UniqueRepresentation, SageObject, CartanType_affine):
    r"""
    A concrete class for affine simple Cartan types

    """

    def __init__(self, letter, n, affine = 1):
        """
        EXAMPLES::

            sage: ct = CartanType(['A',4,1])
            sage: ct == loads(dumps(ct))
            True

        TESTS::

            sage: ct1 = CartanType(['A',3, 1])
            sage: ct2 = CartanType(['B',3, 1])
            sage: ct3 = CartanType(['A',3])
            sage: ct1 == ct1
            True
            sage: ct1 == ct2
            False
            sage: ct1 == ct3
            False

        """
        assert(letter in ['A', 'B', 'C', 'BC', 'D', 'E', 'F', 'G'])
        self.letter = letter
        self.n = n
        self.affine = affine

    def _repr_(self, compact = False):
        """
        TESTS::

            sage: ct = CartanType(['A',3, 1])
            sage: repr(ct)
            "['A', 3, 1]"
        """
        if compact:
            return '%s%s~'%(self.letter, self.n)
        else:
            return "['%s', %s, %s]"%(self.letter, self.n, self.affine)

    def __reduce__(self):
        """
        TESTS::

            sage: T = CartanType(['D', 4, 1])
            sage: T.__reduce__()
            (CartanType, ('D', 4, 1))
            sage: T == loads(dumps(T))
            True

        """
        from sage.combinat.root_system.cartan_type import CartanType
        return (CartanType, (self.letter, self.n, self.affine))

    def __len__(self):
        """
        EXAMPLES::

            sage: len(CartanType(['A',4,1]))
            3
        """
        return 3

    def __getitem__(self, i):
        """
        EXAMPLES::

            sage: t = CartanType(['A', 3, 1])
            sage: t[0]
            'A'
            sage: t[1]
            3
            sage: t[2]
            1
            sage: t[3]
            Traceback (most recent call last):
            ...
            IndexError: index out of range
        """
        if i == 0:
            return self.letter
        elif i==1:
            return self.n
        elif i == 2:
            return self.affine
        else:
            raise IndexError("index out of range")

    def rank(self):
        """
        EXAMPLES::

            sage: CartanType(['A', 4, 1]).rank()
            5
            sage: CartanType(['B', 4, 1]).rank()
            5
            sage: CartanType(['C', 3, 1]).rank()
            4
            sage: CartanType(['D', 4, 1]).rank()
            5
            sage: CartanType(['E', 6, 1]).rank()
            7
            sage: CartanType(['E', 7, 1]).rank()
            8
            sage: CartanType(['F', 4, 1]).rank()
            5
            sage: CartanType(['G', 2, 1]).rank()
            3
            sage: CartanType(['A', 2, 2]).rank()
            2
            sage: CartanType(['A', 6, 2]).rank()
            4
            sage: CartanType(['A', 7, 2]).rank()
            5
            sage: CartanType(['D', 5, 2]).rank()
            5
            sage: CartanType(['E', 6, 2]).rank()
            5
            sage: CartanType(['D', 4, 3]).rank()
            3
        """
        return self.n+1

    def index_set(self):
        """
        Implements :meth:`CartanType_abstract.index_set`.

        The index set for all standard affine Cartan types is of the form `{0,\dots,n}`.

        EXAMPLES::
            sage: CartanType(['A', 5, 1]).index_set()
            [0, 1, 2, 3, 4, 5]

        """
        return range(0,self.n+1)

    def special_node(self):
        r"""
        Implements :meth:'CartanType_abstract.special_node`.

        With the standard labelling conventions, `0` is always a
        special node.

        EXAMPLES::

            sage: CartanType(['A', 3, 1]).special_node()
            0
        """
        return 0

    def type(self):
        """
        Returns the type of self.

        EXAMPLES::

            sage: CartanType(['A', 4, 1]).type()
            'A'
        """
        return self.letter

##########################################################################
class CartanType_standard_untwisted_affine(CartanType_standard_affine):
    r"""
    A concrete class for the standard untwisted affine Cartan types
    """
    def classical(self):
        r"""
        Returns the classical Cartan type associated with self (which
        should be affine)

        EXAMPLES::

            sage: CartanType(['A', 3, 1]).classical()
            ['A', 3]
            sage: CartanType(['B', 3, 1]).classical()
            ['B', 3]
            sage: CartanType(['C', 3, 1]).classical()
            ['C', 3]
            sage: CartanType(['D', 4, 1]).classical()
            ['D', 4]
            sage: CartanType(['E', 6, 1]).classical()
            ['E', 6]
            sage: CartanType(['F', 4, 1]).classical()
            ['F', 4]
            sage: CartanType(['G', 2, 1]).classical()
            ['G', 2]
        """
        return CartanType([self.letter,self.n])

    def is_untwisted_affine(self):
        """
        Implements :meth:'CartanType_affine.is_untwisted_affine` by
        returning True.

        EXAMPLES::

            sage: CartanType(['B', 3, 1]).is_untwisted_affine()
            True

        """
        return True

##############################################################################
# For backward compatibility
class CartanType_simple_finite(object):
    def __setstate__(self, dict):
        """
        Implements the unpickling of Cartan types pickled by Sage <= 4.0.

        EXAMPLES:

        This is the pickle for CartanType(["A", 4])::

            sage: pg_CartanType_simple_finite = unpickle_global('sage.combinat.root_system.cartan_type', 'CartanType_simple_finite')
            sage: si1 = unpickle_newobj(pg_CartanType_simple_finite, ())
            sage: pg_unpickleModule = unpickle_global('twisted.persisted.styles', 'unpickleModule')
            sage: pg_make_integer = unpickle_global('sage.rings.integer', 'make_integer')
            sage: si2 = pg_make_integer('4')
            sage: unpickle_build(si1, {'tools':pg_unpickleModule('sage.combinat.root_system.type_A'), 't':['A', si2], 'letter':'A', 'n':si2})

            sage: si1
            ['A', 4]
            sage: si1.dynkin_diagram()
            O---O---O---O
            1   2   3   4
            A4

        This is quite hacky; in particular unique representation is not preserved::

            sage: si1 == CartanType(["A", 4]) # todo: not implemented
            True
        """
        T = CartanType([dict['letter'], dict['n']])
        self.__class__ = T.__class__
        self.__dict__ = T.__dict__

