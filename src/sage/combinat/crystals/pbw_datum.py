r"""
PBW Data

This contains helper classes and functions which encode PBW data
in finite type.

AUTHORS:

- Dinakar Muthiah (2015-05): initial version
"""

#*****************************************************************************
#       Copyright (C) 2015 Dinakar Muthiah <your email>
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 2 of the License, or
# (at your option) any later version.
#                  http://www.gnu.org/licenses/
#*****************************************************************************

#from sage.misc.lazy_attribute import lazy_attribute
from sage.misc.cachefunc import cached_method
from sage.combinat.root_system.cartan_type import CartanType
from sage.combinat.root_system.coxeter_group import CoxeterGroup
from sage.combinat.root_system.root_system import RootSystem
from sage.combinat.crystals.braid_move_calculator import BraidMoveCalculator

class PBWDatum(object):
    def __init__(self, parent, long_word, lusztig_datum):
        """
        Initialize ``self``.
        """
        self.parent = parent
        self.long_word = tuple(long_word)
        self.lusztig_datum = tuple(lusztig_datum)

    def __repr__(self):
        """
        Return a string representation of ``self``.
        """
        return_str = "PBW Datum element of type {cartan_type} with ".format(
                     cartan_type=self.parent.cartan_type)
        return_str += "long word {long_word} and Lusztig datum {lusztig_datum}".format(
                      long_word=self.long_word,
                      lusztig_datum=self.lusztig_datum) 
        return return_str

    def __eq__(self, other_PBWDatum):
        """
        Check equality.

        EXAMPLES::

            sage: from sage.combinat.crystals.pbw_datum import PBWData, PBWDatum
            sage: P = PBWData("A2")
            sage: L1 = PBWDatum(P, (1,2,1), (1,4,7))
            sage: L2 = PBWDatum(P, (1,2,1), (1,4,7))
            sage: L1 == L2
            True
        """
        return (self.parent == other_PBWDatum.parent and
                self.long_word == other_PBWDatum.long_word and
                self.lusztig_datum == other_PBWDatum.lusztig_datum)

    def is_equivalent_to(self, other_pbw_datum):
        r"""
        Return whether ``self`` is equivalent to ``other_pbw_datum``.
        modulo the tropical Plucker relations.

        EXAMPLES::

            sage: from sage.combinat.crystals.pbw_datum import PBWData, PBWDatum
            sage: P = PBWData("A2")
            sage: L1 = PBWDatum(P, (1,2,1), (1,0,1))
            sage: L2 = PBWDatum(P, (2,1,2), (0,1,0))
            sage: L1.is_equivalent_to(L2)
            True
            sage: L1 == L2
            False
        """
        other_long_word = other_pbw_datum.long_word
        other_lusztig_datum = other_pbw_datum.lusztig_datum
        equiv_pbw_datum = self.convert_to_new_long_word(other_long_word)
        return equiv_pbw_datum.lusztig_datum == other_lusztig_datum
        
    def convert_to_long_word_with_first_letter(self, i):
        r"""
        Return a new PBWDatum equivalent to ``self``
        whose long word begins with ``i``.
        """
        return self.convert_to_new_long_word(self.parent._long_word_begin_with(i))

    def convert_to_new_long_word(self, new_long_word):
        r"""
        Return a new PBWDatum equivalent to ``self``
        whose long word is ``new_long_word``.

        EXAMPLES::

            sage: from sage.combinat.crystals.pbw_datum import PBWData, PBWDatum
            sage: P = PBWData("A2")
            sage: datum = PBWDatum(P, (1,2,1), (1,0,1))
            sage: new_datum = P.convert_to_new_long_word(datum, (2,1,2))
            sage: new_datum.long_word
            (2, 1, 2)
            sage: new_datum.lusztig_datum
            (0, 1, 0)
        """
        return self.parent.convert_to_new_long_word(self, new_long_word)

    def weight(self):
        """
        Return the weight of ``self``.

        EXAMPLES::

            sage: from sage.combinat.crystals.pbw_datum import PBWData, PBWDatum
            sage: P = PBWData("A2")
            sage: L = PBWDatum(P, (1,2,1), (1,1,1))
            sage: L.weight()
            -2*alpha[1] - 2*alpha[2]
        """
        root_list = self.parent._root_list_from(tuple(self.long_word))
        R = self.parent.root_lattice
        return R.linear_combination((root_list[i], -coeff)
                                    for i, coeff in enumerate(self.lusztig_datum))

    def star(self):
        """
        Return the starred version of ``self``, i.e.,
        with reversed `long_word` and `lusztig_datum`

        EXAMPLES::

            sage: from sage.combinat.crystals.pbw_datum import PBWData, PBWDatum
            sage: P = PBWData("A2")
            sage: L1 = PBWDatum(P, (1,2,1), (1,2,3))
            sage: L1.star() == PBWDatum(P, (1,2,1), (3,2,1))
            True
        """
        reversed_long_word = reversed(self.long_word)
        reversed_lusztig_datum = reversed(self.lusztig_datum)
        return PBWDatum(self.parent, reversed_long_word, reversed_lusztig_datum)


class PBWData(object): # UniqueRepresentation?
    """
    Helper class for the set of PBW data.
    """
    def __init__(self, cartan_type):
        """
        Initialize ``self``.
        """
        self.cartan_type = CartanType(cartan_type)
        self.root_system = RootSystem(self.cartan_type)
        self.root_lattice = self.root_system.root_lattice()
        self.weyl_group = self.root_lattice.weyl_group()
        self._braid_move_calc = BraidMoveCalculator(self.weyl_group)

    def convert_to_new_long_word(self, pbw_datum, new_long_word):
        """
        Convert the PBW datum ``pbw_datum`` from its long word to
        ``new_long_word``.

        EXAMPLES::

            sage: from sage.combinat.crystals.pbw_datum import PBWData, PBWDatum
            sage: P = PBWData("A2")
            sage: datum = PBWDatum(P, (1,2,1), (1,0,1))
            sage: new_datum = P.convert_to_new_long_word(datum,(2,1,2))
            sage: new_datum
            PBW Datum element of type ['A', 2] with long word (2, 1, 2)
             and Lusztig datum (0, 1, 0)
            sage: new_datum.long_word
            (2, 1, 2)
            sage: new_datum.lusztig_datum
            (0, 1, 0)
        """
        assert pbw_datum.parent is self
        chain = self._braid_move_calc.chain_of_reduced_words(pbw_datum.long_word,
                                                             new_long_word)
        enhanced_braid_chain = enhance_braid_move_chain(chain, self.cartan_type)
        new_lusztig_datum = compute_new_lusztig_datum(enhanced_braid_chain,
                                                      pbw_datum.lusztig_datum)
        return PBWDatum(self, new_long_word, new_lusztig_datum)

    @cached_method
    def _root_list_from(self, reduced_word):
        """
        Return the list of positive roots in the order determined by 
        ``reduced_word``.

        .. WARNING::

            No error checking is done to verify that ``reduced_word``
            is reduced.

        INPUT:

        - ``reduced_word`` -- a tuple corresponding to a reduced word

        EXAMPLES::

            sage: from sage.combinat.crystals.pbw_datum import PBWData
            sage: P = PBWData(["A",2])
            sage: P._root_list_from((1,2,1))
            [alpha[1], alpha[1] + alpha[2], alpha[2]]
        """
        al = self.root_lattice.simple_roots()
        cur = []
        for i in reversed(reduced_word):
            cur = [al[i]] + [x.simple_reflection(i) for x in cur]
        return cur

    @cached_method
    def _long_word_begin_with(self, i):
        """
        Return a reduced expression of the long word which begins with ``i``.
        """
        si = self.weyl_group.simple_reflection(i)
        w0 = self.weyl_group.long_element()
        return tuple([i] + (si * w0).reduced_word())

#enhanced_braid_chain is an ugly data structure.
def compute_new_lusztig_datum(enhanced_braid_chain, initial_lusztig_datum):
    """
    Return the lusztig datum obtained by applying Tropical Plucker
    relations along ``enhanced_braid_chain`` starting with
    ``initial_lusztig_datum``.

    EXAMPLES::

        sage: from sage.combinat.crystals.braid_move_calculator import BraidMoveCalculator
        sage: from sage.combinat.crystals.pbw_datum import enhance_braid_move_chain
        sage: from sage.combinat.crystals.pbw_datum import compute_new_lusztig_datum
        sage: ct = CartanType(['A', 2])
        sage: W = CoxeterGroup(ct)
        sage: B = BraidMoveCalculator(W)
        sage: chain = B.chain_of_reduced_words((1,2,1),(2,1,2))
        sage: enhanced_braid_chain = enhance_braid_move_chain(chain, ct)
        sage: compute_new_lusztig_datum(enhanced_braid_chain,(1,0,1))    
        (0, 1, 0)

    TESTS::

        sage: from sage.combinat.crystals.braid_move_calculator import BraidMoveCalculator
        sage: from sage.combinat.crystals.pbw_datum import enhance_braid_move_chain
        sage: from sage.combinat.crystals.pbw_datum import compute_new_lusztig_datum
        sage: ct = CartanType(['A', 2])
        sage: W = CoxeterGroup(ct)
        sage: B = BraidMoveCalculator(W)
        sage: chain = B.chain_of_reduced_words((1,2,1), (2,1,2))
        sage: enhanced_braid_chain = enhance_braid_move_chain(chain, ct)
        sage: compute_new_lusztig_datum(enhanced_braid_chain,(1,0,1)) == (0,1,0)
        True
    """
    # Does not currently check that len(initial_lusztig_datum) is appropriate
    new_lusztig_datum = list(initial_lusztig_datum) #shallow copy
    for interval_of_change, type_data in enhanced_braid_chain[1:]:
        a,b = interval_of_change
        old_interval_datum = new_lusztig_datum[a:b]
        new_interval_datum = tropical_plucker_relation(type_data, old_interval_datum)
        new_lusztig_datum[a:b] = new_interval_datum
    return tuple(new_lusztig_datum)

# The tropical plucker relations
def tropical_plucker_relation(a, lusztig_datum):
    r"""
    Apply the tropical Plucker relation of type ``a`` to ``lusztig_datum``.

    INPUT:

    - ``a`` -- a pair ``(x, y)`` of the off-diagonal entries of a
      `2 \times 2` Cartan matrix

    EXAMPLES::

        sage: from sage.combinat.crystals.pbw_datum import tropical_plucker_relation
        sage: tropical_plucker_relation((0,0), (2,3))
        (3, 2)
        # Add more doctests
    """
    n = lusztig_datum
    if a == (0, 0): # A1xA1
        return (n[1], n[0])
    elif a == (-1, -1): # A2
        p = min(n[0], n[2])
        return (n[1]+n[2]-p, p, n[0]+n[1]-p)
    elif a == (-1, -2): # B2
        p1 = min(n[0]+n[1], n[0]+n[3], n[2]+n[3])
        p2 = min(2*n[0]+n[1], 2*n[0]+n[3], 2*n[2]+n[3])
        return (n[1]+2*n[2]+n[3]-p2,
                p2-p1,
                2*p1-p2,
                n[0]+n[1]+n[2]-p1)
    elif a == (-2, -1): # C2
        # I believe this is condition (iii) in Proposition 5.2 of Joel's thesis.
        # (I'm pretty sure this is correct).
        p1 = min(n[0]+n[1], n[0]+n[3], n[2]+n[3])
        p2 = min(n[0]+2*n[1], n[0]+2*n[3], n[2]+2*n[3])
        return (n[1]+n[2]+n[3]-p1,
                2*p1-p2,
                p2-p1,
                n[0]+2*n[1]+n[2]-p2)
    elif a == (-3, -1): # G2
        raise NotImplementedError("type G2 not implemented")

# Maybe we need to be more specific, and pass not the Cartan type, but the root lattice?
# TODO: Move to PBW_data?
def enhance_braid_move_chain(braid_move_chain, cartan_type):
    r"""
    Return a list of tuples that records the data of the long words in
    ``braid_move_chain`` plus the data of the intervals where the braid moves
    occur and the data of the off-diagonal entries of the `2 \times 2` Cartan
    submatrices of each braid move.

    INPUT:

    - ``braid_move_chain`` -- a chain of reduced words in the Weyl group
      of ``cartan_type``
    - ``cartan_type`` -- a finite Cartan type

    OUTPUT:

    A list of 2-tuples
    ``(interval_of_change, cartan_sub_matrix)`` where

    - ``interval_of_change`` is the (half-open) interval of indices where
      the braid move occurs; this is `None` for the first tuple
    - ``cartan_sub_matrix`` is the off-diagonal entries of the `2 \times 2`
      submatrix of the cartan matrix corresponding to the braid move;
      this is `None` for the first tuple

    For a matrix::

        [2 a]
        [b 2]

    the ``cartan_sub_matrix`` is the pair ``(a, b)``.

    TESTS::

        sage: from sage.combinat.crystals.pbw_datum import enhance_braid_move_chain
        sage: braid_chain = [(1, 2, 1, 3, 2, 1),
        ....:                (1, 2, 3, 1, 2, 1),
        ....:                (1, 2, 3, 2, 1, 2),
        ....:                (1, 3, 2, 3, 1, 2),
        ....:                (3, 1, 2, 3, 1, 2),
        ....:                (3, 1, 2, 1, 3, 2),
        ....:                (3, 2, 1, 2, 3, 2),
        ....:                (3, 2, 1, 3, 2, 3)]
        sage: enhanced_chain = enhance_braid_move_chain(braid_chain, CartanType(["A",5])) 
        sage: enhanced_chain[0]
        (None, None)
        sage: enhanced_chain[7]
        ((3, 6), (-1, -1))
    """
    cartan_matrix = cartan_type.cartan_matrix()
    output_list = []
    output_list.append( (None, None) )
    previous_word = braid_move_chain[0]
    # TODO - Optimize this by avoiding calls to diff_interval
    # This likely could be done when performing chain_of_reduced_words
    # Things in here get called the most (about 50x more than enhance_braid_move_chain)
    for current_word in braid_move_chain[1:]:
        interval_of_change = diff_interval(previous_word, current_word)
        i = previous_word[interval_of_change[0]] - 1  # -1 for indexing
        j = current_word[interval_of_change[0]] - 1  # -1 for indexing
        cartan_sub_matrix = (cartan_matrix[i,j], cartan_matrix[j,i])
        output_list.append( (interval_of_change, cartan_sub_matrix) )
        previous_word = current_word
    return output_list

# TODO: Cythonize this if this is still necessary
def diff_interval(list1, list2):
    """
    Return the smallest contiguous half-open interval [a,b) 
    that contains the indices where ``list1`` and ``list2`` differ. 
    Return ``None`` if the lists don't differ.

    INPUT:

    - ``list1``, ``list2`` -- two lists of the same length

    .. NOTE::

        The input is not checked for speed.

    TESTS::

        sage: from sage.combinat.crystals.pbw_datum import diff_interval
        sage: diff_interval([1,2,3,4], [1,2,3,4])
        sage: diff_interval([1,2,3,4], [1,3,2,4])
        (1, 3)
        sage: diff_interval([1,2,4,4,5], [1,3,45,6,3])
        (1, 5)
    """
    L = [i for i,elt in enumerate(list1) if elt != list2[i]]
    if not L:
        return None
    return (L[0], L[-1]+1)

