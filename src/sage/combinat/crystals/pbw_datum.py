r"""
PBW Data

This contains helper classes and functions which encode PBW data
in finite type.

AUTHORS:

- Dinakar Muthiah (2015-05): initial version
"""

#from sage.misc.lazy_attribute import lazy_attribute
from sage.misc.cachefunc import cached_method
#from sage.structure.parent import Parent
#from sage.structure.element import Element
from sage.combinat.root_system.cartan_type import CartanType
from sage.combinat.root_system.coxeter_group import CoxeterGroup
from sage.combinat.root_system.root_system import RootSystem
from sage.combinat.crystals.braid_move_calculator import (BraidMoveCalculator,
                                                          enhance_braid_move_chain)

class PBWDatum(object):
    def __init__(self, parent, long_word, lusztig_datum):
        self.parent = parent
        self.long_word = tuple(long_word)
        self.lusztig_datum = tuple(lusztig_datum)

    def __repr__(self):
        return_str = "PBW Datum element of type {cartan_type} with ".format(
            cartan_type=self.parent.cartan_type)
        return_str += "long word {long_word} and Lusztig datum {lusztig_datum}".format(
            long_word=self.long_word,
            lusztig_datum=self.lusztig_datum) 
        return return_str

    def __eq__(self, other_PBWDatum):
        """
        EXAMPLES::

            sage: P = PBWData(["A2"])
            sage: L1 = P((1,2,1),(1,4,7))
            sage: L2 = P((1,2,1),(1,4,7))
            sage: L1 == L2
            True
        """
        return (self.parent == other_PBWDatum.parent and
                self.long_word == other_PBWDatum.long_word and
                self.lusztig_datum == other_PBWDatum.lusztig_datum)

    def is_equivalent_to(self, other_pbw_datum):
        r"""
        Return whether ``self`` is equivalent to ``other_PBWDatum``.
        modulo the tropical Plucker relations.

        EXAMPLES::

            sage: P = PBWData(["A2"])
            sage: L1 = P((1,2,1),(1,0,1))
            sage: L2 = P((2,1,2),(0,1,0))
            sage: L1.is_equivalent_to(L2)
            True
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

            sage: P = PBWData(["A2"])
            sage: P.convert_to_new_long_word(P((1,2,1),(1,0,1)),(2,1,2)).long_word
            (2, 1, 2)
            sage: P.convert_to_new_long_word(P((1,2,1),(1,0,1)),(2,1,2)).lusztig_datum
            (0, 1, 0)
        """
        return self.parent.convert_to_new_long_word(self, new_long_word)

    def wt(self):
        """
        EXAMPLES::

            sage: P = PBWData(["A",2])
            sage: L = P((1,2,1),(1,1,1))
            sage: L.wt()
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

            sage: P = PBWData(["A2"])
            sage: L1 = P((1,2,1),(1,2,3))
            sage: L1.star() == P((1,2,1),(3,2,1))
            True
        """
        reversed_long_word = reversed(self.long_word)
        reversed_lusztig_datum = reversed(self.lusztig_datum)
        return PBWDatum(self.parent, reversed_long_word, reversed_lusztig_datum)


class PBWData(object): # UniqueRepresentation?
    def __init__(self, cartan_type):
        self.cartan_type = CartanType(cartan_type)
        self.root_system = RootSystem(self.cartan_type)
        self.root_lattice = self.root_system.root_lattice()
        self.weyl_group = self.root_lattice.weyl_group()
        # Is there a more intelligent way to recover the Weyl group
        # from cartan_type?
        self._braid_move_calc = BraidMoveCalculator(self.weyl_group)

    def convert_to_new_long_word(self, pbw_datum, new_long_word):
        assert pbw_datum.parent == self
        chain = self._braid_move_calc.chain_of_reduced_words(pbw_datum.long_word,
                                                             new_long_word)
        enhanced_braid_chain = enhance_braid_move_chain(chain,
                                                        self.cartan_type)
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
        s = self.weyl_group.simple_reflections()
        w0 = self.weyl_group.long_element()
        return tuple([i] + (s[i] * w0).reduced_word())

# Replaced by a doctest
def convert_to_new_long_word_test():
    P = PBWData(["A2"])
    P.convert_to_new_long_word(P((1,2,1),(1,0,1)),(2,1,2))
    assert P.convert_to_new_long_word(P((1,2,1),(1,0,1)),(2,1,2)).long_word == (2,1,2)
    assert P.convert_to_new_long_word(P((1,2,1),(1,0,1)),(2,1,2)).lusztig_datum == (0,1,0)


#enhanced_braid_chain is an ugly data structure.
def compute_new_lusztig_datum(enhanced_braid_chain, initial_lusztig_datum):
    """
    Return the lusztig datum obtained by applying Tropical Plucker relations along
    ``enhanced_braid_chain`` starting with ``initial_lusztig_datum``

    EXAMPLES::

        sage: W = CoxeterGroup(CartanType(["A2"]))
        sage: B = BraidMoveCalculator(W)
        sage: chain = B.chain_of_reduced_words((1,2,1),(2,1,2))
        sage: enhanced_braid_chain = enhance_braid_move_chain(chain,CartanType(["A",2]))
        sage: compute_new_lusztig_datum(enhanced_braid_chain,(1,0,1))    
        (0, 1, 0)
    """
    # Does not currently check that len(initial_lusztig_datum) is appropriate
    new_lusztig_datum = list(initial_lusztig_datum) #shallow copy
    for _, interval_of_change, type_data in enhanced_braid_chain[1:]:
        old_interval_datum = new_lusztig_datum.__getslice__(*interval_of_change)
        new_interval_datum = tropical_plucker_relation(type_data, old_interval_datum)
        new_lusztig_datum.__setslice__(interval_of_change[0], interval_of_change[1],
                                       new_interval_datum)
    return tuple(new_lusztig_datum)

# Replaced by a doctest
def compute_new_lusztig_datum_test():
    W = CoxeterGroup(CartanType(["A2"]))
    B = BraidMoveCalculator(W)
    chain = B.chain_of_reduced_words((1,2,1),(2,1,2))
    enhanced_braid_chain = enhance_braid_move_chain(chain,CartanType(["A",2]))
    #print compute_new_lusztig_datum(enhanced_braid_chain,(1,0,1))
    assert compute_new_lusztig_datum(enhanced_braid_chain,(1,0,1)) == (0,1,0)


# The tropical plucker relations
def tropical_plucker_relation(a, lusztig_datum):
    r"""
    Apply the tropical Plucker relation of type ``a`` to ``lusztig_datum``.

    INPUT:

    - ``a`` -- a pair ``(x, y)`` of the off-diagonal entries of a
      `2 \times 2` Cartan matrix

    EXAMPLES::

        sage: tropical_plucker_relation((0,0), (2,3))
        (3, 2)
    """
    n = lusztig_datum
    if a == (0, 0): # A1xA1
        return (n[1], n[0])
    elif a == (-1, -1): # A2
        p = min(n[0], n[2])
        return (n[1]+n[2]-p, p, n[0]+n[1]-p)
    elif a == (-1, -2): # B2
        # I believe this is condition (iii) in Proposition 5.2 of Joel's thesis.
        # (I'm pretty sure this is correct).
        p1 = min(n[0]+n[1], n[0]+n[3], n[2]+n[3])
        p2 = min(n[0]+2*n[1], n[0]+2*n[3], n[2]+2*n[3])
        return (n[1]+n[2]+n[3]-p1,
                2*p1-p2,
                p2-p1,
                n[0]+2*n[1]+n[2]-p2)
    elif a == (-2, -1): # C2
        p1 = min(n[0]+n[1], n[0]+n[3], n[2]+n[3])
        p2 = min(2*n[0]+n[1], 2*n[0]+n[3], 2*n[2]+n[3])
        return (n[1]+2*n[2]+n[3]-p2,
                p2-p1,
                2*p1-p2,
                n[0]+n[1]+n[2]-p1)
    elif a == (-3, -1): # G2
        raise NotImplementedError("type G2 not implemented")

