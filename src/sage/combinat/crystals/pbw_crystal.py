r""" 
The crystal of PBWData up to equivalence. A model of the crystal 
`B(\infty)` whose elements are PBWDatum up to equivalence by the 
tropical Plucker relations.

AUTHORS:

- Dinakar Muthiah (2015-05-11): initial version
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

from sage.misc.cachefunc import cached_method
from sage.structure.element import Element
from sage.structure.parent import Parent
from sage.structure.unique_representation import UniqueRepresentation
from sage.categories.highest_weight_crystals import HighestWeightCrystals
from sage.categories.infinite_enumerated_sets import InfiniteEnumeratedSets
from sage.combinat.root_system.cartan_type import CartanType
from sage.combinat.crystals.pbw_datum import PBWData, PBWDatum

class PBWCrystalElement(Element):
    """
    A crystal element in the PBW model.
    """
    def __init__(self, parent, long_word, lusztig_datum):
        """
        Initialize ``self``.
        """
        Element.__init__(self, parent)
        self._pbw_datum = PBWDatum(parent._pbw_datum_parent, long_word, lusztig_datum)

    def _repr_(self):
        """
        Return a string representation of ``self``.
        """
        pbw_datum = self._pbw_datum.convert_to_new_long_word(self.parent()._default_word)
        return "PBW monomial with Lusztig datum {}".format(pbw_datum.lusztig_datum)

    def _latex_(self):
        """
        Return a latex representation of ``self``.
        """
        pbw_datum = self._pbw_datum.convert_to_new_long_word(self.parent()._default_word)
        lusztig_datum = list(pbw_datum.lusztig_datum)
        al = self.parent()._pbw_datum_parent._root_list_from(self.parent()._default_word)
        from sage.misc.latex import latex
        ret_str = ' '.join("f_{%s}%s"%(latex(al[i]), "^{%s}"%latex(exp) if i > 1 else "")
                           for i, exp in enumerate(lusztig_datum) if exp)
        if ret_str == '':
            return '1'
        return ret_str

    def lusztig_datum(self, word=None):
        """
        Return the Lusztig datum of ``self`` with respect to the reduced
        expression of the long word ``word``.
        """
        if word is None:
            word = self.parent()._default_word
        else:
            self.parent()._check_is_long_word(word)
        pbw_datum = self._pbw_datum.convert_to_new_long_word(word)
        return tuple(pbw_datum.lusztig_datum)

    def __eq__(self, other):
        """
        Check equality of ``self`` with ``other``.
        """
        other_long_word = other._pbw_datum.long_word
        other_lusztig_datum = other._pbw_datum.lusztig_datum
        equiv_pbw_datum = self._pbw_datum.convert_to_new_long_word(other_long_word)
        return equiv_pbw_datum.lusztig_datum == other_lusztig_datum

    def _cmp_(self, other):
        """
        Return comparison of ``self`` and ``other``.
        """
        i = self.parent().index_set()[0]
        word = self.parent()._pbw_datum_parent._long_word_begin_with(i)
        lusztig_datum = tuple(self._pbw_datum.convert_to_new_long_word(word).lusztig_datum)
        other_lusztig_datum = tuple(other._pbw_datum.convert_to_new_long_word(word).lusztig_datum)
        return cmp(lusztig_datum, other_lusztig_datum)

    @cached_method
    def __hash__(self):
        """
        Return the hash of ``self``.
        """
        i = self.parent().index_set()[0]
        word = self.parent()._pbw_datum_parent._long_word_begin_with(i)
        pbw_datum = self._pbw_datum.convert_to_new_long_word(word)
        return hash(tuple(pbw_datum.lusztig_datum))

    def f(self, i):
        """
        Return the action of `f_i` on ``self``.

        EXAMPLES::

            sage: B = PBWCrystal("D4")
            sage: b = B.module_generators[0]
            sage: c = b.f_string([1,2,3,1,2,3,4])
            sage: c == b.f_string([1,2,4,1,2,3,3])
            True
        """
        equiv_PBWDatum = self._pbw_datum.convert_to_long_word_with_first_letter(i) 
        new_long_word = equiv_PBWDatum.long_word 
        new_lusztig_datum = list(equiv_PBWDatum.lusztig_datum)
        new_lusztig_datum[0] += 1
        return type(self)(self.parent(), new_long_word, tuple(new_lusztig_datum))

    def e(self, i):
        """
        Return the action of `e_i` on ``self``.
        """
        equiv_pbw_datum = self._pbw_datum.convert_to_long_word_with_first_letter(i) 
        new_long_word = equiv_pbw_datum.long_word 
        new_lusztig_datum = list(equiv_pbw_datum.lusztig_datum)
        if new_lusztig_datum[0] == 0:
            return None
        new_lusztig_datum[0] -= 1
        return type(self)(self.parent(), new_long_word, tuple(new_lusztig_datum))

    def epsilon(self, i):
        r"""
        Return `\varepsilon_i` of ``self``.

        EXAMPLES::

            sage: B = PBWCrystal(["A2"])
            sage: s = B((1,2,1),(3,0,0))
            sage: s.epsilon(1)
            3
        """
        equiv_pbw_datum = self._pbw_datum.convert_to_long_word_with_first_letter(i) 
        return equiv_pbw_datum.lusztig_datum[0]

    def phi(self, i):
        r"""
        Return `\varphi_i` of ``self``.

        EXAMPLES::

            sage: B = PBWCrystal(["A2"])
            sage: s = B((1,2,1),(3,0,0))
            sage: s.phi(1)
            -3
            sage: s.phi(2)
            3
        """
        WLR = self.parent().weight_lattice_realization()
        h = WLR.simple_coroots()
        return self.epsilon(i) + self.weight().scalar(h[i])

    def weight(self):
        """
        Return weight of ``self``.

        EXAMPLES::

            sage: B = PBWCrystal(["A",2])
            sage: s = B((1,2,1),(2,2,2)) 
            sage: s.weight()
            -4*alpha[1] - 4*alpha[2]
        """
        WLR = self.parent().weight_lattice_realization()
        al = WLR.simple_roots()
        return WLR.sum(c*al[i] for i,c in self._pbw_datum.wt())

    def star(self):
        r"""
        Return the starred crystal element corresponding
        to ``self``.

        EXAMPLES::

            sage: P = PBWCrystal(["A2"])
            sage: P((1,2,1),(1,2,3)).star() == P((1,2,1),(3,2,1))
            True
        """
        starred_pbw_datum = self._pbw_datum.star()
        return type(self)(self.parent(), starred_pbw_datum.long_word,
                             starred_pbw_datum.lusztig_datum)


class PBWCrystal(Parent, UniqueRepresentation):
    """
    Crystal of `B(\infty)` given by PBW monomials.
    """
    @staticmethod
    def __classcall_private__(cls, cartan_type):
        """
        Normalize input to ensure a unique representation.
        """
        return super(PBWCrystal, cls).__classcall__(cls, CartanType(cartan_type))

    def __init__(self, cartan_type):
        """
        Initialize ``self``.
        """
        self._cartan_type = cartan_type
        self._pbw_datum_parent = PBWData(self._cartan_type)
        category = (HighestWeightCrystals(), InfiniteEnumeratedSets())
        Parent.__init__(self, category=category)

        # There must be a better way to do the following
        i = self._cartan_type.index_set()[0]
        self._default_word = self._pbw_datum_parent._long_word_begin_with(i)
        zero_lusztig_datum = [0]*len(self._default_word)
        self.module_generators = (self.element_class(self, 
                                                     self._default_word,
                                                     zero_lusztig_datum),)

    def _repr_(self):
        """
        Return a string representation of ``self``.
        """
        return "Crystal of PBW data of type {}".format(self._cartan_type)

    def default_long_word(self):
        """
        Return the default long word used to express elements of ``self``.
        """
        return self._default_word

    def _check_is_long_word(self, word):
        """
        Check if ``word`` is a reduced expression of the long of the
        Coxeter group of ``self``.
        """
        W = self._pbw_datum_parent.weyl_group
        if (len(word) != len(self._default_word)
            or W.from_reduced_word(word) != W.long_element()):
            raise ValueError("not a reduced word of the long element")

    def set_default_long_word(self, word):
        """
        Set the default long word used to express elements of ``self``.
        """
        self._check_is_long_word(word)
        self._default_word = tuple(word)

    Element = PBWCrystalElement

