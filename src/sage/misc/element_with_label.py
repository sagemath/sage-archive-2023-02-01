# -*- coding: utf-8 -*-
r"""
Elements with labels.

This module implements a simple wrapper class for pairs consisting of
an "element" and a "label".
For representation purposes (``repr``, ``str``, ``latex``), this pair
behaves like its label, while the element is "silent".
However, these pairs compare like usual pairs (i.e., both element and
label have to be equal for two such pairs to be equal).
This is used for visual representations of graphs and posets.
"""

from sage.misc.latex import latex

class ElementWithLabel:
    """
    Auxiliary class for showing/viewing :class:`Poset`s with
    non-injective labelings.
    For hashing and equality testing the resulting object behaves
    like a tuple ``(element, label)``.
    For any presentation purposes it appears just as ``label`` would.
    """
    def __init__(self, element, label):
        """
        Construct an object that wraps ``element`` but presents itself
        as ``label``.

        TESTS::

            sage: from sage.misc.element_with_label import ElementWithLabel
            sage: e = ElementWithLabel(1, 'a')
            sage: e
            'a'
            sage: e.element
            1
        """
        self.element = element
        self.label = label

    def _latex_(self):
        """
        Return the latex representation of ``self``,
        which is just the latex representation of the label.

        TESTS::

            sage: var('a_1')
            a_1
            sage: from sage.misc.element_with_label import ElementWithLabel
            sage: e = ElementWithLabel(1, a_1)
            sage: latex(e)
            a_{1}
        """
        return latex(self.label)

    def __str__(self):
        """
        Return the string representation of ``self``, which is just
        the string representation of the label.

        TESTS::

            sage: var('a_1')
            a_1
            sage: from sage.misc.element_with_label import ElementWithLabel
            sage: e = ElementWithLabel(1, a_1)
            sage: str(e)
            'a_1'
        """
        return str(self.label)

    def __repr__(self):
        """
        Return the representation of ``self``, which is just
        the representation of the label.

        TESTS::

            sage: var('a_1')
            a_1
            sage: from sage.misc.element_with_label import ElementWithLabel
            sage: e = ElementWithLabel(1, a_1)
            sage: repr(e)
            'a_1'
        """
        return repr(self.label)

    def __hash__(self):
        """
        Return the hash of the labeled element ``self``,
        which is constructed from hashes of both constituents.

        TESTS::

            sage: from sage.misc.element_with_label import ElementWithLabel
            sage: a = ElementWithLabel(1, 'a')
            sage: b = ElementWithLabel(1, 'b')
            sage: d = {}
            sage: d[a] = 'element 1'
            sage: d[b] = 'element 2'
            sage: d
            {'a': 'element 1', 'b': 'element 2'}
            sage: a = ElementWithLabel("a", [2,3])
            sage: hash(a)
            1853891946828512984
        """
        try:
            return hash((hash(self.element), hash(self.label)))
        except TypeError:
            return hash((repr(self.element), repr(self.label)))

    def __eq__(self, other):
        """
        Two labeled elements are equal if and only if both of their
        constituents are equal.

        TESTS::

            sage: from sage.misc.element_with_label import ElementWithLabel
            sage: a = ElementWithLabel(1, 'a')
            sage: b = ElementWithLabel(1, 'b')
            sage: x = ElementWithLabel(1, 'a')
            sage: a == b
            False
            sage: a == x
            True
        """
        return self.element == other.element and self.label == other.label

    def __ne__(self, other):
        """
        Two labeled elements are not equal if and only if first or second
        constituents are not equal.

        TESTS::

            sage: from sage.misc.element_with_label import ElementWithLabel
            sage: a = ElementWithLabel(1, 'a')
            sage: b = ElementWithLabel(1, 'b')
            sage: x = ElementWithLabel(1, 'a')
            sage: a != b
            True
            sage: a != x
            False        
        """
        return not(self == other)

