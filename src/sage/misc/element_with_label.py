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
    like ``element``.
    For any presentation purposes it appears just as ``label`` would.

    TESTS::

    sage: P = Poset({1: [2,3]})
    sage: labs = {i: P.rank(i) for i in range(1, 4)}
    sage: print(labs)
    {1: 0, 2: 1, 3: 1}
    sage: print(P.plot(element_labels=labs))
    Graphics object consisting of 6 graphics primitives
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
        which is just the hash of ``self.element``

        TESTS::

            sage: from sage.misc.element_with_label import ElementWithLabel
            sage: a = ElementWithLabel(1, 'a')
            sage: b = ElementWithLabel(1, 'b')
            sage: d = {}
            sage: d[a] = 'element 1'
            sage: d[b] = 'element 2'
            sage: d
            {'a': 'element 2'}
            sage: a = ElementWithLabel("a", [2,3])
            sage: hash(a)
            12416037344
            sage: hash("a")
            12416037344
        """
        return hash(self.element)

    def __eq__(self, other):
        """
        Two labeled elements are equal if and only if their
        elements are equal.

        TESTS::

            sage: from sage.misc.element_with_label import ElementWithLabel
            sage: a = ElementWithLabel(1, 'a')
            sage: b = ElementWithLabel(1, 'b')
            sage: x = ElementWithLabel(2, 'a')
            sage: a == b
            True
            sage: a == x
            False
        """
        return self.element == other.element

    def __ne__(self, other):
        """
        Two labeled elements are not equal if and only if 
        the elements they are wrapping are not equal.

        TESTS::

            sage: from sage.misc.element_with_label import ElementWithLabel
            sage: a = ElementWithLabel(1, 'a')
            sage: b = ElementWithLabel(1, 'b')
            sage: x = ElementWithLabel(2, 'a')
            sage: a != b
            False
            sage: a != x
            True
        """
        return not(self == other)

