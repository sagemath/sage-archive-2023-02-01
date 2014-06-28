"""
Finitely Presented Groups

Finitely presented groups are constructed as quotients of
:mod:`~sage.groups.free_group`::

    sage: F.<a,b,c> = FreeGroup()
    sage: G = F / [a^2, b^2, c^2, a*b*c*a*b*c]
    sage: G
    Finitely presented group < a, b, c | a^2, b^2, c^2, (a*b*c)^2 >

One can create their elements by mutiplying the generators or by
specifying a Tietze list (see
:meth:`~sage.groups.finitely_presented.FinitelyPresentedGroupElement.Tietze`)
as in the case of free groups::

    sage: G.gen(0) * G.gen(1)
    a*b
    sage: G([1,2,-1])
    a*b*a^-1
    sage: a.parent()
    Free Group on generators {a, b, c}
    sage: G.inject_variables()
    Defining a, b, c
    sage: a.parent()
    Finitely presented group < a, b, c | a^2, b^2, c^2, (a*b*c)^2 >

Notice that, even if they are represented in the same way, the
elements of a finitely presented group and the elements of the
corresponding free group are not the same thing.  However, they can be
converted from one parent to the other::

    sage: F.<a,b,c> = FreeGroup()
    sage: G = F / [a^2,b^2,c^2,a*b*c*a*b*c]
    sage: F([1])
    a
    sage: G([1])
    a
    sage: F([1]) is G([1])
    False
    sage: F([1]) == G([1])
    False
    sage: G(a*b/c)
    a*b*c^-1
    sage: F(G(a*b/c))
    a*b*c^-1

Finitely presented groups are implemented via GAP. You can use the
:meth:`~sage.groups.libgap_wrapper.ParentLibGAP.gap` method to access
the underlying LibGAP object::

    sage: G = FreeGroup(2)
    sage: G.inject_variables()
    Defining x0, x1
    sage: H = G / (x0^2, (x0*x1)^2, x1^2)
    sage: H.gap()
    <fp group on the generators [ x0, x1 ]>

This can be useful, for example, to use GAP functions that are not yet
wrapped in Sage::

    sage: H.gap().LowerCentralSeries()
    [ Group(<fp, no generators known>), Group(<fp, no generators known>) ]

The same holds for the group elements::

    sage: G = FreeGroup(2)
    sage: H = G / (G([1, 1]), G([2, 2, 2]), G([1, 2, -1, -2]));  H
    Finitely presented group < x0, x1 | x0^2, x1^3, x0*x1*x0^-1*x1^-1 >
    sage: a = H([1])
    sage: a
    x0
    sage: a.gap()
    x0
    sage: a.gap().Order()
    2
    sage: type(_)    # note that the above output is not a Sage integer
    <type 'sage.libs.gap.element.GapElement_Integer'>

You can use call syntax to replace the generators with a set of
arbitrary ring elements. For example, take the free abelian group
obtained by modding out the commutator subgroup of the free group::

    sage: G = FreeGroup(2)
    sage: G_ab = G / [G([1, 2, -1, -2])];  G_ab
    Finitely presented group < x0, x1 | x0*x1*x0^-1*x1^-1 >
    sage: a,b = G_ab.gens()
    sage: g =  a * b
    sage: M1 = matrix([[1,0],[0,2]])
    sage: M2 = matrix([[0,1],[1,0]])
    sage: g(3, 5)
    15
    sage: g(M1, M1)
    [1 0]
    [0 4]
    sage: M1*M2 == M2*M1   # matrices do not commute
    False
    sage: g(M1, M2)
    Traceback (most recent call last):
    ...
    ValueError: the values do not satisfy all relations of the group

.. WARNING::

    Some methods are not guaranteed to finish since the word problem
    for finitely presented groups is, in general, undecidable. In
    those cases the process may run unil the available memory is
    exhausted.

REFERENCES:

- :wikipedia:`Presentation_of_a_group`

- :wikipedia:`Word_problem_for_groups`

AUTHOR:

- Miguel Angel Marco Buzunariz
"""

##############################################################################
#       Copyright (C) 2012 Miguel Angel Marco Buzunariz <mmarco@unizar.es>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#
#  The full text of the GPL is available at:
#
#                  http://www.gnu.org/licenses/
##############################################################################


from sage.groups.group import Group
from sage.groups.libgap_wrapper import ParentLibGAP, ElementLibGAP
from sage.groups.libgap_mixin import GroupMixinLibGAP
from sage.structure.unique_representation import UniqueRepresentation
from sage.libs.gap.libgap import libgap
from sage.libs.gap.element import GapElement
from sage.rings.integer import Integer
from sage.rings.integer_ring import IntegerRing
from sage.misc.cachefunc import cached_method
from sage.groups.free_group import FreeGroupElement

from sage.structure.element import Element, MultiplicativeGroupElement
from sage.interfaces.gap import gap
from sage.rings.integer import Integer
from sage.rings.integer_ring import IntegerRing
from sage.functions.generalized import sign
from sage.matrix.constructor import matrix

class FinitelyPresentedGroupElement(FreeGroupElement):
    """
    A wrapper of GAP's Finitely Presented Group elements.

    The elements are created by passing the Tietze list that determines them.

    EXAMPLES::

        sage: G = FreeGroup('a, b')
        sage: H = G / [G([1]), G([2, 2, 2])]
        sage: H([1, 2, 1, -1])
        a*b
        sage: H([1, 2, 1, -2])
        a*b*a*b^-1
        sage: x = H([1, 2, -1, -2])
        sage: x
        a*b*a^-1*b^-1
        sage: y = H([2, 2, 2, 1, -2, -2, -2])
        sage: y
        b^3*a*b^-3
        sage: x*y
        a*b*a^-1*b^2*a*b^-3
        sage: x^(-1)
        b*a*b^-1*a^-1
    """

    def __init__(self, parent, x, check=True):
        """
        The Python constructor.

        See :class:`FinitelyPresentedGroupElement` for details.

        TESTS::

            sage: G = FreeGroup('a, b')
            sage: H = G / [G([1]), G([2, 2, 2])]
            sage: H([1, 2, 1, -1])
            a*b

            sage: TestSuite(G).run()
            sage: TestSuite(H).run()

            sage: G.<a,b> = FreeGroup()
            sage: H = G / (G([1]), G([2, 2, 2]))
            sage: x = H([1, 2, -1, -2])
            sage: TestSuite(x).run()
            sage: TestSuite(G.one()).run()
        """
        if not isinstance(x, GapElement):
            F = parent.free_group()
            free_element = F(x)
            fp_family = parent.gap().Identity().FamilyObj()
            x = libgap.ElementOfFpGroup(fp_family, free_element.gap())
        ElementLibGAP.__init__(self, parent, x)

    def __reduce__(self):
        """
        Used in pickling.

        TESTS::

            sage: F.<a,b> = FreeGroup()
            sage: G = F / [a*b, a^2]
            sage: G.inject_variables()
            Defining a, b
            sage: a.__reduce__()
            (Finitely presented group < a, b | a*b, a^2 >, ((1,),))
            sage: (a*b*a^-1).__reduce__()
            (Finitely presented group < a, b | a*b, a^2 >, ((1, 2, -1),))

            sage: F.<a,b,c> = FreeGroup('a, b, c')
            sage: G = F.quotient([a*b*c/(b*c*a), a*b*c/(c*a*b)])
            sage: G.inject_variables()
            Defining a, b, c
            sage: x = a*b*c
            sage: x.__reduce__()
            (Finitely presented group < a, b, c | a*b*c*a^-1*c^-1*b^-1, a*b*c*b^-1*a^-1*c^-1 >,
             ((1, 2, 3),))
        """
        return (self.parent(), tuple([self.Tietze()]))

    def _repr_(self):
        """
        Return a string representation.

        OUTPUT:

        String.

        EXAMPLES::

            sage: G.<a,b> = FreeGroup()
            sage: H = G / [a^2, b^3]
            sage: H.gen(0)
            a
            sage: H.gen(0)._repr_()
            'a'
            sage: H.one()
            1
        """
        # computing that an element is actually one can be very expensive
        if self.Tietze() == ():
            return '1'
        else:
            return self.gap()._repr_()

    @cached_method
    def Tietze(self):
        """
        Return the Tietze list of the element.

        The Tietze list of a word is a list of integers that represent
        the letters in the word.  A positive integer `i` represents
        the letter corresponding to the `i`-th generator of the group.
        Negative integers represent the inverses of generators.

        OUTPUT:

        A tuple of integers.

        EXAMPLES::

            sage: G = FreeGroup('a, b')
            sage: H = G / (G([1]), G([2, 2, 2]))
            sage: H.inject_variables()
            Defining a, b
            sage: a.Tietze()
            (1,)
            sage: x = a^2*b^(-3)*a^(-2)
            sage: x.Tietze()
            (1, 1, -2, -2, -2, -1, -1)
        """
        tl = self.gap().UnderlyingElement().TietzeWordAbstractWord()
        return tuple(tl.sage())

    def __call__(self, *values, **kwds):
        """
        Replace the generators of the free group with ``values``.

        INPUT:

        - ``*values`` -- a list/tuple/iterable of the same length as
          the number of generators.

        - ``check=True`` -- boolean keyword (default:
          ``True``). Whether to verify that ``values`` satisfy the
          relations in the finitely presented group.

        OUTPUT:

        The product of ``values`` in the order and with exponents
        specified by ``self``.

        EXAMPLES::

            sage: G.<a,b> = FreeGroup()
            sage: H = G / [a/b];  H
            Finitely presented group < a, b | a*b^-1 >
            sage: H.simplified()
            Finitely presented group < a |  >

        The generator `b` can be eliminated using the relation `a=b`. Any
        values that you plug into a word must satisfy this relation::

            sage: A, B = H.gens()
            sage: w = A^2 * B
            sage: w(2,2)
            8
            sage: w(3,3)
            27
            sage: w(1,2)
            Traceback (most recent call last):
            ...
            ValueError: the values do not satisfy all relations of the group
            sage: w(1, 2, check=False)    # result depends on presentation of the group element
            2
        """
        values = list(values)
        if kwds.get('check', True):
            for rel in self.parent().relations():
                rel = rel(values)
                if rel != 1:
                    raise ValueError('the values do not satisfy all relations of the group')
        return super(FinitelyPresentedGroupElement, self).__call__(values)


def wrap_FpGroup(libgap_fpgroup):
    """
    Wrap a GAP finitely presented group.

    This function changes the comparison method of
    ``libgap_free_group`` to comparison by Python ``id``. If you want
    to put the LibGAP free group into a container ``(set, dict)`` then you
    should understand the implications of
    :meth:`~sage.libs.gap.element.GapElement._set_compare_by_id`. To
    be safe, it is recommended that you just work with the resulting
    Sage :class:`FinitelyPresentedGroup`.

    INPUT:

    - ``libgap_fpgroup`` -- a LibGAP finitely presented group

    OUTPUT:

    A Sage :class:`FinitelyPresentedGroup`.

    EXAMPLES:

    First construct a LibGAP finitely presented group::

        sage: F = libgap.FreeGroup(['a', 'b'])
        sage: a_cubed = F.GeneratorsOfGroup()[0] ^ 3
        sage: P = F / libgap([ a_cubed ]);   P
        <fp group of size infinity on the generators [ a, b ]>
        sage: type(P)
        <type 'sage.libs.gap.element.GapElement'>

    Now wrap it::

        sage: from sage.groups.finitely_presented import wrap_FpGroup
        sage: wrap_FpGroup(P)
        Finitely presented group < a, b | a^3 >
    """
    assert libgap_fpgroup.IsFpGroup()
    libgap_fpgroup._set_compare_by_id()
    from sage.groups.free_group import wrap_FreeGroup
    free_group = wrap_FreeGroup(libgap_fpgroup.FreeGroupOfFpGroup())
    relations = tuple( free_group(rel.UnderlyingElement())
                       for rel in libgap_fpgroup.RelatorsOfFpGroup() )
    return FinitelyPresentedGroup(free_group, relations)


class RewritingSystem(object):
    """
    A class that wraps GAP's rewriting systems.

    A rewriting system is a set of rules that allow to transform
    one word in the group to an equivalent one.

    If the rewriting system is confluent, then the transformated
    word is a unique reduced form of the element of the group.

    .. WARNING::

        Note that the process of making a rewriting system confluent
        might not end.

    INPUT:

    - ``G`` -- a group

    REFERENCES:

    - :wikipedia:`Knuth-Bendix_completion_algorithm`

    EXAMPLES::

        sage: F.<a,b> = FreeGroup()
        sage: G = F / [a*b/a/b]
        sage: k = G.rewriting_system()
        sage: k
        Rewriting system of Finitely presented group < a, b | a*b*a^-1*b^-1 >
        with rules:
            a*b*a^-1*b^-1    --->    1

        sage: k.reduce(a*b*a*b)
        (a*b)^2
        sage: k.make_confluent()
        sage: k
        Rewriting system of Finitely presented group < a, b | a*b*a^-1*b^-1 >
        with rules:
            b^-1*a^-1    --->    a^-1*b^-1
            b^-1*a    --->    a*b^-1
            b*a^-1    --->    a^-1*b
            b*a    --->    a*b

        sage: k.reduce(a*b*a*b)
        a^2*b^2

    .. TODO::

        - Include support for different orderings (currently only shortlex
          is used).

        - Include the GAP package kbmag for more functionalities, including
          automatic structures and faster compiled functions.

    AUTHORS:

    - Miguel Angel Marco Buzunariz (2013-12-16)
    """
    def __init__(self, G):
        """
        Initialize ``self``.

        EXAMPLES::

            sage: F.<a,b,c> = FreeGroup()
            sage: G = F / [a^2, b^3, c^5]
            sage: k = G.rewriting_system()
            sage: k
            Rewriting system of Finitely presented group < a, b, c | a^2, b^3, c^5 >
            with rules:
                a^2    --->    1
                b^3    --->    1
                c^5    --->    1
        """
        self._free_group = G.free_group()
        self._fp_group = G
        self._fp_group_gap = G.gap()
        self._monoid_isomorphism = self._fp_group_gap.IsomorphismFpMonoid()
        self._monoid = self._monoid_isomorphism.Image()
        self._gap = self._monoid.KnuthBendixRewritingSystem()

    def __repr__(self):
        """
        Return a string representation.

        EXAMPLES::

            sage: F.<a> = FreeGroup()
            sage: G = F / [a^2]
            sage: k = G.rewriting_system()
            sage: k
            Rewriting system of Finitely presented group < a | a^2 >
            with rules:
                a^2    --->    1
        """
        ret = "Rewriting system of {}\nwith rules:".format(self._fp_group)
        for i in sorted(self.rules().items()): # Make sure they are sorted to the repr is unique
            ret += "\n    {}    --->    {}".format(i[0], i[1])
        return ret

    def free_group(self):
        """
        The free group after which the rewriting system is defined

        EXAMPLES::

            sage: F = FreeGroup(3)
            sage: G = F / [ [1,2,3], [-1,-2,-3] ]
            sage: k = G.rewriting_system()
            sage: k.free_group()
            Free Group on generators {x0, x1, x2}
        """
        return self._free_group

    def finitely_presented_group(self):
        """
        The finitely presented group where the rewriting system is defined.

        EXAMPLES::

            sage: F = FreeGroup(3)
            sage: G = F / [ [1,2,3], [-1,-2,-3], [1,1], [2,2] ]
            sage: k = G.rewriting_system()
            sage: k.make_confluent()
            sage: k
            Rewriting system of Finitely presented group < x0, x1, x2 | x0*x1*x2, x0^-1*x1^-1*x2^-1, x0^2, x1^2 >
            with rules:
                x0^-1    --->    x0
                x1^-1    --->    x1
                x2^-1    --->    x2
                x0^2    --->    1
                x0*x1    --->    x2
                x0*x2    --->    x1
                x1*x0    --->    x2
                x1^2    --->    1
                x1*x2    --->    x0
                x2*x0    --->    x1
                x2*x1    --->    x0
                x2^2    --->    1
            sage: k.finitely_presented_group()
            Finitely presented group < x0, x1, x2 | x0*x1*x2, x0^-1*x1^-1*x2^-1, x0^2, x1^2 >
        """
        return self._fp_group

    def reduce(self, element):
        """
        Applies the rules in the rewriting system to the element, to obtain
        a reduced form.

        If the rewriting system is confluent, this reduced form is unique
        for all words representing the same element.

        EXAMPLES::

            sage: F.<a,b> = FreeGroup()
            sage: G = F/[a^2, b^3, (a*b/a)^3, b*a*b*a]
            sage: k = G.rewriting_system()
            sage: k.reduce(b^4)
            b
            sage: k.reduce(a*b*a)
            a*b*a
        """
        eg = self._fp_group(element).gap()
        egim = self._monoid_isomorphism.Image(eg)
        red = self.gap().ReducedForm(egim.UnderlyingElement())
        redfpmon = self._monoid.One().FamilyObj().ElementOfFpMonoid(red)
        reducfpgr = self._monoid_isomorphism.PreImagesRepresentative(redfpmon)
        tz = reducfpgr.UnderlyingElement().TietzeWordAbstractWord(self._free_group.gap().GeneratorsOfGroup())
        return self._fp_group(tz.sage())

    def gap(self):
        """
        The gap representation of the rewriting system.

        EXAMPLES::

            sage: F.<a,b>=FreeGroup()
            sage: G=F/[a*a,b*b]
            sage: k=G.rewriting_system()
            sage: k.gap()
            Knuth Bendix Rewriting System for Monoid( [ a, A, b, B ] ) with rules
            [ [ a^2, <identity ...> ], [ a*A, <identity ...> ],
              [ A*a, <identity ...> ], [ b^2, <identity ...> ],
              [ b*B, <identity ...> ], [ B*b, <identity ...> ] ]
        """
        return self._gap

    def rules(self):
        """
        Return the rules that form the rewritig system.

        OUTPUT:

        A dictionary containing the rules of the rewriting system.
        Each key is a word in the free group, and its corresponding
        value is the word to which it is reduced.

        EXAMPLES::

            sage: F.<a,b> = FreeGroup()
            sage: G = F / [a*a*a,b*b*a*a]
            sage: k = G.rewriting_system()
            sage: k
            Rewriting system of Finitely presented group < a, b | a^3, b^2*a^2 >
            with rules:
                a^3    --->    1
                b^2*a^2    --->    1

            sage: k.rules()
            {b^2*a^2: 1, a^3: 1}
            sage: k.make_confluent()
            sage: sorted(k.rules().items())
            [(a^-2, a), (a^-1*b^-1, a*b), (a^-1*b, b^-1), (a^2, a^-1),
             (a*b^-1, b), (b^-1*a^-1, a*b), (b^-1*a, b), (b^-2, a^-1),
             (b*a^-1, b^-1), (b*a, a*b), (b^2, a)]
        """
        dic = {}
        grules = self.gap().Rules()
        for i in grules:
            a, b = i
            afpmon = self._monoid.One().FamilyObj().ElementOfFpMonoid(a)
            afg = self._monoid_isomorphism.PreImagesRepresentative(afpmon)
            atz = afg.UnderlyingElement().TietzeWordAbstractWord(self._free_group.gap().GeneratorsOfGroup())
            af = self._free_group(atz.sage())
            if len(af.Tietze()) != 0:
                bfpmon = self._monoid.One().FamilyObj().ElementOfFpMonoid(b)
                bfg = self._monoid_isomorphism.PreImagesRepresentative(bfpmon)
                btz = bfg.UnderlyingElement().TietzeWordAbstractWord(self._free_group.gap().GeneratorsOfGroup())
                bf = self._free_group(btz.sage())
                dic[af]=bf
        return dic

    def is_confluent(self):
        """
        Return ``True`` if the system is confluent and ``False`` otherwise.

        EXAMPLES::

            sage: F = FreeGroup(3)
            sage: G = F / [F([1,2,1,2,1,3,-1]),F([2,2,2,1,1,2]),F([1,2,3])]
            sage: k = G.rewriting_system()
            sage: k.is_confluent()
            False
            sage: k
            Rewriting system of Finitely presented group < x0, x1, x2 | (x0*x1)^2*x0*x2*x0^-1, x1^3*x0^2*x1, x0*x1*x2 >
            with rules:
                x0*x1*x2    --->    1
                x1^3*x0^2*x1    --->    1
                (x0*x1)^2*x0*x2*x0^-1    --->    1

            sage: k.make_confluent()
            sage: k.is_confluent()
            True
            sage: k
            Rewriting system of Finitely presented group < x0, x1, x2 | (x0*x1)^2*x0*x2*x0^-1, x1^3*x0^2*x1, x0*x1*x2 >
            with rules:
                x0^-1    --->    x0
                x1^-1    --->    x1
                x0^2    --->    1
                x0*x1    --->    x2^-1
                x0*x2^-1    --->    x1
                x1*x0    --->    x2
                x1^2    --->    1
                x1*x2^-1    --->    x0*x2
                x1*x2    --->    x0
                x2^-1*x0    --->    x0*x2
                x2^-1*x1    --->    x0
                x2^-2    --->    x2
                x2*x0    --->    x1
                x2*x1    --->    x0*x2
                x2^2    --->    x2^-1
        """
        return self._gap.IsConfluent().sage()

    def make_confluent(self):
        """
        Applies Knuth-Bendix algorithm to try to transform the rewriting
        system into a confluent one.

        Note that this method does not return any object, just changes the
        rewriting sytem internally.

        .. WARNING:

            This algorithm is not granted to finish. Although it may be useful
            in some occasions to run it, interrupt it manually after some time
            and use then the transformed rewriting system. Even if it is not
            confluent, it could be used to reduce some words.

        ALGORITHM:

        Uses GAP's ``MakeConfluent``.

        EXAMPLES::

            sage: F.<a,b> = FreeGroup()
            sage: G = F / [a^2,b^3,(a*b/a)^3,b*a*b*a]
            sage: k = G.rewriting_system()
            sage: k
            Rewriting system of Finitely presented group < a, b | a^2, b^3, a*b^3*a^-1, (b*a)^2 >
            with rules:
                a^2    --->    1
                b^3    --->    1
                (b*a)^2    --->    1
                a*b^3*a^-1    --->    1

            sage: k.make_confluent()
            sage: k
            Rewriting system of Finitely presented group < a, b | a^2, b^3, a*b^3*a^-1, (b*a)^2 >
            with rules:
                a^-1    --->    a
                a^2    --->    1
                b^-1*a    --->    a*b
                b^-2    --->    b
                b*a    --->    a*b^-1
                b^2    --->    b^-1
        """
        try:
            self._gap.MakeConfluent()
        except ValueError:
            raise ValueError('could not make the system confluent')

class FinitelyPresentedGroup(GroupMixinLibGAP, UniqueRepresentation,
    Group, ParentLibGAP):
    """
    A class that wraps GAP's Finitely Presented Groups.

    .. WARNING::

        You should use
        :meth:`~sage.groups.free_group.FreeGroup_class.quotient` to
        construct finitely presented groups as quotients of free
        groups.

    EXAMPLES::

        sage: G.<a,b> = FreeGroup()
        sage: H = G / [a, b^3]
        sage: H
        Finitely presented group < a, b | a, b^3 >
        sage: H.gens()
        (a, b)

        sage: F.<a,b> = FreeGroup('a, b')
        sage: J = F / (F([1]), F([2, 2, 2]))
        sage: J is H
        True

        sage: G = FreeGroup(2)
        sage: H = G / (G([1, 1]), G([2, 2, 2]))
        sage: H.gens()
        (x0, x1)
        sage: H.gen(0)
        x0
        sage: H.ngens()
        2
        sage: H.gap()
        <fp group on the generators [ x0, x1 ]>
        sage: type(_)
        <type 'sage.libs.gap.element.GapElement'>
    """
    Element = FinitelyPresentedGroupElement

    def __init__(self, free_group, relations):
        """
        The Python constructor.

        TESTS::

            sage: G = FreeGroup('a, b')
            sage: H = G / (G([1]), G([2])^3)
            sage: H
            Finitely presented group < a, b | a, b^3 >

            sage: F = FreeGroup('a, b')
            sage: J = F / (F([1]), F([2, 2, 2]))
            sage: J is H
            True

            sage: TestSuite(H).run()
            sage: TestSuite(J).run()
        """
        from sage.groups.free_group import is_FreeGroup
        assert is_FreeGroup(free_group)
        assert isinstance(relations, tuple)
        self._free_group = free_group
        self._relations = relations
        self._assign_names(free_group.variable_names())
        parent_gap = free_group.gap() / libgap([ rel.gap() for rel in relations])
        ParentLibGAP.__init__(self, parent_gap)
        Group.__init__(self)

    def _repr_(self):
        """
        Return a string representation.

        OUTPUT:

        String.

        TESTS::

            sage: G.<a,b> = FreeGroup()
            sage: H = G / (G([1]), G([2])^3)
            sage: H  # indirect doctest
            Finitely presented group < a, b | a, b^3 >
            sage: H._repr_()
            'Finitely presented group < a, b | a, b^3 >'
        """
        gens = ', '.join(self.variable_names())
        rels = ', '.join([ str(r) for r in self.relations() ])
        return 'Finitely presented group ' + '< '+ gens + ' | ' + rels + ' >'

    def _latex_(self):
        """
        Return a LaTeX representation

        OUTPUT:

        String. A valid LaTeX math command sequence.

        TESTS::

            sage: F=FreeGroup(4)
            sage: F.inject_variables()
            Defining x0, x1, x2, x3
            sage: G=F.quotient([x0*x2, x3*x1*x3, x2*x1*x2])
            sage: G._latex_()
            '\\langle x_{0}, x_{1}, x_{2}, x_{3} \\mid x_{0}\\cdot x_{2} , x_{3}\\cdot x_{1}\\cdot x_{3} , x_{2}\\cdot x_{1}\\cdot x_{2}\\rangle'
        """
        r = '\\langle '
        for i in range(self.ngens()):
            r = r+self.gen(i)._latex_()
            if i < self.ngens()-1:
                r = r+', '
        r = r+' \\mid '
        for i in range(len(self._relations)):
            r = r+(self._relations)[i]._latex_()
            if i < len(self.relations())-1:
                r = r+' , '
        r = r+'\\rangle'
        return r

    def free_group(self):
        """
        Return the free group (without relations).

        OUTPUT:

        A :func:`~sage.groups.free_group.FreeGroup`.

        EXAMPLES::

            sage: G.<a,b,c> = FreeGroup()
            sage: H = G / (a^2, b^3, a*b*~a*~b)
            sage: H.free_group()
            Free Group on generators {a, b, c}
            sage: H.free_group() is G
            True
        """
        return self._free_group

    def relations(self):
        """
        Return the relations of the group.

        OUTPUT:

        The relations as a tuple of elements of :meth:`free_group`.

        EXAMPLES::

            sage: F = FreeGroup(5, 'x')
            sage: F.inject_variables()
            Defining x0, x1, x2, x3, x4
            sage: G = F.quotient([x0*x2, x3*x1*x3, x2*x1*x2])
            sage: G.relations()
            (x0*x2, x3*x1*x3, x2*x1*x2)
            sage: all(rel in F for rel in G.relations())
            True
        """
        return self._relations

    @cached_method
    def cardinality(self, limit=4096000):
        """
        Compute the cardinality of ``self``.

        INPUT:

        - ``limit`` -- integer (default: 4096000). The maximal number
          of cosets before the computation is aborted.

        OUTPUT:

        Integer or ``Infinity``. The number of elements in the group.

        EXAMPLES::

            sage: G.<a,b> = FreeGroup('a, b')
            sage: H = G / (a^2, b^3, a*b*~a*~b)
            sage: H.cardinality()
            6

            sage: F.<a,b,c> = FreeGroup()
            sage: J = F / (F([1]), F([2, 2, 2]))
            sage: J.cardinality()
            +Infinity

        ALGORITHM:

            Uses GAP.

        .. WARNING::

            This is in general not a decidable problem, so it is not
            guaranteed to give an answer. If the group is infinite, or
            too big, you should be prepared for a long computation
            that consumes all the memory without finishing if you do
            not set a sensible ``limit``.
        """
        with libgap.global_context('CosetTableDefaultMaxLimit', limit):
            if not libgap.IsFinite(self.gap()):
                from sage.rings.infinity import Infinity
                return Infinity
            try:
                size = self.gap().Size()
            except ValueError:
                raise ValueError('Coset enumeration ran out of memory, is the group finite?')
        return size.sage()

    order = cardinality

    def as_permutation_group(self, limit=4096000):
        """
        Return an isomorphic permutation group.

        The generators of the resulting group correspond to the images
        by the isomorphism of the generators of the given group.

        INPUT:

        - ``limit`` -- integer (default: 4096000). The maximal number
          of cosets before the computation is aborted.

        OUTPUT:

        A Sage
        :func:`~sage.groups.perm_gps.permgroup.PermutationGroup`. If
        the number of cosets exceeds the given ``limit``, a
        ``ValueError`` is returned.

        EXAMPLES::

            sage: G.<a,b> = FreeGroup()
            sage: H = G / (a^2, b^3, a*b*~a*~b)
            sage: H.as_permutation_group()
            Permutation Group with generators [(1,2)(3,5)(4,6), (1,3,4)(2,5,6)]

            sage: G.<a,b> = FreeGroup()
            sage: H = G / [a^3*b]
            sage: H.as_permutation_group(limit=1000)
            Traceback (most recent call last):
            ...
            ValueError: Coset enumeration exceeded limit, is the group finite?

        ALGORITHM:

            Uses GAP's coset enumeration on the trivial subgroup.

        .. WARNING::

            This is in general not a decidable problem (in fact, it is
            not even posible to check if the group is finite or
            not). If the group is infinite, or too big, you should be
            prepared for a long computation that consumes all the
            memory without finishing if you do not set a sensible
            ``limit``.
        """
        with libgap.global_context('CosetTableDefaultMaxLimit', limit):
            try:
                trivial_subgroup = self.gap().TrivialSubgroup()
                coset_table = self.gap().CosetTable(trivial_subgroup).sage()
            except ValueError:
                raise ValueError('Coset enumeration exceeded limit, is the group finite?')
        from sage.combinat.permutation import Permutation
        from sage.groups.perm_gps.permgroup import PermutationGroup
        return PermutationGroup([
                Permutation(coset_table[2*i]) for i in range(len(coset_table)//2)])

    def direct_product(self, H, reduced=False, new_names=True):
        r"""
        Return the direct product of ``self`` with finitely presented
        group ``H``.

        Calls GAP function ``DirectProduct``, which returns the direct
        product of a list of groups of any representation.

        From [JohnsonPG90]_ (pg 45, proposition 4): If `G`, `H` are groups
        presented by `\langle X \mid R \rangle` and `\langle Y \mid S \rangle`
        respectively, then their direct product has the presentation
        `\langle X, Y \mid R, S, [X, Y] \rangle` where `[X, Y]` denotes the
        set of commutators `\{ x^{-1} y^{-1} x y \mid x \in X, y \in Y \}`.

        INPUT:

        - ``H`` -- a finitely presented group

        - ``reduced`` -- (default: ``False``) boolean; if ``True``, then
          attempt to reduce the presentation of the product group

        - ``new_names`` -- (default: ``True``) boolean; If ``True``, then
          lexicographical variable names are assigned to the generators of
          the group to be returned. If ``False``, the group to be returned
          keeps the generator names of the two groups forming the direct
          product. Note that one cannot ask to reduce the output and ask
          to keep the old variable names, as they they may change meaning
          in the output group if its presentation is reduced.

        OUTPUT:

        The direct product of ``self`` with ``H`` as a finitely
        presented group.

        EXAMPLES::

            sage: G = FreeGroup()
            sage: C12 =  ( G / [G([1,1,1,1])] ).direct_product( G / [G([1,1,1])]); C12
            Finitely presented group < a, b | a^4, b^3, a^-1*b^-1*a*b >
            sage: C12.order(), C12.as_permutation_group().is_cyclic()
            (12, True)
            sage: klein = ( G / [G([1,1])] ).direct_product( G / [G([1,1])]); klein
            Finitely presented group < a, b | a^2, b^2, a^-1*b^-1*a*b >
            sage: klein.order(), klein.as_permutation_group().is_cyclic()
            (4, False)

        We can keep the variable names from ``self`` and ``H`` to examine how
        new relations are formed::

            sage: F = FreeGroup("a"); G = FreeGroup("g")
            sage: X = G / [G.0^12]; A = F / [F.0^6]
            sage: X.direct_product(A, new_names=False)
            Finitely presented group < g, a | g^12, a^6, g^-1*a^-1*g*a >
            sage: A.direct_product(X, new_names=False)
            Finitely presented group < a, g | a^6, g^12, a^-1*g^-1*a*g >

        Or we can attempt to reduce the output group presentation::

            sage: F = FreeGroup("a"); G = FreeGroup("g")
            sage: X = G / [G.0]; A = F / [F.0]
            sage: X.direct_product(A, new_names=True)
            Finitely presented group < a, b | a, b, a^-1*b^-1*a*b >
            sage: X.direct_product(A, reduced=True, new_names=True)
            Finitely presented group <  |  >

        But we cannot do both::

            sage: K = FreeGroup(['a','b'])
            sage: D = K / [K.0^5, K.1^8]
            sage: D.direct_product(D, reduced=True, new_names=False)
            Traceback (most recent call last):
            ...
            ValueError: cannot reduce output and keep old variable names

        TESTS::

            sage: G = FreeGroup()
            sage: Dp = (G / [G([1,1])]).direct_product( G / [G([1,1,1,1,1,1])] )
            sage: Dp.as_permutation_group().is_isomorphic(PermutationGroup(['(1,2)','(3,4,5,6,7,8)']))
            True
            sage: C7 = G / [G.0**7]; C6 =  G / [G.0**6]
            sage: C14 = G / [G.0**14]; C3 =  G / [G.0**3]
            sage: C7.direct_product(C6).is_isomorphic(C14.direct_product(C3))
            True
            sage: F = FreeGroup(2); D = F / [F([1,1,1,1,1]),F([2,2]),F([1,2])**2]
            sage: D.direct_product(D).as_permutation_group().is_isomorphic(
            ....: direct_product_permgroups([DihedralGroup(5),DihedralGroup(5)]))
            True

        AUTHORS:

        - Davis Shurbert (2013-07-20): initial version

        REFERENCES:

        .. [JohnsonPG90] D.L. Johnson. *Presentations of Groups*.
           Cambridge University Press. (1990).
        """
        from sage.groups.free_group import FreeGroup, _lexi_gen

        if not isinstance(H, FinitelyPresentedGroup):
            raise TypeError("input must be a finitely presented group")
        if reduced and not new_names:
            raise ValueError("cannot reduce output and keep old variable names")

        fp_product = libgap.DirectProduct([self.gap(), H.gap()])
        GAP_gens = fp_product.FreeGeneratorsOfFpGroup()
        if new_names:
            name_itr = _lexi_gen() # Python generator for lexicographical variable names
            gen_names = [name_itr.next() for i in GAP_gens]
        else:
            gen_names= [str(g) for g in self.gens()] + [str(g) for g in H.gens()]
        # Build the direct product in Sage for better variable names
        ret_F = FreeGroup(gen_names)
        ret_rls = tuple([ret_F(rel_word.TietzeWordAbstractWord(GAP_gens).sage())
            for rel_word in fp_product.RelatorsOfFpGroup()])
        ret_fpg = FinitelyPresentedGroup(ret_F, ret_rls)
        if reduced:
            ret_fpg = ret_fpg.simplified()
        return ret_fpg

    def semidirect_product(self, H, hom, check=True, reduced=False):
        """
        The semidirect product of ``self`` with ``H`` via ``hom``.

        If there exists a homomorphism `\phi` from a group `G` to the
        automorphism group of a group `H`, then we can define the semidirect
        product of `G` with `H` via `\phi` as the cartesian product of `G`
        and `H` with the operation

        .. MATH::

                (g_1, h_1)(g_2, h_2) = (g_1 g_2, \phi(g_2)(h_1) h_2).

        INPUT:

        - ``H`` -- Finitely presented group which is implicitly acted on
          by ``self`` and can be naturally embedded as a normal subgroup
          of the semidirect product.

        - ``hom`` -- Homomorphism from ``self`` to the automorphism group
          of ``H``. Given as a pair, with generators of ``self`` in the
          first slot and the images of the corresponding generators in the
          second. These images must be automorphisms of ``H``, given again
          as a pair of generators and images.

        - ``check`` -- Boolean (default ``True``). If ``False`` the defining
          homomorphism and automorphism images are not tested for validity.
          This test can be costly with large groups, so it can be bypassed
          if the user is confident that his morphisms are valid.

        - ``reduced`` -- Boolean (default ``False``). If ``True`` then the
          method attempts to reduce the presentation of the output group.

        OUTPUT:

        The semidirect product of ``self`` with ``H`` via ``hom`` as a
        finitely presented group. See
        :meth:`PermutationGroup_generic.semidirect_product
        <sage.groups.perm_gps.permgroup.PermutationGroup_generic.semidirect_product>`
        for a more in depth explanation of a semidirect product.

        AUTHORS:

        - Davis Shurbert (8-1-2013)

        EXAMPLES:

        Group of order 12 as two isomorphic semidirect products::

            sage: D4 = groups.presentation.Dihedral(4)
            sage: C3 = groups.presentation.Cyclic(3)
            sage: alpha1 = ([C3.gen(0)],[C3.gen(0)])
            sage: alpha2 = ([C3.gen(0)],[C3([1,1])])
            sage: S1 = D4.semidirect_product(C3, ([D4.gen(1), D4.gen(0)],[alpha1,alpha2]))
            sage: C2 = groups.presentation.Cyclic(2)
            sage: Q = groups.presentation.DiCyclic(3)
            sage: a = Q([1]); b = Q([-2])
            sage: alpha = (Q.gens(), [a,b])
            sage: S2 = C2.semidirect_product(Q, ([C2.0],[alpha]))
            sage: S1.is_isomorphic(S2)
            True

        Dihedral groups can be constructed as semidirect products
        of cyclic groups::

            sage: C2 = groups.presentation.Cyclic(2)
            sage: C8 = groups.presentation.Cyclic(8)
            sage: hom = (C2.gens(), [ ([C8([1])], [C8([-1])]) ])
            sage: D = C2.semidirect_product(C8, hom)
            sage: D.as_permutation_group().is_isomorphic(DihedralGroup(8))
            True

        You can attempt to reduce the presentation of the output group::

            sage: D = C2.semidirect_product(C8, hom); D
            Finitely presented group < a, b, c, d |
             a^2, b^-1*a^-1*b*a*d^-1*c^-1, c^-1*a^-1*c*a*d^-1, d^-1*a^-1*d*a,
             b^2*c^-1, c^-1*b^-1*c*b, d^-1*b^-1*d*b, c^2*d^-1, d^-1*c^-1*d*c, d^2 >
            sage: D = C2.semidirect_product(C8, hom, reduced=True); D
            Finitely presented group < a, b | a^2, (a*b)^2, b^8 >

            sage: C3 = groups.presentation.Cyclic(3)
            sage: C4 = groups.presentation.Cyclic(4)
            sage: hom = (C3.gens(), [(C4.gens(), C4.gens())])
            sage: C3.semidirect_product(C4, hom)
            Finitely presented group < a, b, c |
             a^3, b^-1*a^-1*b*a, c^-1*a^-1*c*a, b^2*c^-1, c^-1*b^-1*c*b, c^2 >
            sage: D = C3.semidirect_product(C4, hom, reduced=True); D
            Finitely presented group < a, b | a^3, b^4, b^-1*a^-1*b*a >
            sage: D.as_permutation_group().is_cyclic()
            True

        You can turn off the checks for the validity of the input morphisms.
        This check is expensive but behavior is unpredictable if inputs are
        invalid and are not caught by these tests::

            sage: C5 = groups.presentation.Cyclic(5)
            sage: C12 = groups.presentation.Cyclic(12)
            sage: hom = (C5.gens(), [(C12.gens(), C12.gens())])
            sage: sp = C5.semidirect_product(C12, hom, check=False); sp
            Finitely presented group < a, b, c, d |
             a^5, b^-1*a^-1*b*a, c^-1*a^-1*c*a, d^-1*a^-1*d*a, b^2*d^-1,
             c^-1*b^-1*c*b, d^-1*b^-1*d*b, c^3, d^-1*c^-1*d*c, d^2 >
            sage: sp.as_permutation_group().is_cyclic(), sp.order()
            (True, 60)

        TESTS:

        The following was fixed in Gap-4.7.2::

            sage: C5.semidirect_product(C12, hom) == sp
            True

        A more complicated semidirect product::

            sage: C = groups.presentation.Cyclic(7)
            sage: D = groups.presentation.Dihedral(5)
            sage: id1 = ([C.0], [(D.gens(),D.gens())])
            sage: Se1 =  C.semidirect_product(D, id1)
            sage: id2 = (D.gens(), [(C.gens(),C.gens()),(C.gens(),C.gens())])
            sage: Se2 =  D.semidirect_product(C ,id2)
            sage: Dp1 = C.direct_product(D);
            sage: Dp1.is_isomorphic(Se1), Dp1.is_isomorphic(Se2)
            (True, True)

        Most checks for validity of input are left to GAP to handle::

            sage: bad_aut = ([C.0], [(D.gens(),[D.0, D.0])])
            sage: C.semidirect_product(D, bad_aut)
            Traceback (most recent call last):
            ...
            ValueError: images of input homomorphism must be automorphisms
            sage: bad_hom = ([D.0, D.1], [(C.gens(),C.gens())])
            sage: D.semidirect_product(C, bad_hom)
            Traceback (most recent call last):
            ...
            ValueError: libGAP: Error, <gens> and <imgs> must be lists of same length
        """
        from sage.groups.free_group import FreeGroup, _lexi_gen

        if not isinstance(H, FinitelyPresentedGroup):
            raise TypeError("input must be a finitely presented group")

        GAP_self = self.gap(); GAP_H = H.gap()
        auto_grp = libgap.AutomorphismGroup(H.gap())
        self_gens = [h.gap() for h in hom[0]]
        # construct image automorphisms in GAP
        GAP_aut_imgs = [ libgap.GroupHomomorphismByImages(GAP_H, GAP_H, [g.gap() for g in gns],
            [i.gap() for i in img]) for (gns, img) in hom[1] ]

        # check for automorphism validity in images of operation defining homomorphism,
        # and construct the defining homomorphism.
        if check:
            if not all([a in libgap.List(libgap.AutomorphismGroup(GAP_H)) for a in GAP_aut_imgs]):
                raise ValueError("images of input homomorphism must be automorphisms")
            GAP_def_hom = libgap.GroupHomomorphismByImages(GAP_self, auto_grp, self_gens, GAP_aut_imgs)
        else:
            GAP_def_hom = GAP_self.GroupHomomorphismByImagesNC( auto_grp, self_gens, GAP_aut_imgs)

        prod = libgap.SemidirectProduct(GAP_self, GAP_def_hom, GAP_H)
        # Convert pc group to fp group
        if prod.IsPcGroup():
            prod = libgap.Image(libgap.IsomorphismFpGroupByPcgs(prod.FamilyPcgs() , 'x'))
        if not prod.IsFpGroup():
            raise NotImplementedError("unable to convert GAP output to equivalent Sage fp group")

        # Convert GAP group object to Sage via Tietze
        # lists for readability of variable names
        GAP_gens = prod.FreeGeneratorsOfFpGroup()
        name_itr = _lexi_gen() # Python generator for lexicographical variable names
        ret_F = FreeGroup([name_itr.next() for i in GAP_gens])
        ret_rls = tuple([ret_F(rel_word.TietzeWordAbstractWord(GAP_gens).sage())
            for rel_word in prod.RelatorsOfFpGroup()])
        ret_fpg = FinitelyPresentedGroup(ret_F, ret_rls)
        if reduced:
            ret_fpg = ret_fpg.simplified()
        return ret_fpg

    def _element_constructor_(self, *args, **kwds):
        """
        Construct an element of ``self``.

        TESTS::

            sage: G.<a,b> = FreeGroup()
            sage: H = G / (G([1]), G([2, 2, 2]))
            sage: H([1, 2, 1, -1]) # indirect doctest
            a*b
            sage: H([1, 2, 1, -2]) # indirect doctest
            a*b*a*b^-1
        """
        if len(args)!=1:
            return self.element_class(self, *args, **kwds)
        x = args[0]
        if x==1:
            return self.one()
        try:
            P = x.parent()
        except AttributeError:
            return self.element_class(self, x, **kwds)
        if P is self._free_group:
            return self.element_class(self, x.Tietze(), **kwds)
        return self.element_class(self, x, **kwds)

    @cached_method
    def abelian_invariants(self):
        r"""
        Return the abelian invariants of ``self``.

        The abelian invariants are given by a list of integers
        `(i_1, \ldots, i_j)`, such that the abelianization of the group is
        isomorphic to `\ZZ / (i_1) \times \cdots \times \ZZ / (i_j)`.

        EXAMPLES::

            sage: G = FreeGroup(4, 'g')
            sage: G.inject_variables()
            Defining g0, g1, g2, g3
            sage: H = G.quotient([g1^2, g2*g1*g2^(-1)*g1^(-1), g1*g3^(-2), g0^4])
            sage: H.abelian_invariants()
            (0, 4, 4)

        ALGORITHM:

        Uses GAP.
        """
        invariants = self.gap().AbelianInvariants()
        return tuple( i.sage() for i in invariants )

    def simplification_isomorphism(self):
        """
        Return an isomorphism from ``self`` to a finitely presented group with
        a (hopefully) simpler presentation.

        EXAMPLES::

            sage: G.<a,b,c> = FreeGroup()
            sage: H = G / [a*b*c, a*b^2, c*b/c^2]
            sage: I = H.simplification_isomorphism()
            sage: I
            Generic morphism:
              From: Finitely presented group < a, b, c | a*b*c, a*b^2, c*b*c^-2 >
              To:   Finitely presented group < b |  >
            sage: I(a)
            b^-2
            sage: I(b)
            b
            sage: I(c)
            b

        TESTS::

            sage: F = FreeGroup(1)
            sage: G = F.quotient([F.0])
            sage: G.simplification_isomorphism()
            Generic morphism:
              From: Finitely presented group < x | x >
              To:   Finitely presented group <  |  >

        ALGORITM:

        Uses GAP.
        """
        I = self.gap().IsomorphismSimplifiedFpGroup()
        domain = self
        codomain = wrap_FpGroup(I.Range())
        phi = lambda x: codomain(I.ImageElm(x.gap()))
        return self.hom(phi, codomain)

    def simplified(self):
        """
        Return an isomorphic group with a (hopefully) simpler presentation.

        OUTPUT:

        A new finitely presented group. Use
        :meth:`simplification_isomorphism` if you want to know the
        isomorphism.

        EXAMPLES::

            sage: G.<x,y> = FreeGroup()
            sage: H = G /  [x ^5, y ^4, y*x*y^3*x ^3]
            sage: H
            Finitely presented group < x, y | x^5, y^4, y*x*y^3*x^3 >
            sage: H.simplified()
            Finitely presented group < x, y | y^4, y*x*y^-1*x^-2, x^5 >

        A more complicate example::

            sage: G.<e0, e1, e2, e3, e4, e5, e6, e7, e8, e9> = FreeGroup()
            sage: rels = [e6, e5, e3, e9, e4*e7^-1*e6, e9*e7^-1*e0,
            ...           e0*e1^-1*e2, e5*e1^-1*e8, e4*e3^-1*e8, e2]
            sage: H = G.quotient(rels);  H
            Finitely presented group < e0, e1, e2, e3, e4, e5, e6, e7, e8, e9 |
            e6, e5, e3, e9, e4*e7^-1*e6, e9*e7^-1*e0, e0*e1^-1*e2, e5*e1^-1*e8, e4*e3^-1*e8, e2 >
            sage: H.simplified()
            Finitely presented group < e0 | e0^2 >
        """
        return self.simplification_isomorphism().codomain()

    def alexander_matrix(self):
        """
        Return the Alexander matrix of the group.

        This matrix is given by the fox derivatives of the relations
        with respect to the generators.

        OUTPUT:

        A group algebra-valued matrix. It depends on the (fixed)
        choice of presentation.

        EXAMPLES::

            sage: G.<a,b,c> = FreeGroup()
            sage: H = G.quotient([a*b/a/b, a*c/a/c, c*b/c/b])
            sage: H.alexander_matrix()
            [     B[1] - B[a*b*a^-1] B[a] - B[a*b*a^-1*b^-1]                       0]
            [     B[1] - B[a*c*a^-1]                       0 B[a] - B[a*c*a^-1*c^-1]]
            [                      0 B[c] - B[c*b*c^-1*b^-1]      B[1] - B[c*b*c^-1]]

            sage: G.<a,b,c,d,e> = FreeGroup()
            sage: H = G.quotient([a*b/a/b, a*c/a/c, a*d/a/d, b*c*d/(c*d*b), b*c*d/(d*b*c)])
            sage: H.alexander_matrix()
            [              B[1] - B[a*b*a^-1]          B[a] - B[a*b*a^-1*b^-1]                                0                                0                                0]
            [              B[1] - B[a*c*a^-1]                                0          B[a] - B[a*c*a^-1*c^-1]                                0                                0]
            [              B[1] - B[a*d*a^-1]                                0                                0          B[a] - B[a*d*a^-1*d^-1]                                0]
            [                               0             B[1] - B[b*c*d*b^-1]   B[b] - B[b*c*d*b^-1*d^-1*c^-1]      B[b*c] - B[b*c*d*b^-1*d^-1]                                0]
            [                               0        B[1] - B[b*c*d*c^-1*b^-1]             B[b] - B[b*c*d*c^-1] B[b*c] - B[b*c*d*c^-1*b^-1*d^-1]                                0]
        """
        rel = self.relations()
        gen = self._free_group.gens()
        return matrix(len(rel), len(gen),
                      lambda i,j: rel[i].fox_derivative(gen[j]))

    def rewriting_system(self):
        """
        Return the rewriting system corresponding to the finitely presented
        group. This rewriting system can be used to reduce words with respect
        to the relations.

        If the rewriting system is transformed into a confluent one, the
        reduction process will give as a result the (unique) reduced form
        of an element.

        EXAMPLES::

            sage: F.<a,b> = FreeGroup()
            sage: G = F / [a^2,b^3,(a*b/a)^3,b*a*b*a]
            sage: k = G.rewriting_system()
            sage: k
            Rewriting system of Finitely presented group < a, b | a^2, b^3, a*b^3*a^-1, (b*a)^2 >
            with rules:
                a^2    --->    1
                b^3    --->    1
                (b*a)^2    --->    1
                a*b^3*a^-1    --->    1

            sage: G([1,1,2,2,2])
            a^2*b^3
            sage: k.reduce(G([1,1,2,2,2]))
            1
            sage: k.reduce(G([2,2,1]))
            b^2*a
            sage: k.make_confluent()
            sage: k.reduce(G([2,2,1]))
            a*b
        """
        return RewritingSystem(self)

