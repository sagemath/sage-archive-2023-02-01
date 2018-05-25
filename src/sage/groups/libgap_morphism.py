r"""
Group homomorphisms for groups with a libGAP backend.

EXAMPLES::

    sage: from sage.groups.abelian_gps.abelian_group_gap import AbelianGroupGap
    sage: A = AbelianGroupGap([2, 4])
    sage: F.<a,b> = FreeGroup()
    sage: f = F.hom([g for g in A.gens()])
    sage: K = f.kernel()
    sage: K
    Group(<free, no generators known>)

AUTHORS:

- Simon Brandhorst (2018-02-08): initial version
"""

# ****************************************************************************
#       Copyright (C) 2018 Simon Brandhorst <sbrandhorst@web.de>
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 2 of the License, or
# (at your option) any later version.
#                  http://www.gnu.org/licenses/
# ****************************************************************************
from __future__ import print_function
from sage.groups.group_homset import GroupHomset_generic
from sage.categories.all import HomsetWithBase, Groups
from sage.rings.integer_ring import ZZ
from sage.groups.libgap_wrapper import ParentLibGAP
from sage.categories.groups import Groups
from sage.categories.morphism import Morphism
from sage.misc.latex import latex


class GroupMorphism_libgap(Morphism):
    r"""
    Group morphism specified by the images of generators.

    This wraps GAP's GroupHomomorphismByImages function.
    Checking if the input defines a group homomorphism can be expensive
    if the group is large.

    INPUT:

    - ``homset`` -- the parent
    - ``imgs`` -- a tuple of generators
    - ``check`` -- (default: ``True``) check if the images define
      a group homomorphism.

    EXAMPLES::

        sage: from sage.groups.abelian_gps.abelian_group_gap import AbelianGroupGap
        sage: A = AbelianGroupGap([2, 4])
        sage: A.hom([g^2 for g in A.gens()])
        Group endomorphism of Abelian group with gap, generator orders (2, 4)

    Homomorphisms can be defined between different kinds of libGAP groups::

        sage: G = MatrixGroup([Matrix(ZZ, 2, [0,1,1,0])])
        sage: f = A.hom([G.0, G(1)])
        sage: f
        Group morphism:
        From: Abelian group with gap, generator orders (2, 4)
        To:   Matrix group over Integer Ring with 1 generators (
        [0 1]
        [1 0]
        )
        sage: G.<a,b> = FreeGroup()
        sage: H = G / (G([1]), G([2])^3)
        sage: f = G.hom(H.gens())
        sage: f
        Group morphism:
          From: Free Group on generators {a, b}
          To:   Finitely presented group < a, b | a, b^3 >

    """

    def __init__(self, homset, imgs, check=True):
        r"""
        Constructor method.

        TESTS:

        The following input does not define a valid group homomorphism::

            sage: from sage.groups.abelian_gps.abelian_group_gap import AbelianGroupGap
            sage: A = AbelianGroupGap([2])
            sage: G = MatrixGroup([Matrix(ZZ, 2, [0,1,-1,0])])
            sage: A.hom([G.0])
            Traceback (most recent call last):
            ...
            ValueError: images do not define a group homomorphism
        """
        from sage.libs.gap.libgap import libgap
        Morphism.__init__(self, homset)
        dom = homset.domain()
        codom = homset.codomain()
        gens = [x.gap() for x in dom.gens()]
        imgs = [codom(x).gap() for x in imgs]
        if check:
            if not len(gens) == len(imgs):
                raise ValueError("provide an image for each generator")
            self._phi = libgap.GroupHomomorphismByImages(dom.gap(), codom.gap(), gens, imgs)
            if not self._phi.IsGroupHomomorphism():
                raise ValueError("images do not define a group homomorphism")
        else:
            self._phi = libgap.GroupHomomorphismByImagesNC(dom.gap(), codom.gap(), gens, imgs)

    def __reduce__(self):
        r"""
        Implements pickling.

        We have to work around the fact that gap does not provide pickling.

        EXAMPLES::

            sage: from sage.groups.abelian_gps.abelian_group_gap import AbelianGroupGap
            sage: G = AbelianGroupGap([3,2,5])
            sage: f = G.hom(G, G.gens())
            sage: f == loads(dumps(f))
            True
        """
        return self.parent(), (tuple(self(g) for g in self.domain().gens()),)

    def _repr_type(self):
        r"""
        Part of the implementation of :meth:`_repr_`.

        EXAMPLES::

            sage: from sage.groups.abelian_gps.abelian_group_gap import AbelianGroupGap
            sage: A = AbelianGroupGap([2, 4])
            sage: A.hom([g^2 for g in A.gens()])
            Group endomorphism of Abelian group with gap, generator orders (2, 4)
        """
        return "Group"

    def gap(self):
        r"""
        Return the underlying LibGAP group homomorphism.

        EXAMPLES::

            sage: from sage.groups.abelian_gps.abelian_group_gap import AbelianGroupGap
            sage: A = AbelianGroupGap([2,4])
            sage: f = A.hom([g^2 for g in A.gens()])
            sage: f.gap()
            [ f1, f2 ] -> [ <identity> of ..., f3 ]
        """
        return self._phi

    def _latex_(self):
        r"""
        Return the latex representation of ``self``.

        EXAMPLES::

            sage: from sage.groups.abelian_gps.abelian_group_gap import AbelianGroupGap
            sage: A = AbelianGroupGap([2,4])
            sage: f = A.hom([g^2 for g in A.gens()])
            sage: f._latex_()
            'Abelian group with gap, generator orders $(2, 4)$ \\rightarrow{} Abelian group with gap, generator orders $(2, 4)$'
        """
        return "%s \\rightarrow{} %s"%(latex(self.domain()), latex(self.codomain()))

    def kernel(self):
        r"""
        Return the kernel of ``self``.

        EXAMPLES::

            sage: from sage.groups.abelian_gps.abelian_group_gap import AbelianGroupGap
            sage: A1 = AbelianGroupGap([6, 6])
            sage: A2 = AbelianGroupGap([3, 3])
            sage: f = A1.hom(A2.gens())
            sage: f.kernel()
            Subgroup of Abelian group with gap, generator orders (6, 6) generated by (f1*f2, f3*f4)
            sage: f.kernel().order()
            4
        """
        dom = self.domain()
        return dom._subgroup_constructor(self.gap().Kernel())

    def pushforward(self, J, *args, **kwds):
        r"""
        The image of an element or a subgroup.

        INPUT:

        ``J`` -- a subgroup or an element of the domain of ``self``

        OUTPUT:

        The image of ``J`` under ``self``.

        .. NOTE::

            ``pushforward`` is the method that is used when a map is called
            on anything that is not an element of its domain. For historical
            reasons, we keep the alias ``image()`` for this method.

        EXAMPLES::

            sage: G.<a,b> = FreeGroup()
            sage: H = G / (G([1]), G([2])^3)
            sage: f = G.hom(H.gens())
            sage: S = G.subgroup([a.gap()])
            sage: f.pushforward(S)
            Group([ a ])
            sage: x = f.image(a)
            sage: x
            a
            sage: x.parent()
            Finitely presented group < a, b | a, b^3 >
        """
        dom = self.domain()
        codom = self.codomain()
        phi = self.gap()
        if isinstance(J, dom.Element) and J in dom:
            return self(J)
        if not isinstance(J, ParentLibGAP):
            raise TypeError("J(=%s) must be a libgap group" %J)
        if dom.gap().IsSubgroup(J.gap()).sage():
            return codom._subgroup_constructor(phi.Image(J.gap()))

    image = pushforward

    def _call_(self, g):
        """
        Call syntax for morphisms.

        Some python code for wrapping GAP's ``Images`` function.

        INPUT:

        - ``g`` -- an element of the domain

        EXAMPLES::

            sage: from sage.groups.abelian_gps.abelian_group_gap import AbelianGroupGap
            sage: A = AbelianGroupGap([2, 4])
            sage: f = A.hom([g^2 for g in A.gens()])
            sage: a = A.gens()[1]
            sage: f(a)
            f3
        """
        return self.codomain()(self.gap().Image(g.gap()))

    def lift(self, h):
        r"""
        Return an element of the domain that maps to ``h``.

        EXAMPLES::

            sage: from sage.groups.abelian_gps.abelian_group_gap import AbelianGroupGap
            sage: A = AbelianGroupGap([2,4])
            sage: f = A.hom([g^2 for g in A.gens()])
            sage: a = A.gens()[1]
            sage: f.lift(a^2)
            f2

        If the element is not in the image, we raise an error::

            sage: f.lift(a)
            Traceback (most recent call last):
            ...
            ValueError: f2 is not an element of the image of Group endomorphism of
            Abelian group with gap, generator orders (2, 4)
        """
        if not h in self.codomain():
            raise TypeError("h (=%s) must be an element of the codomain" %h)
        h = self.codomain()(h)
        phi = self.gap()
        if h.gap() in phi.Image():
            return self.domain()(phi.PreImagesRepresentative(h.gap()))
        else:
            raise ValueError("%s is not an element of the image of %s" %(h, self))

    def preimage(self, S):
        r"""
        Return the preimage of the subgroup ``S``.

        INPUT:

        - ``S`` a subgroup of this group.

        EXAMPLES::

            sage: from sage.groups.abelian_gps.abelian_group_gap import AbelianGroupGap
            sage: A = AbelianGroupGap([2,4])
            sage: f = A.hom([g^2 for g in A.gens()])
            sage: S = A.subgroup(A.gens()[1:])
            sage: f.preimage(S)
            Subgroup of Abelian group with gap, generator orders (2, 4) generated by (f1, f3, 1)
        """
        phi = self.gap()
        if not isinstance(S, ParentLibGAP):
            raise TypeError("%s must be a libGAP group of %s"%(S, self))
        if not self.codomain().gap().IsSubgroup(S.gap()).sage():
            raise ValueError("%s must be a subgroup of %s"%(S, self))
        preimage = phi.PreImage(S.gap())
        return self.codomain()._subgroup_constructor(preimage)

class GroupHomset_libgap(HomsetWithBase):
    r"""
    Homsets of groups with a libgap backend.

    Do not call this directly instead use :meth:`Hom`.

    INPUT:

    - ``G`` -- a libgap group
    - ``H`` -- a libgap group
    - ``category`` -- a category

    OUTPUT:

    The homset of two libgap groups.

    EXAMPLES::

            sage: from sage.groups.abelian_gps.abelian_group_gap import AbelianGroupGap
            sage: A = AbelianGroupGap([2,4])
            sage: H = A.Hom(A)
            sage: H
            Set of Morphisms from Abelian group with gap, generator orders (2, 4)
             to Abelian group with gap, generator orders (2, 4)
             in Category of finite enumerated commutative groups
    """

    def __init__(self, G, H, category=None):
        r"""
        Return the homset of two libgap groups.

        TESTS::

            sage: from sage.groups.abelian_gps.abelian_group_gap import AbelianGroupGap
            sage: A = AbelianGroupGap([2,4])
            sage: H = A.Hom(A)
            sage: TestSuite(H).run(skip="_test_elements")

        """
        if not isinstance(G, ParentLibGAP):
            raise TypeError("G (=%s) must be a ParentLibGAP group." %G)
        if not isinstance(H, ParentLibGAP):
            raise TypeError("H (=%s) must be a ParentLibGAP group." %G)
        category = Groups() & category
        try:
            if G.is_finite() and H.is_finite():
                category = category.Finite()
        except NotImplementedError:
            pass
        from sage.categories.homset import Homset
        HomsetWithBase.__init__(self, G, H, category, ZZ)

    Element = GroupMorphism_libgap

    def _element_constructor_(self, x, check=True):
        r"""
        Handle conversions and coercions.

        EXAMPLES::

            sage: from sage.groups.abelian_gps.abelian_group_gap import AbelianGroupGap
            sage: A = AbelianGroupGap([2,4])
            sage: H = A.Hom(A)
            sage: H([g^2 for g in A.gens()])
            Group endomorphism of Abelian group with gap, generator orders (2, 4)
        """
        if type(x) in [tuple, list]:
            im_gens = x
            im_gens = tuple([self.codomain()(g) for g in im_gens])
            return GroupMorphism_libgap(self, im_gens, check=check)

    def _an_element_(self):
        r"""
        Return an element of ``self``.

        There is only one group homomorphism which is allways defined -
        the trivial one.

        EXAMPLES::

            sage: from sage.groups.abelian_gps.abelian_group_gap import AbelianGroupGap
            sage: A = AbelianGroupGap([2,4])
            sage: G.<a,b> = FreeGroup()
            sage: A.Hom(G).an_element()
            Group morphism:
              From: Abelian group with gap, generator orders (2, 4)
              To:   Free Group on generators {a, b}
        """
        n = len(self.domain().gens())
        one = self.codomain().one()
        return self(n * [one])
