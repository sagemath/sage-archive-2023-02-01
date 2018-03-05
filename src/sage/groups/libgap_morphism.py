r"""
Group homomorphisms for ``LibGAPGroup``.

EXAMPLES::


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
from sage.groups.libgap_wrapper import ParentLibGAP
from sage.categories.groups import Groups
from sage.categories.morphism import Morphism
from sage.misc.latex import latex


class GroupMorphism_libgap(Morphism):
    r"""
    Group morphism specified by the images of generators.

    This wraps GAP's GroupHomomorphismByImages function.
    It can be expensive if the group is large.

    Input:

    - ``homset`` -- the parent
    - ``imgs`` -- a tuple of generators
    - ``check`` -- (default: ``True``) check if the images define a group homomorphism.

    EXAMPLES::

        sage: from sage.groups.abelian_gps.abelian_group_gap import AbelianGroupGap
        sage: A = AbelianGroupGap([2,4])
        sage: A.hom([g^2 for g in A.gens()])
        Group endomorphism of Abelian group with gap, generator orders (2, 4)

    Homomorphisms can be defined between different kinds of libGAP groups::

        sage: G = MatrixGroup([Matrix(ZZ,2,[0,1,1,0])])
        sage: f = A.hom([G.0,G(1)])
        sage: f
        Group morphism:
        From: Abelian group with gap, generator orders (2, 4)
        To:   Matrix group over Integer Ring with 1 generators (
        [0 1]
        [1 0]
        )
    """
    def __init__(self, homset, imgs, check=True):
        r"""
        Constructor method

        TESTS:

        The following input does not define a valid group homomorphism

            sage: from sage.groups.abelian_gps.abelian_group_gap import AbelianGroupGap
            sage: A = AbelianGroupGap([2])
            sage: G = MatrixGroup([Matrix(ZZ,2,[0,1,-1,0])])
            sage: A.hom([G.0])
            Traceback (most recent call last):
            ...
            ValueError: the map (f1,)-->[[ [ 0, 1 ], [ -1, 0 ] ]] is not a homomorphism
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
                raise ValueError('the map {}-->{} is not a homomorphism'.format(dom.gens(), imgs))
        else:
            self._phi = libgap.GroupHomomorphismByImagesNC(dom.gap(), codom.gap(), gens, imgs)

    def _repr_type(self):
        r"""
        Part of the implementation of :meth:`_repr_`.

        EXAMPLES::

            sage: from sage.groups.abelian_gps.abelian_group_gap import AbelianGroupGap
            sage: A = AbelianGroupGap([2,4])
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
            sage: A1 = AbelianGroupGap([6,6])
            sage: A2 = AbelianGroupGap([3,3])
            sage: f = A1.hom(A2.gens())
            sage: f.kernel()
            Subgroup of Abelian group with gap, generator orders (6, 6) generated by (f1*f2, f3*f4)
            sage: f.kernel().order()
            4
        """
        dom = self.domain()
        ker_gen = [dom(g) for g in self.gap().Kernel().GeneratorsOfGroup()]
        return dom.subgroup(ker_gen)

    def pushforward(self, J, *args,**kwds):
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



        """
        dom = self.domain()
        codom = self.codomain()
        phi = self.gap()
        if J in self.domain():
            return self(J)
        if not isinstance(J, ParentLibGAP):
            raise TypeError("J(=%s) must be a libgap group" %J)
        if self.gap().IsSubgroup(J.gap()).sage():
                im_gens = phi.Image(gapJ).GeneratorsOfGroup()
                im_gens = [self.codomain()(g) for g in im_gens]
                return codom.subgroup(im_gens)

    image = pushforward

    def _call_(self, g):
        """
        Call syntax for morphisms.

        Some python code for wrapping GAP's ``Images`` function.
        Returns an error if ``g`` is not in ``G``.

        EXAMPLES::

            sage: from sage.groups.abelian_gps.abelian_group_gap import AbelianGroupGap
            sage: A = AbelianGroupGap([2,4])
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
            sage: f.lift(a)
            Traceback (most recent call last):
            ...
            ValueError: f2 is not an element of the image of Group endomorphism of
            Abelian group with gap, generator orders (2, 4)
        """
        if not h in self.codomain():
            raise ValueError("h (=%s) must be an element of the codomain")
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
        gens = phi.PreImage(S.gap()).GeneratorsOfGroup()
        return self.codomain().subgroup(gens)

class GroupHomset_libgap(GroupHomset_generic):
    r"""

    """
    def __init__(self, G, H, category=None):
        r"""
        Return the homset of two libgap groups.

        INPUT:

        - ``G`` -- a libgap group

        - ``H`` -- a libgap group

        OUTPUT:

        The homset of two libgap groups.

        EXAMPLES::

        """
        if not isinstance(G, ParentLibGAP):
            raise TypeError("G (=%s) must be a ParentLibGAP group." %G)
        if not isinstance(H, ParentLibGAP):
            raise TypeError("H (=%s) must be a ParentLibGAP group." %G)
        category = Groups().or_subcategory(category)
        if G.is_finite() and H.is_finite():
            category = category.Finite()
        from sage.categories.homset import Homset
        Homset.__init__(self, G, H, category)

    def __call__(self, im_gens, check=True):
        """
        Return the homomorphism defined by the images of the generators.

        INPUT:

        - ``im_gens`` -- iterable, the list of images of the
          generators of the domain

        - ``check`` -- bool (optional, default: ``True``), whether to
          check if images define a valid homomorphism

        OUTPUT:

        A Group homomorphism.

        EXAMPLES::


        """
        return GroupMorphism_libgap(self, im_gens, check=check)

    Element = GroupMorphism_libgap
