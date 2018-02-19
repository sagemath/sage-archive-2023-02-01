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
from sage.categories.morphism import Morphism
from sage.misc.latex import latex


class LibGAPGroupMorphism(Morphism):
    r"""
    Group morphism specified by images of generators.

    Some Python code for wrapping GAP's GroupHomomorphismByImages
    function. Can be expensive if G is large.
    """
    def __init__(self, homset, imgs, check=True):
        r"""

        EXAMPLES::


        """
        from sage.libs.gap.libgap import libgap
        Morphism.__init__(self, homset)
        dom = homset.domain()
        codom = homset.codomain()
        gens = [x.gap() for x in dom.gens()]
        imgs = [H(x).gap() for x in imgs]
        if check:
            if not len(gens) == len(imgs):
                raise ValueError("provide an image for each generator")
            self._phi = libgap.GroupHomomorphismByImages(G.gap(), H.gap(), gens, imgs)
            if not self._phi.IsGroupHomomorphism():
                raise ValueError('the map {}-->{} is not a homomorphism'.format(G.gens(), imgsH))
        else:
            self._phi = libgap.GroupHomomorphismByImagesNC(G.gap(), H.gap(), gens, imgs)

    def _repr_type(self):
        r"""
        Part of the implementation of :meth:`_repr_`

        EXAMPLES::

        """
        return "GroupHomomorphism"

    def gap(self):
        r"""
        Return the underlying LibGAP group homomorphism

        OUTPUT:

        A LibGAP element.

        EXAMPLES::


        """
        return self._phi

    def _repr_(self):
        r"""
        EXAMPLES::


        """
        return "Homomorphism : %s --> %s"%(self.domain(),self.codomain())

    def _latex_(self):
        r"""
        EXAMPLES::

        """
        return "%s \\rightarrow{} %s"%(latex(self.domain()), latex(self.codomain()))

    def kernel(self):
        r"""
        Return the kernel of ``self``.

        EXAMPLES::

            sage: from sage.groups.abelian_gps.abelian_group_gap import AbelianGroupGap
            sage: A1 = AbelianGroupGap([6,6])
            sage: A2 = AbelianGroupGap([3,3])
            sage: f = A1.hom(A2)

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
        try:
            if J.is_subgroup(dom):
                im_gens = phi.Image(gapJ).GeneratorsOfGroup()
                im_gens = [self.codomain()(g) for g in im_gens]
                return codom.subgroup(im_gens)
        except AttributeError:
            pass
        # how hard do we want to try?
        try:
            gen = [H(g) for g in J.gens()]
            J = H.subgroup(gen)
            return self.pushforward(gen)
        except AttributeError:
            pass
        J = H(J)
        return self(J)

    image = pushforward

    def _call_(self, g):
        """
        Call syntax for morphisms.

        Some python code for wrapping GAP's ``Images`` function.
        Returns an error if ``g`` is not in ``G``.

        EXAMPLES::


        """
        return self.codomain()(self.gap().Image(g.gap()))

    def lift(self, h):
        r"""
        Return an element of the tomain that maps to ``h``.
        """
        phi = self.gap()
        return self.domain()(phi.PreImagesRepresentative())

    def preimage(self, S):
        r"""
        """
        phi = self.gap()
        try:
            gapS = S.gap()
        except AttributeError:
            gapS = self.codomain(S).gap()

        if self.domain().is_subgroup(S):
            gens = phi.PreImage(gapS)
            return self.codomain().subgroup(gens)
