r"""
Group homomorphisms gap



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

            sage: from sage.groups.matrix_gps.morphism import MatrixGroupMap
            sage: MatrixGroupMap(ZZ.Hom(ZZ))   # mathematical nonsense
            MatrixGroup endomorphism of Integer Ring
        """
        from sage.libs.gap.libgap import libgap
        Morphism.__init__(self, homset)
        G = homset.domain()
        H = homset.codomain()
        gens = [x.gap() for x in G.gens()]
        imgs = [H(x).gap() for x in imgs]
        if check:
            self._phi = libgap.GroupHomomorphismByImages(G.gap(), H.gap(), gens, imgs)
            if not self._phi.IsGroupHomomorphism():
                raise ValueError('the map {}-->{} is not a homomorphism'.format(G.gens(), imgsH))
        else:
            self._phi = libgap.GroupHomomorphismByImagesNC(G.gap(), H.gap(), gens, imgs)

    def _repr_type(self):
        r"""
        Part of the implementation of :meth:`_repr_`

        EXAMPLES::

            sage: from sage.groups.matrix_gps.morphism import MatrixGroupMap
            sage: MatrixGroupMap(ZZ.Hom(ZZ))._repr_type()
            'MatrixGroup'
        """
        return "GroupHomomorphism"

    def gap(self):
        r"""
        Return the underlying LibGAP group homomorphism

        OUTPUT:

        A LibGAP element.

        EXAMPLES::

            sage: F = GF(5); MS = MatrixSpace(F,2,2)
            sage: G = MatrixGroup([MS([1,1,0,1])])
            sage: H = MatrixGroup([MS([1,0,1,1])])
            sage: phi = G.hom(H.gens())
            sage: phi.gap()
            CompositionMapping( [ (6,7,8,10,9)(11,13,14,12,15)(16,19,20,18,17)(21,25,22,24,23) ]
            -> [ [ [ Z(5)^0, 0*Z(5) ], [ Z(5)^0, Z(5)^0 ] ] ], <action isomorphism> )
            sage: type(_)
            <type 'sage.libs.gap.element.GapElement'>
        """
        return self._phi

    def _repr_(self):
        r"""
        EXAMPLES::

            sage: F = GF(5); MS = MatrixSpace(F,2,2)
            sage: G = MatrixGroup([MS([1,1,0,1])])
            sage: H = MatrixGroup([MS([1,0,1,1])])
            sage: phi = G.hom(H.gens())
            sage: phi
            Homomorphism : Matrix group over Finite Field of size 5 with 1 generators (
            [1 1]
            [0 1]
            ) --> Matrix group over Finite Field of size 5 with 1 generators (
            [1 0]
            [1 1]
            )
            sage: phi(MS([1,1,0,1]))
            [1 0]
            [1 1]
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


        """
        gap_ker = self.gap().Kernel()
        G = self.domain()
        return G._subgroup_constructor(gap_ker)

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
        G = self.domain()
        H = self.codomain()
        phi = self.gap()
        C = self.codomain()
        if G.is_subgroup(H):
            return C._subgroup_constructor(phi.Image(gapJ).GeneratorsOfGroup())
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
        phi = self.gap()
        C = self.codomain()
        h = g.gap()
        return C(phi.Image(h))

    def lift(self, h):
        r"""
        """
        phi = self._phi
        return self.domain()(phi.PreImagesRepresentative())

    def preimage(self, S):
        r"""
        """
        phi = self._phi
        gapS = S.gap()
        if self.domain().is_subgroup(S):
            return self.codomain()._subgroup_constructor(phi.PreImage(gapS))
