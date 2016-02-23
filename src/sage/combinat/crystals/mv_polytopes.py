# -*- coding: utf-8 -*-
r"""
Crystal Of Mirković-Vilonen (MV) Polytopes
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

from sage.structure.parent import Parent
from sage.structure.unique_representation import UniqueRepresentation
from sage.categories.highest_weight_crystals import HighestWeightCrystals
from sage.categories.infinite_enumerated_sets import InfiniteEnumeratedSets
from sage.combinat.root_system.cartan_type import CartanType
from sage.combinat.crystals.pbw_crystal import PBWCrystalElement, PBWCrystal

class MVPolytope(PBWCrystalElement):
    """
    A Mirković-Vilonen (MV) polytope.
    """
    def _repr_(self):
        """
        Return a string representation of ``self``.
        """
        pbw_datum = self._pbw_datum.convert_to_new_long_word(self.parent()._default_word)
        return "MV polytope with Lusztig datum {}".format(pbw_datum.lusztig_datum)

    def _latex_(self):
        """
        Return a latex representation of ``self``.
        """
        from sage.misc.latex import latex
        return latex(self.polytope())

    def _polytope_vertices(self, P):
        """
        Return a list of the vertices of ``self`` in ``P``.
        """
        pbw_data = self._pbw_datum.parent
        W = pbw_data.weyl_group
        w0 = W.long_element()
        al = P.simple_roots()

        vertices = set([P.zero()])
        for red in w0.reduced_words():
            cur = P.zero()
            red = tuple(red)
            roots = [P.sum(c*al[a] for a,c in root)
                     for root in pbw_data._root_list_from(red)]
            datum = pbw_data.convert_to_new_long_word(self._pbw_datum, red)
            for i,c in enumerate(datum.lusztig_datum):
                cur = cur + roots[i] * c
                vertices.add(cur)
        return list(vertices)

    def polytope(self, P=None):
        """
        Return a polytope of ``self``.

        INPUT:

        - ``P`` -- (optional) a space to realize the polytope; default is
          the weight lattice realization of the crystal
        """
        if P is None:
            P = self.parent().weight_lattice_realization()

        from sage.geometry.polyhedron.constructor import Polyhedron
        return Polyhedron([v.to_vector() for v in self._polytope_vertices(P)])

    def plot(self, P=None, **options):
        """
        Plot ``self``.

        INPUT:

        - ``P`` -- (optional) a space to realize the polytope; default is
          the weight lattice realization of the crystal

        .. SEEALSO::

            :meth:`~sage.combiant.root_system.root_lattice_realizations.RootLatticeRealizations.ParentMethods.plot_mv_polytope`
        """
        if P is None:
            P = self.parent().weight_lattice_realization()
        return P.plot_mv_polytope(self, **options)

class MVPolytopes(PBWCrystal):
    def _repr_(self):
        """
        Return a string representation of ``self``.
        """
        return "MV polytopes of type {}".format(self._cartan_type)

    Element = MVPolytope

