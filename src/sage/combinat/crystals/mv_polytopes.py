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

    EXAMPLES:

    We can create an animation showing how the MV polytope changes
    under a string of crystal operators::

        sage: MV = crystals.infinity.MVPolytopes(['C', 2])
        sage: u = MV.highest_weight_vector()
        sage: L = RootSystem(['C',2,1]).ambient_space()
        sage: s = [1,2,1,2,2,2,1,1,1,1,2,1,2,2,1,2]
        sage: BB = [[-9, 2], [-10, 2]]
        sage: p = L.plot(reflection_hyperplanes=False, bounding_box=BB)  # long time
        sage: frames = [p + L.plot_mv_polytope(u.f_string(s[:i]),  # long time
        ....:                                  circle_size=0.1,
        ....:                                  wireframe='green',
        ....:                                  fill='purple',
        ....:                                  bounding_box=BB)
        ....:           for i in range(len(s))]
        sage: for f in frames:  # long time
        ....:     f.axes(False)
        sage: animate(frames).show(delay=60) # optional -- ImageMagick # long time
    """
    def _repr_(self):
        """
        Return a string representation of ``self``.

        EXAMPLES::

            sage: MV = crystals.infinity.MVPolytopes(['E', 6])
            sage: b = MV.module_generators[0].f_string([1,2,6,4,3,2,5,2])
            sage: b
            MV polytope with Lusztig datum (0, 1, ..., 1, 0, 0, 0, 0, 0, 0, 3, 1)
        """
        pbw_datum = self._pbw_datum.convert_to_new_long_word(self.parent()._default_word)
        return "MV polytope with Lusztig datum {}".format(pbw_datum.lusztig_datum)

    def _latex_(self):
        r"""
        Return a latex representation of ``self``.

        EXAMPLES::

            sage: MV = crystals.infinity.MVPolytopes(['C', 2])
            sage: b = MV.module_generators[0].f_string([1,2,1,2])
            sage: latex(b)
            \begin{tikzpicture}
            \draw (0, 0) -- (-1, 1) -- (-1, 1) -- (-2, 0) -- (-2, -2);
            \draw (0, 0) -- (0, -2) -- (-1, -3) -- (-1, -3) -- (-2, -2);
            \draw[fill=black] (0, 0) circle (0.1);
            \draw[fill=black] (-2, -2) circle (0.1);
            \end{tikzpicture}

        ::

            sage: MV = crystals.infinity.MVPolytopes(['D',4])
            sage: b = MV.module_generators[0].f_string([1,2,1,2])
            sage: latex(b)
            \text{\texttt{MV{ }polytope{ }...}}
        """
        latex_options = self.parent()._latex_options
        P = latex_options['P']
        plot_options = P.plot_parse_options(projection=latex_options["projection"])
        proj = plot_options.projection
        if proj(P.zero()).parent().dimension() != 2:
            from sage.misc.latex import latex
            return latex(repr(self))

        # We need this to use tikz
        from sage.graphs.graph_latex import setup_latex_preamble
        setup_latex_preamble()

        pbw_data = self._pbw_datum.parent
        W = pbw_data.weyl_group
        w0 = W.long_element()
        al = P.simple_roots()
        ret = "\\begin{tikzpicture}\n"

        final = None
        for red in w0.reduced_words():
            ret += "\\draw "
            cur = proj(P.zero())
            red = tuple(red)
            ret += str(cur)
            roots = [proj(P.sum(c*al[a] for a,c in root))
                     for root in pbw_data._root_list_from(red)]
            datum = pbw_data.convert_to_new_long_word(self._pbw_datum, red)
            for i in reversed(range(len(datum.lusztig_datum))):
                cur -= roots[i] * datum.lusztig_datum[i]
                ret += " -- " + str(cur)
            final = cur
            ret += ";\n"

        if latex_options["mark_endpoints"]:
            circle_size = latex_options["circle_size"]
            ret += "\\draw[fill=black] {} circle ({});\n".format(proj(P.zero()), circle_size)
            ret += "\\draw[fill=black] {} circle ({});\n".format(proj(final), circle_size)
        ret += "\\end{tikzpicture}"
        return ret

    def _polytope_vertices(self, P):
        """
        Return a list of the vertices of ``self`` in ``P``.

        EXAMPLES::

            sage: MV = crystals.infinity.MVPolytopes(['C', 3])
            sage: b = MV.module_generators[0].f_string([1,2,1,2])
            sage: sorted(b._polytope_vertices(MV.weight_lattice_realization()), key=list)
            [(0, 0, 0), (2, 0, -2), (0, 2, -2)]

            sage: MV = crystals.infinity.MVPolytopes(['D', 4])
            sage: b = MV.module_generators[0].f_string([1,2,3,4])
            sage: P = RootSystem(['D',4]).weight_lattice()
            sage: sorted(b._polytope_vertices(P), key=list)  # long time
            [0,
             -Lambda[1] + Lambda[3] + Lambda[4],
             Lambda[1] - Lambda[2] + Lambda[3] + Lambda[4],
             -2*Lambda[2] + 2*Lambda[3] + 2*Lambda[4],
             -Lambda[2] + 2*Lambda[3],
             -Lambda[2] + 2*Lambda[4]]
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

        EXAMPLES::

            sage: MV = crystals.infinity.MVPolytopes(['C', 3])
            sage: b = MV.module_generators[0].f_string([3,2,3,2,1])
            sage: P = b.polytope(); P
            A 3-dimensional polyhedron in QQ^3 defined as the convex hull of 6 vertices
            sage: P.vertices()
            (A vertex at (0, 0, 0),
             A vertex at (0, 1, -1),
             A vertex at (0, 1, 1),
             A vertex at (1, -1, 0),
             A vertex at (1, 1, -2),
             A vertex at (1, 1, 2))
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

            :meth:`~sage.combinat.root_system.root_lattice_realizations.RootLatticeRealizations.ParentMethods.plot_mv_polytope`

        EXAMPLES::

            sage: MV = crystals.infinity.MVPolytopes(['C', 2])
            sage: b = MV.highest_weight_vector().f_string([1,2,1,2,2,2,1,1,1,1,2,1])
            sage: b.plot()
            Graphics object consisting of 12 graphics primitives

        Here is the above example placed inside the ambient space
        of type `C_2`::

        .. PLOT::
            :width: 300 px

            MV = crystals.infinity.MVPolytopes(['C', 2])
            b = MV.highest_weight_vector().f_string([1,2,1,2,2,2,1,1,1,1,2,1])
            L = RootSystem(['C', 2, 1]).ambient_space()
            p = L.plot(reflection_hyperplanes=False, bounding_box=[[-8,2], [-8,2]])
            p += b.plot()
            p.axes(False)
            sphinx_plot(p)
        """
        if P is None:
            P = self.parent().weight_lattice_realization()
        return P.plot_mv_polytope(self, **options)

class MVPolytopes(PBWCrystal):
    """
    The crystal of Mirković-Vilonen (MV) polytopes.

    INPUT:

    - ``cartan_type`` -- a Cartan type
    """
    def __init__(self, cartan_type):
        """
        Initialize ``self``.

        EXAMPLES::

            sage: MV = crystals.infinity.MVPolytopes(['B', 2])
            sage: TestSuite(MV).run()
        """
        PBWCrystal.__init__(self, cartan_type)
        self._latex_options = {"projection": True,
                               "mark_endpoints": True,
                               "P": self.weight_lattice_realization(),
                               "circle_size": 0.1}

    def _repr_(self):
        """
        Return a string representation of ``self``.

        EXAMPLES::

            sage: crystals.infinity.MVPolytopes(['F', 4])
            MV polytopes of type ['F', 4]
        """
        return "MV polytopes of type {}".format(self._cartan_type)

    def set_latex_options(self, **kwds):
        r"""
        Set the latex options for the elements of ``self``.

        INPUT:

        - ``projection`` -- the projection; set to ``True`` to use the
          default projection of the specified weight lattice realization
          (initial: ``True``)
        - ``P`` -- the weight lattice realization to use (initial: the
          weight lattice realization of ``self``)
        - ``mark_endpoints`` -- whether to mark the endpoints (initial: ``True``)
        - ``circle_size`` -- the size of the endpoint circles (initial: 0.1)

        EXAMPLES::

            sage: MV = crystals.infinity.MVPolytopes(['C', 2])
            sage: P = RootSystem(['C', 2]).weight_lattice()
            sage: b = MV.highest_weight_vector().f_string([1,2,1,2])
            sage: latex(b)
            \begin{tikzpicture}
            \draw (0, 0) -- (-1, 1) -- (-1, 1) -- (-2, 0) -- (-2, -2);
            \draw (0, 0) -- (0, -2) -- (-1, -3) -- (-1, -3) -- (-2, -2);
            \draw[fill=black] (0, 0) circle (0.1);
            \draw[fill=black] (-2, -2) circle (0.1);
            \end{tikzpicture}
            sage: MV.set_latex_options(P=P, circle_size=float(0.2))
            sage: latex(b)
            \begin{tikzpicture}
            \draw (0, 0) -- (-2, 1) -- (-2, 1) -- (-2, 0) -- (0, -2);
            \draw (0, 0) -- (2, -2) -- (2, -3) -- (2, -3) -- (0, -2);
            \draw[fill=black] (0, 0) circle (0.2);
            \draw[fill=black] (0, -2) circle (0.2);
            \end{tikzpicture}
            sage: MV.set_latex_options(mark_endpoints=False)
            sage: latex(b)
            \begin{tikzpicture}
            \draw (0, 0) -- (-2, 1) -- (-2, 1) -- (-2, 0) -- (0, -2);
            \draw (0, 0) -- (2, -2) -- (2, -3) -- (2, -3) -- (0, -2);
            \end{tikzpicture}
            sage: MV.set_latex_options(P=MV.weight_lattice_realization(),
            ....:                      circle_size=0.2,
            ....:                      mark_endpoints=True)
        """
        if "projection" in kwds:
            self._latex_options["projection"] = True
            del kwds["projection"]

        if 'P' in kwds:
            self._latex_options['P'] = kwds['P']
            del kwds['P']

        if "mark_endpoints" in kwds:
            self._latex_options["mark_endpoints"] = kwds["mark_endpoints"]
            del kwds["mark_endpoints"]

        if "circle_size" in kwds:
            self._latex_options["circle_size"] = kwds["circle_size"]
            del kwds["circle_size"]

        if kwds:
            raise ValueError("invalid latex option")

    def latex_options(self):
        """
        Return the latex options of ``self``.

        EXAMPLES::

            sage: MV = crystals.infinity.MVPolytopes(['F', 4])
            sage: MV.latex_options()
            {'P': Ambient space of the Root system of type ['F', 4],
             'circle_size': 0.1,
             'mark_endpoints': True,
             'projection': True}
        """
        from copy import copy
        return copy(self._latex_options)

    Element = MVPolytope

