r"""
Associahedron

.. todo::

    - fix adjacency matrix
    - edit graph method to get proper vertex labellings
    - UniqueRepresentation?

AUTHORS:

- Christian Stump
"""
#*****************************************************************************
#       Copyright (C) 2011-2012 Christian Stump <christian.stump@gmail.com>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#  The full text of the GPL is available at:
#
#                  http://www.gnu.org/licenses/
#*****************************************************************************
from sage.geometry.polyhedra import Polyhedron
from sage.combinat.root_system.cartan_type import CartanType
from sage.modules.free_module_element import vector

class Associahedron(Polyhedron):
    r"""
    The generalized associahedron is a polytopal complex with vertices in one-to-one correspondence
    with clusters in the cluster complex, and with edges between two vertices if and only if the associated two
    clusters intersect in codimension 1.

    The associahedron of type `A_n` is one way to realize the classical associahedron as defined in

    http://en.wikipedia.org/wiki/Associahedron.

    A polytopal realization of the associahedron can be found in [CFZ].
    The implementation is based on [CFZ, Theorem 1.5, Remark 1.6, and Corollary 1.9.].

    EXAMPLES::

        sage: Asso = Associahedron(['A',2]); Asso
        Generalized associahedron of type ['A', 2] with 5 vertices
        sage: sorted(Asso.Hrepresentation(), key=repr)
        [An inequality (-1, 0) x + 1 >= 0, An inequality (0, -1) x + 1 >= 0, An inequality (0, 1) x + 1 >= 0, An inequality (1, 0) x + 1 >= 0, An inequality (1, 1) x + 1 >= 0]
        sage: Asso.Vrepresentation()
        [A vertex at (-1, 1), A vertex at (1, 1), A vertex at (1, -1), A vertex at (0, -1), A vertex at (-1, 0)]

        sage: Associahedron(['B',2])
        Generalized associahedron of type ['B', 2] with 6 vertices

   The two pictures of [CFZ] can be recovered with::

        sage: Asso = Associahedron(['A',3]); Asso
        Generalized associahedron of type ['A', 3] with 14 vertices
        sage: Asso.plot()

        sage: Asso = Associahedron(['B',3]); Asso
        Generalized associahedron of type ['B', 3] with 20 vertices
        sage: Asso.plot()

    TESTS::

        sage: sorted(Associahedron(['A',3]).vertices())
        [[-3/2, 0, -1/2], [-3/2, 0, 3/2], [-3/2, 1, -3/2], [-3/2, 2, -3/2], [-3/2, 2, 3/2], [-1/2, -1, -1/2], [-1/2, 0, -3/2], [1/2, -2, 1/2], [1/2, -2, 3/2], [3/2, -2, 1/2], [3/2, -2, 3/2], [3/2, 0, -3/2], [3/2, 2, -3/2], [3/2, 2, 3/2]]
        sage: sorted(Associahedron(['B',3]).vertices())
        [[-3, 0, 0], [-3, 0, 3], [-3, 2, -2], [-3, 4, -3], [-3, 5, -3], [-3, 5, 3], [-2, 1, -2], [-2, 3, -3], [-1, -2, 0], [-1, -1, -1], [1, -4, 1], [1, -3, 0], [2, -5, 2], [2, -5, 3], [3, -5, 2], [3, -5, 3], [3, -3, 0], [3, 3, -3], [3, 5, -3], [3, 5, 3]]

        sage: Associahedron(['A',4]).f_vector()
        (1, 42, 84, 56, 14, 1)
        sage: Associahedron(['B',4]).f_vector()
        (1, 70, 140, 90, 20, 1)

    REFERENCES:

    - [CFZ] Chapoton, Fomin, Zelevinsky - Polytopal realizations of generalized associahedra, arXiv:0202004.
    """
    def __init__(self, cartan_type):
        """
        TESTS::

            sage: Asso = Associahedron(['A',2]); Asso
            Generalized associahedron of type ['A', 2] with 5 vertices
            sage: TestSuite(Asso).run()

        """
        self._cartan_type = CartanType( cartan_type )
        assert self._cartan_type.is_finite()
        root_space = self._cartan_type.root_system().root_space()
        # TODO: generalize this as a method of root lattice realization
        rhocheck = sum( beta.associated_coroot() for beta in root_space.positive_roots() )/2
        I = root_space.index_set()
        inequalities = []
        for orbit in root_space.almost_positive_roots_decomposition():
            c = rhocheck.coefficient(orbit[0].leading_support())
            for beta in orbit:
                inequalities.append( [c] + [ beta.coefficient(i) for i in I ] )
        Polyhedron.__init__(self,ieqs=inequalities)
        # check that there are non non trivial facets
        assert self.n_facets() == len(inequalities)

    def _repr_(self):
        r"""
        Returns a string representation of self.

        EXAMPLES::

            sage: Associahedron(['A',3])._repr_()
            "Generalized associahedron of type ['A', 3] with 14 vertices"
        """
        return 'Generalized associahedron of type %s with %s vertices'%(self._cartan_type,self.n_vertices())

    def cartan_type(self):
        r"""
        Returns the Cartan type of self.

        EXAMPLES::

            sage: Associahedron(['A',3]).cartan_type()
            ['A', 3]
        """
        return self._cartan_type

    def vertices_in_root_space(self):
        r"""
        Returns the vertices of ``self`` as elements in the root space

        EXAMPLES::

            sage: Asso = Associahedron(['A',2])

            sage: Asso.vertices()
            [[-1, 1], [1, 1], [1, -1], [0, -1], [-1, 0]]

            sage: Asso.vertices_in_root_space()
            [-alpha[1] + alpha[2], alpha[1] + alpha[2], alpha[1] - alpha[2], -alpha[2], -alpha[1]]
        """
        root_space = self._cartan_type.root_system().root_space()
        return [ root_space.from_vector(vector(V)) for V in self.vertex_generator() ]
