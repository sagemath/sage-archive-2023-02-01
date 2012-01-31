r"""
Morphisms of simplicial complexes

AUTHORS:

- Benjamin Antieau <d.ben.antieau@gmail.com> (2009.06)

This module implements morphisms of simplicial complexes. The input is given
by a dictionary on the vertex set or the effective vertex set of a simplicial complex.
The initialization checks that faces are sent to faces.

There is also the capability to create the fiber product of two morphisms with the same codomain.

EXAMPLES::

    sage: S = SimplicialComplex(5,[[0,2],[1,5],[3,4]])
    sage: H = Hom(S,S.product(S))
    sage: H.diagonal_morphism()
    Simplicial complex morphism {0: 'L0R0', 1: 'L1R1', 2: 'L2R2', 3: 'L3R3', 4: 'L4R4', 5: 'L5R5'} from Simplicial complex with vertex set (0, 1, 2, 3, 4, 5) and facets {(3, 4), (1, 5), (0, 2)} to Simplicial complex with 36 vertices and 18 facets

    sage: S = SimplicialComplex(5,[[0,2],[1,5],[3,4]])
    sage: T = SimplicialComplex(4,[[0,2],[1,3]])
    sage: f = {0:0,1:1,2:2,3:1,4:3,5:3}
    sage: H = Hom(S,T)
    sage: x = H(f)
    sage: x.image()
    Simplicial complex with vertex set (0, 1, 2, 3) and facets {(1, 3), (0, 2)}
    sage: x.is_surjective()
    False
    sage: x.is_injective()
    False
    sage: x.is_identity()
    False

    sage: S = simplicial_complexes.Sphere(2)
    sage: H = Hom(S,S)
    sage: i = H.identity()
    sage: i.image()
    Simplicial complex with vertex set (0, 1, 2, 3) and facets {(0, 2, 3), (0, 1, 2), (1, 2, 3), (0, 1, 3)}
    sage: i.is_surjective()
    True
    sage: i.is_injective()
    True
    sage: i.is_identity()
    True

    sage: S = simplicial_complexes.Sphere(2)
    sage: H = Hom(S,S)
    sage: i = H.identity()
    sage: j = i.fiber_product(i)
    sage: j
    Simplicial complex morphism {'L1R1': 1, 'L3R3': 3, 'L2R2': 2, 'L0R0': 0} from Simplicial complex with 4 vertices and 4 facets to Simplicial complex with vertex set (0, 1, 2, 3) and facets {(0, 2, 3), (0, 1, 2), (1, 2, 3), (0, 1, 3)}

    sage: S = simplicial_complexes.Sphere(2)
    sage: T = S.product(SimplicialComplex(1,[[0,1]]),rename_vertices = False)
    sage: H = Hom(T,S)
    sage: T
    Simplicial complex with 8 vertices and 12 facets
    sage: T.vertices()
    ((0, 0), (0, 1), (1, 0), (1, 1), (2, 0), (2, 1), (3, 0), (3, 1))
    sage: f = {(0, 0): 0, (0, 1): 0, (1, 0): 1, (1, 1): 1, (2, 0): 2, (2, 1): 2, (3, 0): 3, (3, 1): 3}
    sage: x = H(f)
    sage: U = simplicial_complexes.Sphere(1)
    sage: G = Hom(U,S)
    sage: U
    Simplicial complex with vertex set (0, 1, 2) and facets {(1, 2), (0, 2), (0, 1)}
    sage: g = {0:0,1:1,2:2}
    sage: y = G(g)
    sage: z = y.fiber_product(x)
    sage: z                                     # this is the mapping path space
    Simplicial complex morphism {'L2R(2, 0)': 2, 'L2R(2, 1)': 2, 'L0R(0, 0)': 0, 'L0R(0, 1)': 0, 'L1R(1, 0)': 1, 'L1R(1, 1)': 1} from Simplicial complex with 6 vertices and 6 facets to Simplicial complex with vertex set (0, 1, 2, 3) and facets {(0, 2, 3), (0, 1, 2), (1, 2, 3), (0, 1, 3)}


"""

#*****************************************************************************
# Copyright (C) 2009 D. Benjamin Antieau <d.ben.antieau@gmail.com>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#
#    This code is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty
#    of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
#
#  See the GNU General Public License for more details; the full text
#  is available at:
#
#                  http://www.gnu.org/licenses/
#
#*****************************************************************************

import sage.homology.simplicial_complex as simplicial_complex
import sage.matrix.all as matrix
from sage.structure.sage_object import SageObject
from sage.rings.integer_ring import ZZ
from sage.homology.chain_complex_morphism import ChainComplexMorphism
from sage.combinat.permutation import Permutation
from sage.algebras.steenrod.steenrod_algebra_misc import convert_perm

def is_SimplicialComplexMorphism(x):
    """
    Returns True if and only if x is a morphism of simplicial complexes.

    EXAMPLES::

        sage: from sage.homology.simplicial_complex_morphism import is_SimplicialComplexMorphism
        sage: S = SimplicialComplex(5,[[0,1],[3,4]])
        sage: H = Hom(S,S)
        sage: f = {0:0,1:1,2:2,3:3,4:4,5:5}
        sage: x = H(f)
        sage: is_SimplicialComplexMorphism(x)
        True

    """
    return isinstance(x,SimplicialComplexMorphism)

class SimplicialComplexMorphism(SageObject):
    """
    An element of this class is a morphism of simplicial complexes.
    """
    def __init__(self,f,X,Y):
        """
        Input is a dictionary f, the domain, and the codomain.

        One can define the dictionary either on the vertices of X or on the effective vertices of X (X.effective_vertices()). Note that this
        difference does matter. For instance, it changes the result of the image method, and hence it changes the result of the is_surjective method as well.
        This is because two SimplicialComplexes with the same faces but different vertex sets are not equal.

        EXAMPLES::

            sage: S = SimplicialComplex(5,[[0,1],[3,4]])
            sage: H = Hom(S,S)
            sage: f = {0:0,1:1,2:2,3:3,4:4,5:5}
            sage: g = {0:0,1:1,3:3,4:4}
            sage: x = H(f)
            sage: y = H(g)
            sage: x==y
            False
            sage: x.image()
            Simplicial complex with vertex set (0, 1, 2, 3, 4, 5) and facets {(3, 4), (0, 1)}
            sage: y.image()
            Simplicial complex with vertex set (0, 1, 3, 4) and facets {(3, 4), (0, 1)}
            sage: x.image()==y.image()
            False

        """
        if not isinstance(X,simplicial_complex.SimplicialComplex) or not isinstance(Y,simplicial_complex.SimplicialComplex):
            raise ValueError, "X and Y must be SimplicialComplexes."
        if not set(f.keys())==X._vertex_set.set() and not set(f.keys())==X.effective_vertices().set():
            raise ValueError, "f must be a dictionary from the vertex set of X to single values in the vertex set of Y."
        dim = X.dimension()
        Y_faces = Y.faces()
        for k in range(dim+1):
            for i in X.faces()[k]:
                tup = i.tuple()
                fi = []
                for j in tup:
                    fi.append(f[j])
                v = simplicial_complex.Simplex(set(fi))
            if not v in Y_faces[v.dimension()]:
                raise ValueError, "f must be a dictionary from the vertices of X to the vertices of Y."
        self._vertex_dictionary = f
        self._domain = X
        self._codomain = Y

    def __eq__(self,x):
        """
        Returns True if and only if self == x.

        EXAMPLES::

            sage: S = simplicial_complexes.Sphere(2)
            sage: H = Hom(S,S)
            sage: i = H.identity()
            sage: i
            Simplicial complex morphism {0: 0, 1: 1, 2: 2, 3: 3} from Simplicial complex with vertex set (0, 1, 2, 3) and facets {(0, 2, 3), (0, 1, 2), (1, 2, 3), (0, 1, 3)} to Simplicial complex with vertex set (0, 1, 2, 3) and facets {(0, 2, 3), (0, 1, 2), (1, 2, 3), (0, 1, 3)}
            sage: f = {0:0,1:1,2:2,3:2}
            sage: j = H(f)
            sage: i==j
            False

            sage: T = SimplicialComplex(3,[[1,2]])
            sage: T
            Simplicial complex with vertex set (0, 1, 2, 3) and facets {(1, 2)}
            sage: G = Hom(T,T)
            sage: k = G.identity()
            sage: g = {0:0,1:1,2:2,3:3}
            sage: l = G(g)
            sage: k==l
            True

        """
        if not isinstance(x,SimplicialComplexMorphism) or self._codomain != x._codomain or self._domain != x._domain or self._vertex_dictionary != x._vertex_dictionary:
            return False
        else:
            return True

    def __call__(self,x,orientation=False):
        """
        Input is a simplex of the domain. Output is the image simplex.

        If optional argument ``orientation`` is True, return a pair
        ``(image simplex, oriented)`` where ``oriented`` is 1 or `-1`
        depending on whether the map preserves or reverses the
        orientation of the image simplex.

        EXAMPLES::

            sage: S = simplicial_complexes.Sphere(2)
            sage: T = simplicial_complexes.Sphere(3)
            sage: S
            Simplicial complex with vertex set (0, 1, 2, 3) and facets {(0, 2, 3), (0, 1, 2), (1, 2, 3), (0, 1, 3)}
            sage: T
            Simplicial complex with vertex set (0, 1, 2, 3, 4) and 5 facets
            sage: f = {0:0,1:1,2:2,3:3}
            sage: H = Hom(S,T)
            sage: x = H(f)
            sage: from sage.homology.simplicial_complex import Simplex
            sage: x(Simplex([0,2,3]))
            (0, 2, 3)

        An orientation-reversing example::

            sage: X = SimplicialComplex(1, [[0,1]])
            sage: g = Hom(X,X)({0:1, 1:0})
            sage: g(Simplex([0,1]))
            (0, 1)
            sage: g(Simplex([0,1]), orientation=True)
            ((0, 1), -1)
        """
        dim = self._domain.dimension()
        if not isinstance(x,simplicial_complex.Simplex) or x.dimension() > dim or not x in self._domain.faces()[x.dimension()]:
            raise ValueError, "x must be a simplex of the source of f"
        tup=x.tuple()
        fx=[]
        for j in tup:
            fx.append(self._vertex_dictionary[j])
        if orientation:
            if len(set(fx)) == len(tup):
                oriented = Permutation(convert_perm(fx)).signature()
            else:
                oriented = 1
            return (simplicial_complex.Simplex(set(fx)), oriented)
        else:
            return simplicial_complex.Simplex(set(fx))

    def _repr_(self):
        """
        Print representation

        EXAMPLES::

            sage: S = simplicial_complexes.Sphere(2)
            sage: H = Hom(S,S)
            sage: i = H.identity()
            sage: i
            Simplicial complex morphism {0: 0, 1: 1, 2: 2, 3: 3} from Simplicial complex with vertex set (0, 1, 2, 3) and facets {(0, 2, 3), (0, 1, 2), (1, 2, 3), (0, 1, 3)} to Simplicial complex with vertex set (0, 1, 2, 3) and facets {(0, 2, 3), (0, 1, 2), (1, 2, 3), (0, 1, 3)}
            sage: i._repr_()
            'Simplicial complex morphism {0: 0, 1: 1, 2: 2, 3: 3} from Simplicial complex with vertex set (0, 1, 2, 3) and facets {(0, 2, 3), (0, 1, 2), (1, 2, 3), (0, 1, 3)} to Simplicial complex with vertex set (0, 1, 2, 3) and facets {(0, 2, 3), (0, 1, 2), (1, 2, 3), (0, 1, 3)}'


        """
        return "Simplicial complex morphism " + str(self._vertex_dictionary) + " from " + self._domain._repr_() + " to " + self._codomain._repr_()

    def associated_chain_complex_morphism(self,base_ring=ZZ,augmented=False,cochain=False):
        """
        Returns the associated chain complex morphism of self.

        EXAMPLES::

            sage: S = simplicial_complexes.Sphere(1)
            sage: T = simplicial_complexes.Sphere(2)
            sage: H = Hom(S,T)
            sage: f = {0:0,1:1,2:2}
            sage: x = H(f)
            sage: x
            Simplicial complex morphism {0: 0, 1: 1, 2: 2} from Simplicial complex with vertex set (0, 1, 2) and facets {(1, 2), (0, 2), (0, 1)} to Simplicial complex with vertex set (0, 1, 2, 3) and facets {(0, 2, 3), (0, 1, 2), (1, 2, 3), (0, 1, 3)}
            sage: a = x.associated_chain_complex_morphism()
            sage: a
            Chain complex morphism from Chain complex with at most 2 nonzero terms over Integer Ring to Chain complex with at most 3 nonzero terms over Integer Ring
            sage: a._matrix_dictionary
            {0: [0 0 0]
            [0 1 0]
            [0 0 1]
            [1 0 0],
             1: [0 0 0]
            [0 1 0]
            [0 0 0]
            [1 0 0]
            [0 0 0]
            [0 0 1],
             2: []}
            sage: x.associated_chain_complex_morphism(augmented=True)
            Chain complex morphism from Chain complex with at most 3 nonzero terms over Integer Ring to Chain complex with at most 4 nonzero terms over Integer Ring
            sage: x.associated_chain_complex_morphism(cochain=True)
            Chain complex morphism from Chain complex with at most 3 nonzero terms over Integer Ring to Chain complex with at most 2 nonzero terms over Integer Ring
            sage: x.associated_chain_complex_morphism(augmented=True,cochain=True)
            Chain complex morphism from Chain complex with at most 4 nonzero terms over Integer Ring to Chain complex with at most 3 nonzero terms over Integer Ring
            sage: x.associated_chain_complex_morphism(base_ring=GF(11))
            Chain complex morphism from Chain complex with at most 2 nonzero terms over Finite Field of size 11 to Chain complex with at most 3 nonzero terms over Finite Field of size 11

        Some simplicial maps which reverse the orientation of a few simplices::

            sage: g = {0:1, 1:2, 2:0}
            sage: H(g).associated_chain_complex_morphism()._matrix_dictionary
            {0: [0 0 0]
             [1 0 0]
             [0 1 0]
             [0 0 1], 1: [ 0  0  0]
             [-1  0  0]
             [ 0  0  0]
             [ 0  0  1]
             [ 0  0  0]
             [ 0 -1  0], 2: []}

            sage: X = SimplicialComplex(1, [[0, 1]])
            sage: Hom(X,X)({0:1, 1:0}).associated_chain_complex_morphism()._matrix_dictionary
            {0: [0 1]
             [1 0], 1: [-1]}
        """
        max_dim = max(self._domain.dimension(),self._codomain.dimension())
        min_dim = min(self._domain.dimension(),self._codomain.dimension())
        matrices = {}
        if augmented is True:
            m = matrix.Matrix(base_ring,1,1,1)
            if not cochain:
                matrices[-1] = m
            else:
                matrices[-1] = m.transpose()
        for dim in range(min_dim+1):
#             X_faces = list(self._domain.faces()[dim])
#             Y_faces = list(self._codomain.faces()[dim])
            X_faces = list(self._domain.n_cells(dim))
            Y_faces = list(self._codomain.n_cells(dim))
            num_faces_X = len(X_faces)
            num_faces_Y = len(Y_faces)
            mval = [0 for i in range(num_faces_X*num_faces_Y)]
            for i in X_faces:
                y, oriented = self(i, orientation=True)
                if y.dimension() < dim:
                    pass
                else:
                    mval[X_faces.index(i)+(Y_faces.index(y)*num_faces_X)] = oriented
            m = matrix.Matrix(base_ring,num_faces_Y,num_faces_X,mval,sparse=True)
            if not cochain:
                matrices[dim] = m
            else:
                matrices[dim] = m.transpose()
        for dim in range(min_dim+1,max_dim+1):
            try:
                l1 = len(self._codomain.n_cells(dim))
            except KeyError:
                l1 = 0
            try:
                l2 = len(self._domain.n_cells(dim))
            except KeyError:
                l2 = 0
            m = matrix.zero_matrix(base_ring,l1,l2,sparse=True)
            if not cochain:
                matrices[dim] = m
            else:
                matrices[dim] = m.transpose()
        if not cochain:
            return ChainComplexMorphism(matrices,\
                    self._domain.chain_complex(base_ring=base_ring,augmented=augmented,cochain=cochain),\
                    self._codomain.chain_complex(base_ring=base_ring,augmented=augmented,cochain=cochain))
        else:
            return ChainComplexMorphism(matrices,\
                    self._codomain.chain_complex(base_ring=base_ring,augmented=augmented,cochain=cochain),\
                    self._domain.chain_complex(base_ring=base_ring,augmented=augmented,cochain=cochain))

    def image(self):
        """
        Computes the image simplicial complex of f.

        EXAMPLES::

            sage: S = SimplicialComplex(3,[[0,1],[2,3]])
            sage: T = SimplicialComplex(1,[[0,1]])
            sage: f = {0:0,1:1,2:0,3:1}
            sage: H = Hom(S,T)
            sage: x = H(f)
            sage: x.image()
            Simplicial complex with vertex set (0, 1) and facets {(0, 1)}

            sage: S = SimplicialComplex(2)
            sage: H = Hom(S,S)
            sage: i = H.identity()
            sage: i.image()
            Simplicial complex with vertex set (0, 1, 2) and facets {()}
            sage: i.is_surjective()
            True
            sage: S = SimplicialComplex(5,[[0,1]])
            sage: T = SimplicialComplex(3,[[0,1]])
            sage: f = {0:0,1:1}
            sage: g = {0:0,1:1,2:2,3:3,4:4,5:5}
            sage: H = Hom(S,T)
            sage: x = H(f)
            sage: y = H(g)
            sage: x == y
            False
            sage: x.image()
            Simplicial complex with vertex set (0, 1) and facets {(0, 1)}
            sage: y.image()
            Simplicial complex with vertex set (0, 1, 2, 3, 4, 5) and facets {(0, 1)}

        """
        fa = [self(i) for i in self._domain.facets()]
        return simplicial_complex.SimplicialComplex(set(self._vertex_dictionary.values()),fa,maximality_check=True)

    def domain(self):
        """
        Returns the domain of the morphism.

        EXAMPLES::

            sage: S = SimplicialComplex(3,[[0,1],[2,3]])
            sage: T = SimplicialComplex(1,[[0,1]])
            sage: f = {0:0,1:1,2:0,3:1}
            sage: H = Hom(S,T)
            sage: x = H(f)
            sage: x.domain()
            Simplicial complex with vertex set (0, 1, 2, 3) and facets {(2, 3), (0, 1)}

        """
        return self._domain

    def codomain(self):
        """
        Returns the codomain of the morphism.

        EXAMPLES::

            sage: S = SimplicialComplex(3,[[0,1],[2,3]])
            sage: T = SimplicialComplex(1,[[0,1]])
            sage: f = {0:0,1:1,2:0,3:1}
            sage: H = Hom(S,T)
            sage: x = H(f)
            sage: x.codomain()
            Simplicial complex with vertex set (0, 1) and facets {(0, 1)}

        """
        return self._codomain

    def is_surjective(self):
        """
        Returns True if and only if self is surjective.

        EXAMPLES::

            sage: S = SimplicialComplex(3,[(0,1,2)])
            sage: S
            Simplicial complex with vertex set (0, 1, 2, 3) and facets {(0, 1, 2)}
            sage: T = SimplicialComplex(2,[(0,1)])
            sage: T
            Simplicial complex with vertex set (0, 1, 2) and facets {(0, 1)}
            sage: H = Hom(S,T)
            sage: x = H({0:0,1:1,2:1,3:2})
            sage: x.is_surjective()
            True

            sage: S = SimplicialComplex(3,[[0,1],[2,3]])
            sage: T = SimplicialComplex(1,[[0,1]])
            sage: f = {0:0,1:1,2:0,3:1}
            sage: H = Hom(S,T)
            sage: x = H(f)
            sage: x.is_surjective()
            True

        """
        return self._codomain == self.image()

    def is_injective(self):
        """
        Returns True if and only if self is injective.

        EXAMPLES::

            sage: S = simplicial_complexes.Sphere(1)
            sage: T = simplicial_complexes.Sphere(2)
            sage: U = simplicial_complexes.Sphere(3)
            sage: H = Hom(T,S)
            sage: G = Hom(T,U)
            sage: f = {0:0,1:1,2:0,3:1}
            sage: x = H(f)
            sage: g = {0:0,1:1,2:2,3:3}
            sage: y = G(g)
            sage: x.is_injective()
            False
            sage: y.is_injective()
            True

        """
        v = [self._vertex_dictionary[i[0]] for i in self._domain.faces()[0]]
        for i in v:
            if v.count(i) > 1:
                return False
        return True

    def is_identity(self):
        """
        If x is an identity morphism, returns True. Otherwise, False.

        EXAMPLES::

            sage: T = simplicial_complexes.Sphere(1)
            sage: G = Hom(T,T)
            sage: T
            Simplicial complex with vertex set (0, 1, 2) and facets {(1, 2), (0, 2), (0, 1)}
            sage: j = G({0:0,1:1,2:2})
            sage: j.is_identity()
            True

            sage: S = simplicial_complexes.Sphere(2)
            sage: T = simplicial_complexes.Sphere(3)
            sage: H = Hom(S,T)
            sage: f = {0:0,1:1,2:2,3:3}
            sage: x = H(f)
            sage: x
            Simplicial complex morphism {0: 0, 1: 1, 2: 2, 3: 3} from Simplicial complex with vertex set (0, 1, 2, 3) and facets {(0, 2, 3), (0, 1, 2), (1, 2, 3), (0, 1, 3)} to Simplicial complex with vertex set (0, 1, 2, 3, 4) and 5 facets
            sage: x.is_identity()
            False

        """
        if self._domain != self._codomain:
            return False
        else:
            f = dict()
            for i in self._domain.vertices().set():
                f[i] = i
            if self._vertex_dictionary != f:
                return False
            else:
                return True

    def fiber_product(self,other,rename_vertices = True):
        """
        Fiber product of self and other. Both morphisms should have the same codomain. The method returns a morphism of simplicial complexes, which
        is the morphism from the space of the fiber product to the codomain.

        EXAMPLES::

            sage: S = SimplicialComplex(2,[[0,1],[1,2]])
            sage: T = SimplicialComplex(2,[[0,2]])
            sage: U = SimplicialComplex(1,[[0,1]])
            sage: H = Hom(S,U)
            sage: G = Hom(T,U)
            sage: f = {0:0,1:1,2:0}
            sage: g = {0:0,1:1,2:1}
            sage: x = H(f)
            sage: y = G(g)
            sage: z = x.fiber_product(y)
            sage: z
            Simplicial complex morphism {'L1R2': 1, 'L2R0': 0, 'L0R0': 0} from Simplicial complex with vertex set ('L0R0', 'L1R2', 'L2R0') and facets {('L2R0',), ('L0R0', 'L1R2')} to Simplicial complex with vertex set (0, 1) and facets {(0, 1)}

        """
        if self._codomain != other._codomain:
            raise ValueError, "self and other must have the same codomain."
        X = self._domain.product(other._domain,rename_vertices = rename_vertices)
        v = []
        f = dict()
        eff1 = self._domain.effective_vertices()
        eff2 = other._domain.effective_vertices()
        for i in eff1:
            for j in eff2:
                if self(simplicial_complex.Simplex([i])) == other(simplicial_complex.Simplex([j])):
                    if rename_vertices:
                        v.append("L"+str(i)+"R"+str(j))
                        f["L"+str(i)+"R"+str(j)] = self._vertex_dictionary[i]
                    else:
                        v.append((i,j))
                        f[(i,j)] = self._vertex_dictionary[i]
        return SimplicialComplexMorphism(f,X.generated_subcomplex(v),self._codomain)
