r"""
Homsets between simplicial complexes

EXAMPLES:

::

    sage: S = simplicial_complexes.Sphere(1)
    sage: T = simplicial_complexes.Sphere(2)
    sage: H = Hom(S,T)
    sage: f = {0:0,1:1,2:3}
    sage: x = H(f)
    sage: x
    Simplicial complex morphism {0: 0, 1: 1, 2: 3} from Simplicial complex with vertex set (0, 1, 2) and facets {(1, 2), (0, 2), (0, 1)} to Simplicial complex with vertex set (0, 1, 2, 3) and facets {(0, 2, 3), (0, 1, 2), (1, 2, 3), (0, 1, 3)}
    sage: x.is_injective()
    True
    sage: x.is_surjective()
    False
    sage: x.image()
    Simplicial complex with vertex set (0, 1, 3) and facets {(1, 3), (0, 3), (0, 1)}
    sage: from sage.homology.simplicial_complex import Simplex
    sage: s = Simplex([1,2])
    sage: x(s)
    (1, 3)

TESTS::

    sage: S = simplicial_complexes.Sphere(1)
    sage: T = simplicial_complexes.Sphere(2)
    sage: H = Hom(S,T)
    sage: loads(dumps(H))==H
    True

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

import sage.categories.homset
import sage.homology.simplicial_complex_morphism as simplicial_complex_morphism

def is_SimplicialComplexHomset(x):
    """
    Return True if and only if x is a simplicial complex homspace.

    EXAMPLES::

        sage: H = Hom(SimplicialComplex(2),SimplicialComplex(3))
        sage: H
        Set of Morphisms from Simplicial complex with vertex set (0, 1, 2) and facets {()} to Simplicial complex with vertex set (0, 1, 2, 3) and facets {()} in Category of simplicial complexes
        sage: from sage.homology.simplicial_complex_homset import is_SimplicialComplexHomset
        sage: is_SimplicialComplexHomset(H)
        True
    """
    return isinstance(x, SimplicialComplexHomset)

class SimplicialComplexHomset(sage.categories.homset.Homset):
    def __call__(self, f):
        """
        INPUT:
            f -- a dictionary with keys exactly the vertices of the domain and values vertices of the codomain

        EXAMPLE::

            sage: S = simplicial_complexes.Sphere(3)
            sage: T = simplicial_complexes.Sphere(2)
            sage: f = {0:0,1:1,2:2,3:2,4:2}
            sage: H = Hom(S,T)
            sage: x = H(f)
            sage: x
            Simplicial complex morphism {0: 0, 1: 1, 2: 2, 3: 2, 4: 2} from Simplicial complex with vertex set (0, 1, 2, 3, 4) and 5 facets to Simplicial complex with vertex set (0, 1, 2, 3) and facets {(0, 2, 3), (0, 1, 2), (1, 2, 3), (0, 1, 3)}


        """
        return simplicial_complex_morphism.SimplicialComplexMorphism(f,self.domain(),self.codomain())

    def diagonal_morphism(self,rename_vertices=True):
        """
        Returns the diagonal morphism in Hom(S,SxS).

        EXAMPLES::

            sage: S = simplicial_complexes.Sphere(2)
            sage: H = Hom(S,S.product(S))
            sage: d = H.diagonal_morphism()
            sage: d
            Simplicial complex morphism {0: 'L0R0', 1: 'L1R1', 2: 'L2R2', 3: 'L3R3'} from Simplicial complex with vertex set (0, 1, 2, 3) and facets {(0, 2, 3), (0, 1, 2), (1, 2, 3), (0, 1, 3)} to Simplicial complex with 16 vertices and 96 facets

            sage: T = SimplicialComplex(3)
            sage: U = T.product(T,rename_vertices = False)
            sage: G = Hom(T,U)
            sage: e = G.diagonal_morphism(rename_vertices = False)
            sage: e
            Simplicial complex morphism {0: (0, 0), 1: (1, 1), 2: (2, 2), 3: (3, 3)} from Simplicial complex with vertex set (0, 1, 2, 3) and facets {()} to Simplicial complex with 16 vertices and facets {()}

        """

        if self._codomain == self._domain.product(self._domain,rename_vertices=rename_vertices):
            X = self._domain.product(self._domain,rename_vertices=rename_vertices)
            f = dict()
            if rename_vertices:
                for i in self._domain.vertices().set():
                    f[i] = "L"+str(i)+"R"+str(i)
            else:
                for i in self._domain.vertices().set():
                    f[i] = (i,i)
            return simplicial_complex_morphism.SimplicialComplexMorphism(f,self._domain,X)
        else:
            raise TypeError, "Diagonal morphism is only defined for Hom(X,XxX)."

    def identity(self):
        """
        Returns the identity morphism of Hom(S,S).

        EXAMPLES::

            sage: S = simplicial_complexes.Sphere(2)
            sage: H = Hom(S,S)
            sage: i = H.identity()
            sage: i.is_identity()
            True

            sage: T = SimplicialComplex(3,[[0,1]])
            sage: G = Hom(T,T)
            sage: G.identity()
            Simplicial complex morphism {0: 0, 1: 1, 2: 2, 3: 3} from Simplicial complex with vertex set (0, 1, 2, 3) and facets {(0, 1)} to Simplicial complex with vertex set (0, 1, 2, 3) and facets {(0, 1)}

        """
        if self.is_endomorphism_set():
            f = dict()
            for i in self._domain.vertices().set():
                f[i]=i
            return simplicial_complex_morphism.SimplicialComplexMorphism(f,self._domain,self._codomain)
        else:
            raise TypeError, "Identity map is only defined for endomorphism sets."

    def an_element(self):
        """
        Returns a (non-random) element of self.

        EXAMPLES::

            sage: S = simplicial_complexes.KleinBottle()
            sage: T = simplicial_complexes.Sphere(5)
            sage: H = Hom(S,T)
            sage: x = H.an_element()
            sage: x
            Simplicial complex morphism {0: 0, 1: 0, 2: 0, 'R3': 0, 'L4': 0, 'L5': 0, 'L3': 0, 'R5': 0, 'R4': 0} from Simplicial complex with 9 vertices and 18 facets to Simplicial complex with vertex set (0, 1, 2, 3, 4, 5, 6) and 7 facets

        """
        X_vertices = self._domain.vertices().set()
        try:
            i = self._codomain.vertices().set().__iter__().next()
        except StopIteration:
            if len(X_vertices) == 0:
                return dict()
            else:
                raise TypeError, "There are no morphisms from a non-empty simplicial complex to an empty simplicial comples."
        f = dict()
        for x in X_vertices:
            f[x]=i
        return simplicial_complex_morphism.SimplicialComplexMorphism(f,self._domain,self._codomain)
