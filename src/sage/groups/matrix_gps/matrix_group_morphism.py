"""
Homomorphisms Between Matrix Groups

AUTHORS:

- David Joyner and William Stein (2006-03): initial version

- David Joyner (2006-05): examples
"""

#*****************************************************************************
#       Copyright (C) 2006 David Joyner and William Stein <wstein@gmail.com>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#
#                  http://www.gnu.org/licenses/
#*****************************************************************************

from sage.groups.group import Group
from sage.rings.all import IntegerRing, is_Ring, Integer
from sage.interfaces.gap import gap
from sage.categories.morphism import *
from sage.misc.latex import latex

class MatrixGroupMap(Morphism):
    """
    A set-theoretic map between matrix groups.
    """
    def __init__(self, parent):
        Morphism.__init__(self, parent)

    def _repr_type(self):
        return "MatrixGroup"

class MatrixGroupMorphism(MatrixGroupMap):
    pass

class MatrixGroupMorphism_im_gens(MatrixGroupMorphism):
    """
    Some python code for wrapping GAP's GroupHomomorphismByImages
    function but only for matrix groups. Can be expensive if G is
    large. Returns "fail" if gens does not generate self or if the map
    does not extend to a group homomorphism, self - other.

    TODO: what does it mean to return fail? It's a constructor for a
    class.

    EXAMPLES::

        sage: F = GF(5); MS = MatrixSpace(F,2,2)
        sage: G = MatrixGroup([MS([1,1,0,1])])
        sage: H = MatrixGroup([MS([1,0,1,1])])
        sage: phi = G.hom(H.gens())
        sage: phi
        Homomorphism : Matrix group over Finite Field of size 5 with 1 generators:
         [[[1, 1], [0, 1]]] --> Matrix group over Finite Field of size 5 with 1 generators:
         [[[1, 0], [1, 1]]]
        sage: phi(MS([1,1,0,1]))
        [1 0]
        [1 1]
        sage: F = GF(7); MS = MatrixSpace(F,2,2)
        sage: F.multiplicative_generator()
        3
        sage: G = MatrixGroup([MS([3,0,0,1])])
        sage: a = G.gens()[0]^2
        sage: phi = G.hom([a])
    """
    def __init__(self, homset, imgsH, check=True):
        MatrixGroupMorphism.__init__(self, homset)   # sets the parent
        G = homset.domain()
        H = homset.codomain()
        gaplist_gens = [gap(x) for x in G.gens()]
        gaplist_imgs = [gap(x) for x in imgsH]
        genss = '[%s]'%(','.join(str(v) for v in gaplist_gens))
        imgss = '[%s]'%(','.join(str(v) for v in gaplist_imgs))
        args = '%s, %s, %s, %s'%(G._gap_init_(), H._gap_init_(), genss, imgss)
        self._gap_str = 'GroupHomomorphismByImages(%s)'%args
        self._gap_hom_string = 'phi := %s'%self._gap_str
        phi0 = gap.eval(self._gap_hom_string)
        if phi0=="fail":
            raise ValueError,"The map "+str(gensG)+"-->"+str(imgsH)+" isn't a homomorphism."
        self.hom = gap.eval("phi")

    def _gap_init_(self):
        return self._gap_str

    def kernel(self):
        """
        EXAMPLES::

            sage: F = GF(7); MS = MatrixSpace(F,2,2)
            sage: F.multiplicative_generator()
            3
            sage: G = MatrixGroup([MS([3,0,0,1])])
            sage: a = G.gens()[0]^2
            sage: phi = G.hom([a])
            sage: phi.kernel()
            'Group([ [ [ Z(7)^3, 0*Z(7) ], [ 0*Z(7), Z(7)^0 ] ] ])'
            sage: phi.image(G.gens()[0])
            '[ [ Z(7)^2, 0*Z(7) ], [ 0*Z(7), Z(7)^0 ] ]'
        """
        cmd = self._gap_hom_string
        gap.eval(cmd)
        gap_ker = gap.eval("Kernel(phi)")
        return gap_ker

    def image(self, J):
        """
        J must be a subgroup of G. Computes the subgroup of H which is the
        image of J.

        EXAMPLES::

            sage: F = GF(7); MS = MatrixSpace(F,2,2)
            sage: F.multiplicative_generator()
            3
            sage: G = MatrixGroup([MS([3,0,0,1])])
            sage: a = G.gens()[0]^2
            sage: phi = G.hom([a])
            sage: phi.image(G.gens()[0])
            '[ [ Z(7)^2, 0*Z(7) ], [ 0*Z(7), Z(7)^0 ] ]'
            sage: H = MatrixGroup([MS(a.list())])
            sage: H
            Matrix group over Finite Field of size 7 with 1 generators:
            [[[2, 0], [0, 1]]]
            sage: phi.image(H)
            'Group([ [ [ Z(7)^4, 0*Z(7) ], [ 0*Z(7), Z(7)^0 ] ] ])'
        """
        cmd = self._gap_hom_string
        gap.eval(cmd)
        return gap.eval("Image( phi, "+str(J._gap_init_())+")")

    def _repr_(self):
        """
        EXAMPLES::

            sage: F = GF(5); MS = MatrixSpace(F,2,2)
            sage: G = MatrixGroup([MS([1,1,0,1])])
            sage: H = MatrixGroup([MS([1,0,1,1])])
            sage: phi = G.hom(H.gens())
            sage: phi
            Homomorphism : Matrix group over Finite Field of size 5 with 1 generators:
             [[[1, 1], [0, 1]]] --> Matrix group over Finite Field of size 5 with 1 generators:
             [[[1, 0], [1, 1]]]
            sage: phi(MS([1,1,0,1]))
            [1 0]
            [1 1]
        """
        return "Homomorphism : %s --> %s"%(self.domain(),self.codomain())

    def _latex_(self):
        r"""
        EXAMPLES::

            sage: F = GF(5); MS = MatrixSpace(F,2,2)
            sage: G = MatrixGroup([MS([1,1,0,1])])
            sage: phi = G.hom(G.gens())
            sage: print latex(phi)
            \left\langle \left(\begin{array}{rr}
            1 & 1 \\
            0 & 1
            \end{array}\right) \right\rangle \rightarrow{} \left\langle \left(\begin{array}{rr}
            1 & 1 \\
            0 & 1
            \end{array}\right) \right\rangle
        """
        return "%s \\rightarrow{} %s"%(latex(self.domain()), latex(self.codomain()))

    def __call__( self, g ):
        """
        Some python code for wrapping GAP's Images function for a matrix
        group G. Returns an error if g is not in G.

        EXAMPLES::

            sage: F = GF(5); MS = MatrixSpace(F,2,2)
            sage: g = MS([1,1,0,1])
            sage: G = MatrixGroup([g])
            sage: phi = G.hom(G.gens())
            sage: phi(G.0)
            [1 1]
            [0 1]
            sage: phi(G(g^2))
            [1 2]
            [0 1]

        ::

            sage: F = GF(5); MS = MatrixSpace(F,2,2)
            sage: gens = [MS([1,2,  -1,1]),MS([1,1,  0,1])]
            sage: G = MatrixGroup(gens)
            sage: phi = G.hom(G.gens())
            sage: phi(G.0)
            [1 2]
            [4 1]
            sage: phi(G.1)
            [1 1]
            [0 1]
        """
        cmd = self._gap_hom_string
        gap.eval(cmd)
        G = self.domain()
        F = G.base_ring()
        h = gap(g)
        return gap('Image(phi, %s)'%h.name())._matrix_(F)



