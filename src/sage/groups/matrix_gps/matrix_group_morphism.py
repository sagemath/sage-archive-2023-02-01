"""
Homomorphisms Between Matrix Groups

AUTHORS:

- David Joyner and William Stein (2006-03): initial version

- David Joyner (2006-05): examples

- Simon King (2011-01): cleaning and improving code
"""

#*****************************************************************************
#       Copyright (C) 2006 David Joyner and William Stein <wstein@gmail.com>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#  as published by the Free Software Foundation; either version 2 of
#  the License, or (at your option) any later version.
#                  http://www.gnu.org/licenses/
#*****************************************************************************

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
    large.

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
        phi0 = gap(self)
        if gap.eval("IsGroupHomomorphism(%s)"%phi0.name())!="true":
            raise ValueError,"The map "+str(gensG)+"-->"+str(imgsH)+" isn't a homomorphism."

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

    def _gap_init_(self):
        return self._gap_str

    def kernel(self):
        """
        Return the kernel of ``self``, i.e., a matrix group.

        EXAMPLES::

            sage: F = GF(7); MS = MatrixSpace(F,2,2)
            sage: F.multiplicative_generator()
            3
            sage: G = MatrixGroup([MS([3,0,0,1])])
            sage: a = G.gens()[0]^2
            sage: phi = G.hom([a])
            sage: phi.kernel()
            Matrix group over Finite Field of size 7 with 1 generators:
             [[[6, 0], [0, 1]]]

        """
        gap_ker = gap(self).Kernel()
        from sage.all import MatrixGroup
        F = self.domain().base_ring()
        return MatrixGroup([x._matrix_(F) for x in gap_ker.GeneratorsOfGroup()])

    def pushforward(self, J, *args,**kwds):
        """
        The image of an element or a subgroup.

        INPUT:

        ``J`` -- a subgroup or an element of the domain of ``self``.

        OUTPUT:

        The image of ``J`` under ``self``

        NOTE:

        ``pushforward`` is the method that is used when a map is called on
        anything that is not an element of its domain. For historical reasons,
        we keep the alias ``image()`` for this method.

        EXAMPLES::

            sage: F = GF(7); MS = MatrixSpace(F,2,2)
            sage: F.multiplicative_generator()
            3
            sage: G = MatrixGroup([MS([3,0,0,1])])
            sage: a = G.gens()[0]^2
            sage: phi = G.hom([a])
            sage: phi.image(G.gens()[0]) # indirect doctest
            [2 0]
            [0 1]
            sage: H = MatrixGroup([MS(a.list())])
            sage: H
            Matrix group over Finite Field of size 7 with 1 generators:
             [[[2, 0], [0, 1]]]

        The following tests against trac ticket #10659::

            sage: phi(H)   # indirect doctestest
            Matrix group over Finite Field of size 7 with 1 generators:
             [[[4, 0], [0, 1]]]
        """
        phi = gap(self)
        F = self.codomain().base_ring()
        from sage.all import MatrixGroup
        gapJ = gap(J)
        if gap.eval("IsGroup(%s)"%gapJ.name()) == "true":
            return MatrixGroup([x._matrix_(F) for x in phi.Image(gapJ).GeneratorsOfGroup()])
        return phi.Image(gapJ)._matrix_(F)

    image = pushforward

    def _call_( self, g ):
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

        TEST:

        The following tests that the call method was successfully
        improved in trac ticket #10659::

            sage: O = WeylGroup(['D',6])
            sage: r = prod(O.gens())
            sage: r_ = r^-1
            sage: f = O.hom([r*x*r_ for x in O.gens()])  # long time (19s on sage.math, 2011)
            sage: [f(x) for x in O.gens()]  # long time
            [
            [1 0 0 0 0 0]  [1 0 0 0 0 0]  [1 0 0 0 0 0]  [ 0  0  0  0 -1  0]
            [0 0 1 0 0 0]  [0 1 0 0 0 0]  [0 1 0 0 0 0]  [ 0  1  0  0  0  0]
            [0 1 0 0 0 0]  [0 0 0 1 0 0]  [0 0 1 0 0 0]  [ 0  0  1  0  0  0]
            [0 0 0 1 0 0]  [0 0 1 0 0 0]  [0 0 0 0 1 0]  [ 0  0  0  1  0  0]
            [0 0 0 0 1 0]  [0 0 0 0 1 0]  [0 0 0 1 0 0]  [-1  0  0  0  0  0]
            [0 0 0 0 0 1], [0 0 0 0 0 1], [0 0 0 0 0 1], [ 0  0  0  0  0  1],
            <BLANKLINE>
            [0 0 0 0 0 1]  [ 0  0  0  0  0 -1]
            [0 1 0 0 0 0]  [ 0  1  0  0  0  0]
            [0 0 1 0 0 0]  [ 0  0  1  0  0  0]
            [0 0 0 1 0 0]  [ 0  0  0  1  0  0]
            [0 0 0 0 1 0]  [ 0  0  0  0  1  0]
            [1 0 0 0 0 0], [-1  0  0  0  0  0]
            ]
            sage: f(O)  # long time
            Matrix group over Rational Field with 6 generators:
             [[[1, 0, 0, 0, 0, 0], [0, 0, 1, 0, 0, 0], [0, 1, 0, 0, 0, 0], [0, 0, 0, 1, 0, 0], [0, 0, 0, 0, 1, 0], [0, 0, 0, 0, 0, 1]], [[1, 0, 0, 0, 0, 0], [0, 1, 0, 0, 0, 0], [0, 0, 0, 1, 0, 0], [0, 0, 1, 0, 0, 0], [0, 0, 0, 0, 1, 0], [0, 0, 0, 0, 0, 1]], [[1, 0, 0, 0, 0, 0], [0, 1, 0, 0, 0, 0], [0, 0, 1, 0, 0, 0], [0, 0, 0, 0, 1, 0], [0, 0, 0, 1, 0, 0], [0, 0, 0, 0, 0, 1]], [[0, 0, 0, 0, -1, 0], [0, 1, 0, 0, 0, 0], [0, 0, 1, 0, 0, 0], [0, 0, 0, 1, 0, 0], [-1, 0, 0, 0, 0, 0], [0, 0, 0, 0, 0, 1]], [[0, 0, 0, 0, 0, 1], [0, 1, 0, 0, 0, 0], [0, 0, 1, 0, 0, 0], [0, 0, 0, 1, 0, 0], [0, 0, 0, 0, 1, 0], [1, 0, 0, 0, 0, 0]], [[0, 0, 0, 0, 0, -1], [0, 1, 0, 0, 0, 0], [0, 0, 1, 0, 0, 0], [0, 0, 0, 1, 0, 0], [0, 0, 0, 0, 1, 0], [-1, 0, 0, 0, 0, 0]]]

        """
        phi = gap(self)
        G = self.domain()
        F = G.base_ring()
        h = gap(g)
        return phi.Image(h)._matrix_(F)
