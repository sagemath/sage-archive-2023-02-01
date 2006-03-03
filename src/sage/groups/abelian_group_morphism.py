"""nodoctest
Homomorphisms of abelian groups

TODO:
    * remove nodoctest above to doctest this file
    * there must be a homspace first
    * there should be hom and Hom methods in abelian group
    * using gap for this is probably much more cumbersome slower than doing things directly.

AUTHORS:
    -- David Joyner (2006-03-03): initial version
"""

#*****************************************************************************
#   Copyright (C) 2006 David Joyner and William Stein <wstein@ucsd.edu>
#
#   Distributed under the terms of the GNU General Public License (GPL)
#
#                    http://www.gnu.org/licenses/
#*****************************************************************************

from sage.groups.permgroup import *
from sage.interfaces.gap import *
from sage.categories.morphism import *
from sage.categories.homset import *

from sage.groups.abelian_group import word_problem
from sage.misc.misc import prod

def is_AbelianGroupMorphism(f):
    return isinstance(f, AbelianGroupMorphism);

class AbelianGroupMap(Morphism):
    """
    A set-theoretic map between AbelianGroups.
    """
    def __init__(self, parent):
	Morphism.__init__(self, parent)

    def _repr_type(self):
        return "AbelianGroup"

class AbelianGroupMorphism:
    """
    Some python code for wrapping GAP's GroupHomomorphismByImages
    function for abelian groups. Returns "fail" if
    gens does not generate self or if the map does not
    extend to a group homomorphism, self --> other.

    EXAMPLES:
        sage: from abelian_group import AbelianGroup
        sage: G = AbelianGroup(3,[2,3,4],names="abc"); G
        Abelian group on 3 generators (a, b, c) with invariants [2, 3, 4]
        sage: a,b,c = G.gens()
        sage: H = AbelianGroup(2,[2,3],names="xy"); H
        Abelian group on 2 generators (x, y) with invariants [2, 3]
        sage: x,y = H.gens()

        sage: from abelian_group_morphism import AbelianGroupMorphism
        sage: phi = AbelianGroupMorphism(H,G,[x,y],[a,b])

#    There is a homomorphism from H to G but not from G to H:
#
#        sage: phi = AbelianGroupMorphism_im_gens(G,H,[a*b,a*c],[x,y])
#------------------------------------------------------------
#Traceback (most recent call last):
#  File "<ipython console>", line 1, in ?
#  File ".abeliangp_hom.sage.py", line 737, in __init__
#    raise TypeError, "Sorry, the orders of the corresponding elements in %s, %s must be equal."%(genss,imgss)
#TypeError: Sorry, the orders of the corresponding elements in [a*b, a*c], [x, y] must be equal.
#
#        sage: phi = AbelianGroupMorphism_im_gens(G,H,[a*b,(a*c)^2],[x*y,y])
#------------------------------------------------------------
#Traceback (most recent call last):
#  File "<ipython console>", line 1, in ?
#  File ".abeliangp_hom.sage.py", line 730, in __init__
#    raise TypeError, "Sorry, the list %s must generate G."%genss
#TypeError: Sorry, the list [a*b, c^2] must generate G.


    AUTHOR: David Joyner (2-2006)
    """
    def __init__(self, G, H, genss, imgss ):
        if len(genss) != len(imgss):
            raise TypeError, "Sorry, the lengths of %s, %s must be equal."%(genss,imgss)
        self.domain = G
        self.codomain = H
        if not(G.is_abelian()):
            raise TypeError, "Sorry, the groups must be abelian groups."
    	if not(H.is_abelian()):
            raise TypeError, "Sorry, the groups must be abelian groups."
        G_domain = G.subgroup(genss)
        if G_domain.order() != G.order():
            raise TypeError, "Sorry, the list %s must generate G."%genss
        self.domain_invs = G.invariants()
        self.codomaininvs = H.invariants()
        self.domaingens = genss
        self.codomaingens = imgss
        for i in range(len(self.domaingens)):
            if (self.domaingens[i]).order() != (self.codomaingens[i]).order():
                raise TypeError, "Sorry, the orders of the corresponding elements in %s, %s must be equal."%(genss,imgss)

    def _repr_type(self):
        return "AbelianGroup"

    def domain(self):
        return self._domain

    def range(self):
        return self._codomain

    def codomain(self):
        return self._codomain

    def kernel(self):
        raise NotImplementedError

    def image(self, J):
        """
        J must be a subgroup of G. Computes the subgroup of
        H which is the image of J.
        """
        raise NotImplementedError

    def __call__( self, g ):
        """
    	Some python code for wrapping GAP's Images
    	function but only for permutation groups.
    	Returns an error if g is not in G.

    	EXAMPLES:
            sage: from sage.groups.abelian_group import AbelianGroup
            sage: G = AbelianGroup(3, [2,3,4], names="abc"); G
            Abelian group on 3 generators (a, b, c) with invariants [2, 3, 4]
            sage: a,b,c = G.gens()
            sage: H = AbelianGroup(2, [2,3], names="xy"); H
            Abelian group on 2 generators (x, y) with invariants [2, 3]
            sage: x,y = H.gens()

            sage: from sage.groups.abelian_group_morphism import AbelianGroupMorphism
            sage: phi = AbelianGroupMorphism(H,G,[x,y],[a,b])
            sage: phi(x*y)
            a*b
    	"""
        G = g.parent()
        w = word_problem(self.domaingens,g,display=False)
        n = len(w)
        h = prod([self.codomaingens[i]**(w[i][1]) for i in range(n)])
    	return h





