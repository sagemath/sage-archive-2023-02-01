"""
Homomorphisms of Lie Algebras

AUTHORS:

- Travis Scrimshaw (07-15-2013): Initial implementation
"""
#*****************************************************************************
#  Copyright (C) 2013 Travis Scrimshaw <tscrim at ucdavis.edu>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#                  http://www.gnu.org/licenses/
#******************************************************************************

from sage.categories.morphism import Morphism
from sage.structure.sequence import Sequence, Sequence_generic

class LieAlgebraHomomorphism_im_gens(Morphism):
    r"""
    A homomorphism of Lie algebras.

    Let `\mathfrak{g}` and `\mathfrak{g}^{\prime}` be Lie algebras,
    a linear map `f \colon \mathfrak{g} \to \mathfrak{g}^{\prime}` is a
    homomorphism if `f([x, y]) = [f(x), f(y)]` for all `x, y \in \mathfrak{g}`.
    Thus homomorphisms are completely determined by the image of the
    generators of `\mathfrak{g}`.
    """
    def __init__(self, parent, im_gens, check=True):
        """
        EXAMPLES::

            sage: R.<x,y> = QQ[]
            sage: phi = R.hom([x,x+y]); phi
            Ring endomorphism of Multivariate Polynomial Ring in x, y over Rational Field
              Defn: x |--> x
                    y |--> x + y
            sage: type(phi)
            <type 'sage.rings.morphism.RingHomomorphism_im_gens'>

        Here's another example where the domain isn't free::
        
            sage: S.<xx,yy> = R.quotient(x - y)
            sage: phi = S.hom([xx+1,xx+1])

        Note that one has to specify valid images::
        
            sage: phi = S.hom([xx+1,xx-1])
            Traceback (most recent call last):
            ...
            TypeError: images do not define a valid homomorphism

        There is a check option, but it may be ignored in some cases
        -- it's purpose isn't so you can lie to Sage, but to sometimes
        speed up creation of a homomorphism::

            sage: phi = S.hom([xx+1, xx-1],check=False)
            Traceback (most recent call last):
            ...
            TypeError: images do not define a valid homomorphism
        """
        Morphism.__init__(self, parent)
        if not isinstance(im_gens, Sequence_generic):
            if not isinstance(im_gens, (tuple, list)):
                im_gens = [im_gens]
            im_gens = Sequence(im_gens, parent.codomain(), immutable=True)
        if check:
            if len(im_gens) != len(parent.domain().lie_algebra_generators()):
                raise ValueError("number of images must equal number of generators")
            # TODO: Implement a (meaningful) _is_valid_homomorphism_()
            #if not parent.domain()._is_valid_homomorphism_(parent.codomain(), im_gens):
            #    raise ValueError("relations do not all (canonically) map to 0 under map determined by images of generators.")
        if not im_gens.is_immutable():
            import copy
            im_gens = copy.copy(im_gens)
            im_gens.set_immutable()
        self.__im_gens = im_gens

    def _repr_type(self):
        """
        TESTS::

            sage: f
            sage: type(f)
            <type 'sage.algebras.lie_algebras.morphism.LieAlgebraHomomorphism'>
            sage: f._repr_type()
            'Lie algebra'
            sage: f
            Lie algebra endomorphism of Integer Ring
        """
        return "Lie algebra"

    def im_gens(self):
        """
        Return the images of the generators of the domain.

        OUTPUT:

        - ``list`` -- a copy of the list of gens (it is safe to change this)

        EXAMPLES::

            sage: 
            sage: f = R.hom([x,x+y])
            sage: f.im_gens()
            [x, x + y]
        """
        return list(self.__im_gens)

    def __hash__(self):
        """
        Return the hash of this morphism.

        EXAMPLES::

            sage: R.<x> = ZZ[]
            sage: s = R.hom([x+1])
            sage: type(s)
            <type 'sage.rings.morphism.RingHomomorphism_im_gens'>
            sage: hash(s) == hash(s)
            True
            sage: {s: 1}[s]
            1
        """
        return hash(self.__im_gens)

    def _repr_defn(self):
        """
        Used in constructing string representation of ``self``.

        EXAMPLES::

            sage: R.<x,y> = QQ[]; f = R.hom([x^2,x+y])
            sage: print f._repr_defn()
            x |--> x^2
            y |--> x + y
        """
        D = self.domain()
        ig = self.__im_gens
        return '\n'.join(['%s |--> %s'%(x, ig[i]) for i, x in enumerate(D.gens())])
            
    def _call_(self, x):
        """
        Evaluate this homomorphism at ``x``.

        EXAMPLES::

            sage: R.<x,y,z> = ZZ[]; f = R.hom([2*x,z,y])
            sage: f(x+2*y+3*z)             # indirect doctest
            2*x + 3*y + 2*z
        """
        return x._im_gens_(self.codomain(), self.im_gens())

