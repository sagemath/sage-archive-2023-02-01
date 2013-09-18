"""
Klyachko Bundles.

EXAMPLES::

    sage: X = toric_varieties.dP6xdP6()
    sage: TX = X.tangent_bundle()
    sage: Alt2TX = TX.exterior_power(2)
    sage: # Alt2TX.cohomology(dim=True, weight=(0,0,0,0))

    sage: K = LineBundle(X, X.K())
    sage: antiK = LineBundle(X, -X.K())
    sage: (Alt2TX * K).cohomology(dim=True, weight=(0,0,0,0))

    sage: G_sum = TX + TrivialBundle(X,2)
    sage: V_sum = G_sum.wedge(2) * K
    sage: V_sum.cohomology(dim=True, weight=(0,0,0,0))
    (0, 0, 18, 16, 1)
    sage: Gtilde = G_sum.deformation()
    sage: V = Gtilde.wedge(2) * K
    sage: V.cohomology(dim=True, weight=(0,0,0,0))
    (0, 0, 3, 0, 0)


REFERENCES:

..  [Klyachko]
    Klyachko, Anton Alexander:
    Equivariant Bundles on Toral Varieties
"""

#*****************************************************************************
#       Copyright (C) 2011 Volker Braun <vbraun.name@gmail.com>
#  Distributed under the terms of the GNU General Public License (GPL)
#  as published by the Free Software Foundation; either version 2 of 
#  the License, or (at your option) any later version.  
#                  http://www.gnu.org/licenses/
#*****************************************************************************

from sage.structure.all import SageObject
from sage.rings.all import QQ, ZZ
from sage.misc.all import uniq, cached_method
from sage.modules.free_module import FreeModule, VectorSpace
from sage.matrix.constructor import vector, matrix, block_matrix, zero_matrix
from sage.geometry.cone import is_Cone, IntegralRayCollection

from sage.modules.filtered_vector_space import (
    FilteredVectorSpace, is_FilteredVectorSpace )



def is_KlyachkoBundle(X):
    """
    """
    return isinstance(X, Bundle_class)

##############################################################################


def Bundle(toric_variety, *args, **kwds):
    r"""
    Construct a Klyacho bundle

    EXAMPLES::

        sage: P1 = toric_varieties.P1()
        sage: from sage.schemes.toric.bundle.klyachko import Bundle
        sage: Bundle(P1, [(1,0,0),(0,1,0),(0,0,1)], [(1,[1,2]), (3,[1]), (4,)], [(2,[0]), (3,)])
    """
    base_ring = kwds.pop('base_ring', None)
    check = kwds.pop('check', True)
    assert len(kwds)==0, 'Unkown keyword argument: '+str(kwds.keys()[0])

    if all(is_FilteredVectorSpace(V) for V in args):
        assert base_ring is None
        return Bundle_class(toric_variety, args, check=check)

    if base_ring is None:
        base_ring = QQ
    return Bundle_from_rays_filtrations(toric_variety, args[0], args[1:], 
                                        base_ring=base_ring, check=check)


def Bundle_from_filtered_vector_spaces(toric_variety, filtrations, base_ring=QQ, check=True):
    """
    Construct a Klyachko bundle from filtration data.
    """
    if check:
        assert all(F._rays == filtrations[0]._rays for F in filtrations)
    return Bundle_class(toric_variety, filtrations, check=check)


def Bundle_from_rays_filtrations(toric_variety, rays, filtrations, base_ring=QQ, check=True):
    """
    Construct a Klyachko bundle from filtration data.
    """
    filtrations = [ FilteredVectorSpace(rays,f,base_ring=base_ring, check=check) 
                    for f in filtrations ]
    return Bundle_class(toric_variety, filtrations, check=check)
    


##############################################################################

class Bundle_class(SageObject):
    r"""
    A toric bundle using Klyachko's representation.
    """
    
    def __init__(self, toric_variety, filtrations, check=True):
        r"""
        The Python constructor.
        """
        self._variety = toric_variety
        self._filt = filtrations
        self._base_ring = self._filt[0].base_ring()
        self._rank = self._filt[0].dimension()
        self._vectorspace = VectorSpace(self._base_ring, self._rank)
        
        if not check: return
        if len(filtrations)!=self._variety.fan().nrays():
            raise ValueError('There must be one filtered vector space for each ray of the fan.')
        for i,F in enumerate(self._filt):
            msg = 'The '+str(i)+'-th filtration '+str(F)+' '
            if not F.is_finite():
                raise ValueError(msg+'has infinitely many non-trivial terms.')
            if F.dimension() != self._rank:
                raise ValueError(msg+'has a differnt dimension.')
                

    def variety(self):
        r"""
        Return the base toric variety.

        OUTPUT:

        A toric variety.

        EXAMPLES::

            sage: X = toric_varieties.P2()
            sage: X.bundle.tangent_bundle()
            Rank 2 bundle on 2-d CPR-Fano toric variety covered by 3 affine patches.
            sage: X.tangent_bundle().variety() is X
            True
        """
        return self._variety

    def base_ring(self):
        r"""
        Return the base field.

        EXAMPLES::

            sage: toric_varieties.P2().tangent_bundle().base_ring()
            Rational Field
        """
        return self._base_ring

    def fiber(self):
        r"""
        Return the generic fiber of the vector bundle.

        OUTPUT:

        A vector space over :meth:`base_ring`.

        EXAMPLES::

            sage: toric_varieties.P2().tangent_bundle().fiber()
            Vector space of dimension 2 over Rational Field
        """
        return self._vectorspace
        
    def rank(self):
        r"""
        Return the rank of the vector bundle.

        OUTPUT:

        Integer.

        EXAMPLES::

            sage: toric_varieties.P2().tangent_bundle().rank()
            2
        """
        return self._rank

    def _repr_(self):
        r"""
        Return a string representation.

        OUTPUT:
        
        String.

        EXAMPLES::
        
            sage: toric_varieties.P2().tangent_bundle()
            Rank 2 bundle on 2-d CPR-Fano toric variety covered by 3 affine patches.
        """
        s = 'Rank '+str(self.rank())+' bundle on '+str(self._variety)+'.'
        return s
        
    def filtration(self, alpha):
        r"""
        Return the filtration associated to the ray alpha.

        INPUT:

        - ``alpha`` -- Integer or a `N`-lattice point or a
          one-dimensional cone. Specifies a ray of the fan of the
          toric variety, either via its index or its generator.

        EXAMPLES:
        
            sage: toric_varieties.dP6().tangent_bundle().filtration(0)
            QQ^2 >= QQ^1 >= 0
        """
        try:
            i = ZZ(alpha)
        except TypeError:
            fan = self._variety.fan()
            if is_Cone(alpha):
                cone = fan.embed(alpha)
            else:
                cone = fan.cone_containing(alpha)
            assert cone.dim() == 1
            i = cone.ambient_ray_indices()[0]
        return self._filt[i]

    def filtration_degree(self, alpha, i):
        r"""
        Return the vector subspace ``E^\alpha(i)``.

        EXAMPLES::

            sage: toric_varieties.dP6().tangent_bundle().filtration_degree(0,1)
            Vector space of degree 2 and dimension 1 over Rational Field
            Basis matrix:
            [0 1]
        """
        return self.filtration(alpha).get_degree(i)

    def filtration_intersection(self, sigma, i):
        r"""
        Return the intersection of the filtered subspaces.
        
        INPUT:

        - ``sigma`` -- a cone of the fan of the base toric variety.
        
        - ``i`` -- integer. The filtration degree.

        OUPUT:

        Let the cone be spanned by the rays `\sigma=\langle r_1,\dots,
        r_k\rangle`. This method returns the intersection

        .. math::

            \bigcap_{r\in \{r_1,\dots,r_k\}}
            E^{r}(i)
        
        EXAMPLES::

            sage: X = toric_varieties.P2()
            sage: fan = X.fan()
            sage: V = X.tangent_bundle()
            sage: V.filtration_intersection(fan(1)[0], 1)
            Vector space of degree 2 and dimension 1 over Rational Field
            Basis matrix:
            [1 0]
            sage: V.filtration_intersection(fan(2)[0], 1)
            Vector space of degree 2 and dimension 0 over Rational Field
            Basis matrix:
            []
        """
        sigma = self._variety.fan().embed(sigma)
        V = self.fiber()
        for alpha in sigma.ambient_ray_indices():
            V = V.intersection(self.filtration_degree(alpha,i))
        return V

    def E_degree(self, alpha, m):
        r"""
        Return the vector subspace `E^\alpha(i)`.

        EXAMPLES::

            sage: X = toric_varieties.P2()
            sage: ray = X.fan(1)[0]
            sage: M = X.fan().dual_lattice()
            sage: V = X.tangent_bundle()
            sage: V.E_degree(ray, (1,0))
            Vector space of degree 2 and dimension 1 over Rational Field
            Basis matrix:
            [1 0]
        """
        M = self._variety.fan().dual_lattice()
        m = M(m)
        try:
            ray = self._variety.fan().ray(alpha)
            ray_index = alpha
        except TypeError:
            cone = self._variety.fan().cone_containing(alpha)
            assert cone.dim() == 1
            ray = cone.ray(0)
            ray_index = cone.ambient_ray_indices()[0]
        return self.filtration_degree(ray_index, ray*m)

    def E_intersection(self, sigma, m):
        r"""
        Return the vector subspace ``E^\sigma(m)``.
        
        See [Klyachko]_, equation 4.1.

        INPUT:

        - ``sigma`` -- a cone of the fan of the base toric variety.

        - ``m`` -- a point in the dual lattice of the fan.

        OUPUT:

        The subspace `E^\sigma(m)`
        
        EXAMPLES::

            sage: X = toric_varieties.P2()
            sage: fan = X.fan()
            sage: V = X.tangent_bundle()
            sage: V.E_intersection(fan(1)[0], (1,0))
            Vector space of degree 2 and dimension 1 over Rational Field
            Basis matrix:
            [1 0]
            sage: V.E_intersection(fan(2)[0], (-1,1))
            Vector space of degree 2 and dimension 1 over Rational Field
            Basis matrix:
            [0 1]

        For the empty cone, this is always the whole vector space::
        
            sage: V.E_intersection(fan(0)[0], (1,0))
            Vector space of dimension 2 over Rational Field
        """
        sigma = self._variety.fan().embed(sigma)
        V = self.fiber()
        for alpha in sigma.ambient_ray_indices():
            V = V.intersection(self.E_degree(alpha,m))
        return V

    @cached_method
    def E_quotient(self, sigma, m):
        r"""
        Return the vector space quotient `E_\sigma(m)`.
        
        See [Klyachko]_, equation 4.1.

        INPUT:

        - ``sigma`` -- a cone of the fan of the base toric variety.

        - ``m`` -- a point in the dual lattice of the fan.

        OUPUT:

        The subspace `E_\sigma(m)`
        
        EXAMPLES::

            sage: X = toric_varieties.P2()
            sage: fan = X.fan()
            sage: M = fan.dual_lattice()
            sage: cone = fan(1)[0]
            sage: V = X.tangent_bundle()
            sage: V.E_quotient(cone, M(1,0))
            Vector space quotient V/W of dimension 1 over Rational Field where
            V: Vector space of dimension 2 over Rational Field
            W: Vector space of degree 2 and dimension 1 over Rational Field
            Basis matrix:
            [1 0]
            sage: V.E_quotient(fan(2)[0], (-1,1))
            Vector space quotient V/W of dimension 0 over Rational Field where
            V: Vector space of dimension 2 over Rational Field
            W: Vector space of degree 2 and dimension 2 over Rational Field
            Basis matrix:
            [1 0]
            [0 1]
        """
        sigma = self._variety.fan().embed(sigma)
        V = self.fiber()
        generators = []
        for alpha in sigma.ambient_ray_indices():
            generators.extend(self.E_degree(alpha,m).gens())
        return V.quotient(V.span(generators))

    def E_quotient_projection(self, sigma, tau, m):
        r"""
        Return the projection map `E_\sigma(m) \to E_\tau(m)` where
        `\sigma` is a face of `\tau`.

        INPUT:

        - ``sigma`` -- a cone of the fan of the base toric variety.

        - ``tau`` -- a cone of the fan containing ``sigma``.

        - ``m`` -- a point in the dual lattice of the fan.

        OUTPUT:

        The restriction map

        .. math::
        
            E_\sigma(m) \to E_\tau(m)

        EXAMPLES::

            sage: P3 = toric_varieties.P(3)
            sage: from sage.schemes.toric.bundle.klyachko import Bundle
            sage: rays = [(1,0,0),(0,1,0),(0,0,1)]
            sage: V = Bundle(P3, rays, [(1,[1,2]), (2,[1]), (3,[])], [(1,[0])], [(1,[])], [(1,[])])
            sage: tau = Cone([(1,0,0), (0,1,0)])
            sage: sigma = Cone([(1,0,0)])
            sage: M = P3.fan().dual_lattice()
            sage: m = M(2,1,0)
            sage: V.E_quotient(sigma, m)
            sage: V.E_quotient(tau, m)
            sage: V.E_quotient_projection(sigma, tau, m)
        """
        if not sigma.is_face_of(tau):
            raise ValueError, 'The cone sigma is not a face of the cone tau.'
        E_sigma = self.E_quotient(sigma, m)
        E_tau = self.E_quotient(tau, m)
        images = [ E_tau(E_sigma.lift(g)) for g in E_sigma.gens() ]
        return E_sigma.hom(images, codomain=E_tau)

    def cohomology_complex(self, m):
        r"""
        Return the "cohomology complex" `C^*(m)`

        See [Klyachko]_, equation 4.2.

        EXAMPLES::

            sage: P3 = toric_varieties.P(3)
            sage: from sage.schemes.toric.bundle.klyachko import Bundle
            sage: rays = [(1,0,0),(0,1,0),(0,0,1)]
            sage: V = Bundle(P3, rays, [(1,[1,2]), (2,[1]), (3,[])], [(-1,[1]), (2,[])], [(-1,[])], [(-1,[])])
            sage: tau = Cone([(1,0,0), (0,1,0)])
            sage: sigma = Cone([(1,0,0)])
            sage: M = P3.fan().dual_lattice()
            sage: m = M(1,1,0); m.set_immutable()
            sage: V.cohomology_complex(m)
            Chain complex with at most 4 nonzero terms over Rational Field

            sage: F = CyclotomicField(3)
            sage: P3 = toric_varieties.P(3).change_ring(F)
            sage: V = Bundle(P3, rays, [(1,[1,2]), (2,[1]), (3,[])], [(-1,[1]), (2,[])], [(-1,[])], [(-1,[])])
            sage: V.cohomology_complex(m)
        """
        fan = self._variety.fan()
        C = fan.complex()
        CV = []
        F = self.base_ring()
        for dim in range(1,fan.dim()+1):
            codim = fan.dim() - dim
            d_C = C.differential(codim)
            d_V = []
            for j in range(0, d_C.ncols()):
                tau = fan(dim)[j]
                d_V_row = []
                for i in range(0, d_C.nrows()):
                    sigma = fan(dim-1)[i]
                    if sigma.is_face_of(tau):
                        pr = self.E_quotient_projection(sigma, tau, m)
                        d = d_C[i,j] * pr.matrix().transpose()
                    else:
                        E_sigma = self.E_quotient(sigma, m)
                        E_tau = self.E_quotient(tau, m)
                        d = zero_matrix(F, E_tau.dimension(), E_sigma.dimension())
                    d_V_row.append(d)
                d_V.append(d_V_row)
            #print dim, ':\n', d_V, '\n'
            d_V = block_matrix(d_V, ring=F)
            CV.append(d_V)
            #print dim, ': ', d_V.nrows(), 'x', d_V.ncols(), '\n', d_V
        from sage.homology.chain_complex import ChainComplex
        return ChainComplex(CV, base_ring=self.base_ring())

    def cohomology(self, d=None, weight=None, dim=False):
        r"""
        Return the bundle cohomology groups.

        INPUT:
        
        - ``d`` -- ``None`` (default) or an integer. The degree of
          the cohomology group.

        - ``weight`` -- ``None`` or a point in the dual lattice of the
          fan. The weight of the cohomology group.

        - ``dim`` -- Boolean (default: ``False``). Whether to return
          vector spaces or only their dimension.

        OUTPUT:

        The cohomology group of degree ``d`` and weight ``m``. 

        * If no weight is specified, the unweighted group (sum over all
          weights) is returned. 

        * If no degree is specified, a dictionary whose keys are
          integers and whose values are the cohomology groups is
          returned. If, in addition, ``dim=True``, then an integral
          vector of the dimensions is returned.

        EXAMPLES::

            sage: V = toric_varieties.P2().tangent_bundle()
            sage: V.cohomology(weight=(0,0), dim=True)
            (2, 0, 0)
            sage: for i,j in CartesianProduct(range(-2,3), range(-2,3)):
            ...       HH = V.cohomology(m=(i,j), dim=True)
            ...       if HH.is_zero(): continue
            ...       print 'H^*i(P^2, TP^2)_M('+str(i)+','+str(j)+') =', HH
            H^*i(P^2, TP^2)_M(-1,0) = (1, 0, 0)
            H^*i(P^2, TP^2)_M(-1,1) = (1, 0, 0)
            H^*i(P^2, TP^2)_M(0,-1) = (1, 0, 0)
            H^*i(P^2, TP^2)_M(0,0)  = (2, 0, 0)
            H^*i(P^2, TP^2)_M(0,1)  = (1, 0, 0)
            H^*i(P^2, TP^2)_M(1,-1) = (1, 0, 0)
            H^*i(P^2, TP^2)_M(1,0)  = (1, 0, 0)
        """
        if d:
            return self.cohomology(weight=m, dim=dim)[d]
        C = self.cohomology_complex(weight)
        space_dim = self._variety.dimension()
        C_homology = C.homology()
        HH = dict()
        for d in range(0,space_dim+1):
            try: 
                HH[d] = C_homology[d]
            except KeyError:
                HH[d] = FreeModule(self.base_ring(), 0)
        if dim:
            HH = vector(ZZ, [ HH[i].rank() for i in range(0,space_dim+1) ])
        return HH

    def __cmp__(self, other):
        """
        Compare ``self`` and ``other``
        
        .. warning::
        
            This method tests whether the underlying representation is
            the same. Use :meth:`is_isomorphic` to test for
            mathematical equivalence.

        EXAMPLES
            
            sage: X = toric_varieties.P2()
            sage: cmp(TrivialBundle(X,2), TrivialBundle(X,1))
            1
            sage: cmp(TrivialBundle(X,2), TrivialBundle(X,1) + TrivialBundle(X,1))
            0

            sage: T_X = X.tangent_bundle()
            sage: O_X = TrivialBundle(X,1)
            sage: cmp(T_X+O_X, O_X+T_X)
            -1
            sage: (T_X+O_X).is_isomorphic(O_X+T_X)
            True
        """
        c = cmp(type(self), type(other))
        if c!=0: return c
        c = cmp(self.variety(), other.variety())
        if c!=0: return c
        for i in range(0,len(self._filt)):
            c = cmp(self._filt[i], other._filt[i])
            if c!=0: return c
        return 0

    def is_isomorphic(self, other):
        """
        Test whether ``self`` and ``other`` are isomorphic bundles.
        """
        if not is_KlyachkoBundle(other):
            return False
        if not self.variety().is_isomorphic(other.variety()):
            return False
        # FIXME: this is not enough
        for i in range(0,len(self._filt)):
            if not self._filt[i].is_isomorphic(other._filt[i]):
                return False
        return True
        
    def direct_sum(self, other):
        """
        Return the sum of two vector bundles.
        
        EXAMPLES::

            sage: X = toric_varieties.P2()
            sage: X.tangent_bundle() + TrivialBundle(X,1)
            Rank 3 bundle on 2-d CPR-Fano toric variety covered by 3 affine patches.
            
            sage: TrivialBundle(X,2).is_isomorphic(TrivialBundle(X,1) + TrivialBundle(X,1))
            True
        """
        if not self.variety() == other.variety():
            raise ValueError('The bundles must be over the same base toric variety.')
        filt = [ self._filt[i] + other._filt[i]
                 for i in range(0, self.variety().fan().nrays()) ]
        return Bundle_class(self.variety(), filt, check=True)

    __add__ = direct_sum

    def tensor_product(self, other):
        """
        Return the sum of two vector bundles.
        
        EXAMPLES::

            sage: X = toric_varieties.P2()
            sage: X.tangent_bundle() + TrivialBundle(X,1)
            Rank 3 bundle on 2-d CPR-Fano toric variety covered by 3 affine patches.
            
            sage: TrivialBundle(X,1).is_isomorphic(TrivialBundle(X,1) * TrivialBundle(X,1))
            True
        """
        if not self.variety() == other.variety():
            raise ValueError('The bundles must be over the same base toric variety.')
        filt = [ self._filt[i] * other._filt[i]
                 for i in range(0, self.variety().fan().nrays()) ]
        return Bundle_class(self.variety(), filt, check=True)

    __mul__ = tensor_product

    def exterior_power(self, n):
        """
        EXAMPLES:

            sage: X = toric_varieties.P2_123()
            sage: TX = X.tangent_bundle()
            sage: antiK = LineBundle(X, -X.K() )
            sage: TX.exterior_power(2).is_isomorphic(antiK)
            True
        """
        filt = [ self._filt[i].exterior_power(n)
                 for i in range(0, len(self._filt)) ]
        return Bundle_class(self.variety(), filt, check=True)

    wedge = exterior_power
        
    def symmetric_power(self, n):
        filt = [ self._filt[i].symmetric_power(n)
                 for i in range(0, len(self._filt)) ]
        return Bundle_class(self.variety(), filt, check=True)
        
    def deformation(self, perturbed_rays=None):
        """
        Return a generic deformation of the bundle.

        EXAMPLES::

           sage: P1 = toric_varieties.P1()
           sage: H = P1.divisor(0)
           sage: V = LineBundle(P1, H) + LineBundle(P1, -H)
           sage: V.cohomology(dim=True, weight=(0,))
           (1, 0)
           sage: Vtilde = V.deformation()
           sage: Vtilde.cohomology(dim=True, weight=(0,))
        """
        rays = self._filt[0]._rays
        assert all(F._rays == rays for F in self._filt)
        if perturbed_rays is None:
            from sage.modules.free_module_element import random_vector
            perturbed_rays = [ r + random_vector(self.base_ring(), self.rank()) 
                               for r in rays ]
        else:
            assert len(rays)==len(perturbed_rays)
        filt = [ FilteredVectorSpace(perturbed_rays, F._filt, base_ring=self.base_ring())
                 for F in self._filt ]
        return Bundle_class(self.variety(), filt, check=True)
        
