"""
Klyachko Bundles and Sheaves.

Klyachko bundles are torus-equivariant bundles on toric
varieties. That is, the action of the maximal torus on the toric
variety lifts to an action on the bundle. There is an equivalence of
categories between Klyachko bundles and multiple filtrations (one for
each ray of the fan) of a vector space. The multi-filtrations are
implemented in :mod:`sage.modules.multi_filtered_vector_space`.

EXAMPLES::

    sage: X = toric_varieties.dP6xdP6()
    sage: TX = X.sheaves.tangent_bundle()
    sage: Alt2TX = TX.exterior_power(2)
    sage: # Alt2TX.cohomology(dim=True, weight=(0,0,0,0))

    sage: K = X.sheaves.line_bundle(X.K())
    sage: antiK = X.sheaves.line_bundle(-X.K())
    sage: (Alt2TX * K).cohomology(dim=True, weight=(0,0,0,0))  # long time
    (0, 0, 18, 0, 0)

    sage: G_sum = TX + X.sheaves.trivial_bundle(2)
    sage: V_sum = G_sum.wedge(2) * K                      # long time
    sage: V_sum.cohomology(dim=True, weight=(0,0,0,0))    # long time 
    (0, 0, 18, 16, 1)
    sage: Gtilde = G_sum.deformation()
    sage: V = Gtilde.wedge(2) * K                     # long time
    sage: V.cohomology(dim=True, weight=(0,0,0,0))    # long time
    (0, 0, 3, 0, 0)

REFERENCES:

..  [Klyachko]
    Klyachko, Anton Alexander:
    Equivariant Bundles on Toral Varieties
"""

#*****************************************************************************
#       Copyright (C) 2013 Volker Braun <vbraun.name@gmail.com>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#  as published by the Free Software Foundation; either version 2 of 
#  the License, or (at your option) any later version.  
#                  http://www.gnu.org/licenses/
#*****************************************************************************


from sage.structure.all import SageObject
from sage.rings.all import QQ, ZZ
from sage.misc.all import uniq, cached_method
from sage.matrix.constructor import vector, matrix, block_matrix, zero_matrix
from sage.geometry.cone import is_Cone, IntegralRayCollection

from sage.modules.filtered_vector_space import FilteredVectorSpace, is_FilteredVectorSpace
from sage.modules.multi_filtered_vector_space import MultiFilteredVectorSpace
from sage.schemes.toric.sheaf.util import cached_method_immutable


def is_KlyachkoBundle(X):
    """
    """
    return isinstance(X, KlyachkoBundle_class)


def Bundle(toric_variety, multi_filtration, check=True):
    r"""
    Construct a Klyacho bundle

    INPUT:

    - ``toric_variety`` -- a toric variety. The base space of the bundle.

    - ``multi_filtration`` -- a multi-filtered vectors space with
      multiple filtrations being indexed by the one-dimensional cones
      of the fan. Either an instance of
      :func:`~sage.modules.multi_filtered_vector_space.MultiFilteredVectorSpace`
      or something (like a dictionary of ordinary filtered vector
      spaces).

    EXAMPLES::

        sage: P1 = toric_varieties.P1()
        sage: v1, v2, v3 = [(1,0,0),(0,1,0),(0,0,1)]
        sage: F1 = FilteredVectorSpace({1:[v1, v2, v3], 3:[v1]})
        sage: F2 = FilteredVectorSpace({0:[v1, v2, v3], 2:[v2, v3]})
        sage: P1 = toric_varieties.P1()
        sage: r1, r2 = P1.fan().rays()
        sage: F = MultiFilteredVectorSpace({r1:F1, r2:F2});  F
        Filtrations
            N(-1): QQ^3 >= QQ^2 >= QQ^2 >=  0   >= 0
             N(1): QQ^3 >= QQ^3 >= QQ^1 >= QQ^1 >= 0

    You should use the
    :meth:`~sage.schemes.toric.sheaf.constructor.SheafLibrary.Klyachko`
    method to construct instances::

        sage: P1.sheaves.Klyachko(F)
        Rank 3 bundle on 1-d CPR-Fano toric variety covered by 2 affine patches.

        sage: P1.sheaves.Klyachko({r1:F1, r2:F2})   # alternative
        Rank 3 bundle on 1-d CPR-Fano toric variety covered by 2 affine patches.

    The above is just a shorthand for::

        sage: from sage.schemes.toric.sheaf.klyachko import Bundle
        sage: Bundle(P1, F)
        Rank 3 bundle on 1-d CPR-Fano toric variety covered by 2 affine patches.
    """
    base_ring = toric_variety.base_ring()
    if not hasattr(multi_filtration, 'get_filtration'):
        # try to construct a MultiFilteredVectorSpace
        multi_filtration = MultiFilteredVectorSpace(multi_filtration, base_ring=base_ring, check=check)
    if multi_filtration.base_ring() != base_ring:
        multi_filtration = multi_filtration.change_ring(base_ring)
    return KlyachkoBundle_class(toric_variety, multi_filtration, check=check)



class KlyachkoBundle_class(SageObject):
    r"""
    A toric bundle using Klyachko's representation.

    INPUT:

    - ``toric_variety`` -- a toric variety. The base space of the bundle.

    - ``multi_filtration`` -- a
      :func:`~sage.modules.multi_filtered_vector_space.MultiFilteredVectorSpace`
      with index set the rays of the fan.

    EXAMPLES::


    """
    
    def __init__(self, toric_variety, multi_filtration, check=True):
        r"""
        The Python constructor.
        """
        self._variety = toric_variety
        self._filt = multi_filtration
        
        if not check: return
        from sage.sets.set import Set
        if multi_filtration.index_set() != Set(list(toric_variety.fan().rays())):
            raise ValueError('the index set of the multi-filtration must be'
                             ' all rays of the fan.')
        if not multi_filtration.is_exhaustive():
            raise ValueError('multi-filtration must be exhaustive')
        if not multi_filtration.is_separating():
            raise ValueError('multi-filtration must be separating')

    def variety(self):
        r"""
        Return the base toric variety.

        OUTPUT:

        A toric variety.

        EXAMPLES::

            sage: X = toric_varieties.P2()
            sage: V = X.sheaves.tangent_bundle(); V
            Rank 2 bundle on 2-d CPR-Fano toric variety covered by 3 affine patches.
            sage: V.variety() is X
            True
        """
        return self._variety

    def base_ring(self):
        r"""
        Return the base field.

        EXAMPLES::

            sage: T_P2 = toric_varieties.P2().sheaves.tangent_bundle()
            sage: T_P2.base_ring()
            Rational Field
        """
        return self._filt.base_ring()

    def fiber(self):
        r"""
        Return the generic fiber of the vector bundle.

        OUTPUT:

        A multi-filtered vector space over :meth:`base_ring`.

        EXAMPLES::

            sage: T_P2 = toric_varieties.P2().sheaves.tangent_bundle()
            sage: T_P2.fiber()
            Vector space of dimension 2 over Rational Field
        """
        from sage.modules.all import VectorSpace
        return VectorSpace(self.base_ring(), self.rank())
        
    def rank(self):
        r"""
        Return the rank of the vector bundle.

        OUTPUT:

        Integer.

        EXAMPLES::

            sage: T_P2 = toric_varieties.P2().sheaves.tangent_bundle()
            sage: T_P2.rank()
            2
        """
        return self._filt.dimension()

    def _repr_(self):
        r"""
        Return a string representation.

        OUTPUT:
        
        String.

        EXAMPLES::
        
            sage: toric_varieties.P2().sheaves.tangent_bundle()
            Rank 2 bundle on 2-d CPR-Fano toric variety covered by 3 affine patches.
        """
        s = 'Rank '+str(self.rank())+' bundle on '+str(self._variety)+'.'
        return s
        
    def get_filtration(self, ray=None):
        r"""
        Return the filtration associated to the ``ray``.

        INPUT:

        - ``ray`` -- Integer, a `N`-lattice point, a one-dimensional
          cone, or ``None`` (default). Specifies a ray of the fan of
          the toric variety, either via its index or its generator.

        OUTPUT:

        The filtered vector space associated to the given ``ray``.

        EXAMPLES:
        
            sage: TX = toric_varieties.dP6().sheaves.tangent_bundle()
            sage: TX.get_filtration(0)
            QQ^2 >= QQ^1 >= 0
            sage: TX.get_filtration([-1, -1])
            QQ^2 >= QQ^1 >= 0
            sage: TX.get_filtration(TX.variety().fan(1)[0])
            QQ^2 >= QQ^1 >= 0
            sage: TX.get_filtration()
            Filtrations
                N(-1, -1): QQ^2 >= QQ^1 >= 0
                 N(-1, 0): QQ^2 >= QQ^1 >= 0
                 N(0, -1): QQ^2 >= QQ^1 >= 0
                  N(0, 1): QQ^2 >= QQ^1 >= 0
                  N(1, 0): QQ^2 >= QQ^1 >= 0
                  N(1, 1): QQ^2 >= QQ^1 >= 0
        """
        if ray is None:
            return self._filt
        X = self.variety()
        fan = X.fan()
        if is_Cone(ray):
            if ray.dim() != 1:
                raise ValueError('not a one-dimensional cone')
            ray = ray.ray(0)
        elif ray in ZZ:
            ray = fan.ray(ray)
        else:
            N = fan.lattice()
            ray = N(ray)
            ray.set_immutable()
        return self._filt.get_filtration(ray)

    def get_degree(self, ray, i):
        r"""
        Return the vector subspace ``E^\alpha(i)``.

        EXAMPLES::

            sage: TX = toric_varieties.dP6().sheaves.tangent_bundle()
            sage: TX.get_degree(0, 1)
            Vector space of degree 2 and dimension 1 over Rational Field
            Basis matrix:
            [0 1]
        """
        return self.get_filtration(ray).get_degree(i)

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
            sage: V = X.sheaves.tangent_bundle()
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
            V = V.intersection(self.get_degree(alpha, i))
        return V

    def E_degree(self, alpha, m):
        r"""
        Return the vector subspace `E^\alpha(i)`.

        EXAMPLES::

            sage: X = toric_varieties.P2()
            sage: ray = X.fan(1)[0]
            sage: M = X.fan().dual_lattice()
            sage: V = X.sheaves.tangent_bundle()
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
        return self.get_degree(ray_index, ray*m)

    @cached_method_immutable
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
            sage: V = X.sheaves.tangent_bundle()
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

    @cached_method_immutable
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
            sage: V = X.sheaves.tangent_bundle()
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

    @cached_method_immutable
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
            sage: rays = [(1,0,0), (0,1,0), (0,0,1)]
            sage: F1 = FilteredVectorSpace(rays, {0:[0], 1:[2], 2:[1]})
            sage: F2 = FilteredVectorSpace(3, 0)
            sage: r = P3.fan().rays()
            sage: V = P3.sheaves.Klyachko({r[0]:F1, r[1]:F2, r[2]:F2, r[3]:F2})
            sage: tau = Cone([(1,0,0), (0,1,0)])
            sage: sigma = Cone([(1,0,0)])
            sage: M = P3.fan().dual_lattice()
            sage: m = M(2,1,0)
            sage: V.E_quotient(sigma, m)
            Vector space quotient V/W of dimension 2 over Rational Field where
            V: Vector space of dimension 3 over Rational Field
            W: Vector space of degree 3 and dimension 1 over Rational Field
            Basis matrix:
            [0 1 0]
            sage: V.E_quotient(tau, m)
            Vector space quotient V/W of dimension 2 over Rational Field where
            V: Vector space of dimension 3 over Rational Field
            W: Vector space of degree 3 and dimension 1 over Rational Field
            Basis matrix:
            [0 1 0]
            sage: V.E_quotient_projection(sigma, tau, m)
            Vector space morphism represented by the matrix:
            [1 0]
            [0 1]
            Domain: Vector space quotient V/W of dimension 2 over Rational Field where
            V: Vector space of dimension 3 over Rational Field
            W: Vector space of degree 3 and dimension 1 over Rational Field
            Basis matrix:
            [0 1 0]
            Codomain: Vector space quotient V/W of dimension 2 over Rational Field where
            V: Vector space of dimension 3 over Rational Field
            W: Vector space of degree 3 and dimension 1 over Rational Field
            Basis matrix:
            [0 1 0]
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
            sage: rays = [(1,0,0), (0,1,0), (0,0,1)]
            sage: F1 = FilteredVectorSpace(rays, {0:[0], 1:[2], 2:[1]})
            sage: F2 = FilteredVectorSpace(rays, {0:[1,2], 1:[0]})
            sage: r = P3.fan().rays()
            sage: V = P3.sheaves.Klyachko({r[0]:F1, r[1]:F2, r[2]:F2, r[3]:F2})
            sage: tau = Cone([(1,0,0), (0,1,0)])
            sage: sigma = Cone([(1,0,0)])
            sage: M = P3.fan().dual_lattice()
            sage: m = M(1,1,0); m.set_immutable()
            sage: V.cohomology_complex(m)
            Chain complex with at most 2 nonzero terms over Rational Field

            sage: F = CyclotomicField(3)
            sage: P3 = toric_varieties.P(3).change_ring(F)
            sage: V = P3.sheaves.Klyachko({r[0]:F1, r[1]:F2, r[2]:F2, r[3]:F2})
            sage: V.cohomology_complex(m)
            Chain complex with at most 2 nonzero terms over Cyclotomic
            Field of order 3 and degree 2
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

            sage: V = toric_varieties.P2().sheaves.tangent_bundle()
            sage: V.cohomology(weight=(0,0), dim=True)
            (2, 0, 0)
            sage: for i,j in CartesianProduct(range(-2,3), range(-2,3)):
            ...       HH = V.cohomology(weight=(i,j), dim=True)
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
        from sage.modules.all import FreeModule
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

        EXAMPLES::
            
            sage: X = toric_varieties.P2()
            sage: V1 = X.sheaves.trivial_bundle(1) 
            sage: V2 = X.sheaves.trivial_bundle(2) 
            sage: abs(cmp(V2, V1))
            1
            sage: cmp(V2, V1+V1)
            0

            sage: T_X = X.sheaves.tangent_bundle()
            sage: O_X = X.sheaves.trivial_bundle(1)
            sage: T_X + O_X == O_X + T_X
            False
        """
        c = cmp(type(self), type(other))
        if c!=0: return c
        c = cmp(self.variety(), other.variety())
        if c!=0: return c
        c = cmp(self._filt, other._filt)
        if c!=0: return c
        return 0

    def is_isomorphic(self, other):
        """
        Test whether two bundles are isomorphic.

        EXAMPLES::

            sage: X = toric_varieties.P2()
            sage: T_X = X.sheaves.tangent_bundle()
            sage: O_X = X.sheaves.trivial_bundle(1)
            sage: T_X + O_X == O_X + T_X
            False
            sage: (T_X + O_X).is_isomorphic(O_X + T_X)
            Traceback (most recent call last):
            ...
            NotImplementedError
        """
        raise NotImplementedError

    def direct_sum(self, other):
        """
        Return the sum of two vector bundles.
        
        EXAMPLES::

            sage: X = toric_varieties.P2()
            sage: V1 = X.sheaves.trivial_bundle(1) 
            sage: V2 = X.sheaves.trivial_bundle(2) 
            sage: V2.direct_sum(V1)
            Rank 3 bundle on 2-d CPR-Fano toric variety covered by 3 affine patches.
            
            sage: V1 = X.sheaves.trivial_bundle(1) 
            sage: V2 = X.sheaves.trivial_bundle(2) 
            sage: V2 == V1 + V1
            True
        """
        if not self.variety() == other.variety():
            raise ValueError('The bundles must be over the same base toric variety.')
        filt = self._filt + other._filt
        return self.__class__(self.variety(), filt, check=True)

    __add__ = direct_sum

    def tensor_product(self, other):
        """
        Return the sum of two vector bundles.
        
        EXAMPLES::

            sage: X = toric_varieties.P2()
            sage: X.sheaves.tangent_bundle() + X.sheaves.trivial_bundle(1)
            Rank 3 bundle on 2-d CPR-Fano toric variety covered by 3 affine patches.
            
            sage: OX = X.sheaves.trivial_bundle(1)
            sage: OX == OX * OX
            True
        """
        if not self.variety() == other.variety():
            raise ValueError('The bundles must be over the same base toric variety.')
        filt = self._filt * other._filt
        return self.__class__(self.variety(), filt, check=True)

    __mul__ = tensor_product

    def exterior_power(self, n):
        """
        EXAMPLES:

            sage: X = toric_varieties.P2_123()
            sage: TX = X.sheaves.tangent_bundle()
            sage: antiK = X.sheaves.line_bundle(-X.K())
            sage: TX.exterior_power(2) == antiK
            True
        """
        filt = self._filt.exterior_power(n)
        return self.__class__(self.variety(), filt, check=True)

    wedge = exterior_power
        
    def symmetric_power(self, n):
        filt = self._filt.symmetric_power(n)
        return self.__class__(self.variety(), filt, check=True)

    def dual(self):
        filt = self._filt.dual()
        return self.__class__(self.variety(), filt, check=True)    
        
    def deformation(self, perturbed_rays=None):
        """
        Return a generic torus-equivariant deformation of the bundle.

        EXAMPLES::

           sage: P1 = toric_varieties.P1()
           sage: H = P1.divisor(0)
           sage: V = P1.sheaves.line_bundle(H) + P1.sheaves.line_bundle(-H)
           sage: V.cohomology(dim=True, weight=(0,))
           (1, 0)
           sage: Vtilde = V.deformation()
           sage: Vtilde.cohomology(dim=True, weight=(0,))
           (1, 0)
        """
        return self

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
        return self.__class__(self.variety(), filt, check=True)
        
