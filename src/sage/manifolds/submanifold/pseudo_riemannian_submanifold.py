r"""
Pseudo-Riemannian submanifold of a differentiable manifold

A pseudo-Riemannian submanifold of a differentiable manifold is a differentiable
submanifold which is also pseudo-Riemannian.

AUTHORS:

- Florentin Jaffredo

"""

# *****************************************************************************
#   Copyright (C) 2018 Florentin Jaffredo <florentin.jaffredo@polytechnique.edu>
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 2 of the License, or
# (at your option) any later version.
#                  http://www.gnu.org/licenses/
# *****************************************************************************

from sage.manifolds.differentiable.pseudo_riemannian import \
    PseudoRiemannianManifold
from sage.manifolds.submanifold.differentiable_submanifold import \
    DifferentiableSubmanifold
from sage.rings.infinity import infinity


class PseudoRiemannianSubmanifold(PseudoRiemannianManifold,
                                  DifferentiableSubmanifold):
    r"""
    Pseudo-Riemannian submanifold of a differentiable manifold

    A pseudo-Riemannian submanifold of a differentiable manifold is a
    differentiable submanifold which is also pseudo-Riemannian.

    INPUT:

    - ``n`` -- positive integer; dimension of the manifold
    - ``name`` -- string; name (symbol) given to the manifold
    - ``field`` -- field `K` on which the manifold is
      defined; allowed values are

        - ``'real'`` or an object of type ``RealField`` (e.g., ``RR``) for
           a manifold over `\RR`
        - ``'complex'`` or an object of type ``ComplexField`` (e.g., ``CC``)
           for a manifold over `\CC`
        - an object in the category of topological fields (see
          :class:`~sage.categories.fields.Fields` and
          :class:`~sage.categories.topological_spaces.TopologicalSpaces`)
          for other types of manifolds

    - ``structure`` -- manifold structure (see
      :class:`~sage.manifolds.structure.TopologicalStructure` or
      :class:`~sage.manifolds.structure.RealTopologicalStructure`)
    - ``ambient`` -- (default: ``None``) manifold of destination
      of the immersion. If ``None``, set to ``self``
    - ``base_manifold`` -- (default: ``None``) if not ``None``, must be a
      topological manifold; the created object is then an open subset of
      ``base_manifold``
    - ``latex_name`` -- (default: ``None``) string; LaTeX symbol to
      denote the manifold; if none are provided, it is set to ``name``
    - ``start_index`` -- (default: 0) integer; lower value of the range of
      indices used for "indexed objects" on the manifold, e.g., coordinates
      in a chart
      - ``category`` -- (default: ``None``) to specify the category; if
      ``None``, ``Manifolds(field)`` is assumed (see the category
      :class:`~sage.categories.manifolds.Manifolds`)
    - ``unique_tag`` -- (default: ``None``) tag used to force the construction
      of a new object when all the other arguments have been used previously
      (without ``unique_tag``, the
      :class:`~sage.structure.unique_representation.UniqueRepresentation`
      behavior inherited from
      :class:`~sage.manifolds.subset.ManifoldSubset`
      would return the previously constructed object corresponding to these
      arguments)

    EXAMPLES:

    Let N be a 2-dimensional submanifold of M, 3-dimensional manifold::

        sage: M = Manifold(3, 'M', structure ="pseudo-Riemannian")
        sage: N = Manifold(2, 'N', ambient = M, structure ="pseudo-Riemannian")
        sage: N
        2-dimensional pseudo-Riemannian submanifold N embedded in 3-dimensional
         differentiable manifold M
        sage: CM.<x,y,z> = M.chart()
        sage: CN.<u,v> = N.chart()

    Let's define a 1-dimension foliation indexed by t. The inverse map is needed
    in order to compute the adapted chart in the ambient manifold::

        sage: t = var('t')
        sage: phi = N.diff_map(M, {(CN,CM):[u, v, t+u**2+v**2]}); phi
        Differentiable map from the 2-dimensional pseudo-Riemannian submanifold
         N embedded in 3-dimensional differentiable manifold M to the
         3-dimensional Riemannian manifold M
        sage: phi_inv = M.diff_map(N,{(CM, CN): [x,y]})
        sage: phi_inv_t = M.scalar_field({CM: z-x**2-y**2})

    \phi can then be declared as an embedding from N to M::

        sage: N.set_immersion(phi, phi_inverse = phi_inv, var = t,\
        ....:                 t_inverse = {t: phi_inv_t})
        sage: N.declare_embedding()

    The foliation can also be used to find new charts on the ambient manifold
    that are adapted to the foliation, ie in which the expression of the
    immersion is trivial. At the same time coordinates changes or computed::

        sage: N.adapted_chart()
        [Chart (M, (u_M, v_M, t_M))]
        sage: len(M._coord_changes)
        2

    .. SEEALSO::

        :mod:`sage.manifolds.manifold`
        :mod:`sage.manifolds.submanifold.differentiable_submanifold`
   """
    def __init__(self, n, name, ambient=None, metric_name='g', signature=None,
                 base_manifold=None, diff_degree=infinity, latex_name=None,
                 metric_latex_name=None, start_index=0, category=None,
                 unique_tag=None):
        r"""
        Construct an immersion of a given differentiable manifold.

        EXAMPLES::

            sage: M = Manifold(3, 'M', structure="pseudo-Riemannian")
            sage: N = Manifold(2, 'N', ambient=M, structure="pseudo-Riemannian")
            sage: N
            2-dimensional pseudo-Riemannian submanifold N embedded in
             3-dimensional differentiable manifold M

        """

        PseudoRiemannianManifold.__init__(self, n, name=name,
                                          metric_name=metric_name,
                                          signature=signature,
                                          base_manifold=base_manifold,
                                          diff_degree = diff_degree,
                                          latex_name=latex_name,
                                          metric_latex_name=metric_latex_name,
                                          start_index=start_index,
                                          category=category)
        DifferentiableSubmanifold.__init__(self, n, name, self._field, self._structure,
                                           ambient = ambient, base_manifold=base_manifold,
                                           latex_name=latex_name, start_index=start_index,
                                           category=category)

        self._difft = None          #done
        self._gradt = None          #done
        self._normal = None         #done
        self._lapse = None          #done
        self._shift = None          #done
        self._gamma = None          #done
        self._ambient_gamma = None  #done
        self._K = None              #done
        self._ambient_K = None      #done
        self._ambient_g = None      #done

        self._sgn = 1 if ambient._structure.name == "Riemannian" else -1


    def _repr_(self):
        r"""
        Return a string representation of the submanifold. If no ambient
        manifold is specified, the submanifold is considered as a manifold

        TESTS::

            sage: M = Manifold(3, 'M', structure="pseudo-Riemannian")
            sage: N = Manifold(2, 'N', ambient=M, structure="pseudo-Riemannian")
            sage: N._repr_()
            '2-dimensional pseudo-Riemannian submanifold N embedded in
             3-dimensional differentiable manifold M'

        """
        if self._ambient == None:
            return super(PseudoRiemannianManifold,self).__repr__()
        return "{}-dimensional pseudo-Riemannian submanifold {} embedded " \
               "in {}-dimensional differentiable " \
               "manifold {}".format(self._dim, self._name, self._ambient._dim,
                                    self._ambient._name)

    def ambient_g(self):
        if  not self._embedded or not isinstance(self._ambient,PseudoRiemannianManifold):
            raise ValueError("Submanifold must be embedded in a pseudo-Riemnnian manifold")
        if self._ambient_g is not None:
            return self._ambient_g
        self._ambient_g = self._ambient.metric()
        return self._ambient_g

    def gamma(self):
        if self._gamma is not None:
            return self._gamma
        #self._gamma = PseudoRiemannianManifold.metric(self,r'\gamma_' + self._name)
        self._gamma = self.metric()
        self._gamma.set(self._immersion.pullback(self.ambient_g()))
        return self._gamma

    def difft(self):
        if self._dim_foliation == 0:
            raise ValueError("A foliation is needed to perform this calculation")
        if self._difft is not None:
            return self._difft
        self._difft = self._t_inverse[self._var[0]].differential()
        return self._difft

    def gradt(self):
        if self._dim_foliation == 0:
            raise ValueError("A foliation is needed to perform this calculation")
        if self._gradt is not None:
            return self._gradt
        self._gradt = self._ambient_g.inverse().contract(self.difft())
        return self._gradt

    def normal(self):
        if self._normal is not None:
            return self._normal
        self._normal = self._sgn*self.lapse()*self.gradt()

        # ne marche pas :
        # product = self.ambient_g().contract(self._immersion.pushforward(self.atlas()[0].frame()[0]))
        # for i in range(1,self._dim):
        #     product = product.wedge(self.ambient_g().contract(self._immersion.pushforward(self.atlas()[0].frame()[i])))
        # self._normal = product.hodge_dual(self.ambient_g())
        return self._normal

    def ambient_gamma(self):
        if self._ambient_gamma is not None:
            return self._ambient_gamma
        self._ambient_gamma = self.ambient_g() + \
                              self.ambient_g().contract(self.normal()) * \
                              self.ambient_g().contract(self.normal())
        return self._ambient_gamma

    def lapse(self):
        if self._dim_foliation == 0:
            raise ValueError("A foliation is needed to perform this calculation")
        if self._lapse is not None:
            return self._lapse
        self._lapse = 1/(self._sgn*self.ambient_g()(self.gradt(),self.gradt())).sqrt()
        return self._lapse

    def shift(self):
        if self._dim_foliation == 0:
            raise ValueError("A foliation is needed to perform this calculation")
        if self._shift is not None:
            return self._shift
        self._shift = self._adapted_charts[0].frame()[self._dim]-self.lapse()*self.normal()
        return self._shift

    def ambient_K(self):
        if self._ambient_K is not None:
            return self._ambient_K
        nab = self.ambient_g().connection('nabla', r'\nabla')
        self._ambient_K = -self.ambient_g().contract(nab(self.normal()))-\
                          nab(self.normal()).contract(self.normal()).contract(self.ambient_g())*\
                          self.normal().contract(self.ambient_g())
        return self._ambient_K

    def K(self):
        if self._K is not None:
            return self._K
        self._K = self._immersion.pullback(self.ambient_K())
        return self._K




