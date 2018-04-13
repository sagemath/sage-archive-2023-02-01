r"""
Immersion of a differentiable manifold

Given a differentiable manifold M over the topological field K, an
immersion of M is defined by a differentiable manifold N on the same
field K of dimension smaller, and a differentiable map \phi from N to M
(also called immersion).

In addition, if the differential of \phi is everywhere injective,
then it's called an embedding. \phi can also depend on one or
multiple parameters. As long as the differential of \phi remains
injective in these parameters, it represent a foliation. The dimension
of the foliation is defined as the numbers of parameters.

AUTHORS:

- Florentin Jaffredo

"""

# *****************************************************************************
#  Copyright (C) 2018 Florentin Jaffredo <florentin.jaffredo@polytechnique.edu>
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 2 of the License, or
# (at your option) any later version.
#                  http://www.gnu.org/licenses/
# *****************************************************************************

from sage.manifolds.differentiable.manifold import DifferentiableManifold
from sage.manifolds.submanifold.topological_immersion import \
    TopologicalImmersion
from sage.rings.infinity import infinity


class Immersion(DifferentiableManifold, TopologicalImmersion):
    r"""
    Immersion of a differentiable manifold.

    Given a differentiable manifold M over the topological field K, an
    immersion of M is defined by a differentiable manifold N on the same
    field K of dimension smaller, and a differentiable map \phi from N to M
    (also called immersion).

    In addition, if the differential of \phi is everywhere injective,
    then it's called an embedding. \phi can also depend on one or
    multiple parameters. As long as the differential of \phi remains
    injective in these parameters, it represent a foliation. The dimension
    of the foliation is defined as the numbers of parameters.

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

        sage: M = Manifold(3, 'M')
        sage: N = Manifold(2, 'N', ambient = M)
        sage: N
        2-dimensional differentiable submanifold N embedded in 3-dimensional
         differentiable manifold M
        sage: CM.<x,y,z> = M.chart()
        sage: CN.<u,v> = N.chart()

    Let's define a 1-dimension foliation indexed by t. The inverse map is needed
    in order to compute the adapted chart in the ambient manifold::

        sage: t = var('t')
        sage: phi = N.diff_map(M,{(CN, CM):[u, v, t+u**2+v**2]}); phi
        Differentiable map from the 2-dimensional differentiable submanifold N
         embedded in 3-dimensional differentiable manifold M to the
         3-dimensional differentiable manifold M
        sage: phi_inv = M.diff_map(N, {(CM, CN):[x, y]})
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
       """
    def __init__(self, n, name, field, structure, ambient=None,
                 base_manifold=None, diff_degree=infinity,
                 latex_name=None, start_index=0, category=None,
                 unique_tag=None):
        r"""
        Construct an immersion of a given differentiable manifold.

        EXAMPLES::

            sage: M = Manifold(3,'M')
            sage: N = Manifold(2,'N',ambient = M)
            sage: N
            2-dimensional differentiable submanifold N embedded in 3-dimensional
             differentiable manifold M

        """
        DifferentiableManifold.__init__(self, n, name, field, structure,
                                        base_manifold=base_manifold,
                                        diff_degree=diff_degree,
                                        latex_name=latex_name,
                                        start_index=start_index,
                                        category=category)
        TopologicalImmersion.__init__(self, n, name, field, structure,
                                      ambient=ambient,
                                      base_manifold=base_manifold,
                                      category=category)

    def _repr_(self):
        r"""
        Return a string representation of the submanifold. If no ambient
        manifold is specified, the submanifold is considered as a manifold

        TESTS::

            sage: M = Manifold(3,'M')
            sage: N = Manifold(2,'N',ambient = M)
            sage: N._repr_()
            '2-dimensional differentiable submanifold N embedded in
             3-dimensional differentiable manifold M'

        """
        if self._ambient is None:
            return super(DifferentiableManifold, self).__repr__()
        return "{}-dimensional differentiable submanifold {} embedded in {}-" \
               "dimensional differentiable " \
               "manifold {}".format(self._dim, self._name, self._ambient._dim,
                                    self._ambient._name)
