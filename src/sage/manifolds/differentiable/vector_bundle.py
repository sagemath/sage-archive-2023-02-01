r"""
Differentiable Vector Bundles

Let `K` be a topological field. A `C^k`-differentiable *vector bundle* of rank
`n` over the field `K` and over a `C^k`-differentiable manifold `B` (base space)
is a `C^k`-differentiable manifold `E` (total space) together with a
`C^k` differentiable and surjective projection `\pi: E \to B` such that for
every point `x \in B`
- the set `E_x=\pi^{-1}(x)` has the vector space structure of `K^n`,
- there is a neighborhood `U \subset B` of `x` and a `C^k`-diffeomorphism
  `\varphi: \pi^{-1}(x) \to U \times K^n` such that
  `v \mapsto \varphi^{-1}(y,v)` is a linear isomorphism for any `y \in U`.

AUTHORS:

- Michael Jung (2019) : initial version

"""

#******************************************************************************
#       Copyright (C) 2019 Michael Jung <micjung at uni-potsdam.de>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#  as published by the Free Software Foundation; either version 2 of
#  the License, or (at your option) any later version.
#                  https://www.gnu.org/licenses/
#******************************************************************************

from sage.categories.vector_bundles import VectorBundles
from sage.rings.all import CC
from sage.rings.real_mpfr import RR
from sage.rings.infinity import infinity, minus_infinity
from sage.manifolds.vector_bundle import TopologicalVectorBundle

class DifferentiableVectorBundle(TopologicalVectorBundle):
    r"""
    An instance of this class represents a differentiable vector bundle
    `E \to M`

    INPUT:

    - ``rank`` -- positive integer; rank of the vector bundle
    - ``name```-- string representation given to the total space
    - ``base_space`` -- the base space (differentiable manifold) over which the
      vector bundle is defined
    - ``field`` -- field `K` which gives the fibers the structure of a
     vector space over `K`; allowed values are

      - ``'real'`` or an object of type ``RealField`` (e.g., ``RR``) for
        a vector bundle over `\RR`
      - ``'complex'`` or an object of type ``ComplexField`` (e.g., ``CC``)
        for a vector bundle over `\CC`
      - an object in the category of topological fields (see
        :class:`~sage.categories.fields.Fields` and
        :class:`~sage.categories.topological_spaces.TopologicalSpaces`)
        for other types of topological fields

    - ``latex_name`` -- (default: ``None``) latex representation given to the
      total space
    - ``category`` -- (default: ``None``) to specify the category; if
      ``None``, ``VectorBundles(base_space, c_field).Differentiable()`` (or
      ``Manifolds(field).Smooth()`` if ``diff_degree`` = ``infinity``)`` is
      assumed (see the category
      :class:`~sage.categories.vector_bundles.VectorBundles`)

    EXAMPLES:

    """
    def __init__(self, rank, name, base_space, field='real', latex_name=None,
                 category=None, unique_tag=None):
        r"""
        Construct a differentiable vector bundle over some base space.

        TESTS::



        """
        diff_degree = base_space._diff_degree
        if category is None:
            if field == 'real':
                field_c = RR
            elif field == 'complex':
                field_c = CC
            else:
                field_c = field
            if diff_degree == infinity:
                category = VectorBundles(base_space, field_c).Smooth()
            else:
                category = VectorBundles(base_space, field_c).Differentiable()
        TopologicalVectorBundle.__init__(self, rank, name, base_space,
                                         field=field,
                                         latex_name=latex_name,
                                         category=category)
        self._diff_degree = diff_degree # Override diff degree

    def _repr_(self):
        r"""

        """
        desc = "Differentiable "
        return desc + TopologicalVectorBundle._repr_object_name(self)

# *****************************************************************************

class TensorBundle(DifferentiableVectorBundle):
    r"""

    """
    def __init__(self, manifold, k, l):
        rank = manifold.dim()**(k+l)
        if (k, l) == (1, 0):
            name = "T{}".format(manifold._name)
            latex_name = r'T{}'.format(manifold._latex_name)
        elif (k, l) == (0, 1):
            name = "T*{}".format(manifold._name)
            latex_name = r'T^*{}'.format(manifold._latex_name)
        else:
            name = "T^({}{}){}".format(k, l, manifold._name)
            latex_name = r'T^{{{({}{})}}}{}'.format(k, l, manifold._latex_name)
        DifferentiableVectorBundle.__init__(self, rank, name, manifold,
                                            latex_name=latex_name)
        self._tensor_type = (k, l)
        # TODO: if manifold._atlas not empty, introduce trivializations

    def _repr_(self):
        r"""
        String representation of ``self``.

        TESTS::



        """
        if self._tensor_type == (1,0):
            desc = "Tangent bundle "
        elif self._tensor_type == (0,1):
            desc = "Cotangent bundle "
        else:
            desc = "Tensor bundle "
        desc += self._name + " "
        desc += "over {}".format(self._base_space._name)
        return desc

    def fiber(self, point):
        r"""

        """
        return self._base_space.tangent_space(point).tensor_module(self._tensor_type)
