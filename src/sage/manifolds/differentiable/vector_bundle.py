r"""
Differentiable Vector Bundle

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
from sage.manifolds.vector_bundle import TopologicalVectorBundle

class DifferentiableVectorBundle(TopologicalVectorBundle):
    r"""

    """
    def __init__(self, rank, name, base_space, field='real', latex_name=None,
                 category=None, unique_tag=None):
        r"""

        """
        if category is None:
            if field == 'real':
                field_c = RR
            elif field == 'complex':
                field_c = CC
            else:
                field_c = field
            if diff_degree == infinity:
                category = VectorBundles(field_c).Smooth()
            else:
                category = VectorBundles(field_c).Differentiable()
        TopologicalVectorBundle.__init__(rank, name, base_space, field=field,
                                         latex_name=latex_name,
                                         category=category)
        self._diff_degree = base_space._diff_degree

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
                                            diff_degree=manifold._diff_degree,
                                            latex_name=latex_name)
        self._tensor_type = (k, l)
        # if manifold._atlas not empty, introduce trivializations

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
