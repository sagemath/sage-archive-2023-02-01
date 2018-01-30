r"""
Pseudo-Riemannian Manifolds

A *pseudo-Riemannian manifold* is a pair `(M,g)` where `M` is a real
differentiable manifold `M` (see
:class:`~sage.manifolds.differentiable.manifold.DifferentiableManifold`)
and `g` is a field non-degenerate symmetric bilinear forms on `M`, which is
called the *metric tensor*, or simply the *metric* (see
:class:`~sage.manifolds.differentiable.metric.PseudoRiemannianMetric`).

Two important subcases are

- *Riemannian manifold*: the metric `g` is definite-positive
- *Lorentzian manifold*: the metric `g` has signature `n-2` (positive
  convention) or `2-n` (negative convention), where `n = \dim M`.

All pseudo-Riemannian manifolds are implemented via the class
:class:`PseudoRiemannianManifold`.

AUTHORS:

- Eric Gourgoulhon (2018): initial version

"""

#*****************************************************************************
#       Copyright (C) 2018 Eric Gourgoulhon <eric.gourgoulhon@obspm.fr>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#  as published by the Free Software Foundation; either version 2 of
#  the License, or (at your option) any later version.
#                  http://www.gnu.org/licenses/
#*****************************************************************************

from sage.rings.infinity import infinity
from sage.rings.integer import Integer
from sage.manifolds.structure import (PseudoRiemannianStructure,
                                      RiemannianStructure, LorentzianStructure)
from sage.manifolds.differentiable.manifold import DifferentiableManifold
from sage.manifolds.differentiable.metric import PseudoRiemannianMetric
from sage.manifolds.differentiable.tensorfield import TensorField

###############################################################################

class PseudoRiemannianManifold(DifferentiableManifold):
    r"""
    PseudoRiemannian manifold.

    A *pseudo-Riemannian manifold* is a pair `(M,g)` where `M` is a real
    differentiable manifold `M` (see
    :class:`~sage.manifolds.differentiable.manifold.DifferentiableManifold`)
    and `g` is a field non-degenerate symmetric bilinear forms on `M`, which is
    called the *metric tensor*, or simply the *metric* (see
    :class:`~sage.manifolds.differentiable.metric.PseudoRiemannianMetric`).

    Two important subcases are

    - *Riemannian manifold*: the metric `g` is definite-positive
    - *Lorentzian manifold*: the metric `g` has signature `n-2` (positive
      convention) or `2-n` (negative convention), where `n = \dim M`.

    INPUT:

    - ``n`` -- positive integer; dimension of the manifold
    - ``name`` -- string; name (symbol) given to the manifold
    - ``metric_name`` -- (default: ``'g'``) string; name (symbol) given to the
      metric
    - ``signature`` -- (default: ``None``) signature `S` of the metric as a
      single integer: `S = n_+ - n_-`, where `n_+` (resp. `n_-`) is the
      number of positive terms (resp. number of negative terms) in any
      diagonal writing of the metric components; if ``signature`` is not
      provided, `S` is set to the manifold's dimension (Riemannian
      signature)
    - ``ambient`` -- (default: ``None``) if not ``None``, must be a
      differentiable manifold; the created object is then an open subset of
      ``ambient``
    - ``diff_degree`` -- (default: ``infinity``) degree `k` of
      differentiability
    - ``latex_name`` -- (default: ``None``) string; LaTeX symbol to
      denote the manifold; if none is provided, it is set to ``name``
    - ``metric_latex_name`` -- (default: ``None``) string; LaTeX symbol to
      denote the metric; if none is provided, it is set to ``metric_name``
    - ``start_index`` -- (default: 0) integer; lower value of the range of
      indices used for "indexed objects" on the manifold, e.g. coordinates
      in a chart
    - ``category`` -- (default: ``None``) to specify the category; if ``None``,
      ``Manifolds(RR).Differentiable()`` (or ``Manifolds(RR).Smooth()``
      if ``diff_degree`` = ``infinity``) is assumed (see the category
      :class:`~sage.categories.manifolds.Manifolds`)
    - ``unique_tag`` -- (default: ``None``) tag used to force the construction
      of a new object when all the other arguments have been used previously
      (without ``unique_tag``, the
      :class:`~sage.structure.unique_representation.UniqueRepresentation`
      behavior inherited from
      :class:`~sage.manifolds.subset.ManifoldSubset`, via
      :class:`~sage.manifolds.differentiable.manifold.DifferentiableManifold`
      and :class:`~sage.manifolds.manifold.TopologicalManifold`,
      would return the previously constructed object corresponding to these
      arguments).

    """
    def __init__(self, n, name, metric_name='g', signature=None, ambient=None,
                 diff_degree=infinity, latex_name=None,
                 metric_latex_name=None, start_index=0, category=None,
                 unique_tag=None):
        r"""
        Construct a differentiable manifold.
        """
        if ambient and not isinstance(ambient, PseudoRiemannianManifold):
            raise TypeError("the argument 'ambient' must be a " +
                            "pseudo-Riemannian manifold")
        if signature is None or signature == n:
            structure = RiemannianStructure()
        elif signature == n-2 or signature == 2-n:
            structure = LorentzianStructure()
        else:
            structure = PseudoRiemannianStructure()
        DifferentiableManifold.__init__(self, n, name, 'real', structure,
                                        ambient=ambient,
                                        diff_degree=diff_degree,
                                        latex_name=latex_name,
                                        start_index=start_index,
                                        category=category)
        if ambient:
            self._metric = ambient.metric().restrict(self)
            self._metric_signature = ambient._metric_signature
        else:
            self._metric = None # to be initialized by set_metric() or by metric()
            self._metric_signature = signature
        if not isinstance(metric_name, str):
            raise TypeError("{} is not a string".format(metric_name))
        self._metric_name = metric_name
        if metric_latex_name is None:
            self._metric_latex_name = self._metric_name
        else:
            if not isinstance(metric_latex_name, str):
                raise TypeError("{} is not a string".format(metric_latex_name))
            self._metric_latex_name = metric_latex_name

    def set_metric(self, metric):
        r"""
        Set the metric on ``self``.

        """
        if isinstance(metric, PseudoRiemannianMetric):
            if metric._name != self._metric_name:
                raise ValueError("name of {} does not match ".format(metric) +
                                 "the name '{}' ".format(self._metric_name) +
                                 "declared at the construction of the " +
                                 "{}".format(self))
            if metric.parent()._dest_map is not self.identity_map():
                raise ValueError("{} is not a metric ".format(metric) +
                                 "defined on {}".format(self))
            self._metric = metric
            metric._latex_name = self._metric_latex_name
        elif isinstance(metric, TensorField):
            self._metric = self.metric(self._metric_name,
                                       signature=self._metric_signature,
                                       latex_name=self._metric_latex_name)
            self._metric.set(metric)
        else:
            raise TypeError("{} must be a metric or a ".format(metric) +
                            "of bilinear forms")

    def metric(self, name=None, signature=None, latex_name=None,
               dest_map=None):
        r"""
        Return the metric of the pseudo-Riemannian manifold ``self`` or
        defines a new metric tensor on ``self``.

        INPUT:

        - ``name`` -- (default: ``None``) name given to the metric; if
          ``name`` is ``None`` or matches the name of the metric defining the
          pseudo-Riemannian structure of ``self``, the latter metric is
          returned
        - ``signature`` -- (default: ``None``; ignored if ``name`` is ``None``)
          signature `S` of the metric as a single integer: `S = n_+ - n_-`,
          where `n_+` (resp. `n_-`) is the number of positive terms (resp.
          number of negative terms) in any diagonal writing of the metric
          components; if ``signature`` is not provided, `S` is set to the
          manifold's dimension (Riemannian signature)
        - ``latex_name`` -- (default: ``None``; ignored if ``name`` is ``None``)
          LaTeX symbol to denote the metric; if ``None``, it is formed from
          ``name``
        - ``dest_map`` -- (default: ``None``; ignored if ``name`` is ``None``)
          instance of
          class :class:`~sage.manifolds.differentiable.diff_map.DiffMap`
          representing the destination map `\Phi:\ U \rightarrow M`, where `U`
          is the current manifold; if ``None``, the identity map is assumed
          (case of a metric tensor field *on* `U`)

        OUTPUT:

        - instance of
          :class:`~sage.manifolds.differentiable.metric.PseudoRiemannianMetric`

        EXAMPLES:

        """
        if name is None or name == self._metric_name:
            if self._metric is None:
                self._metric = DifferentiableManifold.metric(self,
                                           self._metric_name,
                                           signature=self._metric_signature,
                                           latex_name=self._metric_latex_name)
            return self._metric
        return DifferentiableManifold.metric(self, name, signature=signature,
                                             latex_name=latex_name,
                                             dest_map=dest_map)
