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

.. RUBRIC:: Example: the sphere as a Riemannian manifold of dimension 2

One starts by declaring `S^2` as a 2-dimensional Riemnnian manifold::

    sage: M = Manifold(2, 'S^2', structure='Riemannian')
    sage: M
    2-dimensional Riemannian manifold S^2
    sage: U = M.open_subset('U')
    sage: stereoN.<x,y> = U.chart()
    sage: V = M.open_subset('V')
    sage: stereoS.<u,v> = V.chart()
    sage: stereoN_to_S = stereoN.transition_map(stereoS,
    ....:                [x/(x^2+y^2), y/(x^2+y^2)], intersection_name='W',
    ....:                restrictions1= x^2+y^2!=0, restrictions2= u^2+v^2!=0)
    sage: stereoN_to_S
    Change of coordinates from Chart (W, (x, y)) to Chart (W, (u, v))
    sage: stereoN_to_S.display()
    u = x/(x^2 + y^2)
    v = y/(x^2 + y^2)
    sage: stereoN_to_S.inverse().display()
    x = u/(u^2 + v^2)
    y = v/(u^2 + v^2)
    sage: M.declare_union(U,V)
    sage: W = U.intersection(V)
    sage: g = M.metric()
    sage: g[stereoN.frame(), 0, 0] = 4/(1 + x^2 + y^2)^2
    sage: g[stereoN.frame(), 1, 1] = 4/(1 + x^2 + y^2)^2
    sage: g.add_comp_by_continuation(stereoS.frame(), W)
    sage: gV = V.metric()
    sage: gV.display()
    g = 4/(u^4 + v^4 + 2*(u^2 + 1)*v^2 + 2*u^2 + 1) du*du + 4/(u^4 + v^4 + 2*(u^2 + 1)*v^2 + 2*u^2 + 1) dv*dv


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
        Construct a pseudo-Riemannian manifold.

        TESTS::

            sage: M = Manifold(4, 'M', structure='pseudo-Riemannian',
            ....:              signature=0)
            sage: M
            4-dimensional pseudo-Riemannian manifold M
            sage: type(M)
            <class 'sage.manifolds.differentiable.pseudo_riemannian.PseudoRiemannianManifold_with_category'>
            sage: X.<w,x,y,z> = M.chart()
            sage: M.metric()
            Pseudo-Riemannian metric g on the 4-dimensional pseudo-Riemannian manifold M
            sage: TestSuite(M).run()

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
            # Default metric associated with the manifold
            if self._metric is None:
                if self._manifold is not self and self._manifold._metric is not None:
                    # case of an open subset with a metric already defined on
                    # the ambient manifold:
                    self._metric = self._manifold._metric.restrict(self)
                else:
                    # creation from scratch:
                    self._metric = DifferentiableManifold.metric(self,
                                           self._metric_name,
                                           signature=self._metric_signature,
                                           latex_name=self._metric_latex_name)
            return self._metric
        # Metric distinct from the default one: it is created by the method
        # metric of the superclass for generic differentiable manifolds:
        return DifferentiableManifold.metric(self, name, signature=signature,
                                             latex_name=latex_name,
                                             dest_map=dest_map)

    def open_subset(self, name, latex_name=None, coord_def={}):
        r"""
        Create an open subset of ``self``.

        An open subset is a set that is (i) included in the manifold and (ii)
        open with respect to the manifold's topology. It is a differentiable
        manifold by itself. Moreover, equipped with the restriction of the
        manifold metric to itself, it is a pseudo-Riemannian manifold. Hence
        the returned object is an instance of
        :class:`PseudoRiemannianManifold`.

        INPUT:

        - ``name`` -- name given to the open subset
        - ``latex_name`` --  (default: ``None``) LaTeX symbol to denote the
          subset; if none is provided, it is set to ``name``
        - ``coord_def`` -- (default: {}) definition of the subset in
          terms of coordinates; ``coord_def`` must a be dictionary with keys
          charts in the manifold's atlas and values the symbolic expressions
          formed by the coordinates to define the subset.

        OUTPUT:

        - instance of :class:`PseudoRiemannianManifold` representing the
          created open subset

        EXAMPLES:

        Open subset of a 2-dimensional Riemannian manifold::

            sage: M = Manifold(2, 'M', structure='Riemannian')
            sage: X.<x,y> = M.chart()
            sage: U = M.open_subset('U', coord_def={X: x>0}); U
            Open subset U of the 2-dimensional Riemannian manifold M
            sage: type(U)
            <class 'sage.manifolds.differentiable.pseudo_riemannian.PseudoRiemannianManifold_with_category'>

        We initialize the metric of ``M``::

            sage: g = M.metric()
            sage: g[0,0], g[1,1] = 1, 1

        Then the metric on ``U`` is determined as the restriction of ``g`` to
        ``U``::

            sage: gU = U.metric(); gU
            Riemannian metric g on the Open subset U of the 2-dimensional Riemannian manifold M
            sage: gU.display()
            g = dx*dx + dy*dy
            sage: gU is g.restrict(U)
            True

        TESTS:

        Open subset created after the initialization of the metric::

            sage: V = M.open_subset('V', coord_def={X: x<0}); V
            Open subset V of the 2-dimensional Riemannian manifold M
            sage: gV = V.metric()
            sage: gV.display()
            g = dx*dx + dy*dy
            sage: gV is g.restrict(V)
            True

        """
        resu = PseudoRiemannianManifold(self._dim, name,
                                        metric_name=self._metric_name,
                                        signature=self._metric_signature,
                                        ambient=self._manifold,
                                        diff_degree=self._diff_degree,
                                        latex_name=latex_name,
                                        metric_latex_name=self._metric_latex_name,
                                        start_index=self._sindex)
        resu._calculus_method = self._calculus_method
        resu._supersets.update(self._supersets)
        for sd in self._supersets:
            sd._subsets.add(resu)
        self._top_subsets.add(resu)
        # Charts on the result from the coordinate definition:
        for chart, restrictions in coord_def.items():
            if chart not in self._atlas:
                raise ValueError("the {} does not belong to ".format(chart) +
                                 "the atlas of {}".format(self))
            chart.restrict(resu, restrictions)
        # Transition maps on the result inferred from those of self:
        for chart1 in coord_def:
            for chart2 in coord_def:
                if chart2 != chart1 and (chart1, chart2) in self._coord_changes:
                    self._coord_changes[(chart1, chart2)].restrict(resu)
        #!# update non-coordinate vector frames and change of frames
        #
        return resu
