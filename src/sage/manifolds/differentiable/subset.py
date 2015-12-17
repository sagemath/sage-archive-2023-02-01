r"""
Open subsets of differentiable manifolds

The class :class:`OpenDifferentiableSubmanifold` implements open subsets of a
differentiable manifold.

AUTHORS:

- Eric Gourgoulhon, Travis Scrimshaw (2015): initial version

REFERENCES:

.. [1] J.M. Lee : *Introduction to Smooth Manifolds*, 2nd ed., Springer
   (New York) (2012); :doi:`10.1007/978-1-4419-9982-5`
.. [2] S. Kobayashi & K. Nomizu : *Foundations of Differential Geometry*,
   vol. 1, Interscience Publishers (New York) (1963)

"""

#*****************************************************************************
#       Copyright (C) 2015 Eric Gourgoulhon <eric.gourgoulhon@obspm.fr>
#       Copyright (C) 2015 Travis Scrimshaw <tscrimsh@umn.edu>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#  as published by the Free Software Foundation; either version 2 of
#  the License, or (at your option) any later version.
#                  http://www.gnu.org/licenses/
#*****************************************************************************

from sage.manifolds.subset import OpenTopologicalSubmanifold
from sage.manifolds.differentiable.manifold import DifferentiableMixin

class OpenDifferentiableSubmanifold(DifferentiableMixin,
                                    OpenTopologicalSubmanifold):
    """
    An open submanifold of a differentiable manifold, which is any open subset
    of the manifold.

    INPUT:

    - ``ambient`` -- differentiable manifold on which the subset is defined
    - ``name`` -- string; name (symbol) given to the subset
    - ``latex_name`` --  (default: ``None``) string; LaTeX symbol to denote the
      subset; if none is provided, it is set to ``name``

    EXAMPLE:

    The unit ball of the Euclidean 2-plane::

        sage: M = Manifold(2, 'M')
        sage: X.<x,y> = M.chart()
        sage: B = M.open_subset('B', coord_def={X: x^2+y^2<1})
        sage: B
        Open subset B of the 2-dimensional differentiable manifold M
        sage: type(B)
        <class 'sage.manifolds.differentiable.subset.OpenDifferentiableSubmanifold_with_category'>
        sage: B.category()
        Join of Category of subobjects of sets and Category of smooth manifolds
         over Real Field with 53 bits of precision

    An open subset of a differentiable manifold being a differentiable manifold
    by itself, ``B`` has all the attributes of a manifold::

        sage: dim(B)
        2
        sage: B.base_field()
        Real Field with 53 bits of precision

    Given the definition of ``B``, it is automatically endowed with a chart,
    which is the restriction of ``X`` to ``B``::

        sage: B.atlas()
        [Chart (B, (x, y))]
        sage: B.default_chart()
        Chart (B, (x, y))
        sage: B.default_chart() is X.restrict(B)
        True

    An point in ``B``::

        sage: p = B.an_element(); p
        Point on the 2-dimensional differentiable manifold M
        sage: X(p)  # the coordinates (x,y) of p
        (0, 0)
        sage: p in B
        True

    Checking whether various points, defined by their coordinates w.r.t.
    chart ``X``,  are in ``B``::

        sage: M((0,1/2)) in B
        True
        sage: M((0,1)) in B
        False
        sage: M((1/2,1)) in B
        False
        sage: M((-1/2,1/3)) in B
        True

    """
    def __init__(self, ambient, name, latex_name=None):
        """
        Initialize ``self``.

        TESTS::

            sage: M = Manifold(2, 'M')
            sage: X.<x,y> = M.chart()
            sage: A = M.open_subset('A', coord_def={X: x^2+y^2<1}); A
            Open subset A of the 2-dimensional differentiable manifold M
            sage: type(A)
            <class 'sage.manifolds.differentiable.subset.OpenDifferentiableSubmanifold_with_category'>
            sage: A.category()
            Join of Category of subobjects of sets and Category of smooth
             manifolds over Real Field with 53 bits of precision
            sage: TestSuite(A).run(skip='_test_not_implemented_methods')

        .. TODO::

            implement method ``lift`` so that ``_test_not_implemented_methods``
            is passed.

        """
        DifferentiableMixin.__init__(self, diff_degree=ambient.diff_degree())
        category = ambient.category().Subobjects()
        OpenTopologicalSubmanifold.__init__(self, ambient, name,
                                            latex_name=latex_name)
