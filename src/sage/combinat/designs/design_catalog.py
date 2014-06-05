r"""
Catalog of designs

This module gathers all designs that can be reached from a Sage sessions
through the ``designs`` objects. In order to create the Witt design on 24 points
it is sufficient to type::

    sage: designs.WittDesign(24) # optional - gap_packages
    Incidence structure with 24 points and 759 blocks

Or a Steiner Triple System on 19 points::

    sage: designs.steiner_triple_system(19)
    Incidence structure with 19 points and 57 blocks

Online database -- Jolla Covering Repository (LJCR) :

There exists an online database of the best known covering designs, the La Jolla
Covering Repository (LJCR), available at [1]_. As it is over 60MB and changes
frequently that database is not included in Sage, but one can obtain individual
coverings and block designs from the LJCR using the method
:meth:`designs.best_known_covering_design_from_LJCR
<sage.combinat.designs.covering_design.best_known_covering_design_www>`::

    sage: C = designs.best_known_covering_design_from_LJCR(7, 3, 2)   # optional - internet
    sage: C                            # optional - internet
    (7,3,2)-covering design of size 7
    Lower bound: 7
    Method: lex covering
    Submitted on: 1996-12-01 00:00:00
    sage: C.incidence_structure()      # optional - internet
    Incidence structure with 7 points and 7 blocks

Currently, this module gathers the following designs :

.. csv-table::
    :class: contentstable
    :widths: 100
    :delim: |

    :meth:`~sage.combinat.designs.block_design.ProjectiveGeometryDesign`
    :meth:`~sage.combinat.designs.block_design.ProjectivePlaneDesign`
    :meth:`~sage.combinat.designs.bibd.balanced_incomplete_block_design`
    :meth:`~sage.combinat.designs.block_design.AffineGeometryDesign`
    :meth:`~sage.combinat.designs.block_design.WittDesign`
    :meth:`~sage.combinat.designs.block_design.HadamardDesign`
    :meth:`~sage.combinat.designs.block_design.Hadamard3Design`
    :meth:`~sage.combinat.designs.latin_squares.mutually_orthogonal_latin_squares`
    :meth:`~sage.combinat.designs.orthogonal_arrays.transversal_design`
    :meth:`~sage.combinat.designs.orthogonal_arrays.orthogonal_array`
    :meth:`~sage.combinat.designs.bibd.steiner_triple_system`
    :meth:`~sage.combinat.designs.steiner_quadruple_systems.steiner_quadruple_system`

And the :meth:`designs.best_known_covering_design_from_LJCR
<sage.combinat.designs.covering_design.best_known_covering_design_www>` function
which queries the LJCR.

.. TODO::

    Implement DerivedDesign and ComplementaryDesign.

REFERENCES:

.. [1] La Jolla Covering Repository,
  http://www.ccrwest.org/cover.html
"""
from sage.combinat.designs.block_design import (ProjectiveGeometryDesign,
                                                DesarguesianProjectivePlaneDesign,
                                                projective_plane,
                                                AffineGeometryDesign,
                                                WittDesign,
                                                HadamardDesign,
                                                Hadamard3Design)

from sage.combinat.designs.steiner_quadruple_systems import steiner_quadruple_system

from sage.combinat.designs.covering_design import best_known_covering_design_www as best_known_covering_design_from_LJCR

from sage.combinat.designs.latin_squares import mutually_orthogonal_latin_squares

from sage.combinat.designs.orthogonal_arrays import transversal_design, orthogonal_array

from sage.combinat.designs.bibd import balanced_incomplete_block_design, steiner_triple_system

# deprecated in june 2014 (#16446)
def BalancedIncompleteBlockDesign(v,k,existence=False,use_LJCR=False):
    r"""
    This function is deprecated.

    TESTS::

        sage: bibd = designs.BalancedIncompleteBlockDesign(21,5)
        doctest:...: DeprecationWarning: designs.BalancedIncompleteBlockDesign is
        deprecated. Please use designs.balanced_incomplete_block_design instead.
        See http://trac.sagemath.org/16446 for details.
        sage: bibd
        Incidence structure with 21 points and 21 blocks
    """
    from sage.misc.superseded import deprecation
    deprecation(16446, "designs.BalancedIncompleteBlockDesign is deprecated. Please use designs.balanced_incomplete_block_design instead.")

    return balanced_incomplete_block_design(v, k, existence=existence, use_LJCR=use_LJCR)

