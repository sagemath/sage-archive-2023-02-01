r"""
Chain complexes

AUTHORS:

- John H. Palmieri (2009-04)

This module implements bounded chain complexes of free `R`-modules,
for any commutative ring `R` (although the interesting things, like
homology, only work if `R` is the integers or a field).

Fix a ring `R`.  A chain complex over `R` is a collection of
`R`-modules `\{C_n\}` indexed by the integers, with `R`-module maps
`d_n : C_n \rightarrow C_{n+1}` such that `d_{n+1} \circ d_n = 0` for
all `n`.  The maps `d_n` are called *differentials*.

One can vary this somewhat: the differentials may decrease degree by
one instead of increasing it: sometimes a chain complex is defined
with `d_n : C_n \rightarrow C_{n-1}` for each `n`. Indeed, the
differentials may change dimension by any fixed integer.

Also, the modules may be indexed over an abelian group other than the
integers, e.g., `\ZZ^{m}` for some integer `m \geq 1`, in which
case the differentials may change the grading by any element of that
grading group.

In this implementation, the ring `R` must be commutative and the
modules `C_n` must be free `R`-modules.  As noted above, homology
calculations will only work if the ring `R` is either `\ZZ` or a field.
The modules may be indexed by any free abelian group.  The
differentials may increase degree by 1 or decrease it, or indeed
change it by any fixed amount: this is controlled by the ``degree``
parameter used in defining the chain complex.
"""


########################################################################
#       Copyright (C) 2013 John H. Palmieri <palmieri@math.washington.edu>
#                          Volker Braun <vbraun.name@gmail.com>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#  as published by the Free Software Foundation; either version 2 of
#  the License, or (at your option) any later version.
#
#                  http://www.gnu.org/licenses/
########################################################################

from copy import copy

from sage.structure.parent import Parent
from sage.structure.element import ModuleElement
from sage.misc.cachefunc import cached_method

from sage.rings.integer_ring import ZZ
from sage.rings.rational_field import QQ
from sage.modules.free_module import FreeModule
from sage.modules.free_module_element import vector
from sage.matrix.matrix0 import Matrix
from sage.matrix.constructor import matrix, prepare_dict
from sage.misc.latex import latex
from sage.rings.all import GF, prime_range
from sage.misc.decorators import rename_keyword


def _latex_module(R, m):
    """
    LaTeX string representing a free module over ``R`` of rank ``m``.

    INPUT:

    - ``R`` -- a commutative ring
    - ``m`` -- non-negative integer

    This is used by the ``_latex_`` method for chain complexes.

    EXAMPLES::

        sage: from sage.homology.chain_complex import _latex_module
        sage: _latex_module(ZZ, 3)
        '\\Bold{Z}^{3}'
        sage: _latex_module(ZZ, 0)
        '0'
        sage: _latex_module(GF(3), 1)
        '\\Bold{F}_{3}^{1}'
    """
    if m == 0:
        return str(latex(0))
    else:
        return str(latex(FreeModule(R, m)))


@rename_keyword(deprecation=15151, check_products='check', check_diffs='check')
def ChainComplex(data=None, **kwds):
    r"""
    Define a chain complex.

    INPUT:

    -  ``data`` -- the data defining the chain complex; see below for
       more details.

    -  ``base_ring`` -- a commutative ring (optional), the ring over
       which the chain complex is defined. If this is not specified,
       it is determined by the data defining the chain complex.

    - ``grading_group`` -- a additive free abelian group (optional,
       default ``ZZ``), the group over which the chain complex is
       indexed.

    -  ``degree`` -- element of grading_group (optional, default 1),
       the degree of the differential.

    - ``check`` -- boolean (optional, default ``True``). If ``True``,
       check that each consecutive pair of differentials are
       composable and have composite equal to zero.

    OUTPUT: a chain complex

    .. WARNING::

       Right now, homology calculations will only work if the base
       ring is either `\ZZ` or a field, so please take this into account
       when defining a chain complex.

    Use data to define the chain complex.  This may be in any of the
    following forms.

    1. a dictionary with integers (or more generally, elements of
       grading_group) for keys, and with ``data[n]`` a matrix representing
       (via left multiplication) the differential coming from degree
       `n`.  (Note that the shape of the matrix then determines the
       rank of the free modules `C_n` and `C_{n+d}`.)

    2. a list/tuple/iterable of the form `[C_0, d_0, C_1, d_1, C_2,
       d_2, ...]`, where each `C_i` is a free module and each `d_i` is
       a matrix, as above.  This only makes sense if ``grading_group``
       is `\ZZ` and ``degree`` is 1.

    3. a list/tuple/iterable of the form `[r_0, d_0, r_1, d_1, r_2,
       d_2, ...]`, where `r_i` is the rank of the free module `C_i`
       and each `d_i` is a matrix, as above.  This only makes sense if
       ``grading_group`` is `\ZZ` and ``degree`` is 1.

    4. a list/tuple/iterable of the form `[d_0, d_1, d_2, ...]` where
       each `d_i` is a matrix, as above.  This only makes sense if
       ``grading_group`` is `\ZZ` and ``degree`` is 1.

    .. NOTE::

       In fact, the free modules `C_i` in case 2 and the ranks `r_i`
       in case 3 are ignored: only the matrices are kept, and from
       their shapes, the ranks of the modules are determined.
       (Indeed, if ``data`` is a list or tuple, then any element which
       is not a matrix is discarded; thus the list may have any number
       of different things in it, and all of the non-matrices will be
       ignored.)  No error checking is done to make sure, for
       instance, that the given modules have the appropriate ranks for
       the given matrices.  However, as long as ``check`` is True, the
       code checks to see if the matrices are composable and that each
       appropriate composite is zero.

    If the base ring is not specified, then the matrices are examined
    to determine a ring over which they are all naturally defined, and
    this becomes the base ring for the complex.  If no such ring can
    be found, an error is raised.  If the base ring is specified, then
    the matrices are converted automatically to this ring when
    defining the chain complex.  If some matrix cannot be converted,
    then an error is raised.

    EXAMPLES::

        sage: ChainComplex()
        Trivial chain complex over Integer Ring

        sage: C = ChainComplex({0: matrix(ZZ, 2, 3, [3, 0, 0, 0, 0, 0])})
        sage: C
        Chain complex with at most 2 nonzero terms over Integer Ring

        sage: m = matrix(ZZ, 2, 2, [0, 1, 0, 0])
        sage: D = ChainComplex([m, m], base_ring=GF(2)); D
        Chain complex with at most 3 nonzero terms over Finite Field of size 2
        sage: D == loads(dumps(D))
        True
        sage: D.differential(0)==m, m.is_immutable(), D.differential(0).is_immutable()
        (True, False, True)

    Note that when a chain complex is defined in Sage, new
    differentials may be created: every nonzero module in the chain
    complex must have a differential coming from it, even if that
    differential is zero::

        sage: IZ = ChainComplex({0: identity_matrix(ZZ, 1)})
        sage: IZ.differential()  # the differentials in the chain complex
        {0: [1], 1: [], -1: []}
        sage: IZ.differential(1).parent()
        Full MatrixSpace of 0 by 1 dense matrices over Integer Ring
        sage: mat = ChainComplex({0: matrix(ZZ, 3, 4)}).differential(1)
        sage: mat.nrows(), mat.ncols()
        (0, 3)

    Defining the base ring implicitly::

        sage: ChainComplex([matrix(QQ, 3, 1), matrix(ZZ, 4, 3)])
        Chain complex with at most 3 nonzero terms over Rational Field
        sage: ChainComplex([matrix(GF(125, 'a'), 3, 1), matrix(ZZ, 4, 3)])
        Chain complex with at most 3 nonzero terms over Finite Field in a of size 5^3

    If the matrices are defined over incompatible rings, an error results::

        sage: ChainComplex([matrix(GF(125, 'a'), 3, 1), matrix(QQ, 4, 3)])
        Traceback (most recent call last):
        ...
        TypeError: unable to find a common ring for all elements

    If the base ring is given explicitly but is not compatible with
    the matrices, an error results::

        sage: ChainComplex([matrix(GF(125, 'a'), 3, 1)], base_ring=QQ)
        Traceback (most recent call last):
        ...
        TypeError: Unable to coerce 0 (<type 
        'sage.rings.finite_rings.element_givaro.FiniteField_givaroElement'>) to Rational
    """
    check = kwds.get('check', True)
    base_ring = kwds.get('base_ring', None)
    grading_group = kwds.get('grading_group', ZZ)
    degree = kwds.get('degree', 1)
    try:
        degree = grading_group(degree)
    except StandardError:
        raise ValueError('degree is not an element of the grading group')

    if data is None or (isinstance(data, (list, tuple)) and len(data) == 0):
        # the zero chain complex
        try:
            zero = grading_group.identity()
        except AttributeError:
            zero = grading_group.zero_element()
        if base_ring is None:
            base_ring = ZZ
        data_dict = {zero: matrix(base_ring, 0, 0, [])}
    elif isinstance(data, dict):  # data is dictionary
        data_dict = data
    else: # data is list/tuple/iterable
        data_matrices = filter(lambda x: isinstance(x, Matrix), data)
        if degree != 1:
            raise ValueError('degree must be +1 if the data argument is a list or tuple')
        if grading_group != ZZ:
            raise ValueError('grading_group must be ZZ if the data argument is a list or tuple')
        data_dict = dict((grading_group(i), m) for i,m in enumerate(data_matrices))

    if base_ring is None:
        _, base_ring = prepare_dict(dict([n, data_dict[n].base_ring()(0)] for n in data_dict))

    # make sure values in data_dict are appropriate matrices
    for n in data_dict.keys():
        if not n in grading_group:
            raise ValueError('one of the dictionary keys is not an element of the grading group')
        mat = data_dict[n]
        if not isinstance(mat, Matrix):
            raise TypeError('One of the differentials in the data is not a matrix')
        if mat.base_ring() is base_ring: 
            if not mat.is_immutable():
                mat = copy(mat)  # do not make any arguments passed immutable
                mat.set_immutable()
        else:
            mat = mat.change_ring(base_ring)
            mat.set_immutable()
        data_dict[n] = mat

    # include any "obvious" zero matrices that are not 0x0 
    for n in data_dict.keys():  # note: data_dict will be mutated in this loop
        mat1 = data_dict[n]
        if (n+degree not in data_dict) and (mat1.nrows() != 0):
            if n+2*degree in data_dict:
                mat2 = matrix(base_ring, data_dict[n+2*degree].ncols(), mat1.nrows())
            else:
                mat2 = matrix(base_ring, 0, mat1.nrows())
            mat2.set_immutable()
            data_dict[n+degree] = mat2
        if (n-degree not in data_dict) and (mat1.ncols() != 0):
            if n-2*degree in data_dict:
                mat0 = matrix(base_ring, mat1.ncols(), data_dict[n-2*degree].nrows())
            else:
                mat0 = matrix(base_ring, mat1.ncols(), 0)
            mat0.set_immutable()
            data_dict[n-degree] = mat0
        
    # check that this is a complex: going twice is zero
    if check:
        for n in data_dict.keys():
            mat0 = data_dict[n]
            try:
                mat1 = data_dict[n+degree]
            except KeyError:
                continue
            try:
                prod = mat1 * mat0
            except TypeError:
                raise TypeError('the differentials d_{%s} and d_{%s} are not compatible: '
                                'their product is not defined' % (n, n+degree))
            if not prod.is_zero():
                raise ValueError('the differentials d_{%s} and d_{%s} are not compatible: '
                                 'their composition is not zero.' % (n, n+degree))

    return ChainComplex_class(grading_group, degree, base_ring, data_dict)
    

class Chain_class(ModuleElement):
    
    def __init__(self, parent, *args, **kwds):
        """
        Parent for all chain complexes over a given ``base_ring``

        EXAMPLES::

            sage: C = ChainComplex({0: matrix(ZZ, 2, 3, [3, 0, 0, 0, 0, 0])}, base_ring=GF(7))
            sage: C.category()
            Category of chain complexes over Finite Field of size 7
        """
        super(Chain_class, self).__init__(parent)

    def is_cycle(self):
        pass  # TODO
    
    def is_boundary(self):
        pass  # TODO



class ChainComplex_class(Parent):

    def __init__(self, grading_group, degree_of_differential, base_ring, differentials):
        """
        See :func:`ChainComplex` for full documentation.

        EXAMPLES::

            sage: C = ChainComplex(); C
            Trivial chain complex over Integer Ring

            sage: D = ChainComplex({0: matrix(ZZ, 2, 3, [3, 0, 0, 0, 0, 0])})
            sage: D
            Chain complex with at most 2 nonzero terms over Integer Ring

            sage: ChainComplex().base_ring()
            Integer Ring
        """
        assert all(d.base_ring() == base_ring and d.is_immutable()
                   for d in differentials.values())
        assert degree_of_differential.parent() is grading_group
        assert grading_group is ZZ or not grading_group.is_multiplicative()
        # all differentials that are not 0x0 must be specified to the constructor
        assert all(dim+degree_of_differential in differentials or d.nrows() == 0
                   for dim, d in differentials.iteritems())
        assert all(dim-degree_of_differential in differentials or d.ncols() == 0
                   for dim, d in differentials.iteritems())
        self._grading_group = grading_group
        self._degree_of_differential = degree_of_differential
        self._diff = differentials

        from sage.categories.all import ChainComplexes
        category = ChainComplexes(base_ring)
        super(ChainComplex_class, self).__init__(base=base_ring, category=category)

    Element = Chain_class

    def _element_constructor_(self, *args):
        """
        The element constructor.

        This is part of the Parent/Element framework. Calling the
        parent uses this method to construct elements.

        EXAMPLES::

            sage: D = ChainComplex({0: matrix(ZZ, 2, 2, [1,0,0,2])})
            sage: D._element_constructor_(0)
        """
        return self.element_class(self, *args)

    @cached_method
    def rank(self, degree, ring=None):
        r"""
        Return the rank of a differential

        INPUT:
        
        - ``degree`` -- an element `\delta` of the grading
          group. Which differential `d_{\delta}` we want to know the
          rank of.
        
        - ``ring`` -- a commutative ring `S` or ``None`` (default). If
          specified, the rank is computed after changing to this ring.
    
        OUTPUT:

        The rank of the differential ``d_{\delta} \otimes_R S`, where
        `R` is the base ring of the chain complex.

        EXAMPLES::

            sage: C = ChainComplex({0:matrix(ZZ, [[2]])})
            sage: C.differential(0)
            [2]
            sage: C.rank(0)
            1
            sage: C.rank(0, ring=GF(2))
            0
        """
        degree = self.grading_group()(degree)
        try:
            d = self._diff[degree]
        except IndexError:
            return ZZ.zero()
        if d.nrows() == 0 or d.ncols() == 0:
            return ZZ.zero()
        if ring is None:
            return d.rank()
        else:
            return d.change_ring(ring).rank()

    def grading_group(self):
        r"""
        Return the grading group

        OUTPUT:

        The discrete abelian group that indexes the individual modules
        of the complex. Usually `\ZZ`.

        EXAMPLES::

            sage: G = AdditiveAbelianGroup([0, 3])
            sage: C = ChainComplex(grading_group=G, degree=G([1,2]))
            sage: C.grading_group()
            Additive abelian group isomorphic to Z/3 + Z
            sage: C.degree_of_differential()
            (2, 1)
        """
        return self._grading_group

    def degree_of_differential(self):
        """
        Return the degree of the differentials of the complex

        OUTPUT:

        An element of the grading group.

        EXAMPLES::

            sage: D = ChainComplex({0: matrix(ZZ, 2, 2, [1,0,0,2])})
            sage: D.degree_of_differential()
            1
        """
        return self._degree_of_differential

    def differential(self, dim=None):
        """
        The differentials which make up the chain complex.

        INPUT:

        - ``dim`` -- element of the grading group (optional, default
           ``None``).  If this is ``None``, return a dictionary of all
           of the differentials.  If this is a single element, return
           the differential starting in that dimension.

        OUTPUT: either a dictionary of all of the differentials or a single
        differential (i.e., a matrix)

        EXAMPLES::

            sage: D = ChainComplex({0: matrix(ZZ, 2, 2, [1,0,0,2])})
            sage: D.differential()
            {0: [1 0]
            [0 2], 1: [], -1: []}
            sage: D.differential(0)
            [1 0]
            [0 2]
            sage: C = ChainComplex({0: identity_matrix(ZZ, 40)})
            sage: C.differential()
            {0: 40 x 40 dense matrix over Integer Ring, 1: [], -1: 40 x 0 dense matrix over Integer Ring}
        """
        if dim is None:
            return copy(self._diff)
        dim = self.grading_group()(dim)
        try:
            return self._diff[dim]
        except KeyError:
            pass
        # all differentials that are not 0x0 are in self._diff
        return matrix(self.base_ring(), 0, 0)

    def dual(self):
        """
        The dual chain complex to ``self``.

        Since all modules in ``self`` are free of finite rank, the
        dual in dimension `n` is isomorphic to the original chain
        complex in dimension `n`, and the corresponding boundary
        matrix is the transpose of the matrix in the original complex.
        This converts a chain complex to a cochain complex and vice
        versa.

        EXAMPLES::

            sage: C = ChainComplex({2: matrix(ZZ, 2, 3, [3, 0, 0, 0, 0, 0])})
            sage: C.degree_of_differential()
            1
            sage: C.differential(2)
            [3 0 0]
            [0 0 0]
            sage: C.dual().degree_of_differential()
            -1
            sage: C.dual().differential(3)
            [3 0]
            [0 0]
            [0 0]
        """
        data = {}
        deg = self.degree_of_differential()
        for d in self.differential():
            data[(d+deg)] = self.differential()[d].transpose()
        return ChainComplex(data, degree=-deg)

    def free_module(self):
        """
        The free module underlying this chain complex.

        EXAMPLES::

            sage: C = ChainComplex({0: matrix(ZZ, 2, 3, [3, 0, 0, 0, 0, 0]), 1: matrix(ZZ, 0, 2)})
            sage: C.free_module()
            Ambient free module of rank 5 over the principal ideal domain Integer Ring

        This defines the forgetful functor from the category of chain
        complexes to the category of free modules::

            sage: FreeModules(ZZ)(C)
            Ambient free module of rank 5 over the principal ideal domain Integer Ring
        """
        rank = sum([mat.ncols() for mat in self.differential().values()])
        return FreeModule(self.base_ring(), rank)

    def __cmp__(self, other):
        """
        Return ``True`` iff this chain complex is the same as other: that
        is, if the base rings and the matrices of the two are the
        same.

        EXAMPLES::

            sage: C = ChainComplex({0: matrix(ZZ, 2, 3, [3, 0, 0, 0, 0, 0])}, base_ring=GF(2))
            sage: D = ChainComplex({0: matrix(GF(2), 2, 3, [1, 0, 0, 0, 0, 0]), 1: matrix(ZZ, 0, 2), 3: matrix(ZZ, 0, 0)})  # base_ring determined from the matrices
            sage: C == D
            True
        """
        if not isinstance(other, ChainComplex_class):
            return cmp(type(other), ChainComplex_class)
        c = cmp(self.base_ring(), other.base_ring())
        if c != 0:
            return c
        R = self.base_ring()
        equal = True
        for d,mat in self.differential().iteritems():
            if d not in other.differential():
                equal = equal and mat.ncols() == 0 and mat.nrows() == 0
            else:
                equal = (equal and
                         other.differential()[d].change_ring(R) == mat.change_ring(R))
        for d,mat in other.differential().iteritems():
            if d not in self.differential():
                equal = equal and mat.ncols() == 0 and mat.nrows() == 0
        if equal:
            return 0
        return -1

    def homology(self, dim=None, **kwds):
        r"""
        The homology of the chain complex in the given dimension.

        INPUT:

        -  ``dim`` -- an element of the grading group for the chain
           complex (optional, default ``None``): the degree in which to
           compute homology. If this is ``None``, return the homology in
           every dimension in which the chain complex is possibly
           nonzero.

        -  ``base_ring`` -- a commutative ring (optional, default is the
           base ring for the chain complex).  Must be either the
           integers `\ZZ` or a field.

        -  ``generators`` -- boolean (optional, default ``False``).  If
           ``True``, return generators for the homology groups along with
           the groups. See :trac:`6100`.

        -  ``verbose`` - boolean (optional, default ``False``).  If
           ``True``, print some messages as the homology is computed.

        -  ``algorithm`` - string (optional, default ``'auto'``).  The
           options are ``'auto'``, ``'dhsw'``, ``'pari'`` or ``'no_chomp'``.
           See below for descriptions.

        OUTPUT:

        If dim is specified, the homology in dimension ``dim``.
        Otherwise, the homology in every dimension as a dictionary
        indexed by dimension.

        ALGORITHM:

        If ``algorithm`` is set to ``'auto'`` (the default), then use
        CHomP if available.  (CHomP is available at the web page
        http://chomp.rutgers.edu/.  It is also an experimental package
        for Sage.)

        CHomP computes homology, not cohomology, and only works over
        the integers or finite prime fields.  Therefore if any of
        these conditions fails, or if CHomP is not present, or if
        ``algorithm`` is set to 'no_chomp', go to plan B: if ``self``
        has a ``_homology`` method -- each simplicial complex has
        this, for example -- then call that.  Such a method implements
        specialized algorithms for the particular type of cell
        complex.

        Otherwise, move on to plan C: compute the chain complex of
        ``self`` and compute its homology groups.  To do this: over a
        field, just compute ranks and nullities, thus obtaining
        dimensions of the homology groups as vector spaces.  Over the
        integers, compute Smith normal form of the boundary matrices
        defining the chain complex according to the value of
        ``algorithm``.  If ``algorithm`` is ``'auto'`` or ``'no_chomp'``,
        then for each relatively small matrix, use the standard Sage
        method, which calls the Pari package.  For any large matrix,
        reduce it using the Dumas, Heckenbach, Saunders, and Welker
        elimination algorithm [DHSW]_: see
        :func:`~sage.homology.matrix_utils.dhsw_snf` for details.

        Finally, ``algorithm`` may also be ``'pari'`` or ``'dhsw'``, which
        forces the named algorithm to be used regardless of the size
        of the matrices and regardless of whether CHomP is available.

        As of this writing, CHomP is by far the fastest option,
        followed by the ``'auto'`` or ``'no_chomp'`` setting of using the
        Dumas, Heckenbach, Saunders, and Welker elimination algorithm
        [DHSW]_ for large matrices and Pari for small ones.

        .. WARNING::

           This only works if the base ring is the integers or a
           field.  Other values will return an error.

        EXAMPLES::

            sage: C = ChainComplex({0: matrix(ZZ, 2, 3, [3, 0, 0, 0, 0, 0])})
            sage: C.homology()
            {0: Z x Z, 1: Z x C3}
            sage: C.homology(dim=1, base_ring = GF(3))
            Vector space of dimension 2 over Finite Field of size 3
            sage: D = ChainComplex({0: identity_matrix(ZZ, 4), 4: identity_matrix(ZZ, 30)})
            sage: D.homology()
            {0: 0, 1: 0, 4: 0, 5: 0}

        Generators: generators are given as
        a list of cycles, each of which is an element in the
        appropriate free module, and hence is represented as a vector::

            sage: C.homology(1, generators=True)  # optional - CHomP
            (Z x C3, [(0, 1), (1, 0)])

        Tests for :trac:`6100`, the Klein bottle with generators::

            sage: d0 = matrix(ZZ, 0,1)
            sage: d1 = matrix(ZZ, 1,3, [[0,0,0]])
            sage: d2 = matrix(ZZ, 3,2, [[1,1], [1,-1], [-1,1]])
            sage: C_k = ChainComplex({0:d0, 1:d1, 2:d2}, degree=-1)
            sage: C_k.homology(generators=true)   # optional - CHomP
            {0: (Z, [(1)]), 1: (Z x C2, [(0, 0, 1), (0, 1, -1)])}

        From a torus using a field::

            sage: T = simplicial_complexes.Torus()
            sage: C_t = T.chain_complex()
            sage: C_t.homology(base_ring=QQ, generators=True)
            {0: [(Vector space of dimension 1 over Rational Field, (0, 0, 0, 0, 0, 0, 1))],
             1: [(Vector space of dimension 1 over Rational Field,
               (0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -1, -1, 0, 0, 0, -1, 0, 1, 0)),
              (Vector space of dimension 1 over Rational Field,
               (0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 1, 0, -1, -1))],
             2: [(Vector space of dimension 1 over Rational Field,
               (1, -1, -1, -1, 1, -1, -1, 1, 1, 1, 1, 1, -1, -1))]}
        """
        from sage.homology.homology_group import HomologyGroup
        from sage.interfaces.chomp import have_chomp, homchain

        if dim is not None and dim not in self.grading_group():
            raise ValueError('dimension is not an element of the grading group')

        algorithm = kwds.get('algorithm', 'auto')
        if algorithm not in ['dhsw', 'pari', 'auto', 'no_chomp']:
            raise NotImplementedError('algorithm not recognized')

        verbose = kwds.get('verbose', False)
        base_ring = kwds.get('base_ring', None)
        generators = kwds.get('generators', False)
        if base_ring is None or base_ring == self.base_ring():
            change_ring = False
            base_ring = self.base_ring()
        else:
            change_ring = True
        if not (base_ring.is_field() or base_ring is ZZ):
            raise NotImplementedError('can only compute homology if the base ring is the integers or a field')

        # try to use CHomP if working over Z or F_p, p a prime.
        if (algorithm == 'auto' and (base_ring == ZZ or
            (base_ring.is_prime_field() and base_ring != QQ))):
            # compute all of homology, then pick off requested dimensions
            H = None
            if have_chomp('homchain'):
                H = homchain(self, **kwds)
                # now pick off the requested dimensions
                if H:
                    if dim is not None:
                        answer = {}
                        if isinstance(dim, (list, tuple)):
                            for d in dim:
                                if d in H:
                                    answer[d] = H[d]
                                else:
                                    answer[d] = HomologyGroup(0, base_ring)
                        else:
                            if dim in H:
                                answer = H[dim]
                            else:
                                answer = HomologyGroup(0, base_ring)
                    else:
                        answer = H
                    return answer
                else:
                    if verbose:
                        print "ran CHomP, but no output."

        # if dim is None, return all of the homology groups
        degree = self.degree_of_differential()
        if dim is None:
            answer = {}
            for n in self._diff.keys():
                if n-degree not in self._diff:
                    continue
                if verbose:
                    print "Computing homology of the chain complex in dimension %s..." % n
                if base_ring == self.base_ring():
                    answer[n] = self.homology(n, verbose=verbose,
                                              generators=generators,
                                              algorithm=algorithm)
                else:
                    answer[n] = self.homology(n, base_ring=base_ring,
                                              verbose=verbose,
                                              generators=generators,
                                              algorithm=algorithm)
            return answer

        # now compute the homology in the given dimension
        if dim in self._diff:
            # d_out is the differential going out of degree dim,
            # d_in is the differential entering degree dim
            d_out_cols = self._diff[dim].ncols()
            d_out_rows = self._diff[dim].nrows()
            if base_ring == ZZ:
                temp_ring = QQ
            else:
                temp_ring = base_ring
            d_out_rank = self.rank(dim, ring=temp_ring)

            if dim - degree in self._diff:
                if change_ring:
                    d_in = self._diff[dim-degree].change_ring(base_ring)
                else:
                    d_in = self._diff[dim-degree]

                if generators:
                    # Find the kernel of the out-going differential.
                    K = self._diff[dim].right_kernel().matrix().transpose().change_ring(base_ring)

                    # Compute the induced map to the kernel
                    S = K.augment(d_in).hermite_form()
                    d_in_induced = S.submatrix(row=0, nrows=d_in.nrows()-d_out_rank,
                                               col=d_in.nrows()-d_out_rank, ncols=d_in.ncols())

                    # Find the SNF of the induced matrix and appropriate generators
                    (N, P, Q) = d_in_induced.smith_form()
                    all_divs = [0]*N.nrows()
                    non_triv = 0
                    for i in range(0, N.nrows()):
                        if i >= N.ncols():
                            break
                        all_divs[i] = N[i][i]
                        if N[i][i] == 1:
                            non_triv = non_triv + 1
                    divisors = filter(lambda x: x != 1, all_divs)
                    gens = (K * P.inverse().submatrix(col=non_triv)).transpose()
                    answer = [(HomologyGroup(1, base_ring, [divisors[i]]), gens[i])
                              for i in range(len(divisors))]
                else:
                    if base_ring.is_field():
                        null = d_out_cols - d_out_rank
                        rk = self.rank(dim-degree, ring=temp_ring)
                        answer = HomologyGroup(null - rk, base_ring)
                    elif base_ring == ZZ:
                        nullity = d_out_cols - d_out_rank
                        if d_in.ncols() == 0:
                            all_divs = [0] * nullity
                        else:
                            if algorithm == 'auto':
                                if ((d_in.ncols() > 300 and d_in.nrows() > 300)
                                    or (min(d_in.ncols(), d_in.nrows()) > 100 and
                                        d_in.ncols() + d_in.nrows() > 600)):
                                    algorithm = 'dhsw'
                                else:
                                    algorithm = 'pari'
                            if algorithm == 'dhsw':
                                from sage.homology.matrix_utils import dhsw_snf
                                all_divs = dhsw_snf(d_in, verbose=verbose)
                            else:
                                algorithm = 'pari'
                                if d_in.is_sparse():
                                    all_divs = d_in.dense_matrix().elementary_divisors(algorithm)
                                else:
                                    all_divs = d_in.elementary_divisors(algorithm)
                        all_divs = all_divs[:nullity]
                        # divisors equal to 1 produce trivial
                        # summands, so filter them out
                        divisors = filter(lambda x: x != 1, all_divs)
                        answer = HomologyGroup(len(divisors), base_ring, divisors)
            else: # no incoming differential: it's zero
                answer = HomologyGroup(d_out_cols - d_out_rank, base_ring)
                if generators: #Include the generators of the nullspace
                    # Find the kernel of the out-going differential.
                    K = self._diff[dim].right_kernel().matrix().change_ring(base_ring)
                    answer = [( answer, vector(base_ring, K.list()) )]
        else:  # chain complex is zero here, so return the zero module
            answer = HomologyGroup(0, base_ring)
            if generators:
                answer = [(answer, vector(base_ring, []))]
        if verbose:
            print "  Homology is %s" % answer
        return answer

    def betti(self, dim=None, **kwds):
        """
        The Betti number of the homology of the chain complex in this
        dimension.

        That is, write the homology in this dimension as a direct sum
        of a free module and a torsion module; the Betti number is the
        rank of the free summand.

        INPUT:

        -  ``dim`` -- an element of the grading group for the chain
           complex or None (optional, default ``None``).  If ``None``,
           then return every Betti number, as a dictionary indexed by
           degree.  If an element of the grading group, then return
           the Betti number in that dimension.

        -  ``base_ring`` -- a commutative ring (optional, default is the
           base ring for the chain complex).  Compute homology with
           these coefficients.  Must be either the integers or a
           field.

        OUTPUT: the Betti number in dimension ``dim`` - the rank of
        the free part of the homology module in this dimension.

        EXAMPLES::

            sage: C = ChainComplex({0: matrix(ZZ, 2, 3, [3, 0, 0, 0, 0, 0])})
            sage: C.betti(0)
            2
            sage: [C.betti(n) for n in range(5)]
            [2, 1, 0, 0, 0]
            sage: C.betti()
            {0: 2, 1: 1}
        """
        base_ring = kwds.get('base_ring', None)

        if base_ring is None:
            base_ring = self.base_ring()
        if base_ring == ZZ:
            base_ring = QQ
        if base_ring.is_field():
            kwds['base_ring'] = base_ring
            H = self.homology(dim=dim, **kwds)
            if isinstance(H, dict):
                return dict([i, H[i].dimension()] for i in H)
            else:
                return H.dimension()
        else:
            raise NotImplementedError, "Not implemented: unable to compute Betti numbers if the base ring is not ZZ or a field."

    def torsion_list(self, max_prime, min_prime=2):
        r"""
        Look for torsion in this chain complex by computing its mod `p`
        homology for a range of primes `p`.

        INPUT:

        -  ``max_prime`` -- prime number: search for torsion mod `p` for
           all `p` strictly less than this number.

        -  ``min_prime`` -- prime (optional, default 2): search for
           torsion mod `p` for primes at least as big as this.

        Return a list of pairs (`p`, ``dims``) where `p` is a prime at
        which there is torsion and ``dims`` is a list of dimensions in
        which this torsion occurs.

        The base ring for the chain complex must be the integers; if
        not, an error is raised.

        Algorithm: let `C` denote the chain complex.  Let `P` equal
        ``max_prime``.  Compute the mod `P` homology of `C`, and use
        this as the base-line computation: the assumption is that this
        is isomorphic to the integral homology tensored with
        `\GF{P}`.  Then compute the mod `p` homology for a range of
        primes `p`, and record whenever the answer differs from the
        base-line answer.

        EXAMPLES::

            sage: C = ChainComplex({0: matrix(ZZ, 2, 3, [3, 0, 0, 0, 0, 0])})
            sage: C.homology()
            {0: Z x Z, 1: Z x C3}
            sage: C.torsion_list(11)
            [(3, [1])]
            sage: C = ChainComplex([matrix(ZZ, 1, 1, [2]), matrix(ZZ, 1, 1), matrix(1, 1, [3])])
            sage: C.homology(1)
            C2
            sage: C.homology(3)
            C3
            sage: C.torsion_list(5)
            [(2, [1]), (3, [3])]
        """
        if self.base_ring() != ZZ:
            raise NotImplementedError('only implemented for base ring the integers')
        answer = []
        torsion_free = self.betti(base_ring=GF(max_prime))
        for p in prime_range(min_prime, max_prime):
            mod_p_betti = self.betti(base_ring=GF(p))
            if mod_p_betti != torsion_free:
                diff_dict = {}
                temp_diff = {}
                D = self.degree_of_differential()
                for i in torsion_free:
                    temp_diff[i] = mod_p_betti.get(i, 0) - torsion_free[i]
                for i in temp_diff:
                    if temp_diff[i] > 0:
                        if i+D in diff_dict:
                            lower = diff_dict[i+D]
                        else:
                            lower = 0
                        current = temp_diff[i]
                        if current > lower:
                            diff_dict[i] = current - lower
                            if i-D in diff_dict:
                                diff_dict[i-D] -= current - lower
                differences = []
                for i in diff_dict:
                    if diff_dict[i] != 0:
                        differences.append(i)
                answer.append((p,differences))
        return answer

    def _Hom_(self, other, category=None):
        """
        Return the set of chain maps between chain complexes ``self``
        and ``other``.

        EXAMPLES::

            sage: S = simplicial_complexes.Sphere(2)
            sage: T = simplicial_complexes.Torus()
            sage: C = S.chain_complex(augmented=True,cochain=True)
            sage: D = T.chain_complex(augmented=True,cochain=True)
            sage: Hom(C,D)  # indirect doctest
            Set of Morphisms from Chain complex with at most 4 nonzero terms over Integer Ring to Chain complex with at most 4 nonzero terms over Integer Ring in Category of chain complexes over Integer Ring
        """
        from sage.homology.chain_complex_homspace import ChainComplexHomspace
        return ChainComplexHomspace(self, other)

    def _flip_(self):
        """
        Flip chain complex upside down (degree `n` gets changed to
        degree `-n`), thus turning a chain complex into a cochain complex
        without changing the homology (except for flipping it, too).

        EXAMPLES::

            sage: C = ChainComplex({2: matrix(ZZ, 2, 3, [3, 0, 0, 0, 0, 0])})
            sage: C.degree_of_differential()
            1
            sage: C.differential(2)
            [3 0 0]
            [0 0 0]
            sage: C._flip_().degree_of_differential()
            -1
            sage: C._flip_().differential(-2)
            [3 0 0]
            [0 0 0]
        """
        data = {}
        deg = self.degree_of_differential()
        for d in self.differential():
            data[-d] = self.differential()[d]
        return ChainComplex(data, degree=-deg)

    def _chomp_repr_(self):
        r"""
        String representation of ``self`` suitable for use by the CHomP
        program.

        Since CHomP can only handle chain complexes, not cochain
        complexes, and since it likes its complexes to start in degree
        0, flip the complex over if necessary, and shift it to start
        in degree 0.  Note also that CHomP only works over the
        integers or a finite prime field.

        EXAMPLES::

            sage: C = ChainComplex({-2: matrix(ZZ, 1, 3, [3, 0, 0])}, degree=-1)
            sage: C._chomp_repr_()
            'chain complex\n\nmax dimension = 1\n\ndimension 0\n   boundary a1 = 0\n\ndimension 1\n   boundary a1 = + 3 * a1 \n   boundary a2 = 0\n   boundary a3 = 0\n\n'
            sage: C = ChainComplex({-2: matrix(ZZ, 1, 3, [3, 0, 0])}, degree=1)
            sage: C._chomp_repr_()
            'chain complex\n\nmax dimension = 1\n\ndimension 0\n   boundary a1 = 0\n\ndimension 1\n   boundary a1 = + 3 * a1 \n   boundary a2 = 0\n   boundary a3 = 0\n\n'
        """
        deg = self.degree_of_differential()
        if (self.grading_group() != ZZ or
            (deg != 1 and deg != -1)):
            raise ValueError('CHomP only works on Z-graded chain complexes with '
                             'differential of degree 1 or -1')
        base_ring = self.base_ring()
        if (base_ring == QQ) or (base_ring != ZZ and not (base_ring.is_prime_field())):
            raise ValueError('CHomP doesn\'t compute over the rationals, only over Z or F_p')
        if deg == -1:
            diffs = self.differential()
        else:
            diffs = self._flip_().differential()

        maxdim = max(diffs)
        mindim = min(diffs)
        # will shift chain complex by subtracting mindim from
        # dimensions, so its bottom dimension is zero.
        s = "chain complex\n\nmax dimension = %s\n\n" % (maxdim - mindim - 1,)

        for i in range(0, maxdim - mindim):
            s += "dimension %s\n" % i
            mat = diffs.get(i + mindim, matrix(base_ring, 0, 0))
            for idx in range(mat.ncols()):
                s += "   boundary a%s = " % (idx + 1)
                # construct list of bdries
                col = mat.column(idx)
                if col.nonzero_positions():
                    for j in col.nonzero_positions():
                        entry = col[j]
                        if entry > 0:
                            sgn = "+"
                        else:
                            sgn = "-"
                            entry = -entry
                        s += "%s %s * a%s " % (sgn, entry, j+1)
                else:
                    s += "0"
                s += "\n"
            s += "\n"
        return s

    def _repr_(self):
        """
        Print representation.

        EXAMPLES::

            sage: C = ChainComplex({0: matrix(ZZ, 2, 3, [3, 0, 0, 0, 0, 0])})
            sage: C._repr_()
            'Chain complex with at most 2 nonzero terms over Integer Ring'
        """
        diffs = filter(lambda mat: mat.nrows() + mat.ncols() > 0,
                       self._diff.values())
        if len(diffs) == 0:
            s = 'Trivial chain complex'
        else:
            s = 'Chain complex with at most {0} nonzero terms'.format(len(diffs)-1)
        s += ' over {0}'.format(self.base_ring())
        return s

    def _ascii_art_(self):
        """
        EXAMPLES::

            sage: C = ChainComplex({0: matrix(ZZ, 2, 3, [3, 0, 0, 0, 0, 0])})
        """
        if self.grading_group() is not ZZ:
            return super(ChainComplex_class, self)._ascii_art_()
        from sage.misc.ascii_art import AsciiArt
        dmin = min(self._diff.keys())
        dmax = max(self._diff.keys())

    def _latex_(self):
        """
        LaTeX print representation.

        EXAMPLES::

            sage: C = ChainComplex({0: matrix(ZZ, 2, 3, [3, 0, 0, 0, 0, 0])})
            sage: C._latex_()
            '\\Bold{Z}^{3} \\xrightarrow{d_{0}} \\Bold{Z}^{2}'
        """
#         Warning: this is likely to screw up if, for example, the
#         degree of the differential is 2 and there are nonzero terms
#         in consecutive dimensions (e.g., in dimensions 0 and 1).  In
#         such cases, the representation might show a differential
#         connecting these terms, although the differential goes from
#         dimension 0 to dimension 2, and from dimension 1 to
#         dimension 3, etc.  I don't know how much effort should be
#         put into trying to fix this.
        string = ""
        dict = self._diff
        deg = self.degree_of_differential()
        ring = self.base_ring()
        if self.grading_group() != ZZ:
            guess = dict.keys()[0]
            if guess-deg in dict:
                string += "\\dots \\xrightarrow{d_{%s}} " % latex(guess-deg)
            string += _latex_module(ring, mat.ncols())
            string += " \\xrightarrow{d_{%s}} \\dots" % latex(guess)
        else:
            backwards = (deg < 0)
            sorted_list = sorted(dict.keys(), reverse=backwards)
            if len(dict) <= 6:
                for n in sorted_list[1:-1]:
                    mat = dict[n]
                    string += _latex_module(ring, mat.ncols())
                    string += " \\xrightarrow{d_{%s}} " % latex(n)
                mat = dict[sorted_list[-1]]
                string += _latex_module(ring, mat.ncols())
            else:
                for n in sorted_list[:2]:
                    mat = dict[n]
                    string += _latex_module(ring, mat.ncols())
                    string += " \\xrightarrow{d_{%s}} " % latex(n)
                string += "\\dots "
                n = sorted_list[-2]
                string += "\\xrightarrow{d_{%s}} " % latex(n)
                mat = dict[sorted_list[-1]]
                string += _latex_module(ring, mat.ncols())
        return string


from sage.structure.sage_object import register_unpickle_override
register_unpickle_override('sage.homology.chain_complex', 'ChainComplex', ChainComplex_class)

