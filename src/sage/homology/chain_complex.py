r"""
Chain complexes

AUTHORS:

- John H. Palmieri (2009-04)

This module implements chain complexes of free `R`-modules, for any
commutative ring `R` (although the interesting things, like homology,
only work if `R` is the integers or a field).

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

from sage.structure.sage_object import SageObject
from sage.rings.integer_ring import ZZ
from sage.rings.rational_field import QQ
from sage.modules.free_module import FreeModule, VectorSpace
from sage.matrix.matrix0 import Matrix
from sage.matrix.constructor import matrix, prepare_dict
from sage.misc.latex import latex
from sage.groups.additive_abelian.additive_abelian_group import AdditiveAbelianGroup_fixed_gens
from sage.modules.free_module import FreeModule
#from sage.groups.abelian_gps.abelian_group import AbelianGroup_class
from sage.rings.all import GF, prime_range

def dhsw_snf(mat, verbose=False):
    """
    Preprocess a matrix using the "Elimination algorithm" described by
    Dumas et al., and then call ``elementary_divisors`` on the
    resulting (smaller) matrix.

    'dhsw' stands for 'Dumas, Heckenbach, Saunders, Welker,' and 'snf'
    stands for 'Smith Normal Form.'

    INPUT:

    -  ``mat`` - an integer matrix, either sparse or dense.

    (They use the transpose of the matrix considered here, so they use
    rows instead of columns.)

    Algorithm: go through ``mat`` one column at a time.  For each
    column, add multiples of previous columns to it until either

      - it's zero, in which case it should be deleted.

      - its first nonzero entry is 1 or -1, in which case it should be kept.

      - its first nonzero entry is something else, in which case it is
        deferred until the second pass.

    Then do a second pass on the deferred columns.

    At this point, the columns with 1 or -1 in the first entry
    contribute to the rank of the matrix, and these can be counted and
    then deleted (after using the 1 or -1 entry to clear out its row).
    Suppose that there were `N` of these.

    The resulting matrix should be much smaller; we then feed it
    to Sage's ``elementary_divisors`` function, and prepend `N` 1's to
    account for the rows deleted in the previous step.

    EXAMPLES::

        sage: from sage.homology.chain_complex import dhsw_snf
        sage: mat = matrix(ZZ, 3, 4, range(12))
        sage: dhsw_snf(mat)
        [1, 4, 0]
        sage: mat = random_matrix(ZZ, 20, 20, x=-1, y=2)
        sage: mat.elementary_divisors() == dhsw_snf(mat)
        True

    REFERENCES:

    - Dumas, Heckenbach, Saunders, Welker, "Computing simplicial
      homology based on efficient Smith normal form algorithms," in
      "Algebra, geometry, and software systems" (2003), 177-206.
    """
    ring = mat.base_ring()
    rows = mat.nrows()
    cols = mat.ncols()
    new_data = {}
    new_mat = matrix(ring, rows, cols, new_data)
    add_to_rank = 0
    zero_cols = 0
    if verbose:
        print "old matrix: %s by %s" % (rows, cols)
    # leading_positions: dictionary of lists indexed by row: if first
    # nonzero entry in column c is in row r, then leading_positions[r]
    # should contain c
    leading_positions = {}
    # pass 1:
    if verbose:
        print "starting pass 1"
    for j in range(cols):
        # new_col is a matrix with one column: sparse matrices seem to
        # be less buggy than sparse vectors (#5184, #5185), and
        # perhaps also faster.
        new_col = mat.matrix_from_columns([j])
        if new_col.is_zero():
            zero_cols += 1
        else:
            check_leading = True
            while check_leading:
                i = new_col.nonzero_positions_in_column(0)[0]
                entry = new_col[i,0]
                check_leading = False
                if i in leading_positions:
                    for c in leading_positions[i]:
                        earlier = new_mat[i,c]
                        # right now we don't check to see if entry divides
                        # earlier, because we don't want to modify the
                        # earlier columns of the matrix.  Deal with this
                        # in pass 2.
                        if entry and earlier.divides(entry):
                            quo = entry.divide_knowing_divisible_by(earlier)
                            new_col = new_col - quo * new_mat.matrix_from_columns([c])
                            entry = 0
                            if not new_col.is_zero():
                                check_leading = True
            if not new_col.is_zero():
                new_mat.set_column(j-zero_cols, new_col.column(0))
                i = new_col.nonzero_positions_in_column(0)[0]
                if i in leading_positions:
                    leading_positions[i].append(j-zero_cols)
                else:
                    leading_positions[i] = [j-zero_cols]
            else:
                zero_cols += 1
    # pass 2:
    # first eliminate the zero columns at the end
    cols = cols - zero_cols
    zero_cols = 0
    new_mat = new_mat.matrix_from_columns(range(cols))
    if verbose:
        print "starting pass 2"
    keep_columns = range(cols)
    check_leading = True
    while check_leading:
        check_leading = False
        new_leading = leading_positions.copy()
        for i in leading_positions:
            if len(leading_positions[i]) > 1:
                j = leading_positions[i][0]
                jth = new_mat[i, j]
                for n in leading_positions[i][1:]:
                    nth = new_mat[i,n]
                    if jth.divides(nth):
                        quo = nth.divide_knowing_divisible_by(jth)
                        new_mat.add_multiple_of_column(n, j, -quo)
                    elif nth.divides(jth):
                        quo = jth.divide_knowing_divisible_by(nth)
                        jth = nth
                        new_mat.swap_columns(n, j)
                        new_mat.add_multiple_of_column(n, j, -quo)
                    else:
                        (g,r,s) = jth.xgcd(nth)
                        (unit,A,B) = r.xgcd(-s)  # unit ought to be 1 here
                        jth_col = new_mat.column(j)
                        nth_col = new_mat.column(n)
                        new_mat.set_column(j, r*jth_col + s*nth_col)
                        new_mat.set_column(n, B*jth_col + A*nth_col)
                        nth = B*jth + A*nth
                        jth = g
                        # at this point, jth should divide nth
                        quo = nth.divide_knowing_divisible_by(jth)
                        new_mat.add_multiple_of_column(n, j, -quo)
                    new_leading[i].remove(n)
                    if new_mat.column(n).is_zero():
                        keep_columns.remove(n)
                        zero_cols += 1
                    else:
                        new_r = new_mat.column(n).nonzero_positions()[0]
                        if new_r in new_leading:
                            new_leading[new_r].append(n)
                        else:
                            new_leading[new_r] = [n]
                        check_leading = True
        leading_positions = new_leading
    # pass 3: get rid of columns which start with 1 or -1
    if verbose:
        print "starting pass 3"
    max_leading = 1
    for i in leading_positions:
        j = leading_positions[i][0]
        entry = new_mat[i,j]
        if entry.abs() == 1:
            add_to_rank += 1
            keep_columns.remove(j)
            for c in new_mat.nonzero_positions_in_row(i):
                if c in keep_columns:
                    new_mat.add_multiple_of_column(c, j, -entry * new_mat[i,c])
        else:
            max_leading = max(max_leading, new_mat[i,j].abs())
    # form the new matrix
    if max_leading != 1:
        new_mat = new_mat.matrix_from_columns(keep_columns)
        if verbose:
            print "new matrix: %s by %s" % (new_mat.nrows(), new_mat.ncols())
        if new_mat.is_sparse():
            ed = [1]*add_to_rank + new_mat.dense_matrix().elementary_divisors()
        else:
            ed = [1]*add_to_rank + new_mat.elementary_divisors()
    else:
        if verbose:
            print "new matrix: all pivots are 1 or -1"
        ed = [1]*add_to_rank

    if len(ed) < rows:
        return ed + [0]*(rows - len(ed))
    else:
        return ed[:rows]

def _latex_module(R, m):
    """
    LaTeX string representing a free module over ``R`` of rank ``m``.

    INPUT:

    - ``R`` - a commutative ring
    - ``m`` - non-negative integer

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

class ChainComplex(SageObject):
    """
    Define a chain complex.

    INPUT:

    -  ``data`` - the data defining the chain complex; see below for
       more details.

    -  ``base_ring`` - a commutative ring (optional), the ring over
       which the chain complex is defined. If this is not specified,
       it is determined by the data defining the chain complex.

    -  ``grading_group`` - a free abelian group (optional, default
       ZZ), the group over which the chain complex is indexed.

    -  ``degree`` - element of grading_group (optional, default 1),
       the degree of the differential.

    -  ``check_products`` - boolean (optional, default True).  If True,
       check that each consecutive pair of differentials are
       composable and have composite equal to zero.

    OUTPUT: a chain complex

    .. warning::

       Right now, homology calculations will only work if the base
       ring is either ZZ or a field, so please take this into account
       when defining a chain complex.

    Use data to define the chain complex.  This may be in any of the
    following forms.

    1. a dictionary with integers (or more generally, elements of
       grading_group) for keys, and with data[n] a matrix representing
       (via left multiplication) the differential coming from degree
       `n`.  (Note that the shape of the matrix then determines the
       rank of the free modules `C_n` and `C_{n+d}`.)

    2. a list or tuple of the form `[C_0, d_0, C_1, d_1, C_2, d_2,
       ...]`, where each `C_i` is a free module and each `d_i` is a
       matrix, as above.  This only makes sense if ``grading_group``
       is `\\ZZ` and ``degree`` is 1.

    3. a list or tuple of the form `[r_0, d_0, r_1, d_1, r_2, d_2,
       ...]`, where `r_i` is the rank of the free module `C_i` and
       each `d_i` is a matrix, as above.  This only makes sense if
       ``grading_group`` is `\\ZZ` and ``degree`` is 1.

    4. a list or tuple of the form `[d_0, d_1, d_2, ...]` where each
       `d_i` is a matrix, as above.  This only makes sense if
       ``grading_group`` is `\\ZZ` and ``degree`` is 1.

    .. note::

       In fact, the free modules `C_i` in case 2 and the ranks `r_i`
       in case 3 are ignored: only the matrices are kept, and from
       their shapes, the ranks of the modules are determined.
       (Indeed, if ``data`` is a list or tuple, then any element which
       is not a matrix is discarded; thus the list may have any number
       of different things in it, and all of the non-matrices will be
       ignored.)  No error checking is done to make sure, for
       instance, that the given modules have the appropriate ranks for
       the given matrices.  However, as long as ``check_products`` is
       True, the code checks to see if the matrices are composable and
       that each appropriate composite is zero.

    If the base ring is not specified, then the matrices are examined
    to determine a ring over which they are all naturally defined, and
    this becomes the base ring for the complex.  If no such ring can
    be found, an error is raised.  If the base ring is specified, then
    the matrices are converted automatically to this ring when
    defining the chain complex.  If some matrix cannot be converted,
    then an error is raised.

    EXAMPLES::

        sage: ChainComplex()
        Chain complex with at most 0 nonzero terms over Integer Ring
        sage: C = ChainComplex({0: matrix(ZZ, 2, 3, [3, 0, 0, 0, 0, 0])})
        sage: C
        Chain complex with at most 2 nonzero terms over Integer Ring
        sage: D = ChainComplex([matrix(ZZ, 2, 2, [0, 1, 0, 0]), matrix(ZZ, 2, 2, [0, 1, 0, 0])], base_ring=GF(2)); D
        Chain complex with at most 3 nonzero terms over Finite Field of size 2
        sage: D == loads(dumps(D))
        True

    Note that when a chain complex is defined in Sage, new
    differentials may be created: every nonzero module in the chain
    complex must have a differential coming from it, even if that
    differential is zero::

        sage: IZ = ChainComplex({0: identity_matrix(ZZ, 1)})
        sage: IZ.differential()  # the differentials in the chain complex
        {0: [1], 1: []}
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
        TypeError: Unable to coerce 0 (<type 'sage.rings.finite_rings.element_givaro.FiniteField_givaroElement'>) to Rational
    """

    def __init__(self, data=None, **kwds):
        """
        See ``ChainComplex`` for full documentation.

        EXAMPLES::

            sage: C = ChainComplex(); C
            Chain complex with at most 0 nonzero terms over Integer Ring
            sage: D = ChainComplex({0: matrix(ZZ, 2, 3, [3, 0, 0, 0, 0, 0])})
            sage: D
            Chain complex with at most 2 nonzero terms over Integer Ring
        """
        check_products = kwds.get('check_products', True) or kwds.get('check_diffs')
        base_ring = kwds.get('base_ring', None)
        grading_group = kwds.get('grading_group', ZZ)
        degree = kwds.get('degree', 1)

        try:
            deg = grading_group(degree)
        except StandardError:
            raise ValueError, "The 'degree' does not appear to be an element of the grading group."
        # check form of data
        new_data = {}
        if data is None or (isinstance(data, (list, tuple)) and len(data) == 0):
            # the zero chain complex
            try:
                zero = grading_group.zero_vector()
            except AttributeError:
                zero = ZZ.zero_element()
            if base_ring is None:
                base_ring = ZZ
            new_data = {zero: matrix(base_ring, 0, 0, [])}
        else:
            if isinstance(data, dict):  # case 1, dictionary
                temp_dict = data
            elif isinstance(data, (list, tuple)):  # cases 2, 3, 4: list or tuple
                if degree != 1:
                    raise ValueError, "The degree must be +1 if the data argument is a list or tuple."
                if grading_group != ZZ:
                    raise ValueError, "The grading_group must be ZZ if the data argument is a list or tuple."
                new_list = filter(lambda x: isinstance(x, Matrix), data)
                temp_dict = dict(zip(range(len(new_list)), new_list))
            else:
                raise TypeError, "The data for a chain complex must be a dictionary, list, or tuple."

            if base_ring is None:
                junk, ring = prepare_dict(dict([n, temp_dict[n].base_ring()(0)] for n in temp_dict))
                base_ring = ring
            # 'dict' is the dictionary of matrices.  check to see that
            # each entry is in fact a matrix, and that the codomain of
            # one matches up with the domain of the next.  if the
            # number of rows in the differential in degree n is
            # nonzero, then make sure that the differential in degree
            # n+degree is present (although it of course may be zero).
            # If it's not present, add a zero matrix of the
            # appropriate shape.  This way, if self._data does not
            # have a key for n, then the complex is zero in degree n.
            for n in temp_dict.keys():
                if not n in grading_group:
                    raise ValueError, "One of the dictionary keys does not appear to be an element of the grading group."
                mat = temp_dict[n]
                if not isinstance(mat, Matrix):
                    raise TypeError, "One of the differentials in the data dictionary indexed by does not appear to be a matrix."
                if mat.base_ring() == base_ring:
                    new_data[n] = mat
                else:
                    new_data[n] = mat.change_ring(base_ring)
            for n in temp_dict.keys():
                mat = temp_dict[n]
                if n+degree in temp_dict:
                    mat2 = temp_dict[n+degree]
                    if check_products:
                        try:
                            prod = mat2 * mat
                        except TypeError:
                            raise TypeError, "The differentials d_{%s} and d_{%s} are not compatible: their product is not defined." % (n, n+degree)
                        if not prod.is_zero():
                            raise ValueError, "The differentials d_{%s} and d_{%s} are not compatible: their composition is not zero." % (n, n+degree)
                else:
                    if not mat.nrows() == 0:
                        if n+2*degree in temp_dict:
                            new_data[n+degree] = matrix(base_ring, temp_dict[n+2*degree].ncols(), mat.nrows())
                        else:
                            new_data[n+degree] = matrix(base_ring, 0, mat.nrows())
        # here ends the initialization/error-checking of the data
        self._grading_group = grading_group
        self._degree = deg
        self._base_ring = base_ring
        self._diff = new_data
        # self._ranks: dictionary for caching the ranks of the
        # differentials: keys are pairs (dim, base_ring)
        self._ranks = {}

    def base_ring(self):
        """
        The base ring for this simplicial complex.

        EXAMPLES::

            sage: ChainComplex().base_ring()
            Integer Ring
        """
        return self._base_ring

    def differential(self, dim=None):
        """
        The differentials which make up the chain complex.

        INPUT:

        -  ``dim`` - element of the grading group (optional, default
           None).  If this is None, return a dictionary of all of the
           differentials.  If this is a single element, return the
           differential starting in that dimension.

        OUTPUT: either a dictionary of all of the differentials or a single
        differential (i.e., a matrix)

        EXAMPLES::

            sage: D = ChainComplex({0: matrix(ZZ, 2, 2, [1,0,0,2])})
            sage: D.differential()
            {0: [1 0]
            [0 2], 1: []}
            sage: D.differential(0)
            [1 0]
            [0 2]
            sage: C = ChainComplex({0: identity_matrix(ZZ, 40)})
            sage: C.differential()
            {0: 40 x 40 dense matrix over Integer Ring, 1: []}
        """
        if dim is None:
            return self._diff
        deg = self._degree
        if dim in self._diff:
            return self._diff[dim]
        if dim+deg in self._diff:
            return matrix(self._base_ring, self._diff[dim+deg].ncols(), 0)
        return matrix(self._base_ring, 0, 0)

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
            sage: C._degree
            1
            sage: C.differential(2)
            [3 0 0]
            [0 0 0]
            sage: C.dual()._degree
            -1
            sage: C.dual().differential(3)
            [3 0]
            [0 0]
            [0 0]
        """
        data = {}
        deg = self._degree
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
        return FreeModule(self._base_ring, rank)

    def __cmp__(self, other):
        """
        Return True iff this chain complex is the same as other: that
        is, if the base rings and the matrices of the two are the
        same.

        EXAMPLES::

            sage: C = ChainComplex({0: matrix(ZZ, 2, 3, [3, 0, 0, 0, 0, 0])}, base_ring=GF(2))
            sage: D = ChainComplex({0: matrix(GF(2), 2, 3, [1, 0, 0, 0, 0, 0]), 1: matrix(ZZ, 0, 2), 3: matrix(ZZ, 0, 0)})  # base_ring determined from the matrices
            sage: C == D
            True
        """
        if self._base_ring != other._base_ring:
            return -1
        R = self._base_ring
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
        """
        The homology of the chain complex in the given dimension.

        INPUT:

        -  ``dim`` - an element of the grading group for the chain
           complex (optional, default None): the degree in which to
           compute homology. If this is None, return the homology in
           every dimension in which the chain complex is possibly
           nonzero.

        -  ``base_ring`` - a commutative ring (optional, default is the
           base ring for the chain complex).  Must be either the
           integers `\\ZZ` or a field.

        -  ``generators`` - boolean (optional, default False).  If
           True, return generators for the homology groups along with
           the groups.  NOTE: this is only implemented if the CHomP
           package is available.

        -  ``verbose`` - boolean (optional, default False).  If True,
           print some messages as the homology is computed.

        -  ``algorithm`` - string (optional, default 'auto').  The
           options are 'auto', 'dhsw', 'pari' or 'no_chomp'.  See
           below for descriptions.

        OUTPUT: if dim is specified, the homology in dimension dim.
        Otherwise, the homology in every dimension as a dictionary
        indexed by dimension.

        ALGORITHM:

        If ``algorithm`` is set to 'auto' (the default), then use
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
        ``algorithm``.  If ``algorithm`` is 'auto' or 'no_chomp', then
        for each relatively small matrix, use the standard Sage
        method, which calls the Pari package.  For any large matrix,
        reduce it using the Dumas, Heckenbach, Saunders, and Welker
        elimination algorithm: see
        :func:`sage.homology.chain_complex.dhsw_snf` for details.

        Finally, ``algorithm`` may also be 'pari' or 'dhsw', which
        forces the named algorithm to be used regardless of the size
        of the matrices and regardless of whether CHomP is available.

        As of this writing, CHomP is by far the fastest option,
        followed by the 'auto' or 'no_chomp' setting of using the
        Dumas, Heckenbach, Saunders, and Welker elimination algorithm
        for large matrices and Pari for small ones.

        .. warning::

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

        Generators (only available via CHomP): generators are given as
        a list of cycles, each of which is an element in the
        appropriate free module, and hence is represented as a vector::

            sage: C.homology(1, generators=True)  # optional - CHomP
            (Z x C3, [(1, 0), (0, 1)])
        """
        from sage.interfaces.chomp import have_chomp, homchain

        algorithm = kwds.get('algorithm', 'auto')
        verbose = kwds.get('verbose', False)
        base_ring = kwds.get('base_ring', None)
        if base_ring is None or base_ring == self._base_ring:
            change_ring = False
            base_ring = self._base_ring
        else:
            change_ring = True
        if (not base_ring.is_field()) and (base_ring != ZZ):
            raise NotImplementedError, "Can't compute homology if the base ring is not the integers or a field."

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
                        field = base_ring.is_prime_field()
                        answer = {}
                        if isinstance(dim, (list, tuple)):
                            for d in dim:
                                if d in H:
                                    answer[d] = H[d]
                                else:
                                    if field:
                                        answer[d] = VectorSpace(base_ring, 0)
                                    else:
                                        answer[d] = HomologyGroup(0)
                        else:
                            if dim in H:
                                answer = H[dim]
                            else:
                                if field:
                                    answer = VectorSpace(base_ring, 0)
                                else:
                                    answer = HomologyGroup(0)
                    else:
                        answer = H
                    return answer
                else:
                    if verbose:
                        print "ran CHomP, but no output."

        if algorithm not in ['dhsw', 'pari', 'auto', 'no_chomp']:
            raise NotImplementedError, "algorithm not recognized"
        # if dim is None, return all of the homology groups
        if dim is None:
            answer = {}
            for n in self._diff.keys():
                if verbose:
                    print "Computing homology of the chain complex in dimension %s..." % n
                if base_ring == self._base_ring:
                    answer[n] = self.homology(n, verbose=verbose,
                                              algorithm=algorithm)
                else:
                    answer[n] = self.homology(n, base_ring=base_ring,
                                              verbose=verbose,
                                              algorithm=algorithm)
            return answer
        # now compute the homology in the given dimension
        degree = self._degree
        # check that dim is in grading_group
        if not dim in self._grading_group:
            raise ValueError, "The dimension does not appear to be an element of the grading group."
        if dim in self._diff:
            # d_out is the differential going out of degree dim,
            # d_in is the differential entering degree dim
            d_out_cols = self._diff[dim].ncols()
            d_out_rows = self._diff[dim].nrows()
            # avoid bugs in rank of sparse mod n matrices (#5099):
            if min(d_out_cols, d_out_rows) == 0:
                d_out_rank = 0
            if base_ring == ZZ:
                temp_ring = QQ
            else:
                temp_ring = base_ring
            if (dim, temp_ring) in self._ranks:
                d_out_rank = self._ranks[(dim, temp_ring)]
            else:
                if min(d_out_cols, d_out_rows) == 0:
                    d_out_rank = 0
                elif change_ring or base_ring == ZZ:
                    d_out_rank = self._diff[dim].change_ring(temp_ring).rank()
                else:
                    d_out_rank = self._diff[dim].rank()
                self._ranks[(dim, temp_ring)] = d_out_rank

            if dim-degree in self._diff:
                if change_ring:
                    d_in = self._diff[dim-degree].change_ring(base_ring)
                else:
                    d_in = self._diff[dim-degree]
                if base_ring.is_field():
                    # avoid bugs in rank of sparse mod n matrices (#5099):
                    if d_out_cols == 0:
                        null = 0
                    elif d_out_rows == 0:
                        null = d_out_cols
                    else:
                        null = d_out_cols - d_out_rank
                    if d_in.nrows() == 0 or d_in.ncols() == 0:
                        rk = 0
                    else:
                        if (dim-degree, temp_ring) in self._ranks:
                            rk = self._ranks[(dim-degree, temp_ring)]
                        else:
                            rk = d_in.rank()
                            self._ranks[(dim-degree, temp_ring)] = rk
                    answer = VectorSpace(base_ring, null - rk)
                elif base_ring == ZZ:
                    nullity = d_out_cols - d_out_rank
                    if d_in.ncols() == 0:
                        all_divs = [0] * nullity
                    else:
                        if algorithm == 'auto' or algorithm == 'no_chomp':
                            if ((d_in.ncols() > 300 and d_in.nrows() > 300)
                                or (min(d_in.ncols(), d_in.nrows()) > 100 and
                                    d_in.ncols() + d_in.nrows() > 600)):
                                algorithm = 'dhsw'
                            else:
                                algorithm = 'pari'
                        if algorithm == 'dhsw':
                            all_divs = dhsw_snf(d_in, verbose=verbose)
                        else:
                            if d_in.is_sparse():
                                all_divs = d_in.dense_matrix().elementary_divisors(algorithm)
                            else:
                                all_divs = d_in.elementary_divisors(algorithm)
                    all_divs = all_divs[:nullity]
                    # divisors equal to 1 produce trivial
                    # summands, so filter them out
                    divisors = filter(lambda x: x != 1, all_divs)
                    answer = HomologyGroup(len(divisors), divisors)
                else:
                    # This code is not in use: base ring isn't a field
                    # or ZZ
                    pass
            else: # no incoming differential: it's zero
                if base_ring.is_field():
                    answer = VectorSpace(base_ring, d_out_cols - d_out_rank)
                elif base_ring == ZZ:
                    nullity = d_out_cols - d_out_rank
                    answer = HomologyGroup(nullity)
                else:
                    # This code is not in use: base ring isn't a field
                    # or ZZ
                    pass
        else:  # chain complex is zero here, so return the zero module
            if base_ring.is_field():
                answer = VectorSpace(base_ring, 0)
            elif base_ring == ZZ:
                answer = HomologyGroup(0)
            else:
                # This code is not in use: base ring isn't a field
                # or ZZ
                answer = FreeModule(self._base_ring, rank=0)
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

        -  ``dim`` - an element of the grading group for the chain
           complex or None (optional, default None).  If None, then
           return every Betti number, as a dictionary indexed by
           degree.  If an element of the grading group, then return
           the Betti number in that dimension.

        -  ``base_ring`` - a commutative ring (optional, default is the
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
            base_ring = self._base_ring
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
        """
        Look for torsion in this chain complex by computing its mod `p`
        homology for a range of primes `p`.

        INPUT:

        -  ``max_prime`` - prime number: search for torsion mod `p` for
           all `p` strictly less than this number.

        -  ``min_prime`` - prime (optional, default 2): search for
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
        `\\GF{P}`.  Then compute the mod `p` homology for a range of
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
        if self._base_ring != ZZ:
            raise ValueError, "This only works if the base ring of the chain complex is the integers"
        else:
            answer = []
            torsion_free = self.betti(base_ring=GF(max_prime))
            for p in prime_range(min_prime, max_prime):
                mod_p_betti = self.betti(base_ring=GF(p))
                if mod_p_betti != torsion_free:
                    diff_dict = {}
                    temp_diff = {}
                    D = self._degree
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

    def category(self):
        """
        Return the category to which this chain complex belongs: the
        category of all chain complexes over the base ring.

        EXAMPLES::

            sage: C = ChainComplex({0: matrix(ZZ, 2, 3, [3, 0, 0, 0, 0, 0])}, base_ring=GF(7))
            sage: C.category()
            Category of chain complexes over Finite Field of size 7
        """
        import sage.categories.all
        return sage.categories.all.ChainComplexes(self.base_ring())

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
        degree `-n`), thus turning a chain complex into a cochain
        complex without changing the homology (except for flipping it,
        too).

        EXAMPLES::

            sage: C = ChainComplex({2: matrix(ZZ, 2, 3, [3, 0, 0, 0, 0, 0])})
            sage: C._degree
            1
            sage: C.differential(2)
            [3 0 0]
            [0 0 0]
            sage: C._flip_()._degree
            -1
            sage: C._flip_().differential(-2)
            [3 0 0]
            [0 0 0]
        """
        data = {}
        deg = self._degree
        for d in self.differential():
            data[-d] = self.differential()[d]
        return ChainComplex(data, degree=-deg)

    def _chomp_repr_(self):
        r"""
        String representation of self suitable for use by the CHomP
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
        if (self._grading_group != ZZ or
            (self._degree != 1 and self._degree != -1)):
            raise ValueError, "CHomP only works on Z-graded chain complexes with differential of degree 1 or -1."
        base_ring = self._base_ring
        if (base_ring == QQ) or (base_ring != ZZ and not (base_ring.is_prime_field())):
            raise ValueError, "CHomP doesn't compute over the rationals, only over Z or F_p."
        if self._degree == -1:
            diffs = self.differential()
        else:
            diffs = self._flip_().differential()

        maxdim = max(diffs)
        mindim = min(diffs)
        # will shift chain complex by subtracting mindim from
        # dimensions, so its bottom dimension is zero.
        s = "chain complex\n\nmax dimension = %s\n\n" % (maxdim - mindim,)

        for i in range(0, maxdim - mindim + 1):
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
        string1 = "Chain complex with at most"
        string2 = " %s nonzero terms over %s" % (len(diffs),
                                                  self._base_ring)
        return string1 + string2

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
        deg = self._degree
        ring = self._base_ring
        if self._grading_group != ZZ:
            guess = dict.keys()[0]
            if guess-deg in dict:
                string += "\\dots \\xrightarrow{d_{%s}} " % latex(guess-deg)
            string += _latex_module(ring, mat.ncols())
            string += " \\xrightarrow{d_{%s}} \\dots" % latex(guess)
        else:
            backwards = (deg < 0)
            sorted_list = sorted(dict.keys(), reverse=backwards)
            if len(dict) <= 6:
                for n in sorted_list[:-1]:
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

class HomologyGroup_class(AdditiveAbelianGroup_fixed_gens):
    """
    Abelian group on `n` generators. This class inherits from
    ``AdditiveAbelianGroup``; see that for more documentation. The main
    difference between the classes is in the print representation.

    EXAMPLES::

        sage: from sage.homology.chain_complex import HomologyGroup
        sage: G = AbelianGroup(5,[5,5,7,8,9]); G
        Multiplicative Abelian Group isomorphic to C5 x C5 x C7 x C8 x C9
        sage: H = HomologyGroup(5,[5,5,7,8,9]); H
        C5 x C5 x C7 x C8 x C9
        sage: G == loads(dumps(G))
        True
        sage: AbelianGroup(4)
        Multiplicative Abelian Group isomorphic to Z x Z x Z x Z
        sage: HomologyGroup(4)
        Z x Z x Z x Z
        sage: HomologyGroup(100)
        Z^100
    """
    def __init__(self, n, invfac):
        """
        See ``HomologyGroup`` for full documentation.

        EXAMPLES::

            sage: from sage.homology.chain_complex import HomologyGroup
            sage: H = HomologyGroup(5,[5,5,7,8,9]); H
            C5 x C5 x C7 x C8 x C9
        """
        n = len(invfac)
        A = ZZ**n
        B = A.span([A.gen(i) * invfac[i] for i in xrange(n)])

        AdditiveAbelianGroup_fixed_gens.__init__(self, A, B, A.gens())
        self._original_invts = invfac

    def _repr_(self):
        """
        Print representation

        EXAMPLES::

            sage: from sage.homology.chain_complex import HomologyGroup
            sage: H = HomologyGroup(7,[4,4,4,4,4,7,7])
            sage: H._repr_()
            'C4^5 x C7 x C7'
            sage: HomologyGroup(6)
            Z^6
        """
        eldv = self._original_invts
        if len(eldv) == 0:
            return "0"
        rank = len(filter(lambda x: x == 0, eldv))
        torsion = sorted(filter(lambda x: x, eldv))
        if rank > 4:
            g = ["Z^%s" % rank]
        else:
            g = ["Z"] * rank
        if len(torsion) != 0:
            printed = []
            for t in torsion:
                numfac = torsion.count(t)
                too_many = (numfac > 4)
                if too_many:
                    if t not in printed:
                        g.append("C%s^%s" % (t, numfac))
                        printed.append(t)
                else:
                    g.append("C%s" % t)
        times = " x "
        return times.join(g)

    def _latex_(self):
        """
        LaTeX representation

        EXAMPLES::

            sage: from sage.homology.chain_complex import HomologyGroup
            sage: H = HomologyGroup(7,[4,4,4,4,4,7,7])
            sage: H._latex_()
            'C_{4}^{5} \\times C_{7} \\times C_{7}'
            sage: latex(HomologyGroup(6))
            \ZZ^{6}
        """
        eldv = self._original_invts
        if len(eldv) == 0:
            return "0"
        rank = len(filter(lambda x: x == 0, eldv))
        torsion = sorted(filter(lambda x: x, eldv))
        if rank > 4:
            g = ["\\ZZ^{%s}" % rank]
        else:
            g = ["\\ZZ"] * rank
        if len(torsion) != 0:
            printed = []
            for t in torsion:
                numfac = torsion.count(t)
                too_many = (numfac > 4)
                if too_many:
                    if t not in printed:
                        g.append("C_{%s}^{%s}" % (t, numfac))
                        printed.append(t)
                else:
                    g.append("C_{%s}" % t)
        times = " \\times "
        return times.join(g)

def HomologyGroup(n, invfac=None):
    """
    Abelian group on `n` generators. This class inherits from
    ``AdditiveAbelianGroup``; see that for more documentation.  The main
    difference between the classes is in the print representation.

    EXAMPLES::

        sage: from sage.homology.chain_complex import HomologyGroup
        sage: G = AbelianGroup(5,[5,5,7,8,9]); G
        Multiplicative Abelian Group isomorphic to C5 x C5 x C7 x C8 x C9
        sage: H = HomologyGroup(5,[5,5,7,8,9]); H
        C5 x C5 x C7 x C8 x C9
        sage: AbelianGroup(4)
        Multiplicative Abelian Group isomorphic to Z x Z x Z x Z
        sage: HomologyGroup(4)
        Z x Z x Z x Z
        sage: HomologyGroup(100)
        Z^100
    """
    # copied from AbelianGroup:
    if invfac is None:
        if isinstance(n, (list, tuple)):
            invfac = n
            n = len(n)
        else:
            invfac = []
    if len(invfac) < n:
        invfac = [0] * (n - len(invfac)) + invfac
    elif len(invfac) > n:
        raise ValueError, "invfac (=%s) must have length n (=%s)"%(invfac, n)
    M = HomologyGroup_class(n, invfac)
    return M
