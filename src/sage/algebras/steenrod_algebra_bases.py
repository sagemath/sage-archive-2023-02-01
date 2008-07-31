r"""
Steenrod algebra bases

AUTHORS:
    - John H. Palmieri (2008-07-30: version 0.9)

This package defines functions for computing various bases of the
Steenrod algebra, and for converting between the Milnor basis and any
other basis.

This packages implements a number of different bases, at least at the
prime 2.  The Milnor and Serre-Cartan bases are the most familiar and
most standard ones, and all of the others are defined in terms of one
of these.  The bases are described in the documentation for the
function \code{steenrod_algebra_basis}; also see the papers by Monks [M]
and Wood [W] for more information about them.  For commutator bases,
see the preprint by Palmieri and Zhang [PZ].

  * 'milnor': Milnor basis.

  * 'serre-cartan' or 'adem' or 'admissible': Serre-Cartan basis.

The other bases are as follows; these are only defined when $p=2$:

  * 'wood_y': Wood's Y basis.

  * 'wood_z': Wood's Z basis.

  * 'wall', 'wall_long': Wall's basis.

  * 'arnon_a', 'arnon_a_long': Arnon's A basis.

  * 'arnon_c': Arnon's C basis.

  * 'pst', 'pst_rlex', 'pst_llex', 'pst_deg', 'pst_revz': various
     $P^s_t$-bases.

  * 'comm', 'comm_rlex', 'comm_llex', 'comm_deg', 'comm_revz', or
     these with '_long' appended: various commutator bases.

EXAMPLES:
    sage: steenrod_algebra_basis(7,'milnor')
    (Sq(0,0,1), Sq(1,2), Sq(4,1), Sq(7))
    sage: steenrod_algebra_basis(5)   # milnor basis is the default
    (Sq(2,1), Sq(5))

The third (optional) argument to \code{steenrod_algebra_basis} is the prime $p$:
    sage: steenrod_algebra_basis(9, 'milnor', p=3)
    (Q_1 P(1), Q_0 P(2))
    sage: steenrod_algebra_basis(9, 'milnor', 3)
    (Q_1 P(1), Q_0 P(2))
    sage: steenrod_algebra_basis(17, 'milnor', 3)
    (Q_2, Q_1 P(3), Q_0 P(0,1), Q_0 P(4))

Other bases:
    sage: steenrod_algebra_basis(7,'admissible')
    (Sq^{7}, Sq^{6} Sq^{1}, Sq^{4} Sq^{2} Sq^{1}, Sq^{5} Sq^{2})
    sage: [x.basis('milnor') for x in steenrod_algebra_basis(7,'admissible')]
    [Sq(7),
    Sq(4,1) + Sq(7),
    Sq(0,0,1) + Sq(1,2) + Sq(4,1) + Sq(7),
    Sq(1,2) + Sq(7)]
    sage: Aw = SteenrodAlgebra(2, basis = 'wall_long')
    sage: [Aw(x) for x in steenrod_algebra_basis(7,'admissible')]
    [Sq^{1} Sq^{2} Sq^{4},
    Sq^{2} Sq^{4} Sq^{1},
    Sq^{4} Sq^{2} Sq^{1},
    Sq^{4} Sq^{1} Sq^{2}]
    sage: steenrod_algebra_basis(13,'admissible',p=3)
    (beta P^{3}, P^{3} beta)
    sage: steenrod_algebra_basis(5,'wall')
    (Q^{2}_{2} Q^{0}_{0}, Q^{1}_{1} Q^{1}_{0})
    sage: steenrod_algebra_basis(5,'wall_long')
    (Sq^{4} Sq^{1}, Sq^{2} Sq^{1} Sq^{2})
    sage: steenrod_algebra_basis(5,'pst-rlex')
    (P^{0}_{1} P^{2}_{1}, P^{1}_{1} P^{0}_{2})

This file also contains a function \code{milnor_convert} which converts
elements from the (default) Milnor basis representation to a
representation in another basis.  The output is a dictionary which
gives the new representation; its form depends on the chosen basis.
For example, in the basis of admissible sequences (a.k.a. the
Serre-Cartan basis), each basis element is of the form $\text{Sq}^a
\text{Sq}^b ...$, and so is represented by a tuple $(a,b,...)$ of
integers.  Thus the dictionary has such tuples as keys, with the
coefficient of the basis element as the associated value:
    sage: from sage.algebras.steenrod_algebra_bases import milnor_convert
    sage: milnor_convert(Sq(2)*Sq(4) + Sq(2)*Sq(5), 'admissible')
    {(5, 1): 1, (6, 1): 1, (6,): 1}
    sage: milnor_convert(Sq(2)*Sq(4) + Sq(2)*Sq(5), 'pst')
    {((1, 1), (2, 1)): 1, ((0, 1), (1, 1), (2, 1)): 1, ((0, 2), (2, 1)): 1}

Users shouldn't need to call \code{milnor_convert}; they should use the
\code{basis} method to view a single element in another basis, or define a
Steenrod algebra with a different default basis and work in that
algebra:
    sage: x = Sq(2)*Sq(4) + Sq(2)*Sq(5)
    sage: x
    Sq(3,1) + Sq(4,1) + Sq(6) + Sq(7)
    sage: x.basis('milnor')  # 'milnor' is the default basis
    Sq(3,1) + Sq(4,1) + Sq(6) + Sq(7)
    sage: x.basis('adem')
    Sq^{5} Sq^{1} + Sq^{6} + Sq^{6} Sq^{1}
    sage: x.basis('pst')
    P^{0}_{1} P^{1}_{1} P^{2}_{1} + P^{0}_{2} P^{2}_{1} + P^{1}_{1} P^{2}_{1}
    sage: A = SteenrodAlgebra(2, basis='pst')
    sage: A(Sq(2) * Sq(4) + Sq(2) * Sq(5))
    P^{0}_{1} P^{1}_{1} P^{2}_{1} + P^{0}_{2} P^{2}_{1} + P^{1}_{1} P^{2}_{1}

**************

INTERNAL DOCUMENTATION:

If you want to implement a new basis for the Steenrod algebra:

    In the file 'steenrod_algebra.py':

        add functionality to \code{get_basis_name}: this should accept as
        input various synonyms for the basis, and its output should be
        an element of \code{_steenrod_basis_unique_names} (see the next
        file).

    In the file 'steenrod_algebra_element.py':

        add name of basis to \code{_steenrod_basis_unique_names}

        add functionality to \code{string_rep}, which probably involves
        adding a \code{BASIS_mono_to_string} function

        add functionality to the \code{_basis_dictionary} method

    In this file ('steenrod_algebra_bases.py'):

        add appropriate lines to \code{steenrod_algebra_basis}

        add a function to compute the basis in a given dimension (to be
        called by \code{steenrod_algebra_basis})

REFERENCES:

    [M]  K. G. Monks, "Change of basis, monomial relations, and $P^s_t$
         bases for the Steenrod algebra," J. Pure Appl. Algebra 125 (1998),
         no. 1-3, 235--260.

    [PZ] J. H. Palmieri and J. J. Zhang, "Commutators in the Steenrod
         algebra," preprint (2008)

    [W]  R. M. W. Wood, "Problems in the Steenrod algebra," Bull. London
         Math. Soc. 30 (1998), no. 5, 449--517.
"""

#*****************************************************************************
#       Copyright (C) 2008 John H. Palmieri <palmieri@math.washington.edu>
#  Distributed under the terms of the GNU General Public License (GPL)
#*****************************************************************************

from sage.rings.integer import Integer
from sage.rings.all import GF

from steenrod_algebra import SteenrodAlgebra, get_basis_name, \
    _steenrod_serre_cartan_basis_names, _steenrod_milnor_basis_names
from steenrod_algebra_element import SteenrodAlgebraElement, Sq, pst


# cache matrices from convert_to_milnor_matrix here:
_conversion_matrices = {}

def convert_to_milnor_matrix(n, basis, p=2):
    r"""
    Change-of-basis matrix, 'basis' to Milnor, in dimension $n$, at the prime $p$.

    INPUT:
        n -- non-negative integer, the dimension
        basis -- string, the basis from which to convert
        p -- positive prime number (optional, default 2)

    OUTPUT:
        matrix -- change-of-basis matrix, a square matrix over GF(p)

    (This is not really intended for casual users, so no error
    checking is made on the integer $n$, the basis name, or the prime.)

    EXAMPLES:
        sage: from sage.algebras.steenrod_algebra_bases import convert_to_milnor_matrix
        sage: convert_to_milnor_matrix(5, 'adem')
        [0 1]
        [1 1]
        sage: convert_to_milnor_matrix(45, 'milnor')
        111 x 111 dense matrix over Finite Field of size 2
        sage: convert_to_milnor_matrix(12,'wall')
        [1 0 0 1 0 0 0]
        [1 1 0 0 0 1 0]
        [0 1 0 1 0 0 0]
        [0 0 0 1 0 0 0]
        [1 1 0 0 1 0 0]
        [0 0 1 1 1 0 1]
        [0 0 0 0 1 0 1]

    The function takes an optional argument, the prime $p$ over which to work:
        sage: convert_to_milnor_matrix(17,'adem',3)
        [0 0 1 1]
        [0 0 0 1]
        [1 1 1 1]
        [0 1 0 1]
        sage: convert_to_milnor_matrix(48,'adem',5)
        [0 1]
        [1 1]
        sage: convert_to_milnor_matrix(36,'adem',3)
        [0 0 1]
        [0 1 0]
        [1 2 0]
    """
    if _conversion_matrices.has_key((n, basis, p)):
        return _conversion_matrices[(n, basis, p)]
    from sage.matrix.constructor import matrix
    milnor_base = steenrod_algebra_basis(n,'milnor',p)
    rows = []
    for poly in steenrod_algebra_basis(n,basis,p):
        for v in milnor_base:
            for m in v._basis_dictionary('milnor'):
                key = m
            if poly._raw['milnor'].has_key(key):
                rows = rows + [poly._raw['milnor'][key]]
            else:
                rows = rows + [0]
    d = len(milnor_base)
    _conversion_matrices[(n, basis, p)] = matrix(GF(p),d,d,rows)
    return matrix(GF(p),d,d,rows)


def convert_from_milnor_matrix(n, basis, p=2):
    r"""
    Change-of-basis matrix, Milnor to 'basis', in dimension $n$.

    INPUT:
        n -- non-negative integer, the dimension
        basis -- string, the basis to which to convert
        p -- positive prime number (optional, default 2)

    OUTPUT:
        matrix -- change-of-basis matrix, a square matrix over GF(p)

    (This is not really intended for casual users, so no error
    checking is made on the integer $n$, the basis name, or the prime.)

    EXAMPLES:
        sage: from sage.algebras.steenrod_algebra_bases import convert_from_milnor_matrix, convert_to_milnor_matrix
        sage: convert_from_milnor_matrix(12,'wall')
        [1 0 0 1 0 0 0]
        [0 0 1 1 0 0 0]
        [0 0 0 1 0 1 1]
        [0 0 0 1 0 0 0]
        [1 0 1 0 1 0 0]
        [1 1 1 0 0 0 0]
        [1 0 1 0 1 0 1]
        sage: convert_from_milnor_matrix(38,'serre_cartan')
        72 x 72 dense matrix over Finite Field of size 2
        sage: x = convert_to_milnor_matrix(20,'wood_y')
        sage: y = convert_from_milnor_matrix(20,'wood_y')
        sage: x*y
        [1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0]
        [0 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0]
        [0 0 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0]
        [0 0 0 1 0 0 0 0 0 0 0 0 0 0 0 0 0]
        [0 0 0 0 1 0 0 0 0 0 0 0 0 0 0 0 0]
        [0 0 0 0 0 1 0 0 0 0 0 0 0 0 0 0 0]
        [0 0 0 0 0 0 1 0 0 0 0 0 0 0 0 0 0]
        [0 0 0 0 0 0 0 1 0 0 0 0 0 0 0 0 0]
        [0 0 0 0 0 0 0 0 1 0 0 0 0 0 0 0 0]
        [0 0 0 0 0 0 0 0 0 1 0 0 0 0 0 0 0]
        [0 0 0 0 0 0 0 0 0 0 1 0 0 0 0 0 0]
        [0 0 0 0 0 0 0 0 0 0 0 1 0 0 0 0 0]
        [0 0 0 0 0 0 0 0 0 0 0 0 1 0 0 0 0]
        [0 0 0 0 0 0 0 0 0 0 0 0 0 1 0 0 0]
        [0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 0 0]
        [0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 0]
        [0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1]

    The function takes an optional argument, the prime $p$ over which to work:
        sage: convert_from_milnor_matrix(17,'adem',3)
        [2 1 1 2]
        [0 2 0 1]
        [1 2 0 0]
        [0 1 0 0]
    """
    mat = convert_to_milnor_matrix(n,basis,p)
    if mat.nrows() != 0:
        return convert_to_milnor_matrix(n,basis,p).inverse()
    else:
        return mat


def make_elt_homogeneous(poly):
    """
    Break element of the Steenrod algebra into a list of homogeneous pieces.

    INPUT:
        poly -- an element of the Steenrod algebra

    OUTPUT:
        list of homogeneous elements of the Steenrod algebra whose sum is poly

    EXAMPLES:
        sage: from sage.algebras.steenrod_algebra_bases import make_elt_homogeneous
        sage: make_elt_homogeneous(Sq(2)*Sq(4) + Sq(2)*Sq(5))
        [Sq(3,1) + Sq(6), Sq(4,1) + Sq(7)]
    """
    degrees = set([m.degree() for m in poly])   # set of degrees in poly
    result = []
    for d in degrees:
        result_d = 0
        for m in poly:
            if m.degree() == d:
                result_d = result_d + m
        result.append(result_d)
    return result


def milnor_convert(poly, basis):
    r"""
    Convert an element of the Steenrod algebra in the Milnor basis
    to its representation in the chosen basis.

    INPUT:
        poly -- element of the Steenrod algebra
        basis -- basis to which to convert

    OUTPUT:
        dict -- dictionary


    This returns a dictionary of terms of the form (mono: coeff),
    where mono is a monomial in 'basis'.  The form of mono depends on
    the chosen basis.

    EXAMPLES:
        sage: from sage.algebras.steenrod_algebra_bases import milnor_convert
        sage: milnor_convert(Sq(2)*Sq(4) + Sq(2)*Sq(5), 'adem')
        {(5, 1): 1, (6, 1): 1, (6,): 1}
        sage: A3 = SteenrodAlgebra(3)
        sage: a = A3.Q(1) * A3.P(2,2); a
        Q_1 P(2,2)
        sage: milnor_convert(a, 'adem')
        {(0, 9, 1, 2, 0): 1, (1, 9, 0, 2, 0): 2}
        sage: milnor_convert(2 * a, 'adem')
        {(0, 9, 1, 2, 0): 2, (1, 9, 0, 2, 0): 1}
        sage: (Sq(2)*Sq(4) + Sq(2)*Sq(5)).basis('adem')
        Sq^{5} Sq^{1} + Sq^{6} + Sq^{6} Sq^{1}
        sage: a.basis('adem')
        P^{9} beta P^{2} + 2 beta P^{9} P^{2}
        sage: milnor_convert(Sq(2)*Sq(4) + Sq(2)*Sq(5), 'pst')
        {((1, 1), (2, 1)): 1, ((0, 1), (1, 1), (2, 1)): 1, ((0, 2), (2, 1)): 1}
        sage: (Sq(2)*Sq(4) + Sq(2)*Sq(5)).basis('pst')
        P^{0}_{1} P^{1}_{1} P^{2}_{1} + P^{0}_{2} P^{2}_{1} + P^{1}_{1}
        P^{2}_{1}
    """
    from sage.matrix.constructor import matrix
    p = poly._prime
    basis_name = get_basis_name(basis, p)
    if basis_name.find('long') >= 0:
        basis_name = basis_name.rsplit('_', 1)[0]
    if basis_name == 'milnor': return poly._raw['milnor']
    dict = {}
    for x in make_elt_homogeneous(poly):
        dim = x.degree()
        vec = []
        for mono in steenrod_algebra_basis(dim,'milnor',p):
            if mono._raw['milnor'].keys()[0] in x._raw['milnor'].keys():
                for entry in mono._raw['milnor']:
                    coeff =  x._raw['milnor'][entry]
                vec = vec + [coeff]
            else:
                vec = vec + [0]
        new_vec = (matrix(GF(p),1,len(vec),vec)
                   * convert_from_milnor_matrix(dim, basis_name, p)).row(0)
        for (mono,coeff) in zip(steenrod_algebra_basis(dim, basis_name, p), new_vec):
            if coeff != 0:
                old_dict = mono._raw[basis_name]
                for entry in old_dict:
                    dict[entry] = coeff
    return dict

def steenrod_algebra_basis(n, basis='milnor', p=2):
    r"""
    Basis for the Steenrod algebra in degree $n$.

    INPUT:
        n -- non-negative integer
        basis -- string, which basis to use  (optional, default = 'milnor')
        p -- positive prime number (optional, default = 2)

    OUTPUT:
        tuple of basis elements for the Steenrod algebra in dimension n

    The choices for the string basis are as follows:

      * 'milnor': Milnor basis.  When $p=2$, the Milnor basis consists
        of symbols of the form $\text{Sq}(m_1, m_2, ..., m_t)$, where each
        $m_i$ is a non-negative integer and if $t>1$, then the last
        entry $m_t > 0$.  When $p$ is odd, the Milnor basis consists
        of symbols of the form $Q_{e_1} Q_{e_2} ... \mathcal{P}(m_1,
        m_2, ..., m_t)$, where $0 \leq e_1 < e_2 < ...$, each $m_i$ is
        a non-negative integer, and if $t>1$, then the last entry $m_t
        > 0$.

      * 'serre-cartan' or 'adem' or 'admissible': Serre-Cartan basis.
        The Serre-Cartan basis consists of 'admissible monomials' in
        the Steenrod operations.  Thus at the prime 2, it consists of
        monomials $\text{Sq}^{m_1} \text{Sq}^{m_2} ... \text{Sq}^{m_t}$
        with $m_i \geq 2m_{i+1}$ for each $i$.  At odd primes, it
        consists of monomials $\beta^{\epsilon_0} \mathcal{P}^{s_1}
        \beta^{\epsilon_1} \mathcal{P}^{s_2} ...  \mathcal{P}^{s_k}
        \beta^{\epsilon_k}$ with each $\epsilon_i$ either 0 or 1, $s_i
        \geq p s_{i+1} + \epsilon_i$, and $s_k \geq 1$.

        When $p=2$, the element $\text{Sq}^a$ equals the Milnor
        element $\text{Sq}(a)$; when $p$ is odd, $\mathcal{P}^a =
        \mathcal{P}(a)$ and $\beta = Q_0$.  Hence for any Serre-Cartan
        basis element, one can represent it in the Milnor basis by
        computing an appropriate product using Milnor multiplication.

    The rest of these bases are only defined when $p=2$.

      * 'wood_y': Wood's Y basis.  For pairs of non-negative integers
        $(m,k)$, let $w(m,k) = \text{Sq}^{2^m (2^{k+1}-1)}$.  Wood's $Y$
        basis consists of monomials $w(m_0,k_0) ... w(m_t, k_t)$ with
        $(m_i,k_i) > (m_{i+1},k_{i+1})$, in left lex order.

      * 'wood_z': Wood's Z basis.  For pairs of non-negative integers
        $(m,k)$, let $w(m,k) = \text{Sq}^{2^m (2^{k+1}-1)}$.  Wood's $Z$
        basis consists of monomials $w(m_0,k_0) ... w(m_t, k_t)$ with
        $(m_i+k_i,m_i) > (m_{i+1}+k_{i+1},m_{i+1})$, in left lex
        order.

      * 'wall' or 'wall_long': Wall's basis.  For any pair of integers
        $(m,k)$ with $m \geq k \geq 0$, let $Q^m_k = \text{Sq}^{2^k}
        \text{Sq}^{2^{k+1}} ... \text{Sq}^{2^m}$.  The elements of
        Wall's basis are monomials $Q^{m_0}_{k_0} ... Q^{m_t}_{k_t}$
        with $(m_i, k_i) > (m_{i+1}, k_{i+1})$, ordered left
        lexicographically.

        (Note that $Q^m_k$ is the reverse of the element $X^m_k$ used in
        defining Arnon's A basis.)

        The standard way of printing elements of the Wall basis is to
        write elements in terms of the $Q^m_k$.  If one sets the basis
        to 'wall_long' instead of 'wall', then each $Q^m_k$ is
        expanded as a product of factors of the form $\text{Sq}^{2^i}$.

      * 'arnon_a' or 'arnon_a_long': Arnon's A basis.  For any pair of
        integers $(m,k)$ with $m \geq k \geq 0$, let $X^m_k =
        \text{Sq}^{2^m} \text{Sq}^{2^{m-1}} ... \text{Sq}^{2^k}$.  The
        elements of Arnon's A basis are monomials $X^{m_0}_{k_0}
        ... X^{m_t}_{k_t}$ with $(m_i, k_i) < (m_{i+1}, k_{i+1})$,
        ordered left lexicographically.

        (Note that $X^m_k$ is the reverse of the element $Q^m_k$ used
        in defining Wall's basis.)

        The standard way of printing elements of Arnon's A basis is to
        write elements in terms of the $X^m_k$.  If one sets the basis
        to 'arnon_a_long' instead of 'arnon_a', then each $X^m_k$ is
        expanded as a product of factors of the form $\text{Sq}^{2^i}$.

      * 'arnon_c': Arnon's C basis.  The elements of Arnon's C basis
        are monomials of the form $\text{Sq}^{t_1} ... \text{Sq}^{t_m}$
        where for each $i$, we have $t_i \leq 2t_{i+1}$ and $2^i | t_{m-i}$.

      * 'pst', 'pst_rlex', 'pst_llex', 'pst_deg', 'pst_revz': various
        $P^s_t$-bases.  For integers $s \geq 0$ and $t > 0$, the element
        $P^s_t$ is the Milnor basis element $\text{Sq}(0, ..., 0, 2^s, 0,
        ...)$, with the nonzero entry in position $t$.  To obtain a
        $P^s_t$-basis, for each set $\{P^{s_1}_{t_1}, ...,
        P^{s_k}_{t_k}\}$ of (distinct) $P^s_t$'s, one chooses an
        ordering and forms the resulting monomial.  The set of all
        such monomials then forms a basis, and so one gets a basis by
        choosing an ordering on each monomial.

        The strings 'rlex', 'llex', etc., correspond to the following
        orderings.  These are all 'global' -- they give a global
        ordering on the $P^s_t$'s, not different orderings depending
        on the monomial.  They order the $P^s_t$'s using the pair of
        integers $(s,t)$ as follows:

          * 'rlex': right lexicographic ordering

          * 'llex': left lexicographic ordering

          * 'deg': ordered by degree, which is the same as left
            lexicographic ordering on the pair $(s+t,t)$

          * 'revz': left lexicographic ordering on the pair $(s+t,s)$,
            which is the reverse of the ordering used (on elements in
            the same degrees as the $P^s_t$'s) in Wood's Z basis:
            'revz' stands for 'reversed Z'.  This is the default:
            'pst' is the same as 'pst_revz'.

      * 'comm', 'comm_rlex', 'comm_llex', 'comm_deg', 'comm_revz', or
        any of these with '_long' appended: various commutator bases.
        Let $c_{i,1} = \text{Sq}^{2^i}$, let $c_{i,2} = [c_{i,1},
        c_{i+1,1}]$, and inductively define $c_{i,k} = [c_{i,k-1},
        c_{i+k-1,1}]$.  Thus $c_{i,k}$ is a $k$-fold iterated
        commutator of the elements $\text{Sq}^{2^i}$, ...,
        $\text{Sq}^{2^{i+k-1}}$.  Note that $\dim c_{i,k} = \dim
        P^i_k$.

        To obtain a commutator basis, for each set $\{c_{s_1,t_1},
        ..., c_{s_k,t_k}\}$ of (distinct) $c_{s,t}$'s, one chooses an
        ordering and forms the resulting monomial.  The set of all
        such monomials then forms a basis, and so one gets a basis by
        choosing an ordering on each monomial.  The strings 'rlex',
        etc., have the same meaning as for the orderings on
        $P^s_t$-bases.  As with the $P^s_t$-bases, 'comm_revz' is the
        default: 'comm' means 'comm_revz'.

        The commutator bases have alternative representations,
        obtained by appending 'long' to their names: instead of, say,
        $c_{2,2}$, the representation is $s_{48}$, indicating the
        commutator of $\text{Sq}^4$ and $\text{Sq}^8$, and $c_{0,3}$,
        which is equal to $[[\text{Sq}^1, \text{Sq}^2], \text{Sq}^4]$,
        is written as $s_{124}$.

    EXAMPLES:
        sage: steenrod_algebra_basis(7,'milnor')
        (Sq(0,0,1), Sq(1,2), Sq(4,1), Sq(7))
        sage: steenrod_algebra_basis(5)   # milnor basis is the default
        (Sq(2,1), Sq(5))

    The third (optional) argument to 'steenrod_algebra_basis' is the prime p:
        sage: steenrod_algebra_basis(9, 'milnor', p=3)
        (Q_1 P(1), Q_0 P(2))
        sage: steenrod_algebra_basis(9, 'milnor', 3)
        (Q_1 P(1), Q_0 P(2))
        sage: steenrod_algebra_basis(17, 'milnor', 3)
        (Q_2, Q_1 P(3), Q_0 P(0,1), Q_0 P(4))

    Other bases:
        sage: steenrod_algebra_basis(7,'admissible')
        (Sq^{7}, Sq^{6} Sq^{1}, Sq^{4} Sq^{2} Sq^{1}, Sq^{5} Sq^{2})
        sage: [x.basis('milnor') for x in steenrod_algebra_basis(7,'admissible')]
        [Sq(7),
        Sq(4,1) + Sq(7),
        Sq(0,0,1) + Sq(1,2) + Sq(4,1) + Sq(7),
        Sq(1,2) + Sq(7)]
        sage: steenrod_algebra_basis(13,'admissible',p=3)
        (beta P^{3}, P^{3} beta)
        sage: steenrod_algebra_basis(5,'wall')
        (Q^{2}_{2} Q^{0}_{0}, Q^{1}_{1} Q^{1}_{0})
        sage: steenrod_algebra_basis(5,'wall_long')
        (Sq^{4} Sq^{1}, Sq^{2} Sq^{1} Sq^{2})
        sage: steenrod_algebra_basis(5,'pst-rlex')
        (P^{0}_{1} P^{2}_{1}, P^{1}_{1} P^{0}_{2})
    """
    if not (isinstance(n, (Integer, int)) and n >= 0):
        raise ValueError, "%s is not a non-negative integer." % n
    basis_name = get_basis_name(basis, p)
    if basis_name.find('long') >= 0:
        long = True
        basis_name = basis_name.rsplit('_', 1)[0]
    else:
        long = False

    ## Milnor basis
    if basis_name == 'milnor':
        return milnor_basis(n,p)
    ## Serre-Cartan basis
    elif basis_name == 'serre-cartan':
        return serre_cartan_basis(n,p)
    ## Atomic bases
    elif p == 2 and (basis_name == 'woody' or basis_name == 'woodz'
                     or basis_name == 'wall' or basis_name == 'arnona'
                     or basis_name.find('pst') >= 0
                     or basis_name.find('comm') >= 0):
        return atomic_basis(n, basis_name, long=long)
    ## Arnon 'C' basis
    elif p == 2 and basis == 'arnonc':
        return arnonC_basis(n)


def restricted_partitions(n, list, no_repeats=False):
    """
    List of 'restricted' partitions of n: partitions with parts taken from list.

    INPUT:
        n -- non-negative integer
        list -- list of positive integers
        no_repeats -- boolean (optional, default = False), if True, only
            return partitions with no repeated parts

    This seems to be faster than RestrictedPartitions, although I
    don't know why.  Maybe because all I want is the list of
    partitions (with each partition represented as a list), not the
    extra stuff provided by RestrictedPartitions.

    EXAMPLES:
        sage: from sage.algebras.steenrod_algebra_bases import restricted_partitions
        sage: restricted_partitions(10, [7,5,1])
        [[7, 1, 1, 1], [5, 5], [5, 1, 1, 1, 1, 1], [1, 1, 1, 1, 1, 1, 1, 1, 1, 1]]
        sage: restricted_partitions(10, [6,5,4,3,2,1], no_repeats=True)
        [[6, 4], [6, 3, 1], [5, 4, 1], [5, 3, 2], [4, 3, 2, 1]]
        sage: restricted_partitions(10, [6,4,2])
        [[6, 4], [6, 2, 2], [4, 4, 2], [4, 2, 2, 2], [2, 2, 2, 2, 2]]
        sage: restricted_partitions(10, [6,4,2], no_repeats=True)
        [[6, 4]]

    'list' may have repeated elements.  If 'no_repeats' is False, this
    has no effect.  If 'no_repeats' is True, and if the repeated
    elements appear consecutively in 'list', then each element may be
    used only as many times as it appears in 'list':
        sage: restricted_partitions(10, [6,4,2,2], no_repeats=True)
        [[6, 4], [6, 2, 2]]
        sage: restricted_partitions(10, [6,4,2,2,2], no_repeats=True)
        [[6, 4], [6, 2, 2], [4, 2, 2, 2]]

    (If the repeated elements don't appear consecutively, the results
    are likely meaningless, containing several partitions more than
    once, for example.)

    In the following examples, 'no_repeats' is False:
        sage: restricted_partitions(10, [6,4,2])
        [[6, 4], [6, 2, 2], [4, 4, 2], [4, 2, 2, 2], [2, 2, 2, 2, 2]]
        sage: restricted_partitions(10, [6,4,2,2,2])
        [[6, 4], [6, 2, 2], [4, 4, 2], [4, 2, 2, 2], [2, 2, 2, 2, 2]]
        sage: restricted_partitions(10, [6,4,4,4,2,2,2,2,2,2])
        [[6, 4], [6, 2, 2], [4, 4, 2], [4, 2, 2, 2], [2, 2, 2, 2, 2]]
    """
    if n < 0:
        return []
    elif n == 0:
        return [[]]
    else:
        results = []
        if no_repeats:
            index = 1
        else:
            index = 0
        old_i = 0
        for i in list:
            if old_i != i:
                for sigma in restricted_partitions(n-i, list[index:], no_repeats):
                    results.append([i] + sigma)
            index += 1
            old_i = i
        return results


def list_to_hist(list):
    """
    Given list of positive integers [a,b,c,...], return corresponding 'histogram'.

    That is, in the output [n1, n2, ...], n1 is the number of 1's in
    the original list, n2 is the number of 2's, etc.

    INPUT:
        list -- list of positive integers
    OUTPUT:
        answer -- list of non-negative integers

    EXAMPLES:
        sage: from sage.algebras.steenrod_algebra_bases import list_to_hist
        sage: list_to_hist([1,2,3,4,2,1,2])
        [2, 3, 1, 1]
        sage: list_to_hist([2,2,2,2])
        [0, 4]
    """
    if len(list) == 0:
        return list
    else:
        answer = [0] * max(list)
        for i in list:
            answer[i-1] += 1
        return answer


########################################################
# Functions for defining bases. They work like this:

# First, check to see if the basis has been computed already: look for
#   the cached result in _steenrod_bases.  If so, filter the result
#   based on the optional argument bound.

# Dictionary in which to cache results:
_steenrod_bases = {}

# Second, if there is no cached result, build directly or recursively,
#   whichever is more convenient.  The details are documented for each
#   basis.  As each basis element is constructed, so is its dictionary
#   for the current basis.  For example, in 'serre_cartan_basis', each
#   new element is obtained by taking an element in a lower-dimensional
#   basis and multiplying, on the right, by Sq(last).  Call this element
#   'new'.  Then set new._raw['serre-cartan'] to be the dictionary of
#   Serre-Cartan monomials describing this element.  This makes
#   conversion of the element to the chosen basis instantaneous; this is
#   useful on its own, and crucial for computing change-of-basis
#   matrices for general conversion.

# Finally, if there is no cached result, cache the result.


def xi_degrees(n,p=2):
    r"""
    Decreasing list of degrees of the xi_i's, starting in degree n.

    INPUT:
        n -- integer

    OUTPUT:
        list -- list of integers

    When p=2: decreasing list of the degrees of the $\xi_i$'s with
    degree at most n.

    At odd primes: decreasing list of these degrees, each divided
    by 2(p-1).

    EXAMPLES:
        sage: sage.algebras.steenrod_algebra_bases.xi_degrees(17)
        [15, 7, 3, 1]
        sage: sage.algebras.steenrod_algebra_bases.xi_degrees(17,p=3)
        [13, 4, 1]
        sage: sage.algebras.steenrod_algebra_bases.xi_degrees(400,p=17)
        [307, 18, 1]
    """
    result = []
    deg = 1
    while deg <= n:
        result.insert(0,deg)
        deg = p*deg + 1
    return result


def milnor_basis(n, p=2):
    r"""
    Milnor basis in dimension $n$.

    INPUT:
        n -- non-negative integer
        p -- positive prime number (optional, default 2)

    OUTPUT:
        tuple of mod p Milnor basis elements in dimension n

    At the prime 2, the Milnor basis consists of symbols of the form
    $\text{Sq}(m_1, m_2, ..., m_t)$, where each $m_i$ is a
    non-negative integer and if $t>1$, then $m_t \neq 0$.
    At odd primes, it consists of symbols of the form
    $Q_{e_1} Q_{e_2} ... P(m_1, m_2, ..., m_t)$, where
    $0 \leq e_1 < e_2 < ...$, each $m_i$ is a non-negative integer,
    and if $t>1$, then $m_t \neq 0$.

    EXAMPLES:
        sage: from sage.algebras.steenrod_algebra_bases import milnor_basis
        sage: milnor_basis(7)
        (Sq(0,0,1), Sq(1,2), Sq(4,1), Sq(7))
        sage: milnor_basis(7, 2)
        (Sq(0,0,1), Sq(1,2), Sq(4,1), Sq(7))
        sage: milnor_basis(9, 3)
        (Q_1 P(1), Q_0 P(2))
        sage: milnor_basis(17, 3)
        (Q_2, Q_1 P(3), Q_0 P(0,1), Q_0 P(4))
        sage: milnor_basis(48, p=5)
        (P(0,1), P(6))
        sage: len(milnor_basis(100,3))
        13
        sage: len(milnor_basis(200,7))
        0
        sage: len(milnor_basis(240,7))
        3
    """
    if n == 0:
        if p == 2:
            return (Sq(0),)
        else:
            return (SteenrodAlgebra(p).P(0),)
    else:
        result = []
        if _steenrod_bases.has_key(('milnor',n,p)):
            result = _steenrod_bases[('milnor',n,p)]
        else:
            if p == 2:
                # build basis from partitions of n, with parts taken
                # from the list xi_degrees(n)
                for sigma in restricted_partitions(n,xi_degrees(n)):
                    sigma_exp = list_to_hist(sigma)
                    deg = 1
                    i = 0
                    exponents = []
                    while deg <= len(sigma_exp):
                        exponents.insert(i,sigma_exp[deg-1])
                        deg = 2*deg + 1
                        i = i + 1
                    result.append(SteenrodAlgebraElement({tuple(exponents): 1}))
                    _steenrod_bases[('milnor',n,p)] = tuple(result)
            else:  # p odd
                # first find the P part of each basis element.
                # in this part of the code (the P part), all dimensions are
                # divided by 2(p-1).
                for dim in range(n/(2*(p-1)) + 1):
                    if dim == 0:
                        P_result = [[0]]
                    else:
                        P_result = []
                    for sigma in restricted_partitions(dim, xi_degrees(dim,p)):
                        sigma_exp = list_to_hist(sigma)
                        deg = 1
                        i = 0
                        p_mono = []
                        while deg <= len(sigma_exp):
                            p_mono.insert(i, sigma_exp[deg-1])
                            deg = p*deg + 1
                            i = i + 1
                        if len(p_mono) > 0:
                            P_result.append(p_mono)
                    # now find the Q part of the basis element.
                    # dimensions here are back to normal.
                    for p_mono in P_result:
                        deg = n - 2*dim*(p-1)
                        q_degrees = [1+2*(p-1)*d for d in
                                     xi_degrees((deg - 1)/(2*(p-1)), p)] + [1]
                        q_degrees_decrease = q_degrees
                        q_degrees.reverse()
                        if deg % (2*(p-1)) <= len(q_degrees):
                            # if this inequality fails, no way to have a partition
                            # with distinct parts.
                            for sigma in restricted_partitions(deg,
                                                               q_degrees_decrease,
                                                               no_repeats = True):
                                index = 0
                                q_mono = []
                                for q in q_degrees:
                                    if q in sigma:
                                        q_mono.append(index)
                                    index += 1
                                result.append(SteenrodAlgebraElement(
                                        {(tuple(q_mono), tuple(p_mono)): 1}, p))
                _steenrod_bases[('milnor',n,p)] = tuple(result)
        return tuple(result)


def serre_cartan_basis(n, p=2, bound=1):
    r"""
    Serre-Cartan basis in dimension $n$.

    INPUT:
        n -- non-negative integer
        bound -- positive integer (optional)
        prime -- positive prime number (optional, default 2)

    OUTPUT:
        tuple of mod p Serre-Cartan basis elements in dimension n

    The Serre-Cartan basis consists of 'admissible monomials in the
    Steenrod squares'.  Thus at the prime 2, it consists of monomials
    $\text{Sq}^{m_1} \text{Sq}^{m_2} ... \text{Sq}^{m_t}$ with $m_i
    \geq 2m_{i+1}$ for each $i$.  At odd primes, it consists of
    monomials $\beta^{e_0} P^{s_1} \beta^{e_1} P^{s_2} ...  P^{s_k}
    \beta^{e_k}$ with each $e_i$ either 0 or 1, $s_i \geq p s_{i+1} +
    e_i$, and $s_k \geq 1$.

    EXAMPLES:
        sage: from sage.algebras.steenrod_algebra_bases import serre_cartan_basis
        sage: serre_cartan_basis(7)
        (Sq^{7}, Sq^{6} Sq^{1}, Sq^{4} Sq^{2} Sq^{1}, Sq^{5} Sq^{2})
        sage: serre_cartan_basis(13,3)
        (beta P^{3}, P^{3} beta)
        sage: serre_cartan_basis(50,5)
        (beta P^{5} P^{1} beta, beta P^{6} beta)

    If optional argument bound is present, include only those
    monomials whose last term is at least bound (when p=2), or those
    for which $s_k - epsilon_k >= bound$ (when p is odd).
        sage: serre_cartan_basis(7, bound=2)
        (Sq^{7}, Sq^{5} Sq^{2})
        sage: serre_cartan_basis(13, 3, bound=3)
        (beta P^{3},)
    """
    if not (isinstance(bound, (Integer, int)) and bound >= 1):
        raise ValueError, "%s is not a positive integer." % bound
    A = SteenrodAlgebra(p, 'serre-cartan')
    if n == 0:
        if p == 2:
            return (Sq(0),)
        else:
            return (A.P(0),)
    else:
        if _steenrod_bases.has_key(('serre-cartan',n,p)):
            result = []
            lookup = _steenrod_bases[('serre-cartan',n,p)]
            for vec in lookup:
                for mono in vec._basis_dictionary('serre-cartan'):
                    if (bound == 1) or \
                            (p == 2 and mono[-1] >= bound) or \
                            (p > 2 and (len(mono)<2 or mono[-2] - mono[-1] >= bound)):
                        result.append(vec)
        else:
            if p == 2:
                # Build basis recursively.  last = last term.
                # last is >= bound, and we will append (last,) to the end of
                # elements from serre_cartan_basis (n - last, bound=2 * last).
                # This means that 2 last <= n - last, or 3 last <= n.
                result = [A(Sq(n))]
                for last in range(bound, 1+n/3):
                    for vec in serre_cartan_basis(n - last, bound = 2*last):
                        new = vec * Sq(last)
                        for m in vec._basis_dictionary('serre-cartan'):
                            sc_mono = m
                        new._raw['serre-cartan'] = {sc_mono + (last,): 1}
                        result.append(A(new))
                if bound == 1: _steenrod_bases[('serre-cartan',n,p)] = tuple(result)
            else: # p odd
                if n % (2 * (p-1)) == 0 and n/(2 * (p-1)) >= bound:
                    a = A.P(int(n/(2 * (p-1))))
                    a._raw['serre-cartan'] = {(0,int(n/(2 * (p-1))),0): 1}
                    result = [a]
                elif n == 1:
                    a = A.Q(0)
                    a._raw['serre-cartan'] = {(1,): 1}
                    result = [a]
                else:
                    result = []
                # 2 cases: append P^{last}, or append P^{last} beta
                # case 1: append P^{last}
                for last in range(bound, 1+n/(2*(p - 1))):
                    if n - 2*(p-1)*last > 0:
                        for vec in serre_cartan_basis(n - 2*(p-1)*last,
                                                      p, p*last):
                            new = vec * A.P(last,)
                            for m in vec._basis_dictionary('serre-cartan'):
                                sc_mono = m
                            new._raw['serre-cartan'] = {sc_mono + (last,0): 1}
                            result.append(A(new))
                # case 2: append P^{last} beta
                if bound == 1:
                    bound = 0
                for last in range(bound+1, 1+n/(2*(p - 1))):
                    basis = serre_cartan_basis(n - 2*(p-1)*last - 1,
                                               p, p*last)
                    for vec in basis:
                        if vec == 1:
                            vec._raw['serre-cartan'] = {(0,): 1}
                        new = vec * A.P(last,) * A.Q(0)
                        for m in vec._basis_dictionary('serre-cartan'):
                            sc_mono = m
                        new._raw['serre-cartan'] = {sc_mono + (last,1): 1}
                        result.append(A(new))
                if bound <= 1: _steenrod_bases[('serre-cartan',n,p)] = tuple(result)
        return tuple(result)


def atomic_basis(n, basis, long=False):
    r"""
    Basis for dimension $n$ made of elements in 'atomic' degrees:
    degrees of the form $2^i (2^j - 1)$.

    INPUT:
        n -- non-negative integer
        basis -- string, the name of the basis

    OUTPUT:
        tuple of basis elements in dimension n

    The atomic bases include Wood's Y and Z bases, Wall's basis,
    Arnon's A basis, the $P^s_t$-bases, and the commutator bases.
    (All of these bases are constructed similarly, hence their
    constructions have been consolidated into a single function.
    Also, See the documentation for 'steenrod_algebra_basis' for
    descriptions of them.)

    EXAMPLES:
        sage: from sage.algebras.steenrod_algebra_bases import atomic_basis
        sage: atomic_basis(6,'woody')
        (Sq^{2} Sq^{3} Sq^{1}, Sq^{4} Sq^{2}, Sq^{6})
        sage: atomic_basis(8,'woodz')
        (Sq^{4} Sq^{3} Sq^{1}, Sq^{7} Sq^{1}, Sq^{6} Sq^{2}, Sq^{8})
        sage: atomic_basis(6,'woodz') == atomic_basis(6, 'woody')
        True
        sage: atomic_basis(9,'woodz') == atomic_basis(9, 'woody')
        False

    Wall's basis:
        sage: atomic_basis(6,'wall')
        (Q^{1}_{1} Q^{1}_{0} Q^{0}_{0}, Q^{2}_{2} Q^{1}_{1}, Q^{2}_{1})

    Elements of the Wall basis have an alternate, 'long'
    representation as monomials in the $\text{Sq}^{2^n}$s:
        sage: atomic_basis(6, 'wall', long=True)
        (Sq^{2} Sq^{1} Sq^{2} Sq^{1}, Sq^{4} Sq^{2}, Sq^{2} Sq^{4})

    Arnon's A basis:
        sage: atomic_basis(7,'arnona')
        (X^{0}_{0} X^{1}_{1} X^{2}_{2},
        X^{0}_{0} X^{2}_{1},
        X^{1}_{0} X^{2}_{2},
        X^{2}_{0})

    These also have a 'long' representation:
        sage: atomic_basis(7,'arnona',long=True)
        (Sq^{1} Sq^{2} Sq^{4},
        Sq^{1} Sq^{4} Sq^{2},
        Sq^{2} Sq^{1} Sq^{4},
        Sq^{4} Sq^{2} Sq^{1})

    $P^s_t$-bases:
        sage: atomic_basis(7,'pst_rlex')
        (P^{0}_{1} P^{1}_{1} P^{2}_{1},
        P^{0}_{1} P^{1}_{2},
        P^{2}_{1} P^{0}_{2},
        P^{0}_{3})
        sage: atomic_basis(7,'pst_llex')
        (P^{0}_{1} P^{1}_{1} P^{2}_{1},
        P^{0}_{1} P^{1}_{2},
        P^{0}_{2} P^{2}_{1},
        P^{0}_{3})
        sage: atomic_basis(7,'pst_deg')
        (P^{0}_{1} P^{1}_{1} P^{2}_{1},
        P^{0}_{1} P^{1}_{2},
        P^{0}_{2} P^{2}_{1},
        P^{0}_{3})
        sage: atomic_basis(7,'pst_revz')
        (P^{0}_{1} P^{1}_{1} P^{2}_{1},
        P^{0}_{1} P^{1}_{2},
        P^{0}_{2} P^{2}_{1},
        P^{0}_{3})

    Commutator bases:
        sage: atomic_basis(7,'comm_rlex')
        (c_{0,1} c_{1,1} c_{2,1}, c_{0,1} c_{1,2}, c_{2,1} c_{0,2}, c_{0,3})
        sage: atomic_basis(7,'comm_llex')
        (c_{0,1} c_{1,1} c_{2,1}, c_{0,1} c_{1,2}, c_{0,2} c_{2,1}, c_{0,3})
        sage: atomic_basis(7,'comm_deg')
        (c_{0,1} c_{1,1} c_{2,1}, c_{0,1} c_{1,2}, c_{0,2} c_{2,1}, c_{0,3})
        sage: atomic_basis(7,'comm_revz')
   	(c_{0,1} c_{1,1} c_{2,1}, c_{0,1} c_{1,2}, c_{0,2} c_{2,1}, c_{0,3})

    Long representations of commutator bases:
        sage: atomic_basis(7,'comm_revz', long=True)
        (s_{1} s_{2} s_{4}, s_{1} s_{24}, s_{12} s_{4}, s_{124})
    """
    def degree_dictionary(n, basis):
        """
        Dictionary of atomic degrees for basis up to degree n.

        The keys for the dictionary are the atomic degrees -- the
        numbers of the form 2^i (2^j - 1) -- which are less than or
        equal to n.  The value associated to such a degree depends on
        basis; it has the form ((s,t), x), where (s,t) is a pair of
        integers which indexes the corresponding element, and x is the
        element in the Milnor basis.
        """
        dict = {}
        if basis.find('wood') >= 0:
            k=0
            m=0
            deg = 2**m * (2**(k+1) - 1)
            while deg <= n:
                dict[deg] = ((m,k), Sq(deg))
                if m>0:
                    m = m - 1
                    k = k + 1
                else:
                    m = k + 1
                    k = 0
                deg = 2**m * (2**(k+1) - 1)
        elif basis.find('wall') >= 0 or basis.find('arnon') >= 0:
            k=0
            m=0
            deg = 2**k * (2**(m-k+1) - 1)
            while deg <= n:
                dict[deg] = ((m,k), Q(m,k,basis))
                if k == 0:
                    m = m + 1
                    k = m
                else:
                    k = k - 1
                deg = 2**k * (2**(m-k+1) - 1)
        elif basis.find('pst') >= 0 or basis.find('comm') >= 0:
            s=0
            t=1
            deg = 2**s * (2**t - 1)
            while deg <= n:
                if basis.find('pst') >= 0:
                    dict[deg] = ((s,t), pst(s,t))
                else:  # comm
                    dict[deg] = ((s,t), commutator(s,t))
                if s == 0:
                    s = t
                    t = 1
                else:
                    s = s - 1
                    t = t + 1
                deg = 2**s * (2**t - 1)
        return dict

    def sorting_pair(s,t,basis):   # pair used for sorting the basis
        if basis.find('wood') >= 0 and basis.find('z') >= 0:
            return (-s-t,-s)
        elif basis.find('wood') >= 0 or basis.find('wall') >= 0 or \
                basis.find('arnon') >= 0:
            return (-s,-t)
        elif basis.find('rlex') >= 0:
            return (t,s)
        elif basis.find('llex') >= 0:
            return (s,t)
        elif basis.find('deg') >= 0:
            return (s+t,t)
        elif basis.find('revz') >= 0:
            return (s+t,s)

    from sage.misc.misc import prod
    if long:
        basis_long_name = basis + "_long"
    else:
        basis_long_name = basis
    A = SteenrodAlgebra(2, basis_long_name)
    if n == 0:
        return (A.Sq(0),)
    else:
        result = []
        if _steenrod_bases.has_key((basis,n)):
            if basis == basis_long_name:
                result = _steenrod_bases[(basis,n)]
            else:
                result = tuple([A(a) for a in _steenrod_bases[(basis,n)]])
        else:
            degrees_etc = degree_dictionary(n, basis)
            degrees = degrees_etc.keys()
            for sigma in restricted_partitions(n, degrees, no_repeats=True):
                big_list = [degrees_etc[part] for part in sigma]
                big_list.sort(key=lambda x: x[0],
                              cmp = lambda x, y: cmp(sorting_pair(x[0], x[1], basis),
                                                     sorting_pair(y[0], y[1], basis)))
                # reverse = True)
                # arnon: sort like wall, then reverse end result
                if basis.find('arnon') >= 0:
                    big_list.reverse()
                list_of_pairs = [d[0] for d in big_list]
                mono_list = [d[1] for d in big_list]
                new = prod(mono_list)
                new._raw[basis] = {tuple(list_of_pairs): 1}
                result.append(A(new))
            _steenrod_bases[(basis,n)] = tuple(result)
        return tuple(result)


def Q(m,k,basis):
    r"""
    Compute $Q^m_k$ (Wall basis) and $X^m_k$ (Arnon's A basis).

    INPUT:
        m,k -- non-negative integers with $m >= k$
        basis -- 'wall' or 'arnona'

    OUTPUT:
        element of Steenrod algebra

    If basis is 'wall', this returns $Q^m_k = Sq(2^k) Sq(2^(k+1))
    ... Sq(2^m)$.  If basis is 'arnona', it returns the reverse of
    this: $X^m_k = Sq(2^m) ... Sq(2^(k+1)) Sq(2^k)$

    EXAMPLES:
        sage: from sage.algebras.steenrod_algebra_bases import Q
        sage: Q(2,2,'wall')
        Sq(4)
        sage: Q(2,2,'arnona')
        Sq(4)
        sage: Q(3,2,'wall')
        Sq(6,2) + Sq(12)
        sage: Q(3,2,'arnona')
        Sq(0,4) + Sq(3,3) + Sq(6,2) + Sq(12)
    """
    exponent = 2**k
    result = Sq(exponent)
    for i in range(m-k):
        exponent = exponent * 2
        if basis == 'wall':
            result = result * Sq(exponent)
        elif basis == 'arnona':
            result = Sq(exponent) * result
    return result


def arnonC_basis(n,bound=1):
    r"""
    Arnon's C basis in dimension $n$.

    INPUT:
        n -- non-negative integer
        bound -- positive integer (optional)

    OUTPUT:
        tuple of basis elements in dimension n

    The elements of Arnon's C basis are monomials of the form
    $\text{Sq}^{t_1} ... \text{Sq}^{t_m}$ where for each $i$, we
    have $t_i \leq 2t_{i+1}$ and $2^i | t_{m-i}$.

    EXAMPLES:
        sage: from sage.algebras.steenrod_algebra_bases import arnonC_basis
        sage: arnonC_basis(7)
        (Sq^{7}, Sq^{2} Sq^{5}, Sq^{4} Sq^{3}, Sq^{4} Sq^{2} Sq^{1})

    If optional argument bound is present, include only those
    monomials whose first term is at least as large as bound:
        sage: arnonC_basis(7,3)
        (Sq^{7}, Sq^{4} Sq^{3}, Sq^{4} Sq^{2} Sq^{1})
    """
    if not (isinstance(bound, (Integer, int)) and bound >= 1):
        raise ValueError, "%s is not a positive integer." % bound
    A = SteenrodAlgebra(2, 'arnonc')
    if n == 0:
        return (A.Sq(0),)
    else:
        if _steenrod_bases.has_key(('arnonc',n)):
            result = []
            lookup = _steenrod_bases[('arnonc',n)]
            for vec in lookup:
                for mono in vec._basis_dictionary('arnonc'):
                    if mono[0] >= bound:
                        result.append(vec)
        else:
            # Build basis recursively.  first = first term.
            # first is >= bound, and we will prepend (first,) to the
            # elements from arnonC_basis (n - first, first / 2).
            # first also must be divisible by 2**(len(old-basis-elt))
            # This means that 3 first <= 2 n.
            result = [A(Sq(n))]
            for first in range(bound,1+2*n/3):
                for vec in arnonC_basis(n - first, max(first/2,1)):
                    if len(vec._basis_dictionary('arnonc')) > 1: print "Error in Arnon's C basis!"
                    for m in vec._basis_dictionary('arnonc'):
                        arnonc_mono = m
                    if first % 2**len(arnonc_mono) == 0:
                        new = Sq(first) * vec
                        new._raw['arnonc'] = {(first,) + arnonc_mono: 1}
                        result.append(A(new))
            if bound == 1: _steenrod_bases[('arnonc',n)] = tuple(result)
        return tuple(result)


# cache commutators here:
_commutators = {}

def commutator(s,t):
    r"""
    $t$th iterated commutator of consecutive $\text{Sq}^{2^i}$'s.

    INPUT:
        s, t: integers

    OUTPUT:
        element of the Steenrod algebra

    If t=1, return $Sq(2^s)$.  Otherwise, return the commutator
    $[commutator(s,t-1), \text{Sq}(2^{s+t-1})]$.

    EXAMPLES:
        sage: from sage.algebras.steenrod_algebra_bases import commutator
        sage: commutator(1,2)
        Sq(0,2)
        sage: commutator(0,4)
        Sq(0,0,0,1)
        sage: commutator(2,2)
        Sq(0,4) + Sq(3,3)

    NOTES:
        commutator(0,n) is equal to $\text{Sq}(0,...,0,1)$, with the
        1 in the $n$th spot.  commutator(i,n) always has
        $\text{Sq}(0,...,0,2^i)$, with $2^i$ in the $n$th spot, as a
        summand, but there may be other terms, as the example of
        commutator(2,2) illustrates.

        That is, commutator(s,t) is equal to $P^s_t$, possibly plus
        other Milnor basis elements.
    """
    if _commutators.has_key((s,t)):
        return _commutators[(s,t)]
    if t == 1:
        answer = Sq(2**s)
    else:
        x = commutator(s,t-1)
        y = Sq(2**(s+t-1))
        answer = x*y + y*x
    _commutators[(s,t)] = answer
    return answer


#############################################################################
def steenrod_basis_error_check(dim,p):
    """
    This performs crude error checking.

    INPUT:
        dim -- non-negative integer
        p -- positive prime number

    OUTPUT:
        None

    This checks to see if the different bases have the same length,
    and if the change-of-basis matrices are invertible.  If something
    goes wrong, an error message is printed.

    This function checks at the prime p as the dimension goes up from
    0, in which case the basis functions use the saved basis
    computations in lower dimensions in the computations.  It also
    checks as the dimension goes down from the top, in which case it
    doesn't have access to the saved computations.  (The saved
    computations are deleted first: the cache _steenrod_bases is set
    to {} before doing the computations.)

    EXAMPLES:
        sage: sage.algebras.steenrod_algebra_bases.steenrod_basis_error_check(12,2)
        p=2, in decreasing order of dimension, starting in dimension 12.
        down to dimension  10
        down to dimension  5
        p=2, now in increasing order of dimension, up to dimension 12
        up to dimension  0
        up to dimension  5
        up to dimension  10
        done checking
        sage: sage.algebras.steenrod_algebra_bases.steenrod_basis_error_check(30,3)
        p=3, in decreasing order of dimension, starting in dimension 30.
        down to dimension  30
        down to dimension  25
        down to dimension  20
        down to dimension  15
        down to dimension  10
        down to dimension  5
        p=3, now in increasing order of dimension, up to dimension 30
        up to dimension  0
        up to dimension  5
        up to dimension  10
        up to dimension  15
        up to dimension  20
        up to dimension  25
        done checking
    """
    global _steenrod_bases

    _steenrod_bases = {}

    if p == 2:
        bases = ('adem','woody', 'woodz', 'wall', 'arnona', 'arnonc',
                 'pst_rlex', 'pst_llex', 'pst_deg', 'pst_revz',
                 'comm_rlex', 'comm_llex', 'comm_deg', 'comm_revz')
    else:
        bases = ('adem', 'milnor')

    print "p=%s, in decreasing order of dimension, starting in dimension %s." % (p, dim)
    for i in range(dim,0,-1):
        if i % 5 == 0: print "down to dimension ", i
        for B in bases:
            if len(steenrod_algebra_basis(i,'milnor')) != len(steenrod_algebra_basis(i,B)):
                print "problem with milnor/" + B + " in dimension ", i
            mat = convert_to_milnor_matrix(i,B,p)
            if mat.nrows() != 0 and not mat.is_invertible():
                print "%s invertibility problem in dim %s at the prime %s" % (B, i, p)

    _steenrod_bases = {}

    print "p=%s, now in increasing order of dimension, up to dimension %s" % (p, dim)
    for i in range(dim):
        if i % 5 == 0: print "up to dimension ", i
        for B in bases:
            if len(steenrod_algebra_basis(i,'milnor')) != len(steenrod_algebra_basis(i,B)):
                print "problem with milnor/" + B + " in dimension ", i
            mat = convert_to_milnor_matrix(i,B,p)
            if mat.nrows() != 0 and not mat.is_invertible():
                print "%s invertibility problem in dim %s at the prime %s" % (B, i, p)

    print "done checking"
