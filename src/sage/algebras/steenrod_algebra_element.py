r"""
Steenrod algebra elements

AUTHORS:
    - John H. Palmieri (2008-07-30: version 0.9)

This package provides for basic algebra with elements in the mod $p$
Steenrod algebra.  In this package, elements in the Steenrod algebra
are represented, by default, using the Milnor basis.

EXAMPLES:
Basic arithmetic, $p=2$.  To construct an element of the mod 2
Steenrod algebra, use the function \code{Sq}:
    sage: a = Sq(1,2)
    sage: b = Sq(4,1)
    sage: z = a + b
    sage: z
    Sq(1,2) + Sq(4,1)
    sage: Sq(4) * Sq(1,2)
    Sq(1,1,1) + Sq(2,3) + Sq(5,2)
    sage: z**2          # non-negative exponents work as they should
    Sq(1,2,1) + Sq(4,1,1)
    sage: z**0
    Sq(0)

Basic arithmetic, $p>2$.  To construct an element of the mod $p$
Steenrod algebra when $p$ is odd, you should first define a Steenrod
algebra, using the \code{SteenrodAlgebra} command:
    sage: SteenrodAlgebra()   # 2 is the default prime
    mod 2 Steenrod algebra
    sage: A3 = SteenrodAlgebra(3)
    sage: A3
    mod 3 Steenrod algebra

Having done this, the newly created algebra \code{A3} has methods
\code{Q} and \code{P} which construct elements of \code{A3}:
    sage: c = A3.Q(1,3,6); c
    Q_1 Q_3 Q_6
    sage: d = A3.P(2,0,1); d
    P(2,0,1)
    sage: c * d
    Q_1 Q_3 Q_6 P(2,0,1)
    sage: e = A3.P(3)
    sage: d * e
    P(5,0,1)
    sage: e * d
    P(1,1,1) + P(5,0,1)
    sage: c * c
    0
    sage: e ** 3
    2 P(1,2)

Note that one can construct an element like \code{c} above in one step,
without first constructing the algebra:
    sage: c = SteenrodAlgebra(3).Q(1,3,6)
    sage: c
    Q_1 Q_3 Q_6

And of course, you can do similar constructions with the mod 2
Steenrod algebra:
    sage: A = SteenrodAlgebra(2); A
    mod 2 Steenrod algebra
    sage: A.Sq(2,3,5)
    Sq(2,3,5)
    sage: A.P(2,3,5)   # when p=2, P = Sq
    Sq(2,3,5)
    sage: A.Q(1,4)     # when p=2, this gives a product of Milnor primitives
    Sq(0,1,0,0,1)

Regardless of the prime, each element has an \code{excess}, and if the
element is homogeneous, a \code{degree}.  The excess of
$\text{Sq}(i_1,i_2,i_3,...)$ is $i_1 + i_2 + i_3 + ...$; when $p$
is odd, the excess of $Q_{0}^{e_0} Q_{1}^{e_1}
... \mathcal{P}(r_1, r_2, ...)$ is $\sum e_i + 2 \sum r_i$.  The
excess of a linear combination of Milnor basis elements is the minimum
of the excesses of those basis elements.

The degree of $\text{Sq}(i_1,i_2,i_3,...)$ is $sum (2^n-1) i_n$, and
when $p$ is odd, the degree of $Q_{0}^{\epsilon_0} Q_{1}^{\epsilon_1}
... \mathcal{P}(r_1, r_2, ...)$ is $\sum \epsilon_i (2p^i - 1) + \sum
r_j (2p^j - 2)$.  The degree of a linear combination of such terms is
only defined if the terms all have the same degree.

Here are some simple examples:
    sage: z = Sq(1,2) + Sq(4,1)
    sage: z.degree()
    7
    sage: (Sq(0,0,1) + Sq(5,3)).degree()
    Element is not homogeneous.
    sage: Sq(7,2,1).excess()
    10
    sage: z.excess()
    3
    sage: B = SteenrodAlgebra(3)
    sage: x = B.Q(1,4)
    sage: y = B.P(1,2,3)
    sage: x.degree()
    166
    sage: x.excess()
    2
    sage: y.excess()
    12

Elements have a \code{weight} in the May filtration, which (when $p=2$) is
related to the \code{height} function defined by Wall:
    sage: Sq(2,1,5).may_weight()
    9
    sage: Sq(2,1,5).wall_height()
    [2, 3, 2, 1, 1]
    sage: b = Sq(4)*Sq(8) + Sq(8)*Sq(4)
    sage: b.may_weight()
    2
    sage: b.wall_height()
    [0, 0, 1, 1]

Odd primary May weights:
    sage: A5 = SteenrodAlgebra(5)
    sage: a = A5.Q(1,2,4)
    sage: b = A5.P(1,2,1)
    sage: a.may_weight()
    10
    sage: b.may_weight()
    8
    sage: (a * b).may_weight()
    18
    sage: A5.P(0,0,1).may_weight()
    3

Since the Steenrod algebra is a Hopf algebra, every element has an
antipode.
    sage: d = Sq(0,0,1); d
    Sq(0,0,1)
    sage: d.antipode()
    Sq(0,0,1)
    sage: Sq(4).antipode()
    Sq(1,1) + Sq(4)
    sage: (Sq(4) * Sq(2)).antipode()
    Sq(6)
    sage: SteenrodAlgebra(7).P(3,1).antipode()
    4 P(3,1) + 5 P(11)

Applying the antipode twice returns the original element:
    sage: y = Sq(8)*Sq(4)
    sage: y == (y.antipode()).antipode()
    True

You can treat elements of the Steenrod algebra like lists of Milnor
basis elements:
    sage: y = Sq(4) * Sq(1,2); y
    Sq(1,1,1) + Sq(2,3) + Sq(5,2)
    sage: for m in y: m
    Sq(1,1,1)
    Sq(2,3)
    Sq(5,2)
    sage: [(m.degree(),m.excess()) for m in y]
    [(11, 3), (11, 5), (11, 7)]

Once you've define a Steenrod algebra, the method \code{pst} is another
way to define elements of it: \code{pst(s,t)} defines the Margolis element
$P^{s}_{t}$, the basis element $\mathcal{P}(0,...,0,p^s)$ with $p^s$
in position $t$:
    sage: A2 = SteenrodAlgebra(2)
    sage: Q2 = A2.pst(0,3)
    sage: Q2
    Sq(0,0,1)
    sage: Q2*Q2
    0
    sage: A2.pst(1,2) == Sq(2)*Sq(4) + Sq(4)*Sq(2)
    True
    sage: A5 = SteenrodAlgebra(5)
    sage: A5.pst(2,2)
    P(0,25)

There are a number of different bases available in which to represent
elements of the Steenrod algebra.  When $p>2$, the choices are the
Milnor basis ('milnor') or the Serre-Cartan basis ('serre-cartan' or
'adem' or 'admissible').  When $p=2$, the choices are those, along
with Wood's Y basis ('wood_y'), Wood's Z basis ('wood_z'), Wall's
basis ('wall' or 'wall_long'), Arnon's A basis ('arnon_a' or
'arnon_a_long'), Arnon's C basis ('arnon_c'), various $P^s_t$ bases
('pst_ORDER' for various values of ORDER), and various commutator
bases ('comm_ORDER' or 'comm_ORDER_long' for various values of ORDER).

See documentation for the function \code{steenrod_algebra_basis} for
full descriptions of these bases.

To access representations of elements with respect to these different
bases, you can either use the \code{basis} method for an element, or define
a Steenrod algebra with respect to a particular basis and then use that:
    sage: c = Sq(2) * Sq(1); c
    Sq(0,1) + Sq(3)
    sage: c.basis('serre-cartan')
    Sq^{2} Sq^{1}
    sage: c.basis('milnor')
    Sq(0,1) + Sq(3)
    sage: adem = SteenrodAlgebra(2, 'serre-cartan')
    sage: x = Sq(7,3,1)   # top class in the subalagebra A(2)
    sage: adem(x)
    Sq^{17} Sq^{5} Sq^{1}
    sage: SteenrodAlgebra(2, 'pst')(x)
    P^{0}_{1} P^{0}_{2} P^{1}_{1} P^{0}_{3} P^{1}_{2} P^{2}_{1}

Multiplication works within bases:
    sage: adem = SteenrodAlgebra(2, 'adem')
    sage: x = adem.Sq(5)
    sage: y = adem.Sq(1)
    sage: x * y
    Sq^{5} Sq^{1}

Multiplying elements defined with respect to two different bases may
have unpredictable results (as far as the basis in which the result is
printed):
    sage: milnor = SteenrodAlgebra(2, 'milnor')
    sage: xm = milnor.Sq(5)
    sage: ym = milnor.Sq(1)
    sage: xm * ym
    Sq(3,1)
    sage: xm * y
    Sq^{5} Sq^{1}
    sage: x * ym
    Sq^{5} Sq^{1}

Several of these bases ('arnon_a', 'wall', 'comm') have alternate,
longer, representations.  These provide ways of expressing elements of
the Steenrod algebra in terms of the $\text{Sq}^{2^n}$.
    sage: Sq(6).basis('arnon_a_long')
    Sq^{1} Sq^{2} Sq^{1} Sq^{2} + Sq^{2} Sq^{4}
    sage: Sq(6).basis('wall_long')
    Sq^{2} Sq^{1} Sq^{2} Sq^{1} + Sq^{2} Sq^{4}
    sage: SteenrodAlgebra(2,'comm_deg_long')(Sq(6))
    s_{1} s_{2} s_{12} + s_{2} s_{4}

**************

INTERNAL DOCUMENTATION:

Here are details on the class \code{SteenrodAlgebraElement} (for people who
want to delve into or extend the code):

Attributes for a \code{SteenrodAlgebraElement} self:

  \code{self._base_field}: $GF(p)$, where $p$ is the associated prime

  \code{self._prime}: $p$

  \code{self._basis}: basis in which to print this element

  \code{self._raw}: dictionary.  keys are basis names, taken from
    \code{_steenrod_basis_unique_names}, and the associated values are
    dictionaries themselves; if the dictionary is nonempty, it gives
    the representation for that element in the given basis.  If it is
    empty, that means that the representation in that basis hasn't
    been computed yet.  The representation of an element with respect
    to a basis (other than the Milnor basis, which is how elements are
    stored internally) isn't computed until requested, either by
    calling the method \code{_basis_dictionary('basis_name')}, or more
    typically, by calling the method \code{basis('basis_name')} or by
    defining a Steenrod algebra at that basis and applying its call
    method to the element.

    The dictionaries are defined as follows.  In the Milnor basis at
    the prime 2, for example, since monomials are of the form
    $\text{Sq}(a,b,c,...)$, then monomials are stored as tuples of
    integers \code{(a,b,c,...)}.  Thus if $y = \text{Sq}(5,3) +
    \text{Sq}(0,0,2)$, then \code{y._raw['milnor']} is \code{\{(0, 0,
    2): 1, (5, 3): 1\}}.  (The 1's following the colons are the
    coefficients of the monomials associated to the tuples.)  Each
    basis has its own representation as a dictionary; Arnon's C basis
    represents basis elements as tuples of integers, just like the
    Milnor basis and the Serre-Cartan basis, while the other bases
    represent basis elements as tuples of pairs of integers.  From the
    descriptions of the bases given in the file
    'steenrod_algebra_bases.py', it should be clear how to associate a
    tuple of pairs of integers to a basis element.  See also the
    function \code{string_rep}.

    When the element is initially defined by calling \code{Sq} or
    \code{SteenrodAlgebraElement}, typically only the 'milnor' dictionary
    is non-empty, while if the element is defined by the function
    \code{steenrod_algebra_basis}, its dictionary for the given basis is
    also initialized correctly.  For example:
        sage: B = steenrod_algebra_basis(6,'adem'); B
        (Sq^{6}, Sq^{5} Sq^{1}, Sq^{4} Sq^{2})
        sage: x = B[1]; x
        Sq^{5} Sq^{1}
        sage: x._raw
        {'milnor': {(3, 1): 1}, 'serre-cartan': {(5, 1): 1}}

    Note that the keys 'milnor' and 'serre-cartan' (a synonym for
    'adem') have nonempty associated values.

    When any element is converted to another basis (by changing the
    basis and then printing the element), its dictionary for that
    basis gets stored, as well:
        sage: x.basis('arnona')
        X^{0}_{0} X^{1}_{0} X^{1}_{1}
        sage: x._raw
        {'arnona': {((0, 0), (1, 0), (1, 1)): 1},
        'milnor': {(3, 1): 1},
        'serre-cartan': {(5, 1): 1}}

Methods for a \code{SteenrodAlgebraElement} self:

Most of these are self-explanatory.

    \code{_mul_}: multiply two elements.  This is done using Milnor
    multiplication, the code for which is in a separate file,
    'steenrod_milnor_multiplication'.  In a long computation, it seems
    that a lot of time is spent here, so one way to speed things up
    would be to optimize the Milnor multiplication routine.

    \code{_basis_dictionary}: compute the dictionary of the element
    with respect to the given basis.  This is basically done by doing
    a basis conversion from the Milnor basis to the given basis.
    There are two parts to this function; first, some elements (e.g.,
    $\text{Sq}(2^n)$) may be easy to convert directly.  This is done
    one basis at a time, and so takes up most of the lines of code.
    If the element is not recognizable as being easy to convert, then
    the function \code{milnor_convert} from the file
    'steenrod_algebra_bases.py' is called.  This does linear algebra:
    it computes the Milnor basis and the new basis in the appropriate
    dimension, computes the change-of-basis matrix, etc.

    \code{basis}: display the element in the given basis.

    \code{milnor}: display the element in the Milnor basis.

    \code{serre_cartan}: display the element in the Serre-Cartan basis.

    \code{adem}: display the element in the Serre-Cartan basis.

    \code{_repr_} and \code{_latex_} call the function
    \code{string_rep}, which has cases depending on the basis.

REFERENCES:

    [Mil] J. W. Milnor, "The Steenrod algebra and its dual," Ann. of Math.
          (2) 67 (1958), 150--171.

    [Mon] K. G. Monks, "Change of basis, monomial relations, and $P^s_t$
          bases for the Steenrod algebra," J. Pure Appl. Algebra 125 (1998),
          no. 1-3, 235--260.

    [Woo] R. M. W. Wood, "Problems in the Steenrod algebra," Bull. London
          Math. Soc. 30 (1998), no. 5, 449--517.
"""

#*****************************************************************************
#       Copyright (C) 2008 John H. Palmieri <palmieri@math.washington.edu>
#  Distributed under the terms of the GNU General Public License (GPL)
#*****************************************************************************

from sage.rings.ring import Algebra
from sage.algebras.algebra_element import AlgebraElement
from sage.structure.parent_gens import ParentWithGens
from sage.structure.element import RingElement
from sage.rings.all import GF
from sage.misc.functional import parent
from sage.rings.integer import Integer

def check_and_trim(nums):
    """
    Check that list or tuple consists of non-negative integers, and
    strip trailing zeroes.

    INPUT:
        nums -- a list or tuple

    OUTPUT:
        new -- a list or tuple

    If nums contains anything other than a non-negative integer, raise
    an exception, identifying the right-most problematic entry.
    Otherwise, return a new list or tuple, obtained from nums by
    omitting any zeroes from the end.

    EXAMPLES:
        sage: from sage.algebras.steenrod_algebra_element import check_and_trim
        sage: check_and_trim([3,4,1])
        [3, 4, 1]
        sage: a=[3,2,1,0,0]
        sage: check_and_trim(a)
        [3, 2, 1]
        sage: a    # check_and_trim doesn't affect its input
        [3, 2, 1, 0, 0]
        sage: check_and_trim([0]*127)
        []
        sage: check_and_trim((1,2,3,4,0,0,0))  # works on tuples, too
        (1, 2, 3, 4)
    """
    index = len(nums)
    for i in range(index-1, -1, -1):
        if nums[i] == 0 and i == index - 1:
            index = i
        if not (isinstance(nums[i], (Integer, int, long)) and nums[i] >= 0):
            print type(nums[i])
            raise ValueError, "%s is not a non-negative integer" % nums[i]
    return nums[:index]


def convert_perm(m):
    """
    Convert tuple m of non-negative integers to a permutation in one-line form.

    INPUT:
        m -- tuple of non-negative integers with no repetitions
    OUTPUT:
        list -- conversion of m to a permutation of the set {1,2,...,len(m)}

    If m=(3,7,4), then one can view m as representing the permutation
    of the set {3,4,7} sending 3 to 3, 4 to 7, and 7 to 4.  This
    function converts m to the list [1,3,2], which represents
    essentially the same permutation, but of the set {1,2,3}.  This
    list can then be passed to Permutation, and its signature can be
    computed.

    EXAMPLES:
        sage: sage.algebras.steenrod_algebra_element.convert_perm((3,7,4))
        [1, 3, 2]
        sage: sage.algebras.steenrod_algebra_element.convert_perm((5,0,6,3))
        [2, 4, 1, 3]
    """
    m2 = list(m)
    m2.sort()
    return [list(m).index(x)+1 for x in m2]


def base_p_expansion(n, p):
    r"""
    Return list of digits in the base p expansion of n.

    INPUT:
        n -- non-negative integer
        p -- positive prime number

    OUTPUT:
        list of digits in the base p expansion of n

    EXAMPLES:
        sage: sage.algebras.steenrod_algebra_element.base_p_expansion(10,2)
        [0, 1, 0, 1]
        sage: sage.algebras.steenrod_algebra_element.base_p_expansion(10,3)
        [1, 0, 1]
        sage: sage.algebras.steenrod_algebra_element.base_p_expansion(10,5)
        [0, 2]
        sage: sage.algebras.steenrod_algebra_element.base_p_expansion(10,7)
        [3, 1]
        sage: sage.algebras.steenrod_algebra_element.base_p_expansion(0,7)
        []
        sage: sage.algebras.steenrod_algebra_element.base_p_expansion(100000,13)
        [4, 9, 6, 6, 3]
    """
    result = []
    while n > 0:
        remainder = n % p
        result.append(remainder)
        n = int((n - remainder)/p)
    return result


def integer_base_2_log(n):
    """
    Largest integer k so that $2^k <= n$

    INPUT:
        n -- positive integer

    OUTPUT:
        k -- integer

    This returns the integer $k$ so that $2^k <= n$ and $2^{k+1} > n$.

    EXAMPLES:
        sage: sage.algebras.steenrod_algebra_element.integer_base_2_log(7)
        2
        sage: sage.algebras.steenrod_algebra_element.integer_base_2_log(8)
        3
        sage: sage.algebras.steenrod_algebra_element.integer_base_2_log(9)
        3
    """
    answer = 0
    while 2**(answer+1) <= n:
        answer += 1
    return answer

# These are for internal use only; they are the names used in the code
# to represent each basis.
_steenrod_basis_unique_names_odd =  ('serre-cartan', 'milnor')
_steenrod_basis_unique_names = _steenrod_basis_unique_names_odd + \
    ('pst_rlex', 'pst_llex', 'pst_deg', 'pst_revz',
     'comm_rlex', 'comm_llex', 'comm_deg', 'comm_revz',
     'comm_rlex_long', 'comm_llex_long', 'comm_deg_long', 'comm_revz_long',
     'woody', 'woodz', 'arnona', 'arnona_long', 'arnonc', 'wall', 'wall_long')

# This dictionary is for caching Steenrod algebras at various primes,
# to make sure that SteenrodAlgebraElements defined at the same prime
# have the same parent.
_steenrod_algebras = {}

class SteenrodAlgebraElement(AlgebraElement):
    r"""
    Element of the mod p Steenrod algebra.

    At the prime 2, use the function 'Sq' to define these, as in
    'w=Sq(4,3,3)' or 'z=Sq(1,2)+Sq(4,1)' or 'q=Sq(8)*Sq(4) + Sq(12)'.

    At odd primes, use the methods 'P' and 'Q' to define these, as
    in 'w=SteenrodAlgebra(3).Q(1,5) * SteenrodAlgebra(3).P(4,3)'.

    EXAMPLES:
        sage: w = Sq(4,3,3)
        sage: w
        Sq(4,3,3)

    The function 'Sq', together with addition, provides an easy way to
    define elements when $p=2$:
        sage: b = Sq(3) + Sq(0,1)
        sage: b
        Sq(0,1) + Sq(3)

    When $p$ is odd, first define a Steenrod algebra to specify the
    prime, and then use the methods 'P' and 'Q', together with
    multiplication and addition:
        sage: A7 = SteenrodAlgebra(7)
        sage: u = A7.Q(0,4); u
        Q_0 Q_4
        sage: v = A7.P(1,2,3); v
        P(1,2,3)
        sage: u * v
        Q_0 Q_4 P(1,2,3)
        sage: 10 * u * v
        3 Q_0 Q_4 P(1,2,3)
        sage: u + v
        P(1,2,3) + Q_0 Q_4
    """

    def __init__(self, poly, p=2, basis='milnor'):
        r"""
        INPUT:
            poly -- dictionary with entries of form (monomial: coefficient)
            Each coefficient is in GF(p), and each monomial is a tuple
            of non-negative integers (a, b, c, ...), corresponding to
            the Milnor basis element Sq(a, b, c, ...).

            At odd primes, the monomials are pairs of tuples: they are
            of the form ((e0, e1, e2, ...), (r1, r2, ...)),
            corresponding to the element $Q_{e_0} Q_{e_1} ... P(r_1, r_2, ...)$.

            Alternatively, poly can be an integer n, in which case it
            is viewed as being in the field GF(p): the resulting
            element is n * Sq(0).

            p -- positive prime number (default 2)

        OUTPUT:
            element of the mod p Steenrod algebra

        EXAMPLES:
            sage: from sage.algebras.steenrod_algebra_element import SteenrodAlgebraElement
            sage: SteenrodAlgebraElement({(1,2,3): 1}, 2)
            Sq(1,2,3)
            sage: SteenrodAlgebraElement({(1,2,3): 4}, 2)
            0
            sage: SteenrodAlgebraElement({((0,3), (1,2)): 5}, 7)
            5 Q_0 Q_3 P(1,2)
            sage: SteenrodAlgebraElement({((0,3), (1,2)): 5}, p=7)
            5 Q_0 Q_3 P(1,2)

        The input can also be an integer, in which case it is treated
        as a multiple, mod p, of the unit element P(0) (a.k.a. Sq(0)
        at the prime 2):
            sage: SteenrodAlgebraElement(3,2)
            Sq(0)
            sage: SteenrodAlgebraElement(6)   # p=2 is the default prime
            0
            sage: SteenrodAlgebraElement(6, p=5)
            P(0)
        """
        from sage.rings.arith import is_prime
        from sage.algebras.steenrod_algebra import SteenrodAlgebra
        if not is_prime(p):
            raise ValueError, "%s is not prime." % p
        else:
            # get cached Steenrod algebra at the prime p
            if _steenrod_algebras.has_key((p,basis)):
                alg = _steenrod_algebras[(p,basis)]
            else:
                _steenrod_algebras[(p,basis)] = SteenrodAlgebra(p,basis)
                alg = _steenrod_algebras[(p,basis)]
            if isinstance(poly, SteenrodAlgebraElement):
                if poly.parent().prime == p:
                    self._raw = poly._raw
                else:
                    raise ValueError, "Mismatch: %s is defined at the prime %s, \
not %s" % (poly, poly.parent().prime, p)
            # basic initializations:
            RingElement.__init__(self, alg)
            F = alg.base_ring()
            self._base_field = F
            self._prime = p
            self._basis = basis
            # now comes most of the work:
            if isinstance(poly, dict):
                new_poly = {}
                for mono in poly:
                    if p == 2:
                        # when p=2, mono is a tuple of integers
                        trimmed = check_and_trim(mono)
                        if new_poly.has_key(trimmed):
                            coeff = F(poly[mono] + new_poly[trimmed])
                        else:
                            coeff = F(poly[mono])
                        if not coeff.is_zero():
                            new_poly[trimmed] = coeff
                    else:
                        # when p is odd, mono is either an empty tuple
                        # or a pair of tuples
                        if len(mono) == 0:
                            if new_poly.has_key(mono):
                                coeff = F(poly[mono] + new_poly[mono])
                            else:
                                coeff = F(poly[mono])
                            if not coeff.is_zero():
                                new_poly[mono] = coeff
                        else:
                            # mono is a pair of tuples
                            mono1, mono2 = mono
                            multiplier = 1
                            # see if there are any repetitions in the Q monomial:
                            if len(mono1) != len(set(mono1)):
                                return self(0)
                            # if not, sort them and introduce the correct sign:
                            if len(mono1) > 0:
                                from sage.combinat.permutation import Permutation
                                multiplier = Permutation(
                                    convert_perm(mono1)).signature()
                                mono1 = tuple(sorted(mono1))
                            trimmed2 = check_and_trim(mono2)
                            if len(mono1) + len(trimmed2) > 0:
                                if new_poly.has_key((mono1,trimmed2)):
                                    coeff = F(multiplier * poly[mono] \
                                                  + new_poly[(mono1, trimmed2)])
                                else:
                                    coeff = F(multiplier * poly[mono])
                                if not coeff.is_zero():
                                    new_poly[(mono1,trimmed2)] = coeff
                            else:
                                if new_poly.has_key(()):
                                    coeff = F(poly[mono] + new_poly[()])
                                else:
                                    coeff = F(poly[mono])
                                if not coeff.is_zero():
                                    new_poly[()] = coeff
                # now define the _raw attribute: a dictionary keyed by the
                # basis names.  set the 'milnor' value to by new_poly
                self._raw = {}
                self._raw['milnor'] = new_poly
            elif isinstance(poly, (Integer, int)) or poly.parent() == F:
                # a scalar, i.e., a scalar multiple of Sq(0) or P(0)
                if F(poly) != 0:
                    new_poly = {(): F(poly)}
                else:
                    new_poly = {}
                self._raw = {}
                if p == 2:
                    basis_list = _steenrod_basis_unique_names
                else:
                    basis_list = _steenrod_basis_unique_names_odd
                for basis in basis_list:
                    self._raw[basis] = new_poly
            elif isinstance(poly, SteenrodAlgebraElement):
                self._raw = poly._raw
            else:
                raise ValueError, "%s is not of the correct form" % poly


    def is_unit(self):
        """
        True if element has a nonzero scalar multiple of P(0) as a
        summand, False otherwise.

        EXAMPLES:
            sage: z = Sq(4,2) + Sq(7,1) + Sq(3,0,1)
            sage: z.is_unit()
            False
            sage: u = 1 + Sq(3,1)
            sage: u == Sq(0) + Sq(3,1)
            True
            sage: u.is_unit()
            True
            sage: A5 = SteenrodAlgebra(5)
            sage: v = A5.P(0)
            sage: (v + v + v).is_unit()
            True
        """
        if self._raw['milnor'].has_key(()):
            return True
        else:
            return False

    def is_nilpotent(self):
        """
        True if element is not a unit, False otherwise.

        EXAMPLES:
            sage: z = Sq(4,2) + Sq(7,1) + Sq(3,0,1)
            sage: z.is_nilpotent()
            True
            sage: u = 1 + Sq(3,1)
            sage: u == Sq(0) + Sq(3,1)
            True
            sage: u.is_nilpotent()
            False
        """
        return not self.is_unit()


    def _add_(self, other):
        """
        Addition for elements of the Steenrod algebra.

        EXAMPLES:
            sage: a = Sq(3,1) + Sq(6)
            sage: b = Sq(0,2)
            sage: c = Sq(6) + Sq(0,2)
            sage: a + b
            Sq(0,2) + Sq(3,1) + Sq(6)
            sage: b + c
            Sq(6)
            sage: A7 = SteenrodAlgebra(7)
            sage: x = A7.P(1,2); y = A7.Q(3,4) * A7.P(1,2)
            sage: x + 3 * y
            P(1,2) + 3 Q_3 Q_4 P(1,2)
            sage: x + 10 * y
            P(1,2) + 3 Q_3 Q_4 P(1,2)
        """
        F = self._base_field
        poly1 = self._raw['milnor']
        poly2 = other._raw['milnor']
        result = poly1.copy()
        for mono in poly2:
            if result.has_key(mono):
                coeff = F(result[mono] + poly2[mono])
            else:
                coeff = F(poly2[mono])
            if coeff == 0:
                del result[mono]
            else:
                result[mono] = coeff
        sum = SteenrodAlgebraElement(result, p=self._prime)
        return sum


    def _neg_(self):
        """
        Multiply every coefficient by the element -1 of the base field.

        EXAMPLES:
            sage: a = Sq(4,2) + Sq(5)
            sage: -a
            Sq(4,2) + Sq(5)
            sage: A5 = SteenrodAlgebra(5)
            sage: b = 2 * A5.P(2,0,1)
            sage: b
            2 P(2,0,1)
            sage: -b
            3 P(2,0,1)
        """
        p = self._prime
        if p == 2:
            return self
        else:
            F = self._base_field
            dict = self._raw['milnor']
            for mono in dict:
                dict[mono] = F(-dict[mono])
            return SteenrodAlgebraElement(dict,p)


    def _sub_(self, other):
        """
        Subtraction for elements of the Steenrod algebra.

        EXAMPLES:
            sage: A7 = SteenrodAlgebra(7)
            sage: A7.P(2,1) - A7.Q(0,3)
            P(2,1) + 6 Q_0 Q_3
        """
        return self._add_(other._neg_())


    def _mul_(self, other):
        """
        Multiplication for elements of the Steenrod algebra.

        EXAMPLES:
            sage: Sq(2) * Sq(1)
            Sq(0,1) + Sq(3)
            sage: Sq(0) * (Sq(6,2) + Sq(9,1))
            Sq(6,2) + Sq(9,1)
            sage: 1 * (Sq(6,2) + Sq(9,1))
            Sq(6,2) + Sq(9,1)
            sage: 4 * (Sq(6,2) + Sq(9,1))
            0
            sage: A5 = SteenrodAlgebra(5)
            sage: A5.P(5) * A5.P(1,2)
            3 P(0,3) + P(6,2)
            sage: A5.Q(1,2,3) * A5.Q(0,5)
            4 Q_0 Q_1 Q_2 Q_3 Q_5
            sage: A5.Q(1,2,3) * A5.Q(0,3)
            0
        """
        from sage.algebras.steenrod_algebra import SteenrodAlgebra
        p = self._prime
        if p == 2:
            from steenrod_milnor_multiplication import milnor_multiplication
        else:
            from steenrod_milnor_multiplication_odd import milnor_multiplication_odd
        F = self._base_field
        poly1 = self._raw['milnor']
        poly2 = other._raw['milnor']
        result = {}
        for mono1 in poly1:
            for mono2 in poly2:
                if len(mono1) == 0:    # multiplying by scalar multiple of one
                    if result.has_key(mono2):
                        result[mono2] = F(result[mono2] + poly1[mono1]
                                          * poly2[mono2])
                    else:
                        result[mono2] = F(poly1[mono1] * poly2[mono2])
                elif len(mono2) == 0:    # multiplying by scalar multiple of one
                    if result.has_key(mono1):
                        result[mono1] = F(result[mono1] + poly1[mono1]
                                          * poly2[mono2])
                    else:
                        result[mono1] = F(poly1[mono1] * poly2[mono2])
                else:
                    if p == 2:
                        new_dict = milnor_multiplication(mono1, mono2)
                    else:
                        new_dict = milnor_multiplication_odd(mono1, mono2, p=p)
                    for new_mono in new_dict:
                        if result.has_key(new_mono):
                            result[new_mono] = F(result[new_mono]
                                                 + new_dict[new_mono]
                                                 * poly1[mono1] * poly2[mono2])
                        else:
                            result[new_mono] = F(new_dict[new_mono]
                                                 * poly1[mono1] * poly2[mono2])
        product = SteenrodAlgebraElement(result, p=p)
        return SteenrodAlgebra(p, basis=self._basis)(product)


    def __cmp__(self,other):
        """
        Two elements are equal iff their difference is zero.

        EXAMPLES:
             sage: A5 = SteenrodAlgebra(5)
             sage: cmp(A5.P(0,1), A5.P(0,2))
             -1
             sage: cmp(A5.P(0,1), A5.pst(0,2))
             0
        """
        difference = self - other
        if len(difference._raw['milnor']) == 0:
            return 0
        else:
            return -1


    def _basis_dictionary(self,basis):
        r"""
        Dictionary of terms of the form (mono: coeff), where mono
        is a monomial in the given basis.

        INPUT:
            basis -- string, basis in which to work

        EXAMPLES:
            sage: c = Sq(2) * Sq(1)
            sage: c._basis_dictionary('milnor')
            {(0, 1): 1, (3,): 1}
            sage: c
            Sq(0,1) + Sq(3)
            sage: c._basis_dictionary('serre-cartan')
            {(2, 1): 1}
            sage: c.basis('serre-cartan')
            Sq^{2} Sq^{1}
            sage: d = Sq(0,0,1)
            sage: d._basis_dictionary('arnonc')
            {(7,): 1, (2, 5): 1, (4, 3): 1, (4, 2, 1): 1}
            sage: d.basis('arnonc')
            Sq^{2} Sq^{5} + Sq^{4} Sq^{2} Sq^{1} + Sq^{4} Sq^{3} + Sq^{7}

        At odd primes:
            sage: e = 2 * SteenrodAlgebra(3).P(1,2)
            sage: e._basis_dictionary('milnor')
            {((), (1, 2)): 2}
            sage: e
            2 P(1,2)
            sage: e._basis_dictionary('serre-cartan')
            {(0, 7, 0, 2, 0): 2, (0, 8, 0, 1, 0): 2}
            sage: e.basis('adem')
            2 P^{7} P^{2} + 2 P^{8} P^{1}

        Implementation: to compute this, take the Milnor representation
        and change bases.  Store the result in self._raw[basis], for later
        use; this way, an element only needs to be converted once.
        """
        def is_power_of_two(n):
            """
            True if and only n is a power of 2
            """
            while n != 0 and n%2 == 0:
                n = n >> 1
            return n == 1

        from steenrod_algebra import _steenrod_serre_cartan_basis_names, \
            _steenrod_milnor_basis_names, get_basis_name
        from steenrod_algebra_bases import milnor_convert

        p = self._prime
        basis_name = get_basis_name(basis, p)
        if basis_name == 'milnor':
            return self._raw['milnor']

        if basis_name.find('long') >= 0:
            basis_name = basis_name.rsplit('_', 1)[0]

        if self._raw.has_key(basis_name) and len(self._raw[basis_name])>0:
            return self._raw[basis_name]
        elif p == 2:
            dict = {}
            for mono in self._raw['milnor']:
                converted = False
                if dict.has_key(mono):
                    old_coeff = dict[mono]
                else:
                    old_coeff = 0
                new_coeff = old_coeff + self._raw['milnor'][mono]
                if len(mono) == 0:   # length 0: no conversion
                    if new_coeff != 0:
                        dict[mono] = new_coeff
                    else:
                        del dict[mono]
                    converted = True
                elif basis_name == 'serre-cartan':
                    if len(mono) == 1:   # length 1: Sq(n) = Sq^{n}, so no conversion
                        if dict.has_key(mono):
                            new_coeff = dict[mono] + self._raw['milnor'][mono]
                            if new_coeff != 0:
                                dict[mono] = new_coeff
                            else:
                                del dict[mono]
                        else:
                            dict[mono] = self._raw['milnor'][mono]
                        converted = True
                elif basis_name == 'woody':
                    if len(mono) == 1 and is_power_of_two(mono[0]):   # no conversion
                        if new_coeff != 0:
                            dict[((integer_base_2_log(mono[0]),0),)] = new_coeff
                        else:
                            del dict[((integer_base_2_log(mono[0]),0),)]
                        converted = True
                elif basis_name == 'woodz':
                    if len(mono) == 1 and is_power_of_two(mono[0]):   # no conversion
                        if new_coeff != 0:
                            dict[((integer_base_2_log(mono[0]),0),)] = new_coeff
                        else:
                            del dict[((integer_base_2_log(mono[0]),0),)]
                        converted = True
                elif basis_name == 'wall':
                    if len(mono) == 1 and is_power_of_two(mono[0]):   # no conversion
                        m = integer_base_2_log(mono[0])
                        if new_coeff != 0:
                            dict[((m,m),)] = new_coeff
                        else:
                            del dict[((m,m),)]
                        converted = True
                elif basis_name == 'arnona':
                    if len(mono) == 1 and is_power_of_two(mono[0]):   # no conversion
                        m = integer_base_2_log(mono[0])
                        if new_coeff != 0:
                            dict[((m,m),)] = new_coeff
                        else:
                            del dict[((m,m),)]
                        converted = True
                elif basis_name == 'arnonc':
                    if len(mono) == 1:   # no conversion
                        if dict.has_key(mono):
                            new_coeff = dict[mono] + self._raw['milnor'][mono]
                            if new_coeff != 0:
                                dict[mono] = new_coeff
                            else:
                                del dict[mono]
                        else:
                            dict[mono] = self._raw['milnor'][mono]
                        converted = True
                if not converted:
                    conversion = milnor_convert(  # conversion required
                        SteenrodAlgebraElement({mono: 1}),
                        basis)
                    for new_mono in conversion:
                        if dict.has_key(new_mono):
                            new_coeff = dict[new_mono] + conversion[new_mono]
                            if new_coeff != 0:
                                dict[new_mono] = new_coeff
                            else:
                                del dict[new_mono]
                        else:
                            dict[new_mono] = conversion[new_mono]
            self._raw[basis] = dict
            return dict
        else:  # p odd
            dict = {}
            for mono in self._raw['milnor']:
                converted = False
                if dict.has_key(mono):
                    old_coeff = dict[mono]
                else:
                    old_coeff = 0
                new_coeff = old_coeff + self._raw['milnor'][mono]
                if len(mono) == 0:   # length 0: no conversion
                    if new_coeff != 0:
                        dict[mono] = new_coeff
                    else:
                        del dict[mono]
                    converted = True
                elif basis in _steenrod_serre_cartan_basis_names:
                    if len(mono[0]) == 0 and len(mono[1]) == 1:
                        # length 1: Sq(n) = Sq^{n}, so no conversion
                        new_mono = (0,mono[1][0], 0)
                        if dict.has_key(new_mono):
                            new_coeff = dict[new_mono] + self._raw['milnor'][mono]
                            if new_coeff != 0:
                                dict[new_mono] = new_coeff
                            else:
                                del dict[new_mono]
                        else:
                            dict[new_mono] = self._raw['milnor'][mono]
                        converted = True
                if not converted:
                    conversion = milnor_convert(  # conversion required
                        SteenrodAlgebraElement({mono: self._raw['milnor'][mono]},
                                               p),
                        basis)
                    for new_mono in conversion:
                        if dict.has_key(new_mono):
                            new_coeff = dict[new_mono] + conversion[new_mono]
                            if new_coeff != 0:
                                dict[new_mono] = new_coeff
                            else:
                                del dict[new_mono]
                        else:
                            dict[new_mono] = conversion[new_mono]
            self._raw[basis] = dict
            return dict


    def basis(self,basis):
        r"""
        Representation of element with respect to basis.

        INPUT:
            basis -- string, basis in which to work.

        OUTPUT:
            Representation of self in given basis

        The choices for basis are:

          * 'milnor' for the Milnor basis.

          * 'serre-cartan', 'serre_cartan', 'sc', 'adem', 'admissible' for
                the Serre-Cartan basis.

          * 'wood_y' for Wood's Y basis.

          * 'wood_z' for Wood's Z basis.

          * 'wall' for Wall's basis.

          * 'wall_long' for Wall's basis, alternate representation

          * 'arnon_a' for Arnon's A basis.

          * 'arnon_a_long' for Arnon's A basis, alternate representation.

          * 'arnon_c' for Arnon's C basis.

          * 'pst', 'pst_rlex', 'pst_llex', 'pst_deg', 'pst_revz' for various
                $P^s_t$-bases.

          * 'comm', 'comm_rlex', 'comm_llex', 'comm_deg', 'comm_revz' for
                various commutator bases.

          * 'comm_long', 'comm_rlex_long', etc., for commutator bases, alternate
                representations.

        See documentation for the function 'steenrod_algebra_basis' for
        descriptions of the different bases.

        EXAMPLES:
            sage: c = Sq(2) * Sq(1)
            sage: c.basis('milnor')
            Sq(0,1) + Sq(3)
            sage: c.basis('serre-cartan')
            Sq^{2} Sq^{1}
            sage: d = Sq(0,0,1)
            sage: d.basis('arnonc')
            Sq^{2} Sq^{5} + Sq^{4} Sq^{2} Sq^{1} + Sq^{4} Sq^{3} + Sq^{7}
        """
        from sage.algebras.steenrod_algebra import SteenrodAlgebra
        return SteenrodAlgebra(p=self._prime, basis=basis)(self)


    def milnor(self):
        r"""
        Milnor representation of self.

        OUTPUT:
            Milnor representation of self.

        EXAMPLES:
            sage: A = SteenrodAlgebra(2, 'adem')
            sage: x = A (Sq(5) * Sq(2) * Sq(1)); x
            Sq^{5} Sq^{2} Sq^{1}
            sage: x.milnor()
            Sq(1,0,1) + Sq(5,1)
            """
        return self.basis('milnor')


    def serre_cartan(self):
        r"""
        Serre-Cartan representation of self.

        OUTPUT:
            Serre-Cartan representation of self.

        EXAMPLES:
            sage: x = Sq(0,1); x
            Sq(0,1)
            sage: x.serre_cartan()
            Sq^{2} Sq^{1} + Sq^{3}
            sage: x.adem()  # 'adem' is a synomym for 'serre_cartan'
            Sq^{2} Sq^{1} + Sq^{3}
            """
        return self.basis('serre-cartan')


    adem = serre_cartan


    def degree(self):
        r"""
        Degree of element.

        OUTPUT:
            degree -- None, or non-negative integer

        The degree of $\text{Sq}(i_1,i_2,i_3,...)$ is
        $i_1 + 3 i_2 + 7 i_3 + ... + (2^n - 1) i_n + ...$.
        When $p$ is odd, the degree of
        $Q_{0}^{e_0} Q_{1}^{e_1} ... P(r_1, r_2, ...)$
        is $\sum e_i (2p^i - 1) + \sum r_j (2p^j - 2)$.

        The degree of a sum is undefined (and this returns 'None'),
        unless each summand has the same degree: that is, unless the
        element is homogeneous.

        EXAMPLES:
            sage: a = Sq(1,2,1)
            sage: a.degree()
            14
            sage: for a in Sq(3) + Sq(5,1): a.degree()
            3
            8
            sage: (Sq(3) + Sq(5,1)).degree()
            Element is not homogeneous.
            sage: B = SteenrodAlgebra(3)
            sage: x = B.Q(1,4)
            sage: y = B.P(1,2,3)
            sage: x.degree()
            166
            sage: y.degree()
            192
        """
        def p_degree(m, mult=1, prime=2):
            """
            For m=(n_1, n_2, n_3, ...), Sum_i 2*n_i*(p^i - 1)
            """
            i = 0
            deg = 0
            for n in m:
                i += 1
                deg += n*mult*(prime**i - 1)
            return deg

        def q_degree(m, prime=3):
            """
            For m=(n_0, n_1, n_2, ...), Sum_i 2*p^{n_i} - 1
            """
            deg = 0
            for n in m:
                deg += 2*prime**n - 1
            return deg

        p = self._prime
        if p == 2:
            degrees = [p_degree(mono) for mono in self._raw['milnor']]
        else:
            degrees = [q_degree(mono1, prime=p)
                           + p_degree(mono2, prime=p, mult=2)
                           for (mono1, mono2) in self._raw['milnor']]
        if min(degrees) == max(degrees):
            return min(degrees)
        else:
            print "Element is not homogeneous."
            return None


    def excess(self):
        r"""
        Excess of element.

        OUTPUT:
            excess -- non-negative integer

        The excess of $\text{Sq}(a,b,c,...)$ is $a + b + c + ...$.
        When $p$ is odd, the excess of $Q_{0}^{e_0}
        Q_{1}^{e_1} ... P(r_1, r_2, ...)$ is $\sum e_i + 2 \sum r_i$.

        The excess of a linear combination of Milnor basis elements is
        the minimum of the excesses of those basis elements.

        See [Kra] for the proofs of these assertions.

        EXAMPLES:
            sage: a = Sq(1,2,3)
            sage: a.excess()
            6
            sage: (Sq(0,0,1) + Sq(4,1) + Sq(7)).excess()
            1
            sage: [m.excess() for m in (Sq(0,0,1) + Sq(4,1) + Sq(7))]
            [1, 5, 7]
            sage: [m for m in (Sq(0,0,1) + Sq(4,1) + Sq(7))]
            [Sq(0,0,1), Sq(4,1), Sq(7)]
            sage: B = SteenrodAlgebra(7)
            sage: a = B.Q(1,2,5)
            sage: b = B.P(2,2,3)
            sage: a.excess()
            3
            sage: b.excess()
            14
            sage: (a + b).excess()
            3
            sage: (a * b).excess()
            17

        REFERENCES:

            [Kra] D. Kraines, "On excess in the Milnor basis," Bull. London
                Math. Soc. 3 (1971), 363-365.
        """
        def excess_odd(mono):
            """
            Excess of mono, where mono has the form ((s0, s1, ...), (r1, r2, ...)).

            Returns the length of the first component, since that is
            the number of factors, plus twice the sum of the terms in
            the second component.
            """
            if len(mono) == 0:
                return 0
            else:
                return len(mono[0]) + 2 * sum(mono[1])

        p = self._prime
        if p == 2:
            excesses = [sum(mono) for mono in self._raw['milnor']]
        else:
            excesses = [excess_odd(mono) for mono in self._raw['milnor']]
        return min(excesses)


    def may_weight(self):
        r"""
        May's 'weight' of element.

        OUTPUT:
            weight -- non-negative integer

        If we let $F_* (A)$ be the May filtration of the Steenrod
        algebra, the weight of an element $x$ is the integer $k$ so
        that $x$ is in $F_k(A)$ and not in $F_{k+1}(A)$.  According to
        Theorem 2.6 in May's thesis [May], the weight of a Milnor
        basis element is computed as follows: first, to compute the
        weight of $P(r_1,r2, ...)$, write each $r_i$ in base
        $p$ as $r_i = \sum_j p^j r_{ij}$.  Then each nonzero binary
        digit $r_{ij}$ contributes $i$ to the weight: the weight is
        $\sum_{i,j} i r_{ij}$.  When $p$ is odd, the weight of $Q_i$
        is $i+1$, so the weight of a product $Q_{i_1} Q_{i_2} ...$ is
        equal $(i_1+1) + (i_2+1) + ...$.  Then the weight of $Q_{i_1}
        Q_{i_2} ...P(r_1,r2, ...)$ is the sum of $(i_1+1) +
        (i_2+1) + ...$ and $\sum_{i,j} i r_{ij}$.

        The weight of a sum of basis elements is the minimum of the
        weights of the summands.

        When $p=2$, we compute the weight on Milnor basis elements by
        adding up the terms in their 'height' -- see the method
        'wall_height' for documentation.  (When $p$ is odd, the height
        of an element is not defined.)

        EXAMPLES:
            sage: Sq(0).may_weight()
            0
            sage: a = Sq(4)
            sage: a.may_weight()
            1
            sage: b = Sq(4)*Sq(8) + Sq(8)*Sq(4)
            sage: b.may_weight()
            2
            sage: Sq(2,1,5).wall_height()
            [2, 3, 2, 1, 1]
            sage: Sq(2,1,5).may_weight()
            9
            sage: A5 = SteenrodAlgebra(5)
            sage: a = A5.Q(1,2,4)
            sage: b = A5.P(1,2,1)
            sage: a.may_weight()
            10
            sage: b.may_weight()
            8
            sage: (a * b).may_weight()
            18
            sage: A5.P(0,0,1).may_weight()
            3

        REFERENCES:

            [May]: J. P. May, "The cohomology of restricted Lie algebras
                and of Hopf algebras; application to the Steenrod
                algebra." Thesis, Princeton Univ., 1964.
        """
        from sage.rings.infinity import Infinity
        p = self._prime
        if self == 0:
            return Infinity
        elif self.is_unit():
            return 0
        elif p == 2:
            wt = Infinity
            for mono in self:
                wt = min(wt, sum(mono.wall_height()))
            return wt
        else: # p odd
            wt = Infinity
            for (mono1, mono2) in self._raw['milnor']:
                P_wt = 0
                index = 1
                for n in mono2:
                    P_wt += index * sum(base_p_expansion(n,p))
                    index += 1
                wt = min(wt, sum(mono1) + len(mono1) + P_wt)
            return wt


    def is_decomposable(self):
        r"""
        return True if element is decomposable, False otherwise.

        OUTPUT:
            decomposable -- boolean

        That is, if element is in the square of the augmentation
        ideal, return True; otherwise, return False.

        EXAMPLES:
            sage: a = Sq(6)
            sage: a.is_decomposable()
            True
            sage: for i in range(9):
            ...       if not Sq(i).is_decomposable():
            ...           print Sq(i)
            Sq(0)
            Sq(1)
            Sq(2)
            Sq(4)
            Sq(8)
        """
        return self.may_weight() > 1


    def wall_height(self):
        r"""
        Wall's 'height' of element.

        OUTPUT:
            height -- list of non-negative integers

        The height of an element of the mod 2 Steenrod algebra is a
        list of non-negative integers, defined as follows: if the
        element is a monomial in the generators $\text{Sq}(2^i)$, then
        the $i$th entry in the list is the number of times
        $\text{Sq}(2^i)$ appears.  For an arbitrary element, write it
        as a sum of such monomials; then its height is the maximum,
        ordered right-lexicographically, of the heights of those
        monomials.

        When $p$ is odd, the height of an element is not defined.

        According to Theorem 3 in [Wall], the height of the Milnor
        basis element $\text{Sq}(r_1, r_2, ...)$ is obtained as
        follows: write each $r_i$ in binary as $r_i = \sum_j 2^j r_{ij}$.
        Then each nonzero binary digit $r_{ij}$ contributes 1 to the
        $k$th entry in the height, for $j \leq k \leq i+j-1$.

        EXAMPLES:
            sage: Sq(0).wall_height()
            []
            sage: a = Sq(4)
            sage: a.wall_height()
            [0, 0, 1]
            sage: b = Sq(4)*Sq(8) + Sq(8)*Sq(4)
            sage: b.wall_height()
            [0, 0, 1, 1]
            sage: Sq(0,0,3).wall_height()
            [1, 2, 2, 1]

        REFERENCES:

            [Wall]: C. T. C. Wall, "Generators and relations for the
                Steenrod algebra," Ann. of Math. (2) \textbf{72} (1960),
                429--444.
        """
        def mono_degree(m):
            """
            If m=(n1,n2,n3, ...), returns the degree of Sq(n1,n2,n3,...).
            That is, it returns sum n_i (2^i - 1).
            """
            i = 0
            deg = 0
            for n in m:
                i += 1
                deg += n*(2**i - 1)
            return deg

        if self._prime > 2:
            raise NotImplementedError, "Wall height is not defined at odd primes."
        if self == 0 or self == 1:
            return []
        result = []
        for r in self._raw['milnor']:
            h = [0]*(1+mono_degree(r))
            i = 1
            for x in r:
                if x > 0:
                    for j in range(1+integer_base_2_log(x)):
                        if (2**j & x) != 0:
                            for k in range(j,i+j):
                                h[k] += 1
                i=i+1
            h.reverse()
            result = max(h, result)
        result.reverse()
        return check_and_trim(result)


    def antipode(self):
        r"""
        Antipode of element.

        OUTPUT:
            antipode -- element of the Steenrod algebra

        Algorithm: according to a result of Milnor's, the antipode of
        $\text{Sq}(n)$ is the sum of all of the Milnor basis elements
        in dimension $n$.  So: convert the element to the Serre-Cartan
        basis and use this formula for the antipode of $\text{Sq}(n)$,
        together with the fact that the antipode is an
        antihomomorphism: if we call the antipode $c$, then $c(ab) =
        c(b) c(a)$.

        At odd primes, a similar method is used: the antipode of
        $P(n)$ is the sum of the Milnor basis elements in
        dimension $n*2(p-1)$, and the antipode of $\beta = Q_0$ is
        $-Q_0$.  So convert to the Serre-Cartan basis, as in the $p=2$
        case.

        EXAMPLES:
            sage: d = Sq(0,0,1); d
            Sq(0,0,1)
            sage: d.antipode()
            Sq(0,0,1)
            sage: Sq(4).antipode()
            Sq(1,1) + Sq(4)
            sage: (Sq(4) * Sq(2)).antipode()
            Sq(6)
            sage: A3 = SteenrodAlgebra(3)
            sage: A3.P(2).antipode()
            P(2)
            sage: A3.P(2,1).antipode()
            2 P(2,1)
            sage: a = SteenrodAlgebra(7).P(3,1)
            sage: a.antipode()
            4 P(3,1) + 5 P(11)

        Applying the antipode twice returns the original element:
            sage: y = Sq(8)*Sq(4)
            sage: y == (y.antipode()).antipode()
            True
        """
        def sum_of_basis(n,p):
            """
            Antipode of P(n) (i.e., of Sq(n) when p=2).

            INPUT:
                n -- integer
                p -- positive prime number

            OUTPUT:
                elt -- element of the Steenrod algebra

            This returns the sum of all of the elements P(...) in the
            Milnor basis in dimension $n$ at the prime p
            """
            from steenrod_algebra_bases import steenrod_algebra_basis
            return sum(steenrod_algebra_basis(n,'milnor',p=p))
        from sage.algebras.steenrod_algebra import SteenrodAlgebra
        from steenrod_algebra_bases import milnor_convert
        result = 0
        p = self._prime
        if p == 2:
            for mono in self._basis_dictionary('serre_cartan'):
                antipode = Sq(0)
                for n in mono:
                    antipode = sum_of_basis(n, p) * antipode
                result = result + antipode
        else:
            from sage.misc.functional import is_even
            for mono in self._basis_dictionary('serre_cartan'):
                antipode = SteenrodAlgebra(p).P(0)
                index = 0
                for n in mono:
                    if is_even(index) and n != 0:
                        antipode = -SteenrodAlgebra(p).Q(0) * antipode
                    else:
                        antipode = sum_of_basis(n*2*(p-1),p) * antipode
                    index += 1
                result = result + antipode
        return result


    def _repr_(self):
        """
        String representation of element.

        OUTPUT:
            string

        The string depends on the basis over which the element is defined.

        EXAMPLES:
            sage: A7 = SteenrodAlgebra(7)
            sage: x = A7.Q(0,3) * A7.P(2,2)
            sage: x._repr_()
            'Q_0 Q_3 P(2,2)'
            sage: x
            Q_0 Q_3 P(2,2)
            sage: a = Sq(0,0,2)
            sage: a
            Sq(0,0,2)
            sage: A2_adem = SteenrodAlgebra(2,'admissible')
            sage: A2_adem(a)
            Sq^{8} Sq^{4} Sq^{2} + Sq^{9} Sq^{4} Sq^{1} + Sq^{10} Sq^{3} Sq^{1} +
            Sq^{10} Sq^{4} + Sq^{11} Sq^{2} Sq^{1} + Sq^{12} Sq^{2} + Sq^{13} Sq^{1}
            + Sq^{14}
            sage: SteenrodAlgebra(2, 'woodz')(a)
            Sq^{6} Sq^{7} Sq^{1} + Sq^{14} + Sq^{4} Sq^{7} Sq^{3} + Sq^{4} Sq^{7}
            Sq^{2} Sq^{1} + Sq^{12} Sq^{2} + Sq^{8} Sq^{6} + Sq^{8} Sq^{4} Sq^{2}
            sage: SteenrodAlgebra(2, 'arnonc')(a)
            Sq^{4} Sq^{2} Sq^{8} + Sq^{4} Sq^{4} Sq^{6} + Sq^{4} Sq^{6} Sq^{4} +
            Sq^{6} Sq^{8} + Sq^{8} Sq^{4} Sq^{2} + Sq^{8} Sq^{6}
            sage: SteenrodAlgebra(2, 'pst_llex')(a)
            P^{1}_{3}
            sage: SteenrodAlgebra(2, 'comm_revz')(a)
            c_{0,1} c_{1,1} c_{0,3} c_{2,1} + c_{0,2} c_{0,3} c_{2,1} + c_{1,3}
        """
        if len(self._raw['milnor']) == 0:
            return "0"
        else:
            return string_rep(self)


    def _latex_(self):
        """
        LaTeX representation of element.

        OUTPUT:
            string

        The string depends on the basis over which the element is defined.

        For any element x in the Steenrod algebra, use 'view(x)' to
        see the typeset LaTeX representation.

        EXAMPLES:
            sage: A7 = SteenrodAlgebra(7)
            sage: x = A7.Q(0,3) * A7.P(2,2)
            sage: x._latex_()
            'Q_{0} Q_{3} \\mathcal{P}(2,2)'
            sage: latex(x)
            Q_{0} Q_{3} \mathcal{P}(2,2)
            sage: b = Sq(0,2)
            sage: b.basis('adem')._latex_()
            '\\text{Sq}^{4} \\text{Sq}^{2} + \\text{Sq}^{5} \\text{Sq}^{1} +
            \\text{Sq}^{6}'
            sage: b.basis('woody')._latex_()
            '\\text{Sq}^{2} \\text{Sq}^{3} \\text{Sq}^{1} + \\text{Sq}^{6} +
            \\text{Sq}^{4} \\text{Sq}^{2}'
            sage: SteenrodAlgebra(2, 'arnona')(b)._latex_()
            'X^{1}_{1} X^{2}_{2}  + X^{2}_{1}'
        """
        if len(self._raw['milnor']) == 0:
            return "0"
        else:
            return string_rep(self,LaTeX=True)


    def __iter__(self):
        """
        Iterator for looping through summands in an element of the
        Steenrod algebra.

        EXAMPLES:
            sage: z = Sq(0,0,1) + Sq(4,1) + Sq(7)
            sage: [m for m in z]
            [Sq(0,0,1), Sq(4,1), Sq(7)]
            sage: [m.excess() for m in z]
            [1, 5, 7]
            sage: for m in z: m * Sq(2)
            Sq(2,0,1)
            Sq(0,3) + Sq(6,1)
            Sq(3,2)
            sage: a = SteenrodAlgebra(5).P(5,5)
            sage: a * a
            P(3,6,1) + 2 P(4,11) + P(9,5,1) + 4 P(10,10)
            sage: for m in a * a: m
            P(3,6,1)
            2 P(4,11)
            P(9,5,1)
            4 P(10,10)

        This loops through the summands in the Milnor basis
        representation of the element.  The element w defined below is
        a single monomial in the Serre-Cartan basis, but a sum of four
        monomials in the Milnor basis:
            sage: w = Sq(4) * Sq(2) * Sq(1)
            sage: A = SteenrodAlgebra(2, 'adem')
            sage: w = A(Sq(4) * Sq(2) * Sq(1)); w
            Sq^{4} Sq^{2} Sq^{1}
            sage: for m in w: m.basis('adem')
            Sq^{4} Sq^{2} Sq^{1} + Sq^{5} Sq^{2} + Sq^{6} Sq^{1} + Sq^{7}
            Sq^{5} Sq^{2} + Sq^{7}
            Sq^{6} Sq^{1} + Sq^{7}
            Sq^{7}
            sage: w.milnor()
            Sq(0,0,1) + Sq(1,2) + Sq(4,1) + Sq(7)
            sage: for m in w: m
            Sq(0,0,1)
            Sq(1,2)
            Sq(4,1)
            Sq(7)
        """
        for m in sorted(self._raw['milnor'].keys()):
            yield SteenrodAlgebraElement({m: self._raw['milnor'][m]},
                                         p = self._prime)


    def additive_order(self):
        """
        The additive order of any element of the mod p Steenrod algebra is p.

        OUTPUT:
            order -- positive prime number

        EXAMPLES:
            sage: z = Sq(4) + Sq(6) + Sq(0)
            sage: z.additive_order()
            2
        """
        return self._prime


def Sq(*nums):
    """
    Milnor element Sq(a,b,c,...).

    INPUT:
        a, b, c, ... -- non-negative integers

    OUTPUT:
        element of the Steenrod algebra

    This returns the Milnor basis element $\text{Sq}(a, b, c, ...)$.

    EXAMPLES:
        sage: Sq(5)
        Sq(5)
        sage: Sq(5) + Sq(2,1) + Sq(5)  # addition is mod 2:
        Sq(2,1)
        sage: (Sq(4,3) + Sq(7,2)).degree()
        13

    Entries must be non-negative integers; otherwise, an error results.

    This function is a good way to define elements of the Steenrod algebra.
    """
    dict = {nums: 1}
    return SteenrodAlgebraElement(dict, p=2)


def pst(s,t,p=2):
    """
    The Margolis element $P^s_t$.

    INPUT:
        s -- non-negative integer
        t -- positive integer
        p -- positive prime number (optional, default 2)

    OUTPUT:
        element of the Steenrod algebra

    This returns the Margolis element $P^s_t$ of the mod p Steenrod
    algebra: the element equal to $P(0,0,...,0,p^s)$, where
    the $p^s$ is in position $t$.

    EXAMPLES:
        sage: from sage.algebras.steenrod_algebra_element import pst
        sage: pst(3,5)
        Sq(0,0,0,0,8)
        sage: pst(1,2) + Sq(4)*Sq(2) + Sq(2)*Sq(4)
        0
        sage: pst(3,5,5)
        P(0,0,0,0,125)
        sage: pst(3,5,p=5)
        P(0,0,0,0,125)
    """
    from sage.algebras.steenrod_algebra import SteenrodAlgebra
    return SteenrodAlgebra(p).pst(s,t)


def degree(x):
    r"""
    Degree of x.

    INPUT:
        x -- element of the Steenrod algebra

    OUTPUT:
        degree -- non-negative integer or None

    The degree of $\text{Sq}(i_1,i_2,i_3,...)$ is
    $i_1 + 3 i_2 + 7 i_3 + ... + (2^n - 1) i_n + ...$.
    When $p$ is odd, the degree of
    $Q_{0}^{e_0} Q_{1}^{e_1} ... P(r_1, r_2, ...)$
    is $\sum e_i (2p^i - 1) + \sum r_j (2p^j - 2)$.

    The degree of a sum is undefined (and this function returns None),
    unless each summand has the same degree: that is, unless the
    element is homogeneous.

    EXAMPLES:
        sage: from sage.algebras.steenrod_algebra_element import degree
        sage: a = Sq(1,2,1)
        sage: degree(a)
        14
        sage: degree(Sq(3) + Sq(5,1))
        Element is not homogeneous.
        sage: B = SteenrodAlgebra(3)
        sage: x = B.Q(1,4)
        sage: y = B.P(1,2,3)
        sage: degree(x)
        166
        sage: degree(y)
        192
    """
    return x.degree()


def excess(x):
    r"""
    Excess of x.

    INPUT:
        x -- element of the Steenrod algebra

    OUTPUT:
        excess -- non-negative integer

    The excess of $\text{Sq}(a,b,c,...)$ is $a + b + c + ...$.
    When $p$ is odd, the excess of
    $Q_{0}^{e_0} Q_{1}^{e_1} ... P(r_1, r_2, ...)$
    is $\sum e_i + 2 \sum r_i$.

    The excess of a linear combination of Milnor basis elements is
    the minimum of the excesses of those basis elements.

    EXAMPLES:
        sage: from sage.algebras.steenrod_algebra_element import excess
        sage: a = Sq(1,2,3)
        sage: excess(a)
        6
        sage: excess(Sq(0,0,1) + Sq(4,1) + Sq(7))
        1
        sage: [excess(m) for m in (Sq(0,0,1) + Sq(4,1) + Sq(7))]
        [1, 5, 7]
        sage: B = SteenrodAlgebra(7)
        sage: a = B.Q(1,2,5)
        sage: b = B.P(2,2,3)
        sage: excess(a)
        3
        sage: excess(b)
        14
        sage: excess(a + b)
        3
        sage: excess(a * b)
        17
    """
    return x.excess()


def milnor(x):
    r"""
    Milnor representation of x.

    INPUT:
        x -- element of the Steenrod algebra

    OUTPUT:
        Milnor representation of x

    EXAMPLES:
        sage: from sage.algebras.steenrod_algebra_element import milnor
        sage: x = Sq(5) * Sq(2) * Sq(1); x.adem()
        Sq^{5} Sq^{2} Sq^{1}
        sage: milnor(x)
        Sq(1,0,1) + Sq(5,1)
    """
    return x.basis('milnor')


def serre_cartan(x):
    r"""
    Serre-Cartan representation of x.

    INPUT:
        x -- element of the Steenrod algebra

    OUTPUT:
        Serre-Cartan representation of x

    EXAMPLES:
        sage: from sage.algebras.steenrod_algebra_element import serre_cartan
        sage: x = Sq(3,2); x
        Sq(3,2)
        sage: serre_cartan(x)
        Sq^{7} Sq^{2}
    """
    return x.basis('adem')

adem = serre_cartan
admissible = serre_cartan

## string representations

def string_rep(element, LaTeX=False, sort=True):
    """
    String representation of element.

    INPUT:
        element -- element of the Steenrod algebra
        LaTeX -- boolean (optional, default False), if True, output LaTeX string
        sort -- boolean (optional, default True), if True, sort output

    OUTPUT:
        string -- string representation of element in current basis

    If LaTeX is True, output a string suitable for LaTeX; otherwise,
    output a plain string.  If sort is True, sort element left
    lexicographically; otherwise, no sorting is done, and so the
    order in which the summands are printed may be unpredictable.

    EXAMPLES:
        sage: from sage.algebras.steenrod_algebra_element import string_rep
        sage: a = Sq(0,0,2)
        sage: A = SteenrodAlgebra(2, 'admissible')
        sage: string_rep(A(a))
        'Sq^{8} Sq^{4} Sq^{2} + Sq^{9} Sq^{4} Sq^{1} + Sq^{10} Sq^{3} Sq^{1} +
        Sq^{10} Sq^{4} + Sq^{11} Sq^{2} Sq^{1} + Sq^{12} Sq^{2} + Sq^{13} Sq^{1}
        + Sq^{14}'
        sage: b = Sq(0,2)
        sage: string_rep(A(b),LaTeX=True)
        '\\text{Sq}^{4} \\text{Sq}^{2} + \\text{Sq}^{5} \\text{Sq}^{1} +
        \\text{Sq}^{6}'
        sage: A_wood_z = SteenrodAlgebra(2, 'woodz')
        sage: string_rep(A_wood_z(a))
        'Sq^{6} Sq^{7} Sq^{1} + Sq^{14} + Sq^{4} Sq^{7} Sq^{3} + Sq^{4} Sq^{7}
        Sq^{2} Sq^{1} + Sq^{12} Sq^{2} + Sq^{8} Sq^{6} + Sq^{8} Sq^{4} Sq^{2}'
        sage: string_rep(SteenrodAlgebra(2, 'arnonc')(a), sort=False)
        'Sq^{4} Sq^{4} Sq^{6} + Sq^{6} Sq^{8} + Sq^{4} Sq^{2} Sq^{8} + Sq^{4}
        Sq^{6} Sq^{4} + Sq^{8} Sq^{4} Sq^{2} + Sq^{8} Sq^{6}'
        sage: string_rep(SteenrodAlgebra(2, 'arnonc')(a))
        'Sq^{4} Sq^{2} Sq^{8} + Sq^{4} Sq^{4} Sq^{6} + Sq^{4} Sq^{6} Sq^{4} +
        Sq^{6} Sq^{8} + Sq^{8} Sq^{4} Sq^{2} + Sq^{8} Sq^{6}'
        sage: string_rep(SteenrodAlgebra(2, 'pst_llex')(a))
        'P^{1}_{3}'
        sage: Ac = SteenrodAlgebra(2, 'comm_revz')
        sage: string_rep(Ac(a),LaTeX=True,sort=False)
        'c_{0,2} c_{0,3} c_{2,1} + c_{1,3} + c_{0,1} c_{1,1} c_{0,3} c_{2,1}'
        sage: string_rep(Ac(a),LaTeX=True)
        'c_{0,1} c_{1,1} c_{0,3} c_{2,1} + c_{0,2} c_{0,3} c_{2,1} + c_{1,3}'
        sage: string_rep(a)
        'Sq(0,0,2)'
        sage: string_rep(a,LaTeX=True)
        '\\text{Sq}(0,0,2)'

    Some odd primary examples:
        sage: A5 = SteenrodAlgebra(5)
        sage: a = A5.P(5,1); b = A5.Q(0,1,3)
        sage: string_rep(b)
        'Q_0 Q_1 Q_3'
        sage: string_rep(a, LaTeX=True)
        '\\mathcal{P}(5,1)'
        sage: A5sc = SteenrodAlgebra(5, 'serre-cartan')
        sage: string_rep(A5sc(a))
        'P^{10} P^{1} + 4 P^{11}'
    """
    if len(element._raw['milnor']) == 0:
        return "0"
    p = element._prime

    basis = element._basis
    dict = element._basis_dictionary(basis)

    if basis == 'milnor':
        mono_to_string = milnor_mono_to_string
    elif basis == 'serre-cartan':
        mono_to_string = serre_cartan_mono_to_string
    elif basis.find('wood') >= 0:
        mono_to_string = wood_mono_to_string
    elif basis == 'wall':
        mono_to_string = wall_mono_to_string
    elif basis == 'wall_long':
        mono_to_string = wall_long_mono_to_string
    elif basis == 'arnona':
        mono_to_string = arnonA_mono_to_string
    elif basis == 'arnona_long':
        mono_to_string = arnonA_long_mono_to_string
    elif basis == 'arnonc':
        mono_to_string = serre_cartan_mono_to_string
    elif basis.find('pst') >= 0:
        mono_to_string = pst_mono_to_string
    elif basis.find('comm') >= 0 and basis.find('long') >= 0:
        mono_to_string = comm_long_mono_to_string
    elif basis.find('comm') >= 0:
        mono_to_string = comm_mono_to_string

    output = ""
    if sort:
        sorted_list = sorted(dict.keys())
    else:
        sorted_list = dict
    for mono in sorted_list:
        if dict[mono] != 1:
            coeff = str(dict[mono]) + " "
        else:
            coeff = ""
        output = output + coeff + mono_to_string(mono, LaTeX,
                                                 p=element._prime) + " + "
    return output.strip(" +")


def milnor_mono_to_string(mono,LaTeX=False,p=2):
    """
    String representation of element of the Milnor basis.

    This is used by the _repr_ and _latex_ methods.

    INPUT:
        mono -- if $p=2$, tuple of non-negative integers (a,b,c,...);
                if $p>2$, pair of tuples of non-negative integers
                ((e0, e1, e2, ...), (r1, r2, ...))
        LaTeX -- boolean (optional, default False), if true, output LaTeX string
        p -- positive prime number (optional, default 2)

    OUTPUT:
        rep -- string

    This returns a string like 'Sq(a,b,c,...)' when p=2, or a
    string like 'Q_e0 Q_e1 Q_e2 ... P(r1, r2, ...)' when p is odd.

    EXAMPLES:
        sage: from sage.algebras.steenrod_algebra_element import milnor_mono_to_string
        sage: milnor_mono_to_string((1,2,3,4))
        'Sq(1,2,3,4)'
        sage: milnor_mono_to_string((1,2,3,4),LaTeX=True)
        '\\text{Sq}(1,2,3,4)'
        sage: milnor_mono_to_string(((1,0), (2,3,1)), p=3)
        'Q_1 Q_0 P(2,3,1)'
        sage: milnor_mono_to_string(((1,0), (2,3,1)), LaTeX=True, p=3)
        'Q_{1} Q_{0} \\mathcal{P}(2,3,1)'

    The empty tuple represents the unit element Sq(0) (or P(0) at an odd prime):
        sage: milnor_mono_to_string(())
        'Sq(0)'
        sage: milnor_mono_to_string((), p=5)
        'P(0)'
    """
    if LaTeX:
        if p == 2:
            sq = "\\text{Sq}"
            P = "\\text{Sq}"
        else:
            P = "\\mathcal{P}"
    else:
        if p == 2:
            sq = "Sq"
            P = "Sq"
        else:
            P = "P"
    if len(mono) == 0 or (p > 2 and len(mono[0]) + len(mono[1]) == 0):
        return P + "(0)"
    else:
        if p == 2:
            string = sq + "(" + str(mono[0])
            for n in mono[1:]:
                string = string + "," + str(n)
            string = string + ")"
        else:
            string = ""
            if len(mono[0]) > 0:
                for e in mono[0]:
                    if LaTeX:
                        string = string + "Q_{" + str(e) + "} "
                    else:
                        string = string + "Q_" + str(e) + " "
            if len(mono[1]) > 0:
                string = string + P + "(" + str(mono[1][0])
                for n in mono[1][1:]:
                    string = string + "," + str(n)
                string = string + ")"
        return string.strip(" ")


def serre_cartan_mono_to_string(mono,LaTeX=False,p=2):
    r"""
    String representation of element of the Serre-Cartan basis.

    This is used by the _repr_ and _latex_ methods.

    INPUT:
        mono -- tuple of positive integers (a,b,c,...) when $p=2$, or tuple
                (e0, n1, e1, n2, ...) when $p>2$, where each ei is 0 or 1, and
                each ni is positive
        LaTeX -- boolean (optional, default False), if true, output LaTeX string
        p -- positive prime number (optional, default 2)

    OUTPUT:
        rep -- string

    This returns a string like '$Sq^{a} Sq^{b} Sq^{c} ...$ when $p=2$,
    or a string like $\beta^{e0} P^{n1} \beta^{e1} P^{n2} ...$ when $p$
    is odd.

    EXAMPLES:
        sage: from sage.algebras.steenrod_algebra_element import serre_cartan_mono_to_string
        sage: serre_cartan_mono_to_string((1,2,3,4))
        'Sq^{1} Sq^{2} Sq^{3} Sq^{4}'
        sage: serre_cartan_mono_to_string((1,2,3,4),LaTeX=True)
        '\\text{Sq}^{1} \\text{Sq}^{2} \\text{Sq}^{3} \\text{Sq}^{4}'
        sage: serre_cartan_mono_to_string((0,5,1,1,0), p=3)
        'P^{5} beta P^{1}'
        sage: serre_cartan_mono_to_string((0,5,1,1,0), p=3, LaTeX=True)
        '\\mathcal{P}^{5} \\beta \\mathcal{P}^{1}'

    The empty tuple represents the unit element $Sq^0$ (or $P^0$ at an odd prime):
        sage: serre_cartan_mono_to_string(())
        'Sq^{0}'
        sage: serre_cartan_mono_to_string((), p=7)
        'P^{0}'
    """
    if LaTeX:
        if p == 2:
            sq = "\\text{Sq}"
            P = "\\text{Sq}"
        else:
            P = "\\mathcal{P}"
    else:
        if p == 2:
            sq = "Sq"
            P = "Sq"
        else:
            P = "P"
    if len(mono) == 0:
        return P + "^{0}"
    else:
        if p == 2:
            string = ""
            for n in mono:
                string = string + sq + "^{" + str(n) + "} "
        else:
            string = ""
            index = 0
            for n in mono:
                from sage.misc.functional import is_even
                if is_even(index):
                    if n == 1:
                        if LaTeX:
                            string = string + "\\beta "
                        else:
                            string = string + "beta "
                else:
                    string = string + P + "^{" + str(n) + "} "
                index += 1
        return string.strip(" ")


def wood_mono_to_string(mono,LaTeX=False,p=2):
    """
    String representation of element of Wood's Y and Z bases.

    This is used by the _repr_ and _latex_ methods.

    INPUT:
        mono -- tuple of pairs of non-negative integers (s,t)

    OUTPUT:
        string -- concatenation of '$Sq^{2^s (2^{t+1}-1)}$' for each pair (s,t)

    EXAMPLES:
        sage: from sage.algebras.steenrod_algebra_element import wood_mono_to_string
        sage: wood_mono_to_string(((1,2),(3,0)))
        'Sq^{14} Sq^{8}'
        sage: wood_mono_to_string(((1,2),(3,0)),LaTeX=True)
        '\\text{Sq}^{14} \\text{Sq}^{8}'

    The empty tuple represents the unit element Sq(0):
        sage: wood_mono_to_string(())
        'Sq(0)'
    """
    if LaTeX:
        sq = "\\text{Sq}"
    else:
        sq = "Sq"
    if len(mono) == 0:
        return sq + "(0)"
    else:
        string = ""
        for (s,t) in mono:
            string = string + sq + "^{" + \
                str(2**s * (2**(t+1)-1)) + "} "
        return string.strip(" ")


def wall_mono_to_string(mono,LaTeX=False,p=2):
    """
    String representation of element of Wall's basis.

    This is used by the _repr_ and _latex_ methods.

    INPUT:
        mono -- tuple of pairs of non-negative integers (m,k) with $m >= k$

    OUTPUT:
        string -- concatenation of '$Q^{m}_{k}$' for each pair (m,k)

    EXAMPLES:
        sage: from sage.algebras.steenrod_algebra_element import wall_mono_to_string
        sage: wall_mono_to_string(((1,2),(3,0)))
        'Q^{1}_{2} Q^{3}_{0}'
        sage: wall_mono_to_string(((1,2),(3,0)),LaTeX=True)
        'Q^{1}_{2} Q^{3}_{0}'

    The empty tuple represents the unit element Sq(0):
        sage: wall_mono_to_string(())
        'Sq(0)'
    """
    if LaTeX:
        sq = "\\text{Sq}"
    else:
        sq = "Sq"
    if len(mono) == 0:
        return sq + "(0)"
    else:
        string = ""
        for (m,k) in mono:
            string = string + "Q^{" + str(m) + "}_{" \
                + str(k) + "} "
        return string.strip(" ")


def wall_long_mono_to_string(mono,LaTeX=False,p=2):
    """
    Alternate string representation of element of Wall's basis.

    This is used by the _repr_ and _latex_ methods.

    INPUT:
        mono -- tuple of pairs of non-negative integers (m,k) with $m >= k$

    OUTPUT:
        string -- concatenation of terms of the form '$Sq^(2^m)$'

    EXAMPLES:
        sage: from sage.algebras.steenrod_algebra_element import wall_long_mono_to_string
        sage: wall_long_mono_to_string(((1,2),(3,0)))
        'Sq^{1} Sq^{2} Sq^{4} Sq^{8}'
        sage: wall_long_mono_to_string(((1,2),(3,0)),LaTeX=True)
        '\\text{Sq}^{1} \\text{Sq}^{2} \\text{Sq}^{4} \\text{Sq}^{8}'

    The empty tuple represents the unit element Sq(0):
        sage: wall_long_mono_to_string(())
        'Sq(0)'
    """
    if LaTeX:
        sq = "\\text{Sq}"
    else:
        sq = "Sq"
    if len(mono) == 0:
        return sq + "(0)"
    else:
        string = ""
        for (m,k) in mono:
            for i in range(k,m+1):
                string = string + sq + "^{" + str(2**i) + "} "
        return string.strip(" ")


def arnonA_mono_to_string(mono,LaTeX=False,p=2):
    """
    String representation of element of Arnon's A basis.

    This is used by the _repr_ and _latex_ methods.

    INPUT:
        mono -- tuple of pairs of non-negative integers (m,k) with $m >= k$

    OUTPUT:
        string -- concatenation of '$X^{m}_{k}$' for each pair (m,k)

    EXAMPLES:
        sage: from sage.algebras.steenrod_algebra_element import arnonA_mono_to_string
        sage: arnonA_mono_to_string(((1,2),(3,0)))
        'X^{1}_{2} X^{3}_{0}'
        sage: arnonA_mono_to_string(((1,2),(3,0)),LaTeX=True)
        'X^{1}_{2} X^{3}_{0}'

    The empty tuple represents the unit element Sq(0):
        sage: arnonA_mono_to_string(())
        'Sq(0)'
    """
    if LaTeX:
        sq = "\\text{Sq}"
    else:
        sq = "Sq"
    if len(mono) == 0:
        return sq + "(0)"
    else:
        string = ""
        for (m,k) in mono:
            string = string + "X^{" + str(m) + "}_{" \
                + str(k) + "} "
        return string.strip(" ")


def arnonA_long_mono_to_string(mono,LaTeX=False,p=2):
    """
    Alternate string representation of element of Arnon's A basis.

    This is used by the _repr_ and _latex_ methods.

    INPUT:
        mono -- tuple of pairs of non-negative integers (m,k) with $m >= k$

    OUTPUT:
        string -- concatenation of strings of the form '$Sq(2^m)$'

    EXAMPLES:
        sage: from sage.algebras.steenrod_algebra_element import arnonA_long_mono_to_string
        sage: arnonA_long_mono_to_string(((1,2),(3,0)))
        'Sq^{8} Sq^{4} Sq^{2} Sq^{1}'
        sage: arnonA_long_mono_to_string(((1,2),(3,0)),LaTeX=True)
        '\\text{Sq}^{8} \\text{Sq}^{4} \\text{Sq}^{2} \\text{Sq}^{1}'

    The empty tuple represents the unit element Sq(0):
        sage: arnonA_long_mono_to_string(())
        'Sq(0)'
    """
    if LaTeX:
        sq = "\\text{Sq}"
    else:
        sq = "Sq"
    if len(mono) == 0:
        return sq + "(0)"
    else:
        string = ""
        for (m,k) in mono:
            for i in range(m,k-1,-1):
                string = string + sq + "^{" + str(2**i) + "} "
        return string.strip(" ")


def pst_mono_to_string(mono,LaTeX=False,p=2):
    r"""
    String representation of element of a $P^s_t$-basis.

    This is used by the _repr_ and _latex_ methods.

    INPUT:
        mono -- tuple of pairs of integers (s,t) with $s >= 0$, $t > 0$

    OUTPUT:
        string -- concatenation of '$P^{s}_{t}$' for each pair (s,t)

    EXAMPLES:
        sage: from sage.algebras.steenrod_algebra_element import pst_mono_to_string
        sage: pst_mono_to_string(((1,2),(0,3)))
        'P^{1}_{2} P^{0}_{3}'
        sage: pst_mono_to_string(((1,2),(0,3)),LaTeX=True)
        'P^{1}_{2} P^{0}_{3}'

    The empty tuple represents the unit element Sq(0):
        sage: pst_mono_to_string(())
        'Sq(0)'
    """
    if LaTeX:
        sq = "\\text{Sq}"
    else:
        sq = "Sq"
    if len(mono) == 0:
        return sq + "(0)"
    else:
        string = ""
        for (s,t) in mono:
            string = string + "P^{" + str(s) + "}_{" \
                + str(t) + "} "
        return string.strip(" ")


def comm_mono_to_string(mono,LaTeX=False,p=2):
    r"""
    String representation of element of a commutator basis.

    This is used by the _repr_ and _latex_ methods.

    INPUT:
        mono -- tuple of pairs of integers (s,t) with $s >= 0$, $t > 0$

    OUTPUT:
        string -- concatenation of '$c_{s,t}$' for each pair (s,t)

    EXAMPLES:
        sage: from sage.algebras.steenrod_algebra_element import comm_mono_to_string
        sage: comm_mono_to_string(((1,2),(0,3)))
        'c_{1,2} c_{0,3}'
        sage: comm_mono_to_string(((1,2),(0,3)),LaTeX=True)
        'c_{1,2} c_{0,3}'

    The empty tuple represents the unit element Sq(0):
        sage: comm_mono_to_string(())
        'Sq(0)'
    """
    if LaTeX:
        sq = "\\text{Sq}"
    else:
        sq = "Sq"
    if len(mono) == 0:
        return sq + "(0)"
    else:
        string = ""
        for (s,t) in mono:
            string = string + "c_{" + str(s) + "," \
                + str(t) + "} "
        return string.strip(" ")


def comm_long_mono_to_string(mono,LaTeX=False,p=2):
    r"""
    Alternate string representation of element of a commutator basis.

    Okay in low dimensions, but gets unwieldy as the dimension increases.

    INPUT:
        mono -- tuple of pairs of integers (s,t) with $s >= 0$, $t > 0$

    OUTPUT:
        string -- concatenation of '$s_{2^s ... 2^(s+t-1)}$' for each pair (s,t)

    EXAMPLES:
        sage: from sage.algebras.steenrod_algebra_element import comm_long_mono_to_string
        sage: comm_long_mono_to_string(((1,2),(0,3)))
        's_{24} s_{124}'
        sage: comm_long_mono_to_string(((1,2),(0,3)),LaTeX=True)
        's_{24} s_{124}'

    The empty tuple represents the unit element Sq(0):
        sage: comm_long_mono_to_string(())
        'Sq(0)'
    """
    if LaTeX:
        sq = "\\text{Sq}"
    else:
        sq = "Sq"
    if len(mono) == 0:
        return sq + "(0)"
    else:
        string = ""
        for (s,t) in mono:
            if s + t > 4:
                comma = ","
            else:
                comma = ""
            string = string + "s_{"
            for i in range(t):
                string = string + str(2**(s+i)) + comma
            string = string.strip(",") + "} "
        return string.strip(" ")
