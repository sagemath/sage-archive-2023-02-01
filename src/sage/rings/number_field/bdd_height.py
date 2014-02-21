import copy
import itertools
import sage.rings.real_mpfr
from sage.rings.number_field.unit_group import UnitGroup
from sage.modules.free_module_element import vector
from sage.matrix.constructor import column_matrix
from sage.rings.rational_field import QQ
from sage.functions.other import ceil
from sage.geometry.polyhedron.constructor import *
from sage.structure.proof.all import number_field

def bdd_norm_pr_gens_iq(K, norm_list):
    r"""
    Compute generators for all principal ideals in an imaginary quadratic field
    ``K`` whose norms are in ``norm_list``.

    The only keys for the output dictionary are integers n appearing in
    ``norm_list``.

    The function will only be called with ``K`` an imaginary quadratic field.
    
    The function will return a dictionary for other number fields, but it may be
    incorrect.

    INPUT:

    - ``K`` - an imaginary quadratic number field
    - ``norm_list`` - a list of positive integers

    OUTPUT:

    - a dictionary, keyed by norm

    EXAMPLES:

    In Q(i), there is one principal ideal of norm 4, two principal ideals of
    norm 5, but no principal ideals of norm 7::

        sage: from sage.rings.number_field.bdd_height import bdd_norm_pr_gens_iq
        sage: K.<g> = NumberField(x^2 + 1)
        sage: L = range(10)
        sage: bdd_pr_ideals = bdd_norm_pr_gens_iq(K, L)
        sage: bdd_pr_ideals[4]
        [2]
        sage: bdd_pr_ideals[5]
        [g + 2, g - 2]
        sage: bdd_pr_ideals[7]
        []

    There are no ideals in the ring of integers with negative norm::

        sage: from sage.rings.number_field.bdd_height import bdd_norm_pr_gens_iq
        sage: K.<g> = NumberField(x^2 + 10)
        sage: L = range(-5,-1)
        sage: bdd_pr_ideals = bdd_norm_pr_gens_iq(K,L)
        sage: bdd_pr_ideals
        {-5: [], -4: [], -3: [], -2: []}

    Calling a key that is not in the input ``norm_list`` raises a KeyError::

        sage: from sage.rings.number_field.bdd_height import bdd_norm_pr_gens_iq
        sage: K.<g> = NumberField(x^2 + 20)
        sage: L = range(100)
        sage: bdd_pr_ideals = bdd_norm_pr_gens_iq(K, L)
        sage: bdd_pr_ideals[100]
        Traceback (most recent call last):
        ...
        KeyError: 100

    """
    gens = dict()
    for n in norm_list:
        gens[n] = K.elements_of_norm(n)
    return gens


def bdd_height_iq(K, height_bound, GRH=False):
    r"""
    Compute all elements in the imaginary quadratic field ``K`` which have
    relative multiplicative height at most ``height_bound``.

    The function will only be called with ``K`` an imaginary quadratic field.
    
    If called with ``K`` not an imaginary quadratic, the function will likely
    yield incorrect output.

    ALGORITHM:

    This is an implementation of Algorithm 5 in John R. Doyle and David Krumm,
    Computing algebraic numbers of bounded height, :arxiv:`1111.4963` (2013).

    INPUT:

    - ``K`` - an imaginary quadratic number field
    - ``height_bound`` - a real number

    OUTPUT:

    - a list of elements of ``K``

    EXAMPLES::

        sage: from sage.rings.number_field.bdd_height import bdd_height_iq
        sage: K.<a> = NumberField(x^2 + 191)
        sage: for t in bdd_height_iq(K,8):
        ....:     print exp(2*t.global_height())
        ....:     
        1.00000000000000
        1.00000000000000
        1.00000000000000
        4.00000000000000
        4.00000000000000
        4.00000000000000
        4.00000000000000
        8.00000000000000
        8.00000000000000
        8.00000000000000
        8.00000000000000
        8.00000000000000
        8.00000000000000
        8.00000000000000
        8.00000000000000

    There are 175 elements of height at most 10 in Q(sqrt(-3))::

        sage: from sage.rings.number_field.bdd_height import bdd_height_iq
        sage: K.<a> = NumberField(x^2 + 3)
        sage: len(list(bdd_height_iq(K,10)))
        175

    The only elements of multiplicative height 1 in a number field are 0 and
    the roots of unity::

        sage: from sage.rings.number_field.bdd_height import bdd_height_iq
        sage: K.<a> = NumberField(x^2 + x + 1)
        sage: list(bdd_height_iq(K,1))
        [0, -a, -a - 1, -1, a, a + 1, 1]

    A number field has no elements of multiplicative height less than 1::

        sage: from sage.rings.number_field.bdd_height import bdd_height_iq
        sage: K.<a> = NumberField(x^2 + 5)
        sage: list(bdd_height_iq(K,0.9))
        []

    """
    if GRH:
        # Assume GRH throughout
        number_field(False)
    if height_bound < 1:
        return
    yield K(0)
    roots_of_unity = K.roots_of_unity()
    for zeta in roots_of_unity:
        yield zeta

    # Get a complete set of ideal class representatives
    class_group_reps = []
    class_group_rep_norms = []
    for c in K.class_group():
        a = c.ideal()
        class_group_reps.append(a)
        class_group_rep_norms.append(a.norm())
    class_number = len(class_group_reps)

    # Find principal ideals of bounded norm
    possible_norm_set = set([])
    for n in xrange(class_number):
        for m in xrange(1, height_bound + 1):
            possible_norm_set.add(m*class_group_rep_norms[n])
    bdd_ideals = bdd_norm_pr_gens_iq(K, possible_norm_set)

    # Distribute the principal ideals
    generator_lists = []
    for n in xrange(class_number):
        this_ideal = class_group_reps[n]
        this_ideal_norm = class_group_rep_norms[n]
        gens = []
        for i in xrange(1, height_bound + 1):
            for g in bdd_ideals[i*this_ideal_norm]:
                if g in this_ideal:
                    gens.append(g)
        generator_lists.append(gens)

    # Build all the output numbers
    for n in xrange(class_number):
        gens = generator_lists[n]
        s = len(gens)
        for i in xrange(s):
            for j in xrange(i + 1, s):
                if K.ideal(gens[i], gens[j]) == class_group_reps[n]:
                    new_number = gens[i]/gens[j]
                    for zeta in roots_of_unity:
                        yield zeta*new_number
                        yield zeta/new_number
    # Stop assuming GRH
    number_field(True)

def bdd_norm_pr_ideal_gens(K, norm_list):
    r"""
    Compute generators for all principal ideals in a number field ``K`` whose
    norms are in ``norm_list``.

    INPUT:

    - ``K`` - a number field
    - ``norm_list`` - a list of positive integers

    OUTPUT:

    - a dictionary, keyed by norm

    EXAMPLES:

    There is only one principal ideal of norm 1, and it is generated by the
    element 1::

        sage: from sage.rings.number_field.bdd_height import bdd_norm_pr_ideal_gens
        sage: K.<g> = QuadraticField(101)
        sage: bdd_norm_pr_ideal_gens(K, [1])
        {1: [1]}

    ::

        sage: from sage.rings.number_field.bdd_height import bdd_norm_pr_ideal_gens
        sage: K.<g> = QuadraticField(123)
        sage: bdd_norm_pr_ideal_gens(K, range(5))
        {0: [], 1: [1], 2: [g - 11], 3: [], 4: [2]}
    
    ::

        sage: from sage.rings.number_field.bdd_height import bdd_norm_pr_ideal_gens
        sage: K.<g> = NumberField(x^5 - x + 19)
        sage: b = bdd_norm_pr_ideal_gens(K, range(100))
        sage: key = ZZ(91)
        sage: b[key]
        [-g^4 - 6*g^3 - 14*g^2 - 22*g - 18,
         -2*g^4 - 6*g^3 + 13*g^2 + 6*g - 47,
         -g^3 + g^2 - g + 7,
         -g^4 - 3*g^3 - 5*g^2 - 6*g + 1]

    """

    negative_norm_units = K.elements_of_norm(-1)
    gens = dict()
    if len(negative_norm_units) == 0:
        for n in norm_list:
            gens[n] = K.elements_of_norm(n) + K.elements_of_norm(-n)
    else:
        for n in norm_list:
            gens[n] = K.elements_of_norm(n)
    return gens

def integer_points_in_polytope(matrix, interval_radius):
    r"""
    Return the set of integer points in the polytope obtained by acting on a
    cube by a linear transformation.
    
    Given an r-by-r matrix ``matrix`` and a real number ``interval_radius``,
    this function finds all integer lattice points in the polytope obtained by
    transforming the cube [-interval_radius,interval_radius]^r via the linear
    map induced by ``matrix``.

    INPUT:

    - ``matrix`` - a square matrix of real numbers
    - ``interval_radius`` - a real number

    OUTPUT:

    - a list of tuples of integers

    EXAMPLES:

    Stretch the interval [-1,1] by a factor of 2 and find the integers in the
    resulting interval::

        sage: from sage.rings.number_field.bdd_height import integer_points_in_polytope
        sage: m = matrix([2])
        sage: r = 1
        sage: integer_points_in_polytope(m,r)
        [(-2), (-1), (0), (1), (2)]

    Integer points inside a parallelogram::

        sage: from sage.rings.number_field.bdd_height import integer_points_in_polytope
        sage: m = matrix([[1, 2],[3, 4]])
        sage: r = RealField()(1.3)
        sage: integer_points_in_polytope(m,r)
        [(-3, -7), (-2, -5), (-2, -4), (-1, -3), (-1, -2), (-1, -1), (0, -1), (0, 0), (0, 1), (1, 1), (1, 2), (1, 3), (2, 4), (2, 5), (3, 7)]

    Integer points inside a parallelepiped::

        sage: from sage.rings.number_field.bdd_height import integer_points_in_polytope
        sage: m = matrix([[1.2,3.7,0.002],[-9.65,-0.00123,5],[1.2,20.7,-107.93]])
        sage: r = 2.7
        sage: L = integer_points_in_polytope(m,r)    # long time (40 s)
        sage: len(L)                                  # long time
        622837

    If ``interval_radius`` is 0, the output should include only the zero tuple::

        sage: from sage.rings.number_field.bdd_height import integer_points_in_polytope
        sage: m = matrix([[1,2,3,7],[4,5,6,2],[7,8,9,3],[0,3,4,5]])
        sage: integer_points_in_polytope(m,0)
        [(0, 0, 0, 0)]

    """

    T = matrix; d = interval_radius; r = T.nrows()

    # Find the vertices of the given box
    box_vertices = [vector(x) for x in itertools.product([-d, d], repeat=r)]

    # Transform the vertices
    T_trans = T.transpose()
    transformed_vertices = [v*T_trans for v in box_vertices]

    # Create polyhedron from transformed vertices and find integer points inside
    return list(Polyhedron(transformed_vertices, base_ring=QQ).integral_points())


def bdd_height(K, height_bound, precision=53, lll=False, GRH=False):
    r"""
    Computes all elements in the number number field ``K`` which have relative
    multiplicative height at most ``height_bound``.

    The algorithm requires arithmetic with floating point numbers;
    ``precision`` gives the user the option to set the precision for such
    computations.
    
    It might be helpful to work with an LLL-reduced system of fundamental 
    units, so the user has the option to perform an LLL reduction for the 
    fundamental units by setting ``lll`` to True.
    
    The algorithm may be made to run calculations assuming the Generalized
    Riemann Hypothesis by setting ``GRH`` equal to True.
    
    The function will only be called for number fields ``K`` with positive unit
    rank. An error will occur if ``K`` is QQ or an imaginary quadratic field.

    ALGORITHM:

    This is an implementation of the main algorithm (Algorithm 3) in John R.
    Doyle and David Krumm, Computing algebraic numbers of bounded
    height, :arxiv:`1111.4963` (2013).

    INPUT:

    - ``height_bound`` - real number
    - ``precision`` - (default: 53) positive integer
    - ``lll`` - (default: False) boolean value
    - ``GRH`` - (default: False) boolean value

    OUTPUT:

    - list of elements of ``K``

    EXAMPLES:

    There are no elements of negative height::

        sage: from sage.rings.number_field.bdd_height import bdd_height
        sage: K.<g> = NumberField(x^5 - x + 7)
        sage: list(bdd_height(K,-3))
        []

    The only nonzero elements of height 1 are the roots of unity::

        sage: from sage.rings.number_field.bdd_height import bdd_height
        sage: K.<g> = QuadraticField(3)
        sage: list(bdd_height(K,1))
        [0, -1, 1]

    ::

        sage: from sage.rings.number_field.bdd_height import bdd_height
        sage: K.<g> = QuadraticField(36865)
        sage: len(list(bdd_height(K,1000))) # long time (50 s)
        54679

    ::

        sage: from sage.rings.number_field.bdd_height import bdd_height
        sage: K.<g> = NumberField(x^6 + 2)
        sage: len(list(bdd_height(K,500))) # long time (290 s)
        124911

    ::

        sage: from sage.rings.number_field.bdd_height import bdd_height
        sage: K.<g> = NumberField(x^3 - 197*x + 39)
        sage: len(list(bdd_height(K, 500))) # long time (25 s)
        3207

    """

    if GRH:
        # Assume GRH throughout
        number_field(False)
    B = height_bound
    r1, r2 = K.signature(); r = r1 + r2 -1
    if B < 1:
        return
    yield K(0)
    roots_of_unity = K.roots_of_unity()
    if B == 1:
        for zeta in roots_of_unity:
            yield zeta
        return
    RF = sage.rings.real_mpfr.RealField(precision)
    embeddings = K.places(prec=precision)
    logB = RF(B).log()

    def log_map(number):
        r"""
        Computes the image of an element of K under the logarithmic map.
        """
        x = number
        x_logs = []
        for i in xrange(r1):
            sigma = embeddings[i]
            x_logs.append(abs(sigma(x)).log())
        for i in xrange(r1, r + 1):
            tau = embeddings[i]
            x_logs.append(2*abs(tau(x)).log())
        return vector(x_logs)

    def log_height_for_generators(n, i, j):
        r"""
        Computes the logarithmic height of elements of the form g_i/g_j.
        """
        gen_logs = generator_logs[n]
        Log_gi = gen_logs[i]; Log_gj = gen_logs[j]
        arch_sum = sum([max(Log_gi[k], Log_gj[k]) for k in range(r + 1)])
        return (arch_sum - class_group_rep_norm_logs[n])

    def packet_height(n, pair, u):
        r"""
        Computes the height of the element of K encoded by a given packet.
        """
        gen_logs = generator_logs[n]
        i = pair[0] ; j = pair[1]
        Log_gi = gen_logs[i]; Log_gj = gen_logs[j]
        Log_u_gi = Log_gi + unit_log_dictionary[u]
        arch_sum = sum([max(Log_u_gi[k], Log_gj[k]) for k in range(r + 1)])
        return (arch_sum - class_group_rep_norm_logs[n])

    class_group_reps = []
    class_group_rep_norms = []
    class_group_rep_norm_logs = []
    for c in K.class_group():
        a = c.ideal()
        a_norm = a.norm()
        class_group_reps.append(a)
        class_group_rep_norms.append(a_norm)
        class_group_rep_norm_logs.append(RF(a_norm).log())
    class_number = len(class_group_reps)

    # Get fundamental units and their images under the log map
    fund_units = UnitGroup(K).fundamental_units()
    fund_unit_logs = [log_map(fund_units[i]) for i in range(r)]
    unit_prec_test = fund_unit_logs
    try:
        [l.change_ring(QQ) for l in unit_prec_test]
    except ValueError:
        raise ValueError('Precision too low.')
    
    # If lll is set to True, find an LLL-reduced system of fundamental units   
    if lll:
        cut_fund_unit_logs = column_matrix(fund_unit_logs).delete_rows([r])
        lll_fund_units = []
        for c in pari(cut_fund_unit_logs).qflll().python().columns():
            new_unit = 1
            for i in xrange(r):
                new_unit *= fund_units[i]**c[i]
            lll_fund_units.append(new_unit)
        fund_units = lll_fund_units
        fund_unit_logs = map(log_map, fund_units)
        unit_prec_test = fund_unit_logs
        try:
            [l.change_ring(QQ) for l in unit_prec_test]
        except ValueError:
            raise ValueError('Precision too low.')

    # Find principal ideals of bounded norm
    possible_norm_set = set([])
    for n in xrange(class_number):
        for m in xrange(1, B + 1):
            possible_norm_set.add(m*class_group_rep_norms[n])
    bdd_ideals = bdd_norm_pr_ideal_gens(K, possible_norm_set)

    # Distribute the principal ideals
    generator_lists = []
    generator_logs = []
    for n in xrange(class_number):
        this_ideal = class_group_reps[n]
        this_ideal_norm = class_group_rep_norms[n]
        gens = []
        gen_logs = []
        for i in xrange(1, B + 1):
            for g in bdd_ideals[i*this_ideal_norm]:
                if g in this_ideal:
                    gens.append(g)
                    gen_logs.append(log_map(g))
        generator_lists.append(gens)
        generator_logs.append(gen_logs)

    # Compute the lists of relevant pairs and corresponding heights
    gen_height_dictionary = dict()
    relevant_pair_lists = [] 
    for n in xrange(class_number):
        relevant_pairs = []
        gens = generator_lists[n]
        s = len(gens)
        for i in xrange(s):
            for j in xrange(i + 1, s):
                if K.ideal(gens[i], gens[j]) == class_group_reps[n]:
                    relevant_pairs.append([i, j])
                    gen_height_dictionary[(n, i, j)] = log_height_for_generators(n, i, j)
        relevant_pair_lists.append(relevant_pairs)

    # Find the bound for units to be computed
    gen_height_list = [gen_height_dictionary[x] for x in gen_height_dictionary.keys()]
    if len(gen_height_list) == 0:
        d = logB
    else:
        d = logB + max(gen_height_list)

    # Create the matrix whose columns are the logs of the fundamental units
    S = column_matrix(fund_unit_logs).delete_rows([r])
    try:
        T = S.inverse()
    except ZeroDivisionError:
        raise ValueError('Precision too low.')

    # Find all integer lattice points in the unit polytope
    U = integer_points_in_polytope(T, ceil(d))

    U0 = []; L0 = []

    # Compute unit heights
    unit_height_dictionary = dict()
    unit_log_dictionary = dict()
    Ucopy = copy.copy(U)

    for u in U:
        u_log = sum([u[j]*fund_unit_logs[j] for j in range(r)])
        unit_log_dictionary[u] = u_log
        u_height = sum([max(u_log[k], 0) for k in range(r + 1)])
        unit_height_dictionary[u] = u_height
        if u_height <= logB:
            U0.append(u)
        if u_height > d:
            Ucopy.remove(u)
    U = Ucopy

    # Sort U by height
    U = sorted(U, key=lambda u: unit_height_dictionary[u])
    U_length = len(U)

    all_unit_tuples  = set(copy.copy(U0))

    # Check candidate heights
    for n in xrange(class_number):
        relevant_pairs = relevant_pair_lists[n]
        for pair in relevant_pairs:
            i = pair[0] ; j = pair[1]
            gen_height = gen_height_dictionary[(n, i, j)]
            u_height_bound = logB + gen_height
            for k in xrange(U_length):
                u = U[k]
                u_height = unit_height_dictionary[u]
                if u_height <= u_height_bound:
                    candidate_height = packet_height(n, pair, u)
                    if candidate_height <= logB:
                        L0.append([n, pair, u])
                        all_unit_tuples.add(u)
                else:
                    break

    # Use previous data to build all necessary units
    units_dictionary = dict()
    for u in all_unit_tuples:
        unit = K(1)
        for j in xrange(r):
            unit *= (fund_units[j])**(u[j])
        units_dictionary[u] = unit

    # Build all the output numbers
    for u in U0:
        unit = units_dictionary[u]
        for zeta in roots_of_unity:
            yield zeta*unit
    for packet in L0:
        n = packet[0] ; pair = packet[1] ; u = packet[2]
        i = pair[0] ; j = pair[1]
        relevant_pairs = relevant_pair_lists[n]
        gens = generator_lists[n]
        unit = units_dictionary[u]
        c = unit*gens[i]/gens[j]
        for zeta in roots_of_unity:
            yield zeta*c
            yield zeta/c

    # Stop assuming GRH
    number_field(True)

