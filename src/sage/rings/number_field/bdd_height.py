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


def _bdd_norm_pr_gens_iq(K, norm_list):
    r"""
    Computes generators for all principal ideals in an imaginary quadratic field ``K`` whose norms are in ``norm_list``.

    The only keys for the output dictionary are integers n appearing in ``norm_list``.

    The function will only be called with ``K`` an imaginary quadratic field. The function will return a dictionary for other number fields, but it may be incorrect.

    INPUT:

    - ``K`` - an imaginary quadratic number field
    - ``norm_list`` - a list of positive integers

    OUTPUT:

    - a dictionary, keyed by norm

    EXAMPLES:

        In Q(i), there is one principal ideal of norm 4, two principal ideals of norm 5, but no principal ideals of norm 7::

        sage: K.<g> = NumberField(x^2 + 1)
        sage: L = range(10)
        sage: bdd_pr_ideals = _bdd_norm_pr_gens_iq(K, L)
        sage: bdd_pr_ideals[4]
        [2]
        sage: bdd_pr_ideals[5]
        [g + 2, g - 2]
        sage: bdd_pr_ideals[7]
        []

        There are no ideals in the ring of integers with negative norm::

        sage: K.<g> = NumberField(x^2 + 10)
        sage: L = range(-10, -1)
        sage: bdd_pr_ideals = _bdd_norm_pr_gens_iq(K, L)
        sage: bdd_pr_ideals
        {-10: [], -9: [], -8: [], -7: [], -6: [], -5: [], -4: [], -3: [], -2: []}


        Calling a key that is not in the input ``norm_list`` raises a KeyError::

        sage: K.<g> = NumberField(x^2 + 20)
        sage: L = range(100)
        sage: bdd_pr_ideals = _bdd_norm_pr_gens_iq(K, L)
        sage: bdd_pr_ideals[100]
        Traceback (most recent call last):
        ...
        KeyError: 100

    """
    gens = dict()
    for n in norm_list:
        gens[n] = K.elements_of_norm(n)
    return gens


def _bdd_height_iq(K, height_bound):
    r"""
    Computes all elements in the imaginary quadratic number field ``K`` which have relative multiplicative height at most
    ``height_bound``.

    The function will only be called with ``K`` an imaginary quadratic field. The function will return a list for other number fields, but it will likely be incorrect.

    ALGORITHM:

    This is an implementation of Algorithm 5 in John R. Doyle and David Krumm, Computing algebraic numbers of bounded height,
    http://arxiv.org/abs/1111.4963 (2013).

    INPUT:

    - ``K`` - an imaginary quadratic number field
    - ``height_bound`` - a real number

    OUTPUT:

    - a list of elements of ``self``

    EXAMPLES::

        sage: K.<a> = NumberField(x^2 + 191)
        sage: _bdd_height_iq(K, 8)
        [0, -1, 1, -1/2, -2, 1/2, 2, 1/16*a + 1/16, -1/12*a + 1/12, -1/16*a - 1/16, 1/12*a - 1/12, 1/16*a - 1/16, -1/12*a - 1/12, -1/16*a + 1/16, 1/12*a + 1/12]

        ::

        sage: K.<a> = NumberField(x^2 + 3)
        sage: L = _bdd_height_iq(K, 10)
        sage: len(L)
        175

        ::

        sage: K.<a> = NumberField(x^2 + 101)
        sage: L = _bdd_height_iq(K, 1000)
        sage: len(L)
        350239

        The only elements of multiplicative height 1 in a number field are 0 and the roots of unity::

        sage: K.<a> = NumberField(x^2 + x + 1)
        sage: _bdd_height_iq(K, 1)
        [0, -a, -a - 1, -1, a, a + 1, 1]

        A number field has no elements of multiplicative height less than 1::

        sage: K.<a> = NumberField(x^2 + 5)
        sage: _bdd_height_iq(K, 0.9)
        []


    """

    if height_bound < 1:
        return []
    roots_of_unity = K.roots_of_unity()
    bdd_numbers = [K(0)] + roots_of_unity # the list of numbers to be returned as final output

    # Find a complete set of ideal class representatives and store their norms.
    class_group_reps = [] # a complete set of class group representatives
    class_group_rep_norms = [] # norms of these ideals
    for c in K.class_group():
        a = c.ideal()
        class_group_reps.append(a)
        class_group_rep_norms.append(a.norm())
    class_number = len(class_group_reps)

    # Find principal ideals of bounded norm.
    possible_norm_set = set([])
    for n in xrange(class_number):
        for m in xrange(1, height_bound + 1):
            possible_norm_set.add(m*class_group_rep_norms[n])
    bdd_ideals = _bdd_norm_pr_gens_iq(K, possible_norm_set)

    # Sort the principal ideals according to which class group representative they are contained in.
    generator_lists = [] # generators for principal ideals of bounded norm
    for n in xrange(class_number):
        this_ideal = class_group_reps[n]
        this_ideal_norm = class_group_rep_norms[n]
        gens = [] #list of generators for principal ideals inside this ideal
        for i in xrange(1, height_bound + 1):
            for g in bdd_ideals[i*this_ideal_norm]:
                if g in this_ideal:
                    gens.append(g)
        generator_lists.append(gens)

    # Build all the numbers of bounded height and include them in the output list.
    for n in xrange(class_number):
        gens = generator_lists[n]
        s = len(gens)
        for i in xrange(s):
            for j in xrange(i + 1, s):
                if K.ideal(gens[i], gens[j]) == class_group_reps[n]:
                    new_number = gens[i]/gens[j]
                    for zeta in roots_of_unity:
                        bdd_numbers += [zeta*new_number, zeta/new_number]

    return bdd_numbers


def _bdd_norm_pr_ideal_gens(K, norm_list):
    r"""
    Computes generators for all principal ideals in a number field ``K`` whose norms are in ``norm_list``.

    INPUT:

    - ``K`` - a number field
    - ``norm_list`` - a list of positive integers

    OUTPUT:

    - a dictionary, keyed by norm

    EXAMPLES::

    There is only one principal ideal of norm 1, and it is generated by the element 1.

    sage: K.<g> = QuadraticField(101)
    sage: _bdd_norm_pr_ideal_gens(K, [1])
    {1: [1]}

    ::

    sage: K.<g> = QuadraticField(123)
    sage: _bdd_norm_pr_ideal_gens(K, range(5))
    {0: [], 1: [1], 2: [g - 11], 3: [], 4: [2]}

    ::

    sage: K.<g> = NumberField(x^5 - x + 19)
    sage: b = _bdd_norm_pr_ideal_gens(K, range(100))
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

def _integer_points_in_polytope(matrix, interval_radius):
    r"""
    Given an r-by-r matrix ``matrix`` and a real number ``interval_radius``, this function finds all integer lattice points in the
    polytope obtained by transforming the region [-interval_radius, interval_radius]^r via the linear map induced by ``matrix``.

    INPUT:

    - ``matrix`` - a square matrix of real numbers
    - ``interval_radius`` - a real number

    OUTPUT:

    - a list of tuples of integers

    EXAMPLES::

        Stretch the interval [-1,1] by a factor of 2 and find the integers in the resulting interval.

        sage: m = matrix([2])
        sage: r = 1
        sage: _integer_points_in_polytope(m,r)
        [(-2,), (-1,), (0,), (1,), (2,)]

    ::

        Integer points inside a parallelogram

        sage: m = matrix([[1, 2],[3, 4]])
        sage: r = RealField()(1.3)
        sage: _integer_points_in_polytope(m,r)
        [(-3, -7), (-2, -5), (-2, -4), (-1, -3), (-1, -2), (-1, -1), (0, -1), (0, 0), (0, 1), (1, 1), (1, 2), (1, 3), (2, 4), (2, 5), (3, 7)]

    ::

        Integer points inside a parallelepiped

        sage: m = matrix([[1.2,3.7,0.002],[-9.65,-0.00123,5],[1.2,20.7,-107.93]])
        sage: r = 2.7
        sage: L = _integer_points_in_polytope(m,r)
        sage: len(L)
        622837

    ::

        If ``interval_radius`` is 0, the output should include only the zero tuple.

        sage: m = matrix([[1,2,3,7],[4,5,6,2],[7,8,9,3],[0,3,4,5]])
        sage: _integer_points_in_polytope(m,0)
        [(0, 0, 0, 0)]

    """

    T = matrix; d = interval_radius
    r = T.nrows()

    # Find the vertices of the given box.
    box_vertices = []
    for x in itertools.product([-d, d], repeat=r):
        box_vertices.append(vector(x))

    # Transform the vertices.
    transformed_vertices = []
    T_trans = T.transpose()
    for v in box_vertices:
        vertex = v*T_trans
        transformed_vertices.append(vertex)

    # Create polyhedron from transformed vertices and find integer points inside
    int_points = Polyhedron(transformed_vertices, base_ring=QQ).integral_points()

    # Return these integer vectors as tuples
    int_point_tuples = []
    for vec in int_points:
        int_point_tuples.append(tuple(vec))
    return int_point_tuples


def bdd_height(self, height_bound, precision=100, GRH=False):
    r"""
    Computes all elements in the number number field ``self`` which have relative multiplicative height at most ``height_bound``.

    The algorithm requires arithmetic with floating point numbers; ``precision`` gives the user the option to set the precision for such
    computations. The algorithm may be made to run calculations assuming the Generalized Riemann Hypothesis by setting ``GRH`` equal to True.

    ALGORITHM:

    This is an implementation of the main algorithm (Algorithm 3) in John R. Doyle and David Krumm, Computing algebraic numbers of bounded
    height, http://arxiv.org/abs/1111.4963 (2013).

    INPUT:

    - ``height_bound`` - real number
    - ``precision`` - (default: 100) positive integer
    - ``GRH`` - (default: False) boolean value

    OUTPUT:

    - list of elements of ``self``

    EXAMPLES::

        There are no elements of negative height

        sage: K.<g> = NumberField(x^5 - x + 7)
        sage: bdd_height(K, -3)
        []

    ::

        The only nonzero elements of height 1 are the roots of unity

        sage: K.<g> = QuadraticField(3)
        sage: bdd_height(K, 1)
        [-1, 1, 0]

    ::

        sage: K.<g> = QuadraticField(-91)
        sage: bdd_height(K, 5)
        [0, -1, 1, -1/2, -2, 1/2, 2, 1/10*g - 3/10, -1/10*g - 3/10, -1/10*g + 3/10, 1/10*g + 3/10]

    ::

        sage: K.<g> = QuadraticField(-107)
        sage: L = bdd_height(K, 1000)
        sage: len(L)
        393775

    ::

        sage: K.<g> = QuadraticField(36865)
        sage: L = bdd_height(K, 1000)
        sage: len(L)
        54679

    ::

        sage: K.<g> = NumberField(x^3 - 197*x + 39)
        sage: L = bdd_height(K, 500)
        sage: len(L)
        3207

    ::

        sage: K.<g> = NumberField(x^6 + 2)
        sage: L = bdd_height(K, 500)
        sage: len(L)
        124911

    """

    if GRH:
        number_field(False) # assume GRH throughout
    B = height_bound
    if B < 1:
        return []
    roots_of_unity = self.roots_of_unity()
    if B == 1:
        return roots_of_unity + [0]
    (r1, r2) = self.signature()
    r = r1 + r2 -1
    if self.degree() == 2 and r == 0:
        return _bdd_height_iq(self, B)

    RF = sage.rings.real_mpfr.RealField(precision)
    embeddings = self.places(prec=precision)
    logB = RF(B).log()

    def _log_map(number):
        r"""
        Computes the image of an element of self under the logarithmic map.
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

    def _log_height_for_generators(n, i, j):
        r"""
        Computes the logarithmic height of elements of the form g_i/g_j.
        """
        gen_logs = generator_logs[n]
        Log_gi = gen_logs[i]; Log_gj = gen_logs[j]
        arch_sum = 0
        for k in xrange(r + 1):
            arch_sum += max(Log_gi[k], Log_gj[k])
        return (arch_sum - class_group_rep_norm_logs[n])

    def _packet_height(n, pair, u):
        r"""
        Computes the height of the element of self encoded by a given packet.
        """
        gen_logs = generator_logs[n]
        i = pair[0] ; j = pair[1]
        Log_gi = gen_logs[i]; Log_gj = gen_logs[j]
        Log_u_gi = Log_gi + unit_log_dictionary[u]
        arch_sum = 0
        for k in xrange(r + 1):
            arch_sum += max(Log_u_gi[k], Log_gj[k])
        return (arch_sum - class_group_rep_norm_logs[n])

    # Find a complete set of ideal class representatives and store the logs of their norms.
    class_group_reps = [] # a complete set of class group representatives
    class_group_rep_norms = [] # norms of these ideals
    class_group_rep_norm_logs = [] # logs of the norms of these ideals
    for c in self.class_group():
        a = c.ideal()
        a_norm = a.norm()
        class_group_reps.append(a)
        class_group_rep_norms.append(a_norm)
        class_group_rep_norm_logs.append(RF(a_norm).log())
    class_number = len(class_group_reps)

    # Find fundamental units and store their images under the log map
    fund_units = UnitGroup(self).fundamental_units()
    fund_unit_logs = [] # images of the fundamental units under the log map
    for i in xrange(r):
        fund_unit_logs.append(_log_map(fund_units[i]))

    # Find principal ideals of bounded norm.
    possible_norm_set = set([])
    for n in xrange(class_number):
        for m in xrange(1, B + 1):
            possible_norm_set.add(m*class_group_rep_norms[n])
    bdd_ideals = _bdd_norm_pr_ideal_gens(self, possible_norm_set)

    # Sort the principal ideals according to which class group representative they are contained in.
    generator_lists = [] # generators for principal ideals of bounded norm
    generator_logs = []
    for n in xrange(class_number):
        this_ideal = class_group_reps[n]
        this_ideal_norm = class_group_rep_norms[n]
        gens = [] # list of generators for principal ideals inside this ideal
        gen_logs = [] # images of the elements in gens under the log map
        for i in xrange(1, B + 1):
            for g in bdd_ideals[i*this_ideal_norm]:
                if g in this_ideal:
                    gens.append(g)
                    gen_logs.append(_log_map(g))
        generator_lists.append(gens)
        generator_logs.append(gen_logs)

    # Compute the lists of relevant pairs and the corresponding heights.
    gen_height_dictionary = dict() # to store heights of g_i/g_j
    relevant_pair_lists = [] # for each ideal class rep A, indices for pairs of generators of A
    for n in xrange(class_number):
        relevant_pairs = [] #list of relevant pairs for this ideal
        gens = generator_lists[n]
        s = len(gens)
        for i in xrange(s):
            for j in xrange(i + 1, s):
                if self.ideal(gens[i], gens[j]) == class_group_reps[n]:
                    relevant_pairs.append([i, j])
                    gen_height_dictionary[(n, i, j)] = _log_height_for_generators(n, i, j)
        relevant_pair_lists.append(relevant_pairs)

    # Find the bound for units to be computed.
    h_max = -1
    for triple in gen_height_dictionary.keys():
        triple_height = gen_height_dictionary[triple]
        if triple_height > h_max:
            h_max = triple_height
    d = logB + h_max

    # Create the matrix S whose columns are the logs of the fundamental units.
    S_columns = []
    for k in xrange(r):
        log_epsilon_k = list(fund_unit_logs[k])
        trash = log_epsilon_k.pop()
        S_columns.append(log_epsilon_k)
    try:
        S = column_matrix(S_columns).change_ring(QQ)
    except ValueError:
        raise ValueError('Please increase the precision.')
    T = S.inverse()

    # Find all integer lattice points in the unit polytope
    U = _integer_points_in_polytope(T, ceil(d))

    L = [self(0)] ; U0 = []; L0 = []

    # Compute unit heights.
    unit_height_dictionary = dict()
    unit_log_dictionary = dict()
    Ucopy = copy.copy(U)
    for u in U:
        # Find the image under the log map of the unit corresponding to u.
        u_log = 0
        for j in xrange(r):
            u_log += u[j]*fund_unit_logs[j]
        unit_log_dictionary[u] = u_log

        # Compute the height of the unit
        u_height = 0
        for k in xrange(r + 1):
            u_height += max(u_log[k], 0)
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

    # Check candidate heights.
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
                    candidate_height = _packet_height(n, pair, u)
                    if candidate_height <= logB:
                        L0.append([n, pair, u])
                        all_unit_tuples.add(u)
                else:
                    break

    # Use previous data to build all units needed to construct the output list.
    units_dictionary = dict()
    for u in all_unit_tuples:
        unit = self(1)
        for j in xrange(r):
            unit *= (fund_units[j])**(u[j])
        units_dictionary[u] = unit

    # Use the computed units and packets to build the final output list.
    for u in U0:
        unit = units_dictionary[u]
        for zeta in roots_of_unity:
            L.append(zeta*unit)
    for packet in L0:
        n = packet[0] ; pair = packet[1] ; u = packet[2]
        i = pair[0] ; j = pair[1]
        relevant_pairs = relevant_pair_lists[n]
        gens = generator_lists[n]
        unit = units_dictionary[u]
        c = unit*gens[i]/gens[j]
        for zeta in roots_of_unity:
            L += [zeta*c, zeta/c]
    number_field(True) # Stop assuming GRH
    return L





