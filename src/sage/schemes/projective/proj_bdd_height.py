r"""
Points of bounded height in projective spaces

This module defines sunctions to compute points of bounded height of a given
number field with height less than a specified bound in projective spaces.

Sage functions to list all elements of a given number field with height less
than a specified bound.

AUTHORS:

- Jing Guo (2022): initial version based on David Krumm's code

REFERENCES:

- [Krumm2016]

"""

import itertools

def QQ_points_of_bounded_height(dim, bound):
    r"""
    ...
    """
    if bound < 1:
        return
    PN = ProjectiveSpace(QQ, dim)
    unit_tuples = list(itertools.product([-1, 1], repeat=dim))
    increasing_tuples = itertools.combinations_with_replacement(range(bound + 1), dim + 1)
    for t in increasing_tuples:
        if gcd(t) == 1:
            for p in itertools.permutations(t):
                for u in unit_tuples:
                    new_point = [a*b for a, b in zip(u, p)] + [p[dim]]
                    yield PN(new_point)

def points_of_bounded_height(K, dim, bound, prec=53):
    r"""
    ...
    """
    if bound < 1:
        return

    r1, r2 = K.signature()
    r = r1 + r2 - 1
    K_degree = K.degree() 
    K_embeddings = K.places(prec=prec)
    roots_of_unity = K.roots_of_unity()
    unit_tuples = list(itertools.product(roots_of_unity, repeat=dim))

    PN = ProjectiveSpace(K, dim)

    Reals = RealField(prec)
    logB = Reals(bound).log()

    class_group_ideals = [c.ideal() for c in K.class_group()]
    class_number = len(class_group_ideals)
    class_group_ideal_norms = [a.norm() for a in class_group_ideals]
    norm_bound = bound*max(class_group_ideal_norms)
    fundamental_units = UnitGroup(K).fundamental_units()
    fund_unit_logs = map(log_map, fundamental_units) # TODO log map

    test_matrix = column_matrix(fund_unit_logs)
    try:
        test_matrix.change_ring(QQ)
    except ValueError:
        raise ValueError('prec too low.')

    cut_fund_unit_logs = column_matrix(fund_unit_logs).delete_rows([r])
    lll_fund_units = []
    for c in pari(cut_fund_unit_logs).qflll().python().columns():
        new_unit = 1
        for i in range(r):
            new_unit *= fundamental_units[i]**c[i]
        lll_fund_units.append(new_unit)
    fundamental_units = lll_fund_units
    fund_unit_logs = map(log_map, fundamental_units)

    possible_norm_set = set([])
    for i in range(class_number):
        for k in range(1, bound + 1):
            possible_norm_set.add(k*class_group_ideal_norms[i])
    
    principal_ideal_gens = dict()
    negative_norm_units = K.elements_of_norm(-1)
    if len(negative_norm_units) == 0:
        for m in possible_norm_set:
            principal_ideal_gens[m] = K.elements_of_norm(m) + K.elements_of_norm(-m)
    else:
        for m in possible_norm_set:
            principal_ideal_gens[m] = K.elements_of_norm(m)

    pr_ideal_gen_logs = dict()
    for key in principal_ideal_gens:
        for y in principal_ideal_gens[key]:
            pr_ideal_gen_logs[y] = log_map(y)
    
    fund_parallelotope_vertices = []
    for coefficient_tuple in itertools.product([-1/2, 1/2], repeat=r):
        vertex = sum([coefficient_tuple[i]*fund_unit_logs[i] for i in range(r)])
        fund_parallelotope_vertices.append(vertex)
      
    D_numbers = []
    for v in range(r + 1):
        D_numbers.append(max([vertex[v] for vertex in fund_parallelotope_vertices]))
    
    A_numbers = []
    for v in range(r + 1):
        A_numbers.append(min([pr_ideal_gen_logs[y][v] for y in pr_ideal_gen_logs]))
    
    aux_constant = (1/K_degree)*Reals(norm_bound).log()
    
    L_numbers = []
    for v in range(r1):
        L_numbers.append(aux_constant + D_numbers[v] - A_numbers[v])
    for v in range(r1, r + 1):
        L_numbers.append(2*aux_constant + D_numbers[v] - A_numbers[v])
    L_numbers = vector(L_numbers).change_ring(QQ)

    T = column_matrix(fund_unit_logs).delete_rows([r]).change_ring(QQ)

    M = ((-1)*matrix.identity(r)).insert_row(r, [1 for i in range(r)])
    M = M.transpose().insert_row(0, [0 for i in range(r + 1)]).transpose()
    M = M.change_ring(QQ)
    M.set_column(0, L_numbers)
    vertices = map(vector, Polyhedron(ieqs=list(M)).vertices())

    T_it = T.inverse().transpose()
    unit_polytope = Polyhedron([v*T_it for v in vertices])

    coordinate_space = dict()
    coordinate_space[0] = [[K(0), log_map(0)]]
    int_points = unit_polytope.integral_points()

    units_with_logs = dict()
    for n in int_points:
        new_unit = 1
        for j in range(r):
            new_unit *= fundamental_units[j]**n[j]
        new_unit_log = sum([n[j]*fund_unit_logs[j] for j in range(r)])
        units_with_logs[n] = [new_unit, new_unit_log]

    for norm in principal_ideal_gens:
        coordinate_list = []
        for y in principal_ideal_gens[norm]:
            for n in int_points:
                unit, unit_log = units_with_logs[n]
                y_log = pr_ideal_gen_logs[y]
                g_log = unit_log + y_log
                bool1 = all(g_log[i] <= aux_constant + D_numbers[i] for i in range(r1))
                bool2 = all(g_log[j] <= 2*aux_constant + D_numbers[j] for j in range(r1, r + 1))
                if bool1 and bool2:
                    g = unit*y
                    coordinate_list.append([g, g_log])
        if len(coordinate_list) > 0:
            coordinate_space[norm] = coordinate_list

    for m in range(class_number):
        a = class_group_ideals[m]
        a_norm = class_group_ideal_norms[m]
        log_a_norm = Reals(a_norm).log()
        a_const = (logB + log_a_norm)/K_degree
        a_coordinates = []

        for k in range(bound + 1):
            norm = k*a_norm
            if coordinate_space.has_key(norm):
                for pair in coordinate_space[norm]:
                    g, g_log = pair
                    if g in a:
                        bool1 = all(g_log[i] <= a_const + D_numbers[i] for i in range(r1))
                        bool2 = all(g_log[j] <= 2*a_const + D_numbers[j] for j in range(r1, r + 1))
                        if bool1 and bool2:
                            a_coordinates.append(pair)    

        t = len(a_coordinates) - 1
        points_in_class_a = set([])
        increasing_tuples = itertools.combinations_with_replacement(range(t + 1), dim + 1)
        log_arch_height_bound = logB + log_a_norm
        for index_tuple in increasing_tuples:
            point_coordinates = [a_coordinates[i][0] for i in index_tuple]
            point_coordinate_logs = [a_coordinates[i][1] for i in index_tuple]
            log_arch_height = sum([max([x[i] for x in point_coordinate_logs]) for i in range(r + 1)])
            if log_arch_height <= log_arch_height_bound and a == K.ideal(point_coordinates):
                for p in itertools.permutations(point_coordinates):
                    for u in unit_tuples:
                        new_point = [i*j for i, j in zip(u, p)] + [p[dim]]
                        yield PN(new_point)
