r"""
Symbolic Equations and Inequalities.

Sage can solve symbolic equations and express inequalities. For
example, we derive the quadratic formula as follows::

    sage: a,b,c = var('a,b,c')
    sage: qe = (a*x^2 + b*x + c == 0)
    sage: qe
    a*x^2 + b*x + c == 0
    sage: print solve(qe, x)
    [
    x == -1/2*(b + sqrt(-4*a*c + b^2))/a,
    x == -1/2*(b - sqrt(-4*a*c + b^2))/a
    ]


The operator, left hand side, and right hand side
--------------------------------------------------

Operators::

    sage: eqn = x^3 + 2/3 >= x - pi
    sage: eqn.operator()
    <built-in function ge>
    sage: (x^3 + 2/3 < x - pi).operator()
    <built-in function lt>
    sage: (x^3 + 2/3 == x - pi).operator()
    <built-in function eq>

Left hand side::

    sage: eqn = x^3 + 2/3 >= x - pi
    sage: eqn.lhs()
    x^3 + 2/3
    sage: eqn.left()
    x^3 + 2/3
    sage: eqn.left_hand_side()
    x^3 + 2/3

Right hand side::

    sage: (x + sqrt(2) >= sqrt(3) + 5/2).right()
    sqrt(3) + 5/2
    sage: (x + sqrt(2) >= sqrt(3) + 5/2).rhs()
    sqrt(3) + 5/2
    sage: (x + sqrt(2) >= sqrt(3) + 5/2).right_hand_side()
    sqrt(3) + 5/2


Arithmetic
----------
Add two symbolic equations::

    sage: var('a,b')
    (a, b)
    sage: m = 144 == -10 * a + b
    sage: n = 136 == 10 * a + b
    sage: m + n
    280 == 2*b
    sage: int(-144) + m
    0 == -10*a + b - 144

Subtract two symbolic equations::

    sage: var('a,b')
    (a, b)
    sage: m = 144 == 20 * a + b
    sage: n = 136 == 10 * a + b
    sage: m - n
    8 == 10*a
    sage: int(144) - m
    0 == -20*a - b + 144

Multiply two symbolic equations::

    sage: x = var('x')
    sage: m = x == 5*x + 1
    sage: n = sin(x) == sin(x+2*pi)
    sage: m * n
    x*sin(x) == (5*x + 1)*sin(2*pi + x)
    sage: m = 2*x == 3*x^2 - 5
    sage: int(-1) * m
    -2*x == -3*x^2 + 5

Divide two symbolic equations::

    sage: x = var('x')
    sage: m = x == 5*x + 1
    sage: n = sin(x) == sin(x+2*pi)
    sage: m/n
    x/sin(x) == (5*x + 1)/sin(2*pi + x)
    sage: m = x != 5*x + 1
    sage: n = sin(x) != sin(x+2*pi)
    sage: m/n
    x/sin(x) != (5*x + 1)/sin(2*pi + x)

Substitution
------------

Substitution into relations::

    sage: x, a = var('x, a')
    sage: eq = (x^3 + a == sin(x/a)); eq
    x^3 + a == sin(x/a)
    sage: eq.substitute(x=5*x)
    125*x^3 + a == sin(5*x/a)
    sage: eq.substitute(a=1)
    x^3 + 1 == sin(x)
    sage: eq.substitute(a=x)
    x^3 + x == sin(1)
    sage: eq.substitute(a=x, x=1)
    x + 1 == sin(1/x)
    sage: eq.substitute({a:x, x:1})
    x + 1 == sin(1/x)

Solving
-------

We can solve equations::

    sage: x = var('x')
    sage: S = solve(x^3 - 1 == 0, x)
    sage: S
    [x == 1/2*I*sqrt(3) - 1/2, x == -1/2*I*sqrt(3) - 1/2, x == 1]
    sage: S[0]
    x == 1/2*I*sqrt(3) - 1/2
    sage: S[0].right()
    1/2*I*sqrt(3) - 1/2
    sage: S = solve(x^3 - 1 == 0, x, solution_dict=True)
    sage: S
    [{x: 1/2*I*sqrt(3) - 1/2}, {x: -1/2*I*sqrt(3) - 1/2}, {x: 1}]
    sage: z = 5
    sage: solve(z^2 == sqrt(3),z)
    Traceback (most recent call last):
    ...
    TypeError: 5 is not a valid variable.

We illustrate finding multiplicities of solutions::

    sage: f = (x-1)^5*(x^2+1)
    sage: solve(f == 0, x)
    [x == -I, x == I, x == 1]
    sage: solve(f == 0, x, multiplicities=True)
    ([x == -I, x == I, x == 1], [1, 1, 5])

We can numerically find roots of equations::

    sage: (x == sin(x)).find_root(-2,2)
    0.0
    sage: (x^5 + 3*x + 2 == 0).find_root(-2,2,x)
    -0.63283452024215225
    sage: (cos(x) == sin(x)).find_root(10,20)
    19.634954084936208

We illustrate some valid error conditions::

    sage: (cos(x) != sin(x)).find_root(10,20)
    Traceback (most recent call last):
    ...
    ValueError: Symbolic equation must be an equality.
    sage: (SR(3)==SR(2)).find_root(-1,1)
    Traceback (most recent call last):
    ...
    RuntimeError: no zero in the interval, since constant expression is not 0.

There must be at most one variable::

    sage: x, y = var('x,y')
    sage: (x == y).find_root(-2,2)
    Traceback (most recent call last):
    ...
    NotImplementedError: root finding currently only implemented in 1 dimension.

Assumptions
-----------

Forgetting assumptions::

    sage: var('x,y')
    (x, y)
    sage: forget() #Clear assumptions
    sage: assume(x>0, y < 2)
    sage: assumptions()
    [x > 0, y < 2]
    sage: (y < 2).forget()
    sage: assumptions()
    [x > 0]


Miscellaneous
-------------

Conversion to Maxima::

    sage: x = var('x')
    sage: eq = (x^(3/5) >= pi^2 + e^i)
    sage: eq._maxima_init_()
    '(x)^(3/5) >= ((%pi)^(2))+(exp(0+%i*1))'
    sage: e1 = x^3 + x == sin(2*x)
    sage: z = e1._maxima_()
    sage: z.parent() is sage.calculus.calculus.maxima
    True
    sage: z = e1._maxima_(maxima)
    sage: z.parent() is maxima
    True
    sage: z = maxima(e1)
    sage: z.parent() is maxima
    True

Conversion to Maple::

    sage: x = var('x')
    sage: eq = (x == 2)
    sage: eq._maple_init_()
    'x = 2'

Comparison::

    sage: x = var('x')
    sage: (x>0) == (x>0)
    True
    sage: (x>0) == (x>1)
    False
    sage: (x>0) != (x>1)
    True

Variables appearing in the relation::

    sage: var('x,y,z,w')
    (x, y, z, w)
    sage: f =  (x+y+w) == (x^2 - y^2 - z^3);   f
    w + x + y == x^2 - y^2 - z^3
    sage: f.variables()
    (w, x, y, z)

LaTeX output::

    sage: latex(x^(3/5) >= pi)
    x^{\frac{3}{5}}  \geq  \pi


More Examples
-------------

::

    sage: x,y,a = var('x,y,a')
    sage: f = x^2 + y^2 == 1
    sage: f.solve(x)
    [x == -sqrt(-y^2 + 1), x == sqrt(-y^2 + 1)]

::

    sage: f = x^5 + a
    sage: solve(f==0,x)
    [x == (-a)^(1/5)*e^(2/5*I*pi), x == (-a)^(1/5)*e^(4/5*I*pi), x == (-a)^(1/5)*e^(-4/5*I*pi), x == (-a)^(1/5)*e^(-2/5*I*pi), x == (-a)^(1/5)]

You can also do arithmetic with inequalities, as illustrated
below::

    sage: var('x y')
    (x, y)
    sage: f = x + 3 == y - 2
    sage: f
    x + 3 == y - 2
    sage: g = f - 3; g
    x == y - 5
    sage: h =  x^3 + sqrt(2) == x*y*sin(x)
    sage: h
    sqrt(2) + x^3 == x*y*sin(x)
    sage: h - sqrt(2)
    x^3 == x*y*sin(x) - sqrt(2)
    sage: h + f
    x + sqrt(2) + x^3 + 3 == x*y*sin(x) + y - 2
    sage: f = x + 3 < y - 2
    sage: g = 2 < x+10
    sage: f - g
    x + 1 < -x + y - 12
    sage: f + g
    x + 5 < x + y + 8
    sage: f*(-1)
    -x - 3 < -y + 2

TESTS:

We test serializing symbolic equations::

    sage: eqn = x^3 + 2/3 >= x
    sage: loads(dumps(eqn))
    x^3 + 2/3 >= x
    sage: loads(dumps(eqn)) == eqn
    True

AUTHORS:

- Bobby Moretti: initial version (based on a trick that Robert
  Bradshaw suggested).

- William Stein: second version

- William Stein (2007-07-16): added arithmetic with symbolic equations

"""
import operator
from sage.calculus.calculus import maxima

def test_relation_maxima(relation):
    """
    Return True if this (in)equality is definitely true. Return False
    if it is false or the algorithm for testing (in)equality is
    inconclusive.

    EXAMPLES::

        sage: from sage.symbolic.relation import test_relation_maxima
        sage: k = var('k')
        sage: pol = 1/(k-1) - 1/k -1/k/(k-1);
        sage: test_relation_maxima(pol == 0)
        True
        sage: f = sin(x)^2 + cos(x)^2 - 1
        sage: test_relation_maxima(f == 0)
        True
        sage: test_relation_maxima( x == x )
        True
        sage: test_relation_maxima( x != x )
        False
        sage: test_relation_maxima( x > x )
        False
        sage: test_relation_maxima( x^2 > x )
        False
        sage: test_relation_maxima( x + 2 > x )
        True
        sage: test_relation_maxima( x - 2 > x )
        False
    """
    m = relation._maxima_()

    #Handle some basic cases first
    if repr(m) in ['0=0']:
        return True
    elif repr(m) in ['0#0', '1#1']:
        return False

    try:
        s = m.parent()._eval_line('is (%s)'%m.name())
    except TypeError, msg:
        raise ValueError, "unable to evaluate the predicate '%s'"%repr(relation)

    if s == 'true':
        return True
    elif s == 'unknown':
        return False

    if relation.operator() != operator.eq:
        return False

    difference = relation.lhs() - relation.rhs()
    if repr(difference) == '0':
        return True

    #Try to apply some simplifications to see if left - right == 0
    simp_list = [difference.simplify_log, difference.simplify_rational, difference.simplify_exp,difference.simplify_radical,difference.simplify_trig]
    for f in simp_list:
        try:
            if repr( f() ).strip() == "0":
                return True
                break
        except:
            pass
    return False


def string_to_list_of_solutions(s):
    r"""
    Used internally by the symbolic solve command to convert the output
    of Maxima's solve command to a list of solutions in Sage's symbolic
    package.

    EXAMPLES:

    We derive the (monic) quadratic formula::

        sage: var('x,a,b')
        (x, a, b)
        sage: solve(x^2 + a*x + b == 0, x)
        [x == -1/2*a - 1/2*sqrt(a^2 - 4*b), x == -1/2*a + 1/2*sqrt(a^2 - 4*b)]

    Behind the scenes when the above is evaluated the function
    :func:`string_to_list_of_solutions` is called with input the
    string `s` below::

        sage: s = '[x=-(sqrt(a^2-4*b)+a)/2,x=(sqrt(a^2-4*b)-a)/2]'
        sage: sage.symbolic.relation.string_to_list_of_solutions(s)
         [x == -1/2*a - 1/2*sqrt(a^2 - 4*b), x == -1/2*a + 1/2*sqrt(a^2 - 4*b)]
    """
    from sage.categories.all import Objects
    from sage.structure.sequence import Sequence
    from sage.calculus.calculus import symbolic_expression_from_maxima_string
    v = symbolic_expression_from_maxima_string(s, equals_sub=True)
    return Sequence(v, universe=Objects(), cr_str=True)

###########
# Solving #
###########

def solve(f, *args, **kwds):
    r"""
    Algebraically solve an equation of system of equations for given
    variables.

    INPUT:

    -  ``f`` - equation or system of equations (given by a
       list or tuple)

    -  ``*args`` - variables to solve for.

    -  ``solution_dict = True`` - return a list of
       dictionaries containing the solutions.

    EXAMPLES::

        sage: x, y = var('x, y')
        sage: solve([x+y==6, x-y==4], x, y)
        [[x == 5, y == 1]]
        sage: solve([x^2+y^2 == 1, y^2 == x^3 + x + 1], x, y)
        [[x == -1/2*I*sqrt(3) - 1/2, y == -1/2*sqrt(-I*sqrt(3) + 3)*sqrt(2)],
         [x == -1/2*I*sqrt(3) - 1/2, y == 1/2*sqrt(-I*sqrt(3) + 3)*sqrt(2)],
         [x == 1/2*I*sqrt(3) - 1/2, y == -1/2*sqrt(I*sqrt(3) + 3)*sqrt(2)],
         [x == 1/2*I*sqrt(3) - 1/2, y == 1/2*sqrt(I*sqrt(3) + 3)*sqrt(2)],
         [x == 0, y == -1],
         [x == 0, y == 1]]
        sage: solutions=solve([x^2+y^2 == 1, y^2 == x^3 + x + 1], x, y, solution_dict=True); solutions
        [{x: -1/2*I*sqrt(3) - 1/2, y: -1/2*sqrt(-I*sqrt(3) + 3)*sqrt(2)},
         {x: -1/2*I*sqrt(3) - 1/2, y: 1/2*sqrt(-I*sqrt(3) + 3)*sqrt(2)},
         {x: 1/2*I*sqrt(3) - 1/2, y: -1/2*sqrt(I*sqrt(3) + 3)*sqrt(2)},
         {x: 1/2*I*sqrt(3) - 1/2, y: 1/2*sqrt(I*sqrt(3) + 3)*sqrt(2)},
         {x: 0, y: -1},
         {x: 0, y: 1}]
        sage: for solution in solutions: print solution[x].n(digits=3), ",", solution[y].n(digits=3)
        -0.500 - 0.866*I , -1.27 + 0.341*I
        -0.500 - 0.866*I , 1.27 - 0.341*I
        -0.500 + 0.866*I , -1.27 - 0.341*I
        -0.500 + 0.866*I , 1.27 + 0.341*I
        0.000 , -1.00
        0.000 , 1.00
        sage: z = 5
        sage: solve([8*z + y == 3, -z +7*y == 0],y,z)
        Traceback (most recent call last):
        ...
        TypeError: 5 is not a valid variable.

    If ``True`` appears in the list of equations it is
    ignored, and if ``False`` appears in the list then no
    solutions are returned. E.g., note that the first
    ``3==3`` evaluates to ``True``, not to a
    symbolic equation.

    ::

        sage: solve([3==3, 1.00000000000000*x^3 == 0], x)
        [x == 0]
        sage: solve([1.00000000000000*x^3 == 0], x)
        [x == 0]

    Here, the first equation evaluates to ``False``, so
    there are no solutions::

        sage: solve([1==3, 1.00000000000000*x^3 == 0], x)
        []

    ::

        sage: var('s,j,b,m,g')
        (s, j, b, m, g)
        sage: sys = [ m*(1-s) - b*s*j, b*s*j-g*j ];
        sage: solve(sys,s,j)
        [[s == 1, j == 0], [s == g/b, j == (b - g)*m/(b*g)]]
        sage: solve(sys,(s,j))
        [[s == 1, j == 0], [s == g/b, j == (b - g)*m/(b*g)]]
        sage: solve(sys,[s,j])
        [[s == 1, j == 0], [s == g/b, j == (b - g)*m/(b*g)]]
    """
    try:
        return f.solve(*args,**kwds)
    except AttributeError:
        from sage.symbolic.ring import is_SymbolicVariable

        if is_SymbolicVariable(args[0]):
            variables = args
        else:
            variables = tuple(args[0])

        for v in variables:
            if not is_SymbolicVariable(v):
                raise TypeError, "%s is not a valid variable."%v

        try:
            f = [s for s in f if s is not True]
        except TypeError:
            raise ValueError, "Unable to solve %s for %s"%(f, args)

        if any(s is False for s in f):
            return []

        m = maxima(f)

        try:
            s = m.solve(variables)
        except:
            raise ValueError, "Unable to solve %s for %s"%(f, args)
        sol_list = string_to_list_of_solutions(repr(s))
        if 'solution_dict' in kwds and kwds['solution_dict']==True:
            sol_dict=[dict([[eq.left(),eq.right()] for eq in solution]) for solution in sol_list]
            return sol_dict
        else:
            return sol_list

def solve_mod(eqns, modulus, solution_dict = False):
    r"""
    Return all solutions to an equation or list of equations modulo the
    given integer modulus. Each equation must involve only polynomials
    in 1 or many variables.

    By default the solutions are returned as `n`-tuples, where `n`
    is the number of variables appearing anywhere in the given
    equations. The variables are in alphabetical order.

    INPUT:


    -  ``eqns`` - equation or list of equations

    -  ``modulus`` - an integer

    -  ``solution_dict`` - (default: False) if True,  return a list of
       dictionaries containing the solutions.

    EXAMPLES::

        sage: var('x,y')
        (x, y)
        sage: solve_mod([x^2 + 2 == x, x^2 + y == y^2], 14)
        [(4, 2), (4, 6), (4, 9), (4, 13)]
        sage: solve_mod([x^2 == 1, 4*x  == 11], 15)
        [(14,)]

    Fermat's equation modulo 3 with exponent 5::

        sage: var('x,y,z')
        (x, y, z)
        sage: solve_mod([x^5 + y^5 == z^5], 3)
        [(0, 0, 0), (0, 1, 1), (0, 2, 2), (1, 0, 1), (1, 1, 2), (1, 2, 0), (2, 0, 2), (2, 1, 0), (2, 2, 1)]

    We can solve with respect to a bigger modulus if it consists only of small prime factors::

        sage: solve_mod([5*x + y == 3, 2*x - 3*y == 9], 3*5*7*11*19*23*29, solution_dict = True)
        [{x: 12915279, y: 8610183}]

    We solve an simple equation modulo 2::

        sage: x,y = var('x,y')
        sage: solve_mod([x == y], 2)
        [(0, 0), (1, 1)]

    .. warning::

       The current implementation splits the modulus into prime powers,
       then naively enumerates all possible solutions and finally combines
       the solution using the Chinese Remainder Theorem.
       The interface is good, but the algorithm is horrible if the modulus
       has some larger prime factors! Sage {does} have the ability to do
       something much faster in certain cases at least by using Groebner
       basis, linear algebra techniques, etc. But for a lot of toy problems
       this function as is might be useful. At least it establishes an interface.
    """
    from sage.rings.all import Integer, Integers, PolynomialRing, factor, crt_basis
    from sage.misc.all import cartesian_product_iterator
    from sage.modules.all import vector
    from sage.matrix.all import matrix

    if not isinstance(eqns, (list, tuple)):
        eqns = [eqns]
    modulus = Integer(modulus)
    if modulus < 1:
         raise ValueError, "the modulus must be a positive integer"
    vars = list(set(sum([list(e.variables()) for e in eqns], [])))
    vars.sort(cmp = lambda x,y: cmp(repr(x), repr(y)))
    n = len(vars)

    factors = [p**i for p,i in factor(modulus)]
    crt_basis = vector(Integers(modulus), crt_basis(factors))
    solutions = [solve_mod_enumerate(eqns, p) for p in factors]

    ans = []
    for solution in cartesian_product_iterator(solutions):
        solution_mat = matrix(Integers(modulus), solution)
        ans.append(tuple(c.dot_product(crt_basis) for c in solution_mat.columns()))

    if solution_dict == True:
        sol_dict = [dict(zip(vars, solution)) for solution in ans]
        return sol_dict
    else:
        return ans

def solve_mod_enumerate(eqns, modulus):
    r"""
    Return all solutions to an equation or list of equations modulo the
    given integer modulus. Each equation must involve only polynomials
    in 1 or many variables.

    The solutions are returned as `n`-tuples, where `n`
    is the number of variables appearing anywhere in the given
    equations. The variables are in alphabetical order.

    INPUT:


    -  ``eqns`` - equation or list of equations

    -  ``modulus`` - an integer


    EXAMPLES::

        sage: var('x,y')
        (x, y)
        sage: solve_mod([x^2 + 2 == x, x^2 + y == y^2], 14)
        [(4, 2), (4, 6), (4, 9), (4, 13)]
        sage: solve_mod([x^2 == 1, 4*x  == 11], 15)
        [(14,)]

    Fermat's equation modulo 3 with exponent 5::

        sage: var('x,y,z')
        (x, y, z)
        sage: solve_mod([x^5 + y^5 == z^5], 3)
        [(0, 0, 0), (0, 1, 1), (0, 2, 2), (1, 0, 1), (1, 1, 2), (1, 2, 0), (2, 0, 2), (2, 1, 0), (2, 2, 1)]

    We solve an simple equation modulo 2::

        sage: x,y = var('x,y')
        sage: solve_mod([x == y], 2)
        [(0, 0), (1, 1)]


    .. warning::

       Currently this naively enumerates all possible solutions.  The
       interface is good, but the algorithm is horrible if the modulus
       is at all large! Sage *does* have the ability to do something
       much faster in certain cases at least by using the Chinese
       Remainder Theorem, Groebner basis, linear algebra techniques,
       etc. But for a lot of toy problems this function as is might be
       useful. At the very least, it establishes an interface.
    """
    from sage.rings.all import Integer, Integers, PolynomialRing
    from sage.symbolic.expression import is_Expression
    from sage.misc.all import cartesian_product_iterator

    if not isinstance(eqns, (list, tuple)):
        eqns = [eqns]
    modulus = Integer(modulus)
    if modulus < 1:
         raise ValueError, "the modulus must be a positive integer"
    vars = list(set(sum([list(e.variables()) for e in eqns], [])))
    vars.sort(cmp = lambda x,y: cmp(repr(x), repr(y)))
    n = len(vars)
    R = Integers(modulus)
    S = PolynomialRing(R, len(vars), vars)
    eqns_mod = [S(eq) if is_Expression(eq) else
                S(eq.lhs() - eq.rhs()) for eq in eqns]
    ans = []
    for t in cartesian_product_iterator([R]*len(vars)):
        is_soln = True
        for e in eqns_mod:
            if e(t) != 0:
                is_soln = False
                break
        if is_soln:
            ans.append(t)

    return ans
