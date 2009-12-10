Linear Programming
==================

Basics
------

What is a Linear Program ?
""""""""""""""""""""""""""

A linear program consists of the following two pieces of information :

    * A linear function, called the objective, which is
      to be maximized or minimized (for example `2 x + y`)
    * Linear constraints on the variables (for example,
      `3 x + y \leq 2` and  `2 x + 3 y \leq 8`)

The solver will then try to find a solution to the system of
constraints such that the objective function is optimized, and
return the values of the variables.



What is a Mixed Integer Linear Program ?
""""""""""""""""""""""""""""""""""""""""""

It is simply a Linear Program such that some variables are forced
to take integer values instead of real values. This difference becomes
very important when one learns that solving a Linear Program can be done
in polynomial time while solving a general Mixed Integer Linear Program
is `NP`-Complete (= there is no polynomial algorithm to solve it,
according to a widely-spread belief that `P\neq NP`)

Why is Linear Programming so useful ?
""""""""""""""""""""""""""""""""""""""

Linear Programming is very useful in many Optimization and
Graph-Theoretical problems because of its wide range of expression.
Most of the time, a natural Linear Program can be easily written
to solve a problem whose solution will be quickly computed thanks
to the wealth of heuristics already contained in Linear Program
Solvers. It is often hard to theoretically find out the execution
time of a Linear Program, though they give very interesting results
in practice.

For more information, you can consult the Wikipedia page dedicated to
Linear Programming : http://en.wikipedia.org/wiki/Linear_programming

How can I solve a linear program using Sage ?
---------------------------------------------

Sage can solve Linear Programs or Mixed Integer Linear Programs through
the class ``MixedIntegerLinearProgram`` defined in ``sage.numerical.mip``. To illustrate how it can be
used, we will try to solve the following problem :

.. MATH::

    \mbox{Maximize : }&2 x_1 + x_2\\
    \mbox{Such that : }&3 x_1 + 4 x_2\leq 2.5\\
    &0.5\leq 1.2 x_1 + 0.5 x_2 \leq 4

First, we need a few informations about ``MIPVariable`` and
how to read the optimal values when the solver has finished its job.

Variables in ``MixedIntegerLinearProgram``
""""""""""""""""""""""""""""""""""""""""""""

A variable linked to an instance of ``MixedIntegerLinearProgram`` behaves exactly as
a dictionary would. It is declared the following way ::

    sage: p=MixedIntegerLinearProgram()
    sage: variable=p.new_variable()

The variable ``variable`` can contain as many keys as you
would like, each of them being formally `unique`. For example
the following constraint (where `P` denotes the pressure,
and `T` the temperature) :

.. MATH::
    2 T_{\mbox{Madrid}} + 3 T_{\mbox{London}} -
    P_{\mbox{Seattle}} + \mbox{flow}_{3,5} +
    8 \mbox{cost}_{(1,3)} + x_3 < 5\\

... can be expressed in Sage (quite naturally, I hope !) this way ::

    sage: p=MixedIntegerLinearProgram()
    sage: temperature=p.new_variable()
    sage: pressure=p.new_variable()
    sage: x=p.new_variable()
    sage: cost=p.new_variable()
    sage: flow=p.new_variable(dim=2)
    sage: p.add_constraint(2*temperature["Madrid"]+3*temperature["London"]-pressure["Seattle"]+flow[3][5]+8*cost[(1,3)]+x[3],max=5)

This example is just meant to show you the different possibilities
offered to you when you use the ``MixedIntegerLinearProgram`` class. You will not need
to declare so many variables in usual applications.

Notice how the variable ``flow`` is defined : you can use any hashable
object as a key for a ``MIPVariable``, but if you think you need
more than one dimension, you need to explicitely say it when
calling ``MixedIntegerLinearProgram.new_variable()``

For the user's convenience, however, there is a default variable
attached to a Linear Program : indeed, the previous implementation
means that each "variable" actually represents a property of
a set of objects (these objects are strings in the case of ``temperature``
or a pair in the case of ``cost``). In some cases, though it is useful to
define an absolute variable which will not be indexed on anything. This can
be done through the following notation ::

    sage: p = MixedIntegerLinearProgram()
    sage: B = p.new_variable()
    sage: p.set_objective( p["first unique variable"] + B[2] + p[-3] )

In this case, two of these "unique" variables are defined through
``p["first unique variable"]`` and ``p[-3]``.

Let us solve this system !
""""""""""""""""""""""""""

Now that we know what are variables,
we are only several lines away from solving our system ::

    sage: # First, we define our MixedIntegerLinearProgram object, setting maximization=True
    sage: p=MixedIntegerLinearProgram( maximization = True )
    sage: x=p.new_variable()
    sage: # Definition of the objective function
    sage: p.set_objective( 2*x[1]+x[2] )
    sage: # Next, the two constraints
    sage: p.add_constraint( 3*x[1]+4*x[2], max=2.5 )
    sage: p.add_constraint( 1.5*x[1]+0.5*x[2], max=4,min=0.5 )
    sage: p.solve() # optional - requires Glpk or COIN-OR/CBC
    1.6666666666666667
    sage: x_sol=p.get_values(x)
    sage: print x_sol
    {1: 0.83333333333333337, 2: 0.0}


The value returned by ``MixedIntegerLinearProgram.solve()`` is the optimal value of
the objective function. To read the values taken by the variables
one needs to call the method ``MixedIntegerLinearProgram.get_values`` which can return
multiple values at the same time if needed (type
``sage: MixedIntegerLinearProgram.get_values?`` for more information on this function)

Some famous examples
---------------------------

Vertex Cover in a graph
""""""""""""""""""""""""

In the Vertex Cover problem, we are given a graph `G` and we want to find
a subset `S` of its vertices of minimal cardinality such that each edge
`e` is incident to at least one vertex of `S`. In order to achieve it, we
define a binary variable `b_v` for each vertex `v`.

.. MATH::

    \mbox{Maximize : }& \sum_{v\in G.vertices()} b_v\\
    \mbox{Such that : }&\forall (u,v)\in G.edges(), b_u + b_v \geq 1\\
    &\forall v,b_v \mbox{ is a binary variable}

In the linear program, the syntax is exactly the same ::

    sage: g=graphs.PetersenGraph()
    sage: p=MixedIntegerLinearProgram(maximization=False)
    sage: b=p.new_variable()
    sage: for (u,v) in g.edges(labels=None):
    ...          p.add_constraint(b[u]+b[v],min=1)
    sage: p.set_binary(b)

And you but have to type ``p.solve()`` to see the result !

Maximum matching in a Graph
""""""""""""""""""""""""""""

In the maximum matching problem, we are given a graph `G`, and we are
looking for a set of edges `M` of maximum cardinality such
that no two edges from `M` are adjacent :

.. MATH::

    \mbox{Maximize : }& \sum_{e\in G.edges()} b_e\\
    \mbox{Such that : }&\forall v\in G.vertices(), \sum_{(v,w)\in G}b_{uv} \leq 1\\
    &\forall e\in G.edges(),b_e \mbox{ is a binary variable}

Here is how this is solved through Sage on a Petersen Graph ::

    sage: g=graphs.PetersenGraph()
    sage: p=MixedIntegerLinearProgram()
    sage: b=p.new_variable(dim=2)
    sage: for u in g.vertices():
    ...    p.add_constraint(sum([b[u][v] for v in g.neighbors(u)]),max=1)
    sage: for (u,v) in g.edges(labels=None):
    ...    p.add_constraint(b[u][v]+b[v][u],min=1,max=1)

And the next step is ``p.solve()`` !


Solvers
-------

Sage solves linear programs by calling specific libraries. Two are available
for the moment :

     * `GLPK <http://www.gnu.org/software/glpk/>`_ : A Linear Program solver
       from `GNU <http://www.gnu.org/>`_
     * `CBC <http://www.coin-or.org/projects/Cbc.xml>`_ : Mixed
       Integer Linear Program solver from
       `COIN-OR <http://www.coin-or.org/>`_

To install them if they are not available on your installation of Sage, type ::

     sage: # To install GLPK
     sage: install_package('glpk') # not tested
     sage: # To install Coin-OR Branch and Cut (CBC)
     sage: install_package('cbc')  # not tested
