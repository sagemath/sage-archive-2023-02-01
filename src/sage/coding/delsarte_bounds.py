r"""
Delsarte, a.k.a. Linear Programming (LP), upper bounds

This module provides LP upper bounds for the parameters of codes.
The exact LP solver PPL is used by default, ensuring that no
rounding/overflow problems occur.

AUTHORS:

- Dmitrii V. (Dima) Pasechnik (2012-10): initial implementation. Minor fixes (2015)
"""

#*****************************************************************************
#       Copyright (C) 2012 Dima Pasechnik <dimpase@gmail.com>
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 2 of the License, or
# (at your option) any later version.
#                  http://www.gnu.org/licenses/
#*****************************************************************************

def Krawtchouk(n,q,l,x,check=True):
    """
    Compute ``K^{n,q}_l(x)``, the Krawtchouk (a.k.a. Kravchuk) polynomial.

    See :wikipedia:`Kravchuk_polynomials`; It is defined by the generating function
    `(1+(q-1)z)^{n-x}(1-z)^x=\sum_{l} K^{n,q}_l(x)z^l` and is equal to

    .. math::

        K^{n,q}_l(x)=\sum_{j=0}^l (-1)^j(q-1)^{(l-j)}{x \choose j}{n-x \choose l-j},

    INPUT:

    - ``n, q, x`` -- arbitrary numbers

    - ``l`` -- a nonnegative integer

    - ``check`` -- check the input for correctness. ``True`` by default. Otherwise, pass it
      as it is. Use ``check=False`` at your own risk.

    EXAMPLES::

        sage: Krawtchouk(24,2,5,4)
        2224
        sage: Krawtchouk(12300,4,5,6)
        567785569973042442072

    TESTS:

    check that the bug reported on :trac:`19561` is fixed::

        sage: Krawtchouk(3,2,3,3)
        -1
        sage: Krawtchouk(int(3),int(2),int(3),int(3))
        -1
        sage: Krawtchouk(int(3),int(2),int(3),int(3),check=False)
        -5
        sage: Kravchuk(24,2,5,4)
        2224

    other unusual inputs ::

        sage: Krawtchouk(sqrt(5),1-I*sqrt(3),3,55.3).n()
        211295.892797... + 1186.42763...*I
        sage: Krawtchouk(-5/2,7*I,3,-1/10)
        480053/250*I - 357231/400
        sage: Krawtchouk(1,1,-1,1)
        Traceback (most recent call last):
        ...
        ValueError: l must be a nonnegative integer
        sage: Krawtchouk(1,1,3/2,1)
        Traceback (most recent call last):
        ...
        TypeError: no conversion of this rational to integer
    """
    from sage.arith.all import binomial
    from sage.arith.srange import srange
    # Use the expression in equation (55) of MacWilliams & Sloane, pg 151
    # We write jth term = some_factor * (j-1)th term
    if check:
        from sage.rings.integer_ring import ZZ
        l0 = ZZ(l)
        if l0 != l or l0<0:
            raise ValueError('l must be a nonnegative integer')
        l = l0
    kraw = jth_term = (q-1)**l * binomial(n, l) # j=0
    for j in srange(1,l+1):
        jth_term *= -q*(l-j+1)*(x-j+1)/((q-1)*j*(n-j+1))
        kraw += jth_term
    return kraw

Kravchuk = Krawtchouk

def _delsarte_LP_building(n, d, d_star, q, isinteger,  solver, maxc = 0):
    """
    LP builder - common for the two functions; not exported.

    EXAMPLES::

        sage: from sage.coding.delsarte_bounds import _delsarte_LP_building
        sage: _,p=_delsarte_LP_building(7, 3, 0, 2, False, "PPL")
        sage: p.show()
        Maximization:
          x_0 + x_1 + x_2 + x_3 + x_4 + x_5 + x_6 + x_7
        Constraints:
          constraint_0: 1 <= x_0 <= 1
          constraint_1: 0 <= x_1 <= 0
          constraint_2: 0 <= x_2 <= 0
          constraint_3: -7 x_0 - 5 x_1 - 3 x_2 - x_3 + x_4 + 3 x_5 + 5 x_6 + 7 x_7 <= 0
          constraint_4: -7 x_0 - 5 x_1 - 3 x_2 - x_3 + x_4 + 3 x_5 + 5 x_6 + 7 x_7 <= 0
          ...
          constraint_16: - x_0 + x_1 - x_2 + x_3 - x_4 + x_5 - x_6 + x_7 <= 0
        Variables:
          x_0 is a continuous variable (min=0, max=+oo)
          ...
          x_7 is a continuous variable (min=0, max=+oo)

    """
    from sage.numerical.mip import MixedIntegerLinearProgram

    p = MixedIntegerLinearProgram(maximization=True, solver=solver)
    A = p.new_variable(integer=isinteger, nonnegative=True)
    p.set_objective(sum([A[r] for r in xrange(n+1)]))
    p.add_constraint(A[0]==1)
    for i in xrange(1,d):
        p.add_constraint(A[i]==0)
    for j in xrange(1,n+1):
        rhs = sum([Krawtchouk(n,q,j,r,check=False)*A[r] for r in xrange(n+1)])
        p.add_constraint(0*A[0] <= rhs)
        if j >= d_star:
          p.add_constraint(0*A[0] <= rhs)
        else: # rhs is proportional to j-th weight of the dual code
          p.add_constraint(0*A[0] == rhs)

    if maxc > 0:
        p.add_constraint(sum([A[r] for r in xrange(n+1)]), max=maxc)
    return A, p

def delsarte_bound_hamming_space(n, d, q, return_data=False, solver="PPL", isinteger=False):
    """
    Find the Delsarte bound [1]_ on codes in Hamming space ``H_q^n`` of minimal distance ``d``


    INPUT:

    - ``n`` -- the code length

    - ``d`` -- the (lower bound on) minimal distance of the code

    - ``q`` -- the size of the alphabet

    - ``return_data`` -- if ``True``, return a triple ``(W,LP,bound)``, where ``W`` is
        a weights vector,  and ``LP`` the Delsarte upper bound LP; both of them are Sage LP
        data.  ``W`` need not be a weight distribution of a code.

    - ``solver`` -- the LP/ILP solver to be used. Defaults to ``PPL``. It is arbitrary
        precision, thus there will be no rounding errors. With other solvers
        (see :class:`MixedIntegerLinearProgram` for the list), you are on your own!

    - ``isinteger`` -- if ``True``, uses an integer programming solver (ILP), rather
        that an LP solver. Can be very slow if set to ``True``.

    EXAMPLES:

    The bound on the size of the `F_2`-codes of length 11 and minimal distance 6::

       sage: delsarte_bound_hamming_space(11, 6, 2)
       12
       sage: a, p, val = delsarte_bound_hamming_space(11, 6, 2, return_data=True)
       sage: [j for i,j in p.get_values(a).iteritems()]
       [1, 0, 0, 0, 0, 0, 11, 0, 0, 0, 0, 0]

    The bound on the size of the `F_2`-codes of length 24 and minimal distance
    8, i.e. parameters of the extened binary Golay code::

       sage: a,p,x=delsarte_bound_hamming_space(24,8,2,return_data=True)
       sage: x
       4096
       sage: [j for i,j in p.get_values(a).iteritems()]
       [1, 0, 0, 0, 0, 0, 0, 0, 759, 0, 0, 0, 2576, 0, 0, 0, 759, 0, 0, 0, 0, 0, 0, 0, 1]

    The bound on the size of `F_4`-codes of length 11 and minimal distance 3::

       sage: delsarte_bound_hamming_space(11,3,4)
       327680/3

    An improvement of a known upper bound (150) from http://www.win.tue.nl/~aeb/codes/binary-1.html ::

       sage: a,p,x= delsarte_bound_hamming_space(23,10,2,return_data=True,isinteger=True); x # long time
       148
       sage: [j for i,j in p.get_values(a).iteritems()]                                      # long time
       [1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 95, 0, 2, 0, 36, 0, 14, 0, 0, 0, 0, 0, 0, 0]

    Note that a usual LP, without integer variables, won't do the trick ::

       sage: delsarte_bound_hamming_space(23,10,2).n(20)
       151.86

    Such an input is invalid::

       sage: delsarte_bound_hamming_space(11,3,-4)
       Solver exception: PPL : There is no feasible solution
       False

    REFERENCES:

    .. [1] \P. Delsarte, An algebraic approach to the association schemes of coding theory,
           Philips Res. Rep., Suppl., vol. 10, 1973.



    """
    from sage.numerical.mip import MIPSolverException
    A, p = _delsarte_LP_building(n, d, 0, q, isinteger, solver)
    try:
        bd=p.solve()
    except MIPSolverException as exc:
        print "Solver exception:", exc
        if return_data:
            return A,p,False
        return False

    if return_data:
        return A,p,bd
    else:
        return bd

def delsarte_bound_additive_hamming_space(n, d, q, d_star=1, q_base=0,
                     return_data=False, solver="PPL", isinteger=False):
   """
   Find a modified Delsarte bound on additive codes in Hamming space ``H_q^n`` of minimal distance ``d``

   Find the Delsarte LP bound on ``F_{q_base}``-dimension of additive codes in
   Hamming space ``H_q^n`` of minimal distance ``d`` with minimal distance of the dual
   code at least ``d_star``.  If ``q_base`` is set to
   non-zero, then  ``q`` is a power of ``q_base``, and the code is, formally, linear over
   ``F_{q_base}``. Otherwise it is assumed that ``q_base==q``.


   INPUT:

   - ``n`` -- the code length

   - ``d`` -- the (lower bound on) minimal distance of the code

   - ``q`` -- the size of the alphabet

   - ``d_star`` -- the (lower bound on) minimal distance of the dual code;
     only makes sense for additive codes.

   - ``q_base`` -- if ``0``, the code is assumed to be linear. Otherwise,
     ``q=q_base^m`` and the code is linear over ``F_{q_base}``.

   - ``return_data`` -- if ``True``, return a triple ``(W,LP,bound)``, where ``W`` is
     a weights vector,  and ``LP`` the Delsarte bound LP; both of them are Sage LP
     data.  ``W`` need not be a weight distribution of a code, or,
     if ``isinteger==False``, even have integer entries.

   - ``solver`` -- the LP/ILP solver to be used. Defaults to ``PPL``. It is arbitrary
     precision, thus there will be no rounding errors. With other solvers
     (see :class:`MixedIntegerLinearProgram` for the list), you are on your own!

   - ``isinteger`` -- if ``True``, uses an integer programming solver (ILP), rather
     that an LP solver. Can be very slow if set to ``True``.

   EXAMPLES:

   The bound on dimension of linear `F_2`-codes of length 11 and minimal distance 6::

       sage: delsarte_bound_additive_hamming_space(11, 6, 2)
       3
       sage: a,p,val=delsarte_bound_additive_hamming_space(11, 6, 2,\
                                      return_data=True)
       sage: [j for i,j in p.get_values(a).iteritems()]
       [1, 0, 0, 0, 0, 0, 5, 2, 0, 0, 0, 0]

   The bound on the dimension of linear `F_4`-codes of length 11 and minimal distance 3::

       sage: delsarte_bound_additive_hamming_space(11,3,4)
       8

   The bound on the `F_2`-dimension of additive `F_4`-codes of length 11 and minimal
   distance 3::

       sage: delsarte_bound_additive_hamming_space(11,3,4,q_base=2)
       16

   Such a d_star is not possible::

       sage: delsarte_bound_additive_hamming_space(11,3,4,d_star=9)
       Solver exception: PPL : There is no feasible solution
       False

   TESTS::

       sage: a,p,x=delsarte_bound_additive_hamming_space(19,15,7,return_data=True,isinteger=True)
       sage: [j for i,j in p.get_values(a).iteritems()]
       [1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 307, 0, 0, 1, 34]
       sage: delsarte_bound_additive_hamming_space(19,15,7,solver='glpk')
       3
       sage: delsarte_bound_additive_hamming_space(19,15,7,isinteger=True,solver='glpk')
       3
   """
   from sage.numerical.mip import MIPSolverException
   if q_base == 0:
      q_base = q

   kk = 0
   while q_base**kk < q:
      kk += 1

   if q_base**kk != q:
      print "Wrong q_base=", q_base, " for q=", q, kk
      return False

   # this implementation assumes that our LP solver to be unable to do a hot
   # restart with an adjusted constraint

   m = kk*n # this is to emulate repeat/until block
   bd = q**n+1

   while q_base**m < bd: # need to solve the LP repeatedly, as this is a new constraint!
                         # we might become infeasible. More precisely, after rounding down
                         # to the closest value of q_base^m, the LP, with the constraint that
                         # the objective function is at most q_base^m,
      A, p = _delsarte_LP_building(n, d, d_star, q, isinteger,  solver, q_base**m)
      try:
        bd=p.solve()
      except MIPSolverException as exc:
        print "Solver exception:", exc
        if return_data:
           return A,p,False
        return False
    # rounding the bound down to the nearest power of q_base, for q=q_base^m
#      bd_r = roundres(log(bd, base=q_base))
      m = -1
      while q_base**(m+1) < bd:
        m += 1
      if q_base**(m+1) == bd:
        m += 1

   if return_data:
      return A, p, m
   else:
      return m
