r"""
Resolvable Balanced Incomplete Block Design (RBIBD)

This module contains everything related to resolvable Balanced Incomplete Block
Designs. The constructions implemented here can be obtained through the
``designs.<tab>`` object::

    designs.resolvable_balanced_incomplete_block_design(15,3)

For Balanced Incomplete Block Design (BIBD) see the module :mod:`bibd
<sage.combinat.designs.bibd>`. A BIBD
is said to be *resolvable* if its blocks can be partitionned into parallel
classes, i.e.  partitions of its ground set.

The main function of this module is
:func:`resolvable_balanced_incomplete_block_design`, which calls all others.

.. csv-table::
    :class: contentstable
    :widths: 30, 70
    :delim: |

    :func:`resolvable_balanced_incomplete_block_design` | Return a resolvable BIBD of parameters `v,k`.
    :func:`kirkman_triple_system` | Return a Kirkman Triple System on `v` points.
    :func:`v_4_1_rbibd` | Return a `(v,4,1)`-RBIBD

References:

.. [Stinson91] D.R. Stinson,
   A survey of Kirkman triple systems and related designs,
   Volume 92, Issues 1-3, 17 November 1991, Pages 371-393,
   Discrete Mathematics,
   http://dx.doi.org/10.1016/0012-365X(91)90294-C.

.. [RCW71] D. K. Ray-Chaudhuri, R. M. Wilson,
   Solution of Kirkman's schoolgirl problem,
   Volume 19, Pages 187-203,
   Proceedings of Symposia in Pure Mathematics

.. [BJL99] T. Beth, D. Jungnickel, H. Lenz,
   Design Theory 2ed.
   Cambridge University Press
   1999

Functions
---------
"""
from sage.rings.arith import is_prime_power
from sage.combinat.designs.bibd import BalancedIncompleteBlockDesign
from sage.categories.sets_cat import EmptySetError

def resolvable_balanced_incomplete_block_design(v,k,existence=False):
    r"""
    Return a resolvable BIBD of parameters `v,k`.

    A BIBD is said to be *resolvable* if its blocks can be partitionned into
    parallel classes, i.e. partitions of the ground set.

    INPUT:

    - ``v,k`` (integers)

    - ``existence`` (boolean) -- instead of building the design, return:

        - ``True`` -- meaning that Sage knows how to build the design

        - ``Unknown`` -- meaning that Sage does not know how to build the
          design, but that the design may exist (see :mod:`sage.misc.unknown`).

        - ``False`` -- meaning that the design does not exist.

    .. SEEALSO::

        - :meth:`IncidenceStructure.is_resolvable`
        - :func:`~sage.combinat.designs.bibd.balanced_incomplete_block_design`

    EXAMPLE::

        sage: KTS15 = designs.resolvable_balanced_incomplete_block_design(15,3); KTS15
        (15,3,1)-Balanced Incomplete Block Design
        sage: KTS15.is_resolvable()
        True
    """
    # Trivial cases
    if v==1 or k==v:
        return balanced_incomplete_block_design(v,k,existence=existence)

    # Non-existence of resolvable BIBD
    if (v < k or
        k < 2 or
        v%k != 0 or
        (v-1) % (k-1) != 0 or
        (v*(v-1)) % (k*(k-1)) != 0 or
        # From the Handbook of combinatorial designs:
        #
        # With lambda>1 the other exceptions is
        # (15,5,2)
        (k==6 and v == 36) or
        # Fisher's inequality
        (v*(v-1))/(k*(k-1)) < v):
        if existence:
            return False
        raise EmptySetError("There exists no ({},{},{})-RBIBD".format(v,k,1))

    if k==2:
        if existence:
            return True
        classes = [[[(c+i)%(v-1),(c+v-i)%(v-1)] for i in range(1,v//2)]
                   for c in range(v-1)]
        for i,classs in enumerate(classes):
            classs.append([v-1,c])

        B = BalancedIncompleteBlockDesign(v,
                                          classes,
                                          check=True,
                                          copy=False)
        B._classes = classes
        return B
    elif k==3:
        return kirkman_triple_system(v,existence=existence)
    elif k==4:
        return v_4_1_rbibd(v,existence=existence)
    else:
        if existence:
            return Unknown
        raise NotImplementedError("I don't know how to build a ({},{},1)-RBIBD!".format(v,3))

def kirkman_triple_system(v,existence=False):
    r"""
    Return a Kirkman Triple System on `v` points.

    A Kirkman Triple System `KTS(v)` is a resolvable Steiner Triple System. It
    exists if and only if `v\equiv 3\pmod{6}`.

    .. NOTE::

        In its current state, this function is **not** able to build a `KTS(v)`
        for all `v` such that `v\equiv 3\pmod{6}`. It implements Theorems 5 and
        6 from [RCW71]_, i.e. constructions 1.1 and 1.2 from [Stinson91]_.

    INPUT:

    - `n` (integer)

    - ``existence`` (boolean; ``False`` by default) -- whether to build the
      `KTS(n)` or only answer whether it exists.

    .. SEEALSO::

        :meth:`IncidenceStructure.is_resolvable`

    EXAMPLES:

    A solution to Kirkmman's original problem::

        sage: kts = designs.kirkman_triple_system(15)
        sage: classes = kts.is_resolvable(1)[1]
        sage: names = '0123456789abcde'
        sage: to_name = lambda (r,s,t): ' '+names[r]+names[s]+names[t]+' '
        sage: rows = [join(('Day {}'.format(i) for i in range(1,8)), '   ')]
        sage: rows.extend(join(map(to_name,row), '   ') for row in zip(*classes))
        sage: print join(rows,'\n')
        Day 1   Day 2   Day 3   Day 4   Day 5   Day 6   Day 7
         07e     18e     29e     3ae     4be     5ce     6de
         139     24a     35b     46c     05d     167     028
         26b     03c     14d     257     368     049     15a
         458     569     06a     01b     12c     23d     347
         acd     7bd     78c     89d     79a     8ab     9bc

    TESTS:

    Construction 1.1 from [Stinson91]_::

        sage: for q in prime_powers(2,40):
        ....:     if q%6 == 1:
        ....:         _ = designs.kirkman_triple_system(2*q+1)

    Construction 1.2 from [Stinson91]_::

        sage: for q in prime_powers(2,40):
        ....:     if q%6 == 1:
        ....:         _ = designs.kirkman_triple_system(3*q)
    """
    if v%6 != 3:
        if existence:
            return False
        raise ValueError("There is no KTS({}) as v!=3 mod(6)".format(v))

    elif v == 3:
        if existence:
            return True
        return BalancedIncompleteBlockDesign(3,[[0,1,2]],k=3,lambd=1)

    elif v == 9:
        if existence:
            return True
        classes = [[[0, 1, 5], [2, 6, 7], [3, 4, 8]],
                   [[1, 6, 8], [3, 5, 7], [0, 2, 4]],
                   [[1, 4, 7], [0, 3, 6], [2, 5, 8]],
                   [[4, 5, 6], [0, 7, 8], [1, 2, 3]]]
        KTS = BalancedIncompleteBlockDesign(v,[tr for cl in classes for tr in cl],k=3,lambd=1,copy=False)
        KTS._classes = classes
        return KTS

    # Construction 1.1 from [Stinson91] (originally Theorem 6 from [RCW71])
    #
    # For all prime powers q=1 mod 6, there exists a KTS(2q+1)
    elif ((v-1)//2)%6 == 1 and is_prime_power((v-1)//2):
        if existence:
            return True
        from sage.rings.finite_rings.constructor import FiniteField as GF
        q = (v-1)//2
        K = GF(q,'x')
        a = K.primitive_element()
        t = (q-1)/6

        # m is the solution of a^m=(a^t+1)/2
        from sage.groups.generic import discrete_log
        m = discrete_log((a**t+1)/2, a)
        assert 2*a**m==a**t+1

        # First parallel class
        first_class = [[(0,1),(0,2),'inf']]
        b0 = K.one(); b1 = a**t; b2 = a**m
        first_class.extend([(b0*a**i,1),(b1*a**i,1),(b2*a**i,2)]
                            for i in range(t)+range(2*t,3*t)+range(4*t,5*t))
        b0 = a**(m+t); b1=a**(m+3*t); b2=a**(m+5*t)
        first_class.extend([[(b0*a**i,2),(b1*a**i,2),(b2*a**i,2)]
                            for i in range(t)])

        # Action of K on the points
        action = lambda v,x : (v+x[0],x[1]) if len(x) == 2 else x

        # relabel to integer
        relabel = {(p,x): i+(x-1)*q
                   for i,p in enumerate(K)
                   for x in [1,2]}
        relabel['inf'] = 2*q

        classes = [[[relabel[action(p,x)] for x in tr] for tr in first_class]
                   for p in K]

        KTS = BalancedIncompleteBlockDesign(v,[tr for cl in classes for tr in cl],k=3,lambd=1,copy=False)

        KTS._classes = classes
        return KTS

    # Construction 1.2 from [Stinson91] (originally Theorem 5 from [RCW71])
    #
    # For all prime powers q=1 mod 6, there exists a KTS(3q)
    elif (v//3)%6 == 1 and is_prime_power(v//3):
        if existence:
            return True
        from sage.rings.finite_rings.constructor import FiniteField as GF
        q = v//3
        K = GF(q,'x')
        a = K.primitive_element()
        t = (q-1)/6
        A0 = [(0,0),(0,1),(0,2)]
        B  = [[(a**i,j),(a**(i+2*t),j),(a**(i+4*t),j)] for j in range(3)
              for i in range(t)]
        A  = [[(a**i,0),(a**(i+2*t),1),(a**(i+4*t),2)] for i in range(6*t)]

        # Action of K on the points
        action = lambda v,x: (v+x[0],x[1])

        # relabel to integer
        relabel = {(p,j): i+j*q
                   for i,p in enumerate(K)
                   for j in range(3)}

        B0  = [A0] + B + A[t:2*t] + A[3*t:4*t] + A[5*t:6*t]

        # Classes
        classes = [[[relabel[action(p,x)] for x in tr] for tr in B0]
                   for p in K]

        for i in range(t)+range(2*t,3*t)+range(4*t,5*t):
            classes.append([[relabel[action(p,x)] for x in A[i]] for p in K])

        KTS = BalancedIncompleteBlockDesign(v,[tr for cl in classes for tr in cl],k=3,lambd=1,copy=False)
        KTS._classes = classes
        return KTS

    else:
        if existence:
            return Unknown
        raise NotImplementedError("I don't know how to build a ({},{},1)-BIBD!".format(v,3))

def v_4_1_rbibd(v,existence=False):
    r"""
    Return a `(v,4,1)`-RBIBD.

    INPUT:

    - `n` (integer)

    - ``existence`` (boolean; ``False`` by default) -- whether to build the
      design or only answer whether it exists.

    .. SEEALSO::

        - :meth:`IncidenceStructure.is_resolvable`
        - :func:`resolvable_balanced_incomplete_block_design`

    .. NOTE::

        A resolvable `(v,4,1)`-BIBD exists whenever `1\equiv 4\pmod(12)`. This
        function, however, only implements a construction of `(v,4,1)`-BIBD such
        that `v=3q+1\equiv 1\pmod{3}` where `q` is a prime power (see VII.7.5.a
        from [BJL99]_).

    EXAMPLE::

        sage: rBIBD = designs.resolvable_balanced_incomplete_block_design(28,4)
        sage: rBIBD.is_resolvable()
        True
        sage: rBIBD.is_t_design(return_parameters=True)
        (True, (2, 28, 4, 1))

    TESTS::

        sage: for q in prime_powers(2,30):
        ....:     if (3*q+1)%12 == 4:
        ....:         _ = designs.resolvable_balanced_incomplete_block_design(3*q+1,4) # indirect doctest
    """
    # Volume 1, VII.7.5.a from [BJL99]_
    if v%3 != 1 or not is_prime_power((v-1)//3):
        if existence:
            return Unknown
        raise NotImplementedError("I don't know how to build a ({},{},1)-RBIBD!".format(v,4))
    from sage.rings.finite_rings.constructor import FiniteField as GF
    q = (v-1)//3
    nn = (q-1)//4
    G = GF(q,'x')
    w = G.primitive_element()
    e = w**(nn)
    assert e**2 == -1

    first_class = [[(w**i,j),(-w**i,j),(e*w**i,j+1),(-e*w**i,j+1)]
                   for i in range(nn) for j in range(3)]

    first_class.append([(0,0),(0,1),(0,2),'inf'])

    label = {p:i for i,p in enumerate(G)}

    classes = [[[v-1 if x=='inf' else (x[1]%3)*q+label[x[0]+g] for x in S]
                for S in first_class]
               for g in G]

    BIBD = BalancedIncompleteBlockDesign(v,
                                         blocks = sum(classes,[]),
                                         k=4,
                                         check=True,
                                         copy=False)
    BIBD._classes = classes
    assert BIBD.is_resolvable()
    return BIBD
