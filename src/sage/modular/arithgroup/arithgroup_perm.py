r"""
Arithmetic subgroups defined by permutations

A theorem of Millington states that an arithmetic subgroup of index `N` is
uniquely determined by two elements generating a transitive subgroup of the
symmetric group `S_N` and satisfying a certain algebraic relation.

These functions are based on Chris Kurth's *KFarey* package.

AUTHORS:

- Chris Kurth (2008): created KFarey package

- David Loeffler (2009): adapted functions from KFarey for inclusion into Sage

"""

################################################################################
#
#       Copyright (C) 2009, The Sage Group -- http://www.sagemath.org/
#
#  Distributed under the terms of the GNU General Public License (GPL)
#
#  The full text of the GPL is available at:
#
#                  http://www.gnu.org/licenses/
#
################################################################################

from all import SL2Z
from arithgroup_generic import ArithmeticSubgroup
from sage.rings.all import Integers, ZZ

Lm = SL2Z([1,1,0,1])
Rm = SL2Z([1,0,1,1])

class ArithmeticSubgroup_Permutation(ArithmeticSubgroup):
    r"""
    An arithmetic subgroup `\Gamma` defined by two permutations, giving the
    action of the parabolic generators

    .. math::

        L = \begin{pmatrix} 1 & 1 \\ 0 & 1\end{pmatrix},\quad R = \begin{pmatrix} 1 & 0 \\ 1 & 1 \end{pmatrix}

    by right multiplication on the coset representatives `\Gamma \backslash {\rm SL}_2(\mathbb{Z})`.


    EXAMPLES:

    We construct a noncongruence subgroup of index 7 (the smallest possible)::

        sage: a2 = SymmetricGroup(7)([(1,2),(3,4),(5,6)]); a3 = SymmetricGroup(7)([(1,3,5),(2,6,7)])
        sage: G = ArithmeticSubgroup_Permutation(a2*a3, ~a2 * ~a3); G
        Arithmetic subgroup corresponding to permutations L=(1,6)(2,3,4,5,7), R=(1,7,6,3,4)(2,5)
        sage: G.index()
        7
        sage: G.dimension_cusp_forms(4)
        1
        sage: G.is_congruence()
        False

    We convert some standard congruence subgroups into permutation form::

        sage: Gamma0(12).as_permutation_group()
        Arithmetic subgroup corresponding to permutations L=(2,3,4,5,6,7,8,9,10,11,12,13)(14,15,16)(17,19,20,18)(21,23,22), R=(1,3,14,17,21,7,24,9,23,20,16,13)(4,18,12)(5,22,11,15)(6,10,19)

    The following is the unique index 2 even subgroup of `{\rm SL}_2(\mathbb{Z})`::

        sage: w = SymmetricGroup(2)([2,1])
        sage: G = ArithmeticSubgroup_Permutation(w, w)
        sage: G.dimension_cusp_forms(6)
        1
        sage: G.genus()
        0

    We test unpickling::

        sage: G == loads(dumps(G))
        True
        sage: G is loads(dumps(G))
        False
    """

    def __init__(self, L, R):
        r"""
        Create an arithmetic subgroup given two permutations.

            sage: w = SymmetricGroup(2)([2,1])
            sage: ArithmeticSubgroup_Permutation(w, w) # indirect doctest
            Arithmetic subgroup corresponding to permutations L=(1,2), R=(1,2)
        """

        from sage.groups.perm_gps.all import PermutationGroup
        self.L = L
        self.R = R
        A1=(L*R**(-1)*L)**2
        A2=(R**(-1)*L)**3
        if A1.order() == 1 and A2.order() == 1:
            # This forces G to be even! -- check this.
            self._permgp = PermutationGroup([L, R])
            if not self._permgp.is_transitive():
                raise ValueError, "Permutations are not transitive"
        else:
            raise ValueError, "Permutations do not satisfy defining relation"

    def __reduce__(self):
        r"""
        Return the data used to construct self. Used in pickling.

        EXAMPLE::

            sage: w = SymmetricGroup(2)([2,1])
            sage: ArithmeticSubgroup_Permutation(w, w).__reduce__()
            (<class 'sage.modular.arithgroup.arithgroup_perm.ArithmeticSubgroup_Permutation'>, ((1,2), (1,2)))
        """
        return (ArithmeticSubgroup_Permutation, (self.L, self.R))

    def perm_group(self):
        r"""
        Return the underlying permutation group.

        EXAMPLE::

            sage: import sage.modular.arithgroup.arithgroup_perm as ap
            sage: ap.HsuExample10().perm_group()
            Permutation Group with generators [(1,4)(2,5,9,10,8)(3,7,6), (1,7,9,10,6)(2,3)(4,5,8)]
        """
        return self._permgp

    def index(self):
        r"""
        Return the index of self in the full modular group.

        EXAMPLE::

            sage: import sage.modular.arithgroup.arithgroup_perm as ap
            sage: ap.HsuExample18().index()
            18
        """
        return self._permgp.degree()

    def _repr_(self):
        r"""
        String representation of self.

        EXAMPLE::

            sage: import sage.modular.arithgroup.arithgroup_perm as ap
            sage: ap.HsuExample18()._repr_()
            'Arithmetic subgroup corresponding to permutations L=(1,2)(3,4)(5,6,7)(8,9,10)(11,12,13,14,15,16,17,18), R=(1,12,18)(2,6,13,9,4,8,17,7)(3,16,14)(5,11)(10,15)'
        """
        return "Arithmetic subgroup corresponding to permutations L=%s, R=%s" % (self.L, self.R)

    def permutation_action(self, x):
        r"""
        Given an element x of `{\rm SL}_2(\mathbb{Z})`, compute the
        permutation of the cosets of self given by right multiplication by x.

        EXAMPLE::

            sage: import sage.modular.arithgroup.arithgroup_perm as ap
            sage: ap.HsuExample10().permutation_action(SL2Z([32, -21, -67, 44]))
            (1,10,5,6,3,8,9,2,7,4)
        """
        x = SL2Z(x)
        w = sl2z_word_problem(x)
        p = LREvalPerm(w, self.L, self.R)
        return p

    def __call__(self, g, check=True):
        r"""
        Create an element of this group from g. If check=True (the default),
        perform the (possibly rather computationally-intensive) check to make
        sure g is in this group.

        EXAMPLE::

            sage: import sage.modular.arithgroup.arithgroup_perm as ap
            sage: ap.HsuExample10()(SL2Z([1,1,0,1]))
            Traceback (most recent call last):
            ...
            TypeError: Not in group
            sage: ap.HsuExample10()(SL2Z([1,1,0,1]), check=False)
            [1 1]
            [0 1]
        """
        g = SL2Z(g)
        if (not check) or (self.permutation_action(g)(1) == 1):
            g._parent = self
            return g
        else:
            raise TypeError, "Not in group"

    def __cmp__(self, other):
        r"""
        Compare self to other.

        EXAMPLES::

            sage: G = ArithmeticSubgroup_Permutation(SymmetricGroup(2)([2,1]), SymmetricGroup(2)([2,1]))
            sage: cmp(G, 1) in [1,-1]
            True
            sage: cmp(G, Gamma0(8))
            1
            sage: cmp(G, G)
            0
        """
        if not isinstance(other, ArithmeticSubgroup_Permutation):
            return cmp(type(self), type(other))
        else:
            return cmp( (self.L, self.R), (other.L, other.R) )

    def is_congruence(self):
        r"""
        Return True if this is a congruence subgroup. Uses Hsu's algorithm, as
        implemented by Chris Kurth in KFarey.

        EXAMPLES:

        This example is congruence -- it's Gamma0(3) in disguise: ::

            sage: import sage.modular.arithgroup.arithgroup_perm as ap
            sage: G=ap.ArithmeticSubgroup_Permutation(SymmetricGroup(4)((2,3,4)), SymmetricGroup(4)((1,3,4))); G
            Arithmetic subgroup corresponding to permutations L=(2,3,4), R=(1,3,4)
            sage: G.is_congruence()
            True

        This one is noncongruence: ::

            sage: import sage.modular.arithgroup.arithgroup_perm as ap
            sage: ap.HsuExample10().is_congruence()
            False
        """
        L = self.L
        R = self.R
        N=L.order()
        e=2**(ZZ(N).valuation(2))
        m=N/e
        if e==1:     #i.e. N is odd
            onehalf=int(Integers(N)(2)**(-1))   #i.e. 2^(-1) mod N
            check=(R**2*L**(-onehalf))**3
            if check.order() == 1:
                return True       #Congruence
            else:
                return False      #Noncongruence
        elif m==1:     #i.e. N is a power of 2
            onefifth=int(Integers(N)(5)**(-1))   #i.e. 2^(-1) mod N
            S=L**20*R**onefifth*L**(-4)*R**(-1)
            #Congruence if these three permutations are trivial:
            Rel1=(L*R**(-1)*L)**(-1) * S * (L*R**(-1)*L) * S
            Rel2=S**(-1)*R*S*R**(-25)
            Rel3=(S*R**5*L*R**(-1)*L)**3
            if Rel1.order()==1 and Rel2.order()==1 and Rel3.order()==1:
                return True       #Congruence
            else:
                return False      #Noncongruence
        else:         #i.e. e>1, M>1
            onehalf=int(Integers(m)(2)**(-1))    #i.e. 2^(-1)  mod m
            onefifth=int(Integers(e)(5)**(-1))   #i.e. 2^(-1)  mod e
            c=int(Integers(m)(e)**(-1))*e        #i.e. e^(-1)e mod m
            d=int(Integers(e)(m)**(-1))*m        #i.e. m^(-1)m mod e
            a=L**c
            b=R**c
            l=L**d
            r=R**d
            s=l**20*r**onefifth*l**(-4)*r**(-1)

            Rel1=a**(-1)*r**(-1)*a*r
            Rel2=(a*b**(-1)*a)**4
            Rel3=(a*b**(-1)*a)**2*(b**(-1)*a)**(-3)
            Rel4=(a*b**(-1)*a)**2*(b**2*a**(-onehalf))**(-3)
            Rel5=(l*r**(-1)*l)**(-1)*s*(l*r**(-1)*l)*s
            Rel6=s**(-1)*r*s*r**(-25)
            Rel7=(l*r**(-1)*l)**2*(s*r**5*l*r**(-1)*l)**(-3)

            if (Rel1.order()==1
                and Rel2.order()==1
                and Rel3.order()==1
                and Rel4.order()==1
                and Rel5.order()==1
                and Rel6.order()==1
                and Rel7.order()==1):
                return True
            else:
                return False


def LREvalPerm(w, L, R):
    r"""
    Given a word w as output by sl2z_word_problem, evaluate the word with the
    given permutations for L and R. Because we are dealing with a right rather
    than a left action, arguments are evaluated back to front.

    EXAMPLE::

        sage: import  sage.modular.arithgroup.arithgroup_perm as ap
        sage: L = SymmetricGroup(4)('(1,2)(3,4)'); R = SymmetricGroup(4)('(1,2,3,4)')
        sage: ap.LREvalPerm([(1,1),(0,1)], L, R) ==  L * R
        True
    """
    from sage.groups.perm_gps.all import PermutationGroup
    G=PermutationGroup([L,R])       #Change this?
    Id=G.identity()
    L=G(L)
    R=G(R)

    checkM=Id

    for i in range(len(w)):
        if w[i][0]==0:
            c=int(w[i][1])
            checkM=L**c*checkM
        elif w[i][0]==1:
            c=int(w[i][1])
            checkM=R**c*checkM
    return checkM

def sl2z_word_problem(A, max_iterations=20):
    r"""
    Given an element of SL2Z, express it as a word in the generators L =
    [1,1,0,1] and R = [1,0,1,1].

    The return format is a list of pairs (a,b), where a = 0 or 1 denoting L or
    R respectively, and b is an integer exponent.

    The parameter iterations (default 20) controls the maximum number of
    iterations to allow in the program's main loop; an error is raised if the
    algorithm has not terminated after this many iterations.

    EXAMPLE::

        sage: import sage.modular.arithgroup.arithgroup_perm as ap
        sage: L = SL2Z([1,1,0,1]); R = SL2Z([1,0,1,1])
        sage: ap.sl2z_word_problem(L)
        [(0, 1)]
        sage: ap.sl2z_word_problem(R**(-1))
        [(1, -1)]
        sage: ap.sl2z_word_problem(L*R)
        [(0, 1), (1, -1), (0, 0), (1, 2)]
    """

    A = SL2Z(A)
    output=[]

    ## If A00 is zero, add A01 to it (probably a better way to do this)
    if A[0,0]==0:
        A=A*Rm
        output=[(1,-1)]+output

    if A[0,0]<0:   # Make sure A00 is positive
        A=SL2Z(-1)*A
        output= [(1,-1), (0,1), (1,-1), (0,1), (1,-1), (0,1)] + output

    if A[0,1]<0:   # if A01 is negative make it positive
        n=(-A[0,1]/A[0,0]).ceil()  #n s.t. 0 <= A[0,1]+n*A[0,0] < A[0,0]
        A=A*Lm**n
        output=[(0, -n)]+output
   ## At this point A00>0 and A01>=0
    timesthrough=0
    while not (A[0,0]==0 or A[0,1]==0):
        if A[0,0]>A[0,1]:
            n=(A[0,0]/A[0,1]).floor()
            A=A*Rm**(-n)
            output=[(1, n)]+output

        else:      #A[0,0]<A[0,1]
            n=(A[0,1]/A[0,0]).floor()
            A=A*Lm**(-n)
            output=[(0, n)]+output
        if timesthrough > max_iterations:
            raise ValueError, "Over-recursion"
        timesthrough=timesthrough+1

    if A==SL2Z(1):
        pass       #done, so don't add R^0
    elif A[0,0]==0:
        c=A[1,1]
        A=A*Lm**(c-1)*Rm*Lm**(-1)
        output=[(0,1),(1,-1),(0, 1-c)]+output
    else:
        c=A[1,0]
        A=A*Rm**(-c)
        output=[(1,c)]+output

    return output

def eval_word(B):
    r"""
    Given a word in the format output by sl2z_word_problem, convert it back
    into an element of SL2(Z).

    EXAMPLES::

        sage: import sage.modular.arithgroup.arithgroup_perm as ap
        sage: ap.eval_word([(0, 1), (1, -1), (0, 0), (1, 3), (0, 2), (1, 9), (0, -1)])
        [ 66 -59]
        [ 47 -42]
    """
    checkM=SL2Z(1)
    for i in range(len(B)):
        if B[i][0]==0:
            c=int(B[i][1])
            checkM=checkM*Lm**c
        elif B[i][0]==1:
            c=int(B[i][1])
            checkM=checkM*Rm**c
    return checkM

def HsuExample10():
    r"""
    An example of an index 10 arithmetic subgroup studied by Tim Hsu.

    EXAMPLE::

        sage: import sage.modular.arithgroup.arithgroup_perm as ap
        sage: ap.HsuExample10()
        Arithmetic subgroup corresponding to permutations L=(1,4)(2,5,9,10,8)(3,7,6), R=(1,7,9,10,6)(2,3)(4,5,8)
    """
    from sage.groups.perm_gps.all import SymmetricGroup
    return ArithmeticSubgroup_Permutation(SymmetricGroup(10)("(1,4)(2,5,9,10,8)(3,7,6)"), SymmetricGroup(10)("(1,7,9,10,6)(2,3)(4,5,8)"))

def HsuExample18():
    r"""
    An example of an index 18 arithmetic subgroup studied by Tim Hsu.

    EXAMPLE::

        sage: import sage.modular.arithgroup.arithgroup_perm as ap
        sage: ap.HsuExample18()
        Arithmetic subgroup corresponding to permutations L=(1,2)(3,4)(5,6,7)(8,9,10)(11,12,13,14,15,16,17,18), R=(1,12,18)(2,6,13,9,4,8,17,7)(3,16,14)(5,11)(10,15)
    """
    from sage.groups.perm_gps.all import SymmetricGroup
    return ArithmeticSubgroup_Permutation(SymmetricGroup(18)("(1,2)(3,4)(5,6,7)(8,9,10)(11,12,13,14,15,16,17,18)"), SymmetricGroup(18)("(1,12,18)(2,6,13,9,4,8,17,7)(3,16,14)(5,11)(10,15)"))

def convert_to_permgroup(G):
    r"""
    Given an arbitrary arithmetic subgroup, convert it to permutation form.

    Note that the permutation representation is not always unique, so if G is
    already of permutation type, then the return value won't necessarily be
    identical to G, but it will represent the same subgroup.

    EXAMPLES::

        sage: import sage.modular.arithgroup.arithgroup_perm as ap
        sage: ap.convert_to_permgroup(Gamma0(5))
        Arithmetic subgroup corresponding to permutations L=(2,3,4,5,6), R=(1,3,5,4,6)
        sage: ap.convert_to_permgroup(ap.HsuExample10())
        Arithmetic subgroup corresponding to permutations L=(1,2)(3,5,6,7,8)(4,9,10), R=(1,9,6,7,10)(2,5,8)(3,4)
    """

    from all import is_ArithmeticSubgroup
    from sage.groups.perm_gps.all import SymmetricGroup
    if not is_ArithmeticSubgroup(G):
        raise TypeError, "%s is not an arithmetic subgroup of SL2(Z)" % G

    if not G.is_even():
        raise NotImplementedError, "Only subgroups containing -1 are currently implemented"

    X = list(G.coset_reps())

    if not (X[0] in G):
        raise Exception, "Coset reps wrongly ordered!"

    # find permutation action of L
    perm = []
    for i in xrange(len(X)):
        for j in xrange(len(X)):
            if (X[i] * Lm) * X[j]**(-1) in G:
                perm.append(j)
                break
        if len(perm) != i+1:
            raise ArithmeticError
    Lperm = SymmetricGroup(len(X))([i+1 for i in perm])

    # same for R
    perm = []
    for i in xrange(len(X)):
        for j in xrange(len(X)):
            if (X[i] * Rm) * X[j]**(-1) in G:
                perm.append(j)
                break
        if len(perm) != i+1:
            raise ArithmeticError
    Rperm = SymmetricGroup(len(X))([i+1 for i in perm])
    return ArithmeticSubgroup_Permutation(Lperm, Rperm)
