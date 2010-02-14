# -*- coding: utf-8 -*-
#import ell_point
#import formal_group
#import ell_torsion
#from ell_generic import EllipticCurve_generic, is_EllipticCurve
#from ell_number_field import EllipticCurve_number_field

#import sage.groups.all
import sage.rings.arith as arith
import sage.rings.all as rings
from sage.rings.all import ZZ

def prove_BSD(self, verbosity=0, simon=False, proof=None, secs_hi=30):
    r"""
    Attempts to prove the Birch and Swinnerton-Dyer conjectural
    formula for `E`, returning a list of primes `p` for which this
    function fails to prove BSD(E,p).  Here, BSD(E,p) is the
    statement: "the Birch and Swinnerton-Dyer formula holds up to a
    rational number coprime to `p`."

    INPUT:

        - ``verbosity`` - int, how much information about the proof to print.

            - 0 - print nothing
            - 1 - print sketch of proof
            - 2 - print information about remaining primes

        - ``simon`` - bool (default False), whether to use two_descent or
          simon_two_descent at p=2.

        - ``proof`` - bool or None (default: None, see
          proof.elliptic_curve or sage.structure.proof). If False, this
          function just immediately returns the empty list.

        - ``secs_hi`` - maximum number of seconds to try to compute the
          Heegner index before switching over to trying to compute the
          Heegner index bound. (Rank 0 only!)

    NOTE:

    When printing verbose output, phrases such as "by Mazur" are referring
    to the following list of papers:

    REFERENCES:

    .. [Cha] B. Cha. Vanishing of some cohomology goups and bounds for the
       Shafarevich-Tate groups of elliptic curves. J. Number Theory, 111:154-
       178, 2005.
    .. [Jetchev] D. Jetchev. Global divisibility of Heegner points and
       Tamagawa numbers. Compos. Math. 144 (2008), no. 4, 811--826.
    .. [Kato] K. Kato. p-adic Hodge theory and values of zeta functions of
       modular forms. Astérisque, (295):ix, 117-290, 2004.
    .. [Kolyvagin] V. A. Kolyvagin. On the structure of Shafarevich-Tate
       groups. Algebraic geometry, 94--121, Lecture Notes in Math., 1479,
       Springer, Berlin, 1991.
    .. [LumStein] A. Lum, W. Stein. Verification of the Birch and
       Swinnerton-Dyer Conjecture for Elliptic Curves with Complex
       Multiplication (unpublished)
    .. [Mazur] B. Mazur. Modular curves and the Eisenstein ideal. Inst.
       Hautes Études Sci. Publ. Math. No. 47 (1977), 33--186 (1978).
    .. [Rubin] K. Rubin. The "main conjectures" of Iwasawa theory for
       imaginary quadratic fields. Invent. Math. 103 (1991), no. 1, 25--68.
    .. [SteinWuthrich] W. Stein and C. Wuthrich. Computations about
       Tate-Shafarevich groups using Iwasawa theory.
       http://wstein.org/papers/shark, February 2008.
    .. [SteinEtAl] G. Grigorov, A. Jorza, S. Patrikis, W. Stein,
       C. Tarniţǎ. Computational verification of the Birch and
       Swinnerton-Dyer conjecture for individual elliptic curves.
       Math. Comp. 78 (2009), no. 268, 2397--2425.


    EXAMPLES::

        sage: EllipticCurve('11a').prove_BSD(verbosity=2)
        p = 2: True by 2-descent
        True for p not in {2, 5} by Kolyvagin.
        True for p=5 by Mazur
        []

        sage: EllipticCurve('14a').prove_BSD(verbosity=2)
        p = 2: True by 2-descent
        True for p not in {2, 3} by Kolyvagin.
        Remaining primes:
        p = 3: reducible, not surjective, good ordinary, divides a Tamagawa number
        [3]
        sage: EllipticCurve('14a').prove_BSD(simon=True)
        [3]

    A rank two curve::

        sage: E = EllipticCurve('389a')

    We know nothing with proof=True::

        sage: E.prove_BSD()
        Set of all prime numbers: 2, 3, 5, 7, ...

    We (think we) know everything with proof=False::

        sage: E.prove_BSD(proof=False)
        []

    A curve of rank 0 and prime conductor::

        sage: E = EllipticCurve('19a')
        sage: E.prove_BSD(verbosity=2)
        p = 2: True by 2-descent
        True for p not in {2, 3} by Kolyvagin.
        True for p=3 by Mazur
        []

        sage: E = EllipticCurve('37a')
        sage: E.rank()
        1
        sage: E._EllipticCurve_rational_field__rank
        {True: 1}
        sage: E.analytic_rank = lambda : 0
        sage: E.prove_BSD()
        Traceback (most recent call last):
        ...
        RuntimeError: It seems that the rank conjecture does not hold for this curve (Elliptic Curve defined by y^2 + y = x^3 - x over Rational Field)! This may be a counterexample to BSD, but is more likely a bug.

    We test the consistency check for the 2-part of Sha::

        sage: E = EllipticCurve('37a')
        sage: S = E.sha(); S
        Shafarevich-Tate group for the Elliptic Curve defined by y^2 + y = x^3 - x over Rational Field
        sage: S.an = lambda : 4
        sage: E.prove_BSD()
        Traceback (most recent call last):
        ...
        RuntimeError: Two descent => ord_2(#Sha[2]) = 0, but ord_2(#sha_an) = 2.

    An example with a Tamagawa number at 5::

        sage: E = EllipticCurve('123a1')
        sage: E.prove_BSD(verbosity=2)
        p = 2: True by 2-descent
        True for p not in {2, 5} by Kolyvagin.
        Remaining primes:
        p = 5: reducible, not surjective, good ordinary, divides a Tamagawa number
        [5]

    A curve for which 3 divides the order of the Shafarevich-Tate group::

        sage: E = EllipticCurve('681b')
        sage: E.prove_BSD(verbosity=2)               # long time
        p = 2: True by 2-descent
        True for p not in {2, 3} by Kolyvagin.
        ALERT: p = 3 left in Kolyvagin bound
            0 <= ord_p(#Sha) <= 2
            ord_p(#Sha_an) = 2
        Remaining primes:
        p = 3: irreducible, surjective, non-split multiplicative
        [3]

    A curve for which we need to use ``heegner_index_bound``::

        sage: E = EllipticCurve('198b')
        sage: E.prove_BSD(verbosity=1, secs_hi=1)
        p = 2: True by 2-descent
        Timeout stopped Heegner index computation...
        Proceeding to use heegner_index_bound instead.
        True for p not in {2, 3} by Kolyvagin.
        [3]

    TESTS:

    This was fixed by trac #8184 and #7575::

        sage: EllipticCurve('438e1').prove_BSD(verbosity=1)
        p = 2: True by 2-descent
        True for p not in {2} by Kolyvagin.
        []

    ::

        sage: E = EllipticCurve('960d1')
        sage: E.prove_BSD(verbosity=1)
        p = 2: True by 2-descent
        Timeout stopped Heegner index computation...
        Proceeding to use heegner_index_bound instead.
        True for p not in {2} by Kolyvagin.
        []

    """
    if proof is None:
        from sage.structure.proof.proof import get_flag
        proof = get_flag(proof, "elliptic_curve")
    else:
        proof = bool(proof)
    if not proof:
        return []
    two_tor_rk = self.two_torsion_rank()
    Sha = self.sha()
    sha_an = Sha.an()
    if simon:
        rank_lower_bd, two_sel_rk_bd = self.simon_two_descent()[:2]
        two_sel_rk_bd -= two_tor_rk
        sha2_lower_bd = 0
    else:
        MWRC = self.mwrank_curve()
        two_sel_rk_bd = MWRC.selmer_rank()
        two_sel_rk_bd -= two_tor_rk
        rank_lower_bd = MWRC.rank()
        sha2_lower_bd = two_sel_rk_bd - MWRC.rank_bound()
    sha2_upper_bd = two_sel_rk_bd - rank_lower_bd
    # note: two_sel_rk_bd will include sha
    rank = None
    if rank_lower_bd > 1:
        # We do not know BSD(E,p) for even a single p, since it's
        # an open problem to show that L^r(E,1)/(Reg*Omega) is
        # rational for any curve with r >= 2.
        from sage.sets.all import Primes
        return Primes()
    if sha2_upper_bd == sha2_lower_bd:
        rank = rank_lower_bd
        if sha_an.ord(2) != sha2_lower_bd:
            raise RuntimeError("Two descent => ord_2(#Sha[2]) = %d, but ord_2(#sha_an) = %d."%(sha2_lower_bd,sha_an.ord(2)))
        if verbosity > 0:
            print 'p = 2: True by 2-descent'
        two_proven = True
    else:
        if sha2_upper_bd < sha2_lower_bd:
            raise RuntimeError("Apparent contradiction: ord_2(#Sha[2]_an) == %d, rank(2-Sel)-rank(2-tor) = %d, rank >= %d, ord_2(#Sha[2]) >= %d"%(sha_an.ord(2),two_sel_rk_bd,rank_lower_bd,sha2_lower_bd))
        if two_sel_rk_bd - rank_lower_bd == sha_an.ord(2):
            if verbosity > 0:
                print 'p = 2: ord_2(#Sha[2]_an) == %d >= ord_2(#Sha[2]).'%sha_an.ord(2)
        elif two_sel_rk_bd - rank_lower_bd > sha_an.ord(2):
            if verbosity > 0:
                print 'p = 2: ord_2(#Sha[2]_an) == %d, and ord_2(#Sha[2]) <= %d.'%(sha_an.ord(2),two_sel_rk_bd - rank_lower_bd)
        else:
            raise RuntimeError("Apparent contradiction: ord_2(#Sha[2]_an) == %d, rank(2-Sel)-rank(2-tor) = %d, rank >= %d"%(sha_an.ord(2),two_sel_rk_bd,rank_lower_bd))
        two_proven = False
        if verbosity > 0:
            print 'Looking rank up in database...'
        rank = self.rank(use_database=True, only_use_mwrank=False)

    if rank != self.analytic_rank():
        raise RuntimeError("It seems that the rank conjecture does not hold for this curve (%s)! This may be a counterexample to BSD, but is more likely a bug."%(self))

    # We replace self by the optimal curve, which we can do since
    # truth of BSD(E,p) is invariant under isogeny.
    self = self.optimal_curve()

    N = self.conductor()

    # reduce set of remaining primes to a finite set
    import signal
    remaining_primes = []
    kolyvagin_primes = []
    heegner_index = None
    if self.rank() == 0:
        try:
            old_alarm = signal.alarm(secs_hi)
            old_alarm_set = (old_alarm != 0)
            for D in self.heegner_discriminants_list(10):
                I = None
                while I is None:
                    dsl=15
                    try:
                        I = self.heegner_index(D, descent_second_limit=dsl)
                    except RuntimeError as err:
                        if err.args[0][-33:] == 'Generators not provably computed.':
                            dsl += 1
                        else: raise RuntimeError(err)
                J = I.is_int()
                if J[0] and J[1]>0:
                    I = J[1]
                else:
                    J = (2*I).is_int()
                    if J[0] and J[1]>0:
                        I = J[1]
                    else:
                        I = None
                if I is not None:
                    if heegner_index is None:
                        heegner_index = I
                        break # no big long loops just yet...
            old_alarm_sub = signal.alarm(0)
            if old_alarm_set:
                old_alarm -= old_alarm_sub
        except KeyboardInterrupt:
            if signal.alarm(0)==0:
                print 'Timeout stopped Heegner index computation...'
                print 'Proceeding to use heegner_index_bound instead.'
            else:
                raise KeyboardInterrupt
        except RuntimeError as err:
            if err.args[0][:37] == 'End Of File (EOF) in read_nonblocking':
                print 'Computing Heegner index failed due to mwrank interface: see #7575.'
            else: raise RuntimeError(err)
            old_alarm_sub = signal.alarm(0)
            if old_alarm_set:
                old_alarm -= old_alarm_sub
                if old_alarm <= 0:
                    raise KeyboardInterrupt
                signal.alarm(old_alarm)
        if old_alarm_set: # in case alarm was already set...
            if old_alarm <= 0:
                raise KeyboardInterrupt
            signal.alarm(old_alarm)
        if heegner_index is None:
            for D in self.heegner_discriminants_list(100):
                max_height = 12
                heegner_primes = -1
                while heegner_primes == -1:
                    max_height += 1
                    heegner_primes, _ = self.heegner_index_bound(D, max_height=max_height)
                if isinstance(heegner_primes, list):
                    break
            if not isinstance(heegner_primes, list):
                raise RuntimeError("Tried 100 Heegner discriminants, and heegner_index_bound failed each time.")
            if 2 in heegner_primes:
                heegner_primes.remove(2)
        else:
            heegner_primes = [p for p in arith.prime_divisors(heegner_index) if p!=2]
    else: # rank 1
        for D in self.heegner_discriminants_list(10):
            I = self.heegner_index(D)
            J = I.is_int()
            if J[0] and J[1]>0:
                I = J[1]
            else:
                J = (2*I).is_int()
                if J[0] and J[1]>0:
                    I = J[1]
                else:
                    continue
            heegner_index = I
            break
        heegner_primes = [p for p in arith.prime_divisors(heegner_index) if p!=2]

    if self.has_cm():
        # ensure that CM is by a maximal order
        non_max_j_invs = [-12288000, 54000, 287496, 16581375]
        if self.j_invariant() in non_max_j_invs:
            for E in self.isogeny_class()[0]:
                if E.j_invariant() not in non_max_j_invs:
                    Sha = E.sha()
                    sha_an = Sha.an()
                    if verbosity > 0:
                        print 'CM by non maximal order: switching curves'
                    break
        else:
            E = self
        if E.analytic_rank() == 0:
            if verbosity > 0:
                print 'p >= 5: true by Rubin'
            remaining_primes.append(3)
        else:
            K = rings.QuadraticField(E.cm_discriminant(), 'a')
            D_K = K.disc()
            D_E = E.discriminant()
            if len(K.factor(3)) == 1: # 3 does not split in K
                remaining_primes.append(3)
            for p in arith.prime_divisors(D_K):
                if p >= 5:
                    remaining_primes.append(p)
            for p in arith.prime_divisors(D_E):
                if p >= 5 and D_K%p and len(K.factor(p)) == 1: # p is inert in K
                    remaining_primes.append(p)
            for p in heegner_primes:
                if p >= 5 and D_E%p != 0 and D_K%p != 0 and len(K.factor(p)) == 1: # p is good for E and inert in K
                    kolyvagin_primes.append(p)
            assert sha_an in ZZ and sha_an > 0
            for p in arith.prime_divisors(sha_an):
                if p >= 5 and D_K%p != 0 and len(K.factor(p)) == 1:
                    if E.is_good(p):
                        if verbosity > 2 and p in heegner_primes and heegner_index is None:
                            print 'ALERT: Prime p (%d) >= 5 dividing sha_an, good for E, inert in K, in heegner_primes, should not divide the actual Heegner index'
                        # Note that the following check is not entirely
                        # exhaustive, in case there is a p not dividing
                        # the Heegner index in heegner_primes,
                        # for which only an outer bound was computed
                        if p not in heegner_primes:
                            raise RuntimeError("p = %d divides sha_an, is of good reduction for E, inert in K, and does not divide the Heegner index. This may be a counterexample to BSD, but is more likely a bug. %s"%(p,self))
            if verbosity > 0:
                print 'True for p not in {%s} by Kolyvagin (via Stein & Lum -- unpublished) and Rubin.'%str(list(set(remaining_primes).union(set(kolyvagin_primes))))[1:-1]
    else: # no CM
        E = self
        # do some tricks to get to a finite set without calling bound_kolyvagin
        remaining_primes = E.galois_representation().non_surjective()
        for p in heegner_primes:
            if p not in remaining_primes:
                remaining_primes.append(p)
        assert sha_an in ZZ and sha_an > 0
        for p in arith.prime_divisors(sha_an):
            if p not in remaining_primes:
                remaining_primes.append(p)
        if 2 in remaining_primes: remaining_primes.remove(2)
        if verbosity > 0:
            print 'True for p not in {' + str([2]+list(remaining_primes))[1:-1] + '} by Kolyvagin.'
        primes_to_remove = []
        for p in remaining_primes:
            if E.galois_representation().is_surjective(p) and not E.has_additive_reduction(p):
                if E.has_nonsplit_multiplicative_reduction(p):
                    if E.rank() > 0:
                        continue
                if p==3:
                    if (not (E.is_ordinary(p) and E.is_good(p))) and (not E.has_split_multiplicative_reduction(p)):
                        continue
                    if E.rank() > 0:
                        continue
                if verbosity > 1:
                    print 'p = %d: Trying p_primary_bound'%p
                p_bound = Sha.p_primary_bound(p)
                if sha_an.ord(p) == 0 and p_bound == 0:
                    if verbosity > 0:
                        print 'True for p=%d by Stein-Wuthrich.'%p
                    primes_to_remove.append(p)
                else:
                    print 'Analytic %d-rank is '%p + str(sha_an.ord(p)) + ', actual %d-rank is at most %d.'%(p, p_bound)
                    print '    by Stein-Wuthrich.\n'
        for p in primes_to_remove:
            remaining_primes.remove(p)
        kolyvagin_primes = []
        for p in remaining_primes:
            if E.galois_representation().is_surjective(p):
                kolyvagin_primes.append(p)
        for p in kolyvagin_primes:
            remaining_primes.remove(p)
    # apply other hypotheses which imply Kolyvagin's bound holds
    bounded_primes = []
    D_K = rings.QuadraticField(D, 'a').disc()
    assert 2 not in remaining_primes
    # Cha's hypothesis
    for p in remaining_primes:
        if D_K%p != 0 and N%(p**2) != 0 and E.galois_representation().is_irreducible(p):
            if verbosity > 0:
                print 'Kolyvagin\'s bound for p = %d applies by Cha.'%p
            kolyvagin_primes.append(p)
    # Stein et al.
    if not E.has_cm():
        L = arith.lcm([F.torsion_order() for F in E.isogeny_class()[0]])
        for p in remaining_primes:
            if p in kolyvagin_primes: continue
            if L%p != 0:
                if len(arith.prime_divisors(D_K)) == 1:
                    if D_K%p != 0:
                        if verbosity > 0:
                            print 'Kolyvagin\'s bound for p = %d applies by Stein et al.'%p
                        kolyvagin_primes.append(p)
                else:
                    if verbosity > 0:
                        print 'Kolyvagin\'s bound for p = %d applies by Stein et al.'%p
                    kolyvagin_primes.append(p)
    for p in kolyvagin_primes:
        if p in remaining_primes:
            remaining_primes.remove(p)

    prime_bounds = []
    # apply Kolyvagin's bound
    primes_to_remove = []
    for p in kolyvagin_primes:
        if sha_an.ord(p) == 0 and p not in heegner_primes:
                if verbosity > 0:
                    print 'True for p = %d by Kolyvagin bound.'%p
                primes_to_remove.append(p)
                continue
        if heegner_index is not None: # p must divide heegner_index
            ord_p_bound = 2*heegner_index.ord(p)
            # Here Jetchev's results apply.
            m_max = max([E.tamagawa_number(q).ord(p) for q in N.prime_divisors()])
            if m_max > 0 and verbosity > 0:
                print 'Jetchev\'s results apply (at p = %d) with m_max ='%p, m_max
            ord_p_bound -= 2*m_max
            if ord_p_bound == 0:
                if sha_an.ord(p) != 0:
                    raise RuntimeError("p = %d: ord_p_bound == 0, but sha_an.ord(p) == %d. This appears to be a counterexample to BSD, but is more likely a bug."%(p,sha_an.ord(p)))
                if verbosity > 0:
                    print 'True for p = %d by Kolyvagin bound.'%p
                primes_to_remove.append(p)
                continue
        elif p not in heegner_primes:
            ord_p_bound = 0
        else:
            from sage.rings.infinity import Infinity
            ord_p_bound = Infinity
            if verbosity > 0:
                print 'p = %d may divide the Heegner index, for which only a bound was computed.'%p
        if verbosity > 0:
            print 'ALERT: p = %d left in Kolyvagin bound'%p
            print '    0 <= ord_p(#Sha) <=', ord_p_bound
            print '    ord_p(#Sha_an) =', sha_an.ord(p)
    for p in primes_to_remove:
        kolyvagin_primes.remove(p)
    remaining_primes = list( set(remaining_primes).union(set(kolyvagin_primes)) )

    # Kato's bound
    if rank == 0 and not E.has_cm():
        L_over_Omega = E.lseries().L_ratio()
        kato_primes = Sha.bound_kato()
        primes_to_remove = []
        for p in remaining_primes:
            if p not in kato_primes:
                if verbosity > 0:
                    print 'Kato further implies that #Sha[%d] is trivial.'%p
                primes_to_remove.append(p)
            if p not in [2,3] and N%p != 0:
                if E.galois_representation().is_surjective(p):
                    if verbosity > 1:
                        print 'Kato might apply nontrivially for %d'%p
                    # ordp(sha) <= ordp(L_over_omega)
        for p in primes_to_remove:
            remaining_primes.remove(p)

    # Mazur
    if N.is_prime():
        for p in remaining_primes:
            if E.galois_representation().is_reducible(p):
                remaining_primes.remove(p)
                if verbosity > 0:
                    print 'True for p=%s by Mazur'%p

    if two_proven is False:
        remaining_primes.append(2)

    # print some extra information
    remaining_primes.sort()
    if verbosity > 1:
        if len(remaining_primes) > 0:
            print 'Remaining primes:'
        for p in remaining_primes:
            s = 'p = ' + str(p) + ': '
            if not E.galois_representation().is_reducible(p):
                s += 'ir'
            s += 'reducible, '
            if not E.galois_representation().is_surjective(p):
                s += 'not '
            s += 'surjective, '
            a_p = E.an(p)
            if E.is_good(p):
                if a_p%p != 0:
                    s += 'good ordinary'
                else:
                    s += 'good supersingular'
            else:
                assert E.is_minimal()
                if a_p == 0:
                    s += 'additive'
                elif a_p == 1:
                    s += 'split multiplicative'
                elif a_p == -1:
                    s += 'non-split multiplicative'
            if E.tamagawa_product()%p==0:
                s += ', divides a Tamagawa number'
            print s

    return remaining_primes

