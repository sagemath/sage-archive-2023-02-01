# -*- coding: utf-8 -*-
"""
This file contains a bunch of tests extracted from the published book
'Modular Forms: a Computational Approach' by William Stein, AMS 2007.

I made some minor changes below from the book. Also, in about 10
places Sage is broken which breaks examples below -- this are commented
out and trac tickets have been opened for these issues.
"""

from sage.all_cmdline import *;
import sage.plot.plot; sage.plot.plot.DOCTEST_MODE=True

def warning_function(f):
    import warnings

    def doctest_showwarning(message, category, filename, lineno, file=f):
        try:
            file.write(warnings.formatwarning(message, category, 'doctest', lineno))
        except IOError:
            pass # the file (probably stdout) is invalid
    return doctest_showwarning

def change_warning_output(file):
    import warnings
    warnings.showwarning = warning_function(file)
def example_0():	r""">>> set_random_seed(0L)

>>> change_warning_output(sys.stdout)


>>> G = SL(Integer(2),ZZ); G###line 72:_sage_    : G = SL(2,ZZ); G
Special Linear Group of degree 2 over Integer Ring
>>> S, T = G.gens()###line 74:_sage_    : S, T = G.gens()
>>> S###line 75:_sage_    : S
[ 0  1]
[-1  0]
>>> T###line 78:_sage_    : T
[1 1]
[0 1]
"""

def example_1():	r""">>> set_random_seed(0L)

>>> change_warning_output(sys.stdout)


>>> delta_qexp(Integer(6))###line 862:_sage_    : delta_qexp(6)
q - 24*q^2 + 252*q^3 - 1472*q^4 + 4830*q^5 + O(q^6)
"""

def example_2():	r""">>> set_random_seed(0L)

>>> change_warning_output(sys.stdout)


>>> bernoulli(Integer(12))###line 948:_sage_    : bernoulli(12)
-691/2730
>>> bernoulli(Integer(50))###line 950:_sage_    : bernoulli(50)
495057205241079648212477525/66
>>> len(str(bernoulli(Integer(10000))))###line 952:_sage_    : len(str(bernoulli(10000)))
27706
"""

def example_3():	r""">>> set_random_seed(0L)

>>> change_warning_output(sys.stdout)


>>> E4 = eisenstein_series_qexp(Integer(4), Integer(3))###line 1192:_sage_    : E4 = eisenstein_series_qexp(4, 3)
>>> E6 = eisenstein_series_qexp(Integer(6), Integer(3))###line 1193:_sage_    : E6 = eisenstein_series_qexp(6, 3)
>>> E4**Integer(6)###line 1194:_sage_    : E4^6
1/191102976000000 + 1/132710400000*q
                  + 203/44236800000*q^2 + O(q^3)
>>> E4**Integer(3)*E6**Integer(2)###line 1197:_sage_    : E4^3*E6^2
1/3511517184000 - 1/12192768000*q
                - 377/4064256000*q^2 + O(q^3)
>>> E6**Integer(4)###line 1200:_sage_    : E6^4
1/64524128256 - 1/32006016*q + 241/10668672*q^2 + O(q^3)
"""

def example_4():	r""">>> set_random_seed(0L)

>>> change_warning_output(sys.stdout)


>>> victor_miller_basis(Integer(28),Integer(5))###line 1306:_sage_    : victor_miller_basis(28,5)
[
1 + 15590400*q^3 + 36957286800*q^4 + O(q^5),
q + 151740*q^3 + 61032448*q^4 + O(q^5),
q^2 + 192*q^3 - 8280*q^4 + O(q^5)
]
"""

def example_5():	r""">>> set_random_seed(0L)

>>> change_warning_output(sys.stdout)


>>> R = QQ[['q']]; (q,) = R._first_ngens(Integer(1))###line 1630:_sage_    : R.<q> = QQ[['q']]
>>> F4 =  Integer(240) * eisenstein_series_qexp(Integer(4),Integer(3))###line 1631:_sage_    : F4 =  240 * eisenstein_series_qexp(4,3)
>>> F6 = -Integer(504) * eisenstein_series_qexp(Integer(6),Integer(3))###line 1632:_sage_    : F6 = -504 * eisenstein_series_qexp(6,3)
>>> F4**Integer(3)###line 1633:_sage_    : F4^3
1 + 720*q + 179280*q^2 + O(q^3)
>>> Delta = (F4**Integer(3) - F6**Integer(2))/Integer(1728); Delta###line 1635:_sage_    : Delta = (F4^3 - F6^2)/1728; Delta
q - 24*q^2 + O(q^3)
>>> F4**Integer(3) - Integer(720)*Delta###line 1637:_sage_    : F4^3 - 720*Delta
1 + 196560*q^2 + O(q^3)
"""

def example_6():	r""">>> set_random_seed(0L)

>>> change_warning_output(sys.stdout)


>>> # M = ModularForms(Integer(1),Integer(36), prec=Integer(6)).echelon_form()   # not tested -- used to work; trac #4356
>>> # M.basis()  # not tested
[
1 + 6218175600*q^4 + 15281788354560*q^5 + O(q^6),
q + 57093088*q^4 + 37927345230*q^5 + O(q^6),
q^2 + 194184*q^4 + 7442432*q^5 + O(q^6),
q^3 - 72*q^4 + 2484*q^5 + O(q^6)
]


Next we compute the matrix of the Hecke operator $T_2$.

>>> # T2 = M.hecke_matrix(Integer(2)); T2   # not tested
[34359738369   0       6218175600 9026867482214400]
[          0   0      34416831456    5681332472832]
[          0   1           194184       -197264484]
[          0   0              -72           -54528]


Finally we compute and factor its characteristic polynomial.


>>> # T2.charpoly().factor()         # not tested
(x - 34359738369) *
   (x^3 - 139656*x^2 - 59208339456*x - 1467625047588864)
"""

def example_7():	r""">>> set_random_seed(0L)

>>> change_warning_output(sys.stdout)


>>> bernoulli_mod_p(Integer(23))###line 1853:_sage_    : bernoulli_mod_p(23)
[1, 4, 13, 17, 13, 6, 10, 5, 10, 9, 15]
"""

def example_8():	r""">>> set_random_seed(0L)

>>> change_warning_output(sys.stdout)


>>> set_modsym_print_mode ('modular')###line 2476:_sage_    : set_modsym_print_mode ('modular')
>>> M = ModularSymbols(Integer(11), Integer(2))###line 2477:_sage_    : M = ModularSymbols(11, 2)
>>> M.basis()###line 2478:_sage_    : M.basis()
({Infinity,0}, {7/8,1}, {8/9,1})
>>> S = M.cuspidal_submodule()###line 2480:_sage_    : S = M.cuspidal_submodule()
>>> S.integral_basis()     # basis over ZZ.###line 2481:_sage_    : S.integral_basis()     # basis over ZZ.
({7/8,1}, {8/9,1})
>>> set_modsym_print_mode ('manin')    # set it back###line 2483:_sage_    : set_modsym_print_mode ('manin')    # set it back
"""

def example_9():	r""">>> set_random_seed(0L)

>>> change_warning_output(sys.stdout)


>>> convergents(Integer(4)/Integer(7))###line 2642:_sage_    : convergents(4/7)
[0, 1, 1/2, 4/7]
"""

def example_10():	r""">>> set_random_seed(0L)

>>> change_warning_output(sys.stdout)


>>> M = ModularSymbols(Integer(2),Integer(2))###line 2728:_sage_    : M = ModularSymbols(2,2)
>>> M###line 2729:_sage_    : M
Modular Symbols space of dimension 1 for Gamma_0(2)
of weight 2 with sign 0 over Rational Field
>>> M.manin_generators()###line 2732:_sage_    : M.manin_generators()
[(0,1), (1,0), (1,1)]

>>> M = ModularSymbols(Integer(3),Integer(2))###line 2735:_sage_    : M = ModularSymbols(3,2)
>>> M.manin_generators()###line 2736:_sage_    : M.manin_generators()
[(0,1), (1,0), (1,1), (1,2)]

>>> M = ModularSymbols(Integer(6),Integer(2))###line 2739:_sage_    : M = ModularSymbols(6,2)
>>> M.manin_generators()###line 2740:_sage_    : M.manin_generators()
[(0,1), (1,0), (1,1), (1,2), (1,3), (1,4), (1,5), (2,1),
 (2,3), (2,5), (3,1), (3,2)]
"""

def example_11():	r""">>> set_random_seed(0L)

>>> change_warning_output(sys.stdout)


>>> M = ModularSymbols(Integer(2),Integer(2))###line 2749:_sage_    : M = ModularSymbols(2,2)
>>> [x.lift_to_sl2z(Integer(2)) for x in M.manin_generators()]###line 2750:_sage_    : [x.lift_to_sl2z(2) for x in M.manin_generators()]
[[1, 0, 0, 1], [0, -1, 1, 0], [1, 0, 1, 1]]
>>> M = ModularSymbols(Integer(6),Integer(2))###line 2752:_sage_    : M = ModularSymbols(6,2)
>>> x = M.manin_generators()[Integer(9)]###line 2753:_sage_    : x = M.manin_generators()[9]
>>> x###line 2754:_sage_    : x
(2,5)
>>> x.lift_to_sl2z(Integer(6))###line 2756:_sage_    : x.lift_to_sl2z(6)
[1, 2, 2, 5]
"""

def example_12():	r""">>> set_random_seed(0L)

>>> change_warning_output(sys.stdout)


>>> M = ModularSymbols(Integer(2),Integer(2))###line 2766:_sage_    : M = ModularSymbols(2,2)
>>> M.manin_basis()###line 2767:_sage_    : M.manin_basis()
[1]
>>> [M.manin_generators()[i] for i in M.manin_basis()]###line 2769:_sage_    : [M.manin_generators()[i] for i in M.manin_basis()]
[(1,0)]
>>> M = ModularSymbols(Integer(6),Integer(2))###line 2771:_sage_    : M = ModularSymbols(6,2)
>>> M.manin_basis()###line 2772:_sage_    : M.manin_basis()
[1, 10, 11]
>>> [M.manin_generators()[i] for i in M.manin_basis()]###line 2774:_sage_    : [M.manin_generators()[i] for i in M.manin_basis()]
[(1,0), (3,1), (3,2)]


Thus, e.g., every element of $\sM_2(\Gamma_0(6))$ is a $\Q$-linear
combination of the three symbols $[(1,0)]$, $[(3,1)]$, and $[(3,2)]$.
We can write each of these as a modular symbol using the
{\tt modular\_symbol\_rep} function.

\sageindex{modular symbols printing}

>>> M.basis()###line 2786:_sage_    : M.basis()
((1,0), (3,1), (3,2))
>>> [x.modular_symbol_rep() for x in M.basis()]###line 2788:_sage_    : [x.modular_symbol_rep() for x in M.basis()]
[{Infinity,0}, {-1,-2/3}, {-1/2,-1/3}]
"""

def example_13():	r""">>> set_random_seed(0L)

>>> change_warning_output(sys.stdout)


>>> M = ModularSymbols(Integer(2),Integer(2))###line 2797:_sage_    : M = ModularSymbols(2,2)
>>> M.manin_gens_to_basis()###line 2798:_sage_    : M.manin_gens_to_basis()
[-1]
[ 1]
[ 0]
"""

def example_14():	r""">>> set_random_seed(0L)

>>> change_warning_output(sys.stdout)


>>> M = ModularSymbols(Integer(2),Integer(2))###line 2812:_sage_    : M = ModularSymbols(2,2)
>>> x = (Integer(1),Integer(0)); M(x)###line 2813:_sage_    : x = (1,0); M(x)
(1,0)
"""

def example_15():	r""">>> set_random_seed(0L)

>>> change_warning_output(sys.stdout)


>>> M = ModularSymbols(Integer(6),Integer(2))###line 2823:_sage_    : M = ModularSymbols(6,2)
>>> M.manin_gens_to_basis()###line 2824:_sage_    : M.manin_gens_to_basis()
[-1  0  0]
[ 1  0  0]
[ 0  0  0]
[ 0 -1  1]
[ 0 -1  0]
[ 0 -1  1]
[ 0  0  0]
[ 0  1 -1]
[ 0  0 -1]
[ 0  1 -1]
[ 0  1  0]
[ 0  0  1]
"""

def example_16():	r""">>> set_random_seed(0L)

>>> change_warning_output(sys.stdout)


>>> M = ModularSymbols(Integer(6),Integer(2))###line 2845:_sage_    : M = ModularSymbols(6,2)
>>> M((Integer(0),Integer(1)))###line 2846:_sage_    : M((0,1))
-(1,0)
>>> M((Integer(1),Integer(2)))###line 2848:_sage_    : M((1,2))
-(3,1) + (3,2)
"""

def example_17():	r""">>> set_random_seed(0L)

>>> change_warning_output(sys.stdout)


>>> list(HeilbronnCremona(Integer(2)))###line 2950:_sage_    : HeilbronnCremonaList(2)
[[1, 0, 0, 2], [2, 0, 0, 1], [2, 1, 0, 1], [1, 0, 1, 2]]
>>> list(HeilbronnCremona(Integer(3)))###line 2952:_sage_    : HeilbronnCremonaList(3)
[[1, 0, 0, 3], [3, 1, 0, 1], [1, 0, 1, 3], [3, 0, 0, 1],
 [3, -1, 0, 1], [-1, 0, 1, -3]]
>>> list(HeilbronnCremona(Integer(5)))###line 2955:_sage_    : HeilbronnCremonaList(5)
[[1, 0, 0, 5], [5, 2, 0, 1], [2, 1, 1, 3], [1, 0, 3, 5],
 [5, 1, 0, 1], [1, 0, 1, 5], [5, 0, 0, 1], [5, -1, 0, 1],
 [-1, 0, 1, -5], [5, -2, 0, 1], [-2, 1, 1, -3],
 [1, 0, -3, 5]]
>>> len(HeilbronnCremona(Integer(37)))###line 2960:_sage_    : len(HeilbronnCremonaList(37))
128
>>> len(HeilbronnCremona(Integer(389)))###line 2962:_sage_    : len(HeilbronnCremonaList(389))
1892
>>> len(HeilbronnCremona(Integer(2003)))###line 2964:_sage_    : len(HeilbronnCremonaList(2003))
11662
"""

def example_18():	r""">>> set_random_seed(0L)

>>> change_warning_output(sys.stdout)


>>> M = ModularSymbols(Integer(2),Integer(2))###line 2975:_sage_    : M = ModularSymbols(2,2)
>>> M.T(Integer(2)).matrix()###line 2976:_sage_    : M.T(2).matrix()
[1]
"""

def example_19():	r""">>> set_random_seed(0L)

>>> change_warning_output(sys.stdout)


>>> M = ModularSymbols(Integer(6), Integer(2))###line 2988:_sage_    : M = ModularSymbols(6, 2)
>>> M.T(Integer(2)).matrix()###line 2989:_sage_    : M.T(2).matrix()
[ 2  1 -1]
[-1  0  1]
[-1 -1  2]
>>> M.T(Integer(3)).matrix()###line 2993:_sage_    : M.T(3).matrix()
[3 2 0]
[0 1 0]
[2 2 1]
>>> M.T(Integer(3)).fcp()  # factored characteristic polynomial###line 2997:_sage_    : M.T(3).fcp()  # factored characteristic polynomial
(x - 3) * (x - 1)^2
"""

def example_20():	r""">>> set_random_seed(0L)

>>> change_warning_output(sys.stdout)


>>> M = ModularSymbols(Integer(39), Integer(2))###line 3010:_sage_    : M = ModularSymbols(39, 2)
>>> T2 = M.T(Integer(2))###line 3011:_sage_    : T2 = M.T(2)
>>> T2.matrix()###line 3012:_sage_    : T2.matrix()
[ 3  0 -1  0  0  1  1 -1  0]
[ 0  0  2  0 -1  1  0  1 -1]
[ 0  1  0 -1  1  1  0  1 -1]
[ 0  0  1  0  0  1  0  1 -1]
[ 0 -1  2  0  0  1  0  1 -1]
[ 0  0  1  1  0  1  1 -1  0]
[ 0  0  0 -1  0  1  1  2  0]
[ 0  0  0  1  0  0  2  0  1]
[ 0  0 -1  0  0  0  1  0  2]
>>> T2.fcp()     # factored characteristic polynomial###line 3022:_sage_    : T2.fcp()     # factored characteristic polynomial
(x - 1)^2 * (x - 3)^3 * (x^2 + 2*x - 1)^2


The Hecke operators commute, so their
eigenspace structures are related.


>>> T2 = M.T(Integer(2)).matrix()###line 3031:_sage_    : T2 = M.T(2).matrix()
>>> T5 = M.T(Integer(5)).matrix()###line 3032:_sage_    : T5 = M.T(5).matrix()
>>> T2*T5 - T5*T2 == Integer(0)###line 3033:_sage_    : T2*T5 - T5*T2 == 0
True
>>> T5.charpoly().factor()###line 3035:_sage_    : T5.charpoly().factor()
(x - 2)^2 * (x - 6)^3 * (x^2 - 8)^2
"""

def example_21():	r""">>> set_random_seed(0L)

>>> change_warning_output(sys.stdout)


>>> M = ModularSymbols(Integer(39), Integer(2))###line 3044:_sage_    : M = ModularSymbols(39, 2)
>>> M.T(Integer(2)).decomposition()###line 3045:_sage_    : M.T(2).decomposition()
[
Modular Symbols subspace of dimension 2 of Modular Symbols space of dimension 9 for Gamma_0(39) of weight 2 with sign 0 over Rational Field,
Modular Symbols subspace of dimension 3 of Modular Symbols space of dimension 9 for Gamma_0(39) of weight 2 with sign 0 over Rational Field,
Modular Symbols subspace of dimension 4 of Modular Symbols space of dimension 9 for Gamma_0(39) of weight 2 with sign 0 over Rational Field
]
"""

def example_22():	r""">>> set_random_seed(0L)

>>> change_warning_output(sys.stdout)


>>> M = ModularSymbols(Integer(2), Integer(2))###line 3108:_sage_    : M = ModularSymbols(2, 2)
>>> M.boundary_map()###line 3109:_sage_    : M.boundary_map()
Hecke module morphism boundary map defined by the matrix
[ 1 -1]
Domain: Modular Symbols space of dimension 1 for
Gamma_0(2) of weight ...
Codomain: Space of Boundary Modular Symbols for
Congruence Subgroup Gamma0(2) ...
>>> M.cuspidal_submodule()###line 3116:_sage_    : M.cuspidal_submodule()
Modular Symbols subspace of dimension 0 of Modular
Symbols space of dimension 1 for Gamma_0(2) of weight
2 with sign 0 over Rational Field
"""

def example_23():	r""">>> set_random_seed(0L)

>>> change_warning_output(sys.stdout)


>>> M = ModularSymbols(Integer(11), Integer(2))###line 3125:_sage_    : M = ModularSymbols(11, 2)
>>> M.boundary_map().matrix()###line 3126:_sage_    : M.boundary_map().matrix()
[ 1 -1]
[ 0  0]
[ 0  0]
>>> M.cuspidal_submodule()###line 3130:_sage_    : M.cuspidal_submodule()
Modular Symbols subspace of dimension 2 of Modular
Symbols space of dimension 3 for Gamma_0(11) of weight
2 with sign 0 over Rational Field
>>> S = M.cuspidal_submodule(); S###line 3134:_sage_    : S = M.cuspidal_submodule(); S
Modular Symbols subspace of dimension 2 of Modular
Symbols space of dimension 3 for Gamma_0(11) of weight
2 with sign 0 over Rational Field
>>> S.basis()###line 3138:_sage_    : S.basis()
((1,8), (1,9))


The following illustrates that the
Hecke operators preserve $\sS_2(\Gamma_0(N))$:

>>> S.T(Integer(2)).matrix()###line 3146:_sage_    : S.T(2).matrix()
[-2  0]
[ 0 -2]
>>> S.T(Integer(3)).matrix()###line 3149:_sage_    : S.T(3).matrix()
[-1  0]
[ 0 -1]
>>> S.T(Integer(5)).matrix()###line 3152:_sage_    : S.T(5).matrix()
[1 0]
[0 1]


A nontrivial fact is that for $p$ prime the eigenvalue of each
of these matrices is $p+1 - \#E(\F_p)$, where $E$ is the
elliptic curve $X_0(11)$ defined by the (affine) equation
$y^2 + y = x^3 - x^2 - 10x - 20.$
For example, we have

>>> E = EllipticCurve([Integer(0),-Integer(1),Integer(1),-Integer(10),-Integer(20)])###line 3164:_sage_    : E = EllipticCurve([0,-1,1,-10,-20])
>>> Integer(2) + Integer(1) - E.Np(Integer(2))###line 3165:_sage_    : 2 + 1 - E.Np(2)
-2
>>> Integer(3) + Integer(1) - E.Np(Integer(3))###line 3167:_sage_    : 3 + 1 - E.Np(3)
-1
>>> Integer(5) + Integer(1) - E.Np(Integer(5))###line 3169:_sage_    : 5 + 1 - E.Np(5)
1
>>> Integer(7) + Integer(1) - E.Np(Integer(7))###line 3171:_sage_    : 7 + 1 - E.Np(7)
-2


The same numbers appear as the eigenvalues
of Hecke operators:

>>> [S.T(p).matrix()[Integer(0),Integer(0)] for p in [Integer(2),Integer(3),Integer(5),Integer(7)]]###line 3179:_sage_    : [S.T(p).matrix()[0,0] for p in [2,3,5,7]]
[-2, -1, 1, -2]
"""

def example_24():	r""">>> set_random_seed(0L)

>>> change_warning_output(sys.stdout)


>>> M = ModularSymbols(Integer(11)); M.basis()###line 3372:_sage_    : M = ModularSymbols(11); M.basis()
((1,0), (1,8), (1,9))
>>> S = M.cuspidal_submodule(); S###line 3374:_sage_    : S = M.cuspidal_submodule(); S
Modular Symbols subspace of dimension 2 of Modular
Symbols space of dimension 3 for Gamma_0(11) of weight
2 with sign 0 over Rational Field


We compute a few Hecke operators, and then read off a
nonzero cusp form, which forms a basis for $S_2(\Gamma_0(11))$:

>>> S.T(Integer(2)).matrix()###line 3384:_sage_    : S.T(2).matrix()
[-2  0]
[ 0 -2]
>>> S.T(Integer(3)).matrix()###line 3387:_sage_    : S.T(3).matrix()
[-1  0]
[ 0 -1]
"""

def example_25():	r""">>> set_random_seed(0L)

>>> change_warning_output(sys.stdout)


>>> M = ModularSymbols(Integer(33))###line 3401:_sage_    : M = ModularSymbols(33)
>>> S = M.cuspidal_submodule(); S###line 3402:_sage_    : S = M.cuspidal_submodule(); S
Modular Symbols subspace of dimension 6 of Modular
Symbols space of dimension 9 for Gamma_0(33) of weight
2 with sign 0 over Rational Field


Thus $\dim S_2(\Gamma_0(33)) = 3$.

>>> R = PowerSeriesRing(QQ,names=('q',)); (q,) = R._first_ngens(Integer(1))###line 3411:_sage_    : R.<q> = PowerSeriesRing(QQ)
>>> v = [S.T(n).matrix()[Integer(0),Integer(0)] for n in range(Integer(1),Integer(6))]###line 3412:_sage_    : v = [S.T(n).matrix()[0,0] for n in range(1,6)]
>>> f00 = sum(v[n-Integer(1)]*q**n for n in range(Integer(1),Integer(6))) + O(q**Integer(6))###line 3413:_sage_    : f00 = sum(v[n-1]*q^n for n in range(1,6)) + O(q^6)
>>> f00###line 3414:_sage_    : f00
q - q^2 - q^3 + q^4 + O(q^6)


This gives us one basis element of $S_2(\Gamma_0(33))$.  It remains
to find two others.   We find

>>> v = [S.T(n).matrix()[Integer(0),Integer(1)] for n in range(Integer(1),Integer(6))]###line 3422:_sage_    : v = [S.T(n).matrix()[0,1] for n in range(1,6)]
>>> f01 = sum(v[n-Integer(1)]*q**n for n in range(Integer(1),Integer(6))) + O(q**Integer(6))###line 3423:_sage_    : f01 = sum(v[n-1]*q^n for n in range(1,6)) + O(q^6)
>>> f01###line 3424:_sage_    : f01
-2*q^3 + O(q^6)


and

>>> v = [S.T(n).matrix()[Integer(1),Integer(0)] for n in range(Integer(1),Integer(6))]###line 3431:_sage_    : v = [S.T(n).matrix()[1,0] for n in range(1,6)]
>>> f10 = sum(v[n-Integer(1)]*q**n for n in range(Integer(1),Integer(6))) + O(q**Integer(6))###line 3432:_sage_    : f10 = sum(v[n-1]*q^n for n in range(1,6)) + O(q^6)
>>> f10###line 3433:_sage_    : f10
q^3 + O(q^6)


This third one is (to our precision) a scalar multiple of the second,
so we look further.

>>> v = [S.T(n).matrix()[Integer(1),Integer(1)] for n in range(Integer(1),Integer(6))]###line 3441:_sage_    : v = [S.T(n).matrix()[1,1] for n in range(1,6)]
>>> f11 = sum(v[n-Integer(1)]*q**n for n in range(Integer(1),Integer(6))) + O(q**Integer(6))###line 3442:_sage_    : f11 = sum(v[n-1]*q^n for n in range(1,6)) + O(q^6)
>>> f11###line 3443:_sage_    : f11
q - 2*q^2 + 2*q^4 + q^5 + O(q^6)
"""

def example_26():	r""">>> set_random_seed(0L)

>>> change_warning_output(sys.stdout)


>>> M = ModularSymbols(Integer(23))###line 3467:_sage_    : M = ModularSymbols(23)
>>> S = M.cuspidal_submodule()###line 3468:_sage_    : S = M.cuspidal_submodule()
>>> S###line 3469:_sage_    : S
Modular Symbols subspace of dimension 4 of Modular
Symbols space of dimension 5 for Gamma_0(23) of weight
2 with sign 0 over Rational Field
>>> f = S.q_expansion_cuspforms(Integer(6))###line 3473:_sage_    : f = S.q_expansion_cuspforms(6)
>>> f(Integer(0),Integer(0))###line 3474:_sage_    : f(0,0)
q - 2/3*q^2 + 1/3*q^3 - 1/3*q^4 - 4/3*q^5 + O(q^6)
>>> f(Integer(0),Integer(1))###line 3476:_sage_    : f(0,1)
O(q^6)
>>> f(Integer(1),Integer(0))###line 3478:_sage_    : f(1,0)
-1/3*q^2 + 2/3*q^3 + 1/3*q^4 - 2/3*q^5 + O(q^6)


Thus a basis for $S_2(\Gamma_0(23))$ is
\begin{align*}
f_{0,0} &= q - \frac{2}{3}q^{2} + \frac{1}{3}q^{3} - \frac{1}{3}q^{4} - \frac{4}{3}q^{5} + \cdots, \\
f_{1,0} &= -\frac{1}{3}q^{2} + \frac{2}{3}q^{3} + \frac{1}{3}q^{4} - \frac{2}{3}q^{5} + \cdots.
\end{align*}
Or, in echelon form,
\begin{align*}
 & q - q^{3} - q^{4} + \cdots\\
 & \,\,\,\,\,q^{2} - 2q^{3} - q^{4} + 2q^{5} + \cdots
\end{align*}
which we computed using

>>> S.q_expansion_basis(Integer(6))###line 3495:_sage_    : S.q_expansion_basis(6)
[
q - q^3 - q^4 + O(q^6),
q^2 - 2*q^3 - q^4 + 2*q^5 + O(q^6)
]
"""

def example_27():	r""">>> set_random_seed(0L)

>>> change_warning_output(sys.stdout)


>>> R = Integers(Integer(49))###line 3845:_sage_    : R = Integers(49)
>>> R###line 3846:_sage_    : R
Ring of integers modulo 49


The {\tt unit\_gens} command computes the minimal generators
for $(\Z/N\Z)^*$, as defined above.

>>> R.unit_gens()###line 3854:_sage_    : R.unit_gens()
[3]
>>> Integers(Integer(25)).unit_gens()###line 3856:_sage_    : Integers(25).unit_gens()
[2]
>>> Integers(Integer(100)).unit_gens()###line 3858:_sage_    : Integers(100).unit_gens()
[51, 77]
>>> Integers(Integer(200)).unit_gens()###line 3860:_sage_    : Integers(200).unit_gens()
[151, 101, 177]
>>> Integers(Integer(2005)).unit_gens()###line 3862:_sage_    : Integers(2005).unit_gens()
[402, 1206]
>>> Integers(Integer(200000000)).unit_gens()###line 3864:_sage_    : Integers(200000000).unit_gens()
[174218751, 51562501, 187109377]
"""

def example_28():	r""">>> set_random_seed(0L)

>>> change_warning_output(sys.stdout)


>>> list(DirichletGroup(Integer(5)))###line 3914:_sage_    : list(DirichletGroup(5))
[[1], [zeta4], [-1], [-zeta4]]
>>> list(DirichletGroup(Integer(5), QQ))###line 3916:_sage_    : list(DirichletGroup(5, QQ))
[[1], [-1]]
"""

def example_29():	r""">>> set_random_seed(0L)

>>> change_warning_output(sys.stdout)


>>> G = DirichletGroup(Integer(200))###line 3962:_sage_    : G = DirichletGroup(200)
>>> G###line 3963:_sage_    : G
Group of Dirichlet characters of modulus 200 over
Cyclotomic Field of order 20 and degree 8
>>> G.exponent()###line 3966:_sage_    : G.exponent()
20
>>> G.gens()###line 3968:_sage_    : G.gens()
([-1, 1, 1], [1, -1, 1], [1, 1, zeta20])


We construct $\eps$.

>>> K = G.base_ring()###line 3975:_sage_    : K = G.base_ring()
>>> zeta = K.gen(0)###line 3976:_sage_    : zeta = K.0
>>> eps = G([Integer(1),-Integer(1),zeta**Integer(5)])###line 3977:_sage_    : eps = G([1,-1,zeta^5])
>>> eps###line 3978:_sage_    : eps
[1, -1, zeta20^5]


Finally, we evaluate $\eps$ at $3$.

>>> eps(Integer(3))###line 3985:_sage_    : eps(3)
zeta20^5
>>> -zeta**Integer(15)###line 3987:_sage_    : -zeta^15
zeta20^5
"""

def example_30():	r""">>> set_random_seed(0L)

>>> change_warning_output(sys.stdout)


>>> kronecker(Integer(151),Integer(200))###line 4365:_sage_    : kronecker(151,200)
1
>>> kronecker(Integer(101),Integer(200))###line 4367:_sage_    : kronecker(101,200)
-1
>>> kronecker(Integer(177),Integer(200))###line 4369:_sage_    : kronecker(177,200)
1
"""

def example_31():	r""">>> set_random_seed(0L)

>>> change_warning_output(sys.stdout)


>>> G = DirichletGroup(Integer(20))###line 4556:_sage_    : G = DirichletGroup(20)
>>> G.galois_orbits()###line 4557:_sage_    : G.galois_orbits()
[
[[1, 1]],
[[1, zeta4], [1, -zeta4]],
[[1, -1]],
[[-1, 1]],
[[-1, zeta4], [-1, -zeta4]],
[[-1, -1]]
]
"""

def example_32():	r""">>> set_random_seed(0L)

>>> change_warning_output(sys.stdout)


>>> G = DirichletGroup(Integer(11), QQ); G###line 4629:_sage_    : G = DirichletGroup(11, QQ); G
Group of Dirichlet characters of modulus 11 over
Rational Field


A Dirichlet character prints as a matrix that gives
the values of the character on canonical generators
of $(\Z/N\Z)^*$ (as discussed below).

>>> list(G)###line 4639:_sage_    : list(G)
[[1], [-1]]
>>> eps = G.gen(0)      # 0th generator for Dirichlet group###line 4641:_sage_    : eps = G.0      # 0th generator for Dirichlet group
>>> eps###line 4642:_sage_    : eps
[-1]


The character $\varepsilon$ takes the value $-1$ on the
unit generator.

>>> G.unit_gens()###line 4650:_sage_    : G.unit_gens()
[2]
>>> eps(Integer(2))###line 4652:_sage_    : eps(2)
-1
>>> eps(Integer(3))###line 4654:_sage_    : eps(3)
1


It is $0$ on any integer not coprime to $11$:

>>> [eps(Integer(11)*n) for n in range(Integer(10))]###line 4661:_sage_    : [eps(11*n) for n in range(10)]
[0, 0, 0, 0, 0, 0, 0, 0, 0, 0]
"""

def example_33():	r""">>> set_random_seed(0L)

>>> change_warning_output(sys.stdout)


>>> R = CyclotomicField(Integer(4))###line 4669:_sage_    : R = CyclotomicField(4)
>>> CyclotomicField(Integer(4))###line 4670:_sage_    : CyclotomicField(4)
Cyclotomic Field of order 4 and degree 2


Then we define $G = D(15,\Q({\zeta_4}))$.

>>> G = DirichletGroup(Integer(15), R)###line 4677:_sage_    : G = DirichletGroup(15, R)
>>> G###line 4678:_sage_    : G
Group of Dirichlet characters of modulus 15 over
Cyclotomic Field of order 4 and degree 2


Next we list each of its elements.

>>> list(G)###line 4686:_sage_    : list(G)
[[1, 1], [-1, 1], [1, zeta4], [-1, zeta4], [1, -1],
[-1, -1], [1, -zeta4], [-1, -zeta4]]


Now we evaluate the second generator of $G$ on various
integers:

>>> e = G.gen(1)###line 4695:_sage_    : e = G.1
>>> e(Integer(4))###line 4696:_sage_    : e(4)
-1
>>> e(-Integer(1))###line 4698:_sage_    : e(-1)
-1
>>> e(Integer(5))###line 4700:_sage_    : e(5)
0


Finally we list all the values of $e$.

>>> [e(n) for n in range(Integer(15))]###line 4707:_sage_    : [e(n) for n in range(15)]
[0, 1, zeta4, 0, -1, 0, 0, zeta4, -zeta4,
    0, 0, 1, 0, -zeta4, -1]
"""

def example_34():	r""">>> set_random_seed(0L)

>>> change_warning_output(sys.stdout)


>>> G = DirichletGroup(Integer(15), GF(Integer(5))); G###line 4715:_sage_    : G = DirichletGroup(15, GF(5)); G
Group of Dirichlet characters of modulus 15
      over Finite Field of size 5


We list all the elements of $G$, again represented by lists
that give the images of each unit generator, as an element of $\F_5$.

>>> list(G)###line 4724:_sage_    : list(G)
  [[1, 1], [4, 1], [1, 2], [4, 2], [1, 4], [4, 4],
   [1, 3], [4, 3]]


We evaluate the second generator of $G$ on several integers.

>>> e = G.gen(1)###line 4732:_sage_    : e = G.1
>>> e(-Integer(1))###line 4733:_sage_    : e(-1)
4
>>> e(Integer(2))###line 4735:_sage_    : e(2)
2
>>> e(Integer(5))###line 4737:_sage_    : e(5)
0
>>> print [e(n) for n in range(Integer(15))]###line 4739:_sage_    : print [e(n) for n in range(15)]
[0, 1, 2, 0, 4, 0, 0, 2, 3, 0, 0, 1, 0, 3, 4]
"""

def example_35():	r""">>> set_random_seed(0L)

>>> change_warning_output(sys.stdout)


>>> G = DirichletGroup(Integer(5))###line 4961:_sage_    : G = DirichletGroup(5)
>>> e = G.gen(0)###line 4962:_sage_    : e = G.0
>>> e(Integer(2))###line 4963:_sage_    : e(2)
zeta4


We compute the Bernoulli number $B_{1,\eps}$.

>>> e.bernoulli(Integer(1))###line 4970:_sage_    : e.bernoulli(1)
-1/5*zeta4 - 3/5


We compute $B_{9,\eps}$.

>>> e.bernoulli(Integer(9))###line 4977:_sage_    : e.bernoulli(9)
-108846/5*zeta4 - 176868/5
"""

def example_36():	r""">>> set_random_seed(0L)

>>> change_warning_output(sys.stdout)


>>> E = EisensteinForms(Gamma1(Integer(13)),Integer(2))###line 5205:_sage_    : E = EisensteinForms(Gamma1(13),2)
>>> E.eisenstein_series()###line 5206:_sage_    : E.eisenstein_series()
[
1/2 + q + 3*q^2 + 4*q^3 + 7*q^4 + 6*q^5 + O(q^6),
-7/13*zeta6 - 11/13 + q + (2*zeta6 + 1)*q^2 + (-3*zeta6 + 1)*q^3 + (6*zeta6 - 3)*q^4 - 4*q^5 + O(q^6),
q + (zeta6 + 2)*q^2 + (-zeta6 + 3)*q^3 + (3*zeta6 + 3)*q^4 + 4*q^5 + O(q^6),
-zeta6 + q + (2*zeta6 - 1)*q^2 + (3*zeta6 - 2)*q^3 + (-2*zeta6 - 1)*q^4 + 6*q^5 + O(q^6),
q + (zeta6 + 1)*q^2 + (zeta6 + 2)*q^3 + (zeta6 + 2)*q^4 + 6*q^5 + O(q^6),
-1 + q - q^2 + 4*q^3 + 3*q^4 - 4*q^5 + O(q^6),
q + q^2 + 4*q^3 + 3*q^4 + 4*q^5 + O(q^6),
zeta6 - 1 + q + (-2*zeta6 + 1)*q^2 + (-3*zeta6 + 1)*q^3 + (2*zeta6 - 3)*q^4 + 6*q^5 + O(q^6),
q + (-zeta6 + 2)*q^2 + (-zeta6 + 3)*q^3 + (-zeta6 + 3)*q^4 + 6*q^5 + O(q^6),
7/13*zeta6 - 18/13 + q + (-2*zeta6 + 3)*q^2 + (3*zeta6 - 2)*q^3 + (-6*zeta6 + 3)*q^4 - 4*q^5 + O(q^6),
q + (-zeta6 + 3)*q^2 + (zeta6 + 2)*q^3 + (-3*zeta6 + 6)*q^4 + 4*q^5 + O(q^6)
]


We can also compute the parameters $\chi,\psi,t$ that
define each series:

>>> e = E.eisenstein_series()###line 5213:_sage_    : e = E.eisenstein_series()
>>> for e in E.eisenstein_series():###line 5214:_sage_    : for e in E.eisenstein_series():
...       print e.parameters()
...
([1], [1], 13)
([1], [zeta6], 1)
([zeta6], [1], 1)
([1], [zeta6 - 1], 1)
([zeta6 - 1], [1], 1)
([1], [-1], 1)
([-1], [1], 1)
([1], [-zeta6], 1)
([-zeta6], [1], 1)
([1], [-zeta6 + 1], 1)
([-zeta6 + 1], [1], 1)
"""

def example_37():	r""">>> set_random_seed(0L)

>>> change_warning_output(sys.stdout)


>>> dimension_cusp_forms(Gamma0(Integer(2007)),Integer(2))###line 5413:_sage_    : dimension_cusp_forms(Gamma0(2007),2)
221
>>> dimension_eis(Gamma0(Integer(2007)),Integer(2))###line 5415:_sage_    : dimension_eis(Gamma0(2007),2)
7
>>> dimension_modular_forms(Gamma0(Integer(2007)),Integer(2))###line 5417:_sage_    : dimension_modular_forms(Gamma0(2007),2)
228
"""

def example_38():	r""">>> set_random_seed(0L)

>>> change_warning_output(sys.stdout)


>>> dimension_new_cusp_forms(Gamma0(Integer(11)),Integer(12))###line 5494:_sage_    : dimension_new_cusp_forms(Gamma0(11),12)
8
>>> dimension_cusp_forms(Gamma0(Integer(11)),Integer(12))###line 5496:_sage_    : dimension_cusp_forms(Gamma0(11),12)
10
>>> dimension_new_cusp_forms(Gamma0(Integer(2007)),Integer(12))###line 5498:_sage_    : dimension_new_cusp_forms(Gamma0(2007),12)
1017
>>> dimension_cusp_forms(Gamma0(Integer(2007)),Integer(12))###line 5500:_sage_    : dimension_cusp_forms(Gamma0(2007),12)
2460
"""

def example_39():	r""">>> set_random_seed(0L)

>>> change_warning_output(sys.stdout)


>>> dimension_cusp_forms(Gamma1(Integer(2007)),Integer(2))###line 5617:_sage_    : dimension_cusp_forms(Gamma1(2007),2)
147409
>>> dimension_eis(Gamma1(Integer(2007)),Integer(2))###line 5619:_sage_    : dimension_eis(Gamma1(2007),2)
3551
>>> dimension_modular_forms(Gamma1(Integer(2007)),Integer(2))###line 5621:_sage_    : dimension_modular_forms(Gamma1(2007),2)
150960
"""

def example_40():	r""">>> set_random_seed(0L)

>>> change_warning_output(sys.stdout)


>>> G = DirichletGroup(Integer(2007))###line 5808:_sage_    : G = DirichletGroup(2007)
>>> e = prod(G.gens(), G(Integer(1)))###line 5809:_sage_    : e = prod(G.gens(), G(1))


Next we compute the dimension of the four spaces.

>>> dimension_cusp_forms(e,Integer(2))###line 5815:_sage_    : dimension_cusp_forms(e,2)
222
>>> dimension_cusp_forms(e,Integer(3))###line 5817:_sage_    : dimension_cusp_forms(e,3)
0
>>> dimension_cusp_forms(e,Integer(4))###line 5819:_sage_    : dimension_cusp_forms(e,4)
670
>>> dimension_cusp_forms(e,Integer(24))###line 5821:_sage_    : dimension_cusp_forms(e,24)
5150


We can also compute dimensions of the corresponding
spaces of Eisenstein series.

>>> dimension_eis(e,Integer(2))###line 5829:_sage_    : dimension_eis(e,2)
4
>>> dimension_eis(e,Integer(3))###line 5831:_sage_    : dimension_eis(e,3)
0
>>> dimension_eis(e,Integer(4))###line 5833:_sage_    : dimension_eis(e,4)
4
>>> dimension_eis(e,Integer(24))###line 5835:_sage_    : dimension_eis(e,24)
4
"""

def example_41():	r""">>> set_random_seed(0L)

>>> change_warning_output(sys.stdout)


>>> G = DirichletGroup(Integer(2007), QQ)###line 5874:_sage_    : G = DirichletGroup(2007, QQ)
>>> e = prod(G.gens(), G(Integer(1)))###line 5875:_sage_    : e = prod(G.gens(), G(1))
>>> dimension_new_cusp_forms(e,Integer(2))###line 5876:_sage_    : dimension_new_cusp_forms(e,2)
76
"""

def example_42():	r""">>> set_random_seed(0L)

>>> change_warning_output(sys.stdout)


>>> p = Integer(389)###line 6131:_sage_    : p = 389
>>> k = GF(p)###line 6132:_sage_    : k = GF(p)
>>> a = k(Integer(7)/Integer(13)); a###line 6133:_sage_    : a = k(7/13); a
210
>>> a.rational_reconstruction()###line 6135:_sage_    : a.rational_reconstruction()
7/13
"""

def example_43():	r""">>> set_random_seed(0L)

>>> change_warning_output(sys.stdout)


>>> M = MatrixSpace(QQ,Integer(4),Integer(8))###line 6513:_sage_    : M = MatrixSpace(QQ,4,8)
>>> A = M([[-9,6,7,3,1,0,0,0],[-10,3,8,2,0,1,0,0],[3,-6,2,8,0,0,1,0],[-8,-6,-8,6,0,0,0,1]])


%end_sage
First choose the ``random'' prime $p=41$, which does not
divide any of the entries of $A$, and compute the
echelon form of the reduction of $A$ modulo $41$.\sageindex{echelon form}
%begin_sage

>>> A41 = MatrixSpace(GF(Integer(41)),Integer(4),Integer(8))(A)###line 6525:_sage_    : A41 = MatrixSpace(GF(41),4,8)(A)
>>> E41 = A41.echelon_form()###line 6526:_sage_    : E41 = A41.echelon_form()


%end_sage
The echelon form of $A\pmod{41}$ is
$$
\left(\begin{array}{rrrrrrrr}
1&0&0&2&0&20&33&18\\
0&1&0&40&0&30&7&1\\
0&0&1&39&0&19&13&17\\
0&0&0&0&1&31&0&37
\end{array}\right).
$$
Thus we take $c_0=0$, $c_1=1$, $c_2=2$, and $c_3 = 4$.
Also $r_i=i$ for $i=0,1,2,3$.
%begin_sage
Next extract the submatrix $B$.

>>> B = A.matrix_from_columns([Integer(0),Integer(1),Integer(2),Integer(4)])###line 6545:_sage_    : B = A.matrix_from_columns([0,1,2,4])


%end_sage
The submatrix $B$ is
$$
  B = \left(\begin{array}{rrrr}
-9&6&7&1\\
-10&3&8&0\\
3&-6&2&0\\
-8&-6&-8&0
\end{array}\right).
$$
The inverse of $B$ is
$$
 B^{-1} = \left(\begin{array}{rrrr}
0&-\frac{5}{92}&\frac{1}{46}&-\frac{9}{184}\vspace{2pt}\\
0&-\frac{1}{138}&-\frac{3}{23}&-\frac{11}{276}\vspace{2pt}\\
0&\frac{11}{184}&\frac{7}{92}&-\frac{17}{368}\vspace{2pt}\\
1&-\frac{159}{184}&\frac{41}{92}&\frac{45}{368}
\end{array}\right).
$$
Multiplying by $A$ yields
$$
E = B^{-1} A =
\left(\begin{array}{rrrrrrrr}
1&0&0&-\frac{21}{92}&0&-\frac{5}{92}&\frac{1}{46}&-\frac{9}{184}\vspace{2pt}\\
0&1&0&-\frac{179}{138}&0&-\frac{1}{138}&-\frac{3}{23}&-\frac{11}{276}\vspace{2pt}\\
0&0&1&\frac{83}{184}&0&\frac{11}{184}&\frac{7}{92}&-\frac{17}{368}\vspace{2pt}\\
0&0&0&\frac{1025}{184}&1&-\frac{159}{184}&\frac{41}{92}&\frac{45}{368}
\end{array}\right).
$$
%begin_sage

>>> E = B**(-Integer(1))*A###line 6580:_sage_    : E = B^(-1)*A


%end_sage
This is {\em not} the echelon form of $A$.  Indeed, it is not
even in echelon form, since the last row is not
normalized so the leftmost nonzero entry is $1$.
We thus choose another random prime, say $p=43$.
The echelon form mod $43$ has columns $0,1,2,3$ as pivot columns.
We thus extract the matrix
$$
 B = \left(\begin{array}{rrrr}
-9&6&7&3\\
-10&3&8&2\\
3&-6&2&8\\
-8&-6&-8&6
\end{array}\right).
$$
%begin_sage

>>> B = A.matrix_from_columns([Integer(0),Integer(1),Integer(2),Integer(3)])###line 6601:_sage_    : B = A.matrix_from_columns([0,1,2,3])
"""

def example_44():	r""">>> set_random_seed(0L)

>>> change_warning_output(sys.stdout)


>>> M = ModularSymbols(Gamma0(Integer(6))); M###line 8524:_sage_    : M = ModularSymbols(Gamma0(6)); M
Modular Symbols space of dimension 3 for Gamma_0(6)
of weight 2 with sign 0 over Rational Field
>>> M.new_subspace()###line 8527:_sage_    : M.new_subspace()
Modular Symbols subspace of dimension 1 of Modular
Symbols space of dimension 3 for Gamma_0(6) of weight
2 with sign 0 over Rational Field
>>> M.old_subspace()###line 8531:_sage_    : M.old_subspace()
Modular Symbols subspace of dimension 3 of Modular
Symbols space of dimension 3 for Gamma_0(6) of weight
2 with sign 0 over Rational Field
"""

def example_45():	r""">>> set_random_seed(0L)

>>> change_warning_output(sys.stdout)


>>> G = DirichletGroup(Integer(13))###line 9526:_sage_    : G = DirichletGroup(13)
>>> G###line 9527:_sage_    : G
Group of Dirichlet characters of modulus 13 over
Cyclotomic Field of order 12 and degree 4
>>> dimension_modular_forms(Gamma1(Integer(13)),Integer(2))###line 9530:_sage_    : dimension_modular_forms(Gamma1(13),2)
13
>>> [dimension_modular_forms(e,Integer(2)) for e in G]###line 9532:_sage_    : [dimension_modular_forms(e,2) for e in G]
[1, 0, 3, 0, 2, 0, 2, 0, 2, 0, 3, 0]
"""

def example_46():	r""">>> set_random_seed(0L)

>>> change_warning_output(sys.stdout)


>>> G = DirichletGroup(Integer(100))###line 9537:_sage_    : G = DirichletGroup(100)
>>> G###line 9538:_sage_    : G
Group of Dirichlet characters of modulus 100 over
Cyclotomic Field of order 20 and degree 8
>>> dimension_modular_forms(Gamma1(Integer(100)),Integer(2))###line 9541:_sage_    : dimension_modular_forms(Gamma1(100),2)
370
>>> v = [dimension_modular_forms(e,Integer(2)) for e in G]; v###line 9543:_sage_    : v = [dimension_modular_forms(e,2) for e in G]; v
[24, 0, 0, 17, 18, 0, 0, 17, 18, 0, 0, 21, 18, 0, 0, 17,
  18, 0, 0, 17, 24, 0, 0, 17, 18, 0, 0, 17, 18, 0, 0, 21,
  18, 0, 0, 17, 18, 0, 0, 17]
>>> sum(v)###line 9547:_sage_    : sum(v)
370
"""

def example_47():	r""">>> set_random_seed(0L)

>>> change_warning_output(sys.stdout)


>>> S = CuspForms(Gamma0(Integer(45)), Integer(2), prec=Integer(14)); S###line 9655:_sage_    : S = CuspForms(Gamma0(45), 2, prec=14); S
Cuspidal subspace of dimension 3 of Modular Forms space
of dimension 10 for Congruence Subgroup Gamma0(45) of
weight 2 over Rational Field
>>> S.basis()###line 9659:_sage_    : S.basis()
[
q - q^4 - q^10 - 2*q^13 + O(q^14),
q^2 - q^5 - 3*q^8 + 4*q^11 + O(q^14),
q^3 - q^6 - q^9 - q^12 + O(q^14)
]
>>> #  (BUG FIX THIS!  trac 4357)
>>> # S.new_subspace().basis()   # not tested
(q - q^4 - q^10 - 2*q^13 + O(q^14),)
>>> CuspForms(Gamma0(Integer(9)),Integer(2))###line 9667:_sage_    : CuspForms(Gamma0(9),2)
Cuspidal subspace of dimension 0 of Modular Forms space
of dimension 3 for Congruence Subgroup Gamma0(9) of
weight 2 over Rational Field
>>> CuspForms(Gamma0(Integer(15)),Integer(2), prec=Integer(10)).basis()###line 9671:_sage_    : CuspForms(Gamma0(15),2, prec=10).basis()
[
q - q^2 - q^3 - q^4 + q^5 + q^6 + 3*q^8 + q^9 + O(q^10)
]
"""


if __name__ ==  '__main__':
    import doctest, sys
    s = doctest.testmod(sys.modules[__name__],
                   optionflags=doctest.NORMALIZE_WHITESPACE
                              |doctest.ELLIPSIS,
                   verbose=False,
                   globs=globals())
    quit_sage(verbose=False)
    sys.exit(s[0])
