r"""
This file contains doctests for the Chapter "k-Schur function primer"
for the book "k-Schur functions and affine Schubert calculus" for code
written by Anne Schilling and Mike Zabrocki, 2012
"""

"""
Sage example in ./kschurnotes/notes-mike-anne.tex, line 178::

    sage: P = Partitions(4); P
    Partitions of the integer 4
    sage: P.list()
    [[4], [3, 1], [2, 2], [2, 1, 1], [1, 1, 1, 1]]

Sage example in ./kschurnotes/notes-mike-anne.tex, line 185::

    sage: la=Partition([2,2]); mu=Partition([3,1])
    sage: mu.dominates(la)
    True

Sage example in ./kschurnotes/notes-mike-anne.tex, line 191::

    sage: ord = lambda x,y: y.dominates(x)
    sage: P = Poset([Partitions(6), ord], facade=True)
    sage: H = P.hasse_diagram()

Sage example in ./kschurnotes/notes-mike-anne.tex, line 202::

    sage: la=Partition([4,3,3,3,2,2,1])
    sage: la.conjugate()
    [7, 6, 4, 1]
    sage: la.k_split(4)
    [[4], [3, 3], [3, 2], [2, 1]]

Sage example in ./kschurnotes/notes-mike-anne.tex, line 210::

    sage: p = SkewPartition([[2,1],[1]])
    sage: p.is_connected()
    False

Sage example in ./kschurnotes/notes-mike-anne.tex, line 484::

    sage: la = Partition([4,3,3,3,2,2,1])
    sage: kappa = la.k_skew(4); kappa
    [[12, 8, 5, 5, 2, 2, 1], [8, 5, 2, 2]]

Sage example in ./kschurnotes/notes-mike-anne.tex, line 490::

    sage: kappa.row_lengths()
    [4, 3, 3, 3, 2, 2, 1]

Sage example in ./kschurnotes/notes-mike-anne.tex, line 495::

    sage: tau = Core([12,8,5,5,2,2,1],5)
    sage: mu = tau.to_bounded_partition(); mu
    [4, 3, 3, 3, 2, 2, 1]
    sage: mu.to_core(4)
    [12, 8, 5, 5, 2, 2, 1]

Sage example in ./kschurnotes/notes-mike-anne.tex, line 503::

    sage: Cores(3,6).list()
    [[6, 4, 2], [5, 3, 1, 1], [4, 2, 2, 1, 1], [3, 3, 2, 2, 1, 1]]

Sage example in ./kschurnotes/notes-mike-anne.tex, line 572::

    sage: W = WeylGroup(['A',4,1])
    sage: S = W.simple_reflections()
    sage: [s.reduced_word() for s in S]
    [[0], [1], [2], [3], [4]]

Sage example in ./kschurnotes/notes-mike-anne.tex, line 580::

    sage: w = W.an_element(); w
    [ 2  0  0  1 -2]
    [ 2  0  0  0 -1]
    [ 1  1  0  0 -1]
    [ 1  0  1  0 -1]
    [ 1  0  0  1 -1]
    sage: w.reduced_word()
    [0, 1, 2, 3, 4]
    sage: w = W.from_reduced_word([2,1,0])
    sage: w.is_affine_grassmannian()
    True

Sage example in ./kschurnotes/notes-mike-anne.tex, line 651::

    sage: c = Core([7,3,1],5)
    sage: c.affine_symmetric_group_simple_action(2)
    [8, 4, 1, 1]
    sage: c.affine_symmetric_group_simple_action(0)
    [7, 3, 1]

Sage example in ./kschurnotes/notes-mike-anne.tex, line 661::

    sage: k=4; length=3
    sage: W = WeylGroup(['A',k,1])
    sage: G = W.affine_grassmannian_elements_of_given_length(length)
    sage: [w.reduced_word() for w in G]
    [[2, 1, 0], [4, 1, 0], [3, 4, 0]]

    sage: C = Cores(k+1,length)
    sage: [c.to_grassmannian().reduced_word() for c in C]
    [[2, 1, 0], [4, 1, 0], [3, 4, 0]]

Sage example in ./kschurnotes/notes-mike-anne.tex, line 725::

    sage: la = Partition([4,3,3,3,2,2,1])
    sage: c = la.to_core(4); c
    [12, 8, 5, 5, 2, 2, 1]
    sage: W = WeylGroup(['A',4,1])
    sage: w = W.from_reduced_word([4,1,0,2,1,4,3,2,0,4,3,1,0,4,3,2,1,0])
    sage: c.to_grassmannian() == w
    True

Sage example in ./kschurnotes/notes-mike-anne.tex, line 769::

    sage: la = Partition([4,3,3,3,2,2,1])
    sage: la.k_conjugate(4)
    [3, 2, 2, 2, 2, 1, 1, 1, 1, 1, 1, 1]

Sage example in ./kschurnotes/notes-mike-anne.tex, line 1080::

    sage: c = Core([3,1,1],3)
    sage: c.weak_covers()
    [[4, 2, 1, 1]]
    sage: c.strong_covers()
    [[5, 3, 1], [4, 2, 1, 1], [3, 2, 2, 1, 1]]

Sage example in ./kschurnotes/notes-mike-anne.tex, line 1088::

    sage: kappa = Core([4,1],4)
    sage: tau = Core([2,1],4)
    sage: tau.weak_le(kappa)
    False
    sage: tau.strong_le(kappa)
    True

Sage example in ./kschurnotes/notes-mike-anne.tex, line 1097::

    sage: C = sum(([c for c in Cores(4,m)] for m in range(7)),[])
    sage: ord = lambda x,y: x.weak_le(y)
    sage: P = Poset([C, ord], cover_relations = False)
    sage: H = P.hasse_diagram()
    sage: view(H)        #optional

Sage example in ./kschurnotes/notes-mike-anne.tex, line 1262::

    sage: Sym = SymmetricFunctions(QQ)
    sage: h = Sym.homogeneous()
    sage: m = Sym.monomial()

Sage example in ./kschurnotes/notes-mike-anne.tex, line 1268::

    sage: f = h[3,1]+h[2,2]
    sage: m(f)
    10*m[1, 1, 1, 1] + 7*m[2, 1, 1] + 5*m[2, 2] + 4*m[3, 1] + 2*m[4]

Sage example in ./kschurnotes/notes-mike-anne.tex, line 1275::

    sage: f.scalar(h[2,1,1])
    7
    sage: m(f).coefficient([2,1,1])
    7

Sage example in ./kschurnotes/notes-mike-anne.tex, line 1283::

    sage: p = Sym.power()
    sage: e = Sym.elementary()
    sage: sum( (-1)**(i-1)*e[4-i]*p[i] for i in range(1,4) ) - p[4]
    4*e[4]
    sage: sum( (-1)**(i-1)*p[i]*e[4-i] for i in range(1,4) ) - p[4]
    1/6*p[1, 1, 1, 1] - p[2, 1, 1] + 1/2*p[2, 2] + 4/3*p[3, 1] - p[4]

Sage example in ./kschurnotes/notes-mike-anne.tex, line 1366::

    sage: Sym = SymmetricFunctions(QQ)
    sage: s = Sym.schur()
    sage: m = Sym.monomial()
    sage: h = Sym.homogeneous()
    sage: m(s[1,1,1])
    m[1, 1, 1]
    sage: h(s[1,1,1])
    h[1, 1, 1] - 2*h[2, 1] + h[3]

Sage example in ./kschurnotes/notes-mike-anne.tex, line 1377::

    sage: p = Sym.power()
    sage: s = Sym.schur()
    sage: p(s[1,1,1])
    1/6*p[1, 1, 1] - 1/2*p[2, 1] + 1/3*p[3]
    sage: p(s[2,1])
    1/3*p[1, 1, 1] - 1/3*p[3]
    sage: p(s[3])
    1/6*p[1, 1, 1] + 1/2*p[2, 1] + 1/3*p[3]

Sage example in ./kschurnotes/notes-mike-anne.tex, line 1388::

    sage: s[2,1].scalar(s[1,1,1])
    0
    sage: s[2,1].scalar(s[2,1])
    1

Sage example in ./kschurnotes/notes-mike-anne.tex, line 1571::

  sage: Sym = SymmetricFunctions(FractionField(QQ['t']))
  sage: Qp = Sym.hall_littlewood().Qp()
  sage: Qp.base_ring()
  Fraction Field of Univariate Polynomial Ring in t over Rational Field
  sage: s = Sym.schur()
  sage: s(Qp[1,1,1])
  s[1, 1, 1] + (t^2+t)*s[2, 1] + t^3*s[3]

Sage example in ./kschurnotes/notes-mike-anne.tex, line 1584::

  sage: t = Qp.t
  sage: s[2,1].scalar(s[3].theta_qt(t,0))
  t^2 - t

Sage example in ./kschurnotes/notes-mike-anne.tex, line 1590::

  sage: s(Qp([1,1])).hl_creation_operator([3])
  s[3, 1, 1] + t*s[3, 2] + (t^2+t)*s[4, 1] + t^3*s[5]
  sage: s(Qp([3,1,1]))
  s[3, 1, 1] + t*s[3, 2] + (t^2+t)*s[4, 1] + t^3*s[5]

Sage example in ./kschurnotes/notes-mike-anne.tex, line 1619::

  sage: Sym = SymmetricFunctions(FractionField(QQ['q,t']))
  sage: Mac = Sym.macdonald()
  sage: H = Mac.H()
  sage: s = Sym.schur()
  sage: for la in Partitions(3):
  ...     print "H", la, "=", s(H(la))
  H [3] = q^3*s[1, 1, 1] + (q^2+q)*s[2, 1] + s[3]
  H [2, 1] = q*s[1, 1, 1] + (q*t+1)*s[2, 1] + t*s[3]
  H [1, 1, 1] = s[1, 1, 1] + (t^2+t)*s[2, 1] + t^3*s[3]

Sage example in ./kschurnotes/notes-mike-anne.tex, line 1632::

  sage: Sym = SymmetricFunctions(FractionField(QQ['t']))
  sage: Mac = Sym.macdonald(q=0)
  sage: H = Mac.H()
  sage: s = Sym.schur()
  sage: for la in Partitions(3):
  ...    print "H",la, "=", s(H(la))
  H [3] = s[3]
  H [2, 1] = s[2, 1] + t*s[3]
  H [1, 1, 1] = s[1, 1, 1] + (t^2+t)*s[2, 1] + t^3*s[3]
  sage: Qp = Sym.hall_littlewood().Qp()
  sage: s(Qp[1, 1, 1])
  s[1, 1, 1] + (t^2+t)*s[2, 1] + t^3*s[3]

Sage example in ./kschurnotes/notes-mike-anne.tex, line 1647::

  sage: Sym = SymmetricFunctions(FractionField(QQ['q']))
  sage: Mac = Sym.macdonald(t=0)
  sage: H = Mac.H()
  sage: s = Sym.schur()
  sage: for la in Partitions(3):
  ...    print "H",la, "=", s(H(la))
  H [3] = q^3*s[1, 1, 1] + (q^2+q)*s[2, 1] + s[3]
  H [2, 1] = q*s[1, 1, 1] + s[2, 1]
  H [1, 1, 1] = s[1, 1, 1]

Sage example in ./kschurnotes/notes-mike-anne.tex, line 1852::

    sage: SemistandardTableaux([5,2],[4,2,1]).list()
    [[[1, 1, 1, 1, 2], [2, 3]], [[1, 1, 1, 1, 3], [2, 2]]]

Sage example in ./kschurnotes/notes-mike-anne.tex, line 1857::

    sage: P = Partitions(4)
    sage: P.list()
    [[4], [3, 1], [2, 2], [2, 1, 1], [1, 1, 1, 1]]
    sage: n = P.cardinality(); n
    5
    sage: K = matrix(QQ,n,n,
    ...           [[SemistandardTableaux(la,mu).cardinality()
    ...            for mu in P] for la in P])
    sage: K
    [1 1 1 1 1]
    [0 1 1 2 3]
    [0 0 1 1 2]
    [0 0 0 1 3]
    [0 0 0 0 1]

Sage example in ./kschurnotes/notes-mike-anne.tex, line 1966::

    sage: Sym = SymmetricFunctions(QQ)
    sage: ks4 = Sym.kschur(4,t=1)
    sage: h = Sym.homogeneous()
    sage: ks4(h[4,3,1])
    ks4[4, 3, 1] + ks4[4, 4]

    sage: ks6 = Sym.kschur(6,t=1)
    sage: ks6(h[4,3,1])
    ks6[4, 3, 1] + ks6[4, 4] + ks6[5, 2, 1] + 2*ks6[5, 3]
      + ks6[6, 1, 1] + ks6[6, 2]

Sage example in ./kschurnotes/notes-mike-anne.tex, line 1980::

    sage: Sym = SymmetricFunctions(QQ)
    sage: ks = Sym.kschur(3,t=1)
    sage: ks.realization_of()
    3-bounded Symmetric Functions over Rational Field with t=1
    sage: s = Sym.schur()
    sage: s.realization_of()
    Symmetric Functions over Rational Field

Sage example in ./kschurnotes/notes-mike-anne.tex, line 2077::

    sage: Sym = SymmetricFunctions(QQ)
    sage: ks = Sym.kschur(3,t=1)
    sage: h = Sym.homogeneous()
    sage: for mu in Partitions(7, max_part =3):
    ...     print h(ks(mu))
    ...
    h[3, 3, 1]
    h[3, 2, 2] - h[3, 3, 1]
    h[3, 2, 1, 1] - h[3, 2, 2]
    h[3, 1, 1, 1, 1] - 2*h[3, 2, 1, 1] + h[3, 3, 1]
    h[2, 2, 2, 1] - h[3, 2, 1, 1] - h[3, 2, 2] + h[3, 3, 1]
    h[2, 2, 1, 1, 1] - 2*h[2, 2, 2, 1] - h[3, 1, 1, 1, 1]
        + 2*h[3, 2, 1, 1] + h[3, 2, 2] - h[3, 3, 1]
    h[2, 1, 1, 1, 1, 1] - 3*h[2, 2, 1, 1, 1] + 2*h[2, 2, 2, 1]
        + h[3, 2, 1, 1] - h[3, 2, 2]
    h[1, 1, 1, 1, 1, 1, 1] - 4*h[2, 1, 1, 1, 1, 1] + 4*h[2, 2, 1, 1, 1]
        + 2*h[3, 1, 1, 1, 1] - 4*h[3, 2, 1, 1] + h[3, 3, 1]

Sage example in ./kschurnotes/notes-mike-anne.tex, line 2153::

    sage: mu = Partition([3,2,1])
    sage: c = mu.to_core(3)
    sage: w = c.to_grassmannian()
    sage: w.stanley_symmetric_function()
    4*m[1, 1, 1, 1, 1, 1] + 3*m[2, 1, 1, 1, 1] + 2*m[2, 2, 1, 1]
        + m[2, 2, 2] + 2*m[3, 1, 1, 1] + m[3, 2, 1]
    sage: w.reduced_words()
    [[2, 0, 3, 2, 1, 0], [0, 2, 3, 2, 1, 0], [0, 3, 2, 3, 1, 0],
     [0, 3, 2, 1, 3, 0]]

Sage example in ./kschurnotes/notes-mike-anne.tex, line 2202::

    sage: Sym = SymmetricFunctions(QQ)
    sage: ks = Sym.kschur(3,t=1)
    sage: h = Sym.homogeneous()
    sage: c = Partition([3,2,1]).to_core(3)
    sage: f = c.to_grassmannian().stanley_symmetric_function()*h([1]); f
    28*m[1, 1, 1, 1, 1, 1, 1] + 19*m[2, 1, 1, 1, 1, 1] + 12*m[2, 2, 1, 1, 1]
    + 7*m[2, 2, 2, 1] + 11*m[3, 1, 1, 1, 1] + 6*m[3, 2, 1, 1] + 3*m[3, 2, 2]
    + 2*m[3, 3, 1] + 2*m[4, 1, 1, 1] + m[4, 2, 1]
    sage: for la in Partitions(7, max_part=3):
    ...    print la, f.scalar( ks(la) )
    [3, 3, 1] 2
    [3, 2, 2] 1
    [3, 2, 1, 1] 3
    [3, 1, 1, 1, 1] 1
    [2, 2, 2, 1] 0
    [2, 2, 1, 1, 1] 0
    [2, 1, 1, 1, 1, 1] 0
    [1, 1, 1, 1, 1, 1, 1] 0

Sage example in ./kschurnotes/notes-mike-anne.tex, line 2225::

    sage: S = [la for la in Partitions(7, max_part=3) if f.scalar(ks(la))>0]
    sage: for p in S:
    ...     print p, SkewPartition([p.to_core(3).to_partition(),c.to_partition()])
    ...
    [3, 3, 1] [[7, 4, 1], [5, 2, 1]]
    [3, 2, 2] [[5, 2, 2], [5, 2, 1]]
    [3, 2, 1, 1] [[6, 3, 1, 1], [5, 2, 1]]
    [3, 1, 1, 1, 1] [[5, 2, 1, 1, 1], [5, 2, 1]]

Sage example in ./kschurnotes/notes-mike-anne.tex, line 2533::

  sage: W = WeylGroup(['A',3,1])
  sage: [w.reduced_word() for w in W.pieri_factors()]
  [[], [0], [1], [2], [3], [1, 0], [2, 0], [0, 3], [2, 1], [3, 1], [3, 2],
   [2, 1, 0], [1, 0, 3], [0, 3, 2], [3, 2, 1]]

Sage example in ./kschurnotes/notes-mike-anne.tex, line 2541::

  sage: A = NilCoxeterAlgebra(WeylGroup(['A',3,1]), prefix = 'A')
  sage: A.homogeneous_noncommutative_variables([2])
  A1*A0 + A2*A0 + A0*A3 + A3*A2 + A3*A1 + A2*A1

Sage example in ./kschurnotes/notes-mike-anne.tex, line 2548::

  sage: A.k_schur_noncommutative_variables([2,2])
  A0*A3*A1*A0 + A3*A1*A2*A0 + A1*A2*A0*A1 + A3*A2*A0*A3 + A2*A0*A3*A1
  + A2*A3*A1*A2

Sage example in ./kschurnotes/notes-mike-anne.tex, line 2555::

  sage: Sym = SymmetricFunctions(ZZ)
  sage: ks = Sym.kschur(5,t=1)
  sage: ks[2,1]*ks[2,1]
  ks5[2, 2, 1, 1] + ks5[2, 2, 2] + ks5[3, 1, 1, 1] + 2*ks5[3, 2, 1]
  + ks5[3, 3] + ks5[4, 2]

Sage example in ./kschurnotes/notes-mike-anne.tex, line 2683::

    sage: la = Partition([2,2])
    sage: la.k_conjugate(2).conjugate()
    [4]
    sage: la = Partition([2,1,1])
    sage: la.k_conjugate(2).conjugate()
    [3, 1]
    sage: la = Partition([1,1,1,1])
    sage: la.k_conjugate(2).conjugate()
    [2, 2]

Sage example in ./kschurnotes/notes-mike-anne.tex, line 2695::

    sage: Sym = SymmetricFunctions(FractionField(QQ['t']))
    sage: ks = Sym.kschur(2)
    sage: ks[2,2].omega_t_inverse()
    1/t^2*ks2[1, 1, 1, 1]
    sage: ks[2,1,1].omega_t_inverse()
    1/t*ks2[2, 1, 1]
    sage: ks[1,1,1,1].omega_t_inverse()
    1/t^2*ks2[2, 2]

Sage example in ./kschurnotes/notes-mike-anne.tex, line 2706::

    sage: Sym = SymmetricFunctions(FractionField(QQ['q,t']))
    sage: H = Sym.macdonald().H()
    sage: ks = Sym.kschur(2)
    sage: ks(H[2,2])
    q^2*ks2[1, 1, 1, 1] + (q*t+q)*ks2[2, 1, 1] + ks2[2, 2]
    sage: ks(H[2,1,1])
    q*ks2[1, 1, 1, 1] + (q*t^2+1)*ks2[2, 1, 1] + t*ks2[2, 2]
    sage: ks(H[1,1,1,1])
    ks2[1, 1, 1, 1] + (t^3+t^2)*ks2[2, 1, 1] + t^4*ks2[2, 2]

Sage example in ./kschurnotes/notes-mike-anne.tex, line 2858::

    sage: t = Tableau([[1,1,1,2,3,7],[2,2,3,5],[3,4],[4,5],[6]])
    sage: t.charge()
    9

Sage example in ./kschurnotes/notes-mike-anne.tex, line 3064::

    sage: la = Partition([3,2,1,1])
    sage: la.k_atom(4)
    [[[1, 1, 1], [2, 2], [3], [4]], [[1, 1, 1, 4], [2, 2], [3]]]

Sage example in ./kschurnotes/notes-mike-anne.tex, line 3163::

    sage: s = SymmetricFunctions(QQ['t']).schur()
    sage: G1 = s[1]
    sage: G211 = G1.hl_creation_operator([2,1]); G211
    s[2, 1, 1] + t*s[2, 2] + t*s[3, 1]
    sage: G3211 = G211.hl_creation_operator([3]); G3211
    s[3, 2, 1, 1] + t*s[3, 2, 2] + t*s[3, 3, 1] + t*s[4, 1, 1, 1]
     + (2*t^2+t)*s[4, 2, 1] + t^2*s[4, 3] + (t^3+t^2)*s[5, 1, 1]
     + 2*t^3*s[5, 2] + t^4*s[6, 1]

Sage example in ./kschurnotes/notes-mike-anne.tex, line 3529::

    sage: Sym = SymmetricFunctions(FractionField(QQ['t']))
    sage: ks4 = Sym.kschur(4)
    sage: ks4([3, 1, 1]).hl_creation_operator([1])
    (t-1)*ks4[2, 2, 1, 1] + t^2*ks4[3, 1, 1, 1] + t^3*ks4[3, 2, 1]
     + (t^3-t^2)*ks4[3, 3] + t^4*ks4[4, 1, 1]
    sage: ks4([3, 1, 1]).hl_creation_operator([2])
    t*ks4[3, 2, 1, 1] + t^2*ks4[3, 3, 1] + t^2*ks4[4, 1, 1, 1]
     + t^3*ks4[4, 2, 1]
    sage: ks4([3, 1, 1]).hl_creation_operator([3])
    ks4[3, 3, 1, 1] + t*ks4[4, 2, 1, 1] + t^2*ks4[4, 3, 1]
    sage: ks4([3, 1, 1]).hl_creation_operator([4])
    ks4[4, 3, 1, 1]

Sage example in ./kschurnotes/notes-mike-anne.tex, line 3596::

    sage: Sym = SymmetricFunctions(FractionField(QQ['t']))
    sage: ks3 = Sym.kschur(3)
    sage: ks3([3,2]).omega()
    Traceback (most recent call last):
    ...
    ValueError: t^2*s[1, 1, 1, 1, 1] + t*s[2, 1, 1, 1] + s[2, 2, 1] is not
    in the image of Generic morphism:
    From: 3-bounded Symmetric Functions over Fraction Field of Univariate
    Polynomial Ring in t over Rational Field in the 3-Schur basis
    To:   Symmetric Functions over Fraction Field of Univariate Polynomial Ring
    in t over Rational Field in the Schur basis

    sage: s = Sym.schur()
    sage: s(ks3[3,2])
    s[3, 2] + t*s[4, 1] + t^2*s[5]
    sage: t = s.base_ring().gen()
    sage: invert = lambda x: s.base_ring()(x.subs(t=1/t))
    sage: ks3(s(ks3([3,2])).omega().map_coefficients(invert))
    1/t^2*ks3[1, 1, 1, 1, 1]

Sage example in ./kschurnotes/notes-mike-anne.tex, line 3618::

    sage: ks3[3,2].omega_t_inverse()
    1/t^2*ks3[1, 1, 1, 1, 1]

Sage example in ./kschurnotes/notes-mike-anne.tex, line 3786::

    sage: Sym = SymmetricFunctions(FractionField(QQ['t']))
    sage: ks3 = Sym.kschur(3)
    sage: ks3[3,1].coproduct()
    ks3[] # ks3[3, 1] + ks3[1] # ks3[2, 1] + (t+1)*ks3[1] # ks3[3]
    + ks3[1, 1] # ks3[2] + ks3[2] # ks3[1, 1] + (t+1)*ks3[2] # ks3[2]
    + ks3[2, 1] # ks3[1] + (t+1)*ks3[3] # ks3[1]  + ks3[3, 1] # ks3[]

Sage example in ./kschurnotes/notes-mike-anne.tex, line 3820::

    sage: Sym = SymmetricFunctions(FractionField(QQ['t']))
    sage: ks2 = Sym.kschur(2)
    sage: ks3 = Sym.kschur(3)
    sage: ks5 = Sym.kschur(5)
    sage: ks5(ks3[2])*ks5(ks2[1])
    ks5[2, 1] + ks5[3]
    sage: ks5(ks3[2])*ks5(ks2[2,1])
    ks5[2, 2, 1] + ks5[3, 1, 1] + (t+1)*ks5[3, 2] + (t+1)*ks5[4, 1]
      + t*ks5[5]

Sage example in ./kschurnotes/notes-mike-anne.tex, line 3875::

    sage: Sym = SymmetricFunctions(FractionField(QQ['t']))
    sage: ks3 = Sym.kschur(3)
    sage: ks4 = Sym.kschur(4)
    sage: ks5 = Sym.kschur(5)
    sage: ks4(ks3[3,2,1,1])
    ks4[3, 2, 1, 1] + t*ks4[3, 3, 1] + t*ks4[4, 1, 1, 1] + t^2*ks4[4, 2, 1]
    sage: ks5(ks3[3,2,1,1])
    ks5[3, 2, 1, 1] + t*ks5[3, 3, 1] + t*ks5[4, 1, 1, 1] + t^2*ks5[4, 2, 1]
     + t^2*ks5[4, 3] + t^3*ks5[5, 1, 1]

    sage: ks5(ks4[3,2,1,1])
    ks5[3, 2, 1, 1]
    sage: ks5(ks4[4,3,3,2,1,1])
    ks5[4, 3, 3, 2, 1, 1] + t*ks5[4, 4, 3, 1, 1, 1]
     + t^2*ks5[5, 3, 3, 1, 1, 1]
    sage: ks5(ks4[4,3,3,2,1,1,1])
    ks5[4, 3, 3, 2, 1, 1, 1] + t*ks5[4, 3, 3, 3, 1, 1]
     + t*ks5[4, 4, 3, 1, 1, 1, 1] + t^2*ks5[4, 4, 3, 2, 1, 1]
     + t^2*ks5[5, 3, 3, 1, 1, 1, 1] + t^3*ks5[5, 3, 3, 2, 1, 1]
     + t^4*ks5[5, 4, 3, 1, 1, 1]

Sage example in ./kschurnotes/notes-mike-anne.tex, line 3942::

    sage: Sym = SymmetricFunctions(FractionField(QQ['q,t']))
    sage: H = Sym.macdonald().H()
    sage: ks = Sym.kschur(3)
    sage: ks(H[3])
    q^3*ks3[1, 1, 1] + (q^2+q)*ks3[2, 1] + ks3[3]
    sage: ks(H[3,2])
    q^4*ks3[1, 1, 1, 1, 1] + (q^3*t+q^3+q^2)*ks3[2, 1, 1, 1]
     + (q^3*t+q^2*t+q^2+q)*ks3[2, 2, 1]
     + (q^2*t+q*t+q)*ks3[3, 1, 1] + ks3[3, 2]
    sage: ks(H[3,1,1])
    q^3*ks3[1, 1, 1, 1, 1] + (q^3*t^2+q^2+q)*ks3[2, 1, 1, 1]
     + (q^2*t^2+q^2*t+q*t+q)*ks3[2, 2, 1]
     + (q^2*t^2+q*t^2+1)*ks3[3, 1, 1] + t*ks3[3, 2]

Sage example in ./kschurnotes/notes-mike-anne.tex, line 4047::

    sage: la = Partition([2,1])
    sage: w = la.to_core(3).to_grassmannian()
    sage: f = w.stanley_symmetric_function(); f
    2*m[1, 1, 1] + m[2, 1]

Sage example in ./kschurnotes/notes-mike-anne.tex, line 4055::

    sage: Sym = SymmetricFunctions(QQ)
    sage: p = Sym.power()
    sage: ks3 = Sym.kschur(3,1)
    sage: for mu in Partitions(5, max_part=3):
    ...    print mu, (f*p[2]).scalar( ks3(mu) )
    [3, 2] 1
    [3, 1, 1] 1
    [2, 2, 1] 0
    [2, 1, 1, 1] -1
    [1, 1, 1, 1, 1] -1

Sage example in ./kschurnotes/notes-mike-anne.tex, line 4153::

    sage: R = QQ[I]; z4 = R.zeta(4)
    sage: Sym = SymmetricFunctions(R)
    sage: ks3z = Sym.kschur(3,t=z4)
    sage: ks3 = Sym.kschur(3,t=1)
    sage: p = Sym.p()
    sage: p(ks3z[2, 2, 2, 2, 2, 2, 2, 2])
    1/12*p[4, 4, 4, 4] + 1/4*p[8, 8] + (-1/3)*p[12, 4]
    sage: p(ks3[2,2])
    1/12*p[1, 1, 1, 1] + 1/4*p[2, 2] + (-1/3)*p[3, 1]
    sage: p(ks3[2,2]).plethysm(p[4])
    1/12*p[4, 4, 4, 4] + 1/4*p[8, 8] + (-1/3)*p[12, 4]

Sage example in ./kschurnotes/notes-mike-anne.tex, line 4167::

    sage: ks3z[3, 3, 3, 3]*ks3z[2, 1]
    ks3[3, 3, 3, 3, 2, 1]

Sage example in ./kschurnotes/notes-mike-anne.tex, line 4328::

    sage: Sym = SymmetricFunctions(QQ['t'].fraction_field())
    sage: HLQp = Sym.hall_littlewood().Qp()
    sage: HLP = Sym.hall_littlewood().P()
    sage: def dual_k_schur(k, la):
    ...     ks = Sym.kschur(k)
    ...     return sum( ks(HLQp(mu)).coefficient(la)*HLP(mu)
    ...         for mu in Partitions(add(la), max_part=k) )
    sage: ks3 = Sym.kschur(3)
    sage: f = dual_k_schur(3,[2,1,1])*dual_k_schur(3,[3,2,1])
    sage: for mu in Partitions(10,max_part=3):
    ...     print mu, f.scalar(ks3(mu))
    [3, 3, 3, 1] t^5 + 3*t^4 + 3*t^3 + 4*t^2 + 2*t + 1
    [3, 3, 2, 2] t^4 + t^3 + 3*t^2 + 2*t + 1
    [3, 3, 2, 1, 1] 2*t^5 + 3*t^4 + 5*t^3 + 6*t^2 + 4*t + 2
    [3, 3, 1, 1, 1, 1] t^5 + t^4 + 4*t^3 + 4*t^2 + 3*t + 1
    [3, 2, 2, 2, 1] t^4 + 3*t^3 + 4*t^2 + 3*t + 1
    [3, 2, 2, 1, 1, 1] 2*t^2 + t + 1
    [3, 2, 1, 1, 1, 1, 1] 2*t^5 + 3*t^4 + 4*t^3 + 3*t^2 + t
    [3, 1, 1, 1, 1, 1, 1, 1] t^5 + 2*t^4 + t^3
    [2, 2, 2, 2, 2] t^5 + 2*t^4 + 2*t^3 + t^2
    [2, 2, 2, 2, 1, 1] t^3 + t^2
    [2, 2, 2, 1, 1, 1, 1] t^4 + t^3 + t^2
    [2, 2, 1, 1, 1, 1, 1, 1] 0
    [2, 1, 1, 1, 1, 1, 1, 1, 1] t^7 + t^6
    [1, 1, 1, 1, 1, 1, 1, 1, 1, 1] 0

Sage example in ./kschurnotes/notes-mike-anne.tex, line 4433::

    sage: ks2 = SymmetricFunctions(QQ['t']).kschur(2)
    sage: HLQp = SymmetricFunctions(QQ['t']).hall_littlewood().Qp()
    sage: ks2( (HLQp(ks2[1,1])*HLQp(ks2[1])).restrict_parts(2) )
    ks2[1, 1, 1] + (-t+1)*ks2[2, 1]

"""
