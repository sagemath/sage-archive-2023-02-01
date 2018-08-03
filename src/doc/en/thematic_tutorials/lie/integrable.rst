Integrable Highest Weight Representations of Affine Lie algebras
================================================================

In this section `\mathfrak{g}` can be an arbitrary Kac-Moody Lie Algebra
made with a symmetrizable, indecomposable Cartan matrix.

Suppose that `V` is a representation with a weight decomposition as in
:ref:`roots_and_weights`.  Let `\alpha` be a real root, and let
`\mathfrak{g}_\alpha` be the corresponding root space, that is,
the one-dimensional weight space for `\alpha` in the adjoint
representation of `\mathfrak{g}` on itself. Then `-\alpha` is also a root. The
two one-dimensional spaces `\mathfrak{g}_\alpha` and `\mathfrak{g}_{-\alpha}`
generate a Lie algebra isomorphic to `\mathfrak{sl}_2`. The module `V` is
called *integrable* if for each such `\alpha` the representation of
`\mathfrak{sl}_2` obtained this way integrates to a representation of the Lie
group `\text{SL}_2`. Since this group contains an element that stabilizes
`\mathfrak{h}` and induces the corresponding simple reflection on the weight
lattice, integrability implies that the weight multiplicities are invariant
under the action of the Weyl group.

If the Kac-Moody Lie algebra `\mathfrak{g}` is finite-dimensional
then the integrable highest weight representations are
just the irreducible finite-dimensional ones. For a general
Kac-Moody Lie algebra the integrable highest weight representations
are the analogs of the finite-dimensional ones,
discussed in :file:`weyl_character_ring`, and their
theory has many aspects in common with the finite-dimensional
representations of finite-dimensional simple Lie algebras,
such as the parametrization by dominant weights, and
generalizations of the Weyl denominator and character
formulas, due to Macdonald and Kac respectively.

If `\lambda` is a dominant weight, then the irreducible
highest weight module `L(\lambda)` defined in :ref:`roots_and_weights`
is integrable. Moreover every highest weight integrable representation arises
this way, so these representations are in bijection with the cone `P^+` of
dominant weights.

The affine case
---------------

Now we assume that `\mathfrak{g}` is affine. The integrable
highest weight representations (and their crystals) are
extremely interesting. Integrable highest weight representations of
`\mathfrak{g}` arise in a variety of contexts, from string
theory to the modular representation theory of the symmetric
groups, and the theory of modular forms. One particular
representation, `L(\Lambda_0)` is particularly ubiquitous,
and this is called the *basic representation*.

Therefore in [KMPS]_ (published in 1990) tabulated data for
many of these representations. They wrote

    We present a vast quantity of numerical data in tabular form, this
    being the only source for such information. The computations are tedious
    and not particularly straightforward when it is necessary to carry them
    out individually. We hope the appearance of this book will spur interest
    in a field that has become, in barely 20 years, deeply rewarding and
    full of promise for the future. It would indeed be gratifying if these
    tables were to appear to the scientists of 2040 as obsolete as the
    dust-gathering compilations of transcendental functions appear for us
    today because of their availability on every pocket calculator.

As we will explain, Sage can reproduce the contents of these tables. 
Moreover the tables in [KMPS]_ are limited to the untwisted types,
but Sage also implements the twisted types.

Although Sage can reproduce the tables in the second volume of [KMPS]_, the
work remains very useful. The first volume is a down-to-earth
and very helpful exposition of the theory of integrable representations of
affine Lie algebras with explicit examples and explanations of the
connections with mathematical physics and vertex operators.

The support of an integrable highest weight representation
----------------------------------------------------------

Let `\Lambda\in P^+` and let `V=L(\lambda)` be the integrable representation
with highest weight `\Lambda`. If `\mu` is another weight, let `\text{mult}(\mu)` denote the
multiplicity of the weight `\mu` in `L(\lambda)`. Define the
*support* of the representation `\text{supp}(V)` to be the set
of `\mu` such that `\text{mult}(\mu) > 0`.

If `\text{mult}(\mu)>0` then `\lambda-\mu` is a linear combination
of the simple roots with nonnegative integer coefficients.
Moreover `\text{supp}(V)` is contained in the paraboloid

.. MATH::

    (\Lambda+\rho | \Lambda+\rho) - (\mu+\rho | \mu+\rho) \geq 0

where `(\, | \,)` is the invariant inner product on the weight
lattice and `\rho` is the Weyl vector (:ref:`untwisted_affine`).
Moreover if `\mu\in\text{supp}(V)` then `\Lambda-\mu`
is an element of the root lattice `Q` ([Kac]_, Propositions 11.3 and 11.4).
    
We organize the weight multiplicities into sequences called *string functions*
or *strings* as follows. By [Kac]_, Proposition 11.3 or Corollary 11.9, for fixed `\mu`
the function `\text{mult}(\mu - k\delta)` of `k` is a increasing sequence.
We adjust `\mu` by a multiple of `\delta` to the beginning
of the positive part of the sequence. Thus we define
`\mu` to be *maximal* if `\text{mult}(\mu) \neq 0` but `\text{mult}(\mu + \delta) = 0`.

Since `\delta` is fixed under the action of the affine Weyl group, and since
the weight multiplicities are Weyl group invariant, the function
`k \mapsto \text{mult}(\mu - k \delta)` is unchanged if `\mu` is replaced by `w(\mu)`
for some Weyl group element `w`. Now every Weyl orbit contains a dominant
weight.  Therefore in enumerating the string we may assume that the weight
`\mu` is dominant. There are only a finite number of dominant maximal
weights. Thus there are only a finite number of such strings to be computed.

Remarkably, [KacPeterson]_ showed that each string is the set of Fourier
coefficients of a modular form; see also [Kac]_ Chapters 12 and 13. To this end
we define the *modular characteristic*

.. MATH::

    m_\Lambda = \frac{|\Lambda+\rho|^2}{2(k+h^\vee)} - \frac{|\rho|^2}{2h^\vee}

Here `k=(\Lambda|\delta)` is the *level* of the representation and
`h^\vee` is the dual Coxeter number (:ref:`coxeternumber`).
If `\mu` is a weight, define

.. MATH::

    m_{\Lambda,\mu} = m_\Lambda - \frac{|\mu|^2}{2k}.

Let `\lambda` be a weight, which we may assume maximal. Then Kac and Peterson
defined the *string function*

.. MATH::

    c_\mu^\Lambda = q^{m_{\Lambda,\mu}}\sum_{n\in\ZZ}\text{mult}(\mu-n\delta)q^n.

It is a modular form. See [Kac]_, Corollary 13.10. Although these
do arise as partition functions in string theory, the term "string" here does
not refer to physical strings.

Sage methods for integrable representations
-------------------------------------------

In this section we will show how to use Sage to compute with
integrable highest weight Lie algebras.
For further documentation, see the reference manual
:class:`~sage.combinat.root_system.integrable_representations.IntegrableRepresentation`

In the following example, we work with the integrable representation
with highest weight `2\Lambda_0` for `\widehat{\mathfrak{sl}}_2`,
that is, `A_1^{(1)}`. First we create a dominant weight in
the extended weight lattice, then create the ``IntegrableRepresentation``
class. We compute the strings. There are
two, since there are two dominant multiple weights. One of them
is the highest weight `2\Lambda_0`, and the other is `2\Lambda_1-\delta`.
We apply the simple reflection `s_0` to the second, giving
`2\Lambda_1-\delta`, a maximal weight that is not dominant.
Then we compute the string function at this weight, which we see
agrees with the string function for the corresponding dominant
maximal weight::

    sage: L = RootSystem("A1~").weight_lattice(extended=True)
    sage: Lambda = L.fundamental_weights()
    sage: delta = L.null_root()
    sage: W = L.weyl_group(prefix="s")
    sage: [s0,s1]=W.simple_reflections()
    sage: V = IntegrableRepresentation(2*Lambda[0])
    sage: V.strings()
    {2*Lambda[0]: [1, 1, 3, 5, 10, 16, 28, 43, 70, 105, 161, 236],
     2*Lambda[1] - delta: [1, 2, 4, 7, 13, 21, 35, 55, 86, 130, 196, 287]}
     sage: [mw1,mw2] = V.dominant_maximal_weights(); mw1,mw2
     (2*Lambda[0], 2*Lambda[1] - delta)
     sage: s0.action(mw2)
     2*Lambda[1] - delta
     sage: [V.mult(s0.action(mw2)-k*delta) for k in [0..10]]
     [1, 2, 4, 7, 13, 21, 35, 55, 86, 130, 196]







