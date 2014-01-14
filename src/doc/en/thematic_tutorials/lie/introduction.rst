--------------------------
The Scope of this Document
--------------------------


Lie groups and algebras
-----------------------

Sage can be used to do standard computations for Lie groups and Lie
algebras. The following categories of representations are equivalent:

- Complex representations of a compact, semisimple simply connected
  Lie group `G`.

- Complex representations of its Lie algebra `\mathfrak{g}`. This is a
  real Lie algebra, so representations are not required to be complex
  linear maps.

- Complex representations of its complexified Lie algebra
  `\mathfrak{g}_{\mathbf{C}} = \mathbf{C} \otimes \mathfrak{g}`. This
  is a complex Lie algebra and representations are required to be
  complex linear transformations.

- The complex analytic representations of the semisimple
  simply-connected complex analytic group `G_{\mathbf{C}}` having
  `\mathfrak{g}_{\mathbf{C}}` as its Lie algebra.

- Modules of the universal enveloping algebra
  `U(\mathfrak{g}_{\mathbf{C}})`.

- Modules of the quantized enveloping algebra
  `U_q(\mathfrak{g}_{\mathbf{C}})`.

For example, we could take `G = SU(n)`,
`\mathfrak{g} = \mathfrak{sl}(n, \mathbf{R})`,
`\mathfrak{g}_{\mathbf{C}} = \mathfrak{sl}(n, \mathbf{C})` and
`G_{\mathbf{C}} = SL(n, \mathbf{C})`. Because these categories are the same, their
representations may be studied simultaneously. The above equivalences
may be expanded to include reductive groups like `U(n)` and `GL(n)`
with a bit of care.

Here are some typical problems that can be solved using Sage:

- Decompose a module in any one of these categories into irreducibles.

- Compute the Frobenius-Schur indicator of an irreducible module.

- Compute the tensor product of two modules.

- If `H` is a subgroup of `G`, study the restriction of modules for
  `G` to `H`. The solution to this problem is called a *branching rule*.

- Find the multiplicities of the weights of the representation.

In addition to its representations, which we may study as above, a Lie
group has various related structures. These include:

- The Weyl Group `W`.

- The Weight Lattice.

- The Root System

- The Cartan Type.

- The Dynkin diagram.

- The extended Dynkin diagram.

Sage contains methods for working with these structures.

If there is something you need that is not implemented, getting it
added to Sage will likely be possible. You may write your own
algorithm for an unimplemented task, and if it is something others
will be interested in, it is probably possible to get it added to
Sage.


Combinatorics
-------------

Sage supports a great many related mathematical objects. Some of these
properly belong to combinatorics. It is beyond the scope of these
notes to cover all the combinatorics in Sage, but we will try to touch
on those combinatorial methods which have some connection with Lie
groups and representation theory. These include:

- The affine Weyl group, an infinite group containing `W`.

- Kashiwara crystals, which are combinatorial analogs of modules in
  the above categories.

- Coxeter group methods applicable to Weyl groups and the affine Weyl
  group, such as Bruhat order.

- The Iwahori Hecke algebras, which are deformations of the group
  algebras of `W` and the affine Weyl group.

- Kazhdan-Lusztig polynomials.
