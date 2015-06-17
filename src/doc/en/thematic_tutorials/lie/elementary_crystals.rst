-------------------
Elementary crystals
-------------------

T-crystal
---------

Let `\lambda` be a weight. As defined in [Kashiwara1993]_ the crystal
`T_{\lambda} = \{ t_{\lambda} \}` is a single element crystal with the
crystal structure defined by

.. MATH::

    \mathrm{wt}(t_\lambda) = \lambda, \quad
    e_i t_{\lambda} = f_i t_{\lambda} = 0, \quad
    \varepsilon_i(t_{\lambda}) = \varphi_i(t_{\lambda}) = -\infty.

The crystal `T_{\lambda}` shifts the weights of the vertices in a crystal
`B` by `\lambda` when tensored with `B`, but leaves the graph structure of
`B` unchanged. That is, for all `b \in B`, we have `\mathrm{wt}(t_\lambda
\otimes b) = \mathrm{wt}(b) + \lambda`::

    sage: B = crystals.Tableaux(['A',2],shape=[2,1])
    sage: T = crystals.elementary.T(['A',2], B.Lambda()[1] + B.Lambda()[2])
    sage: V = crystals.TensorProduct(T,B)
    sage: for x in V:
    ....:     print x.weight()
    ....:
    (4, 2, 0)
    (3, 3, 0)
    (3, 2, 1)
    (3, 1, 2)
    (2, 2, 2)
    (4, 1, 1)
    (3, 2, 1)
    (2, 3, 1)
    sage: for x in B:
        print x.weight() + T[0].weight()
    ....:
    (4, 2, 0)
    (3, 3, 0)
    (3, 2, 1)
    (3, 1, 2)
    (2, 2, 2)
    (4, 1, 1)
    (3, 2, 1)
    (2, 3, 1)


C-crystal
---------

Defined in [Kashiwara1993]_, the component crystal `C = \{c\}` is the single
element crystal whose crystal structure is defined by

.. MATH::

    \mathrm{wt}(c) = 0, \quad
    e_i c = f_i c = 0, \quad
    \varepsilon_i(c) = \varphi_i(c) = 0.

Note `C \cong B(0)`, where `B(0)` is the highest weight crystal of highest
weight `0`.


R-crystal
---------

For a fixed weight `\lambda`, the crystal `R_{\lambda} = \{ r_{\lambda} \}`
is a single element crystal with the crystal structure defined by

.. MATH::

    \mathrm{wt}(r_{\lambda}) = \lambda, \quad
    e_i r_{\lambda} = f_i r_{\lambda} = 0, \quad
    \varepsilon_i(r_{\lambda}) = -\langle h_i, \lambda\rangle, \quad
    \varphi_i(r_{\lambda}) = 0,

where `\{h_i\}` are the simple coroots.

Tensoring `R_{\lambda}` with a crystal `B` results in shifting the weights
of the vertices in `B` by `\lambda` and may also cut a subset out of the
original graph of `B`.  That is, `\mathrm{wt}(r_{\lambda} \otimes b) =
\mathrm{wt}(b) + \lambda`, where `b \in B`, provided `r_{\lambda} \otimes
b \neq 0`. For example, the crystal graph of `B(\lambda)` is the same as
the crystal graph of `R_{\lambda} \otimes B(\infty)` generated from the
component `r_{\lambda} \otimes u_{\infty}`.


`i`-th elementary crystal
-------------------------

For `i` an element of the index set of type `X`, the crystal `B_i` of type
`X` is the set

.. MATH::

    B_i = \{ b_i(m) : m \in \ZZ \},

where the crystal stucture is given by

.. MATH::

    \begin{aligned}
    \mathrm{wt}\bigl(b_i(m)\bigr) &= m\alpha_i \\
    \varphi_j\bigl(b_i(m)\bigr) &= \begin{cases}
        m & \text{ if } j=i, \\
        -\infty & \text{ if } j\neq i,
    \end{cases} \\
    \varepsilon_j\bigl(b_i(m)\bigr) &= \begin{cases}
        -m & \text{ if } j=i, \\
        -\infty & \text{ if } j\neq i,
    \end{cases} \\
    e_j b_i(m) &= \begin{cases}
        b_i(m+1) & \text{ if } j=i, \\
        0 & \text{ if } j\neq i,
    \end{cases} \\
    f_j b_i(m) &= \begin{cases}
        b_i(m-1) & \text{ if } j=i, \\
        0 & \text{ if } j\neq i.
    \end{cases}
    \end{aligned}
