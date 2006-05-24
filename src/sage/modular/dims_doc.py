"""
Dimensions of spaces of modular forms

The dimension formulas and implementations in this module grew out of
a program that Bruce Caskel wrote (around 1996) in PARI, which Kevin
Buzzard subsequently extended.  I (William Stein) then implemented it
in C++ for HECKE.  I also implemented it in MAGMA.  Also, the
functions for dimensions of spaces with nontrivial character are based
on a paper (that has no proofs) by Cohen and Oesterle (Springer
Lecture notes in math, volume 627, pages 69--78).

The formulas here are more complete than in HECKE or MAGMA.

Currently the input to each function below is an integer and either a
Dirichlet character $\eps$ or a congruence subgroup, which must be
either $\Gamma_0(N)$ or $\Gamma_1(N)$.  If the input is a Dirichlet
character $\eps$, the dimensions are for subspaces of
$M_k(\Gamma_1(N), \eps)$, where $N$ is the modulus of $\eps$.
"""

__doc_exclude = []


from dims import (dimension_cusp_forms,
                  dimension_eis,
                  dimension_modular_forms,
                  dimension_new_cusp_forms)

