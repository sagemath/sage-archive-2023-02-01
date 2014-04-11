.. linkall

===================
The Tsetlin library
===================

Introduction
------------

In this section, we study a simple random walk (or Markov chain),
called the *Tsetlin library*. It will give us the opportunity to see
the interplay between combinatorics, linear algebra, representation
theory and computer exploration, without requiring heavy theoretical
background. I hope this encourages everyone to play around with
this or similar systems and investigate their properties!

It has been known for several years that the theory of group
representations can facilitate the study of systems whose evolution
is random (Markov chains), breaking them down into simpler systems.
More recently it has been realized that generalizing this (namely
replacing the invertibility axiom for groups by other axioms)
could explain the behavior of other particularly simple Markov chains
such as the Tsetlin library.

The Tsetlin library
-------------------

Consider a bookshelf in a library containing `n` distinct books. When
a person borrows a book and then returns it, it gets placed back on the
shelf to the right of all books. This is what we naturally do with our
pile of shirts in the closet: after use and cleaning, the shirt is placed
on the top of its pile. Hence the most popular books/shirts will more
likely appear on the right/top of the shelf/pile.

This type of organization has the advantage of being self-adaptive:

- The books most often used accumulate on the right and thus can easily
  be found.
- If the use changes over time, the system adapts.

In fact, this type of strategy is used not only in everyday life, but also
in computer science. The natural questions that arise are:

- *Steady state*: To which state(s) does the system converge to? This,
  among other things, is used to evaluate the average access time to
  a book.
- *The rate of convergence*: How fast does the system adapt to a changing
  environment .

Let us formalize the description. The Tsetlin library is a discrete Markov
chain (discrete time, discrete state space) described by:

- The state space `\Omega_n` is given by the set of all permutations of the
  `n` books.
- The transition operators are denoted by `\tau_i \colon \Omega_n \to \Omega_n`.
  When `\tau_i` is applied to a permutation `\sigma`, the number `i` is moved
  to the end of the permutation.
- We assign parameters `x_i \ge 0` for all `1\le i\le n` with
  `\sum_{i=1}^n x_i = 1`. The parameter `x_i` indicates the probability of
  choosing the operator `\tau_i`.
