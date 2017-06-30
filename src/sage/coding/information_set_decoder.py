# -*- coding: utf-8 -*-
r"""
Information set decoding for linear codes

Information-set decoding is a probabilistic decoding strategy that
essentially tries to guess `k` correct positions in the received word,
where `k` is the dimension of the code. A codeword agreeing with the
received word on the guessed position can easily be computed, and their
difference is one possible error vector. A "correct" guess is assumed when
this error vector has low Hamming weight.

This simple algorithm is not very efficient in itself, but there are numerous
refinements to the strategy that make it very capable over rather large codes.
Still, the decoding algorithm is exponential in dimension of the code and the
log of the field size.

The ISD strategy requires choosing how many errors is deemed acceptable. One
choice could be `d/2`, where `d` is the minimum distance of the code, but
sometimes `d` is not known, or sometimes more errors are expected. If one
chooses anything above `d/2`, the algorithm does not guarantee to return a
nearest codeword.

AUTHORS:

- David Lucas,Yann Laigle-Chapuy, Johan Rosenkilde (2016-02, 2017-06): initial
  version
"""

#******************************************************************************
#       Copyright (C) 2005 David Joyner <wdjoyner@gmail.com>
#                     2006 William Stein <wstein@gmail.com>
#
#  Distributed under the terms of the GNU General Public License (GPL),
#  version 2 or later (at your preference).
#
#                  http://www.gnu.org/licenses/
#******************************************************************************
# python3
from __future__ import division, print_function, absolute_import
from six.moves import range
from six import iteritems

from sage.rings.integer_ring import ZZ
from sage.rings.integer import Integer
from sage.modules.free_module_element import vector
from .decoder import Decoder, DecodingError
from sage.arith.all import binomial

class AbstractInformationSetDecoder(Decoder):
    r"""
    Abstract base class for Information-set decoder for linear codes.

    For a description of the information-set decoding paradigm (ISD), see
    :func:`sage.coding.LinearCodeInformationSetDecoder`. There are many
    variants of information-set decoding, and this class serves as a base
    class for all such implementations.

    To sub-class this class, override ``decode_to_code`` and call the super
    constructor from ``__init__``.

    INPUT:

    - ``code`` -- A linear code for which to decode.

    - ``number_errors`` -- an integer, the maximal number of errors to accept as
      correct decoding. An interval can also be specified by giving a pair of
      integers, where both end values are taken to be in the interval.

    - ``algorithm_name`` -- A name for the specific ISD algorithm used.

    EXAMPLES::

    Sage ships with some information set decoders available::

        sage: C = codes.GolayCode(GF(2))
        sage: D = C.decoder("InformationSet", 2)
        sage: D
        Information set decoder (LeeBrickell) for [24, 12, 8] Extended Golay code over GF(2) decoding up to 2 errors

    It is straightforward to define your own::

        sage: from sage.coding.isd_decoding import AbstractInformationSetDecoder
        sage: class MinimalISD(AbstractInformationSetDecoder):
        ....:   def __init__(self, code, number_errors):
        ....:       super(MinimalISD, self).__init__(code, number_errors, "MyISD")
        ....:   def decode_to_code(self, r):
        ....:       # Here goes your ISD algorithm
        ....:       raise NotImplementedError
        sage: D = MinimalISD(C, 3); D
        Information set decoder (MyISD) for [24, 12, 8] Extended Golay code over GF(2) decoding up to 3 errors
    """

    def __init__(self, code, number_errors, algorithm_name):
        r"""
        TESTS:

        ``number_errors`` has to be either a list of Integers/ints, a tuple of Integers/ints,
        or an Integer/int::

            sage: C = codes.GolayCode(GF(2))
            sage: D = C.decoder("InformationSet", "aa")
            Traceback (most recent call last):
            ...
            ValueError: number_errors must be an integer or a pair of integers

        If ``number_errors`` is passed as a list/tuple, it has to contain only two values,
        the first one being at most the second one::

            sage: C = codes.GolayCode(GF(2))
            sage: D = C.decoder("InformationSet", (4, 2))
            Traceback (most recent call last):
            ...
            ValueError: number_errors should be a positive integer or a valid interval within the positive integers
        """
        if isinstance(number_errors, (Integer, int)):
            number_errors = (0, number_errors)
        if isinstance(number_errors, (tuple, list)):
            if not len(number_errors) == 2:
                raise ValueError("number_errors should be either an integer"
                                 " or a pair of integers")
            if not (number_errors[0] in ZZ and number_errors[1] in ZZ):
                raise ValueError("All elements of number_errors have to be"
                                 " positive integers")
            if 0 > number_errors[0] or number_errors[0] > number_errors[1]:
                raise ValueError(
                        "number_errors should be a positive integer or"
                        " a valid interval within the positive integers")
            if number_errors[1] > code.length():
                raise ValueError("The provided number of errors should be at"
                                 " most the code's length")
        else:
            raise ValueError("number_errors must be an integer or a pair of integers")

        self._number_errors = number_errors
        self._algorithm_name = algorithm_name
        super(AbstractInformationSetDecoder, self).__init__(
            code, code.ambient_space(), code._default_encoder_name)

    def _format_number_errors(self):
        r"""
        Format the number of errors when calling ``_repr_`` or ``_latex_``.

        EXAMPLES::

            sage: C = codes.GolayCode(GF(2))
            sage: D = C.decoder("InformationSet", 3)
            sage: D._format_number_errors()
            'up to 3'
            sage: D = C.decoder("InformationSet", (2,3))
            sage: D._format_number_errors()
            'between 2 and 3'
            sage: D = C.decoder("InformationSet", (3,3))
            sage: D._format_number_errors()
            'exactly 3'
        """
        if self._number_errors[0] == 0:
            return "up to {0}".format(self._number_errors[1])
        if self._number_errors[0] == self._number_errors[1]:
            return "exactly {0}".format(self._number_errors[0])
        return "between {0} and {1}".format(self._number_errors[0], self._number_errors[1])

    def _repr_(self):
        r"""
        Returns a string representation of this decoder.

        EXAMPLES::

            sage: C = codes.GolayCode(GF(2))
            sage: D = C.decoder("InformationSet", 2)
            sage: D
            Information set decoder (LeeBrickell) for [24, 12, 8] Extended Golay code over GF(2) decoding up to 2 errors
        """
        return "Information set decoder ({}) for {} decoding {} errors ".format(self._algorithm_name, self.code(), self._format_number_errors())

    def _latex_(self):
        r"""
        Returns a latex representation of this decoder.

        EXAMPLES::

            sage: C = codes.GolayCode(GF(2))
            sage: D = C.decoder("InformationSet", 2)
            sage: latex(D)
            \textnormal{Information set decoder (LeeBrickell) for }[24, 12, 8] \textnormal{ Extended Golay Code over } \Bold{F}_{2} \textnormal{decoding up to 2 errors}
        """
        return "\\textnormal{{Information set decoder ({}) for }}{} \\textnormal{{decoding {} errors}}".format(self._algorithm_name, self.code()._latex_(), self._format_number_errors())

    def decoding_radius(self):
        r"""
        Return the maximal number of errors this decoder can decode.

        EXAMPLES::

            sage: C = codes.GolayCode(GF(2))
            sage: D = C.decoder("InformationSet", 2)
            sage: D.decoding_radius()
            2
        """
        return self._number_errors[1]

    def decoding_interval(self):
        r"""
         A pair of integers specifying the interval of number of errors this
         decoder will attempt to correct.

         The interval includes both end values.


        EXAMPLES::

            sage: C = codes.GolayCode(GF(2))
            sage: D = C.decoder("InformationSet", 2)
            sage: D.decoding_interval()
            (0, 2)
        """
        return self._number_errors




class LinearCodeISD_LeeBrickell(AbstractInformationSetDecoder):
    r"""
    Information-set decoder for any linear code using the Lee--Brickell
    algorithm.

    For a description of the information-set decoding paradigm (ISD), see
    :func:`sage.coding.LinearCodeInformationSetDecoder`.

    This implements the Lee--Brickell variant of ISD, see [LB1988] for the
    original binary case, and [Pet10] for the `q`-ary extension.

    Let `C` be a `[n, k]`-linear code over `GF(q)`, and let `r \in GF(q)^{n}` be
    a received word in a transmission. We seek the codeword whose Hamming
    distance from `r` is minimal. Let `p` and `w` be integers, such that `0\leq
    p\leq w`, Let `G` be a generator matrix of `C`, and for any set of indices
    `I`, we write `G_{I}` for the matrix formed by the columns of `G` indexed by
    `I`. The Lee--Brickell ISD loops the following until it is successful:

        1. Choose an information set `I` of `C`.
        2. Compute `r' = r - r_{I}\times G_I^{-1} \times G`
        3. Consider every size-`p` subset of `I`, `\{a_1, \dots, a_p\}`.
           For each `m = (m_1, \dots, m_p) \in GF(q)^{p}`, compute
           the error vector `e = r' - \sum_{i=1}^{p} m_i\times g_{a_i}`,
        4. If `e` has a Hamming weight at most `w`, return `y-e`.

    INPUT:

    - ``code`` -- A linear code for which to decode.

    - ``number_errors`` -- an integer, the maximal number of errors to accept as
      correct decoding. An interval can also be specified by giving a pair of
      integers, where both end values are taken to be in the interval.

    - ``search_size`` -- (optional) the size of subsets to use on step 3 of the
      algorithm as described above. Usually a small number. It has to be at most
      the largest allowed number of errors. A good choice will be approximated
      if this option is not set; see
      :meth:`sage.coding.LinearCodeISD_LeeBrickell._calibrate_search_size`
      for details.

    EXAMPLES::

        sage: C = codes.GolayCode(GF(2))
        sage: D = C.decoder("InformationSet", 2, algorithm="LeeBrickell")
        sage: D
        Information set decoder (LeeBrickell) for [24, 12, 8] Extended Golay code over GF(2) decoding up to 2 errors

        sage: C = codes.GolayCode(GF(2))
        sage: D = C.decoder("InformationSet", (2,3), algorithm="LeeBrickell")
        sage: D
        Information set decoder (LeeBrickell) for [24, 12, 8] Extended Golay code over GF(2) decoding between 2 and 3 errors

    We can also information set decode non-binary codes::

        sage: C = codes.GolayCode(GF(3))
        sage: D = C.decoder("InformationSet", 2, algorithm="LeeBrickell")
        sage: D
        Information set decoder (LeeBrickell) for [12, 6, 6] Extended Golay code over GF(3) decoding up to 2 errors

    The decoder class can be invoked in two other ways as well::

        sage: C = codes.GolayCode(GF(2))
        sage: D = codes.decoders.LinearCodeInformationSetDecoder(C, 3, algorithm="LeeBrickell")
        sage: D
        Information set decoder (LeeBrickell) for [24, 12, 8] Extended Golay code over GF(2) decoding up to 3 errors

        sage: from sage.coding.isd_decoding import LinearCodeISD_LeeBrickell
        sage: C = codes.GolayCode(GF(2))
        sage: D = LinearCodeISD_LeeBrickell(C, 3)
        sage: D
        Information set decoder (LeeBrickell) for [24, 12, 8] Extended Golay code over GF(2) decoding up to 3 errors

    """
    def __init__(self, code, number_errors, search_size = None):
        r"""
        TESTS:

        If ``search_size`` is bigger than a possible value for ``number_errors``, an error
        will be raised::

            sage: C = codes.GolayCode(GF(2))
            sage: D = C.decoder("InformationSet", (1, 3), algorithm="LeeBrickell", search_size=5)
            Traceback (most recent call last):
            ...
            ValueError: The search size parameter has to be at most the maximal number of allowed errors
        """
        super(LinearCodeISD_LeeBrickell, self).__init__(code, number_errors, "LeeBrickell")
        if not search_size is None:
            if not isinstance(search_size, (Integer, int)) or search_size < 0:
                raise ValueError("The search size parameter has to be a positive integer")
            if search_size > self.decoding_interval()[1]:
                raise ValueError("The search size parameter has to be at most"
                                 " the maximal number of allowed errors")
            self._search_size = search_size
            self._search_size_specified = True
        else:
            self._search_size_specified = False

    def __eq__(self, other):
        r"""
        Tests equality between information set decoder objects.

        EXAMPLES::

            sage: C = codes.GolayCode(GF(2))
            sage: from sage.coding.isd_decoding import LinearCodeISD_LeeBrickell
            sage: D = LinearCodeISD_LeeBrickell(C, 4)
            sage: D == LinearCodeISD_LeeBrickell(C, 4)
            True
            sage: D == LinearCodeISD_LeeBrickell(C, 5)
            False
            sage: other_search = 1 if D.search_size() != 1 else 2
            sage: D == LinearCodeISD_LeeBrickell(C, 4, search_size=other_search)
            False

        Lee-Brickell ISD objects can be equal only if they have both calibrated
        the search size, or if they both had it set and to the same value::

            sage: D2 = LinearCodeISD_LeeBrickell(C, 4, search_size=D.search_size())
            sage: D == D2
            False
            sage: D2 == LinearCodeISD_LeeBrickell(C, 4, search_size=D.search_size())
            True
        """
        return isinstance(other, LinearCodeISD_LeeBrickell)\
                and self.code() == other.code()\
                and self.decoding_interval() == other.decoding_interval()\
                and self._search_size_specified == other._search_size_specified\
                and (not self._search_size_specified or self.search_size() == other.search_size())

    def _lee_brickell_algorithm(self, r, tau, p):
        r"""
        The Lee-Brickell algorithm as described in the class doc.

        INPUT:

        - `r` -- a received word, i.e. a vector in the ambient space of
          :meth:`decoder.Decoder.code`.

        - `tau` -- an integer interval of acceptable Hamming weights for the error.

        - `p` -- an integer, the search size.

        OUTPUT: A codeword whose distance to `r` satisfies `tau`.

        EXAMPLES::

            sage: M = matrix(GF(2), [[1, 0, 0, 0, 0, 1, 0, 0, 1, 0, 1, 0],\
                                     [0, 0, 1, 0, 0, 1, 0, 1, 0, 0, 1, 1],\
                                     [0, 0, 0, 1, 0, 1, 0, 1, 1, 0, 0, 0],\
                                     [0, 0, 0, 0, 1, 1, 0, 0, 1, 1, 0, 1],\
                                     [0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 0, 1]])
            sage: C = codes.LinearCode(M)
            sage: D = C.decoder('InformationSet', (2,2))
            sage: c = C.random_element()
            sage: Chan = channels.StaticErrorRateChannel(C.ambient_space(), 2)
            sage: y = Chan(c)
            sage: c_out = D._lee_brickell_algorithm(y, (2, 2), 2)
            sage: (y - c).hamming_weight() == 2
            True
        """
        import itertools
        from sage.misc.prandom import sample
        C = self.code()
        n, k = C.length(), C.dimension()
        F = C.base_ring()
        G = C.generator_matrix()
        Fstar = F.list()[1:]
        while True:
            # step 1.
            I = sample(range(n), k)
            Gi = G.matrix_from_columns(I)
            try:
                Gi_inv = Gi.inverse()
            except ZeroDivisionError:
                # I was not an information set
                continue
            Gt = Gi_inv * G
            #step 2.
            y = r - vector([r[i] for i in I]) * Gt
            g = Gt.rows()
            #step 3.
            for pi in range(p+1):
                for A in itertools.combinations(range(k), pi):
                    for m in itertools.product(Fstar, repeat=pi):
                        e = y - sum(m[i]*g[A[i]] for i in range(pi))
                        errs = e.hamming_weight()
                        if  errs >= tau[0] and errs <= tau[1]:
                            return r - e

    def decode_to_code(self, r):
        r"""
        Decodes a received word with respect to the associated code of this decoder.

        WARNING:

        If there is no codeword within the decoding radius of this decoder, this
        method will never terminate.

        INPUT:

        - ``r`` -- a vector in the ambient space of :meth:`decoder.Decoder.code`.

        OUTPUT: a codeword of :meth:`decoder.Decoder.code`.

        EXAMPLES::

            sage: M = matrix(GF(2), [[1, 0, 0, 0, 0, 1, 0, 0, 1, 0, 1, 0],\
                                     [0, 0, 1, 0, 0, 1, 0, 1, 0, 0, 1, 1],\
                                     [0, 0, 0, 1, 0, 1, 0, 1, 1, 0, 0, 0],\
                                     [0, 0, 0, 0, 1, 1, 0, 0, 1, 1, 0, 1],\
                                     [0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 0, 1]])
            sage: C = codes.LinearCode(M)
            sage: D = C.decoder('InformationSet', 2, algorithm="LeeBrickell")
            sage: c = C.random_element()
            sage: Chan = channels.StaticErrorRateChannel(C.ambient_space(), 2)
            sage: y = Chan(c)
            sage: D.decode_to_code(y) in C
            True

        Information-set decoding a non-binary code::

            sage: C = codes.GolayCode(GF(3)); C
            [12, 6, 6] Extended Golay code over GF(3)
            sage: D = C.decoder('InformationSet', 2, algorithm="LeeBrickell")
            sage: c = C.random_element()
            sage: Chan = channels.StaticErrorRateChannel(C.ambient_space(), 2)
            sage: y = Chan(c)
            sage: D.decode_to_code(y) in C
            True
        """
        C = self.code()
        if r in C:
            return r
        return self._lee_brickell_algorithm(r, self.decoding_interval(),
                                                self.search_size())

    def _calibrate_search_size(self):
        r"""
        Run some test computations to estimate the optimal search size.

        We should simply choose `p` such that the average expected time is
        minimal. The algorithm succeeds when it chooses an information set with
        at least `k - p` correct positions, where `k` is the dimension of the
        code and `p` the search size. The expected number of trials we need
        before this occurs is::

            binom{n}{k}/(\rho \sum_{i=0}^p \binom{n-\tau}{k-i}\binom{\tau}{i})

        Here `\rho` is the fraction of `k` subsets of indices which are
        information sets. If `T` is the average time for steps 1 and 2
        (including selecting `I` until an information set is found), while `P(i)`
        is the time for the body of the ``for``-loop in step 3 for `m` of weight
        `i`, then each information set trial takes roughly time `T +
        \sum_{i=0}^{p} P(i) \binom{k}{i} (q-1)^i`, where `\GF{q}` is the base
        field.

        The values `T` and `P` are here estimated by running a few test
        computations similar to those done by the decoding algorithm.
        We don't explicitly estimate `rho`.

        OUTPUT: A list of floats representing the estimated decoding times for
        each possible value of `p`.

        EXAMPLES::

            sage: from sage.coding.isd_decoding import LinearCodeISD_LeeBrickell
            sage: C = codes.GolayCode(GF(2))
            sage: D = LinearCodeISD_LeeBrickell(C, 3); D
            Information set decoder (LeeBrickell) for [24, 12, 8] Extended Golay code over GF(2) decoding up to 3 errors
            sage: D._calibrate_search_size() #random
            [0.002337092727272496, 0.0008162108571427874, 0.0010929299999999607]
        """
        from sage.misc.prandom import sample
        from sage.stats.basic_stats import mean
        from sage.modules.free_module_element import random_vector
        from sage.matrix.special import random_matrix
        from sage.misc.prandom import randint
        import time
        C = self.code()
        G = C.generator_matrix()
        n, k = C.length(), C.dimension()
        tau = self.decoding_radius()
        F = C.base_ring()
        q = F.cardinality()
        Fstar = F.list()[1:]
        def time_information_set_steps():
            before = time.clock()
            while True:
                I = sample(range(n), k)
                Gi = G.matrix_from_columns(I)
                try:
                    Gi_inv = Gi.inverse()
                except ZeroDivisionError:
                    continue
                return time.clock() - before
        def time_search_loop(p):
            y = random_vector(F, n)
            g = random_matrix(F, p, n).rows()
            scalars = [  [ Fstar[randint(0,q-2)] for i in range(p) ]
                             for s in range(100) ]
            before = time.clock()
            for m in scalars:
                e = y - sum(m[i]*g[i] for i in range(p))
                errs = e.hamming_weight()
            return (time.clock() - before)/100.
        T = mean([ time_information_set_steps() for s in range(5) ])
        P = [ time_search_loop(p) for p in range(tau+1) ]

        estimates = []
        for p in range(tau+1):
            iters = 1.* binomial(n, k)/ \
                sum( binomial(n-tau, k-i)*binomial(tau,i) for i in range(p+1) )
            estimate = iters*(T + \
                sum(P[pi] * (q-1)**pi * binomial(k, pi) for pi in range(p+1) ))
            estimates.append(estimate)
        return estimates

    def search_size(self):
        r"""
        The search-size parameter for this Lee--Brickel information-set decoder.

        EXAMPLES::

            sage: C = codes.GolayCode(GF(2))
            sage: D = C.decoder("InformationSet", 2, algorithm="LeeBrickell", search_size=2)
            sage: D.search_size()
            2

        If not set, calibration will determine a sensible value::

            sage: C = codes.GolayCode(GF(2))
            sage: D = C.decoder("InformationSet", 2, algorithm="LeeBrickell")
            sage: D.search_size()
            2
        """
        if not hasattr(self, "_search_size"):
            estimates = self._calibrate_search_size()
            self._search_size = 0
            for p in range(1, len(estimates)):
                if estimates[p] < estimates[self._search_size]:
                    self._search_size = p
        return self._search_size


def LinearCodeInformationSetDecoder(code, number_errors, algorithm=None, **kwargs):
    r"""
    Information-set decoder for any linear code.

    This is a factory function that selects an appropriate information set
    decoding (ISD) algorithm.

    Information-set decoding is a probabilistic decoding strategy that
    essentially tries to guess `k` correct positions in the received word,
    where `k` is the dimension of the code. A codeword agreeing with the
    received word on the guessed position can easily be computed, and their
    difference is one possible error vector. A "correct" guess is assumed when
    this error vector has low Hamming weight.

    This simple algorithm is not very efficient in itself, but there are numerous
    refinements to the strategy. This factory function allows access to the ones
    implemented in Sage, and selects an appropriate one if the user does not
    specify a preference.

    The ISD strategy requires choosing how many errors is deemed acceptable. One
    choice could be `d/2`, where `d` is the minimum distance of the code, but
    sometimes `d` is not known, or sometimes more errors are expected. If one
    chooses anything above `d/2`, the algorithm does not guarantee to return a
    nearest codeword.

    WARNING::

        If there is no codeword within the specified decoding distance, then the
        decoding algorithm is not promised to terminate.

    INPUT:

    - ``code`` -- A linear code for which to decode.

    - ``number_errors`` -- an integer, the maximal number of errors to accept as
      correct decoding. An interval can also be specified by giving a pair of
      integers, where both end values are taken to be in the interval.

    - ``algorithm`` -- (optional) the ISD algorithm to employ. If this is not
      set, an appropriate one will be chosen.

    - ``**kwargs`` -- (optional) any number of named arguments passed on to the
      ISD algorithm. Such are usually not required, and they can only be set if
      ``algorithm`` is set to a specific algorithm. See the documentation for
      each individual ISD algorithm class for information on any named arguments
      they may accept. The easiest way to access this documentation is to first
      construct the decoder without passing any named arguments, and then using
      the `?` help on the constructed object.

    EXAMPLES::

    The principal way to access this function is through the
    :meth:`sage.code.linear_code.AbstractLinearCode.decoder` method:

        sage: C = codes.GolayCode(GF(3))
        sage: D = C.decoder("InformationSet", 2); D
        Information set decoder (LeeBrickell) for [12, 6, 6] Extended Golay code over GF(3) decoding up to 2 errors

    You can specify which algorithm you wish to use, and you should do so in
    order to pass special parameters to it::

        sage: C = codes.GolayCode(GF(3))
        sage: D2 = C.decoder("InformationSet", 2, algorithm="LeeBrickell", search_size=2); D2
        Information set decoder (LeeBrickell) for [12, 6, 6] Extended Golay code over GF(3) decoding up to 2 errors

    If you specify an algorithm which is not known, you get a friendly error message::

        sage: from sage.coding.isd_decoding import LinearCodeInformationSetDecoder
        sage: D = LinearCodeInformationSetDecoder(C, 2, algorithm="NoSuchThing"); D
        Traceback (most recent call last):
        ...
        ValueError: Unknown ISD algorithm 'NoSuchThing'. The known algorithms are ['LeeBrickell'].

    There are two other ways to access this function::

        sage: D = codes.decoders.LinearCodeInformationSetDecoder(C, 2); D
        Information set decoder (LeeBrickell) for [12, 6, 6] Extended Golay code over GF(3) decoding up to 2 errors

        sage: from sage.coding.isd_decoding import LinearCodeInformationSetDecoder
        sage: D = LinearCodeInformationSetDecoder(C, 2); D
        Information set decoder (LeeBrickell) for [12, 6, 6] Extended Golay code over GF(3) decoding up to 2 errors
    """
    if algorithm is None:
        if kwargs:
            raise ValueError("Additional arguments to an information-set decoder"
                            " algorithm are only allowed if a specific"
                            " algorithm is selected by setting the algorithm"
                            " keyword")
        algorithm = "LeeBrickell"
    knowns = {
        "LeeBrickell": LinearCodeISD_LeeBrickell
        }
    if algorithm in knowns:
        return knowns[algorithm](code, number_errors, **kwargs)
    else:
        raise ValueError("Unknown ISD algorithm '{}'."
                         " The known algorithms are {}."\
                         .format(algorithm, knowns.keys()))

AbstractInformationSetDecoder._decoder_type = {"hard-decision",
    "probabilistic", "not-always-closest", "bounded-distance", "might-fail"}
