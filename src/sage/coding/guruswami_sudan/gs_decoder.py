r"""
Guruswami-Sudan decoder for Generalized Reed-Solomon codes

REFERENCES:

    .. [GS99] Venkatesan Guruswami and Madhu Sudan, Improved Decoding of
       Reed-Solomon Codes and Algebraic-Geometric Codes, 1999

    .. [N13] Johan S. R. Nielsen, List Decoding of Algebraic Codes, 2013
"""

#*****************************************************************************
#       Copyright (C) 2015 David Lucas <david.lucas@inria.fr>
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 2 of the License, or
# (at your option) any later version.
#                  http://www.gnu.org/licenses/
#*****************************************************************************

from sage.coding.grs import GeneralizedReedSolomonCode
from sage.modules.free_module_element import vector
from sage.coding.decoder import Decoder
from sage.coding.guruswami_sudan.interpolation import construct_Q_linalg
from sage.coding.guruswami_sudan.rootfinding import rootfind_roth_ruckenstein
from sage.coding.guruswami_sudan.utils import (list_decoding_range,
                                               gilt,
                                               solve_degree2_to_integer_range,
                                               find_minimal_satisfiable)
from sage.functions.other import binomial, floor, sqrt

class GRSGuruswamiSudanDecoder(Decoder):
    r"""
    A decoder based on the Guruswami-Sudan list-decoding algorithm.

    One can read [GS99]_ to learn more about this algorithm.

    INPUT:

    - ``code`` -- A code associated to this decoder.

    - ``tau`` -- (default: ``None``) an integer, number of errors one wants Guruswami-Sudan algorithm
      to correct.

    - ``parameters`` -- (default: ``None``) a pair of integers, where:
        - the first integer is the multiplicity parameter of Guruswami-Sudan algorithm and
        - the second integer is the list size parameter.

    - ``interpolation_alg`` -- (default: ``None``) the name of the interpolation algorithm that will be used,
      or a the method which performs the interpolation. See NOTES section for details on signature.
      One can use the following names:

        * ``LinearAlgebra`` -- uses a linear system solver.

    - ``root_finder`` -- (default: ``None``) the name of the rootfinding algorithm that will be used,
      or a the method which performs the rootfinding. See NOTES section for details on signature.
      One can use the following names:

        * ``RothRuckenstein`` -- uses Roth-Ruckenstein algorithm.

    .. NOTE::

        One has to provide either ``C`` or ``parameters``. If none is
        given, an exception will be raised.

        If one wants to provide a method as ``rootfinder``, its signature
        has to be: ``my_rootfinder(Q, maxd=default_value, precision=default_value)``.
        See :meth:`sage.coding.guruswami_sudan.rootfinding.rootfind_roth_ruckenstein`
        for an example.

        If one wants to provide a method as ``interpolation_alg``, its signature
        has to be: ``my_inter(interpolation_points, tau, s_and_l, wy)``.
        See :meth:`sage.coding.guruswami_sudan.interpolation.construct_Q_linalg`
        for an example.

    EXAMPLES::

        sage: C = codes.GeneralizedReedSolomonCode(GF(251).list()[:250], 70)
        sage: D = codes.decoders.GRSGuruswamiSudanDecoder(C, tau = 97)
        sage: D
        Guruswami-Sudan decoder for [250, 70, 181] Generalized Reed-Solomon Code over Finite Field of size 251

    One can specify ``s`` and ``l`` instead of ``tau``::

        sage: D = codes.decoders.GRSGuruswamiSudanDecoder(C, parameters = (1,2))
        sage: D
        Guruswami-Sudan decoder for [250, 70, 181] Generalized Reed-Solomon Code over Finite Field of size 251

    One can pass a method as ``root_finder`` (works also for ``interpolation_alg``):

        sage: from sage.coding.guruswami_sudan.rootfinding import rootfind_roth_ruckenstein
        sage: rf = rootfind_roth_ruckenstein
        sage: D = codes.decoders.GRSGuruswamiSudanDecoder(C, parameters = (1,2), root_finder = rf)
        sage: D
        Guruswami-Sudan decoder for [250, 70, 181] Generalized Reed-Solomon Code over Finite Field of size 251
    """

    ####################### static methods ###############################
    IMPOSSIBLE_PARAMS = "Impossible parameters for the Guruswami-Sudan algorithm"

    @staticmethod
    def best_parameters_given_tau(tau, C = None, n_k = None):
        r"""
        Returns the best ``s`` and ``l`` possible according to input parameters.

        INPUT:

        - ``tau`` -- an integer, number of errrors one expects Guruswami-Sudan algorithm
          to correct
        - ``C`` -- (default: ``None``) a :class:`GeneralizedReedSolomonCode`
        - ``n_k`` -- (default: ``None``) a tuple of integers, respectively the
          length and the dimension of the :class:`GeneralizedReedSolomonCode`

        OUTPUT:

        - ``(s, l)`` -- a couple of integers, where:
            - ``s`` is the multiplicity parameter of Guruswami-Sudan algorithm and
            - ``l`` is the list size parameter

        .. NOTE::

            One has to provide either ``C`` or ``(n, k)``. If none or both are
            given, an exception will be raised.

        EXAMPLES::

            sage: tau, n, k = 97, 250, 70
            sage: codes.decoders.GRSGuruswamiSudanDecoder.best_parameters_given_tau(tau, n_k = (n, k))
            (1, 2)

        Another one with a bigger tau::

            sage: tau, n, k = 118, 250, 70
            sage: codes.decoders.GRSGuruswamiSudanDecoder.best_parameters_given_tau(tau, n_k = (n, k))
            (47, 89)
        """
        if C is not None and n_k is not None:
            raise ValueError("Please provide only the code or its length and dimension")
        elif C is None and n_k is None:
            raise ValueError("Please provide either the code or its length and dimension")
        elif C is not None:
            n, k = C.length(), C.dimension()
        elif n_k is not None and not isinstance(n_k, tuple):
            raise ValueError("n_k has to be a tuple")
        elif n_k is not None:
            n, k = n_k[0], n_k[1]
        (firsts, firstl) = GRSGuruswamiSudanDecoder._suitable_parameters_from_tau(tau, n_k = (n, k))
        def try_l(l):
            (mins,maxs) = solve_degree2_to_integer_range(n, n-2*(l+1)*(n-tau), (k-1)*l*(l+1))
            if maxs > 0 and maxs >= mins:
                return max(1, mins)
            else:
                return None
        l = find_minimal_satisfiable(try_l, firstl)
        s = try_l(l)
        assert GRSGuruswamiSudanDecoder.gs_satisfactory(tau, s, l, n_k = (n, k)) , IMPOSSIBLE_PARAMS
        return (s, l)

    @staticmethod
    def guruswami_sudan_decoding_radius(C = None, n_k = None, l = None, s = None):
        r"""
        Returns the maximal decoding radius of the Guruswami-Sudan decoder.

        If ``s`` is set but ``l`` is not it will return the best decoding radius using this ``s``
        alongside with the required ``l``, same for ``l``.

        INPUT:

        - ``C`` -- (default: ``None``) a :class:`GeneralizedReedSolomonCode`
        - ``n_k`` -- (default: ``None``) a tuple of integers, respectively the
          length and the dimension of the :class:`GeneralizedReedSolomonCode`
        - ``s`` -- (default: ``None``) an integer, the multiplicity parameter of Guruswami-Sudan algorithm
        - ``l`` -- (default: ``None``) an integer, the list size parameter

        .. NOTE::

            One has to provide either ``C`` or ``n_k``. If none or both are
            given, an exception will be raised.

        EXAMPLES::

            sage: n, k = 250, 70
            sage: codes.decoders.GRSGuruswamiSudanDecoder.guruswami_sudan_decoding_radius(n_k = (n, k))
            (118, (47, 89))
        """
        if C is not None and n_k is not None:
            raise ValueError("Please provide only the code or its length and dimension")
        elif C is None and n_k is None:
            raise ValueError("Please provide either the code or its length and dimension")
        elif C is not None:
            n, k = C.length(), C.dimension()
        elif n_k is not None and not isinstance(n_k, tuple):
            raise ValueError("n_k has to be a tuple")
        elif n_k is not None:
            n, k = n_k[0], n_k[1]

        def get_tau(s,l):
            if s<=0 or l<=0:
                return -1
            return gilt(n - n/2*(s+1)/(l+1) - (k-1)/2*l/s)
        if l==None and s==None:
            tau = list_decoding_range(n,n-k+1)[1]
            return (tau, GRSGuruswamiSudanDecoder.best_parameters_given_tau(tau, n_k = (n, k)))
        if l!=None and s!=None:
            return (get_tau(s,l), (s,l))
        if s!= None:
            # maximising tau under condition
            # n*(s+1 choose 2) < (ell+1)*s*(n-tau) - (ell+1 choose 2)*(k-1)
            # knowing n and s, we can just minimise
            # ( n*(s+1 choose 2) + (ell+1 choose 2)*(k-1) )/(ell+1)
            # Differentiating and setting to zero yields ell best choice:
            lmax = sqrt(n*s*(s+1)/(k-1)) - 1
            #the best integral value will be
            (l,tau) = find_integral_max(lmax, lambda l: get_tau(s,l))
            assert GRSGuruswamiSudanDecoder.gs_satisfactory(tau,s,l, n_k = (n, k)), IMPOSSIBLE_PARAMS
            #Note that we have not proven that this ell is minimial in integral
            #sense! It just seems that this most often happens
            return (tau,(s,l))
        if l!= None:
            # Acquired similarly to when restricting s
            smax = sqrt((k-1)/n*l*(l+1))
            (s,tau) = find_integral_max(smax, lambda s: get_tau(s,l))
            assert GRSGuruswamiSudanDecoder.gs_satisfactory(tau,s,l, n_k = (n, k)), IMPOSSIBLE_PARAMS
            return (get_tau(s,l), (s,l))

    @staticmethod
    def _suitable_parameters_from_tau(tau, C = None, n_k = None):
        r"""
        Returns ``s`` and ``l`` according to input parameters.

        These parameters are not guaranteed to be the best ones possible
        for the provided ``tau``, but are a good approximation of the best ones.

        See [N13]_ pages 53-54, proposition 3.11 for details.

        INPUT:

        - ``tau`` -- an integer, number of errrors one expects Guruswami-Sudan algorithm
          to correct
        - ``C`` -- (default: ``None``) a :class:`GeneralizedReedSolomonCode`
        - ``n_k`` -- (default: ``None``) a tuple of integers, respectively the
          length and the dimension of the :class:`GeneralizedReedSolomonCode`

        OUTPUT:

        - ``(s, l)`` -- a couple of integers, where:
            - ``s`` is the multiplicity parameter of Guruswami-Sudan algorithm and
            - ``l`` is the list size parameter

        ..NOTE::

            One has to provide either ``C`` or ``(n, k)``. If none or both are
            given, an exception will be raised.

        EXAMPLES::

            sage: tau = 97
            sage: n, k = 250, 70
            sage: codes.decoders.GRSGuruswamiSudanDecoder._suitable_parameters_from_tau(tau, n_k = (n, k))
            (2, 3)

        Same one with a GRS code::

            sage: C = codes.GeneralizedReedSolomonCode(GF(251).list()[:250], 70)
            sage: codes.decoders.GRSGuruswamiSudanDecoder._suitable_parameters_from_tau(tau, C = C)
            (2, 3)

        Another one with a bigger ``tau``::

            sage: codes.decoders.GRSGuruswamiSudanDecoder._suitable_parameters_from_tau(118, C = C)
            (47, 89)
        """
        if C is not None and n_k is not None:
            raise ValueError("Please provide only the code or its length and dimension")
        elif C is None and n_k is None:
            raise ValueError("Please provide either the code or its length and dimension")
        elif C is not None:
            n, k = C.length(), C.dimension()
        elif n_k is not None and not isinstance(n_k, tuple):
            raise ValueError("n_k has to be a tuple")
        elif n_k is not None:
            n, k = n_k[0], n_k[1]

        w = k - 1
        atau = n - tau
        smin = tau * w / (atau ** 2 - n * w)
        s = floor(1 + smin)
        D = (s - smin) * (atau ** 2 - n * w) * s + (w**2) /4
        l = floor(atau / w * s + 0.5 - sqrt(D)/w)
        assert GRSGuruswamiSudanDecoder.gs_satisfactory(tau,s,l, n_k = (n, k)) , IMPOSSIBLE_PARAMS
        return (s, l)

    @staticmethod
    def gs_satisfactory(tau, s, l, C = None, n_k = None):
        r"""
        Returns whether input parameters satisfy the governing equation of
        Guruswami-Sudan.

        See [N13]_ page 49, definition 3.3 and proposition 3.4 for details.

        INPUT:

        - ``tau`` -- an integer, number of errrors one expects Guruswami-Sudan algorithm
          to correct
        - ``s`` -- an integer, multiplicity parameter of Guruswami-Sudan algorithm
        - ``l`` -- an integer, list size parameter
        - ``C`` -- (default: ``None``) a :class:`GeneralizedReedSolomonCode`
        - ``n_k`` -- (default: ``None``) a tuple of integers, respectively the
          length and the dimension of the :class:`GeneralizedReedSolomonCode`

        ..NOTE::

            One has to provide either ``C`` or ``(n, k)``. If none or both are
            given, an exception will be raised.

        EXAMPLES::

            sage: tau, s, l = 97, 1, 2
            sage: n, k = 250, 70
            sage: codes.decoders.GRSGuruswamiSudanDecoder.gs_satisfactory(tau, s, l, n_k = (n, k))
            True

        One can also pass a GRS code::

            sage: C = codes.GeneralizedReedSolomonCode(GF(251).list()[:250], 70)
            sage: codes.decoders.GRSGuruswamiSudanDecoder.gs_satisfactory(tau, s, l, C = C)
            True

        Another example where ``s`` and ``l`` does not satisfy the equation::

            sage: tau, s, l = 118, 47, 80
            sage: codes.decoders.GRSGuruswamiSudanDecoder.gs_satisfactory(tau, s, l, n_k = (n, k))
            False

        If one provides both ``C`` and ``n_k`` an exception is returned::

            sage: tau, s, l = 97, 1, 2
            sage: n, k = 250, 70
            sage: C = codes.GeneralizedReedSolomonCode(GF(251).list()[:250], 70)
            sage: codes.decoders.GRSGuruswamiSudanDecoder.gs_satisfactory(tau, s, l, C = C, n_k = (n, k))
            Traceback (most recent call last):
            ...
            ValueError: Please provide only the code or its length and dimension

        Same if one provides none of these::

            sage: codes.decoders.GRSGuruswamiSudanDecoder.gs_satisfactory(tau, s, l)
            Traceback (most recent call last):
            ...
            ValueError: Please provide either the code or its length and dimension
        """
        if C is not None and n_k is not None:
            raise ValueError("Please provide only the code or its length and dimension")
        elif C is None and n_k is None:
            raise ValueError("Please provide either the code or its length and dimension")
        elif C is not None:
            n, k = C.length(), C.dimension()
        elif n_k is not None and not isinstance(n_k, tuple):
            raise ValueError("n_k has to be a tuple")
        elif n_k is not None:
            n, k = n_k[0], n_k[1]
        return l > 0 and s > 0 and n * s * (s+1) < (l+1) * (2*s*(n-tau) - (k-1) * l)



    ####################### decoder itself ###############################
    def __init__(self, code, tau = None, parameters = None, interpolation_alg = None, root_finder = None):
        r"""
        TESTS:

        If neither ``tau`` nor ``parameters`` is given, an exception is returned::

            sage: C = codes.GeneralizedReedSolomonCode(GF(251).list()[:250], 70)
            sage: D = codes.decoders.GRSGuruswamiSudanDecoder(C)
            Traceback (most recent call last):
            ...
            ValueError: Specify either tau or parameters

        """
        n, k = code.length(), code.dimension()
        if tau and parameters:
            assert GRSGuruswamiSudanDecoder.gs_satisfactory(tau, parameters[0], parameters[1], C = code), IMPOSSIBLE_PARAMS
            self._tau, self._s, self._ell = tau, parameters[0], parameters[1]
        elif tau:
            self._tau = tau
            self._s, self._ell = GRSGuruswamiSudanDecoder.best_parameters_given_tau(tau, n_k = (n, k))
        elif parameters:
            self._s = parameters[0]
            self._ell = parameters[1]
            (self._tau,_) = GRSGuruswamiSudanDecoder.guruswami_sudan_decoding_radius(C = code, s=self._s, l=self._ell)
        else:
            raise ValueError("Specify either tau or parameters")
        if hasattr(interpolation_alg, '__call__'):
            self.interpolation_alg = interpolation_alg
        elif interpolation_alg == None or interpolation_alg == "LinearAlgebra":
            self.interpolation_alg = construct_Q_linalg
        else:
            raise ValueError("Please provide a method or one of the allowed strings for interpolation_alg")
        if hasattr(root_finder, '__call__'):
            self.root_finder = root_finder
        elif root_finder == None or interpolation_alg == "RothRuckenstein":
            self.root_finder = rootfind_roth_ruckenstein
        else:
            raise ValueError("Please provide a method or one of the allowed strings for root_finder")
        super(GRSGuruswamiSudanDecoder, self).__init__(code, code.ambient_space(), "EvaluationPolynomial")

    def _repr_(self):
        r"""
        Returns a string representation of ``self``.

        EXAMPLES::

            sage: C = codes.GeneralizedReedSolomonCode(GF(251).list()[:250], 70)
            sage: D = C.decoder("GuruswamiSudan", tau = 97)
            sage: D
            Guruswami-Sudan decoder for [250, 70, 181] Generalized Reed-Solomon Code over Finite Field of size 251
        """
        return "Guruswami-Sudan decoder for %s" % self.code()

    def _latex_(self):
        r"""
        Returns a string representation of ``self``.

        EXAMPLES::

            sage: C = codes.GeneralizedReedSolomonCode(GF(251).list()[:250], 70)
            sage: D = C.decoder("GuruswamiSudan", tau = 97)
            sage: latex(D)
            \textnormal{Guruswami-Sudan decoder for } [250, 70, 181] \textnormal{ Generalized Reed-Solomon Code over } \Bold{F}_{251}
        """
        return "\\textnormal{Guruswami-Sudan decoder for } %s" % self.code()._latex_()

    def __eq__(self, other):
        r"""
        Tests equality between GRSGuruswamiSudanDecoder objects.

        EXAMPLES::

            sage: C = codes.GeneralizedReedSolomonCode(GF(251).list()[:250], 70)
            sage: D1 = C.decoder("GuruswamiSudan", tau = 97)
            sage: D2 = C.decoder("GuruswamiSudan", tau = 97)
            sage: D1.__eq__(D2)
            True
        """
        return isinstance(other, GRSGuruswamiSudanDecoder)\
                and self.code() == other.code()\
                and self.multiplicity() == other.multiplicity()\
                and self.list_size() == other.list_size()\
                and self.interpolation_algorithm() == other.interpolation_algorithm()\
                and self.rootfinding_algorithm() == other.rootfinding_algorithm()

    def __ne__(self, other):
        r"""
        Tests inequality between GRSGuruswamiSudanDecoder objects.

        EXAMPLES::

            sage: C = codes.GeneralizedReedSolomonCode(GF(251).list()[:250], 70)
            sage: D1 = C.decoder("GuruswamiSudan", tau = 97)
            sage: D2 = C.decoder("GuruswamiSudan", tau = 98)
            sage: D1.__ne__(D2)
            True
        """
        return not self.__eq__(other)

    def interpolation_algorithm(self):
        r"""
        Returns the interpolation algorithm that will be used.

        Remember that its signature has to be:
        ``my_inter(interpolation_points, tau, s_and_l, wy)``.
        See :meth:`sage.coding.guruswami_sudan.interpolation.construct_Q_linalg`
        for an example.

        EXAMPLES::

            sage: C = codes.GeneralizedReedSolomonCode(GF(251).list()[:250], 70)
            sage: D = C.decoder("GuruswamiSudan", tau = 97)
            sage: D.interpolation_algorithm() #random
            <function construct_Q_linalg at 0x7f9d55753500>
        """
        return self.interpolation_alg

    def rootfinding_algorithm(self):
        r"""
        Returns the rootfinding algorithm that will be used.

        Remember that its signature has to be:
        ``my_rootfinder(Q, maxd=default_value, precision=default_value)``.
        See :meth:`sage.coding.guruswami_sudan.rootfinding.rootfind_roth_ruckenstein`
        for an example.

        EXAMPLES::

            sage: C = codes.GeneralizedReedSolomonCode(GF(251).list()[:250], 70)
            sage: D = C.decoder("GuruswamiSudan", tau = 97)
            sage: D.rootfinding_algorithm() #random
            <function rootfind_roth_ruckenstein at 0x7fea00618848>
        """
        return self.root_finder

    def multiplicity(self):
        r"""
        Returns the multiplicity parameter of ``self``.

        EXAMPLES::

            sage: C = codes.GeneralizedReedSolomonCode(GF(251).list()[:250], 70)
            sage: D = C.decoder("GuruswamiSudan", tau = 97)
            sage: D.multiplicity()
            1
        """
        return self._s

    def list_size(self):
        r"""
        Returns the list size parameter of ``self``.

        EXAMPLES::

            sage: C = codes.GeneralizedReedSolomonCode(GF(251).list()[:250], 70)
            sage: D = C.decoder("GuruswamiSudan", tau = 97)
            sage: D.list_size()
            2
        """
        return self._ell

    def decode_to_message(self, r):
        r"""
        Decodes ``r`` to the message space of ``self``'s connected encoder.

        INPUT:

        - ``r`` -- a element of the input space of ``self``.

        EXAMPLES::

            sage: C = codes.GeneralizedReedSolomonCode(GF(17).list()[:15], 6)
            sage: D = codes.decoders.GRSGuruswamiSudanDecoder(C, tau=5)
            sage: F.<x> = GF(17)[]
            sage: m = 9*x^5 + 10*x^4 + 9*x^3 + 7*x^2 + 15*x + 2
            sage: c = D.connected_encoder().encode(m)
            sage: r = vector(GF(17), [3,1,4,2,14,1,0,4,13,12,1,16,1,13,15])
            sage: (c-r).hamming_weight()
            5
            sage: m in D.decode_to_message(r)
            True

        TESTS:

        If one has provided a method as a ``root_finder`` or a ``interpolation_alg`` which
        does not fits the allowed signature, an exception will be raised::

            sage: C = codes.GeneralizedReedSolomonCode(GF(17).list()[:15], 6)
            sage: D = codes.decoders.GRSGuruswamiSudanDecoder(C, tau=5, root_finder=next_prime)
            sage: F.<x> = GF(17)[]
            sage: m = 9*x^5 + 10*x^4 + 9*x^3 + 7*x^2 + 15*x + 2
            sage: c = D.connected_encoder().encode(m)
            sage: r = vector(GF(17), [3,1,4,2,14,1,0,4,13,12,1,16,1,13,15])
            sage: m in D.decode_to_message(r)
            Traceback (most recent call last):
            ...
            TypeError: next_prime() got an unexpected keyword argument 'maxd'
        """
        C = self.code()
        n,k,d,alphas,colmults, s, l = C.length(), C.dimension(), C.minimum_distance(),\
                C.evaluation_points(), C.column_multipliers(), self.multiplicity(), self.list_size()
        tau = self.decoding_radius()
        ## SETUP INTERPOLATION PROBLEM
        wy = k-1
        points = [ (alphas[i], r[i]/colmults[i]) for i in range(0,len(alphas)) ]
        ## SOLVE INTERPOLATION
        try:
            Q = self.interpolation_algorithm()(points, tau, (s,l), wy)
        except TypeError, e:
            raise e
        ## EXAMINE THE FACTORS AND CONVERT TO CODEWORDS
        try:
            polynomials = self.rootfinding_algorithm()(Q, maxd = None)
        except TypeError, e:
            raise e
        if not polynomials:
            return None
        return [f for f in polynomials]

    def decode_to_code(self, r):
        r"""
        Corrects the errors in ``r`` and returns a codeword.

        INPUT:

        - ``r`` -- a element of the input space of ``self``.

        EXAMPLES::

            sage: C = codes.GeneralizedReedSolomonCode(GF(17).list()[:15], 6)
            sage: D = codes.decoders.GRSGuruswamiSudanDecoder(C, tau=5)
            sage: c = vector(GF(17), [2,1,2,1,14,1,11,4,13,0,1,16,1,13,15])
            sage: c in C
            True
            sage: r = vector(GF(17), [3,1,4,2,14,1,0,4,13,12,1,16,1,13,15])
            sage: r in C
            False
            sage: c in D.decode_to_code(r)
            True
        """
        return [ self.connected_encoder().encode(i) for i in self.decode_to_message(r) ]

    def decoding_radius(self):
        r"""
        Returns the maximal number of errors that ``self`` is able to correct.

        EXAMPLES::

            sage: C = codes.GeneralizedReedSolomonCode(GF(251).list()[:250], 70)
            sage: D = C.decoder("GuruswamiSudan", tau = 97)
            sage: D.decoding_radius()
            97

        Another one where tau is not a given input::

            sage: C = codes.GeneralizedReedSolomonCode(GF(251).list()[:250], 70)
            sage: D = C.decoder("GuruswamiSudan", parameters = (2,4))
            sage: D.decoding_radius()
            105
        """
        return self._tau


####################### registration ###############################

GeneralizedReedSolomonCode._registered_decoders["GuruswamiSudan"] = GRSGuruswamiSudanDecoder
GRSGuruswamiSudanDecoder._decoder_type = {"list-decoder", "always-succeed", "hard-decision"}
