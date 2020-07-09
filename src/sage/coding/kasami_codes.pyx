# -*- coding: utf-8 -*-
r"""
Kasami code

This module implements a construction for the extended Kasami codes.
The "regular" Kasami codes are obtained from truncating the extended version.
The coset graphs of the Kasami codes are distance-regular.

In particular, the extended Kasami codes result in distance-regular graphs with intersection arrays:
    * `[q^{2j+1}, q^{2j+1} - 1, q^{2j+1} - q, q^{2j+1} - q^{2j} + 1; 1, q, q^{2j} -1, q^{2j+1}]`
    * `[q^2, q^2 - 1, q^2 - q, 1; 1, q, q^2 - 1, q^2]`

The Kasami codes result in distance-regular graphs with intersection arrays:
    * `[q^{2j+1} - 1, q^{2j+1} - q, q^{2j+1} - q^{2j} + 1; 1, q, q^{2j} -1]`
    * `[q^2 - 1, q^2 - q, 1; 1, q, q^2 - 1]`

REFERENCES:

- [BCN1989]_ p. 358 for a definition.

AUTHORS:

- Ivo Maffei (2020-07-09): initial version
"""

#*****************************************************************************
#       Copyright (C) 2020 Ivo Maffei <ivomaffei@gmail.com>
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 2 of the License, or
# (at your option) any later version.
#                  http://www.gnu.org/licenses/
#*****************************************************************************

from sage.rings.finite_rings.finite_field_constructor import GF
from sage.modules.free_module import VectorSpace
from sage.modules.free_module_element import vector
from sage.matrix.constructor import matrix
from sage.coding.linear_code import LinearCode
from sage.coding.linear_code import AbstractLinearCode, LinearCodeGeneratorMatrixEncoder
from sage.arith.misc import is_prime_power, gcd

class KasamiCode(AbstractLinearCode):
    r"""
    Representation of a Kasami Code.

    INPUT:

    - ``s,t`` -- (integer) the parameters of the Kasami code

    - ``extended`` -- (default: ``True``) if set to ``True``, creates an extended Kasami
      code.

    EXAMPLES::

        sage: codes.KasamiCode(16,4)
        (16, 4) Extended Kasami code

        sage: codes.KasamiCode(8,4)
        Traceback (most recent call last):
        ...
        ValueError: The parameters(=8,4) are invalid. Check the documentation

    The extended Kasami code is the extension of the Kasami code::

        sage: C = codes.KasamiCode(16, 4, extended=False)
        sage: Cext = C.extended_code()
        sage: D = codes.KasamiCode(16, 4, extended=True)
        sage: D.generator_matrix() == Cext.generator_matrix()
        True

    .. SEEALSO::

        :mod:`sage.coding.linear_code`.

    REFERENCES:

    For more information on Kasami codes and their use see [BCN1989]_.

    TESTS:

         sage: C1 = codes.KasamiCode(16, 4)
         sage: C2 = codes.KasamiCode(16, 4, extended=False)
         sage: C1.parameters() == C2.parameters()
         True
         sage: C1 == C2
         False
         sage: C1.minimum_distance() == C2.minimum_distance()+1 
         True

         sage: C = codes.KasamiCode(4,2)
         sage: C.dimension()
         0
         sage: C.generator_matrix()
         []

    """

    _registered_encoders = {}
    _registered_decoders = {}

    def __init__(self, s, t, extended=True):
        r"""
        Constructor for the ``KasamiCode`` class.

        TESTS::

            sage: codes.KasamiCode(64,8)
            (64, 8) Extended Kasami code

            sage: codes.KasamiCode(64,8, extended=False)
            (64, 8) Kasami code

            sage: codes.KasamiCode(3,5)
            Traceback (most recent call last):
            ...
            ValueError: The parameter t(=5) must be a power of 2
        """
        # Check validity of s and t
        (p,i) = is_prime_power(t,get_data=True)
        if p != 2:
            raise ValueError("The parameter t(={}) must be a power of 2".format(t))

        if s != t*t:
            # then we must have s=q^{2j+1} and t = q^m
            (p,k) = is_prime_power(s,get_data=True)
            if p != 2:
                raise ValueError("The parameter s(={}) must be a power of 2".format(s))

            # q= 2^l here l = gcd(k,i)
            l = gcd(i,k)
            q = 2**l
            m = i // l

            if (k//l) % 2 == 0:
                raise ValueError("The parameter s(={}) is invalid. Check the documentation".format(s))

            j = ((k//l) - 1) // 2

            # gcd(m,2*j+1) = gcd( i/l, k/l) = 1
            if m > j:
                raise ValueError("The parameters(={},{}) are invalid. Check the documentation".format(s,t))

        # s and t are valid!!!
        self._s = s
        self._t = t

        length = s-1
        if extended:
            length += 1
        
        super(KasamiCode, self).__init__(GF(2), length, "GeneratorMatrix", "Syndrome")

    def parameters(self):
        r"""
        Return the parameters `s,t` of ``self``.

        Examples::

            sage: C = codes.KasamiCode(16, 4, extended=True)
            sage: C.parameters()
            (16, 4)
            sage: D = codes.KasamiCode(16, 4, extended=False)
            sage: D.parameters()
            (16, 4)

            sage: C = codes.KasamiCode(8,2)
            sage: C.parameters()
            (8, 2)
        """
        return (self._s,self._t)

    def __eq__(self, other):
        r"""
        Test equality between Kasami Code objects.

        EXAMPLES::

            sage: C1 = codes.KasamiCode(8,2)
            sage: C2 = codes.KasamiCode(8,2)
            sage: C1.__eq__(C2)
            True
        """
        return isinstance(other, KasamiCode) \
                and self.parameters() == other.parameters() \
                and self.length() == other.length()

    def _repr_(self):
        r"""
        Return a string representation of ``self``.

        EXAMPLES::

            sage: codes.KasamiCode(4,2,extended=True)
            (4, 2) Extended Kasami code
        """
        ext = ""
        if self.length() == self._s:
            ext = " Extended"
        return "(%s, %s)%s Kasami code"\
                % (self._s, self._t, ext)

    def _latex_(self):
        r"""
        Return a latex representation of ``self``.

        EXAMPLES::

            sage: C = codes.KasamiCode(16,4)
            sage: latex(C)
            (16, 4) \textnormal{ Extended Kasami code }
        """
        ext = ""
        if self.length() == self._s:
            ext = "Extended"
        return "(%s, %s) \\textnormal{ %s Kasami code }"\
                % (self._s, self._t, ext)

    def generator_matrix(self):
        r"""
        Return a generator matrix of ``self``.

        EXAMPLES::

            sage: C = codes.KasamiCode(16, 4, extended=False)
            sage: C.generator_matrix()
            [1 0 0 0 0 0 0 0 0 1 0 0 1 1 1]
            [0 1 0 0 0 0 0 0 0 1 1 0 1 0 0]
            [0 0 1 0 0 0 0 0 0 0 1 1 0 1 0]
            [0 0 0 1 0 0 0 0 0 0 0 1 1 0 1]
            [0 0 0 0 1 0 0 0 0 0 0 0 1 1 0]
            [0 0 0 0 0 1 0 0 0 1 1 0 1 1 1]
            [0 0 0 0 0 0 1 0 0 0 1 1 0 1 1]
            [0 0 0 0 0 0 0 1 0 1 1 1 0 0 1]
            [0 0 0 0 0 0 0 0 1 1 0 1 0 0 0]
        """
        if self.length() == self._s:
            return _extended_Kasami_code(self._s,self._t)
        return _Kasami_code(self._s,self._t)

def _extended_Kasami_code(const int s, const int t):
    r"""
    Return a generator matrix for the  extended Kasami code with parameters `s,t`.

    The extended Kasami code with parameters `(s,t)` is defined as

    .. MATH::

        \{ v \in GF(2)^s \mid
        \sum_{a \in GF(s)} v_a =
        \sum_{a \in GF(s)} a v_a =
        \sum_{a \in GF(s)} a^{t+1} v_a = 0 \}

    The only valid parameters `s,t` are given by the below, where `q` is a power of 2:
        * `s = q^{2j+1}`, `t = q^m` with `m \leq j` and `\gcd(m,2j+1) = 1`
        * `s = q^2`, `t=q`

    INPUT:

    - ``s,t`` -- (integers); powers of 2; parameters of the code


    EXAMPLES::

        sage: C = codes.KasamiCode(16,4)
        sage: C.minimum_distance()
        4

        sage: codes.KasamiCode(8,4)
        Traceback (most recent call last):
        ...
        ValueError: The parameters(=8,4) are invalid. Check the documentation

    .. SEEALSO::

        :class:`sage.coding.kasami_codes.KasamiCode`.

    ALGORITHM:

    We generate spanning sets for the subspaces:
        * `\{v \in GF(2)^s \mid \sum_{a \in GF(s)} v_a = 0\}`
        * `\{v \in GF(2)^s \mid \sum_{a \in GF(s)} a v_a = 0\}`
        * `\{v \in GF(2)^s \mid \sum_{a \in GF(s)} a^{t+1} v_a = 0\}`

    Then we compute our codebook by taking the intersection of the above subspaces.

    TESTS::

        sage: C = codes.KasamiCode(4,2)
        sage: C.generator_matrix()
        []

        sage: C = codes.KasamiCode(8,2)
        sage: C.generator_matrix()
        [1 1 1 1 1 1 1 1]
        sage: C.minimum_distance()
        8

        sage: C = codes.KasamiCode(16,4)
        sage: C.minimum_distance()
        4

        sage: C = codes.KasamiCode(64,4)
        sage: C.minimum_distance()  # long time
        4
    """
    F2 = GF(2)
    V = VectorSpace(F2, s)
    elemsFs = [x for x in GF(s)]

    #we ensure that 0 is the first element of elemsFs
    if not elemsFs[0].is_zero():
        for i in range(s):
            if elemsFs[i].is_zero:
                a = elemsFs[0]
                elemsFs[0] = elemsFs[i]
                elemsFs[i] = a
                break

    FsToInt = { x : i for i,x in enumerate(elemsFs)}
    elemsFsT = [x**(t+1) for x in elemsFs]
    FsTToInt = { x: i for i,x in enumerate(elemsFsT)}

    e1 = [0]*s
    e1[0] = 1
    e1 = vector(F2,e1,immutable=True)

    W1_basis = []
    for i in range(s-1):
        v = [0]*s
        v[i] = 1
        v[s-1] = 1
        W1_basis.append(v)
    W1 = V.span(W1_basis) #W1 satisfies \sum v[i] = 0

    W2_basis = set([e1])#not really a basis...
    for i in range(1,s):#avoid x = 0
        x = elemsFs[i]
        for j in range(i+1,s):
            y = elemsFs[j]
            v = [0]*s
            v[i] = 1
            v[j] = 1
            v[ FsToInt[(x+y)] ] = 1
            v = vector(F2,v,immutable=True)
            W2_basis.add(v)
    W2 = V.span(W2_basis) #W2 satisfies \sum v[i]elemsFs[i] = 0


    W3_basis = set([e1]) #again not really a basis
    for i in range(1,s): #avoid x = 0^(t+1) = 0
        x = elemsFsT[i]
        for j in range(i+1,s):
            y = elemsFsT[j]
            v = [0]*s
            v[i] = 1
            v[j] = 1
            v[ FsTToInt[(x+y)] ] = 1
            v=vector(F2,v,immutable=True)
            W3_basis.add(v)
    W3 = V.span(W3_basis)

    W = W2.intersection(W3)
    codebook = W.intersection(W1)

    return codebook.basis_matrix()

def _Kasami_code(const int s, const int t):
    r"""
    Return the generator matrix of the Kasami code with parameters `s,t`.

    The Kasami code `(s,t)` is obtained from the extended
    Kasami code `(s,t)`, via truncation of all words.

    INPUT:

    - ``s,t`` -- (integers); powers of 2; parameters of the code

    EXAMPLES::

        sage: codes.KasamiCode(8, 2, extended=False)
        (8, 2) Kasami code

        sage: codes.KasamiCode(4, 2, extended=False)
        (4, 2) Kasami code

    .. SEEALSO::

        :class:`sage.coding.kasami_codes.KasamiCode`,
        :meth:`sage.coding.kasami_codes._extended_Kasami_code`.

    TESTS::

        sage: C = codes.KasamiCode(8, 2, extended=False)
        sage: C.generator_matrix()
        [1 1 1 1 1 1 1]
        sage: C.minimum_distance()
        7

        sage: C = codes.KasamiCode(4, 2, extended=False)
        sage: C.generator_matrix()
        []

        sage: C = codes.KasamiCode(16, 4, extended=False)
        sage: C.minimum_distance()
        3

        sage: C = codes.KasamiCode(64,4, extended=False)
        sage: C.minimum_distance()  # long time
        3
    """
    M = _extended_Kasami_code(s,t)
    newM = [v[:-1] for v in M]
    
    return matrix(GF(2), newM)



KasamiCode._registered_encoders["GeneratorMatrix"] = LinearCodeGeneratorMatrixEncoder
