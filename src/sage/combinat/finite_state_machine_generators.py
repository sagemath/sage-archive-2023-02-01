r"""
Common Transducers (Finite State Machines Generators)

Transducers in Sage can be built through the ``transducers``
object. It contains generators for common finite state machines. For example,

::

    sage: I = transducers.Identity([0, 1, 2])

generates an identity transducer on the alphabet `\{0, 1, 2\}`.

To construct transducers manually, you can use the class
:class:`Transducer`. See :mod:`~sage.combinat.finite_state_machine`
for more details and a lot of examples.

**Transducers**

.. csv-table::
    :class: contentstable
    :widths: 30, 70
    :delim: |

    :meth:`~TransducerGenerators.Identity` | Returns a transducer realizing the identity map.
    :meth:`~TransducerGenerators.abs` | Returns a transducer realizing absolute value.
    :meth:`~TransducerGenerators.operator` | Returns a transducer realizing a binary operation.
    :meth:`~TransducerGenerators.add` | Returns a transducer realizing addition.
    :meth:`~TransducerGenerators.sub` | Returns a transducer realizing subtraction.
    :meth:`~TransducerGenerators.CountSubblockOccurrences` | Returns a transducer counting the occurrences of a subblock.
    :meth:`~TransducerGenerators.weight` | Returns a transducer realizing the Hamming weight
    :meth:`~TransducerGenerators.GrayCode` | Returns a transducer realizing binary Gray code.

AUTHORS:

- Clemens Heuberger (2014-04-07): initial version
- Sara Kropf (2014-04-10): some changes in TransducerGenerator
- Daniel Krenn (2014-04-15): improved common docstring during review
- Sara Kropf (2014-04-29): weight transducer

ACKNOWLEDGEMENT:

- Daniel Krenn, Clemens Heuberger and Sara Kropf are supported by the
  Austrian Science Fund (FWF): P 24644-N26.

Functions and methods
---------------------

"""
#*****************************************************************************
#       Copyright (C) 2014 Clemens Heuberger <clemens.heuberger@aau.at>
#                     2014 Sara Kropf <sara.kropf@aau.at>
#                     2014 Daniel Krenn <devel@danielkrenn.at>
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#                  http://www.gnu.org/licenses/
#*****************************************************************************

from sage.combinat.finite_state_machine import Transducer
from sage.rings.integer_ring import ZZ

class TransducerGenerators(object):
    r"""
    A class consisting of constructors for several common transducers.

    A list of all transducers in this database is available via tab
    completion. Type "``transducers.``" and then hit tab to see which
    transducers are available.

    The transducers currently in this class include:

    - :meth:`~Identity`
    - :meth:`~abs`
    - :meth:`~TransducerGenerators.operator`
    - :meth:`~add`
    - :meth:`~sub`
    - :meth:`~CountSubblockOccurrences`
    - :meth:`~GrayCode`

    """

    def Identity(self, input_alphabet):
        """
        Returns the identity transducer realizing the identity map.

        INPUT:

        - ``input_alphabet`` -- a list or other iterable.

        OUTPUT:

        A transducer mapping each word over ``input_alphabet`` to
        itself.

        EXAMPLES::

            sage: T = transducers.Identity([0, 1])
            sage: sorted(T.transitions())
            [Transition from 0 to 0: 0|0,
             Transition from 0 to 0: 1|1]
            sage: T.initial_states()
            [0]
            sage: T.final_states()
            [0]
            sage: T.input_alphabet
            [0, 1]
            sage: T.output_alphabet
            [0, 1]
            sage: sage.combinat.finite_state_machine.FSMOldProcessOutput = False
            sage: T([0, 1, 0, 1, 1])
            [0, 1, 0, 1, 1]

        """
        return Transducer(
            [(0, 0, d, d) for d in input_alphabet],
            input_alphabet=input_alphabet,
            output_alphabet=input_alphabet,
            initial_states=[0],
            final_states=[0])

    def CountSubblockOccurrences(self, block, input_alphabet):
        """
        Returns a transducer counting the number of (possibly
        overlapping) occurrences of a block in the input.

        INPUT:

        - ``block`` -- a list (or other iterable) of letters.

        - ``input_alphabet`` -- a list or other iterable.

        OUTPUT:

        A transducer counting (in unary) the number of occurrences of the given
        block in the input.  Overlapping occurrences are counted several
        times.

        Denoting the block by `b_0\ldots b_{k-1}`, the input word by
        `i_0\ldots i_L` and the output word by `o_0\ldots o_L`, we
        have `o_j = 1` if and only if `i_{j-k+1}\ldots i_{j} = b_0\ldots
        b_{k-1}`. Otherwise, `o_j = 0`.

        EXAMPLES:

        #.  Counting the number of ``10`` blocks over the alphabet
            ``[0, 1]``::

                sage: T = transducers.CountSubblockOccurrences(
                ....:     [1, 0],
                ....:     [0, 1])
                sage: sorted(T.transitions())
                [Transition from () to (): 0|0,
                 Transition from () to (1,): 1|0,
                 Transition from (1,) to (): 0|1,
                 Transition from (1,) to (1,): 1|0]
                sage: T.input_alphabet
                [0, 1]
                sage: T.output_alphabet
                [0, 1]
                sage: T.initial_states()
                [()]
                sage: T.final_states()
                [(), (1,)]

            Check some sequence::

                sage: sage.combinat.finite_state_machine.FSMOldProcessOutput = False
                sage: T([0, 1, 0, 1, 1, 0])
                [0, 0, 1, 0, 0, 1]

        #.  Counting the number of ``11`` blocks over the alphabet
            ``[0, 1]``::

                sage: T = transducers.CountSubblockOccurrences(
                ....:     [1, 1],
                ....:     [0, 1])
                sage: sorted(T.transitions())
                [Transition from () to (): 0|0,
                 Transition from () to (1,): 1|0,
                 Transition from (1,) to (): 0|0,
                 Transition from (1,) to (1,): 1|1]

            Check some sequence::

                sage: sage.combinat.finite_state_machine.FSMOldProcessOutput = False
                sage: T([0, 1, 0, 1, 1, 0])
                [0, 0, 0, 0, 1, 0]

        #.  Counting the number of ``1010`` blocks over the
            alphabet ``[0, 1, 2]``::

                sage: T = transducers.CountSubblockOccurrences(
                ....:     [1, 0, 1, 0],
                ....:     [0, 1, 2])
                sage: sorted(T.transitions())
                [Transition from () to (): 0|0,
                 Transition from () to (1,): 1|0,
                 Transition from () to (): 2|0,
                 Transition from (1,) to (1, 0): 0|0,
                 Transition from (1,) to (1,): 1|0,
                 Transition from (1,) to (): 2|0,
                 Transition from (1, 0) to (): 0|0,
                 Transition from (1, 0) to (1, 0, 1): 1|0,
                 Transition from (1, 0) to (): 2|0,
                 Transition from (1, 0, 1) to (1, 0): 0|1,
                 Transition from (1, 0, 1) to (1,): 1|0,
                 Transition from (1, 0, 1) to (): 2|0]
                sage: input =  [0, 1, 0, 1, 0, 1, 0, 0, 1, 0, 1, 0, 2]
                sage: output = [0, 0, 0, 0, 1, 0, 1, 0, 0, 0, 0, 1, 0]
                sage: sage.combinat.finite_state_machine.FSMOldProcessOutput = False
                sage: T(input) == output
                True

        """
        block_as_tuple = tuple(block)

        def starts_with(what, pattern):
            return len(what) >= len(pattern) \
                and what[:len(pattern)] == pattern

        def transition_function(read, input):
            current = read + (input, )
            if starts_with(block_as_tuple, current) \
                    and len(block_as_tuple) > len(current):
                return (current, 0)
            else:
                k = 1
                while not starts_with(block_as_tuple, current[k:]):
                    k += 1
                return (current[k:], int(block_as_tuple == current))

        T = Transducer(
            transition_function,
            input_alphabet=input_alphabet,
            output_alphabet=[0, 1],
            initial_states=[()])
        for s in T.iter_states():
            s.is_final = True
        return T


    def operator(self, operator, input_alphabet, number_of_operands=2):
        r"""
        Returns a transducer which realizes an operation
        on tuples over the given input alphabet.

        INPUT:

        - ``operator`` -- operator to realize. It is a function which
          takes ``number_of_operands`` input arguments (each out of
          ``input_alphabet``).

        - ``input_alphabet``  -- a list or other iterable.

        - ``number_of_operands`` -- (default: `2`) it specifies the number
          of input arguments the operator takes.

        OUTPUT:

        A transducer mapping an input letter `(i_1, \dots, i_n)` to
        `\mathrm{operator}(i_1, \dots, i_n)`. Here, `n` equals
        ``number_of_operands``.

        The input alphabet of the generated transducer is the cartesian
        product of ``number_of_operands`` copies of ``input_alphabet``.

        EXAMPLE:

        The following binary transducer realizes component-wise
        addition (this transducer is also available as :meth:`.add`)::

            sage: import operator
            sage: T = transducers.operator(operator.add, [0, 1])
            sage: T.transitions()
            [Transition from 0 to 0: (0, 0)|0,
             Transition from 0 to 0: (0, 1)|1,
             Transition from 0 to 0: (1, 0)|1,
             Transition from 0 to 0: (1, 1)|2]
            sage: T.input_alphabet
            [(0, 0), (0, 1), (1, 0), (1, 1)]
            sage: T.initial_states()
            [0]
            sage: T.final_states()
            [0]
            sage: T([(0, 0), (0, 1), (1, 0), (1, 1)])
            [0, 1, 1, 2]

        Note that for a unary operator the input letters of the
        new transducer are tuples of length `1`::

            sage: T = transducers.operator(abs,
            ....:                          [-1, 0, 1],
            ....:                          number_of_operands=1)
            sage: T([-1, 1, 0])
            Traceback (most recent call last):
            ...
            ValueError: Invalid input sequence.
            sage: T([(-1,), (1,), (0,)])
            [1, 1, 0]

        Compare this with the transducer generated by :meth:`.abs`::

            sage: T = transducers.abs([-1, 0, 1])
            sage: T([-1, 1, 0])
            [1, 1, 0]
        """
        from itertools import product

        def transition_function(state, operands):
            return (0, operator(*operands))
        pairs = list(product(input_alphabet, repeat=number_of_operands))
        return Transducer(transition_function,
                          input_alphabet=pairs,
                          initial_states=[0],
                          final_states=[0])


    def add(self, input_alphabet):
        """
        Returns a transducer which realizes addition on pairs over the
        given input alphabet.

        INPUT:

        - ``input_alphabet``  -- a list or other iterable.

        OUTPUT:

        A transducer mapping an input word `(i_0, i'_0)\ldots (i_k, i'_k)`
        to the word `(i_0 + i'_0)\ldots (i_k + i'_k)`.

        The input alphabet of the generated transducer is the cartesian
        product of two copies of ``input_alphabet``.

        EXAMPLE:

        The following transducer realizes letter-wise
        addition::

            sage: T = transducers.add([0, 1])
            sage: T.transitions()
            [Transition from 0 to 0: (0, 0)|0,
             Transition from 0 to 0: (0, 1)|1,
             Transition from 0 to 0: (1, 0)|1,
             Transition from 0 to 0: (1, 1)|2]
            sage: T.input_alphabet
            [(0, 0), (0, 1), (1, 0), (1, 1)]
            sage: T.initial_states()
            [0]
            sage: T.final_states()
            [0]
            sage: T([(0, 0), (0, 1), (1, 0), (1, 1)])
            [0, 1, 1, 2]
        """
        import operator
        return self.operator(operator.add, input_alphabet)


    def sub(self, input_alphabet):
        """
        Returns a transducer which realizes subtraction on pairs over
        the given input alphabet.

        INPUT:

        - ``input_alphabet``  -- a list or other iterable.

        OUTPUT:

        A transducer mapping an input word `(i_0, i'_0)\ldots (i_k, i'_k)`
        to the word `(i_0 - i'_0)\ldots (i_k - i'_k)`.

        The input alphabet of the generated transducer is the cartesian
        product of two copies of ``input_alphabet``.

        EXAMPLE:

        The following transducer realizes letter-wise
        subtraction::

            sage: T = transducers.sub([0, 1])
            sage: T.transitions()
            [Transition from 0 to 0: (0, 0)|0,
             Transition from 0 to 0: (0, 1)|-1,
             Transition from 0 to 0: (1, 0)|1,
             Transition from 0 to 0: (1, 1)|0]
            sage: T.input_alphabet
            [(0, 0), (0, 1), (1, 0), (1, 1)]
            sage: T.initial_states()
            [0]
            sage: T.final_states()
            [0]
            sage: T([(0, 0), (0, 1), (1, 0), (1, 1)])
            [0, -1, 1, 0]
        """
        import operator
        return self.operator(operator.sub, input_alphabet)

    def weight(self, input_alphabet, zero=0):
        r"""
        Returns a transducer which realizes the Hamming weight of the input
        over the given input alphabet.

        INPUT:

        - ``input_alphabet`` -- a list or other iterable.

        - ``zero`` -- the zero symbol in the alphabet used

        OUTPUT:

        A transducer mapping `i_0\ldots i_k` to `(i_0\neq 0)\ldots(i_k\neq 0)`.

        The Hamming weight is defined as the number of non-zero digits in the
        input sequence over the alphabet ``input_alphabet`` (see
        :wikipedia:`Hamming_weight`). The output sequence of the transducer is
        a unary encoding of the Hamming weight. Thus the sum of the output
        sequence is the Hamming weight of the input.

        EXAMPLES::

            sage: W = transducers.weight([-1, 0, 2])
            sage: W.transitions()
            [Transition from 0 to 0: -1|1,
             Transition from 0 to 0: 0|0,
             Transition from 0 to 0: 2|1]
            sage: unary_weight = W([-1, 0, 0, 2, -1])
            sage: unary_weight
            [1, 0, 0, 1, 1]
            sage: weight = add(unary_weight)
            sage: weight
            3

        Also the joint Hamming weight can be computed::

            sage: v1 = vector([-1, 0])
            sage: v0 = vector([0, 0])
            sage: W = transducers.weight([v1, v0])
            sage: unary_weight = W([v1, v0, v1, v0])
            sage: add(unary_weight)
            2

        For the input alphabet ``[-1, 0, 1]`` the weight transducer is the
        same as the absolute value transducer
        :meth:`~TransducerGenerators.abs`::

            sage: W = transducers.weight([-1, 0, 1])
            sage: A = transducers.abs([-1, 0, 1])
            sage: W == A
            True

        For other input alphabets, we can specify the zero symbol::

            sage: W = transducers.weight(['a', 'b'], zero='a')
            sage: add(W(['a', 'b', 'b']))
            2
        """
        def weight(state, input):
            weight = int(input != zero)
            return (0, weight)
        return Transducer(weight, input_alphabet=input_alphabet,
                          initial_states=[0],
                          final_states=[0])


    def abs(self, input_alphabet):
        """
        Returns a transducer which realizes the letter-wise
        absolute value of an input word over the given input alphabet.

        INPUT:

        - ``input_alphabet``  -- a list or other iterable.

        OUTPUT:

        A transducer mapping `i_0\ldots i_k`
        to `|i_0|\ldots |i_k|`.

        EXAMPLE:

        The following transducer realizes letter-wise
        absolute value::

            sage: T = transducers.abs([-1, 0, 1])
            sage: T.transitions()
            [Transition from 0 to 0: -1|1,
             Transition from 0 to 0: 0|0,
             Transition from 0 to 0: 1|1]
            sage: T.initial_states()
            [0]
            sage: T.final_states()
            [0]
            sage: T([-1, -1, 0, 1])
            [1, 1, 0, 1]

        """
        return Transducer(lambda state, input: (0, abs(input)),
                          input_alphabet=input_alphabet,
                          initial_states=[0],
                          final_states=[0])


    def GrayCode(self):
        """
        Returns a transducer converting the standard binary
        expansion to Gray code.

        INPUT:

        Nothing.

        OUTPUT:

        A transducer.

        Cf. the :wikipedia:`Gray_code` for a description of the Gray code.

        EXAMPLE::

            sage: G = transducers.GrayCode()
            sage: G
            Transducer with 3 states
            sage: sage.combinat.finite_state_machine.FSMOldProcessOutput = False
            sage: for v in srange(0, 10):
            ....:     print v, G(v.digits(base=2) + [0])
            0 []
            1 [1]
            2 [1, 1]
            3 [0, 1]
            4 [0, 1, 1]
            5 [1, 1, 1]
            6 [1, 0, 1]
            7 [0, 0, 1]
            8 [0, 0, 1, 1]
            9 [1, 0, 1, 1]

        In the example :ref:`Gray Code <finite_state_machine_gray_code_example>`
        in the documentation of the
        :mod:`~sage.combinat.finite_state_machine` module, the Gray code
        transducer is derived from the algorithm converting the binary
        expansion to the Gray code. The result is the same as the one
        given here.
        """
        z = ZZ(0)
        o = ZZ(1)
        return Transducer([[0, 1, z, None],
                           [0, 2, o, None],
                           [1, 1, z, z],
                           [1, 2, o, o],
                           [2, 1, z, o],
                           [2, 2, o, z]],
                          initial_states=[0],
                          final_states=[1])


# Easy access to the transducer generators from the command line:
transducers = TransducerGenerators()
