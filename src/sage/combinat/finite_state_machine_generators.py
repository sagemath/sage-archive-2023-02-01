r"""
Common Finite State Machines

Generators for common finite state machines.

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
    :meth:`~TransducerGenerators.GrayCode` | Returns a transducer realizing binary Gray code.

Functions and methods
---------------------

"""

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

        - ``input_alphabet`` -- the input alphabet.

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
        - ``input_alphabet`` -- input alphabet.

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

    def operator(self, operator, input_alphabet):
        r"""
        Returns a transducer which realizes the binary operator over
        an input alphabet.

        INPUT:

        - ``operator`` -- binary operator to realize (a map defined on
          ``input_alphabet`` `\times` ``input_alphabet``).
        - ``input_alphabet``  -- input alphabet.

        OUTPUT:

        A transducer mapping `(i_0, i'_0)\ldots (i_k, i'_k)`
        to `(i_0\odot i'_0)\ldots (i_k \odot i'_k)` where
        `\odot` stands for the operator given.

        EXAMPLE:

        The following transducer realizes component-wise
        addition::

            sage: import operator
            sage: T = transducers.operator(operator.add,
            ....:                           [0, 1])
            sage: T.transitions()
            [Transition from 0 to 0: (0, 0)|0,
             Transition from 0 to 0: (0, 1)|1,
             Transition from 0 to 0: (1, 0)|1,
             Transition from 0 to 0: (1, 1)|2]
            sage: T.input_alphabet
            [(0, 0), (0, 1), (1, 0), (1, 1)]
        """
        from itertools import product

        def transition_function(state, (i, o)):
            return(0, operator(i, o))
        pairs = list(product(input_alphabet, repeat=2))
        return Transducer(transition_function,
                          input_alphabet=pairs,
                          initial_states=[0])

    def add(self, input_alphabet):
        """
        Returns a transducer which realizes the component-wise
        addition over an input alphabet.

        INPUT:

        - ``input_alphabet``  -- input alphabet.

        OUTPUT:

        A transducer mapping `(i_0, i'_0)\ldots (i_k, i'_k)`
        to `(i_0 + i'_0)\ldots (i_k + i'_k)`.

        EXAMPLE:

        The following transducer realizes component-wise
        addition::

            sage: T = transducers.add([0, 1])
            sage: T.transitions()
            [Transition from 0 to 0: (0, 0)|0,
             Transition from 0 to 0: (0, 1)|1,
             Transition from 0 to 0: (1, 0)|1,
             Transition from 0 to 0: (1, 1)|2]
            sage: T.input_alphabet
            [(0, 0), (0, 1), (1, 0), (1, 1)]
        """
        import operator
        return self.operator(operator.add, input_alphabet)

    def sub(self, input_alphabet):
        """
        Returns a transducer which realizes the component-wise
        addition over an input alphabet.

        INPUT:

        - ``input_alphabet``  -- input alphabet.

        OUTPUT:

        A transducer mapping `(i_0, i'_0)\ldots (i_k, i'_k)`
        to `(i_0 - i'_0)\ldots (i_k - i'_k)`.

        EXAMPLE:

        The following transducer realizes component-wise
        subtraction::

            sage: T = transducers.sub([0, 1])
            sage: T.transitions()
            [Transition from 0 to 0: (0, 0)|0,
             Transition from 0 to 0: (0, 1)|-1,
             Transition from 0 to 0: (1, 0)|1,
             Transition from 0 to 0: (1, 1)|0]
            sage: T.input_alphabet
            [(0, 0), (0, 1), (1, 0), (1, 1)]

        """
        import operator
        return self.operator(operator.sub, input_alphabet)


    def abs(self, input_alphabet):
        """
        Returns a transducer which realizes the component-wise
        absolute value over an input alphabet.

        INPUT:

        - ``input_alphabet``  -- input alphabet.

        OUTPUT:

        A transducer mapping `i_0\ldots i_k`
        to `|i_0|\ldots |i_k|`.

        EXAMPLE:

        The following transducer realizes component-wise
        absolute value::

            sage: T = transducers.abs([-1, 0, 1])
            sage: T.transitions()
            [Transition from 0 to 0: -1|1,
             Transition from 0 to 0: 0|0,
             Transition from 0 to 0: 1|1]
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
            sage: for v in srange(0,10):
            ....:     print v, G(v.digits(base=2)+[0])
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
