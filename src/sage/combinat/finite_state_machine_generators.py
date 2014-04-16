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
    :meth:`~TransducerGenerators.CountSubblockOccurrences` | Returns a transducer counting the occurrences of a subblock.

AUTHORS:

- Clemens Heuberger (2014-04-07): initial version
- Sara Kropf (2014-04-10): some changes in TransducerGenerator
- Daniel Krenn (2014-04-15): improved common docstring during review

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

class TransducerGenerators(object):
    r"""
    A class consisting of constructors for several common transducers.

    A list of all transducers in this database is available via tab
    completion. Type "``transducers.``" and then hit tab to see which
    transducers are available.

    The transducers currently in this class include:

    - :meth:`~Identity`
    - :meth:`~CountSubblockOccurrences`

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

# Easy access to the transducer generators from the command line:
transducers = TransducerGenerators()
