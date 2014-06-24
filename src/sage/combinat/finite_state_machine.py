# -*- coding: utf-8 -*-
r"""
Finite State Machines, Automata, Transducers

This module adds support for finite state machines, automata and
transducers. See class :class:`FiniteStateMachine` and the examples
below for details creating one.

Examples
========


A simple finite state machine
-----------------------------

We can easily create a finite state machine by

::

    sage: fsm = FiniteStateMachine()
    sage: fsm
    Finite state machine with 0 states

By default this is the empty finite state machine, so not very
interesting. Let's create some states and transitions::

    sage: from sage.combinat.finite_state_machine import FSMState, FSMTransition
    sage: day = FSMState('day')
    sage: night = FSMState('night')
    sage: sunrise = FSMTransition(night, day)
    sage: sunset = FSMTransition(day, night)

And now let's add those states and transitions to our finite state machine::

    sage: fsm.add_transition(sunrise)
    Transition from 'night' to 'day': -|-
    sage: fsm.add_transition(sunset)
    Transition from 'day' to 'night': -|-

Note that the states are added automatically, since they are present
in the transitions. We could add the states manually by

::

    sage: fsm.add_state(day)
    'day'
    sage: fsm.add_state(night)
    'night'

Anyhow, we got the following finite state machine::

    sage: fsm
    Finite state machine with 2 states

We can also obtain the underlying directed graph by

::

    sage: fsm.graph()
    Digraph on 2 vertices

To visualize a finite state machine, we can use
:func:`~sage.misc.latex.latex` and run the result through LaTeX,
see the section on :ref:`finite_state_machine_LaTeX_output`
below.

Alternatively, we could have created the finite state machine above
simply by

::

    sage: FiniteStateMachine([('night', 'day'), ('day', 'night')])
    Finite state machine with 2 states

or by

::

    sage: fsm = FiniteStateMachine()
    sage: day = fsm.add_state('day')
    sage: night = fsm.add_state('night')
    sage: sunrise = fsm.add_transition(night, day)
    sage: sunset = fsm.add_transition(day, night)
    sage: fsm
    Finite state machine with 2 states

.. _finite_state_machine_recognizing_NAFs_example:

A simple Automaton (recognizing NAFs)
---------------------------------------

We want to build an automaton which recognizes non-adjacent forms
(NAFs), i.e., sequences which have no adjacent non-zeros.
We use `0`, `1`, and `-1` as digits::

    sage: NAF = Automaton(
    ....:     {'A': [('A', 0), ('B', 1), ('B', -1)], 'B': [('A', 0)]})
    sage: NAF.state('A').is_initial = True
    sage: NAF.state('A').is_final = True
    sage: NAF.state('B').is_final = True
    sage: NAF
    Automaton with 2 states

Of course, we could have specified the initial and final states
directly in the definition of ``NAF`` by ``initial_states=['A']`` and
``final_states=['A', 'B']``.

So let's test the automaton with some input::

    sage: sage.combinat.finite_state_machine.FSMOldProcessOutput = False  # activate new output behavior
    sage: NAF([0])
    True
    sage: NAF([0, 1])
    True
    sage: NAF([1, -1])
    False
    sage: NAF([0, -1, 0, 1])
    True
    sage: NAF([0, -1, -1, -1, 0])
    False
    sage: NAF([-1, 0, 0, 1, 1])
    False

Alternatively, we could call that by

::

    sage: NAF.process([0, -1, 0, 1])
    (True, 'B')

which gives additionally the state in which we arrived.


.. _finite_state_machine_LaTeX_output:

LaTeX output
------------

We can visualize a finite state machine by converting it to LaTeX by
using the usual function :func:`~sage.misc.latex.latex`. Within LaTeX,
TikZ is used for typesetting the graphics, see the
:wikipedia:`PGF/TikZ`.

::

    sage: print latex(NAF)
    \begin{tikzpicture}[auto, initial text=, >=latex]
    \node[state, accepting, initial] (v0) at (3.000000, 0.000000) {$\text{\texttt{A}}$};
    \node[state, accepting] (v1) at (-3.000000, 0.000000) {$\text{\texttt{B}}$};
    \path[->] (v1.5.00) edge node[rotate=0.00, anchor=south] {$0$} (v0.175.00);
    \path[->] (v0.185.00) edge node[rotate=360.00, anchor=north] {$1, -1$} (v1.355.00);
    \path[->] (v0) edge[loop above] node {$0$} ();
    \end{tikzpicture}

We can turn this into a graphical representation. Before doing this,
we have to :func:`setup the latex preamble <setup_latex_preamble>` and
make sure that TikZ pictures are not rendered by mathjax, but by
actually running LaTeX.

::

    sage: sage.combinat.finite_state_machine.setup_latex_preamble()
    sage: latex.mathjax_avoid_list('tikzpicture')
    sage: view(NAF) # not tested

To actually see this, use the live documentation in the Sage notebook
and execute the cells in this and the previous section.

Several options can be set to customize the output, see
:meth:`~FiniteStateMachine.latex_options` for details. In particular,
we use :meth:`~FiniteStateMachine.format_letter_negative` to format
`-1` as `\overline{1}`.

::

    sage: NAF.latex_options(
    ....:     coordinates={'A': (0, 0),
    ....:                  'B': (6, 0)},
    ....:     initial_where={'A': 'below'},
    ....:     format_letter=NAF.format_letter_negative,
    ....:     format_state_label=lambda x:
    ....:         r'\mathcal{%s}' % x.label()
    ....: )
    sage: print latex(NAF)
    \begin{tikzpicture}[auto, initial text=, >=latex]
    \node[state, accepting, initial, initial where=below] (v0) at (0.000000, 0.000000) {$\mathcal{A}$};
    \node[state, accepting] (v1) at (6.000000, 0.000000) {$\mathcal{B}$};
    \path[->] (v1.185.00) edge node[rotate=360.00, anchor=north] {$0$} (v0.355.00);
    \path[->] (v0.5.00) edge node[rotate=0.00, anchor=south] {$1, \overline{1}$} (v1.175.00);
    \path[->] (v0) edge[loop above] node {$0$} ();
    \end{tikzpicture}
    sage: view(NAF) # not tested

A simple transducer (binary inverter)
-------------------------------------

Let's build a simple transducer, which rewrites a binary word by
iverting each bit::

    sage: inverter = Transducer({'A': [('A', 0, 1), ('A', 1, 0)]},
    ....:     initial_states=['A'], final_states=['A'])

We can look at the states and transitions::

    sage: inverter.states()
    ['A']
    sage: for t in inverter.transitions():
    ....:     print t
    Transition from 'A' to 'A': 0|1
    Transition from 'A' to 'A': 1|0

Now we apply a word to it and see what the transducer does::

    sage: inverter([0, 1, 0, 0, 1, 1, 0, 0, 0, 1, 1, 1])
    [1, 0, 1, 1, 0, 0, 1, 1, 1, 0, 0, 0]

``True`` means, that we landed in a final state, that state is labeled
``'A'``, and we also got an output.


A transducer which performs division by `3` in binary
-----------------------------------------------------

Now we build a transducer, which divides a binary number by `3`.
The labels of the states are the remainder of the division.
The transition function is

::

    sage: def f(state_from, read):
    ....:     if state_from + read <= 1:
    ....:         state_to = 2*state_from + read
    ....:         write = 0
    ....:     else:
    ....:         state_to = 2*state_from + read - 3
    ....:         write = 1
    ....:     return (state_to, write)

which assumes reading a binary number from left to right.
We get the transducer with

::

    sage: D = Transducer(f, initial_states=[0], final_states=[0],
    ....:                input_alphabet=[0, 1])

Let us try to divide `12` by `3`::

    sage: D([1, 1, 0, 0])
    [0, 1, 0, 0]

Now we want to divide `13` by `3`::

    sage: D([1, 1, 0, 1])
    Traceback (most recent call last):
    ...
    ValueError: Invalid input sequence.

The raised ``ValueError``
means `13` is not divisible by `3`.

.. _finite_state_machine_gray_code_example:

Gray Code
---------

The Gray code is a binary :wikipedia:`numeral system <Numeral_system>`
where two successive values differ in only one bit, cf. the
:wikipedia:`Gray_code`. The Gray code of an integer `n` is obtained by
a bitwise xor between the binary expansion of `n` and the binary
expansion of `\\lfloor n/2\\rfloor`; the latter corresponds to a
shift by one position in binary.

The purpose of this example is to construct a transducer converting the
standard binary expansion to the Gray code by translating this
construction into operations with transducers.

For this construction, the least significant digit is at
the left-most position.
Note that it is easier to shift everything to
the right first, i.e., multiply by `2` instead of building
`\\lfloor n/2\\rfloor`. Then, we take the input xor with the right
shift of the input and forget the first letter.

We first construct a transducer shifting the binary expansion to the
right. This requires storing the previously read digit in a state.

::

    sage: def shift_right_transition(state, digit):
    ....:     if state == 'I':
    ....:         return (digit, None)
    ....:     else:
    ....:         return (digit, state)
    sage: shift_right_transducer = Transducer(
    ....:     shift_right_transition,
    ....:     initial_states=['I'],
    ....:     input_alphabet=[0, 1],
    ....:     final_states=[0])
    sage: shift_right_transducer.transitions()
    [Transition from 'I' to 0: 0|-,
     Transition from 'I' to 1: 1|-,
     Transition from 0 to 0: 0|0,
     Transition from 0 to 1: 1|0,
     Transition from 1 to 0: 0|1,
     Transition from 1 to 1: 1|1]
    sage: sage.combinat.finite_state_machine.FSMOldProcessOutput = False
    sage: shift_right_transducer([0, 1, 1, 0])
    [0, 1, 1]
    sage: shift_right_transducer([1, 0, 0])
    [1, 0]

The output of the shifts above look a bit weird (from a right-shift
transducer, we would expect, for example, that ``[1, 0, 0]`` was
mapped to ``[0, 1, 0]``), since we write ``None`` instead of the zero
at the left. Further, note that only `0` is listed as a final state
as we have to enforce that a most significant zero is read as the last
input letter in order to flush the last digit::

    sage: shift_right_transducer([0, 1, 0, 1])
    Traceback (most recent call last):
    ...
    ValueError: Invalid input sequence.

Next, we construct the transducer performing the xor operation. We also
have to take ``None`` into account as our ``shift_right_transducer``
waits one iteration until it starts writing output. This corresponds
with our intention to forget the first letter.

::

    sage: def xor_transition(state, digits):
    ....:    if digits[0] is None or digits[1] is None:
    ....:        return (0, None)
    ....:    else:
    ....:        return (0, digits[0].__xor__(digits[1]))
    sage: from itertools import product
    sage: xor_transducer = Transducer(
    ....:    xor_transition,
    ....:    initial_states=[0],
    ....:    final_states=[0],
    ....:    input_alphabet=list(product([None, 0, 1], [0, 1])))
    sage: xor_transducer.transitions()
    [Transition from 0 to 0: (None, 0)|-,
     Transition from 0 to 0: (None, 1)|-,
     Transition from 0 to 0: (0, 0)|0,
     Transition from 0 to 0: (0, 1)|1,
     Transition from 0 to 0: (1, 0)|1,
     Transition from 0 to 0: (1, 1)|0]
    sage: xor_transducer([(None, 0), (None, 1), (0, 0), (0, 1), (1, 0), (1, 1)])
    [0, 1, 1, 0]
    sage: xor_transducer([(0, None)])
    Traceback (most recent call last):
    ...
    ValueError: Invalid input sequence.

The transducer computing the Gray code is then constructed as a
:meth:`cartesian product <Transducer.cartesian_product>` between the
shifted version and the original input (represented here by the
``shift_right_transducer`` and the :meth:`identity transducer
<sage.combinat.finite_state_machine_generators.TransducerGenerators.Identity>`,
respectively). This cartesian product is then fed into the
``xor_transducer`` as a :meth:`composition
<FiniteStateMachine.composition>` of transducers.

As described in :meth:`Transducer.cartesian_product`, we have to
temporarily set
``finite_state_machine.FSMOldCodeTransducerCartesianProduct`` to
``False`` in order to disable backwards compatible code.

::

    sage: sage.combinat.finite_state_machine.FSMOldCodeTransducerCartesianProduct = False
    sage: product_transducer = shift_right_transducer.cartesian_product(transducers.Identity([0, 1]))
    sage: sage.combinat.finite_state_machine.FSMOldCodeTransducerCartesianProduct = True
    sage: Gray_transducer = xor_transducer(product_transducer)
    sage: Gray_transducer.transitions()
    [Transition from (('I', 0), 0) to ((0, 0), 0): 0|-,
     Transition from (('I', 0), 0) to ((1, 0), 0): 1|-,
     Transition from ((0, 0), 0) to ((0, 0), 0): 0|0,
     Transition from ((0, 0), 0) to ((1, 0), 0): 1|1,
     Transition from ((1, 0), 0) to ((0, 0), 0): 0|1,
     Transition from ((1, 0), 0) to ((1, 0), 0): 1|0]

There is a :meth:`prepackaged transducer
<sage.combinat.finite_state_machine_generators.TransducerGenerators.GrayCode>`
for Gray code, let's see whether they agree. We have to use
:meth:`~FiniteStateMachine.relabeled` to relabel our states with
integers. Note that equality of finite state machines does not compare
initial and final states, so we have to do it on our own here.

::

    sage: constructed = Gray_transducer.relabeled()
    sage: packaged = transducers.GrayCode()
    sage: constructed == packaged
    True
    sage: constructed.initial_states() == packaged.initial_states()
    True
    sage: constructed.final_states() == packaged.final_states()
    True

Finally, we check that this indeed computes the Gray code of the first
10 non-negative integers. Note that we add a trailing zero at the most
significant position of the input in order to flush all output digits.
This is due to the left shift which delays its output.

::

    sage: for n in srange(10):
    ....:     Gray_transducer(n.bits() + [0])
    []
    [1]
    [1, 1]
    [0, 1]
    [0, 1, 1]
    [1, 1, 1]
    [1, 0, 1]
    [0, 0, 1]
    [0, 0, 1, 1]
    [1, 0, 1, 1]


Using the hook-functions
------------------------

Let's use the previous example "divison by `3`" to demonstrate the
optional state and transition parameters ``hook``.

First, we define, what those functions should do. In our case, this is
just saying in which state we are and which transition we take

::

    sage: def state_hook(state, process):
    ....:     print "We are now in State %s." % (state.label(),)
    sage: from sage.combinat.finite_state_machine import FSMWordSymbol
    sage: def transition_hook(transition, process):
    ....:     print ("Currently we go from %s to %s, "
    ....:            "reading %s and writing %s." % (
    ....:                transition.from_state, transition.to_state,
    ....:                FSMWordSymbol(transition.word_in),
    ....:                FSMWordSymbol(transition.word_out)))

Now, let's add these hook-functions to the existing transducer::

    sage: for s in D.iter_states():
    ....:     s.hook = state_hook
    sage: for t in D.iter_transitions():
    ....:     t.hook = transition_hook

Rerunning the process again now gives the following output::

    sage: D.process([1, 1, 0, 1])
    We are now in State 0.
    Currently we go from 0 to 1, reading 1 and writing 0.
    We are now in State 1.
    Currently we go from 1 to 0, reading 1 and writing 1.
    We are now in State 0.
    Currently we go from 0 to 0, reading 0 and writing 0.
    We are now in State 0.
    Currently we go from 0 to 1, reading 1 and writing 0.
    We are now in State 1.
    (False, 1, [0, 1, 0, 0])

The example above just explains the basic idea of using
hook-functions. In the following, we will use those hooks more seriously.


Detecting sequences with same number of `0` and `1`
---------------------------------------------------

Suppose we have a binary input and want to accept all sequences with
the same number of `0` and `1`. This cannot be done with a finite
automaton. Anyhow, we can make usage of the hook functions to extend
our finite automaton by a counter::

    sage: from sage.combinat.finite_state_machine import FSMState, FSMTransition
    sage: C = FiniteStateMachine()
    sage: def update_counter(state, process):
    ....:     l = process.read_letter()
    ....:     process.fsm.counter += 1 if l == 1 else -1
    ....:     if process.fsm.counter > 0:
    ....:         next_state = 'positive'
    ....:     elif process.fsm.counter < 0:
    ....:         next_state = 'negative'
    ....:     else:
    ....:         next_state = 'zero'
    ....:     return FSMTransition(state, process.fsm.state(next_state),
    ....:                          l, process.fsm.counter)
    sage: C.add_state(FSMState('zero', hook=update_counter,
    ....:             is_initial=True, is_final=True))
    'zero'
    sage: C.add_state(FSMState('positive', hook=update_counter))
    'positive'
    sage: C.add_state(FSMState('negative', hook=update_counter))
    'negative'

Now, let's input some sequence::

    sage: C.counter = 0; C([1, 1, 1, 1, 0, 0])
    (False, 'positive', [1, 2, 3, 4, 3, 2])

The result is False, since there are four `1` but only two `0`. We
land in the state ``positive`` and we can also see the values of the
counter in each step.

Let's try some other examples::

    sage: C.counter = 0; C([1, 1, 0, 0])
    (True, 'zero', [1, 2, 1, 0])
    sage: C.counter = 0; C([0, 1, 0, 0])
    (False, 'negative', [-1, 0, -1, -2])


AUTHORS:

- Daniel Krenn (2012-03-27): initial version
- Clemens Heuberger (2012-04-05): initial version
- Sara Kropf (2012-04-17): initial version
- Clemens Heuberger (2013-08-21): release candidate for Sage patch
- Daniel Krenn (2013-08-21): release candidate for Sage patch
- Sara Kropf (2013-08-21): release candidate for Sage patch
- Clemens Heuberger (2013-09-02): documentation improved
- Daniel Krenn (2013-09-13): comments from trac worked in
- Clemens Heuberger (2013-11-03): output (labels) of determinisation,
    product, composition, etc. changed (for consistency),
    representation of state changed, documentation improved
- Daniel Krenn (2013-11-04): whitespaces in documentation corrected
- Clemens Heuberger (2013-11-04): full_group_by added
- Daniel Krenn (2013-11-04): next release candidate for Sage patch
- Sara Kropf (2013-11-08): fix for adjacency matrix
- Clemens Heuberger (2013-11-11): fix for prepone_output
- Daniel Krenn (2013-11-11): comments from trac 15078 included:
    docstring of FiniteStateMachine rewritten, Automaton and Transducer
    inherited from FiniteStateMachine
- Daniel Krenn (2013-11-25): documentation improved according to
    comments from trac 15078

ACKNOWLEDGEMENT:

- Daniel Krenn, Clemens Heuberger and Sara Kropf are supported by the
  Austrian Science Fund (FWF): P 24644-N26.

"""

#*****************************************************************************
#  Copyright (C) 2012, 2013 Daniel Krenn <math+sage@danielkrenn.at>
#                2012, 2013 Clemens Heuberger <clemens.heuberger@aau.at>
#                2012, 2013 Sara Kropf <sara.kropf@aau.at>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#  as published by the Free Software Foundation; either version 2 of
#  the License, or (at your option) any later version.
#                http://www.gnu.org/licenses/
#*****************************************************************************

from sage.structure.sage_object import SageObject
from sage.graphs.digraph import DiGraph
from sage.matrix.constructor import matrix
from sage.rings.integer_ring import ZZ
from sage.rings.real_mpfr import RR
from sage.symbolic.ring import SR
from sage.calculus.var import var
from sage.misc.latex import latex
from sage.misc.misc import verbose
from sage.functions.trig import cos, sin, atan2
from sage.symbolic.constants import pi

from copy import copy
from copy import deepcopy

import itertools
from itertools import imap, ifilter
from collections import defaultdict


def full_group_by(l, key=lambda x: x):
    """
    Group iterable ``l`` by values of ``key``.

    INPUT:

    - iterable ``l``
    - key function ``key``

    OUTPUT:

    A list of pairs ``(k, elements)`` such that ``key(e)=k`` for all
    ``e`` in ``elements``.

    This is similar to ``itertools.groupby`` except that lists are
    returned instead of iterables and no prior sorting is required.

    We do not require

    - that the keys are sortable (in contrast to the
      approach via ``sorted`` and ``itertools.groupby``) and
    - that the keys are hashable (in contrast to the
      implementation proposed in `<http://stackoverflow.com/a/15250161>`_).

    However, it is required

    - that distinct keys have distinct ``str``-representations.

    The implementation is inspired by
    `<http://stackoverflow.com/a/15250161>`_, but non-hashable keys are
    allowed.

    EXAMPLES::

        sage: from sage.combinat.finite_state_machine import full_group_by
        sage: t = [2/x, 1/x, 2/x]
        sage: r = full_group_by([0,1,2], key=lambda i:t[i])
        sage: sorted(r, key=lambda p:p[1])
        [(2/x, [0, 2]), (1/x, [1])]
        sage: from itertools import groupby
        sage: for k, elements in groupby(sorted([0,1,2],
        ....:                            key=lambda i:t[i]),
        ....:                            key=lambda i:t[i]):
        ....:     print k, list(elements)
        2/x [0]
        1/x [1]
        2/x [2]

    Note that the behavior is different from ``itertools.groupby``
    because neither `1/x<2/x` nor `2/x<1/x` does hold.

    Here, the result ``r`` has been sorted in order to guarantee a
    consistent order for the doctest suite.
    """
    elements = defaultdict(list)
    original_keys = {}
    for item in l:
        k = key(item)
        s = str(k)
        if s in original_keys:
            if original_keys[s]!=k:
                raise ValueError("Two distinct elements with representation "
                                 "%s " % s)
        else:
            original_keys[s]=k
        elements[s].append(item)
    return [(original_keys[s], values ) for (s, values) in elements.items()]

#*****************************************************************************

FSMEmptyWordSymbol = '-'
EmptyWordLaTeX = r'\varepsilon'
EndOfWordLaTeX = r'\$'
FSMOldCodeTransducerCartesianProduct = True
FSMOldProcessOutput = True  # See trac #16132 (deprecation).
tikz_automata_where = {"right": 0,
                       "above": 90,
                       "left": 180,
                       "below": 270}

def FSMLetterSymbol(letter):
    """
    Returns a string associated to the input letter.

    INPUT:

    - ``letter`` -- the input letter or ``None`` (representing the
      empty word).

    OUTPUT:

    If ``letter`` is ``None`` the symbol for the empty word
    ``FSMEmptyWordSymbol`` is returned, otherwise the string
    associated to the letter.

    EXAMPLES::

        sage: from sage.combinat.finite_state_machine import FSMLetterSymbol
        sage: FSMLetterSymbol(0)
        '0'
        sage: FSMLetterSymbol(None)
        '-'
    """
    return FSMEmptyWordSymbol if letter is None else repr(letter)


def FSMWordSymbol(word):
    """
    Returns a string of ``word``. It may returns the symbol of the
    empty word ``FSMEmptyWordSymbol``.

    INPUT:

    - ``word`` -- the input word.

    OUTPUT:

    A string of ``word``.

    EXAMPLES::

        sage: from sage.combinat.finite_state_machine import FSMWordSymbol
        sage: FSMWordSymbol([0, 1, 1])
        '0,1,1'
    """
    if not isinstance(word, list):
        return FSMLetterSymbol(word)
    if len(word) == 0:
        return FSMEmptyWordSymbol
    s = ''
    for letter in word:
        s += (',' if len(s) > 0 else '') + FSMLetterSymbol(letter)
    return s


#*****************************************************************************


def is_FSMState(S):
    """
    Tests whether or not ``S`` inherits from :class:`FSMState`.

    TESTS::

        sage: from sage.combinat.finite_state_machine import is_FSMState, FSMState
        sage: is_FSMState(FSMState('A'))
        True
    """
    return isinstance(S, FSMState)


class FSMState(SageObject):
    """
    Class for a state of a finite state machine.

    INPUT:

    - ``label`` -- the label of the state.

    - ``word_out`` -- (default: ``None``) a word that is written when
      the state is reached.

    - ``is_initial`` -- (default: ``False``)

    - ``is_final`` -- (default: ``False``)

    - ``final_word_out`` -- (default: ``None``) a word that is written when
      the state is reached as the last state of some input; only for final
      states

      .. WARNING::

          ``final_word_out`` is not implemented everywhere. Currently it is
          implemented in :meth:`FiniteStateMachine.process`,
          in the LaTeX output,
          :meth:`FiniteStateMachine.composition`,
          :meth:`FiniteStateMachine.output_projection`,
          :meth:`FiniteStateMachine.prepone_output`,
          :meth:`Transducer.cartesian_product`, but not in
          :meth:`FiniteStateMachine.determine_alphabets`,
          :meth:`FiniteStateMachine.equivalence_classes`,
          :meth:`FiniteStateMachine.transposition`,
          :meth:`Transducer.intersection` and
          :meth:`Transducer.simplification`.

    - ``hook`` -- (default: ``None``) A function which is called when
      the state is reached during processing input.

    - ``color`` -- (default: ``None``) In order to distinguish states,
      they can be given an arbitrary "color" (an arbitrary object).
      This is used in :meth:`FiniteStateMachine.equivalence_classes`:
      states of different colors are never considered to be
      equivalent. Note that :meth:`Automaton.determinisation` requires
      that ``color`` is hashable.

    - ``allow_label_None`` -- (default: ``False``) If ``True`` allows also
      ``None`` as label. Note that a state with label ``None`` is used in
      :class:`FSMProcessIterator`.

    OUTPUT:

    Returns a state of a finite state machine.

    EXAMPLES::

        sage: from sage.combinat.finite_state_machine import FSMState
        sage: A = FSMState('state 1', word_out=0, is_initial=True)
        sage: A
        'state 1'
        sage: A.label()
        'state 1'
        sage: B = FSMState('state 2')
        sage: A == B
        False

    We can also define a final output word of a final state which is
    used if the input of a transducer leads to this state. Such final
    output words are used in subsequential transducers. ::

        sage: C = FSMState('state 3', is_final=True, final_word_out='end')
        sage: C.final_word_out
        ['end']

    The final output word can be a single letter, ``None`` or a list of
    letters::

        sage: A = FSMState('A')
        sage: A.is_final = True
        sage: A.final_word_out = 2
        sage: A.final_word_out
        [2]
        sage: A.final_word_out = [2, 3]
        sage: A.final_word_out
        [2, 3]

    Only final states can have a final output word which is not
    ``None``::

        sage: B = FSMState('B')
        sage: B.final_word_out is None
        True
        sage: B.final_word_out = 2
        Traceback (most recent call last):
        ...
        ValueError: Only final states can have a final output word,
        but state B is not final.

    Setting the ``final_word_out`` of a final state to ``None`` is the same as
    setting it to ``[]`` and is also the default for a final state::

        sage: C = FSMState('C', is_final=True)
        sage: C.final_word_out
        []
        sage: C.final_word_out = None
        sage: C.final_word_out
        []
        sage: C.final_word_out = []
        sage: C.final_word_out
        []

    It is not allowed to use ``None`` as a label::

        sage: from sage.combinat.finite_state_machine import FSMState
        sage: FSMState(None)
        Traceback (most recent call last):
        ...
        ValueError: Label None reserved for a special state, choose another label.

    This can be overridden by::

        sage: FSMState(None, allow_label_None=True)
        None

    Note that :meth:`Automaton.determinisation` requires that ``color``
    is hashable::

        sage: A = Automaton([[0, 0, 0]], initial_states=[0])
        sage: A.state(0).color = []
        sage: A.determinisation()
        Traceback (most recent call last):
        ...
        TypeError: unhashable type: 'list'
        sage: A.state(0).color = ()
        sage: A.determinisation()
        Automaton with 1 states
    """
    def __init__(self, label, word_out=None,
                 is_initial=False, is_final=False, final_word_out=None,
                 hook=None, color=None, allow_label_None=False):
        """
        See :class:`FSMState` for more information.

        EXAMPLES::

            sage: from sage.combinat.finite_state_machine import FSMState
            sage: FSMState('final', is_final=True)
            'final'

        TESTS::

            sage: A = FSMState('A', is_final=True)
            sage: A.final_word_out
            []
            sage: A.is_final = True
            sage: A = FSMState('A', is_final=True, final_word_out='end')
            sage: A.final_word_out
            ['end']
            sage: A = FSMState('A', is_final=True,
            ....:              final_word_out=['e', 'n', 'd'])
            sage: A.final_word_out
            ['e', 'n', 'd']
            sage: A = FSMState('A', is_final=True, final_word_out=[])
            sage: A.final_word_out
            []
            sage: A = FSMState('A', is_final=True, final_word_out=None)
            sage: A.final_word_out
            []
            sage: A = FSMState('A', is_final=False)
            sage: A.final_word_out is None
            True
            sage: A.is_final = False
            sage: A = FSMState('A', is_final=False, final_word_out='end')
            Traceback (most recent call last):
            ...
            ValueError: Only final states can have a final output word,
            but state A is not final.
            sage: A = FSMState('A', is_final=False,
            ....:              final_word_out=['e', 'n', 'd'])
            Traceback (most recent call last):
            ...
            ValueError: Only final states can have a final output word,
            but state A is not final.
            sage: A = FSMState('A', is_final=False, final_word_out=None)
            sage: A.final_word_out is None
            True
            sage: A = FSMState('A', is_final=False, final_word_out=[])
            Traceback (most recent call last):
            ...
            ValueError: Only final states can have a final output word,
            but state A is not final.
        """
        if not allow_label_None and label is None:
            raise ValueError("Label None reserved for a special state, "
                             "choose another label.")
        self._label_ = label

        if isinstance(word_out, list):
            self.word_out = word_out
        elif word_out is not None:
            self.word_out = [word_out]
        else:
            self.word_out = []

        self.is_initial = is_initial
        self._final_word_out_ = None
        self.is_final = is_final
        self.final_word_out = final_word_out

        if hook is not None:
            if hasattr(hook, '__call__'):
                self.hook = hook
            else:
                raise TypeError('Wrong argument for hook.')

        self.color = color


    def __lt__(self, other):
        """
        Returns True if label of ``self`` is less than label of
        ``other``.

        INPUT:

        - `other` -- a state.

        OUTPUT:

        True or False.

        EXAMPLE::

            sage: from sage.combinat.finite_state_machine import FSMState
            sage: FSMState(0) < FSMState(1)
            True
        """
        return self.label() < other.label()


    @property
    def final_word_out(self):
        """
        The final output word of a final state which is written if the
        state is reached as the last state of the input of the finite
        state machine. For a non-final state, the value is ``None``.

        ``final_word_out`` can be a single letter, a list or ``None``,
        but for a final-state, it is always saved as a list.

        EXAMPLES::

            sage: from sage.combinat.finite_state_machine import FSMState
            sage: A = FSMState('A', is_final=True, final_word_out=2)
            sage: A.final_word_out
            [2]
            sage: A.final_word_out = 3
            sage: A.final_word_out
            [3]
            sage: A.final_word_out = [3, 4]
            sage: A.final_word_out
            [3, 4]
            sage: A.final_word_out = None
            sage: A.final_word_out
            []
            sage: B = FSMState('B')
            sage: B.final_word_out is None
            True

        A non-final state cannot have a final output word::

            sage: B.final_word_out = [3, 4]
            Traceback (most recent call last):
            ...
            ValueError: Only final states can have a final
            output word, but state B is not final.
        """
        return self._final_word_out_


    @final_word_out.setter
    def final_word_out(self, final_word_out):
        """
        Sets the value of the final output word of a final state.

        INPUT:

        - ``final_word_out`` -- a list, any element or ``None``.

        OUTPUT:

        Nothing.

        TESTS::

            sage: from sage.combinat.finite_state_machine import FSMState
            sage: B = FSMState('B')
            sage: B.final_word_out = []
            Traceback (most recent call last):
            ...
            ValueError: Only final states can have a final
            output word, but state B is not final.
            sage: B.final_word_out = None
            sage: B.final_word_out is None
            True
        """
        if not self.is_final:
            if final_word_out is not None:
                raise ValueError("Only final states can have a "
                                 "final output word, but state %s is not final."
                                 % (self.label()))
            else:
                self._final_word_out_ = None
        elif isinstance(final_word_out, list):
            self._final_word_out_ = final_word_out
        elif final_word_out is not None:
            self._final_word_out_ = [final_word_out]
        else:
            self._final_word_out_ = []


    @property
    def is_final(self):
        """
        Describes whether the state is final or not.

        ``True`` if the state is final and ``False`` otherwise.

        EXAMPLES::

            sage: from sage.combinat.finite_state_machine import FSMState
            sage: A = FSMState('A', is_final=True, final_word_out=3)
            sage: A.is_final
            True
            sage: A.is_final = False
            Traceback (most recent call last):
            ...
            ValueError: State A cannot be non-final, because it has a
            final output word. Only final states can have a final output
            word.
            sage: A.final_word_out = None
            sage: A.is_final = False
            sage: A.is_final
            False
        """
        return (self.final_word_out is not None)


    @is_final.setter
    def is_final(self, is_final):
        """
        Defines the state as a final state or a non-final state.

        INPUT:

        - ``is_final`` -- ``True`` if the state should be final and
          ``False`` otherwise.

        OUTPUT:

        Nothing.

        TESTS::

            sage: from sage.combinat.finite_state_machine import FSMState
            sage: A = FSMState('A', is_final=True)
            sage: A.final_word_out
            []
            sage: A.is_final = False
            sage: A.final_word_out is None
            True
            sage: A = FSMState('A', is_final=True, final_word_out='a')
            sage: A.is_final = False
            Traceback (most recent call last):
            ...
            ValueError: State A cannot be non-final, because it has a
            final output word. Only final states can have a final output
            word.
            sage: A = FSMState('A', is_final=True, final_word_out=[])
            sage: A.is_final = False
            sage: A.final_word_out is None
            True
        """
        if is_final and self.final_word_out is None:
            self._final_word_out_ = []
        elif not is_final:
            if not self.final_word_out:
                self._final_word_out_ = None
            else:
                raise ValueError("State %s cannot be non-final, because it "
                                 "has a final output word. Only final states "
                                 "can have a final output word. "
                                 % self.label())


    def label(self):
        """
        Returns the label of the state.

        INPUT:

        Nothing.

        OUTPUT:

        The label of the state.

        EXAMPLES::

            sage: from sage.combinat.finite_state_machine import FSMState
            sage: A = FSMState('state')
            sage: A.label()
            'state'
        """
        return self._label_


    def __copy__(self):
        """
        Returns a (shallow) copy of the state.

        INPUT:

        Nothing.

        OUTPUT:

        A new state.

        EXAMPLES::

            sage: from sage.combinat.finite_state_machine import FSMState
            sage: A = FSMState('A')
            sage: copy(A)
            'A'
        """
        new = FSMState(self.label(), self.word_out,
                       self.is_initial, self.is_final)
        if hasattr(self, 'hook'):
            new.hook = self.hook
        return new


    copy = __copy__


    def __deepcopy__(self, memo):
        """
        Returns a deep copy of the state.

        INPUT:

        - ``memo`` -- a dictionary storing already processed elements.

        OUTPUT:

        A new state.

        EXAMPLES::

            sage: from sage.combinat.finite_state_machine import FSMState
            sage: A = FSMState('A')
            sage: deepcopy(A)
            'A'
        """
        try:
            label = self._deepcopy_relabel_
        except AttributeError:
            label = deepcopy(self.label(), memo)
        new = FSMState(label, deepcopy(self.word_out, memo),
                       self.is_initial, self.is_final)
        if hasattr(self, 'hook'):
            new.hook = deepcopy(self.hook, memo)
        new.color = deepcopy(self.color, memo)
        new.final_word_out = deepcopy(self.final_word_out, memo)
        return new


    def deepcopy(self, memo=None):
        """
        Returns a deep copy of the state.

        INPUT:

        - ``memo`` -- (default: ``None``) a dictionary storing already
          processed elements.

        OUTPUT:

        A new state.

        EXAMPLES::

            sage: from sage.combinat.finite_state_machine import FSMState
            sage: A = FSMState((1, 3), color=[1, 2],
            ....:              is_final=True, final_word_out=3)
            sage: B = deepcopy(A)
            sage: B
            (1, 3)
            sage: B.label == A.label
            True
            sage: B.label is A.label
            False
            sage: B.color == A.color
            True
            sage: B.color is A.color
            False
            sage: B.is_final == A.is_final
            True
            sage: B.is_final is A.is_final
            True
            sage: B.final_word_out == A.final_word_out
            True
            sage: B.final_word_out is A.final_word_out
            False
        """
        return deepcopy(self, memo)


    def relabeled(self, label, memo=None):
        """
        Returns a deep copy of the state with a new label.

        INPUT:

        - ``label`` -- the label of new state.

        - ``memo`` -- (default: ``None``) a dictionary storing already
          processed elements.

        OUTPUT:

        A new state.

        EXAMPLES::

            sage: from sage.combinat.finite_state_machine import FSMState
            sage: A = FSMState('A')
            sage: A.relabeled('B')
            'B'

        """
        self._deepcopy_relabel_ = label
        new = deepcopy(self, memo)
        del self._deepcopy_relabel_
        return new


    def __hash__(self):
        """
        Returns a hash value for the object.

        INPUT:

        Nothing.

        OUTPUT:

        The hash of this state.

        TESTS::

            sage: from sage.combinat.finite_state_machine import FSMState
            sage: A = FSMState('A')
            sage: hash(A) #random
            -269909568
        """
        return hash(self.label())


    def _repr_(self):
        """
        Returns the string "label".

        INPUT:

        Nothing.

        OUTPUT:

        A string.

        TESTS:

            sage: from sage.combinat.finite_state_machine import FSMState
            sage: FSMState('A')._repr_()
            "'A'"
        """
        return repr(self.label())


    def __eq__(left, right):
        """
        Returns True if two states are the same, i.e., if they have
        the same labels.

        Note that the hooks and whether the states are initial or
        final are not checked.

        INPUT:

        - ``left`` -- a state.

        - ``right`` -- a state.

        OUTPUT:

        True or False.

        EXAMPLES::

            sage: from sage.combinat.finite_state_machine import FSMState
            sage: A = FSMState('A')
            sage: B = FSMState('A', is_initial=True)
            sage: A == B
            True
        """
        if not is_FSMState(right):
            return False
        return left.label() == right.label()


    def __ne__(left, right):
        """
        Tests for inequality, complement of __eq__.

        INPUT:

        - ``left`` -- a state.

        - ``right`` -- a state.

        OUTPUT:

        True or False.

        EXAMPLES::

            sage: from sage.combinat.finite_state_machine import FSMState
            sage: A = FSMState('A', is_initial=True)
            sage: B = FSMState('A', is_final=True)
            sage: A != B
            False
        """
        return (not (left == right))


    def __nonzero__(self):
        """
        Returns True.

        INPUT:

        Nothing.

        OUTPUT:

        True or False.

        TESTS::

            sage: from sage.combinat.finite_state_machine import FSMState
            sage: FSMState('A').__nonzero__()
            True
        """
        return True  # A state cannot be zero (see __init__)


#*****************************************************************************


def is_FSMTransition(T):
    """
    Tests whether or not ``T`` inherits from :class:`FSMTransition`.

    TESTS::

        sage: from sage.combinat.finite_state_machine import is_FSMTransition, FSMTransition
        sage: is_FSMTransition(FSMTransition('A', 'B'))
        True
    """
    return isinstance(T, FSMTransition)


class FSMTransition(SageObject):
    """
    Class for a transition of a finite state machine.

    INPUT:

    - ``from_state`` -- state from which transition starts.

    - ``to_state`` -- state in which transition ends.

    - ``word_in`` -- the input word of the transitions (when the
      finite state machine is used as automaton)

    - ``word_out`` -- the output word of the transitions (when the
      finite state machine is used as transducer)

    OUTPUT:

    A transition of a finite state machine.

    EXAMPLES::

        sage: from sage.combinat.finite_state_machine import FSMState, FSMTransition
        sage: A = FSMState('A')
        sage: B = FSMState('B')
        sage: S = FSMTransition(A, B, 0, 1)
        sage: T = FSMTransition('A', 'B', 0, 1)
        sage: T == S
        True
        sage: U = FSMTransition('A', 'B', 0)
        sage: U == T
        False

    """
    def __init__(self, from_state, to_state,
                 word_in=None, word_out=None,
                 hook=None):
        """
        See :class:`FSMTransition` for more information.

        EXAMPLES::

            sage: from sage.combinat.finite_state_machine import FSMTransition
            sage: FSMTransition('A', 'B', 0, 1)
            Transition from 'A' to 'B': 0|1
        """
        if is_FSMState(from_state):
            self.from_state = from_state
        else:
            self.from_state = FSMState(from_state)
        if is_FSMState(to_state):
            self.to_state = to_state
        else:
            self.to_state = FSMState(to_state)

        if isinstance(word_in, list):
            self.word_in = word_in
        elif word_in is not None:
            self.word_in = [word_in]
        else:
            self.word_in = []

        if isinstance(word_out, list):
            self.word_out = word_out
        elif word_out is not None:
            self.word_out = [word_out]
        else:
            self.word_out = []

        if hook is not None:
            if hasattr(hook, '__call__'):
                self.hook = hook
            else:
                raise TypeError('Wrong argument for hook.')


    def __lt__(self, other):
        """
        Returns True if ``self`` is less than ``other`` with respect to the
        key ``(self.from_state, self.word_in, self.to_state, self.word_out)``.

        INPUT:

        - `other` -- a transition.

        OUTPUT:

        True or False.

        EXAMPLE::

            sage: from sage.combinat.finite_state_machine import FSMTransition
            sage: FSMTransition(0,1,0,0) < FSMTransition(1,0,0,0)
            True
        """
        return (self.from_state, self.word_in, self.to_state, self.word_out) < \
            (other.from_state, other.word_in, other.to_state, other.word_out)


    def __copy__(self):
        """
        Returns a (shallow) copy of the transition.

        INPUT:

        Nothing.

        OUTPUT:

        A new transition.

        EXAMPLES::

            sage: from sage.combinat.finite_state_machine import FSMTransition
            sage: t = FSMTransition('A', 'B', 0)
            sage: copy(t)
            Transition from 'A' to 'B': 0|-
        """
        new = FSMTransition(self.from_state, self.to_state,
                            self.word_in, self.word_out)
        if hasattr(self, 'hook'):
            new.hook = self.hook
        return new


    copy = __copy__

    def __deepcopy__(self, memo):
        """
        Returns a deep copy of the transition.

        INPUT:

        - ``memo`` -- a dictionary storing already processed elements.

        OUTPUT:

        A new transition.

        EXAMPLES::

            sage: from sage.combinat.finite_state_machine import FSMTransition
            sage: t = FSMTransition('A', 'B', 0)
            sage: deepcopy(t)
            Transition from 'A' to 'B': 0|-
        """
        new = FSMTransition(deepcopy(self.from_state, memo),
                            deepcopy(self.to_state, memo),
                            deepcopy(self.word_in, memo),
                            deepcopy(self.word_out, memo))
        if hasattr(self, 'hook'):
            new.hook = deepcopy(self.hook, memo)
        return new


    def deepcopy(self, memo=None):
        """
        Returns a deep copy of the transition.

        INPUT:

        - ``memo`` -- (default: ``None``) a dictionary storing already
          processed elements.

        OUTPUT:

        A new transition.

        EXAMPLES::

            sage: from sage.combinat.finite_state_machine import FSMTransition
            sage: t = FSMTransition('A', 'B', 0)
            sage: deepcopy(t)
            Transition from 'A' to 'B': 0|-
        """
        return deepcopy(self, memo)


    def __hash__(self):
        """
        Since transitions are mutable, they should not be hashable, so
        we return a type error.

        INPUT:

        Nothing.

        OUTPUT:

        The hash of this transition.

        EXAMPLES::

            sage: from sage.combinat.finite_state_machine import FSMTransition
            sage: hash(FSMTransition('A', 'B'))
            Traceback (most recent call last):
            ...
            TypeError: Transitions are mutable, and thus not hashable.

        """
        raise TypeError("Transitions are mutable, and thus not hashable.")


    def _repr_(self):
        """
        Represents a transitions as from state to state and input, output.

        INPUT:

        Nothing.

        OUTPUT:

        A string.

        EXAMPLES::

            sage: from sage.combinat.finite_state_machine import FSMTransition
            sage: FSMTransition('A', 'B', 0, 0)._repr_()
            "Transition from 'A' to 'B': 0|0"

        """
        return "Transition from %s to %s: %s" % (repr(self.from_state),
                                                 repr(self.to_state),
                                                 self._in_out_label_())


    def _in_out_label_(self):
        """
        Returns the input and output of a transition as
        "word_in|word_out".

        INPUT:

        Nothing.

        OUTPUT:

        A string of the input and output labels.

        EXAMPLES::

            sage: from sage.combinat.finite_state_machine import FSMTransition
            sage: FSMTransition('A', 'B', 0, 1)._in_out_label_()
            '0|1'
        """
        return "%s|%s" % (FSMWordSymbol(self.word_in),
                          FSMWordSymbol(self.word_out))


    def __eq__(left, right):
        """
        Returns True if the two transitions are the same, i.e., if the
        both go from the same states to the same states and read and
        write the same words.

        Note that the hooks are not checked.

        INPUT:

        - ``left`` -- a transition.

        - ``right`` -- a transition.

        OUTPUT:

        True or False.

        EXAMPLES::

            sage: from sage.combinat.finite_state_machine import FSMState, FSMTransition
            sage: A = FSMState('A', is_initial=True)
            sage: t1 = FSMTransition('A', 'B', 0, 1)
            sage: t2 = FSMTransition(A, 'B', 0, 1)
            sage: t1 == t2
            True
        """
        if not is_FSMTransition(right):
            raise TypeError('Only instances of FSMTransition ' \
                'can be compared.')
        return left.from_state == right.from_state \
            and left.to_state == right.to_state \
            and left.word_in == right.word_in \
            and left.word_out == right.word_out


    def __ne__(left, right):
        """

        INPUT:

        - ``left`` -- a transition.

        - ``right`` -- a transition.

        OUTPUT:

        True or False.
        Tests for inequality, complement of __eq__.

        EXAMPLES::

            sage: from sage.combinat.finite_state_machine import FSMState, FSMTransition
            sage: A = FSMState('A', is_initial=True)
            sage: t1 = FSMTransition('A', 'B', 0, 1)
            sage: t2 = FSMTransition(A, 'B', 0, 1)
            sage: t1 != t2
            False
        """
        return (not (left == right))


    def __nonzero__(self):
        """
        Returns True.

        INPUT:

        Nothing.

        OUTPUT:

        True or False.

        EXAMPLES::

            sage: from sage.combinat.finite_state_machine import FSMTransition
            sage: FSMTransition('A', 'B', 0).__nonzero__()
            True
        """
        return True  # A transition cannot be zero (see __init__)


#*****************************************************************************


def is_FiniteStateMachine(FSM):
    """
    Tests whether or not ``FSM`` inherits from :class:`FiniteStateMachine`.

    TESTS::

        sage: from sage.combinat.finite_state_machine import is_FiniteStateMachine
        sage: is_FiniteStateMachine(FiniteStateMachine())
        True
        sage: is_FiniteStateMachine(Automaton())
        True
        sage: is_FiniteStateMachine(Transducer())
        True
    """
    return isinstance(FSM, FiniteStateMachine)


def duplicate_transition_ignore(old_transition, new_transition):
    """
    Default function for handling duplicate transitions in finite
    state machines. This implementation ignores the occurrence.

    See the documentation of the ``on_duplicate_transition`` parameter
    of :class:`FiniteStateMachine`.

    INPUT:

    - ``old_transition`` -- A transition in a finite state machine.

    - ``new_transition`` -- A transition, identical to ``old_transition``,
      which is to be inserted into the finite state machine.

    OUTPUT:

    The same transition, unchanged.

    EXAMPLES::

        sage: from sage.combinat.finite_state_machine import duplicate_transition_ignore
        sage: from sage.combinat.finite_state_machine import FSMTransition
        sage: duplicate_transition_ignore(FSMTransition(0, 0, 1),
        ....:                             FSMTransition(0, 0, 1))
        Transition from 0 to 0: 1|-
    """
    return old_transition


def duplicate_transition_raise_error(old_transition, new_transition):
    """
    Alternative function for handling duplicate transitions in finite
    state machines. This implementation raises a ``ValueError``.

    See the documentation of the ``on_duplicate_transition`` parameter
    of :class:`FiniteStateMachine`.

    INPUT:

    - ``old_transition`` -- A transition in a finite state machine.

    - ``new_transition`` -- A transition, identical to ``old_transition``,
      which is to be inserted into the finite state machine.

    OUTPUT:

    Nothing. A ``ValueError`` is raised.

    EXAMPLES::

        sage: from sage.combinat.finite_state_machine import duplicate_transition_raise_error
        sage: from sage.combinat.finite_state_machine import FSMTransition
        sage: duplicate_transition_raise_error(FSMTransition(0, 0, 1),
        ....:                                  FSMTransition(0, 0, 1))
        Traceback (most recent call last):
        ...
        ValueError: Attempting to re-insert transition Transition from 0 to 0: 1|-
    """
    raise ValueError("Attempting to re-insert transition %s" % old_transition)


def duplicate_transition_add_input(old_transition, new_transition):
    """
    Alternative function for handling duplicate transitions in finite
    state machines. This implementation adds the input label of the
    new transition to the input label of the old transition.  This is
    intended for the case where a Markov chain is modelled by a finite
    state machine using the input labels as transition probabilities.

    See the documentation of the ``on_duplicate_transition`` parameter
    of :class:`FiniteStateMachine`.

    INPUT:

    - ``old_transition`` -- A transition in a finite state machine.

    - ``new_transition`` -- A transition, identical to ``old_transition``,
      which is to be inserted into the finite state machine.

    OUTPUT:

    A transition whose input weight is the sum of the input
    weights of ``old_transition`` and ``new_transition``.

    EXAMPLES::

        sage: from sage.combinat.finite_state_machine import duplicate_transition_add_input
        sage: from sage.combinat.finite_state_machine import FSMTransition
        sage: duplicate_transition_add_input(FSMTransition('a', 'a', 1/2),
        ....:                                FSMTransition('a', 'a', 1/2))
        Transition from 'a' to 'a': 1|-

    Input labels must be lists of length 1::

        sage: duplicate_transition_add_input(FSMTransition('a', 'a', [1, 1]),
        ....:                                FSMTransition('a', 'a', [1, 1]))
        Traceback (most recent call last):
        ...
        TypeError: Trying to use duplicate_transition_add_input on
        "Transition from 'a' to 'a': 1,1|-" and
        "Transition from 'a' to 'a': 1,1|-",
        but input words are assumed to be lists of length 1
    """
    if (hasattr(old_transition.word_in, '__iter__')
        and len(old_transition.word_in) == 1
        and hasattr(new_transition.word_in, '__iter__')
        and len(new_transition.word_in) == 1):
        old_transition.word_in = [old_transition.word_in[0]
                                  + new_transition.word_in[0]]
    else:
        raise TypeError('Trying to use duplicate_transition_add_input on ' +
                        '"%s" and "%s", ' % (old_transition, new_transition) +
                        'but input words are assumed to be lists of length 1')
    return old_transition


class FiniteStateMachine(SageObject):
    """
    Class for a finite state machine.

    A finite state machine is a finite set of states connected by
    transitions.

    INPUT:

    - ``data`` -- can be any of the following:

      #. a dictionary of dictionaries (of transitions),

      #. a dictionary of lists (of states or transitions),

      #. a list (of transitions),

      #. a function (transition function),

      #. an other instance of a finite state machine.

    - ``initial_states`` and ``final_states`` -- the initial and
      final states of this machine

    - ``input_alphabet`` and ``output_alphabet`` -- the input and
      output alphabets of this machine

    - ``determine_alphabets`` -- If ``True``, then the function
      :meth:`.determine_alphabets` is called after ``data`` was read and
      processed, if ``False``, then not. If it is ``None``, then it is
      decided during the construction of the finite state machine
      whether :meth:`.determine_alphabets` should be called.

    - ``with_final_word_out`` -- If given (not ``None``), then the
      function :meth:`.with_final_word_out` (more precisely, its inplace
      pendant :meth:`.construct_final_word_out`) is called with input
      ``letters=with_final_word_out`` at the end of the creation
      process.

    - ``store_states_dict`` -- If ``True``, then additionally the states
      are stored in an interal dictionary for speed up.

    - ``on_duplicate_transition`` -- A function which is called when a
      transition is inserted into ``self`` which already existed (same
      ``from_state``, same ``to_state``, same ``word_in``, same ``word_out``).

      This function is assumed to take two arguments, the first being
      the already existing transition, the second being the new
      transition (as an ``FSMTransition``). The function must return the
      (possibly modified) original transition.

      By default, we have ``on_duplicate_transition=None``, which is
      interpreted as
      ``on_duplicate_transition=duplicate_transition_ignore``, where
      ``duplicate_transition_ignore`` is a predefined function
      ignoring the occurrence. Other such predefined functions are
      ``duplicate_transition_raise_error`` and
      ``duplicate_transition_add_input``.

    OUTPUT:

    A finite state machine.

    The object creation of :class:`Automaton` and :class:`Transducer`
    is the same as the one described here (i.e. just replace the word
    ``FiniteStateMachine`` by ``Automaton`` or ``Transducer``).

    Each transition of an automaton has an input label. Automata can,
    for example, be determinised (see
    :meth:`Automaton.determinisation`) and minimized (see
    :meth:`Automaton.minimization`). Each transition of a transducer
    has an input and an output label. Transducers can, for example, be
    simplified (see :meth:`Transducer.simplification`).

    EXAMPLES::

        sage: from sage.combinat.finite_state_machine import FSMState, FSMTransition

    See documentation for more examples.

    We illustrate the different input formats:

    #.  The input-data can be a dictionary of dictionaries, where

        - the keys of the outer dictionary are state-labels (from-states of
          transitions),
        - the keys of the inner dictionaries are state-labels (to-states of
          transitions),
        - the values of the inner dictionaries specify the transition
          more precisely.

        The easiest is to use a tuple consisting of an input and an
        output word::

            sage: FiniteStateMachine({'a':{'b':(0, 1), 'c':(1, 1)}})
            Finite state machine with 3 states

        Instead of the tuple anything iterable (e.g. a list) can be
        used as well.

        If you want to use the arguments of ``FSMTransition``
        directly, you can use a dictionary::

            sage: FiniteStateMachine({'a':{'b':{'word_in':0, 'word_out':1},
            ....:                          'c':{'word_in':1, 'word_out':1}}})
            Finite state machine with 3 states

        In the case you already have instances of ``FSMTransition``, it is
        possible to use them directly::

            sage: FiniteStateMachine({'a':{'b':FSMTransition('a', 'b', 0, 1),
            ....:                          'c':FSMTransition('a', 'c', 1, 1)}})
            Finite state machine with 3 states

    #.  The input-data can be a dictionary of lists, where the keys
        are states or label of states.

        The list-elements can be states::

            sage: a = FSMState('a')
            sage: b = FSMState('b')
            sage: c = FSMState('c')
            sage: FiniteStateMachine({a:[b, c]})
            Finite state machine with 3 states

        Or the list-elements can simply be labels of states::

            sage: FiniteStateMachine({'a':['b', 'c']})
            Finite state machine with 3 states

        The list-elements can also be transitions::

            sage: FiniteStateMachine({'a':[FSMTransition('a', 'b', 0, 1),
            ....:                          FSMTransition('a', 'c', 1, 1)]})
            Finite state machine with 3 states

        Or they can be tuples of a label, an input word and an output
        word specifying a transition::

            sage: FiniteStateMachine({'a':[('b', 0, 1), ('c', 1, 1)]})
            Finite state machine with 3 states

    #.  The input-data can be a list, where its elements specify
        transitions::

            sage: FiniteStateMachine([FSMTransition('a', 'b', 0, 1),
            ....:                     FSMTransition('a', 'c', 1, 1)])
            Finite state machine with 3 states

        It is possible to skip ``FSMTransition`` in the example above::

            sage: FiniteStateMachine([('a', 'b', 0, 1), ('a', 'c', 1, 1)])
            Finite state machine with 3 states

        The parameters of the transition are given in tuples. Anyhow,
        anything iterable (e.g. a list) is possible.

        You can also name the parameters of the transition. For this
        purpose you take a dictionary::

            sage: FiniteStateMachine([{'from_state':'a', 'to_state':'b',
            ....:                      'word_in':0, 'word_out':1},
            ....:                     {'from_state':'a', 'to_state':'c',
            ....:                      'word_in':1, 'word_out':1}])
            Finite state machine with 3 states

        Other arguments, which :class:`FSMTransition` accepts, can be
        added, too.

    #.  The input-data can also be function acting as transition
        function:

        This function has two input arguments:

        #. a label of a state (from which the transition starts),

        #. a letter of the (input-)alphabet (as input-label of the transition).

        It returns a tuple with the following entries:

        #. a label of a state (to which state the transition goes),

        #. a letter of or a word over the (output-)alphabet (as
           output-label of the transition).

        It may also output a list of such tuples if several
        transitions from the from-state and the input letter exist
        (this means that the finite state machine is
        non-deterministic).

        If the transition does not exist, the function should raise a
        ``LookupError`` or return an empty list.

        When constructing a finite state machine in this way, some
        inital states and an input alphabet have to be specified.

        ::

            sage: def f(state_from, read):
            ....:     if int(state_from) + read <= 2:
            ....:         state_to = 2*int(state_from)+read
            ....:         write = 0
            ....:     else:
            ....:         state_to = 2*int(state_from) + read - 5
            ....:         write = 1
            ....:     return (str(state_to), write)
            sage: F = FiniteStateMachine(f, input_alphabet=[0, 1],
            ....:                        initial_states=['0'],
            ....:                        final_states=['0'])
            sage: F([1, 0, 1])
            (True, '0', [0, 0, 1])

    #.  The input-data can be an other instance of a finite state machine::

            sage: FiniteStateMachine(FiniteStateMachine([]))
            Traceback (most recent call last):
            ...
            NotImplementedError


    The following examples demonstrate the use of ``on_duplicate_transition``::

        sage: F = FiniteStateMachine([['a', 'a', 1/2], ['a', 'a', 1/2]])
        sage: F.transitions()
        [Transition from 'a' to 'a': 1/2|-]

    ::

        sage: from sage.combinat.finite_state_machine import duplicate_transition_raise_error
        sage: F1 = FiniteStateMachine([['a', 'a', 1/2], ['a', 'a', 1/2]],
        ....:                         on_duplicate_transition=duplicate_transition_raise_error)
        Traceback (most recent call last):
        ...
        ValueError: Attempting to re-insert transition Transition from 'a' to 'a': 1/2|-

    Use ``duplicate_transition_add_input`` to emulate a Markov chain,
    the input labels are considered as transition probabilities::

        sage: from sage.combinat.finite_state_machine import duplicate_transition_add_input
        sage: F = FiniteStateMachine([['a', 'a', 1/2], ['a', 'a', 1/2]],
        ....:                        on_duplicate_transition=duplicate_transition_add_input)
        sage: F.transitions()
        [Transition from 'a' to 'a': 1|-]

    Use ``with_final_word_out`` to construct final output::

        sage: T = Transducer([(0, 1, 0, 0), (1, 0, 0, 0)],
        ....:                initial_states=[0],
        ....:                final_states=[0],
        ....:                with_final_word_out=0)
        sage: for s in T.iter_final_states():
        ....:     print s, s.final_word_out
        0 []
        1 [0]

    TESTS::

        sage: a = FSMState('S_a', 'a')
        sage: b = FSMState('S_b', 'b')
        sage: c = FSMState('S_c', 'c')
        sage: d = FSMState('S_d', 'd')
        sage: FiniteStateMachine({a:[b, c], b:[b, c, d],
        ....:                     c:[a, b], d:[a, c]})
        Finite state machine with 4 states

    We have several constructions which lead to the same finite
    state machine::

        sage: A = FSMState('A')
        sage: B = FSMState('B')
        sage: C = FSMState('C')
        sage: FSM1 = FiniteStateMachine(
        ....:  {A:{B:{'word_in':0, 'word_out':1},
        ....:   C:{'word_in':1, 'word_out':1}}})
        sage: FSM2 = FiniteStateMachine({A:{B:(0, 1), C:(1, 1)}})
        sage: FSM3 = FiniteStateMachine(
        ....:  {A:{B:FSMTransition(A, B, 0, 1),
        ....:      C:FSMTransition(A, C, 1, 1)}})
        sage: FSM4 = FiniteStateMachine({A:[(B, 0, 1), (C, 1, 1)]})
        sage: FSM5 = FiniteStateMachine(
        ....:  {A:[FSMTransition(A, B, 0, 1), FSMTransition(A, C, 1, 1)]})
        sage: FSM6 = FiniteStateMachine(
        ....:  [{'from_state':A, 'to_state':B, 'word_in':0, 'word_out':1},
        ....:   {'from_state':A, 'to_state':C, 'word_in':1, 'word_out':1}])
        sage: FSM7 = FiniteStateMachine([(A, B, 0, 1), (A, C, 1, 1)])
        sage: FSM8 = FiniteStateMachine(
        ....:  [FSMTransition(A, B, 0, 1), FSMTransition(A, C, 1, 1)])

        sage: FSM1 == FSM2 == FSM3 == FSM4 == FSM5 == FSM6 == FSM7 == FSM8
        True

    It is possible to skip ``FSMTransition`` in the example above.

    Some more tests for different input-data::

        sage: FiniteStateMachine({'a':{'a':[0, 0], 'b':[1, 1]},
        ....:                     'b':{'b':[1, 0]}})
        Finite state machine with 2 states

        sage: a = FSMState('S_a', 'a')
        sage: b = FSMState('S_b', 'b')
        sage: c = FSMState('S_c', 'c')
        sage: d = FSMState('S_d', 'd')
        sage: t1 = FSMTransition(a, b)
        sage: t2 = FSMTransition(b, c)
        sage: t3 = FSMTransition(b, d)
        sage: t4 = FSMTransition(c, d)
        sage: FiniteStateMachine([t1, t2, t3, t4])
        Finite state machine with 4 states
    """

    #*************************************************************************
    # init
    #*************************************************************************


    def __init__(self,
                 data=None,
                 initial_states=None, final_states=None,
                 input_alphabet=None, output_alphabet=None,
                 determine_alphabets=None,
                 with_final_word_out=None,
                 store_states_dict=True,
                 on_duplicate_transition=None):
        """
        See :class:`FiniteStateMachine` for more information.

        TEST::

            sage: FiniteStateMachine()
            Finite state machine with 0 states
        """
        self._states_ = []  # List of states in the finite state
                            # machine.  Each state stores a list of
                            # outgoing transitions.
        if store_states_dict:
            self._states_dict_ = {}

        if initial_states is not None:
            if not hasattr(initial_states, '__iter__'):
                raise TypeError('Initial states must be iterable ' \
                    '(e.g. a list of states).')
            for s in initial_states:
                state = self.add_state(s)
                state.is_initial = True

        if final_states is not None:
            if not hasattr(final_states, '__iter__'):
                raise TypeError('Final states must be iterable ' \
                    '(e.g. a list of states).')
            for s in final_states:
                state = self.add_state(s)
                state.is_final = True

        self.input_alphabet = input_alphabet
        self.output_alphabet = output_alphabet

        if on_duplicate_transition is None:
            on_duplicate_transition = duplicate_transition_ignore
        if hasattr(on_duplicate_transition, '__call__'):
            self.on_duplicate_transition=on_duplicate_transition
        else:
            raise TypeError, 'on_duplicate_transition must be callable'

        if data is None:
            pass
        elif is_FiniteStateMachine(data):
            raise NotImplementedError
        elif hasattr(data, 'iteritems'):
            # data is a dict (or something similar),
            # format: key = from_state, value = iterator of transitions
            for (sf, iter_transitions) in data.iteritems():
                self.add_state(sf)
                if hasattr(iter_transitions, 'iteritems'):
                    for (st, transition) in iter_transitions.iteritems():
                        self.add_state(st)
                        if is_FSMTransition(transition):
                            self.add_transition(transition)
                        elif hasattr(transition, 'iteritems'):
                            self.add_transition(sf, st, **transition)
                        elif hasattr(transition, '__iter__'):
                            self.add_transition(sf, st, *transition)
                        else:
                            self.add_transition(sf, st, transition)
                elif hasattr(iter_transitions, '__iter__'):
                    for transition in iter_transitions:
                        if hasattr(transition, '__iter__'):
                            L = [sf]
                            L.extend(transition)
                        elif is_FSMTransition(transition):
                            L = transition
                        else:
                            L = [sf, transition]
                        self.add_transition(L)
                else:
                    raise TypeError('Wrong input data for transition.')
            if determine_alphabets is None and input_alphabet is None \
                    and output_alphabet is None:
                determine_alphabets = True
        elif hasattr(data, '__iter__'):
            # data is a something that is iterable,
            # items are transitions
            for transition in data:
                if is_FSMTransition(transition):
                    self.add_transition(transition)
                elif hasattr(transition, 'iteritems'):
                    self.add_transition(transition)
                elif hasattr(transition, '__iter__'):
                    self.add_transition(transition)
                else:
                    raise TypeError('Wrong input data for transition.')
            if determine_alphabets is None and input_alphabet is None \
                    and output_alphabet is None:
                determine_alphabets = True
        elif hasattr(data, '__call__'):
            self.add_from_transition_function(data)
        else:
            raise TypeError('Cannot decide what to do with data.')

        if determine_alphabets:
            self.determine_alphabets()

        if with_final_word_out is not None:
            self.construct_final_word_out(with_final_word_out)


    #*************************************************************************
    # copy and hash
    #*************************************************************************


    def __copy__(self):
        """
        Returns a (shallow) copy of the finite state machine.

        INPUT:

        Nothing.

        OUTPUT:

        A new finite state machine.

        TESTS::

            sage: copy(FiniteStateMachine())
            Traceback (most recent call last):
            ...
            NotImplementedError
        """
        raise NotImplementedError


    copy = __copy__


    def empty_copy(self, memo=None):
        """
        Returns an empty deep copy of the finite state machine, i.e.,
        ``input_alphabet``, ``output_alphabet``, ``on_duplicate_transition``
        are preserved, but states and transitions are not.

        INPUT:

        - ``memo`` -- a dictionary storing already processed elements.

        OUTPUT:

        A new finite state machine.

        EXAMPLES::

            sage: from sage.combinat.finite_state_machine import duplicate_transition_raise_error
            sage: F = FiniteStateMachine([('A', 'A', 0, 2), ('A', 'A', 1, 3)],
            ....:                        input_alphabet=[0, 1],
            ....:                        output_alphabet=[2, 3],
            ....:                        on_duplicate_transition=duplicate_transition_raise_error)
            sage: FE = F.empty_copy(); FE
            Finite state machine with 0 states
            sage: FE.input_alphabet
            [0, 1]
            sage: FE.output_alphabet
            [2, 3]
            sage: FE.on_duplicate_transition == duplicate_transition_raise_error
            True
        """
        new = self.__class__()
        new.input_alphabet = deepcopy(self.input_alphabet, memo)
        new.output_alphabet = deepcopy(self.output_alphabet, memo)
        new.on_duplicate_transition = self.on_duplicate_transition
        return new


    def __deepcopy__(self, memo):
        """
        Returns a deep copy of the finite state machine.

        INPUT:

        - ``memo`` -- a dictionary storing already processed elements.

        OUTPUT:

        A new finite state machine.

        EXAMPLES::

            sage: F = FiniteStateMachine([('A', 'A', 0, 1), ('A', 'A', 1, 0)])
            sage: deepcopy(F)
            Finite state machine with 1 states
        """
        relabel = hasattr(self, '_deepcopy_relabel_')
        new = self.empty_copy(memo=memo)
        relabel_iter = itertools.count(0)
        for state in self.iter_states():
            if relabel:
                if self._deepcopy_labels_ is None:
                    state._deepcopy_relabel_ = relabel_iter.next()
                elif hasattr(self._deepcopy_labels_, '__call__'):
                    state._deepcopy_relabel_ = self._deepcopy_labels_(state.label())
                elif hasattr(self._deepcopy_labels_, '__getitem__'):
                    state._deepcopy_relabel_ = self._deepcopy_labels_[state.label()]
                else:
                    raise TypeError("labels must be None, a callable "
                                    "or a dictionary.")
            s = deepcopy(state, memo)
            if relabel:
                del state._deepcopy_relabel_
            new.add_state(s)
        for transition in self.iter_transitions():
            new.add_transition(deepcopy(transition, memo))
        return new


    def deepcopy(self, memo=None):
        """
        Returns a deep copy of the finite state machine.

        INPUT:

        - ``memo`` -- (default: ``None``) a dictionary storing already
          processed elements.

        OUTPUT:

        A new finite state machine.

        EXAMPLES::

            sage: F = FiniteStateMachine([('A', 'A', 0, 1), ('A', 'A', 1, 0)])
            sage: deepcopy(F)
            Finite state machine with 1 states
        """
        return deepcopy(self, memo)


    def relabeled(self, memo=None, labels=None):
        """
        Returns a deep copy of the finite state machine, but the
        states are relabeled.

        INPUT:

        - ``memo`` -- (default: ``None``) a dictionary storing already
          processed elements.

        - ``labels`` -- (default: ``None``) a dictionary or callable
          mapping old labels to new labels. If ``None``, then the new
          labels are integers starting with 0.

        OUTPUT:

        A new finite state machine.

        EXAMPLES::

            sage: FSM1 = FiniteStateMachine([('A', 'B'), ('B', 'C'), ('C', 'A')])
            sage: FSM1.states()
            ['A', 'B', 'C']
            sage: FSM2 = FSM1.relabeled()
            sage: FSM2.states()
            [0, 1, 2]
            sage: FSM3 = FSM1.relabeled(labels={'A': 'a', 'B': 'b', 'C': 'c'})
            sage: FSM3.states()
            ['a', 'b', 'c']
            sage: FSM4 = FSM2.relabeled(labels=lambda x: 2*x)
            sage: FSM4.states()
            [0, 2, 4]

        TESTS::

            sage: FSM2.relabeled(labels=1)
            Traceback (most recent call last):
            ...
            TypeError: labels must be None, a callable or a dictionary.
        """
        self._deepcopy_relabel_ = True
        self._deepcopy_labels_ = labels
        new = deepcopy(self, memo)
        del self._deepcopy_relabel_
        del self._deepcopy_labels_
        return new


    def induced_sub_finite_state_machine(self, states):
        """
        Returns a sub-finite-state-machine of the finite state machine
        induced by the given states.

        INPUT:

        - ``states`` -- a list (or an iterator) of states (either labels or
          instances of :class:`FSMState`) of the sub-finite-state-machine.

        OUTPUT:

        A new finite state machine. It consists (of deep copies) of
        the given states and (deep copies) of all transitions of ``self``
        between these states.

        EXAMPLE::

            sage: FSM = FiniteStateMachine([(0, 1, 0), (0, 2, 0),
            ....:                           (1, 2, 0), (2, 0, 0)])
            sage: sub_FSM = FSM.induced_sub_finite_state_machine([0, 1])
            sage: sub_FSM.states()
            [0, 1]
            sage: sub_FSM.transitions()
            [Transition from 0 to 1: 0|-]
            sage: FSM.induced_sub_finite_state_machine([3])
            Traceback (most recent call last):
            ...
            ValueError: 3 is not a state of this finite state machine.

        TESTS:

        Make sure that the links between transitions and states
        are still intact::

            sage: sub_FSM.transitions()[0].from_state is sub_FSM.state(0)
            True

        """
        good_states = set()
        for state in states:
            if not self.has_state(state):
                raise ValueError("%s is not a state of this finite state machine." % state)
            good_states.add(self.state(state))

        memo = {}
        new = self.empty_copy(memo=memo)
        for state in good_states:
            s = deepcopy(state, memo)
            new.add_state(s)

        for state in good_states:
            for transition in self.iter_transitions(state):
                if transition.to_state in good_states:
                    new.add_transition(deepcopy(transition, memo))

        return new


    def __hash__(self):
        """
        Since finite state machines are mutable, they should not be
        hashable, so we return a type error.

        INPUT:

        Nothing.

        OUTPUT:

        The hash of this finite state machine.

        EXAMPLES::

            sage: hash(FiniteStateMachine())
            Traceback (most recent call last):
            ...
            TypeError: Finite state machines are mutable, and thus not hashable.
        """
        if getattr(self, "_immutable", False):
            return hash((tuple(self.states()), tuple(self.transitions())))
        raise TypeError("Finite state machines are mutable, " \
            "and thus not hashable.")


    #*************************************************************************
    # operators
    #*************************************************************************


    def __or__(self, other):
        """
        Returns the disjoint union of the finite state machines self and other.

        INPUT:

        - ``other`` -- a finite state machine.

        OUTPUT:

        A new finite state machine.

        TESTS::

            sage: FiniteStateMachine() | FiniteStateMachine([('A', 'B')])
            Traceback (most recent call last):
            ...
            NotImplementedError
        """
        if is_FiniteStateMachine(other):
            return self.disjoint_union(other)

    __add__ = __or__

    def __iadd__(self, other):
        """
        TESTS::

            sage: F = FiniteStateMachine()
            sage: F += FiniteStateMachine()
            Traceback (most recent call last):
            ...
            NotImplementedError
        """
        raise NotImplementedError


    def __and__(self, other):
        """
        Returns the intersection of ``self`` with ``other``.

        TESTS::

            sage: FiniteStateMachine() & FiniteStateMachine([('A', 'B')])
            Traceback (most recent call last):
            ...
            NotImplementedError
        """
        if is_FiniteStateMachine(other):
            return self.intersection(other)


    def __imul__(self, other):
        """
        TESTS::

            sage: F = FiniteStateMachine()
            sage: F *= FiniteStateMachine()
            Traceback (most recent call last):
            ...
            NotImplementedError
        """
        raise NotImplementedError


    def __call__(self, *args, **kwargs):
        """
        .. WARNING::

            The default output of this method is scheduled to change.
            This docstring describes the new default behaviour, which can
            already be achieved by setting
            ``FSMOldProcessOutput`` to ``False``.

        Calls either method :meth:`.composition` or :meth:`.process`
        (with ``full_output=False``).

        By setting ``FSMOldProcessOutput`` to ``False``
        the new desired output is produced.

        EXAMPLES::

            sage: sage.combinat.finite_state_machine.FSMOldProcessOutput = False  # activate new output behavior
            sage: from sage.combinat.finite_state_machine import FSMState
            sage: A = FSMState('A', is_initial=True, is_final=True)
            sage: binary_inverter = Transducer({A:[(A, 0, 1), (A, 1, 0)]})
            sage: binary_inverter([0, 1, 0, 0, 1, 1])
            [1, 0, 1, 1, 0, 0]

        ::

            sage: F = Transducer([('A', 'B', 1, 0), ('B', 'B', 1, 1),
            ....:                 ('B', 'B', 0, 0)],
            ....:                initial_states=['A'], final_states=['B'])
            sage: G = Transducer([(1, 1, 0, 0), (1, 2, 1, 0),
            ....:                 (2, 2, 0, 1), (2, 1, 1, 1)],
            ....:                initial_states=[1], final_states=[1])
            sage: H = G(F)
            sage: H.states()
            [('A', 1), ('B', 1), ('B', 2)]
        """
        if len(args) == 0:
            raise TypeError("Called with too few arguments.")
        if is_FiniteStateMachine(args[0]):
            return self.composition(*args, **kwargs)
        if hasattr(args[0], '__iter__'):
            if not kwargs.has_key('full_output'):
                kwargs['full_output'] = False
            return self.process(*args, **kwargs)
        raise TypeError("Do not know what to do with that arguments.")


    #*************************************************************************
    # tests
    #*************************************************************************


    def __nonzero__(self):
        """
        Returns True if the finite state machine consists of at least
        one state.

        INPUT:

        Nothing.

        OUTPUT:

        True or False.

        TESTS::

            sage: FiniteStateMachine().__nonzero__()
            False
        """
        return len(self._states_) > 0


    def __eq__(left, right):
        """
        Returns True if the two finite state machines are equal, i.e.,
        if they have the same states and the same transitions.

        INPUT:

        - ``left`` -- a finite state machine.

        - ``right`` -- a finite state machine.

        OUTPUT:

        True or False.

        EXAMPLES::

            sage: F = FiniteStateMachine([('A', 'B', 1)])
            sage: F == FiniteStateMachine()
            False
        """
        if not is_FiniteStateMachine(right):
            raise TypeError('Only instances of FiniteStateMachine ' \
                'can be compared.')
        if len(left._states_) != len(right._states_):
            return False
        for state in left.iter_states():
            if state not in right._states_:
                return False
            left_transitions = state.transitions
            right_transitions = right.state(state).transitions
            if len(left_transitions) != len(right_transitions):
                return False
            for t in left_transitions:
                if t not in right_transitions:
                    return False
        return True


    def __ne__(left, right):
        """
        Tests for inequality, complement of :meth:`.__eq__`.

        INPUT:

        - ``left`` -- a finite state machine.

        - ``right`` -- a finite state machine.

        OUTPUT:

        True or False.

        EXAMPLES::

            sage: E = FiniteStateMachine([('A', 'B', 0)])
            sage: F = Automaton([('A', 'B', 0)])
            sage: G = Transducer([('A', 'B', 0, 1)])
            sage: E == F
            True
            sage: E == G
            False
        """
        return (not (left == right))


    def __contains__(self, item):
        """
        Returns true, if the finite state machine contains the
        state or transition item. Note that only the labels of the
        states and the input and output words are tested.

        INPUT:

        - ``item`` -- a state or a transition.

        OUTPUT:

        True or False.

        EXAMPLES::

            sage: from sage.combinat.finite_state_machine import FSMState, FSMTransition
            sage: F = FiniteStateMachine([('A', 'B', 0), ('B', 'A', 1)])
            sage: FSMState('A', is_initial=True) in F
            True
            sage: 'A' in F
            False
            sage: FSMTransition('A', 'B', 0) in F
            True
        """
        if is_FSMState(item):
            return self.has_state(item)
        if is_FSMTransition(item):
            return self.has_transition(item)
        return False


    def is_Markov_chain(self):
        """
        Checks whether ``self`` is a Markov chain where the transition
        probabilities are modeled as input labels.

        INPUT:

        Nothing.

        OUTPUT:

        True or False.

        ``on_duplicate_transition`` must be
        ``duplicate_transition_add_input`` and the sum of the input
        weights of the transitions leaving a state must add up to 1.

        EXAMPLES::

            sage: from sage.combinat.finite_state_machine import duplicate_transition_add_input
            sage: F = Transducer([[0, 0, 1/4, 0], [0, 1, 3/4, 1],
            ....:                 [1, 0, 1/2, 0], [1, 1, 1/2, 1]],
            ....:                on_duplicate_transition=duplicate_transition_add_input)
            sage: F.is_Markov_chain()
            True

        ``on_duplicate_transition`` must be ``duplicate_transition_add_input``::

            sage: F = Transducer([[0, 0, 1/4, 0], [0, 1, 3/4, 1],
            ....:                 [1, 0, 1/2, 0], [1, 1, 1/2, 1]])
            sage: F.is_Markov_chain()
            False

        Sum of input labels of the transitions leaving states must be 1::

            sage: F = Transducer([[0, 0, 1/4, 0], [0, 1, 3/4, 1],
            ....:                 [1, 0, 1/2, 0]],
            ....:                on_duplicate_transition=duplicate_transition_add_input)
            sage: F.is_Markov_chain()
            False
        """

        if self.on_duplicate_transition != duplicate_transition_add_input:
            return False

        return all((sum(t.word_in[0] for t in state.transitions) - 1).is_zero()
                   for state in self.states())


    #*************************************************************************
    # representations / LaTeX
    #*************************************************************************


    def _repr_(self):
        """
        Represents the finite state machine as "Finite state machine
        with n states" where n is the number of states.

        INPUT:

        Nothing.

        OUTPUT:

        A string.

        EXAMPLES::

            sage: FiniteStateMachine()._repr_()
            'Finite state machine with 0 states'
        """
        return "Finite state machine with %s states" % len(self._states_)

    default_format_letter = latex
    format_letter = default_format_letter


    def format_letter_negative(self, letter):
        r"""
        Format negative numbers as overlined numbers, everything
        else by standard LaTeX formatting.

        INPUT:

        ``letter`` -- anything.

        OUTPUT:

        Overlined absolute value if letter is a negative integer,
        :func:`latex(letter) <sage.misc.latex.latex>` otherwise.

        EXAMPLES::

            sage: A = Automaton([(0, 0, -1)])
            sage: map(A.format_letter_negative, [-1, 0, 1, 'a', None])
             ['\\overline{1}', 0, 1, \text{\texttt{a}}, \mbox{\rm None}]
            sage: A.latex_options(format_letter=A.format_letter_negative)
            sage: print(latex(A))
            \begin{tikzpicture}[auto, initial text=, >=latex]
            \node[state] (v0) at (3.000000, 0.000000) {$0$};
            \path[->] (v0) edge[loop above] node {$\overline{1}$} ();
            \end{tikzpicture}
        """
        if letter in ZZ and letter < 0:
            return r'\overline{%d}' % -letter
        else:
            return latex(letter)


    def format_transition_label_reversed(self, word):
        r"""
        Format words in transition labels in reversed order.

        INPUT:

        ``word`` -- list of letters.

        OUTPUT:

        String representation of ``word`` suitable to be typeset in
        mathematical mode, letters are written in reversed order.

        This is the reversed version of
        :meth:`.default_format_transition_label`.

        In digit expansions, digits are frequently processed from the
        least significant to the most significant position, but it is
        customary to write the least significant digit at the
        right-most position. Therefore, the labels have to be
        reversed.

        EXAMPLE::

            sage: T = Transducer([(0, 0, 0, [1, 2, 3])])
            sage: T.format_transition_label_reversed([1, 2, 3])
            '3 2 1'
            sage: T.latex_options(format_transition_label=T.format_transition_label_reversed)
            sage: print latex(T)
            \begin{tikzpicture}[auto, initial text=, >=latex]
            \node[state] (v0) at (3.000000, 0.000000) {$0$};
            \path[->] (v0) edge[loop above] node {$0\mid 3 2 1$} ();
            \end{tikzpicture}

        TEST:

            Check that #16357 is fixed::

                sage: T = Transducer()
                sage: T.format_transition_label_reversed([])
                '\\varepsilon'
        """
        return self.default_format_transition_label(reversed(word))


    def default_format_transition_label(self, word):
        r"""
        Default formatting of words in transition labels for LaTeX output.

        INPUT:

        ``word`` -- list of letters

        OUTPUT:

        String representation of ``word`` suitable to be typeset in
        mathematical mode.

        -   For a non-empty word: Concatenation of the letters, piped through
            ``self.format_letter`` and separated by blanks.
        -   For an empty word:
            ``sage.combinat.finite_state_machine.EmptyWordLaTeX``.

        There is also a variant :meth:`.format_transition_label_reversed`
        writing the words in reversed order.

        EXAMPLES:

        #.  Example of a non-empty word::

                sage: T = Transducer()
                sage: print T.default_format_transition_label(
                ....:    ['a', 'alpha', 'a_1', '0', 0, (0, 1)])
                \text{\texttt{a}} \text{\texttt{alpha}}
                \text{\texttt{a{\char`\_}1}} 0 0 \left(0, 1\right)

        #.  In the example above, ``'a'`` and ``'alpha'`` should perhaps
            be symbols::

                sage: var('a alpha a_1')
                (a, alpha, a_1)
                sage: print T.default_format_transition_label([a, alpha, a_1])
                a \alpha a_{1}

        #.  Example of an empty word::

                sage: print T.default_format_transition_label([])
                \varepsilon

            We can change this by setting
            ``sage.combinat.finite_state_machine.EmptyWordLaTeX``::

                sage: sage.combinat.finite_state_machine.EmptyWordLaTeX = ''
                sage: T.default_format_transition_label([])
                ''

            Finally, we restore the default value::

                sage: sage.combinat.finite_state_machine.EmptyWordLaTeX = r'\varepsilon'

        #.  This method is the default value for
            ``FiniteStateMachine.format_transition_label``. That can be changed to be
            any other function::

                sage: A = Automaton([(0, 1, 0)])
                sage: def custom_format_transition_label(word):
                ....:     return "t"
                sage: A.latex_options(format_transition_label=custom_format_transition_label)
                sage: print latex(A)
                \begin{tikzpicture}[auto, initial text=, >=latex]
                \node[state] (v0) at (3.000000, 0.000000) {$0$};
                \node[state] (v1) at (-3.000000, 0.000000) {$1$};
                \path[->] (v0) edge node[rotate=360.00, anchor=south] {$t$} (v1);
                \end{tikzpicture}

        TEST:

            Check that #16357 is fixed::

                sage: T = Transducer()
                sage: T.default_format_transition_label([])
                '\\varepsilon'
                sage: T.default_format_transition_label(iter([]))
                '\\varepsilon'
        """
        result = " ".join(imap(self.format_letter, word))
        if result:
            return result
        else:
            return EmptyWordLaTeX

    format_transition_label = default_format_transition_label


    def latex_options(self,
                      coordinates=None,
                      format_state_label=None,
                      format_letter=None,
                      format_transition_label=None,
                      loop_where=None,
                      initial_where=None,
                      accepting_style=None,
                      accepting_distance=None,
                      accepting_where=None,
                      accepting_show_empty=None):
        r"""
        Set options for LaTeX output via
        :func:`~sage.misc.latex.latex` and therefore
        :func:`~sage.misc.latex.view`.

        INPUT:

        - ``coordinates`` -- a dictionary or a function mapping labels
          of states to pairs interpreted as coordinates. If no
          coordinates are given, states a placed equidistantly on a
          circle of radius `3`. See also :meth:`.set_coordinates`.

        - ``format_state_label`` -- a function mapping labels of
          states to a string suitable for typesetting in LaTeX's
          mathematics mode. If not given, :func:`~sage.misc.latex.latex`
          is used.

        - ``format_letter`` -- a function mapping letters of the input
          and output alphabets to a string suitable for typesetting in
          LaTeX's mathematics mode. If not given,
          :meth:`.default_format_transition_label` uses
          :func:`~sage.misc.latex.latex`.

        - ``format_transition_label`` -- a function mapping words over
          the input and output alphabets to a string suitable for
          typesetting in LaTeX's mathematics mode. If not given,
          :meth:`.default_format_transition_label` is used.

        - ``loop_where`` -- a dictionary or a function mapping labels of
          initial states to one of ``'above'``, ``'left'``, ``'below'``,
          ``'right'``. If not given, ``'above'`` is used.

        - ``initial_where`` -- a dictionary or a function mapping
          labels of initial states to one of ``'above'``, ``'left'``,
          ``'below'``, ``'right'``. If not given, TikZ' default
          (currently ``'left'``) is used.

        - ``accepting_style`` -- one of ``'accepting by double'`` and
          ``'accepting by arrow'``. If not given, ``'accepting by
          double'`` is used unless there are non-empty final output
          words.

        - ``accepting_distance`` -- a string giving a LaTeX length
          used for the length of the arrow leading from a final state.
          If not given, TikZ' default (currently ``'3ex'``) is used
          unless there are non-empty final output words, in which case
          ``'7ex'`` is used.

        - ``accepting_where`` -- a dictionary or a function mapping
          labels of final states to one of ``'above'``, ``'left'``,
          ``'below'``, ``'right'``. If not given, TikZ' default
          (currently ``'right'``) is used. If the final state has a
          final output word, it is also possible to give an angle
          in degrees.

        - ``accepting_show_empty`` -- if ``True`` the arrow of an
          empty final output word is labeled as well. Note that this
          implicitly implies ``accepting_style='accepting by
          arrow'``. If not given, the default ``False`` is used.

        OUTPUT:

        Nothing.

        As TikZ (cf. the :wikipedia:`PGF/TikZ`) is used to typeset
        the graphics, the syntax is oriented on TikZ' syntax.

        This is a convenience function collecting all options for
        LaTeX output. All of its functionality can also be achieved by
        directly setting the attributes

        - ``coordinates``, ``format_label``, ``loop_where``,
          ``initial_where``, and ``accepting_where`` of
          :class:`FSMState` (here, ``format_label`` is a callable
          without arguments, everything else is a specific value);

        - ``format_label`` of :class:`FSMTransition` (``format_label``
          is a callable without arguments);

        - ``format_state_label``, ``format_letter``,
          ``format_transition_label``, ``accepting_style``,
          ``accepting_distance``, and ``accepting_show_empty``
          of :class:`FiniteStateMachine`.

        This function, however, also (somewhat) checks its input and
        serves to collect documentation on all these options.

        The function can be called several times, only those arguments
        which are not ``None`` are taken into account. By the same
        means, it can be combined with directly setting some
        attributes as outlined above.

        See also :func:`setup_latex_preamble` or the example below on
        how to setup the LaTeX environment.

        EXAMPLES:

        See also the section on :ref:`finite_state_machine_LaTeX_output`
        in the introductory examples of this module.

        ::

            sage: from sage.combinat.finite_state_machine import setup_latex_preamble
            sage: setup_latex_preamble()
            sage: latex.mathjax_avoid_list('tikzpicture')
            sage: T = Transducer(initial_states=['I'],
            ....:     final_states=[0, 3])
            sage: for j in srange(4):
            ....:     T.add_transition('I', j, 0, [0, j])
            ....:     T.add_transition(j, 'I', 0, [0, -j])
            ....:     T.add_transition(j, j, 0, 0)
            Transition from 'I' to 0: 0|0,0
            Transition from 0 to 'I': 0|0,0
            Transition from 0 to 0: 0|0
            Transition from 'I' to 1: 0|0,1
            Transition from 1 to 'I': 0|0,-1
            Transition from 1 to 1: 0|0
            Transition from 'I' to 2: 0|0,2
            Transition from 2 to 'I': 0|0,-2
            Transition from 2 to 2: 0|0
            Transition from 'I' to 3: 0|0,3
            Transition from 3 to 'I': 0|0,-3
            Transition from 3 to 3: 0|0
            sage: T.add_transition('I', 'I', 0, 0)
            Transition from 'I' to 'I': 0|0
            sage: T.state(3).final_word_out = [0, 0]
            sage: T.latex_options(
            ....:     coordinates={'I': (0, 0),
            ....:                  0: (-6, 3),
            ....:                  1: (-2, 3),
            ....:                  2: (2, 3),
            ....:                  3: (6, 3)},
            ....:     format_state_label=lambda x: r'\mathbf{%s}' % x,
            ....:     format_letter=lambda x: r'w_{%s}' % x,
            ....:     format_transition_label=lambda x:
            ....:         r"{\scriptstyle %s}" % T.default_format_transition_label(x),
            ....:     loop_where={'I': 'below', 0: 'left', 1: 'above',
            ....:                 2: 'right', 3:'below'},
            ....:     initial_where=lambda x: 'above',
            ....:     accepting_style='accepting by double',
            ....:     accepting_distance='10ex',
            ....:     accepting_where={0: 'left', 3: 45}
            ....:     )
            sage: T.state('I').format_label=lambda: r'\mathcal{I}'
            sage: latex(T)
            \begin{tikzpicture}[auto, initial text=, >=latex]
            \node[state, initial, initial where=above] (v0) at (0.000000, 0.000000) {$\mathcal{I}$};
            \node[state, accepting, accepting where=left] (v1) at (-6.000000, 3.000000) {$\mathbf{0}$};
            \node[state, accepting, accepting where=45] (v2) at (6.000000, 3.000000) {$\mathbf{3}$};
            \path[->] (v2.45.00) edge node[rotate=45.00, anchor=south] {$\$ \mid {\scriptstyle w_{0} w_{0}}$} ++(45.00:10ex);
            \node[state] (v3) at (-2.000000, 3.000000) {$\mathbf{1}$};
            \node[state] (v4) at (2.000000, 3.000000) {$\mathbf{2}$};
            \path[->] (v0.61.31) edge node[rotate=56.31, anchor=south] {${\scriptstyle w_{0}}\mid {\scriptstyle w_{0} w_{2}}$} (v4.231.31);
            \path[->] (v1.-21.57) edge node[rotate=-26.57, anchor=south] {${\scriptstyle w_{0}}\mid {\scriptstyle w_{0} w_{0}}$} (v0.148.43);
            \path[->] (v0.31.57) edge node[rotate=26.57, anchor=south] {${\scriptstyle w_{0}}\mid {\scriptstyle w_{0} w_{3}}$} (v2.201.57);
            \path[->] (v2) edge[loop below] node {${\scriptstyle w_{0}}\mid {\scriptstyle w_{0}}$} ();
            \path[->] (v2.-148.43) edge node[rotate=26.57, anchor=north] {${\scriptstyle w_{0}}\mid {\scriptstyle w_{0} w_{-3}}$} (v0.21.57);
            \path[->] (v3.-51.31) edge node[rotate=-56.31, anchor=south] {${\scriptstyle w_{0}}\mid {\scriptstyle w_{0} w_{-1}}$} (v0.118.69);
            \path[->] (v4) edge[loop right] node[rotate=90, anchor=north] {${\scriptstyle w_{0}}\mid {\scriptstyle w_{0}}$} ();
            \path[->] (v3) edge[loop above] node {${\scriptstyle w_{0}}\mid {\scriptstyle w_{0}}$} ();
            \path[->] (v1) edge[loop left] node[rotate=90, anchor=south] {${\scriptstyle w_{0}}\mid {\scriptstyle w_{0}}$} ();
            \path[->] (v4.-118.69) edge node[rotate=56.31, anchor=north] {${\scriptstyle w_{0}}\mid {\scriptstyle w_{0} w_{-2}}$} (v0.51.31);
            \path[->] (v0.158.43) edge node[rotate=333.43, anchor=north] {${\scriptstyle w_{0}}\mid {\scriptstyle w_{0} w_{0}}$} (v1.328.43);
            \path[->] (v0.128.69) edge node[rotate=303.69, anchor=north] {${\scriptstyle w_{0}}\mid {\scriptstyle w_{0} w_{1}}$} (v3.298.69);
            \path[->] (v0) edge[loop below] node {${\scriptstyle w_{0}}\mid {\scriptstyle w_{0}}$} ();
            \end{tikzpicture}
            sage: view(T) # not tested

        To actually see this, use the live documentation in the Sage notebook
        and execute the cells.

        By changing some of the options, we get the following output::

            sage: T.latex_options(
            ....:     format_transition_label=T.default_format_transition_label,
            ....:     accepting_style='accepting by arrow',
            ....:     accepting_show_empty=True
            ....:     )
            sage: latex(T)
            \begin{tikzpicture}[auto, initial text=, >=latex, accepting text=, accepting/.style=accepting by arrow, accepting distance=10ex]
            \node[state, initial, initial where=above] (v0) at (0.000000, 0.000000) {$\mathcal{I}$};
            \node[state] (v1) at (-6.000000, 3.000000) {$\mathbf{0}$};
            \path[->] (v1.180.00) edge node[rotate=360.00, anchor=south] {$\$ \mid \varepsilon$} ++(180.00:10ex);
            \node[state] (v2) at (6.000000, 3.000000) {$\mathbf{3}$};
            \path[->] (v2.45.00) edge node[rotate=45.00, anchor=south] {$\$ \mid w_{0} w_{0}$} ++(45.00:10ex);
            \node[state] (v3) at (-2.000000, 3.000000) {$\mathbf{1}$};
            \node[state] (v4) at (2.000000, 3.000000) {$\mathbf{2}$};
            \path[->] (v0.61.31) edge node[rotate=56.31, anchor=south] {$w_{0}\mid w_{0} w_{2}$} (v4.231.31);
            \path[->] (v1.-21.57) edge node[rotate=-26.57, anchor=south] {$w_{0}\mid w_{0} w_{0}$} (v0.148.43);
            \path[->] (v0.31.57) edge node[rotate=26.57, anchor=south] {$w_{0}\mid w_{0} w_{3}$} (v2.201.57);
            \path[->] (v2) edge[loop below] node {$w_{0}\mid w_{0}$} ();
            \path[->] (v2.-148.43) edge node[rotate=26.57, anchor=north] {$w_{0}\mid w_{0} w_{-3}$} (v0.21.57);
            \path[->] (v3.-51.31) edge node[rotate=-56.31, anchor=south] {$w_{0}\mid w_{0} w_{-1}$} (v0.118.69);
            \path[->] (v4) edge[loop right] node[rotate=90, anchor=north] {$w_{0}\mid w_{0}$} ();
            \path[->] (v3) edge[loop above] node {$w_{0}\mid w_{0}$} ();
            \path[->] (v1) edge[loop left] node[rotate=90, anchor=south] {$w_{0}\mid w_{0}$} ();
            \path[->] (v4.-118.69) edge node[rotate=56.31, anchor=north] {$w_{0}\mid w_{0} w_{-2}$} (v0.51.31);
            \path[->] (v0.158.43) edge node[rotate=333.43, anchor=north] {$w_{0}\mid w_{0} w_{0}$} (v1.328.43);
            \path[->] (v0.128.69) edge node[rotate=303.69, anchor=north] {$w_{0}\mid w_{0} w_{1}$} (v3.298.69);
            \path[->] (v0) edge[loop below] node {$w_{0}\mid w_{0}$} ();
            \end{tikzpicture}
            sage: view(T) # not tested

        TESTS::

            sage: T.latex_options(format_state_label='Nothing')
            Traceback (most recent call last):
            ...
            TypeError: format_state_label must be callable.
            sage: T.latex_options(format_letter='')
            Traceback (most recent call last):
            ...
            TypeError: format_letter must be callable.
            sage: T.latex_options(format_transition_label='')
            Traceback (most recent call last):
            ...
            TypeError: format_transition_label must be callable.
            sage: T.latex_options(loop_where=37)
            Traceback (most recent call last):
            ...
            TypeError: loop_where must be a callable or a
            dictionary.
            sage: T.latex_options(loop_where=lambda x: 'top')
            Traceback (most recent call last):
            ...
            ValueError: loop_where for I must be in ['below',
            'right', 'above', 'left'].
            sage: T.latex_options(initial_where=90)
            Traceback (most recent call last):
            ...
            TypeError: initial_where must be a callable or a
            dictionary.
            sage: T.latex_options(initial_where=lambda x: 'top')
            Traceback (most recent call last):
            ...
            ValueError: initial_where for I must be in ['below',
            'right', 'above', 'left'].
            sage: T.latex_options(accepting_style='fancy')
            Traceback (most recent call last):
            ...
            ValueError: accepting_style must be in ['accepting by
            double', 'accepting by arrow'].
            sage: T.latex_options(accepting_where=90)
            Traceback (most recent call last):
            ...
            TypeError: accepting_where must be a callable or a
            dictionary.
            sage: T.latex_options(accepting_where=lambda x: 'top')
            Traceback (most recent call last):
            ...
            ValueError: accepting_where for 0 must be in ['below',
            'right', 'above', 'left'].
            sage: T.latex_options(accepting_where={0: 'above', 3: 'top'})
            Traceback (most recent call last):
            ...
            ValueError: accepting_where for 3 must be a real number or
            be in ['below', 'right', 'above', 'left'].
        """
        if coordinates is not None:
            self.set_coordinates(coordinates)

        if format_state_label is not None:
            if not hasattr(format_state_label, '__call__'):
                raise TypeError('format_state_label must be callable.')
            self.format_state_label = format_state_label

        if format_letter is not None:
            if not hasattr(format_letter, '__call__'):
                raise TypeError('format_letter must be callable.')
            self.format_letter = format_letter

        if format_transition_label is not None:
            if not hasattr(format_transition_label, '__call__'):
                raise TypeError('format_transition_label must be callable.')
            self.format_transition_label = format_transition_label

        if loop_where is not None:
            permissible = list(tikz_automata_where.iterkeys())
            for state in self.states():
                if hasattr(loop_where, '__call__'):
                    where = loop_where(state.label())
                else:
                    try:
                        where = loop_where[state.label()]
                    except TypeError:
                        raise TypeError("loop_where must be a "
                                        "callable or a dictionary.")
                    except KeyError:
                        continue
                if where in permissible:
                    state.loop_where = where
                else:
                    raise ValueError('loop_where for %s must be in %s.' %
                                     (state.label(), permissible))

        if initial_where is not None:
            permissible = list(tikz_automata_where.iterkeys())
            for state in self.iter_initial_states():
                if hasattr(initial_where, '__call__'):
                    where = initial_where(state.label())
                else:
                    try:
                        where = initial_where[state.label()]
                    except TypeError:
                        raise TypeError("initial_where must be a "
                                        "callable or a dictionary.")
                    except KeyError:
                        continue
                if where in permissible:
                    state.initial_where = where
                else:
                    raise ValueError('initial_where for %s must be in %s.' %
                                     (state.label(), permissible))

        if accepting_style is not None:
            permissible = ['accepting by double',
                           'accepting by arrow']
            if accepting_style in permissible:
                self.accepting_style = accepting_style
            else:
                raise ValueError('accepting_style must be in %s.' %
                    permissible)

        if accepting_distance is not None:
            self.accepting_distance = accepting_distance

        if accepting_where is not None:
            permissible = list(tikz_automata_where.iterkeys())
            for state in self.iter_final_states():
                if hasattr(accepting_where, '__call__'):
                    where = accepting_where(state.label())
                else:
                    try:
                        where = accepting_where[state.label()]
                    except TypeError:
                        raise TypeError("accepting_where must be a "
                                        "callable or a dictionary.")
                    except KeyError:
                        continue
                if where in permissible:
                    state.accepting_where = where
                elif hasattr(state, 'final_word_out') \
                        and state.final_word_out:
                    if where in RR:
                        state.accepting_where = where
                    else:
                        raise ValueError('accepting_where for %s must '
                                         'be a real number or be in %s.' %
                                         (state.label(), permissible))

                else:
                    raise ValueError('accepting_where for %s must be in %s.' %
                                     (state.label(), permissible))

        if accepting_show_empty is not None:
            self.accepting_show_empty = accepting_show_empty


    def _latex_(self):
        r"""
        Returns a LaTeX code for the graph of the finite state machine.

        INPUT:

        Nothing.

        OUTPUT:

        A string.

        EXAMPLES::

            sage: F = FiniteStateMachine([('A', 'B', 1, 2)],
            ....:                        initial_states=['A'],
            ....:                        final_states=['B'])
            sage: F.state('A').initial_where='below'
            sage: print latex(F)  # indirect doctest
            \begin{tikzpicture}[auto, initial text=, >=latex]
            \node[state, initial, initial where=below] (v0) at (3.000000, 0.000000) {$\text{\texttt{A}}$};
            \node[state, accepting] (v1) at (-3.000000, 0.000000) {$\text{\texttt{B}}$};
            \path[->] (v0) edge node[rotate=360.00, anchor=south] {$ $} (v1);
            \end{tikzpicture}
        """
        def label_rotation(angle, both_directions):
            """
            Given an angle of a transition, compute the TikZ string to
            rotate the label.
            """
            angle_label = angle
            anchor_label = "south"
            if angle > 90 or angle <= -90:
                angle_label = angle + 180
                if both_directions:
                    # if transitions in both directions, the transition to the
                    # left has its label below the transition, otherwise above
                    anchor_label = "north"
            return "rotate=%.2f, anchor=%s" % (angle_label, anchor_label)

        options = ["auto", "initial text=", ">=latex"]

        nonempty_final_word_out = False
        for state in self.iter_final_states():
            if state.final_word_out:
                nonempty_final_word_out = True
                break

        if hasattr(self, "accepting_style"):
            accepting_style = self.accepting_style
        elif nonempty_final_word_out:
            accepting_style = "accepting by arrow"
        else:
            accepting_style = "accepting by double"

        if accepting_style == "accepting by arrow":
            options.append("accepting text=")
            options.append("accepting/.style=%s" % accepting_style)

        if hasattr(self, "accepting_distance"):
            accepting_distance = self.accepting_distance
        elif nonempty_final_word_out:
            accepting_distance = "7ex"
        else:
            accepting_distance = None
        if accepting_style == "accepting by arrow" and accepting_distance:
            options.append("accepting distance=%s"
                           % accepting_distance)

        if hasattr(self, "accepting_show_empty"):
            accepting_show_empty = self.accepting_show_empty
        else:
            accepting_show_empty = False

        result = "\\begin{tikzpicture}[%s]\n" % ", ".join(options)
        j = 0;
        for vertex in self.iter_states():
            if not hasattr(vertex, "coordinates"):
                vertex.coordinates = (3*cos(2*pi*j/len(self.states())),
                                      3*sin(2*pi*j/len(self.states())))
            options = ""
            if vertex.is_final:
                if not (vertex.final_word_out
                        and accepting_style == "accepting by arrow") \
                        and not accepting_show_empty:
                    # otherwise, we draw a custom made accepting path
                    # with label below
                    options += ", accepting"
                    if hasattr(vertex, "accepting_where"):
                        options += ", accepting where=%s" % (
                            vertex.accepting_where,)
            if vertex.is_initial:
                options += ", initial"
            if hasattr(vertex, "initial_where"):
                options += ", initial where=%s" % vertex.initial_where
            if hasattr(vertex, "format_label"):
                label = vertex.format_label()
            elif hasattr(self, "format_state_label"):
                label = self.format_state_label(vertex)
            else:
                label = latex(vertex.label())
            result += "\\node[state%s] (v%d) at (%f, %f) {$%s$};\n" % (
                options, j, vertex.coordinates[0],
                vertex.coordinates[1], label)
            vertex._number_ = j
            if vertex.is_final and (vertex.final_word_out or accepting_show_empty):
                angle = 0
                if hasattr(vertex, "accepting_where"):
                    angle = tikz_automata_where.get(vertex.accepting_where,
                                                    vertex.accepting_where)
                result += "\\path[->] (v%d.%.2f) edge node[%s] {$%s \mid %s$} ++(%.2f:%s);\n" % (
                    j, angle,
                    label_rotation(angle, False),
                    EndOfWordLaTeX,
                    self.format_transition_label(vertex.final_word_out),
                    angle, accepting_distance)

            j += 1
        adjacent = {}
        for source in self.iter_states():
            for target in self.iter_states():
                transitions = [transition for transition in source.transitions if transition.to_state == target]
                adjacent[source, target] = transitions

        for ((source, target), transitions) in adjacent.iteritems():
            if len(transitions) > 0:
                labels = []
                for transition in transitions:
                    if hasattr(transition, "format_label"):
                        labels.append(transition.format_label())
                    else:
                        labels.append(self._latex_transition_label_(
                                transition, self.format_transition_label))
                label = ", ".join(labels)
                if source != target:
                    angle = atan2(
                        target.coordinates[1] - source.coordinates[1],
                        target.coordinates[0] - source.coordinates[0]) * 180/pi
                    both_directions = len(adjacent[target, source]) > 0
                    if both_directions:
                        angle_source = ".%.2f" % ((angle + 5).n(),)
                        angle_target = ".%.2f" % ((angle + 175).n(),)
                    else:
                        angle_source = ""
                        angle_target = ""
                    result += "\\path[->] (v%d%s) edge node[%s] {$%s$} (v%d%s);\n" % (
                        source._number_, angle_source,
                        label_rotation(angle, both_directions),
                        label,
                        target._number_, angle_target)
                else:
                    loop_where = "above"
                    if hasattr(source, "loop_where"):
                        loop_where = source.loop_where
                    rotation = {'left': '[rotate=90, anchor=south]',
                                'right': '[rotate=90, anchor=north]'}
                    result += "\\path[->] (v%d) edge[loop %s] node%s {$%s$} ();\n" % (
                        source._number_,
                        loop_where, rotation.get(loop_where, ''),
                        label)

        result += "\\end{tikzpicture}"
        return result


    def _latex_transition_label_(self, transition, format_function=latex):
        r"""
        Returns the proper transition label.

        INPUT:

        - ``transition`` - a transition

        - ``format_function`` - a function formatting the labels

        OUTPUT:

        A string.

        TESTS::

            sage: F = FiniteStateMachine([('A', 'B', 0, 1)])
            sage: t = F.transitions()[0]
            sage: F._latex_transition_label_(t)
            ' '
        """
        return ' '

    def set_coordinates(self, coordinates, default=True):
        """
        Set coordinates of the states for the LaTeX representation by
        a dictionary or a function mapping labels to coordinates.

        INPUT:

        - ``coordinates`` -- a dictionary or a function mapping labels
          of states to pairs interpreted as coordinates.

        - ``default`` -- If ``True``, then states not given by
          ``coordinates`` get a default position on a circle of
          radius 3.

        OUTPUT:

        Nothing.

        EXAMPLES::

            sage: F = Automaton([[0, 1, 1], [1, 2, 2], [2, 0, 0]])
            sage: F.set_coordinates({0: (0, 0), 1: (2, 0), 2: (1, 1)})
            sage: F.state(0).coordinates
            (0, 0)

        We can also use a function to determine the coordinates::

            sage: F = Automaton([[0, 1, 1], [1, 2, 2], [2, 0, 0]])
            sage: F.set_coordinates(lambda l: (l, 3/(l+1)))
            sage: F.state(2).coordinates
            (2, 1)
        """
        states_without_coordinates = []
        for state in self.iter_states():
            try:
                state.coordinates = coordinates[state.label()]
                continue
            except (KeyError, TypeError):
                pass

            try:
                state.coordinates = coordinates(state.label())
                continue
            except TypeError:
                pass

            states_without_coordinates.append(state)

        if default:
            n = len(states_without_coordinates)
            for j, state in enumerate(states_without_coordinates):
                state.coordinates = (3*cos(2*pi*j/n),
                                     3*sin(2*pi*j/n))


    #*************************************************************************
    # other
    #*************************************************************************


    def _matrix_(self, R=None):
        """
        Returns the adjacency matrix of the finite state machine.
        See :meth:`.adjacency_matrix` for more information.

        EXAMPLES::

            sage: B = FiniteStateMachine({0: {0: (0, 0), 'a': (1, 0)},
            ....:                         'a': {2: (0, 0), 3: (1, 0)},
            ....:                         2:{0:(1, 1), 4:(0, 0)},
            ....:                         3:{'a':(0, 1), 2:(1, 1)},
            ....:                         4:{4:(1, 1), 3:(0, 1)}},
            ....:                        initial_states=[0])
            sage: B._matrix_()
            [1 1 0 0 0]
            [0 0 1 1 0]
            [x 0 0 0 1]
            [0 x x 0 0]
            [0 0 0 x x]
        """
        return self.adjacency_matrix()


    def adjacency_matrix(self, input=None,
                         entry=None):
        """
        Returns the adjacency matrix of the underlying graph.

        INPUT:

        - ``input`` -- Only transitions with input label ``input`` are
          respected.

        - ``entry`` -- The function ``entry`` takes a transition and the
          return value is written in the matrix as the entry
          ``(transition.from_state, transition.to_state)``. The default
          value (``None``) of entry takes the variable ``x`` to the
          power of the sum of the output word of the transition.

        OUTPUT:

        A matrix.

        If any label of a state is not an integer, the finite state
        machine is relabeled at the beginning.  If there are more than
        one transitions between two states, then the different return
        values of ``entry`` are added up.

        EXAMPLES::

            sage: B = FiniteStateMachine({0:{0:(0, 0), 'a':(1, 0)},
            ....:                         'a':{2:(0, 0), 3:(1, 0)},
            ....:                         2:{0:(1, 1), 4:(0, 0)},
            ....:                         3:{'a':(0, 1), 2:(1, 1)},
            ....:                         4:{4:(1, 1), 3:(0, 1)}},
            ....:                        initial_states=[0])
            sage: B.adjacency_matrix()
            [1 1 0 0 0]
            [0 0 1 1 0]
            [x 0 0 0 1]
            [0 x x 0 0]
            [0 0 0 x x]

        This is equivalent to::

            sage: matrix(B)
            [1 1 0 0 0]
            [0 0 1 1 0]
            [x 0 0 0 1]
            [0 x x 0 0]
            [0 0 0 x x]

        It is also possible to use other entries in the adjacency matrix::

            sage: B.adjacency_matrix(entry=(lambda transition: 1))
            [1 1 0 0 0]
            [0 0 1 1 0]
            [1 0 0 0 1]
            [0 1 1 0 0]
            [0 0 0 1 1]
            sage: B.adjacency_matrix(1, entry=(lambda transition:
            ....:     exp(I*transition.word_out[0]*var('t'))))
            [      0       1       0       0       0]
            [      0       0       0       1       0]
            [e^(I*t)       0       0       0       0]
            [      0       0 e^(I*t)       0       0]
            [      0       0       0       0 e^(I*t)]
            sage: a = Automaton([(0, 1, 0),
            ....:                (1, 2, 0),
            ....:                (2, 0, 1),
            ....:                (2, 1, 0)],
            ....:               initial_states=[0],
            ....:               final_states=[0])
            sage: a.adjacency_matrix()
            [0 1 0]
            [0 0 1]
            [1 1 0]

        """
        def default_function(transitions):
            var('x')
            return x**sum(transition.word_out)

        if entry is None:
            entry = default_function

        relabeledFSM = self
        l = len(relabeledFSM.states())
        for state in self.iter_states():
            if state.label() not in ZZ or state.label() >= l \
                     or state.label() < 0:
                relabeledFSM = self.relabeled()
                break
        dictionary = {}
        for transition in relabeledFSM.iter_transitions():
            if input is None or transition.word_in == [input]:
                if (transition.from_state.label(),
                    transition.to_state.label()) in dictionary:
                    dictionary[(transition.from_state.label(),
                                transition.to_state.label())] \
                                += entry(transition)
                else:
                    dictionary[(transition.from_state.label(),
                                transition.to_state.label())] \
                                = entry(transition)
        return matrix(len(relabeledFSM.states()), dictionary)


    def determine_alphabets(self, reset=True):
        """
        Determines the input and output alphabet according to the
        transitions in self.

        INPUT:

        - ``reset`` -- If reset is ``True``, then the existing input
          and output alphabets are erased, otherwise new letters are
          appended to the existing alphabets.

        OUTPUT:

        Nothing.

        After this operation the input alphabet and the output
        alphabet of self are a list of letters.

        EXAMPLES::

            sage: T = Transducer([(1, 1, 1, 0), (1, 2, 2, 1),
            ....:                 (2, 2, 1, 1), (2, 2, 0, 0)],
            ....:                determine_alphabets=False)
            sage: (T.input_alphabet, T.output_alphabet)
            (None, None)
            sage: T.determine_alphabets()
            sage: (T.input_alphabet, T.output_alphabet)
            ([0, 1, 2], [0, 1])
       """
        if reset:
            ain = set()
            aout = set()
        else:
            ain = set(self.input_alphabet)
            aout = set(self.output_alphabet)

        for t in self.iter_transitions():
            for letter in t.word_in:
                ain.add(letter)
            for letter in t.word_out:
                aout.add(letter)
        self.input_alphabet = list(ain)
        self.output_alphabet = list(aout)


    #*************************************************************************
    # get states and transitions
    #*************************************************************************


    def states(self):
        """
        Returns the states of the finite state machine.

        INPUT:

        Nothing.

        OUTPUT:

        The states of the finite state machine as list.

        EXAMPLES::

            sage: FSM = Automaton([('1', '2', 1), ('2', '2', 0)])
            sage: FSM.states()
            ['1', '2']

        """
        return copy(self._states_)


    def iter_states(self):
        """
        Returns an iterator of the states.

        INPUT:

        Nothing.

        OUTPUT:

        An iterator of the states of the finite state machine.

        EXAMPLES::

            sage: FSM = Automaton([('1', '2', 1), ('2', '2', 0)])
            sage: [s.label() for s in FSM.iter_states()]
            ['1', '2']
        """
        return iter(self._states_)


    def transitions(self, from_state=None):
        """
        Returns a list of all transitions.

        INPUT:

        - ``from_state`` -- (default: ``None``) If ``from_state`` is
          given, then a list of transitions starting there is given.

        OUTPUT:

        A list of all transitions.

        EXAMPLES::

            sage: FSM = Automaton([('1', '2', 1), ('2', '2', 0)])
            sage: FSM.transitions()
            [Transition from '1' to '2': 1|-,
             Transition from '2' to '2': 0|-]
        """
        return list(self.iter_transitions(from_state))


    def iter_transitions(self, from_state=None):
        """
        Returns an iterator of all transitions.

        INPUT:

        - ``from_state`` -- (default: ``None``) If ``from_state`` is
          given, then a list of transitions starting there is given.

        OUTPUT:

        An iterator of all transitions.

        EXAMPLES::

            sage: FSM = Automaton([('1', '2', 1), ('2', '2', 0)])
            sage: [(t.from_state.label(), t.to_state.label())
            ....:     for t in FSM.iter_transitions('1')]
            [('1', '2')]
            sage: [(t.from_state.label(), t.to_state.label())
            ....:     for t in FSM.iter_transitions('2')]
            [('2', '2')]
            sage: [(t.from_state.label(), t.to_state.label())
            ....:     for t in FSM.iter_transitions()]
            [('1', '2'), ('2', '2')]
        """
        if from_state is None:
            return self._iter_transitions_all_()
        else:
            return iter(self.state(from_state).transitions)


    def _iter_transitions_all_(self):
        """
        Returns an iterator over all transitions.

        INPUT:

        Nothing.

        OUTPUT:

        An iterator over all transitions.

        EXAMPLES::

            sage: FSM = Automaton([('1', '2', 1), ('2', '2', 0)])
            sage: [(t.from_state.label(), t.to_state.label())
            ....:     for t in FSM._iter_transitions_all_()]
            [('1', '2'), ('2', '2')]
        """
        for state in self.iter_states():
            for t in state.transitions:
                yield t


    def initial_states(self):
        """
        Returns a list of all initial states.

        INPUT:

        Nothing.

        OUTPUT:

        A list of all initial states.

        EXAMPLES::

            sage: from sage.combinat.finite_state_machine import FSMState
            sage: A = FSMState('A', is_initial=True)
            sage: B = FSMState('B')
            sage: F = FiniteStateMachine([(A, B, 1, 0)])
            sage: F.initial_states()
            ['A']
        """
        return list(self.iter_initial_states())


    def iter_initial_states(self):
        """
        Returns an iterator of the initial states.

        INPUT:

        Nothing.

        OUTPUT:

        An iterator over all initial states.

        EXAMPLES::

            sage: from sage.combinat.finite_state_machine import FSMState
            sage: A = FSMState('A', is_initial=True)
            sage: B = FSMState('B')
            sage: F = FiniteStateMachine([(A, B, 1, 0)])
            sage: [s.label() for s in F.iter_initial_states()]
            ['A']
        """
        return itertools.ifilter(lambda s:s.is_initial, self.iter_states())


    def final_states(self):
        """
        Returns a list of all final states.

        INPUT:

        Nothing.

        OUTPUT:

        A list of all final states.

        EXAMPLES::

            sage: from sage.combinat.finite_state_machine import FSMState
            sage: A = FSMState('A', is_final=True)
            sage: B = FSMState('B', is_initial=True)
            sage: C = FSMState('C', is_final=True)
            sage: F = FiniteStateMachine([(A, B), (A, C)])
            sage: F.final_states()
            ['A', 'C']
        """
        return list(self.iter_final_states())


    def iter_final_states(self):
        """
        Returns an iterator of the final states.

        INPUT:

        Nothing.

        OUTPUT:

        An iterator over all initial states.

        EXAMPLES::

            sage: from sage.combinat.finite_state_machine import FSMState
            sage: A = FSMState('A', is_final=True)
            sage: B = FSMState('B', is_initial=True)
            sage: C = FSMState('C', is_final=True)
            sage: F = FiniteStateMachine([(A, B), (A, C)])
            sage: [s.label() for s in F.iter_final_states()]
            ['A', 'C']
        """
        return itertools.ifilter(lambda s:s.is_final, self.iter_states())


    def state(self, state):
        """
        Returns the state of the finite state machine.

        INPUT:

        - ``state`` -- If ``state`` is not an instance of
          :class:`FSMState`, then it is assumed that it is the label
          of a state.

        OUTPUT:

        Returns the state of the finite state machine corresponding to
        ``state``.

        If no state is found, then a ``LookupError`` is thrown.

        EXAMPLES::

            sage: from sage.combinat.finite_state_machine import FSMState
            sage: A = FSMState('A')
            sage: FSM = FiniteStateMachine([(A, 'B'), ('C', A)])
            sage: FSM.state('A') == A
            True
            sage: FSM.state('xyz')
            Traceback (most recent call last):
            ...
            LookupError: No state with label xyz found.
        """
        def what(s, switch):
            if switch:
                return s.label()
            else:
                return s
        switch = is_FSMState(state)

        try:
            return self._states_dict_[what(state, switch)]
        except AttributeError:
            for s in self.iter_states():
                if what(s, not switch) == state:
                    return s
        except KeyError:
            pass
        raise LookupError("No state with label %s found." % (what(state, switch),))


    def transition(self, transition):
        """
        Returns the transition of the finite state machine.

        INPUT:

        - ``transition`` -- If ``transition`` is not an instance of
          :class:`FSMTransition`, then it is assumed that it is a
          tuple ``(from_state, to_state, word_in, word_out)``.

        OUTPUT:

        Returns the transition of the finite state machine
        corresponding to ``transition``.

        If no transition is found, then a ``LookupError`` is thrown.

        EXAMPLES::

            sage: from sage.combinat.finite_state_machine import FSMTransition
            sage: t = FSMTransition('A', 'B', 0)
            sage: F = FiniteStateMachine([t])
            sage: F.transition(('A', 'B', 0))
            Transition from 'A' to 'B': 0|-
            sage: id(t) == id(F.transition(('A', 'B', 0)))
            True
       """
        if not is_FSMTransition(transition):
            transition = FSMTransition(*transition)
        for s in self.iter_transitions(transition.from_state):
            if s == transition:
                return s
        raise LookupError("No transition found.")


    #*************************************************************************
    # properties (state and transitions)
    #*************************************************************************


    def has_state(self, state):
        """
        Returns whether ``state`` is one of the states of the finite
        state machine.

        INPUT:

        - ``state`` can be a :class:`FSMState` or a label of a state.

        OUTPUT:

        True or False.

        EXAMPLES::

            sage: FiniteStateMachine().has_state('A')
            False
        """
        try:
            self.state(state)
            return True
        except LookupError:
            return False


    def has_transition(self, transition):
        """
        Returns whether ``transition`` is one of the transitions of
        the finite state machine.

        INPUT:

        - ``transition`` has to be a :class:`FSMTransition`.

        OUTPUT:

        True or False.

        EXAMPLES::

            sage: from sage.combinat.finite_state_machine import FSMTransition
            sage: t = FSMTransition('A', 'A', 0, 1)
            sage: FiniteStateMachine().has_transition(t)
            False
            sage: FiniteStateMachine().has_transition(('A', 'A', 0, 1))
            Traceback (most recent call last):
            ...
            TypeError: Transition is not an instance of FSMTransition.
        """
        if is_FSMTransition(transition):
            return transition in self.iter_transitions()
        raise TypeError("Transition is not an instance of FSMTransition.")


    def has_initial_state(self, state):
        """
        Returns whether ``state`` is one of the initial states of the
        finite state machine.

        INPUT:

        - ``state`` can be a :class:`FSMState` or a label.

        OUTPUT:

        True or False.

        EXAMPLES::

            sage: F = FiniteStateMachine([('A', 'A')], initial_states=['A'])
            sage: F.has_initial_state('A')
            True
        """
        try:
            return self.state(state).is_initial
        except LookupError:
            return False


    def has_initial_states(self):
        """
        Returns whether the finite state machine has an initial state.

        INPUT:

        Nothing.

        OUTPUT:

        True or False.

        EXAMPLES::

            sage: FiniteStateMachine().has_initial_states()
            False
        """
        return len(self.initial_states()) > 0


    def has_final_state(self, state):
        """
        Returns whether ``state`` is one of the final states of the
        finite state machine.

        INPUT:

        - ``state`` can be a :class:`FSMState` or a label.

        OUTPUT:

        True or False.

        EXAMPLES::

            sage: FiniteStateMachine(final_states=['A']).has_final_state('A')
            True
        """
        try:
            return self.state(state).is_final
        except LookupError:
            return False


    def has_final_states(self):
        """
        Returns whether the finite state machine has a final state.

        INPUT:

        Nothing.

        OUTPUT:

        True or False.

        EXAMPLES::

            sage: FiniteStateMachine().has_final_states()
            False
        """
        return len(self.final_states()) > 0


    #*************************************************************************
    # properties
    #*************************************************************************


    def is_deterministic(self):
        """
        Returns whether the finite finite state machine is deterministic.

        INPUT:

        Nothing.

        OUTPUT:

        ``True`` or ``False``.

        A finite state machine is considered to be deterministic if
        each transition has input label of length one and for each
        pair `(q,a)` where `q` is a state and `a` is an element of the
        input alphabet, there is at most one transition from `q` with
        input label `a`.

        TESTS::

            sage: fsm = FiniteStateMachine()
            sage: fsm.add_transition(('A', 'B', 0, []))
            Transition from 'A' to 'B': 0|-
            sage: fsm.is_deterministic()
            True
            sage: fsm.add_transition(('A', 'C', 0, []))
            Transition from 'A' to 'C': 0|-
            sage: fsm.is_deterministic()
            False
            sage: fsm.add_transition(('A', 'B', [0,1], []))
            Transition from 'A' to 'B': 0,1|-
            sage: fsm.is_deterministic()
            False
        """
        for state in self.iter_states():
            for transition in state.transitions:
                if len(transition.word_in) != 1:
                    return False

            transition_classes_by_word_in = full_group_by(
                state.transitions,
                key=lambda t: t.word_in)

            for key,transition_class in transition_classes_by_word_in:
                if len(transition_class) > 1:
                    return False
        return True


    def is_complete(self):
        """
        Returns whether the finite state machine is complete.

        INPUT:

        Nothing.

        OUTPUT:

        ``True`` or ``False``.

        A finite state machine is considered to be complete if
        each transition has an input label of length one and for each
        pair `(q, a)` where `q` is a state and `a` is an element of the
        input alphabet, there is exactly one transition from `q` with
        input label `a`.

        EXAMPLES::

            sage: fsm = FiniteStateMachine([(0, 0, 0, 0),
            ....:                           (0, 1, 1, 1),
            ....:                           (1, 1, 0, 0)],
            ....:                          determine_alphabets=False)
            sage: fsm.is_complete()
            Traceback (most recent call last):
            ...
            ValueError: No input alphabet is given. Try calling determine_alphabets().
            sage: fsm.input_alphabet = [0, 1]
            sage: fsm.is_complete()
            False
            sage: fsm.add_transition((1, 1, 1, 1))
            Transition from 1 to 1: 1|1
            sage: fsm.is_complete()
            True
            sage: fsm.add_transition((0, 0, 1, 0))
            Transition from 0 to 0: 1|0
            sage: fsm.is_complete()
            False
        """
        if self.input_alphabet is None:
            raise ValueError("No input alphabet is given. "
                             "Try calling determine_alphabets().")

        for state in self.iter_states():
            for transition in state.transitions:
                if len(transition.word_in) != 1:
                    return False

            transition_classes_by_word_in = full_group_by(
                state.transitions,
                key=lambda t: t.word_in)

            for key, transition_class in transition_classes_by_word_in:
                if len(transition_class) > 1:
                    return False

            # all input labels are lists, extract the only element
            outgoing_alphabet = [key[0] for key, transition_class in
                                 transition_classes_by_word_in]
            if not sorted(self.input_alphabet) == sorted(outgoing_alphabet):
                return False

        return True


    def is_connected(self):
        """
        TESTS::

            sage: FiniteStateMachine().is_connected()
            Traceback (most recent call last):
            ...
            NotImplementedError
        """
        raise NotImplementedError


    #*************************************************************************
    # let the finite state machine work
    #*************************************************************************


    def process(self, *args, **kwargs):
        """
        Returns whether the finite state machine accepts the input, the state
        where the computation stops and which output is generated.

        INPUT:

        - ``input_tape`` -- The input tape can be a list with entries from
          the input alphabet.

        - ``initial_state`` -- (default: ``None``) The state in which
          to start. If this parameter is ``None`` and there is only
          one initial state in the machine, then this state is taken.

        OUTPUT:

        A triple, where

        - the first entry is ``True`` if the input string is accepted,

        - the second gives the reached state after processing the
          input tape (This is a state with label ``None`` if the input
          could not be processed, i.e., when at one point no
          transition to go could be found.), and

        - the third gives a list of the output labels used during
          processing (in the case the finite state machine runs as
          transducer).

        Note that in the case the finite state machine is not
        deterministic, one possible path is gone. This means, that in
        this case the output can be wrong.

        EXAMPLES::

            sage: from sage.combinat.finite_state_machine import FSMState
            sage: A = FSMState('A', is_initial = True, is_final = True)
            sage: binary_inverter = FiniteStateMachine({A:[(A, 0, 1), (A, 1, 0)]})
            sage: binary_inverter.process([0, 1, 0, 0, 1, 1])
            (True, 'A', [1, 0, 1, 1, 0, 0])

        Alternatively, we can invoke this function by::

            sage: binary_inverter([0, 1, 0, 0, 1, 1])
            (True, 'A', [1, 0, 1, 1, 0, 0])

        ::

            sage: NAF_ = FSMState('_', is_initial = True, is_final = True)
            sage: NAF1 = FSMState('1', is_final = True)
            sage: NAF = FiniteStateMachine(
            ....:     {NAF_: [(NAF_, 0), (NAF1, 1)], NAF1: [(NAF_, 0)]})
            sage: [NAF.process(w)[0] for w in [[0], [0, 1], [1, 1], [0, 1, 0, 1],
            ....:                           [0, 1, 1, 1, 0], [1, 0, 0, 1, 1]]]
            [True, True, False, True, False, False]
        """
        it = self.iter_process(*args, **kwargs)
        for _ in it:
            pass
        return (it.accept_input, it.current_state, it.output_tape)


    def iter_process(self, input_tape=None, initial_state=None, **kwargs):
        """
        See :meth:`.process` for more informations.

        EXAMPLES::

            sage: inverter = Transducer({'A': [('A', 0, 1), ('A', 1, 0)]},
            ....:     initial_states=['A'], final_states=['A'])
            sage: it = inverter.iter_process(input_tape=[0, 1, 1])
            sage: for _ in it:
            ....:     pass
            sage: it.output_tape
            [1, 0, 0]
        """
        return FSMProcessIterator(self, input_tape, initial_state, **kwargs)


    #*************************************************************************
    # change finite state machine (add/remove state/transitions)
    #*************************************************************************


    def add_state(self, state):
        """
        Adds a state to the finite state machine and returns the new
        state. If the state already exists, that existing state is
        returned.

        INPUT:

        - ``state`` is either an instance of FSMState or, otherwise, a
          label of a state.

        OUTPUT:

        The new or existing state.

        EXAMPLES::

            sage: from sage.combinat.finite_state_machine import FSMState
            sage: F = FiniteStateMachine()
            sage: A = FSMState('A', is_initial=True)
            sage: F.add_state(A)
            'A'
        """
        try:
            return self.state(state)
        except LookupError:
            pass
        # at this point we know that we have a new state
        if is_FSMState(state):
            s = state
        else:
            s = FSMState(state)
        s.transitions = list()
        self._states_.append(s)
        try:
            self._states_dict_[s.label()] = s
        except AttributeError:
            pass
        return s


    def add_states(self, states):
        """
        Adds several states. See add_state for more information.

        INPUT:

        - ``states`` -- a list of states or iterator over states.

        OUTPUT:

        Nothing.

        EXAMPLES::

            sage: F = FiniteStateMachine()
            sage: F.add_states(['A', 'B'])
            sage: F.states()
            ['A', 'B']
        """
        for state in states:
            self.add_state(state)


    def add_transition(self, *args, **kwargs):
        """
        Adds a transition to the finite state machine and returns the
        new transition.

        If the transition already exists, the return value of
        ``self.on_duplicate_transition`` is returned. See the
        documentation of :class:`FiniteStateMachine`.

        INPUT:

        The following forms are all accepted:

        ::

            sage: from sage.combinat.finite_state_machine import FSMState, FSMTransition
            sage: A = FSMState('A')
            sage: B = FSMState('B')

            sage: FSM = FiniteStateMachine()
            sage: FSM.add_transition(FSMTransition(A, B, 0, 1))
            Transition from 'A' to 'B': 0|1

            sage: FSM = FiniteStateMachine()
            sage: FSM.add_transition(A, B, 0, 1)
            Transition from 'A' to 'B': 0|1

            sage: FSM = FiniteStateMachine()
            sage: FSM.add_transition(A, B, word_in=0, word_out=1)
            Transition from 'A' to 'B': 0|1

            sage: FSM = FiniteStateMachine()
            sage: FSM.add_transition('A', 'B', {'word_in': 0, 'word_out': 1})
            Transition from 'A' to 'B': {'word_in': 0, 'word_out': 1}|-

            sage: FSM = FiniteStateMachine()
            sage: FSM.add_transition(from_state=A, to_state=B,
            ....:                    word_in=0, word_out=1)
            Transition from 'A' to 'B': 0|1

            sage: FSM = FiniteStateMachine()
            sage: FSM.add_transition({'from_state': A, 'to_state': B,
            ....:                    'word_in': 0, 'word_out': 1})
            Transition from 'A' to 'B': 0|1

            sage: FSM = FiniteStateMachine()
            sage: FSM.add_transition((A, B, 0, 1))
            Transition from 'A' to 'B': 0|1

            sage: FSM = FiniteStateMachine()
            sage: FSM.add_transition([A, B, 0, 1])
            Transition from 'A' to 'B': 0|1

        If the states ``A`` and ``B`` are not instances of FSMState, then
        it is assumed that they are labels of states.

        OUTPUT:

        The new transition.
        """
        if len(args) + len(kwargs) == 0:
            return
        if len(args) + len(kwargs) == 1:
            if len(args) == 1:
                d = args[0]
                if is_FSMTransition(d):
                    return self._add_fsm_transition_(d)
            else:
                d = kwargs.itervalues().next()
            if hasattr(d, 'iteritems'):
                args = []
                kwargs = d
            elif hasattr(d, '__iter__'):
                args = d
                kwargs = {}
            else:
                raise TypeError("Cannot decide what to do with input.")

        data = dict(zip(
                ('from_state', 'to_state', 'word_in', 'word_out', 'hook'),
                args))
        data.update(kwargs)

        data['from_state'] = self.add_state(data['from_state'])
        data['to_state'] = self.add_state(data['to_state'])

        return self._add_fsm_transition_(FSMTransition(**data))


    def _add_fsm_transition_(self, t):
        """
        Adds a transition.

        INPUT:

        - ``t`` -- an instance of ``FSMTransition``.

        OUTPUT:

        The new transition.

        TESTS::

            sage: from sage.combinat.finite_state_machine import FSMTransition
            sage: F = FiniteStateMachine()
            sage: F._add_fsm_transition_(FSMTransition('A', 'B'))
            Transition from 'A' to 'B': -|-
        """
        try:
            existing_transition = self.transition(t)
        except LookupError:
            pass
        else:
            return self.on_duplicate_transition(existing_transition, t)
        from_state = self.add_state(t.from_state)
        self.add_state(t.to_state)
        from_state.transitions.append(t)
        return t


    def add_from_transition_function(self, function, initial_states=None,
                                     explore_existing_states=True):
        """
        Constructs a finite state machine from a transition function.

        INPUT:

        - ``function`` may return a tuple (new_state, output_word) or a
          list of such tuples.

        - ``initial_states`` -- If no initial states are given, the
          already existing initial states of self are taken.

        - If ``explore_existing_states`` is True (default), then
          already existing states in self (e.g. already given final
          states) will also be processed if they are reachable from
          the initial states.

        OUTPUT:

        Nothing.

        EXAMPLES::

            sage: F = FiniteStateMachine(initial_states=['A'],
            ....:                        input_alphabet=[0, 1])
            sage: def f(state, input):
            ....:     return [('A', input), ('B', 1-input)]
            sage: F.add_from_transition_function(f)
            sage: F.transitions()
            [Transition from 'A' to 'A': 0|0,
            Transition from 'A' to 'B': 0|1,
            Transition from 'A' to 'A': 1|1,
            Transition from 'A' to 'B': 1|0,
            Transition from 'B' to 'A': 0|0,
            Transition from 'B' to 'B': 0|1,
            Transition from 'B' to 'A': 1|1,
            Transition from 'B' to 'B': 1|0]

        Initial states can also be given as a parameter::

            sage: F = FiniteStateMachine(input_alphabet=[0,1])
            sage: def f(state, input):
            ....:     return [('A', input), ('B', 1-input)]
            sage: F.add_from_transition_function(f,initial_states=['A'])
            sage: F.initial_states()
            ['A']

        Already existing states in the finite state machine (the final
        states in the example below) are also explored::

            sage: F = FiniteStateMachine(initial_states=[0],
            ....:                        final_states=[1],
            ....:                        input_alphabet=[0])
            sage: def transition_function(state, letter):
            ....:     return(1-state, [])
            sage: F.add_from_transition_function(transition_function)
            sage: F.transitions()
            [Transition from 0 to 1: 0|-,
             Transition from 1 to 0: 0|-]

        If ``explore_existing_states=False``, however, this behavior
        is turned off, i.e., already existing states are not
        explored::

            sage: F = FiniteStateMachine(initial_states=[0],
            ....:                        final_states=[1],
            ....:                        input_alphabet=[0])
            sage: def transition_function(state, letter):
            ....:     return(1-state, [])
            sage: F.add_from_transition_function(transition_function,
            ....:                                explore_existing_states=False)
            sage: F.transitions()
            [Transition from 0 to 1: 0|-]

        TEST::

            sage: F = FiniteStateMachine(initial_states=['A'])
            sage: def f(state, input):
            ....:     return [('A', input), ('B', 1-input)]
            sage: F.add_from_transition_function(f)
            Traceback (most recent call last):
            ...
            ValueError: No input alphabet is given.
            Try calling determine_alphabets().

        ::

            sage: def transition(state, where):
            ....:     return (vector([0, 0]), 1)
            sage: Transducer(transition, input_alphabet=[0], initial_states=[0])
            Traceback (most recent call last):
            ...
            TypeError: mutable vectors are unhashable
        """
        if self.input_alphabet is None:
            raise ValueError("No input alphabet is given. "
                               "Try calling determine_alphabets().")

        if initial_states is None:
            not_done = self.initial_states()
        elif hasattr(initial_states, '__iter__'):
            not_done = []
            for s in initial_states:
                state = self.add_state(s)
                state.is_initial = True
                not_done.append(state)
        else:
            raise TypeError('Initial states must be iterable ' \
                '(e.g. a list of states).')
        if len(not_done) == 0:
            raise ValueError("No state is initial.")
        if explore_existing_states:
            ignore_done = self.states()
            for s in not_done:
                try:
                    ignore_done.remove(s)
                except ValueError:
                    pass
        else:
            ignore_done = []
        while len(not_done) > 0:
            s = not_done.pop(0)
            for letter in self.input_alphabet:
                try:
                    return_value = function(s.label(), letter)
                except LookupError:
                    continue
                if not hasattr(return_value, "pop"):
                    return_value = [return_value]
                try:
                    for (st_label, word) in return_value:
                        pass
                except TypeError:
                    raise ValueError("The callback function for "
                                     "add_from_transition is expected "
                                     "to return a pair (new_state, "
                                     "output_label) or a list of such pairs. "
                                     "For the state %s and the input "
                                     "letter %s, it however returned %s, "
                                     "which is not acceptable."
                                     % (s.label(), letter, return_value))
                for (st_label, word) in return_value:
                    if not self.has_state(st_label):
                        not_done.append(self.add_state(st_label))
                    elif len(ignore_done) > 0:
                        u = self.state(st_label)
                        if u in ignore_done:
                            not_done.append(u)
                            ignore_done.remove(u)
                    self.add_transition(s, st_label,
                                        word_in=letter, word_out=word)


    def add_transitions_from_function(self, function, labels_as_input=True):
        """
        Adds one or more transitions if ``function(state, state)``
        says that there are some.

        INPUT:

        - ``function`` -- a transition function. Given two states
          ``from_state`` and ``to_state`` (or their labels if
          ``label_as_input`` is true), this function shall return a
          tuple ``(word_in, word_out)`` to add a transition from
          ``from_state`` to ``to_state`` with input and output labels
          ``word_in`` and ``word_out``, respectively. If no such
          addition is to be added, the transition function shall
          return ``None``. The transition function may also return
          a list of such tuples in order to add multiple transitions
          between the pair of states.

        - ``label_as_input`` -- (default: ``True``)

        OUTPUT:

        Nothing.

        EXAMPLES::

            sage: F = FiniteStateMachine()
            sage: F.add_states(['A', 'B', 'C'])
            sage: def f(state1, state2):
            ....:     if state1 == 'C':
            ....:         return None
            ....:     return (0, 1)
            sage: F.add_transitions_from_function(f)
            sage: len(F.transitions())
            6

        Multiple transitions are also possible::

            sage: F = FiniteStateMachine()
            sage: F.add_states([0, 1])
            sage: def f(state1, state2):
            ....:     if state1 != state2:
            ....:          return [(0, 1), (1, 0)]
            ....:     else:
            ....:          return None
            sage: F.add_transitions_from_function(f)
            sage: F.transitions()
            [Transition from 0 to 1: 0|1,
             Transition from 0 to 1: 1|0,
             Transition from 1 to 0: 0|1,
             Transition from 1 to 0: 1|0]

        TESTS::

            sage: F = FiniteStateMachine()
            sage: F.add_state(0)
            0
            sage: def f(state1, state2):
            ....:     return 1
            sage: F.add_transitions_from_function(f)
            Traceback (most recent call last):
            ...
            ValueError: The callback function for add_transitions_from_function
            is expected to return a pair (word_in, word_out) or a list of such
            pairs. For states 0 and 0 however, it returned 1,
            which is not acceptable.

        """
        for s_from in self.iter_states():
            for s_to in self.iter_states():
                try:
                    if labels_as_input:
                        return_value = function(s_from.label(), s_to.label())
                    else:
                        return_value = function(s_from, s_to)
                except LookupError:
                    continue
                if return_value is None:
                    continue
                if not hasattr(return_value, "pop"):
                    transitions = [return_value]
                else:
                    transitions = return_value
                for t in transitions:
                    if not hasattr(t, '__getitem__'):
                         raise ValueError("The callback function for "
                                          "add_transitions_from_function "
                                          "is expected to return a "
                                          "pair (word_in, word_out) or a "
                                          "list of such pairs. For "
                                          "states %s and %s however, it "
                                          "returned %s, which is not "
                                          "acceptable." % (s_from, s_to, return_value))
                    label_in = t[0]
                    try:
                        label_out = t[1]
                    except LookupError:
                        label_out = None
                    self.add_transition(s_from, s_to, label_in, label_out)


    def delete_transition(self, t):
        """
        Deletes a transition by removing it from the list of transitions of
        the state, where the transition starts.

        INPUT:

        - ``t`` -- a transition.

        OUTPUT:

        Nothing.

        EXAMPLES::

            sage: F = FiniteStateMachine([('A', 'B', 0), ('B', 'A', 1)])
            sage: F.delete_transition(('A', 'B', 0))
            sage: F.transitions()
            [Transition from 'B' to 'A': 1|-]
        """
        transition = self.transition(t)
        transition.from_state.transitions.remove(transition)


    def delete_state(self, s):
        """
        Deletes a state and all transitions coming or going to this state.

        INPUT:

        - ``s`` -- a label of a state or an :class:`FSMState`.

        OUTPUT:

        Nothing.

        EXAMPLES::

            sage: from sage.combinat.finite_state_machine import FSMTransition
            sage: t1 = FSMTransition('A', 'B', 0)
            sage: t2 = FSMTransition('B', 'B', 1)
            sage: F = FiniteStateMachine([t1, t2])
            sage: F.delete_state('A')
            sage: F.transitions()
            [Transition from 'B' to 'B': 1|-]

        TESTS::

            sage: F._states_
            ['B']
            sage: F._states_dict_  # This shows that #16024 is fixed.
            {'B': 'B'}
        """
        state = self.state(s)
        for transition in self.transitions():
            if transition.to_state == state:
                self.delete_transition(transition)
        self._states_.remove(state)
        try:
            del self._states_dict_[state.label()]
        except AttributeError:
            pass


    def remove_epsilon_transitions(self):
        """
        TESTS::

            sage: FiniteStateMachine().remove_epsilon_transitions()
            Traceback (most recent call last):
            ...
            NotImplementedError
        """
        raise NotImplementedError


    def accessible_components(self):
        """
        Returns a new finite state machine with the accessible states
        of self and all transitions between those states.

        INPUT:

        Nothing.

        OUTPUT:

        A finite state machine with the accessible states of self and
        all transitions between those states.

        A state is accessible if there is a directed path from an
        initial state to the state. If self has no initial states then
        a copy of the finite state machine self is returned.

        EXAMPLES::

            sage: F = Automaton([(0, 0, 0), (0, 1, 1), (1, 1, 0), (1, 0, 1)],
            ....:               initial_states=[0])
            sage: F.accessible_components()
            Automaton with 2 states

        ::

            sage: F = Automaton([(0, 0, 1), (0, 0, 1), (1, 1, 0), (1, 0, 1)],
            ....:               initial_states=[0])
            sage: F.accessible_components()
            Automaton with 1 states
        """
        if len(self.initial_states()) == 0:
            return deepcopy(self)

        memo = {}
        def accessible(sf, read):
            trans = [x for x in self.transitions(sf) if x.word_in[0] == read]
            return map(lambda x: (deepcopy(x.to_state, memo), x.word_out),
                       trans)

        new_initial_states=map(lambda x: deepcopy(x, memo),
                               self.initial_states())
        result = self.empty_copy()
        result.add_from_transition_function(accessible,
                                            initial_states=new_initial_states)
        for final_state in self.iter_final_states():
            try:
                new_final_state=result.state(final_state.label)
                new_final_state.is_final=True
            except LookupError:
                pass
        return result


    # *************************************************************************
    # creating new finite state machines
    # *************************************************************************


    def disjoint_union(self, other):
        """
        TESTS::

            sage: F = FiniteStateMachine([('A', 'A')])
            sage: FiniteStateMachine().disjoint_union(F)
            Traceback (most recent call last):
            ...
            NotImplementedError
        """
        raise NotImplementedError


    def concatenation(self, other):
        """
        TESTS::

            sage: F = FiniteStateMachine([('A', 'A')])
            sage: FiniteStateMachine().concatenation(F)
            Traceback (most recent call last):
            ...
            NotImplementedError
        """
        raise NotImplementedError


    def Kleene_closure(self):
        """
        TESTS::

            sage: FiniteStateMachine().Kleene_closure()
            Traceback (most recent call last):
            ...
            NotImplementedError
        """
        raise NotImplementedError


    def intersection(self, other):
        """
        TESTS::

            sage: FiniteStateMachine().intersection(FiniteStateMachine())
            Traceback (most recent call last):
            ...
            NotImplementedError
        """
        raise NotImplementedError


    def product_FiniteStateMachine(self, other, function,
                                   new_input_alphabet=None,
                                   only_accessible_components=True,
                                   final_function=None):
        """
        Returns a new finite state machine whose states are
        pairs of states of the original finite state machines.

        INPUT:

        - ``other`` -- a finite state machine

        - ``function`` has to accept two transitions from `A` to `B`
          and `C` to `D` and returns a pair ``(word_in, word_out)``
          which is the label of the transition `(A, C)` to `(B,
          D)`. If there is no transition from `(A, C)` to `(B, D)`,
          then ``function`` should raise a ``LookupError``.

        - ``new_input_alphabet`` (optional)-- the new input alphabet
          as a list.

        - ``only_accessible_components`` -- If true (default), then
          the result is piped through ``accessible_components``. If no
          ``new_input_alphabet`` is given, it is determined by
          :meth:`.determine_alphabets`.

        - ``final_function`` -- A function mapping two final states of
          the original finite state machines to the final output of
          the corresponding state in the new finite state machine. By
          default, the final output is the empty word if both final
          outputs of the constituent states are empty; otherwise, a
          ``ValueError`` is raised.

        OUTPUT:

        A finite state machine whose states are pairs of states of the
        original finite state machines. A state is initial or
        final if both constituent states are initial or final,
        respectively.

        The labels of the transitions are defined by ``function``.

        The final output of a final state is determined by calling
        ``final_function`` on the constituent states.

        The color of a new state is the tuple of colors of the
        constituent states of ``self`` and ``other``.

        EXAMPLES::

            sage: F = Automaton([('A', 'B', 1), ('A', 'A', 0), ('B', 'A', 2)],
            ....:               initial_states=['A'], final_states=['B'],
            ....:               determine_alphabets=True)
            sage: G = Automaton([(1, 1, 1)], initial_states=[1], final_states=[1])
            sage: def addition(transition1, transition2):
            ....:     return (transition1.word_in[0] + transition2.word_in[0],
            ....:             None)
            sage: H = F.product_FiniteStateMachine(G, addition, [0, 1, 2, 3], only_accessible_components=False)
            sage: H.transitions()
            [Transition from ('A', 1) to ('B', 1): 2|-,
             Transition from ('A', 1) to ('A', 1): 1|-,
             Transition from ('B', 1) to ('A', 1): 3|-]
            sage: H1 = F.product_FiniteStateMachine(G, addition, [0, 1, 2, 3], only_accessible_components=False)
            sage: H1.states()[0].label()[0] is F.states()[0]
            True
            sage: H1.states()[0].label()[1] is G.states()[0]
            True

        ::

            sage: F = Automaton([(0,1,1/4), (0,0,3/4), (1,1,3/4), (1,0,1/4)],
            ....:                initial_states=[0] )
            sage: G = Automaton([(0,0,1), (1,1,3/4), (1,0,1/4)],
            ....:                initial_states=[0] )
            sage: H = F.product_FiniteStateMachine(
            ....:         G, lambda t1,t2: (t1.word_in[0]*t2.word_in[0], None))
            sage: H.states()
            [(0, 0), (1, 0)]

        ::

            sage: F = Automaton([(0,1,1/4), (0,0,3/4), (1,1,3/4), (1,0,1/4)],
            ....:                initial_states=[0] )
            sage: G = Automaton([(0,0,1), (1,1,3/4), (1,0,1/4)],
            ....:                initial_states=[0] )
            sage: H = F.product_FiniteStateMachine(G,
            ....:                                  lambda t1,t2: (t1.word_in[0]*t2.word_in[0], None),
            ....:                                  only_accessible_components=False)
            sage: H.states()
            [(0, 0), (1, 0), (0, 1), (1, 1)]

        Also final output words are considered according to the function
        ``final_function``::

            sage: F = Transducer([(0, 1, 0, 1), (1, 1, 1, 1), (1, 1, 0, 1)],
            ....:                final_states=[1])
            sage: F.state(1).final_word_out = 1
            sage: G = Transducer([(0, 0, 0, 1), (0, 0, 1, 0)], final_states=[0])
            sage: G.state(0).final_word_out = 1
            sage: def minus(t1, t2):
            ....:     return (t1.word_in[0] - t2.word_in[0],
            ....:                t1.word_out[0] - t2.word_out[0])
            sage: H = F.product_FiniteStateMachine(G, minus)
            Traceback (most recent call last):
            ...
            ValueError: A final function must be given.
            sage: def plus(s1, s2):
            ....:     return s1.final_word_out[0] + s2.final_word_out[0]
            sage: H = F.product_FiniteStateMachine(G, minus,
            ....:                                  final_function=plus)
            sage: H.final_states()
            [(1, 0)]
            sage: H.final_states()[0].final_word_out
            [2]

        TESTS:

        Check that colors are correctly dealt with. In particular, the
        new colors have to be hashable such that
        :meth:`Automaton.determinisation` does not fail::

            sage: A = Automaton([[0, 0, 0]], initial_states=[0])
            sage: B = A.product_FiniteStateMachine(A,
            ....:                                  lambda t1, t2: (0, None))
            sage: B.states()[0].color
            (None, None)
            sage: B.determinisation()
            Automaton with 1 states

        """
        def default_final_function(s1, s2):
            if s1.final_word_out or s2.final_word_out:
                raise ValueError("A final function must be given.")
            return []

        if final_function is None:
            final_function = default_final_function

        result = self.empty_copy()
        if new_input_alphabet is not None:
            result.input_alphabet = new_input_alphabet
        else:
            result.input_alphabet = None

        for transition1 in self.transitions():
            for transition2 in other.transitions():
                try:
                    word = function(transition1, transition2)
                except LookupError:
                    continue
                result.add_transition((transition1.from_state,
                                       transition2.from_state),
                                      (transition1.to_state,
                                       transition2.to_state),
                                      word[0], word[1])
        for state in result.states():
            if all(map(lambda s: s.is_initial, state.label())):
                state.is_initial = True
            if all(map(lambda s: s.is_final, state.label())):
                state.is_final = True
                state.final_word_out = final_function(*state.label())
            state.color = tuple(map(lambda s: s.color, state.label()))

        if only_accessible_components:
            if result.input_alphabet is None:
                result.determine_alphabets()
            return result.accessible_components()
        else:
            return result


    def composition(self, other, algorithm=None,
                    only_accessible_components=True):
        """
        Returns a new transducer which is the composition of ``self``
        and ``other``.

        INPUT:

        - ``other`` -- a transducer

        - ``algorithm`` -- can be one of the following

          - ``direct`` -- The composition is calculated directly.

            There can be arbitrarily many initial and final states,
            but the input and output labels must have length 1.

            WARNING: The output of other is fed into self.

          - ``explorative`` -- An explorative algorithm is used.

            At least the following restrictions apply, but are not
            checked:
            - both self and other have exactly one initial state
            - all input labels of transitions have length exactly 1

            The input alphabet of self has to be specified.

            This is a very limited implementation of composition.
            WARNING: The output of ``other`` is fed into ``self``.

          If algorithm is ``None``, then the algorithm is chosen
          automatically (at the moment always ``direct``).

        OUTPUT:

        A new transducer.

        The labels of the new finite state machine are pairs of states
        of the original finite state machines. The color of a new
        state is the tuple of colors of the constituent states.


        EXAMPLES::

            sage: F = Transducer([('A', 'B', 1, 0), ('B', 'A', 0, 1)],
            ....:                initial_states=['A', 'B'], final_states=['B'],
            ....:                determine_alphabets=True)
            sage: G = Transducer([(1, 1, 1, 0), (1, 2, 0, 1),
            ....:                 (2, 2, 1, 1), (2, 2, 0, 0)],
            ....:                initial_states=[1], final_states=[2],
            ....:                determine_alphabets=True)
            sage: Hd = F.composition(G, algorithm='direct')
            sage: Hd.initial_states()
            [(1, 'B'), (1, 'A')]
            sage: Hd.transitions()
            [Transition from (1, 'B') to (1, 'A'): 1|1,
             Transition from (1, 'A') to (2, 'B'): 0|0,
             Transition from (2, 'B') to (2, 'A'): 0|1,
             Transition from (2, 'A') to (2, 'B'): 1|0]

        ::

            sage: F = Transducer([('A', 'B', 1, [1, 0]), ('B', 'B', 1, 1),
            ....:                 ('B', 'B', 0, 0)],
            ....:                initial_states=['A'], final_states=['B'])
            sage: G = Transducer([(1, 1, 0, 0), (1, 2, 1, 0),
            ....:                 (2, 2, 0, 1), (2, 1, 1, 1)],
            ....:                initial_states=[1], final_states=[1])
            sage: He = G.composition(F, algorithm='explorative')
            sage: He.transitions()
            [Transition from ('A', 1) to ('B', 2): 1|0,1,
             Transition from ('B', 2) to ('B', 2): 0|1,
             Transition from ('B', 2) to ('B', 1): 1|1,
             Transition from ('B', 1) to ('B', 1): 0|0,
             Transition from ('B', 1) to ('B', 2): 1|0]

        Also final output words are considered if ``algorithm='direct'`` or
        ``None``::

            sage: F = Transducer([('A', 'B', 1, 0), ('B', 'A', 0, 1)],
            ....:                initial_states=['A', 'B'],
            ....:                final_states=['A', 'B'])
            sage: F.state('A').final_word_out = 0
            sage: F.state('B').final_word_out = 1
            sage: G = Transducer([(1, 1, 1, 0), (1, 2, 0, 1),
            ....:                 (2, 2, 1, 1), (2, 2, 0, 0)],
            ....:                initial_states=[1], final_states=[2])
            sage: G.state(2).final_word_out = 0
            sage: Hd = F.composition(G, algorithm='direct')
            sage: Hd.final_states()
            [(2, 'B')]

        Note that ``(2, 'A')`` is not final, as the final output `0`
        of state `2` of `G` cannot be processed in state ``'A'`` of
        `F`.

        ::

            sage: [s.final_word_out for s in Hd.final_states()]
            [[1, 0]]

        Be aware that after composition, different transitions may
        share the same output label (same python object)::

            sage: F = Transducer([ ('A','B',0,0), ('B','A',0,0)],
            ....:                initial_states=['A'],
            ....:                final_states=['A'])
            sage: F.transitions()[0].word_out is F.transitions()[1].word_out
            False
            sage: G = Transducer([('C','C',0,1)],)
            ....:                initial_states=['C'],
            ....:                final_states=['C'])
            sage: H = G.composition(F)
            sage: H.transitions()[0].word_out is H.transitions()[1].word_out
            True

        TESTS:

        Due to the limitations of the two algorithms the following
        (examples from above, but different algorithm used) does not
        give a full answer or does not work

        In the following, ``algorithm='explorative'`` is inadequate,
        as ``F`` has more than one initial state::

            sage: F = Transducer([('A', 'B', 1, 0), ('B', 'A', 0, 1)],
            ....:                initial_states=['A', 'B'], final_states=['B'],
            ....:                determine_alphabets=True)
            sage: G = Transducer([(1, 1, 1, 0), (1, 2, 0, 1),
            ....:                 (2, 2, 1, 1), (2, 2, 0, 0)],
            ....:                initial_states=[1], final_states=[2],
            ....:                determine_alphabets=True)
            sage: He = F.composition(G, algorithm='explorative')
            sage: He.initial_states()
            [(1, 'A')]
            sage: He.transitions()
            [Transition from (1, 'A') to (2, 'B'): 0|0,
             Transition from (2, 'B') to (2, 'A'): 0|1,
             Transition from (2, 'A') to (2, 'B'): 1|0]

        In the following example, ``algorithm='direct'`` is inappropriate
        as there are edges with output labels of length greater than 1::

            sage: F = Transducer([('A', 'B', 1, [1, 0]), ('B', 'B', 1, 1),
            ....:                 ('B', 'B', 0, 0)],
            ....:                initial_states=['A'], final_states=['B'])
            sage: G = Transducer([(1, 1, 0, 0), (1, 2, 1, 0),
            ....:                 (2, 2, 0, 1), (2, 1, 1, 1)],
            ....:                initial_states=[1], final_states=[1])
            sage: Hd = G.composition(F, algorithm='direct')
        """
        if algorithm is None:
            algorithm = 'direct'
        if algorithm == 'direct':
            return self._composition_direct_(other, only_accessible_components)
        elif algorithm == 'explorative':
            return self._composition_explorative_(other)
        else:
            raise ValueError("Unknown algorithm %s." % (algorithm,))


    def _composition_direct_(self, other, only_accessible_components=True):
        """
        See :meth:`.composition` for details.

        TESTS::

            sage: F = Transducer([('A', 'B', 1, 0), ('B', 'A', 0, 1)],
            ....:                initial_states=['A', 'B'], final_states=['B'],
            ....:                determine_alphabets=True)
            sage: G = Transducer([(1, 1, 1, 0), (1, 2, 0, 1),
            ....:                 (2, 2, 1, 1), (2, 2, 0, 0)],
            ....:                initial_states=[1], final_states=[2],
            ....:                determine_alphabets=True)
            sage: Hd = F._composition_direct_(G)
            sage: Hd.initial_states()
            [(1, 'B'), (1, 'A')]
            sage: Hd.transitions()
            [Transition from (1, 'B') to (1, 'A'): 1|1,
             Transition from (1, 'A') to (2, 'B'): 0|0,
             Transition from (2, 'B') to (2, 'A'): 0|1,
             Transition from (2, 'A') to (2, 'B'): 1|0]
        """
        def function(transition1, transition2):
            if transition1.word_out == transition2.word_in:
                return (transition1.word_in, transition2.word_out)
            else:
                raise LookupError

        result = other.product_FiniteStateMachine(
            self, function,
            only_accessible_components=only_accessible_components,
            final_function=lambda s1, s2: [])

        for state_result in result.iter_states():
            state = state_result.label()[0]
            if state.is_final:
                accept, state_to, output = self.process(
                    state.final_word_out,
                    initial_state=self.state(state_result.label()[1]))
                if not accept:
                    state_result.is_final = False
                else:
                    state_result.is_final = True
                    state_result.final_word_out = output

        return result


    def _composition_explorative_(self, other):
        """
        See :meth:`.composition` for details.

        TESTS::

            sage: F = Transducer([('A', 'B', 1, [1, 0]), ('B', 'B', 1, 1),
            ....:                 ('B', 'B', 0, 0)],
            ....:                initial_states=['A'], final_states=['B'])
            sage: G = Transducer([(1, 1, 0, 0), (1, 2, 1, 0),
            ....:                 (2, 2, 0, 1), (2, 1, 1, 1)],
            ....:                initial_states=[1], final_states=[1])
            sage: He = G._composition_explorative_(F)
            sage: He.transitions()
            [Transition from ('A', 1) to ('B', 2): 1|0,1,
             Transition from ('B', 2) to ('B', 2): 0|1,
             Transition from ('B', 2) to ('B', 1): 1|1,
             Transition from ('B', 1) to ('B', 1): 0|0,
             Transition from ('B', 1) to ('B', 2): 1|0]

        Check that colors are correctly dealt with. In particular, the
        new colors have to be hashable such that
        :meth:`Automaton.determinisation` does not fail::

            sage: A = Automaton([[0, 0, 0]], initial_states=[0])
            sage: B = A.composition(A, algorithm='explorative')
            sage: B.states()[0].color
            (None, None)
            sage: B.determinisation()
            Automaton with 1 states

        TODO:

        The explorative algorithm should be re-implemented using the
        process iterators of both finite state machines.
        """
        def composition_transition(state, input):
            (state1, state2) = state
            transition1 = None
            for transition in other.iter_transitions(state1):
                if transition.word_in == [input]:
                    transition1 = transition
                    break
            if transition1 is None:
                raise LookupError
            new_state1 = transition1.to_state
            new_state2 = state2
            output = []
            for o in transition1.word_out:
                transition2 = None
                for transition in self.iter_transitions(new_state2):
                    if transition.word_in == [o]:
                        transition2 = transition
                        break
                if transition2 is None:
                    raise LookupError
                new_state2 = transition2.to_state
                output += transition2.word_out
            return ((new_state1, new_state2), output)

        F = other.empty_copy()
        new_initial_states = [(other.initial_states()[0], self.initial_states()[0])]
        F.add_from_transition_function(composition_transition,
                                       initial_states=new_initial_states)

        for state in F.states():
            if all(map(lambda s: s.is_final, state.label())):
                state.is_final = True
            state.color = tuple(map(lambda s: s.color, state.label()))

        return F


    def input_projection(self):
        """
        Returns an automaton where the output of each transition of
        self is deleted.

        INPUT:

        Nothing

        OUTPUT:

        An automaton.

        EXAMPLES::

            sage: F = FiniteStateMachine([('A', 'B', 0, 1), ('A', 'A', 1, 1),
            ....:                         ('B', 'B', 1, 0)])
            sage: G = F.input_projection()
            sage: G.transitions()
            [Transition from 'A' to 'B': 0|-,
             Transition from 'A' to 'A': 1|-,
             Transition from 'B' to 'B': 1|-]
        """
        return self.projection(what='input')


    def output_projection(self):
        """
        Returns a automaton where the input of each transition of self
        is deleted and the new input is the original output.

        INPUT:

        Nothing

        OUTPUT:

        An automaton.

        EXAMPLES::

            sage: F = FiniteStateMachine([('A', 'B', 0, 1), ('A', 'A', 1, 1),
            ....:                         ('B', 'B', 1, 0)])
            sage: G = F.output_projection()
            sage: G.transitions()
            [Transition from 'A' to 'B': 1|-,
             Transition from 'A' to 'A': 1|-,
             Transition from 'B' to 'B': 0|-]

        Final output words are also considered correctly::

            sage: H = Transducer([('A', 'B', 0, 1), ('A', 'A', 1, 1),
            ....:                 ('B', 'B', 1, 0), ('A', ('final', 0), 0, 0)],
            ....:                final_states=['A', 'B'])
            sage: H.state('B').final_word_out = 2
            sage: J = H.output_projection()
            sage: J.states()
            ['A', 'B', ('final', 0), ('final', 1)]
            sage: J.transitions()
            [Transition from 'A' to 'B': 1|-,
             Transition from 'A' to 'A': 1|-,
             Transition from 'A' to ('final', 0): 0|-,
             Transition from 'B' to 'B': 0|-,
             Transition from 'B' to ('final', 1): 2|-]
            sage: J.final_states()
            ['A', ('final', 1)]
        """
        return self.projection(what='output')


    def projection(self, what='input'):
        """
        Returns an Automaton which transition labels are the projection
        of the transition labels of the input.

        INPUT:

        - ``what`` -- (default: ``input``) either ``input`` or ``output``.

        OUTPUT:

        An automaton.

        EXAMPLES::

            sage: F = FiniteStateMachine([('A', 'B', 0, 1), ('A', 'A', 1, 1),
            ....:                         ('B', 'B', 1, 0)])
            sage: G = F.projection(what='output')
            sage: G.transitions()
            [Transition from 'A' to 'B': 1|-,
             Transition from 'A' to 'A': 1|-,
             Transition from 'B' to 'B': 0|-]
        """
        new = Automaton()
        # TODO: use empty_copy() in order to
        # preserve on_duplicate_transition and future extensions.
        # for this, empty_copy would need a new optional argument
        # use_class=None ?

        if what == 'input':
            new.input_alphabet = copy(self.input_alphabet)
        elif what == 'output':
            new.input_alphabet = copy(self.output_alphabet)
        else:
            raise NotImplementedError

        state_mapping = {}
        for state in self.iter_states():
            state_mapping[state] = new.add_state(deepcopy(state))
        for transition in self.iter_transitions():
            if what == 'input':
                new_word_in = transition.word_in
            elif what == 'output':
                new_word_in = transition.word_out
            else:
                raise NotImplementedError
            new.add_transition((state_mapping[transition.from_state],
                                state_mapping[transition.to_state],
                                new_word_in, None))

        if what == 'output':
            states = [s for s in self.iter_final_states() if s.final_word_out]
            if not states:
                return new
            number = 0
            while new.has_state(('final', number)):
                number += 1
            final = new.add_state(('final', number))
            final.is_final = True
            for state in states:
                output = state.final_word_out
                new.state(state_mapping[state]).final_word_out = []
                new.state(state_mapping[state]).is_final = False
                new.add_transition((state_mapping[state], final, output, None))

        return new


    def transposition(self):
        """
        Returns a new finite state machine, where all transitions of the
        input finite state machine are reversed.

        INPUT:

        Nothing.

        OUTPUT:

        A new finite state machine.

        EXAMPLES::

            sage: aut = Automaton([('A', 'A', 0), ('A', 'A', 1), ('A', 'B', 0)],
            ....:                 initial_states=['A'], final_states=['B'])
            sage: aut.transposition().transitions('B')
            [Transition from 'B' to 'A': 0|-]

        ::

            sage: aut = Automaton([('1', '1', 1), ('1', '2', 0), ('2', '2', 0)],
            ....:                 initial_states=['1'], final_states=['1', '2'])
            sage: aut.transposition().initial_states()
            ['1', '2']
        """
        transposition = self.empty_copy()

        for state in self.states():
            transposition.add_state(deepcopy(state))

        for transition in self.transitions():
            transposition.add_transition(
                transition.to_state.label(), transition.from_state.label(),
                transition.word_in, transition.word_out)

        for initial in self.initial_states():
            state = transposition.state(initial.label())
            if not initial.is_final:
                state.is_final = True
                state.is_initial = False

        for final in self.final_states():
            state = transposition.state(final.label())
            if not final.is_initial:
                state.is_final = False
                state.is_initial = True

        return transposition


    def split_transitions(self):
        """
        Returns a new transducer, where all transitions in self with input
        labels consisting of more than one letter
        are replaced by a path of the corresponding length.

        INPUT:

        Nothing.

        OUTPUT:

        A new transducer.

        EXAMPLES::

            sage: A = Transducer([('A', 'B', [1, 2, 3], 0)],
            ....:                initial_states=['A'], final_states=['B'])
            sage: A.split_transitions().states()
            [('A', ()), ('B', ()),
             ('A', (1,)), ('A', (1, 2))]
        """
        new = self.empty_copy()
        for state in self.states():
            new.add_state(FSMState((state, ()), is_initial=state.is_initial,
                                   is_final=state.is_final))
        for transition in self.transitions():
            for j in range(len(transition.word_in)-1):
                new.add_transition((
                        (transition.from_state, tuple(transition.word_in[:j])),
                        (transition.from_state, tuple(transition.word_in[:j+1])),
                        transition.word_in[j],
                        []))
            new.add_transition((
                    (transition.from_state, tuple(transition.word_in[:-1])),
                    (transition.to_state, ()),
                    transition.word_in[-1:],
                    transition.word_out))
        return new


    def final_components(self):
        """
        Returns the final components of a finite state machine as finite
        state machines.

        INPUT:

        Nothing.

        OUTPUT:

        A list of finite state machines, each representing a final
        component of ``self``.

        A final component of a transducer ``T`` is a strongly connected
        component ``C`` such that there are no transitions of ``T``
        leaving ``C``.

        The final components are the only parts of a transducer which
        influence the main terms of the asympotic behaviour of the sum
        of output labels of a transducer, see [HKP2014]_ and [HKW2014]_.

        EXAMPLES::

            sage: T = Transducer([['A', 'B', 0, 0], ['B', 'C', 0, 1],
            ....:                 ['C', 'B', 0, 1], ['A', 'D', 1, 0],
            ....:                 ['D', 'D', 0, 0], ['D', 'B', 1, 0],
            ....:                 ['A', 'E', 2, 0], ['E', 'E', 0, 0]])
            sage: FC = T.final_components()
            sage: sorted(FC[0].transitions())
            [Transition from 'B' to 'C': 0|1,
             Transition from 'C' to 'B': 0|1]
            sage: FC[1].transitions()
            [Transition from 'E' to 'E': 0|0]

        Another example (cycle of length 2)::

            sage: T = Automaton([[0, 1, 0], [1, 0, 0]])
            sage: len(T.final_components()) == 1
            True
            sage: T.final_components()[0].transitions()
            [Transition from 0 to 1: 0|-,
             Transition from 1 to 0: 0|-]

        REFERENCES:

        .. [HKP2014] Clemens Heuberger, Sara Kropf, and Helmut
           Prodinger, *Asymptotic analysis of the sum of the output of
           transducer*, in preparation.
        """
        DG = self.digraph()
        condensation = DG.strongly_connected_components_digraph()
        final_labels = [v for v in condensation.vertices() if condensation.out_degree(v) == 0]
        return [self.induced_sub_finite_state_machine(map(self.state, component))
                for component in final_labels]


    # *************************************************************************
    # simplifications
    # *************************************************************************


    def prepone_output(self):
        """
        For all paths, shift the output of the path from one
        transition to the earliest possible preceeding transition of
        the path.

        INPUT:

        Nothing.

        OUTPUT:

        Nothing.

        Apply the following to each state `s` (except initial states) of the
        finite state machine as often as possible:

        If the letter `a` is a prefix of the output label of all transitions from
        `s` (including the final output of `s`), then remove it from all these
        labels and append it to all output labels of all transitions leading
        to `s`.

        We assume that the states have no output labels, but final outputs are
        allowed.

        EXAMPLES::

            sage: A = Transducer([('A', 'B', 1, 1),
            ....:                 ('B', 'B', 0, 0),
            ....:                 ('B', 'C', 1, 0)],
            ....:                initial_states=['A'],
            ....:                final_states=['C'])
            sage: A.prepone_output()
            sage: A.transitions()
            [Transition from 'A' to 'B': 1|1,0,
             Transition from 'B' to 'B': 0|0,
             Transition from 'B' to 'C': 1|-]

        ::

            sage: B = Transducer([('A', 'B', 0, 1),
            ....:                 ('B', 'C', 1, [1, 1]),
            ....:                 ('B', 'C', 0, 1)],
            ....:                initial_states=['A'],
            ....:                final_states=['C'])
            sage: B.prepone_output()
            sage: B.transitions()
            [Transition from 'A' to 'B': 0|1,1,
             Transition from 'B' to 'C': 1|1,
             Transition from 'B' to 'C': 0|-]

        If initial states are not labeled as such, unexpected results may be
        obtained::

            sage: C = Transducer([(0,1,0,0)])
            sage: C.prepone_output() # doctest: +ELLIPSIS
            verbose 0 (...: finite_state_machine.py, prepone_output)
            All transitions leaving state 0 have an output label with
            prefix 0.  However, there is no inbound transition and it
            is not an initial state. This routine (possibly called by
            simplification) therefore erased this prefix from all
            outbound transitions.
            sage: C.transitions()
            [Transition from 0 to 1: 0|-]

        Also the final output of final states can be changed::

            sage: T = Transducer([('A', 'B', 0, 1),
            ....:                 ('B', 'C', 1, [1, 1]),
            ....:                 ('B', 'C', 0, 1)],
            ....:                initial_states=['A'],
            ....:                final_states=['B'])
            sage: T.state('B').final_word_out = [1]
            sage: T.prepone_output()
            sage: T.transitions()
            [Transition from 'A' to 'B': 0|1,1,
             Transition from 'B' to 'C': 1|1,
             Transition from 'B' to 'C': 0|-]
            sage: T.state('B').final_word_out
            []

        ::

            sage: S = Transducer([('A', 'B', 0, 1),
            ....:                 ('B', 'C', 1, [1, 1]),
            ....:                 ('B', 'C', 0, 1)],
            ....:                initial_states=['A'],
            ....:                final_states=['B'])
            sage: S.state('B').final_word_out = [0]
            sage: S.prepone_output()
            sage: S.transitions()
            [Transition from 'A' to 'B': 0|1,
             Transition from 'B' to 'C': 1|1,1,
             Transition from 'B' to 'C': 0|1]
            sage: S.state('B').final_word_out
            [0]

        Output labels do not have to be hashable::

            sage: C = Transducer([(0, 1, 0, []),
            ....:                 (1, 0, 0, [vector([0, 0]), 0]),
            ....:                 (1, 1, 1, [vector([0, 0]), 1]),
            ....:                 (0, 0, 1, 0)],
            ....:                 determine_alphabets=False,
            ....:                 initial_states=[0])
            sage: C.prepone_output()
            sage: sorted(C.transitions())
            [Transition from 0 to 1: 0|(0, 0),
             Transition from 0 to 0: 1|0,
             Transition from 1 to 0: 0|0,
             Transition from 1 to 1: 1|1,(0, 0)]
        """
        def find_common_output(state):
            if any(itertools.ifilter(
                    lambda transition: not transition.word_out,
                    self.transitions(state))) \
                   or state.is_final and not state.final_word_out:
                return tuple()
            first_letters = map(lambda transition: transition.word_out[0],
                                self.transitions(state))
            if state.is_final:
                first_letters = first_letters + [state.final_word_out[0]]
            if not first_letters:
                return tuple()
            first_item = first_letters.pop()
            if all([item == first_item for item in first_letters]):
                return (first_item,)
            return tuple()

        changed = 1
        iteration = 0
        while changed > 0:
            changed = 0
            iteration += 1
            for state in self.iter_states():
                if state.is_initial:
                    continue
                if state.word_out:
                    raise NotImplementedError(
                        "prepone_output assumes that all states have "
                        "empty output word, but state %s has output "
                        "word %s" % (state, state.word_out))
                common_output = find_common_output(state)
                if common_output:
                    changed += 1
                    if state.is_final:
                        assert state.final_word_out[0] == common_output[0]
                        state.final_word_out = state.final_word_out[1:]
                    for transition in self.transitions(state):
                        assert transition.word_out[0] == common_output[0]
                        transition.word_out = transition.word_out[1:]
                    found_inbound_transition = False
                    for transition in self.iter_transitions():
                        if transition.to_state == state:
                            transition.word_out = transition.word_out \
                                + [common_output[0]]
                            found_inbound_transition = True
                    if not found_inbound_transition:
                        verbose(
                            "All transitions leaving state %s have an "
                            "output label with prefix %s. However, "
                            "there is no inbound transition and it is "
                            "not an initial state. This routine "
                            "(possibly called by simplification) "
                            "therefore erased this prefix from all "
                            "outbound transitions." %
                            (state, common_output[0]),
                            level=0)


    def equivalence_classes(self):
        r"""
        Returns a list of equivalence classes of states.

        INPUT:

        Nothing.

        OUTPUT:

        A list of equivalence classes of states.

        Two states `a` and `b` are equivalent if and only if there is
        a bijection `\varphi` between paths starting at `a` and paths
        starting at `b` with the following properties: Let `p_a` be a
        path from `a` to `a'` and `p_b` a path from `b` to `b'` such
        that `\varphi(p_a)=p_b`, then

        - `p_a.\mathit{word}_\mathit{in}=p_b.\mathit{word}_\mathit{in}`,
        - `p_a.\mathit{word}_\mathit{out}=p_b.\mathit{word}_\mathit{out}`,
        - `a'` and `b'` have the same output label, and
        - `a'` and `b'` are both final or both non-final.

        The function :meth:`.equivalence_classes` returns a list of
        the equivalence classes to this equivalence relation.

        This is one step of Moore's minimization algorithm.

        .. SEEALSO::

            :meth:`.minimization`

        EXAMPLES::

            sage: fsm = FiniteStateMachine([("A", "B", 0, 1), ("A", "B", 1, 0),
            ....:                           ("B", "C", 0, 0), ("B", "C", 1, 1),
            ....:                           ("C", "D", 0, 1), ("C", "D", 1, 0),
            ....:                           ("D", "A", 0, 0), ("D", "A", 1, 1)])
            sage: fsm.equivalence_classes()
            [['A', 'C'], ['B', 'D']]
        """

        # Two states `a` and `b` are j-equivalent if and only if there
        # is a bijection `\varphi` between paths of length <= j
        # starting at `a` and paths starting at `b` with the following
        # properties: Let `p_a` be a path from `a` to `a'` and `p_b` a
        # path from `b` to `b'` such that `\varphi(p_a)=p_b`, then
        #
        # - `p_a.\mathit{word}_{in}=p_b.\mathit{word}_{in}`,
        # - `p_a.\mathit{word}_{out}=p_b.\mathit{word}_{out}`,
        # - `a'` and `b'` have the same output label, and
        # - `a'` and `b'` are both final or both non-final.

        # If for some j the relations j-1 equivalent and j-equivalent
        # coincide, then they are equal to the equivalence relation
        # described in the docstring.

        # classes_current holds the equivalence classes of
        # j-equivalence, classes_previous holds the equivalence
        # classes of j-1 equivalence.

        # initialize with 0-equivalence
        classes_previous = []
        key_0 = lambda state: (state.is_final, state.color, state.word_out)
        states_grouped = full_group_by(self.states(), key=key_0)
        classes_current = [equivalence_class for
                           (key,equivalence_class) in states_grouped]

        while len(classes_current) != len(classes_previous):
            class_of = {}
            classes_previous = classes_current
            classes_current = []

            for k in range(len(classes_previous)):
                for state in classes_previous[k]:
                    class_of[state] = k

            key_current = lambda state: sorted(
                [(transition.word_in,
                  transition.word_out,
                  class_of[transition.to_state])
                 for transition in state.transitions])

            for class_previous in classes_previous:
                states_grouped = full_group_by(class_previous, key=key_current)
                classes_current.extend([equivalence_class for
                                       (key,equivalence_class) in states_grouped])

        return classes_current


    def quotient(self, classes):
        r"""
        Constructs the quotient with respect to the equivalence
        classes.

        INPUT:

        - ``classes`` is a list of equivalence classes of states.

        OUTPUT:

        A finite state machine.

        The labels of the new states are tuples of states of the
        ``self``, corresponding to ``classes``.

        Assume that `c` is a class, and `a` and `b` are states in
        `c`. Then there is a bijection `\varphi` between the
        transitions from `a` and the transitions from `b` with the
        following properties: if `\varphi(t_a)=t_b`, then

        - `t_a.\mathit{word}_\mathit{in}=t_b.\mathit{word}_\mathit{in}`,
        - `t_a.\mathit{word}_\mathit{out}=t_b.\mathit{word}_\mathit{out}`, and
        - `t_a` and `t_b` lead to some equivalent states `a'` and `b'`.

        Non-initial states may be merged with initial states, the
        resulting state is an initial state.

        All states in a class must have the same ``is_final`` and
        ``word_out`` values.

        EXAMPLES::

            sage: fsm = FiniteStateMachine([("A", "B", 0, 1), ("A", "B", 1, 0),
            ....:                           ("B", "C", 0, 0), ("B", "C", 1, 1),
            ....:                           ("C", "D", 0, 1), ("C", "D", 1, 0),
            ....:                           ("D", "A", 0, 0), ("D", "A", 1, 1)])
            sage: fsmq = fsm.quotient([[fsm.state("A"), fsm.state("C")],
            ....:                      [fsm.state("B"), fsm.state("D")]])
            sage: fsmq.transitions()
            [Transition from ('A', 'C')
                          to ('B', 'D'): 0|1,
             Transition from ('A', 'C')
                          to ('B', 'D'): 1|0,
             Transition from ('B', 'D')
                          to ('A', 'C'): 0|0,
             Transition from ('B', 'D')
                          to ('A', 'C'): 1|1]
            sage: fsmq.relabeled().transitions()
            [Transition from 0 to 1: 0|1,
             Transition from 0 to 1: 1|0,
             Transition from 1 to 0: 0|0,
             Transition from 1 to 0: 1|1]
            sage: fsmq1 = fsm.quotient(fsm.equivalence_classes())
            sage: fsmq1 == fsmq
            True
            sage: fsm.quotient([[fsm.state("A"), fsm.state("B"), fsm.state("C"), fsm.state("D")]])
            Traceback (most recent call last):
                ...
            AssertionError: Transitions of state 'A' and 'B' are incompatible.
        """
        new = self.empty_copy()
        state_mapping = {}

        # Create new states and build state_mapping
        for c in classes:
            new_label = tuple(c)
            new_state = c[0].relabeled(new_label)
            new.add_state(new_state)
            for state in c:
                state_mapping[state] = new_state

        # Copy data from old transducer
        for c in classes:
            new_state = state_mapping[c[0]]
            sorted_transitions = sorted(
                [(state_mapping[t.to_state], t.word_in, t.word_out)
                 for t in c[0].transitions])
            for transition in self.iter_transitions(c[0]):
                new.add_transition(
                    from_state = new_state,
                    to_state = state_mapping[transition.to_state],
                    word_in = transition.word_in,
                    word_out = transition.word_out)

            # check that all class members have the same information (modulo classes)
            for state in c:
                new_state.is_initial = new_state.is_initial or state.is_initial
                assert new_state.is_final == state.is_final, \
                    "Class %s mixes final and non-final states" % (c,)
                assert new_state.word_out == state.word_out, \
                    "Class %s mixes different word_out" % (c,)
                assert new_state.color == state.color, \
                    "Class %s mixes different colors" % (c,)
                assert sorted_transitions == sorted(
                    [(state_mapping[t.to_state], t.word_in, t.word_out)
                     for t in state.transitions]), \
                    "Transitions of state %s and %s are incompatible." % (c[0], state)
        return new


    def merged_transitions(self):
        """
        Merges transitions which have the same ``from_state``,
        ``to_state`` and ``word_out`` while adding their ``word_in``.

        INPUT:

        Nothing.

        OUTPUT:

        A finite state machine with merged transitions. If no mergers occur,
        return ``self``.

        EXAMPLE::

            sage: from sage.combinat.finite_state_machine import duplicate_transition_add_input
            sage: T = Transducer([[1, 2, 1/4, 1], [1, -2, 1/4, 1], [1, -2, 1/2, 1],
            ....:                 [2, 2, 1/4, 1], [2, -2, 1/4, 1], [-2, -2, 1/4, 1],
            ....:                 [-2, 2, 1/4, 1], [2, 3, 1/2, 1], [-2, 3, 1/2, 1]],
            ....:                on_duplicate_transition=duplicate_transition_add_input)
            sage: T1 = T.merged_transitions()
            sage: T1 is T
            False
            sage: sorted(T1.transitions())
            [Transition from -2 to -2: 1/4|1,
             Transition from -2 to 2: 1/4|1,
             Transition from -2 to 3: 1/2|1,
             Transition from 1 to 2: 1/4|1,
             Transition from 1 to -2: 3/4|1,
             Transition from 2 to -2: 1/4|1,
             Transition from 2 to 2: 1/4|1,
             Transition from 2 to 3: 1/2|1]

        Applying the function again does not change the result::

            sage: T2 = T1.merged_transitions()
            sage: T2 is T1
            True
        """

        def key(transition):
            return (transition.to_state, transition.word_out)

        new = self.empty_copy()
        changed = False
        state_dict = {}
        memo = {}

        for state in self.states():
            new_state = deepcopy(state,memo)
            state_dict[state] = new_state
            new.add_state(new_state)

        for state in self.states():
            grouped_transitions = itertools.groupby(sorted(state.transitions, key=key), key=key)
            for (to_state, word_out), transitions in grouped_transitions:
                transition_list = list(transitions)
                changed = changed or len(transition_list) > 1
                word_in = 0
                for transition in transition_list:
                    if hasattr(transition.word_in, '__iter__') and len(transition.word_in) == 1:
                        word_in += transition.word_in[0]
                    else:
                        raise TypeError('%s does not have a list of length 1 as word_in' % transition)
                new.add_transition((state, to_state, word_in, word_out))

        if changed:
            return new
        else:
            return self


    def markov_chain_simplification(self):
        """
        Consider ``self`` as Markov chain with probabilities as input labels
        and simplify it.

        INPUT:

        Nothing.

        OUTPUT:

        Simplified version of ``self``.

        EXAMPLE::

            sage: from sage.combinat.finite_state_machine import duplicate_transition_add_input
            sage: T = Transducer([[1, 2, 1/4, 0], [1, -2, 1/4, 0], [1, -2, 1/2, 0],
            ....:                 [2, 2, 1/4, 1], [2, -2, 1/4, 1], [-2, -2, 1/4, 1],
            ....:                 [-2, 2, 1/4, 1], [2, 3, 1/2, 2], [-2, 3, 1/2, 2]],
            ....:                initial_states=[1],
            ....:                final_states=[3],
            ....:                on_duplicate_transition=duplicate_transition_add_input)
            sage: T1 = T.markov_chain_simplification()
            sage: sorted(T1.transitions())
            [Transition from ((1,),) to ((2, -2),): 1|0,
             Transition from ((2, -2),) to ((2, -2),): 1/2|1,
             Transition from ((2, -2),) to ((3,),): 1/2|2]
        """
        current = self.merged_transitions()
        number_states = len(current.states())

        while True:
            current = current.simplification()
            new_number_states = len(current.states())
            new = current.merged_transitions()
            if new is current and number_states == new_number_states:
                return new
            current = new
            number_states = new_number_states


    def with_final_word_out(self, letters, allow_non_final=True):
        """
        Constructs a new finite state machine with final output words
        for all states by implicitly reading trailing letters until a
        final state is reached.

        INPUT:

        - ``letters`` -- either an element of the input alphabet or a
          list of such elements. This is repeated cyclically when
          needed.

        - ``allow_non_final`` -- a boolean (default: ``True``) which
          indicates whether we allow that some states may be non-final
          in the resulting finite state machine. I.e., if ``False`` then
          each state has to have a path to a final state with input
          label matching ``letters``.

        OUTPUT:

        A finite state machine.

        The inplace version of this function is
        :meth:`.construct_final_word_out`.

        Suppose for the moment a single element ``letter`` as input
        for ``letters``. This is equivalent to ``letters = [letter]``.
        We will discuss the general case below.

        Let ``word_in`` be a word over the input alphabet and assume
        that the original finite state machine transforms ``word_in`` to
        ``word_out`` reaching a possibly non-final state ``s``. Let
        further `k` be the minimum number of letters ``letter`` such
        that there is a path from ``s`` to some final state ``f`` whose
        input label consists of `k` copies of ``letter`` and whose
        output label is ``path_word_out``. Then the state ``s`` of the
        resulting finite state machine is a final state with final
        output ``path_word_out + f.final_word_out``. Therefore, the new
        finite state machine transforms ``word_in`` to ``word_out +
        path_word_out + f.final_word_out``.

        This is e.g. useful for finite state machines operating on digit
        expansions: there, it is sometimes required to read a sufficient
        number of trailing zeros (at the most significant positions) in
        order to reach a final state and to flush all carries. In this
        case, this method constructs an essentially equivalent finite
        state machine in the sense that it not longer requires adding
        sufficiently many trailing zeros. However, it is the
        responsibility of the user to make sure that if adding trailing
        zeros to the input anyway, the output is equivalent.

        If ``letters`` consists of more than one letter, then it is
        assumed that (not necessarily complete) cycles of ``letters``
        are appended as trailing input.

        .. SEEALSO::

            :ref:`example on Gray code <finite_state_machine_gray_code_example>`

        EXAMPLES:

            #.  A simple transducer transforming `00` blocks to `01`
                blocks::

                    sage: T = Transducer([(0, 1, 0, 0), (1, 0, 0, 1)],
                    ....:                initial_states=[0],
                    ....:                final_states=[0])
                    sage: T.process([0, 0, 0])
                    (False, 1, [0, 1, 0])
                    sage: T.process([0, 0, 0, 0])
                    (True, 0, [0, 1, 0, 1])
                    sage: F = T.with_final_word_out(0)
                    sage: for f in F.iter_final_states():
                    ....:     print f, f.final_word_out
                    0 []
                    1 [1]
                    sage: F.process([0, 0, 0])
                    (True, 1, [0, 1, 0, 1])
                    sage: F.process([0, 0, 0, 0])
                    (True, 0, [0, 1, 0, 1])

            #.  A more realistic example: Addition of `1` in binary. We
                construct a transition function transforming the input
                to its binary expansion::

                    sage: def binary_transition(carry, input):
                    ....:     value = carry + input
                    ....:     if value.mod(2) == 0:
                    ....:         return (value/2, 0)
                    ....:     else:
                    ....:         return ((value-1)/2, 1)

                Now, we only have to start with a carry of `1` to
                get the required transducer::

                    sage: T = Transducer(binary_transition,
                    ....:                input_alphabet=[0, 1],
                    ....:                initial_states=[1],
                    ....:                final_states=[0])

                We test this for the binary expansion of `7`::

                    sage: T.process([1, 1, 1])
                    (False, 1, [0, 0, 0])

                The final carry `1` has not be flushed yet, we have to add a
                trailing zero::

                    sage: T.process([1, 1, 1, 0])
                    (True, 0, [0, 0, 0, 1])

                We check that with this trailing zero, the transducer
                performs as advertised::

                    sage: all(ZZ(T(k.bits()+[0]), base=2) == k + 1
                    ....:     for k in srange(16))
                    True

                However, most of the time, we produce superfluous trailing
                zeros::

                    sage: T(11.bits()+[0])
                    [0, 0, 1, 1, 0]

                We now use this method::

                    sage: F = T.with_final_word_out(0)
                    sage: for f in F.iter_final_states():
                    ....:     print f, f.final_word_out
                    1 [1]
                    0 []

                The same tests as above, but we do not have to pad with
                trailing zeros anymore::

                    sage: F.process([1, 1, 1])
                    (True, 1, [0, 0, 0, 1])
                    sage: all(ZZ(F(k.bits()), base=2) == k + 1
                    ....:     for k in srange(16))
                    True

                No more trailing zero in the output::

                    sage: F(11.bits())
                    [0, 0, 1, 1]
                    sage: all(F(k.bits())[-1] == 1
                    ....:     for k in srange(16))
                    True

            #.  Here is an example, where we allow trailing repeated `10`::

                    sage: T = Transducer([(0, 1, 0, 'a'),
                    ....:                 (1, 2, 1, 'b'),
                    ....:                 (2, 0, 0, 'c')],
                    ....:                initial_states=[0],
                    ....:                final_states=[0])
                    sage: F = T.with_final_word_out([1, 0])
                    sage: for f in F.iter_final_states():
                    ....:     print f, ''.join(f.final_word_out)
                    0
                    1 bc

                Trying this with trailing repeated `01` does not produce
                a ``final_word_out`` for state ``1``, but for state ``2``::

                    sage: F = T.with_final_word_out([0, 1])
                    sage: for f in F.iter_final_states():
                    ....:     print f, ''.join(f.final_word_out)
                    0
                    2 c

            #.  Here another example with a more-letter trailing input::

                    sage: T = Transducer([(0, 1, 0, 'a'),
                    ....:                 (1, 2, 0, 'b'), (1, 2, 1, 'b'),
                    ....:                 (2, 3, 0, 'c'), (2, 0, 1, 'e'),
                    ....:                 (3, 1, 0, 'd'), (3, 1, 1, 'd')],
                    ....:                initial_states=[0],
                    ....:                final_states=[0],
                    ....:                with_final_word_out=[0, 0, 1, 1])
                    sage: for f in T.iter_final_states():
                    ....:     print f, ''.join(f.final_word_out)
                    0
                    1 bcdbcdbe
                    2 cdbe
                    3 dbe

        TESTS:

            #.  Reading copies of ``letter`` may result in a cycle. In
                this simple example, we have no final state at all::

                    sage: T = Transducer([(0, 1, 0, 0), (1, 0, 0, 0)],
                    ....:                initial_states=[0])
                    sage: T.with_final_word_out(0)
                    Traceback (most recent call last):
                    ...
                    ValueError: The finite state machine contains
                    a cycle starting at state 0 with input label 0
                    and no final state.

            #.  A unique transition with input word ``letter`` is
                required::

                    sage: T = Transducer([(0, 1, 0, 0), (0, 2, 0, 0)])
                    sage: T.with_final_word_out(0)
                    Traceback (most recent call last):
                    ...
                    ValueError: No unique transition leaving state 0
                    with input label 0.

                It is not a problem if there is no transition starting
                at state ``1`` with input word ``letter``::

                    sage: T = Transducer([(0, 1, 0, 0)])
                    sage: F = T.with_final_word_out(0)
                    sage: for f in F.iter_final_states():
                    ....:     print f, f.final_word_out

                Anyhow, you can override this by::

                    sage: T = Transducer([(0, 1, 0, 0)])
                    sage: T.with_final_word_out(0, allow_non_final=False)
                    Traceback (most recent call last):
                    ...
                    ValueError: No unique transition leaving state 1
                    with input label 0.

            #.  All transitions must have input labels of length `1`::

                    sage: T = Transducer([(0, 0, [], 0)])
                    sage: T.with_final_word_out(0)
                    Traceback (most recent call last):
                    ...
                    NotImplementedError: All transitions must have input
                    labels of length 1. Consider calling split_transitions().
                    sage: T = Transducer([(0, 0, [0, 1], 0)])
                    sage: T.with_final_word_out(0)
                    Traceback (most recent call last):
                    ...
                    NotImplementedError: All transitions must have input
                    labels of length 1. Consider calling split_transitions().

            #.  An empty list as input is not allowed::

                    sage: T = Transducer([(0, 0, [], 0)])
                    sage: T.with_final_word_out([])
                    Traceback (most recent call last):
                    ...
                    ValueError: letters is not allowed to be an empty list.
        """
        new = deepcopy(self)
        new.construct_final_word_out(letters, allow_non_final)
        return new


    def construct_final_word_out(self, letters, allow_non_final=True):
        """
        This is an inplace version of :meth:`.with_final_word_out`. See
        :meth:`.with_final_word_out` for documentation and examples.

        TESTS::

            sage: T = Transducer([(0, 1, 0, 0), (1, 0, 0, 1)],
            ....:                initial_states=[0],
            ....:                final_states=[0])
            sage: F = T.with_final_word_out(0)
            sage: T.construct_final_word_out(0)
            sage: T == F  # indirect doctest
            True
            sage: T = Transducer([(0, 1, 0, None)],
            ....:                final_states=[1])
            sage: F = T.with_final_word_out(0)
            sage: F.state(0).final_word_out
            []
        """
        from itertools import cycle, izip_longest

        if not isinstance(letters, list):
            letters = [letters]
        elif not letters:
            raise ValueError(
                "letters is not allowed to be an empty list.")

        in_progress = set()
        cache = {}

        def find_final_word_out(state):
            # The return value is the output which is produced when
            # reading the given letters until a final state is reached.
            # If no final state can be reached, then None is returned.
            # For final states, the final word out is returned.
            # For final states with empty final output, that is [].
            position, letter = trailing_letters.next()
            if state.is_final:
                return state.final_word_out

            if (state, position) in cache:
                return cache[state, position]

            if (state, position) in in_progress:
                raise ValueError(
                    "The finite state machine contains a cycle "
                    "starting at state %s with input label %s "
                    "and no final state." % (state, letter))

            if any(len(t.word_in) != 1 for t in state.transitions):
                raise NotImplementedError(
                    "All transitions must have input labels of length "
                    "1. Consider calling split_transitions().")

            transitions = [t for t in state.transitions
                           if t.word_in == [letter]]
            if allow_non_final and not transitions:
                final_word_out = None
            elif len(transitions) != 1:
                raise ValueError(
                    "No unique transition leaving state %s with input "
                    "label %s." % (state, letter))
            else:
                in_progress.add((state, position))
                next_word = find_final_word_out(transitions[0].to_state)
                if next_word is not None:
                    final_word_out = transitions[0].word_out + next_word
                else:
                    final_word_out = None
                in_progress.remove((state, position))

            cache[state, position] = final_word_out
            return final_word_out

        for state in self.iter_states():
            assert(not in_progress)
            # trailing_letters is an infinite iterator additionally
            # marking positions
            trailing_letters = cycle(enumerate(letters))
            find_final_word_out(state)

        # actual modifications can only be carried out after all final words
        # have been computed as it may not be permissible to stop at a
        # formerly non-final state unless a cycle has been completed.

        for (state, position), final_word_out in cache.iteritems():
            if position == 0 and final_word_out is not None:
                state.is_final = True
                state.final_word_out = final_word_out


    # *************************************************************************
    # other
    # *************************************************************************


    def graph(self, edge_labels='words_in_out'):
        """
        Returns the graph of the finite state machine with labeled
        vertices and labeled edges.

        INPUT:

        - ``edge_label``: (default: ``'words_in_out'``) can be
             - ``'words_in_out'`` (labels will be strings ``'i|o'``)
             - a function with which takes as input a transition
               and outputs (returns) the label

        OUTPUT:

        A graph.

        EXAMPLES::

            sage: from sage.combinat.finite_state_machine import FSMState
            sage: A = FSMState('A')
            sage: T = Transducer()
            sage: T.graph()
            Digraph on 0 vertices
            sage: T.add_state(A)
            'A'
            sage: T.graph()
            Digraph on 1 vertex
            sage: T.add_transition(('A', 'A', 0, 1))
            Transition from 'A' to 'A': 0|1
            sage: T.graph()
            Looped digraph on 1 vertex
        """
        if edge_labels == 'words_in_out':
            label_fct = lambda t:t._in_out_label_()
        elif hasattr(edge_labels, '__call__'):
            label_fct = edge_labels
        else:
            raise TypeError('Wrong argument for edge_labels.')

        graph_data = []
        isolated_vertices = []
        for state in self.iter_states():
            transitions = state.transitions
            if len(transitions) == 0:
                isolated_vertices.append(state.label())
            for t in transitions:
                graph_data.append((t.from_state.label(), t.to_state.label(),
                                   label_fct(t)))

        G = DiGraph(graph_data)
        G.add_vertices(isolated_vertices)
        return G


    digraph = graph


    def plot(self):
        """
        Plots a graph of the finite state machine with labeled
        vertices and labeled edges.

        INPUT:

        Nothing.

        OUTPUT:

        A plot of the graph of the finite state machine.

        TESTS::

            sage: FiniteStateMachine([('A', 'A', 0)]).plot()
        """
        return self.graph(edge_labels='words_in_out').plot()


    def predecessors(self, state, valid_input=None):
        """
        Lists all predecessors of a state.

        INPUT:

        - ``state`` -- the state from which the predecessors should be
          listed.

        - ``valid_input`` -- If ``valid_input`` is a list, then we
          only consider transitions whose input labels are contained
          in ``valid_input``. ``state`` has to be a :class:`FSMState`
          (not a label of a state). If input labels of length larger
          than `1` are used, then ``valid_input`` has to be a list of
          lists.

        OUTPUT:

        A list of states.

        EXAMPLES::

            sage: A = Transducer([('I', 'A', 'a', 'b'), ('I', 'B', 'b', 'c'),
            ....:                 ('I', 'C', 'c', 'a'), ('A', 'F', 'b', 'a'),
            ....:                 ('B', 'F', ['c', 'b'], 'b'), ('C', 'F', 'a', 'c')],
            ....:                initial_states=['I'], final_states=['F'])
            sage: A.predecessors(A.state('A'))
            ['A', 'I']
            sage: A.predecessors(A.state('F'), valid_input=['b', 'a'])
            ['F', 'C', 'A', 'I']
            sage: A.predecessors(A.state('F'), valid_input=[['c', 'b'], 'a'])
            ['F', 'C', 'B']
        """
        if valid_input is not None:
            valid_list = list()
            for input in valid_input:
                input_list = input
                if not isinstance(input_list, list):
                    input_list = [input]
                valid_list.append(input_list)
            valid_input = valid_list

        unhandeled_direct_predecessors = {s:[] for s in self.states() }
        for t in self.transitions():
            if valid_input is None or t.word_in in valid_input:
                unhandeled_direct_predecessors[t.to_state].append(t.from_state)
        done = []
        open = [state]
        while len(open) > 0:
            s = open.pop()
            candidates = unhandeled_direct_predecessors[s]
            if candidates is not None:
                open.extend(candidates)
                unhandeled_direct_predecessors[s] = None
                done.append(s)
        return(done)

    def asymptotic_moments(self, variable=SR.symbol('n')):
        r"""
        Returns the main terms of expectation and variance of the sum
        of output labels and its covariance with the sum of input
        labels.

        INPUT:

        - ``variable`` -- a symbol denoting the length of the input,
          by default `n`.

        OUTPUT:

        A dictionary consisting of

        - ``expectation`` -- `e n + \operatorname{Order}(1)`,
        - ``variance`` -- `v n + \operatorname{Order}(1)`,
        - ``covariance`` -- `c n + \operatorname{Order}(1)`

        for suitable constants `e`, `v` and `c`.

        Assume that all input and output labels are numbers and that
        ``self`` is complete and has only one final component. Assume
        further that this final component is aperiodic. Furthermore,
        assume that there is exactly one initial state and that all
        states are final.

        Denote by `X_n` the sum of output labels written by the
        finite state machine when reading a random input word of
        length `n` over the input alphabet (assuming
        equidistribution).

        Then the expectation of `X_n` is `en+O(1)`, the variance
        of `X_n` is `vn+O(1)` and the covariance of `X_n` and
        the sum of input labels is `cn+O(1)`, cf. [HKW2014]_,
        Theorem 2.

        In the case of non-integer input or output labels, performance
        degrades significantly. For rational input and output labels,
        consider rescaling to integers. This limitation comes from the
        fact that determinants over polynomial rings can be computed
        much more efficiently than over the symbolic ring. In fact, we
        compute (parts) of a trivariate generating function where the
        input and output labels are exponents of some indeterminates,
        see [HKW2014]_, Theorem 2 for details. If those exponents are
        integers, we can use a polynomial ring.

        .. WARNING::

            If not all states are final, we only print a warning. This is
            for a transitional period to accomodate subsequential
            transducers while those are not yet implemented
            (cf. :trac:`16191`).

        EXAMPLES:

        #.  A trivial example: write the negative of the input::

                sage: T = Transducer([(0, 0, 0, 0), (0, 0, 1, -1)],
                ....:                initial_states=[0],
                ....:                final_states=[0])
                sage: T([0, 1, 1])
                [0, -1, -1]
                sage: moments = T.asymptotic_moments()
                sage: moments['expectation']
                -1/2*n + Order(1)
                sage: moments['variance']
                1/4*n + Order(1)
                sage: moments['covariance']
                -1/4*n + Order(1)

        #.  For the case of the Hamming weight of the non-adjacent-form
            (NAF) of integers, cf. the :wikipedia:`Non-adjacent_form`
            and the :ref:`example on recognizing NAFs
            <finite_state_machine_recognizing_NAFs_example>`, the
            following agrees with the results in [HP2007]_.

            We first use the transducer to convert the standard binary
            expansion to the NAF given in [HP2007]_::

                sage: NAF = Transducer([(0, 0, 0, 0),
                ....:                   (0, '.1', 1, None),
                ....:                   ('.1', 0, 0, [1, 0]),
                ....:                   ('.1', 1, 1, [-1, 0]),
                ....:                   (1, 1, 1, 0),
                ....:                   (1, '.1', 0, None)],
                ....:                  initial_states=[0],
                ....:                  final_states=[0])

            As an example, we compute the NAF of `27` by this
            transducer. Note that we have to add two trailing (at the
            most significant positions) digits `0` in order to be sure
            to reach the final state.

            ::

                sage: binary_27 = 27.bits()
                sage: binary_27
                [1, 1, 0, 1, 1]
                sage: NAF_27 = NAF(binary_27+[0, 0])
                sage: NAF_27
                [-1, 0, -1, 0, 0, 1, 0]
                sage: ZZ(NAF_27, base=2)
                27

            Next, we are only interested in the Hamming weight::

                sage: def weight(state, input):
                ....:     if input is None:
                ....:         result = 0
                ....:     else:
                ....:         result = ZZ(input != 0)
                ....:     return (0, result)
                sage: weight_transducer = Transducer(weight,
                ....:                                input_alphabet=[-1, 0, 1],
                ....:                                initial_states=[0],
                ....:                                final_states=[0])
                sage: NAFweight = weight_transducer.composition(
                ....:     NAF,
                ....:     algorithm='explorative').relabeled()
                sage: sorted(NAFweight.transitions())
                [Transition from 0 to 0: 0|0,
                 Transition from 0 to 1: 1|-,
                 Transition from 1 to 0: 0|1,0,
                 Transition from 1 to 2: 1|1,0,
                 Transition from 2 to 1: 0|-,
                 Transition from 2 to 2: 1|0]
                sage: NAFweight(binary_27 + [0, 0])
                [1, 0, 1, 0, 0, 1, 0]

            Now, we actually compute the asymptotic moments::

                sage: moments = NAFweight.asymptotic_moments()
                verbose 0 (...) Not all states are final. Proceeding
                under the assumption that you know what you are doing.
                sage: moments['expectation']
                1/3*n + Order(1)
                sage: moments['variance']
                2/27*n + Order(1)
                sage: moments['covariance']
                Order(1)

            In this case, we can ignore the warning: we could have all
            states as final states if we would have a final output
            label, i.e., a subsequential transducer. However, this is
            not yet implemented in this package, cf. :trac:`16191`.

        #.  This is Example 3.1 in [HKW2014]_, where a transducer with
            variable output labels is given. There, the aim was to
            choose the output labels of this very simple transducer such
            that the input and output sum are asymptotically
            independent, i.e., the constant `c` vanishes.

            ::

                sage: var('a_1, a_2, a_3, a_4')
                (a_1, a_2, a_3, a_4)
                sage: T = Transducer([[0, 0, 0, a_1], [0, 1, 1, a_3],
                ....:                 [1, 0, 0, a_4], [1, 1, 1, a_2]],
                ....:                initial_states=[0], final_states=[0, 1])
                sage: moments = T.asymptotic_moments() # doctest: +ELLIPSIS
                verbose 0 (...) Non-integer output weights lead to
                significant performance degradation.
                sage: moments['expectation']
                1/4*(a_1 + a_2 + a_3 + a_4)*n + Order(1)
                sage: moments['covariance']
                -1/4*(a_1 - a_2)*n + Order(1)

            Therefore, the asymptotic covariance vanishes if and only if
            `a_2=a_1`.

        #.  This is Example 6.2 in [HKW2014]_, dealing with the
            transducer converting the binary expansion of an integer
            into Gray code (cf. the :wikipedia:`Gray_code` and the
            :ref:`example on Gray code
            <finite_state_machine_gray_code_example>`)::

                sage: moments = transducers.GrayCode().asymptotic_moments()
                verbose 0 (...) Not all states are final. Proceeding
                under the assumption that you know what you are doing.
                sage: moments['expectation']
                1/2*n + Order(1)
                sage: moments['variance']
                1/4*n + Order(1)
                sage: moments['covariance']
                Order(1)

            Also in this case, we can ignore the warning: we could have
            all states as final states if we would have a final output
            label, i.e., a subsequential transducer. However, this is
            not yet implemented in this package, cf. :trac:`16191`.

        #.  This is the first part of Example 6.3 in [HKW2014]_,
            counting the number of 10 blocks in the standard binary
            expansion. The least significant digit is at the left-most
            position::

                sage: block10 = transducers.CountSubblockOccurrences(
                ....:     [1, 0],
                ....:     input_alphabet=[0, 1])
                sage: sorted(block10.transitions())
                [Transition from () to (): 0|0,
                 Transition from () to (1,): 1|0,
                 Transition from (1,) to (): 0|1,
                 Transition from (1,) to (1,): 1|0]
                sage: moments = block10.asymptotic_moments()
                sage: moments['expectation']
                1/4*n + Order(1)
                sage: moments['variance']
                1/16*n + Order(1)
                sage: moments['covariance']
                Order(1)

        #.  This is the second part of Example 6.3 in [HKW2014]_,
            counting the number of 11 blocks in the standard binary
            expansion. The least significant digit is at the left-most
            position::

                sage: block11 = transducers.CountSubblockOccurrences(
                ....:     [1, 1],
                ....:     input_alphabet=[0, 1])
                sage: sorted(block11.transitions())
                [Transition from () to (): 0|0,
                 Transition from () to (1,): 1|0,
                 Transition from (1,) to (): 0|0,
                 Transition from (1,) to (1,): 1|1]
                sage: var('N')
                N
                sage: moments = block11.asymptotic_moments(N)
                sage: moments['expectation']
                1/4*N + Order(1)
                sage: moments['variance']
                5/16*N + Order(1)
                sage: correlation = (moments['covariance'].coefficient(N) /
                ....:                (1/2 * sqrt(moments['variance'].coefficient(N))))
                sage: correlation
                2/5*sqrt(5)

        #.  This is Example 6.4 in [HKW2014]_, counting the number of
            01 blocks minus the number of 10 blocks in the standard binary
            expansion. The least significant digit is at the left-most
            position::

                sage: block01 = transducers.CountSubblockOccurrences(
                ....:     [0, 1],
                ....:     input_alphabet=[0, 1])
                sage: sage.combinat.finite_state_machine.FSMOldCodeTransducerCartesianProduct = False
                sage: product_01x10 = block01.cartesian_product(block10)
                sage: block_difference = transducers.sub([0, 1])(product_01x10)
                sage: T = block_difference.simplification().relabeled()
                sage: sage.combinat.finite_state_machine.FSMOldCodeTransducerCartesianProduct = True
                sage: T.transitions()
                [Transition from 0 to 1: 0|-1,
                 Transition from 0 to 0: 1|0,
                 Transition from 1 to 1: 0|0,
                 Transition from 1 to 0: 1|1,
                 Transition from 2 to 1: 0|0,
                 Transition from 2 to 0: 1|0]
                sage: moments = T.asymptotic_moments()
                sage: moments['expectation']
                Order(1)
                sage: moments['variance']
                Order(1)
                sage: moments['covariance']
                Order(1)

        #.  The finite state machine must have a unique final component::

                sage: T = Transducer([(0, -1, -1, -1), (0, 1, 1, 1),
                ....:                 (-1, -1, -1, -1), (-1, -1, 1, -1),
                ....:                 (1, 1, -1, 1), (1, 1, 1, 1)],
                ....:                initial_states=[0],
                ....:                final_states=[0, 1, -1])
                sage: T.asymptotic_moments()
                Traceback (most recent call last):
                ...
                NotImplementedError: asymptotic_moments is only
                implemented for finite state machines with one final
                component.

            In this particular example, the first letter of the input
            decides whether we reach the loop at `-1` or the loop at
            `1`. In the first case, we have `X_n = -n`, while we have
            `X_n = n` in the second case. Therefore, the expectation
            `E(X_n)` of `X_n` is `E(X_n) = 0`. We get `(X_n-E(X_n))^2 =
            n^2` in all cases, which results in a variance of `n^2`.

            So this example shows that the variance may be non-linear if
            there is more than one final component.

        TESTS:

        #.  An input alphabet must be given::

                sage: T = Transducer([[0, 0, 0, 0]],
                ....:                initial_states=[0], final_states=[0],
                ....:                determine_alphabets=False)
                sage: T.asymptotic_moments()
                Traceback (most recent call last):
                ...
                ValueError: No input alphabet is given.
                Try calling determine_alphabets().

        #.  The finite state machine must have a unique initial state::

                sage: T = Transducer([(0, 0, 0, 0)])
                sage: T.asymptotic_moments()
                Traceback (most recent call last):
                ...
                ValueError: A unique initial state is required.

        #.  The finite state machine must be complete::

                sage: T = Transducer([[0, 0, 0, 0]],
                ....:                initial_states=[0], final_states=[0],
                ....:                input_alphabet=[0, 1])
                sage: T.asymptotic_moments()
                Traceback (most recent call last):
                ...
                NotImplementedError: This finite state machine is
                not complete.

        #.  The final component of the finite state machine must be
            aperiodic::

                sage: T = Transducer([(0, 1, 0, 0), (1, 0, 0, 0)],
                ....:                initial_states=[0], final_states=[0, 1])
                sage: T.asymptotic_moments()
                Traceback (most recent call last):
                ...
                NotImplementedError: asymptotic_moments is only
                implemented for finite state machines whose unique final
                component is aperiodic.

        #.  Non-integer input or output labels lead to a warning::

                sage: T = Transducer([[0, 0, 0, 0], [0, 0, 1, -1/2]],
                ....:                initial_states=[0], final_states=[0])
                sage: moments = T.asymptotic_moments() # doctest: +ELLIPSIS
                verbose 0 (...) Non-integer output weights lead to
                significant performance degradation.
                sage: moments['expectation']
                -1/4*n + Order(1)
                sage: moments['variance']
                1/16*n + Order(1)
                sage: moments['covariance']
                -1/8*n + Order(1)

            This warning can be silenced by :func:`~sage.misc.misc.set_verbose`::

                sage: set_verbose(-1, "finite_state_machine.py")
                sage: moments = T.asymptotic_moments()
                sage: moments['expectation']
                -1/4*n + Order(1)
                sage: moments['variance']
                1/16*n + Order(1)
                sage: moments['covariance']
                -1/8*n + Order(1)
                sage: set_verbose(0, "finite_state_machine.py")

        #.  Check whether ``word_out`` of ``FSMState`` are correctly
            dealt with::

                sage: from sage.combinat.finite_state_machine import FSMState
                sage: s = FSMState(0, word_out=2,
                ....:              is_initial=True,
                ....:              is_final=True)
                sage: T = Transducer([(s, s, 0, 1)],
                ....:                initial_states=[s], final_states=[s])
                sage: T([0, 0])
                [2, 1, 2, 1, 2]
                sage: T.asymptotic_moments()['expectation']
                3*n + Order(1)

            The same test for non-integer output::

                sage: from sage.combinat.finite_state_machine import FSMState
                sage: s = FSMState(0, word_out=2/3)
                sage: T = Transducer([(s, s, 0, 1/2)],
                ....:                initial_states=[s], final_states=[s])
                sage: T.asymptotic_moments()['expectation'] # doctest: +ELLIPSIS
                verbose 0 (...) Non-integer output weights lead to
                significant performance degradation.
                7/6*n + Order(1)

        ALGORITHM:

        See [HKW2014]_, Theorem 2.

        REFERENCES:

        .. [HKW2014] Clemens Heuberger, Sara Kropf and Stephan Wagner,
           *Combinatorial Characterization of Independent Transducers via
           Functional Digraphs*, :arxiv:`1404.3680`.

        .. [HP2007] Clemens Heuberger and Helmut Prodinger, *The Hamming
           Weight of the Non-Adjacent-Form under Various Input Statistics*,
           Periodica Mathematica Hungarica Vol. 55 (1), 2007, pp. 81–96,
           :doi:`10.1007/s10998-007-3081-z`.
        """
        from sage.calculus.functional import derivative
        from sage.rings.polynomial.polynomial_ring_constructor import PolynomialRing
        from sage.rings.rational_field import QQ

        if self.input_alphabet is None:
            raise ValueError("No input alphabet is given. "
                             "Try calling determine_alphabets().")

        if len(self.initial_states()) != 1:
            raise ValueError("A unique initial state is required.")

        if len(self.final_states()) != len(self.states()):
            verbose("Not all states are final. Proceeding under the "
                    "assumption that you know what you are doing.",
                    level=0)

        if not self.is_complete():
            raise NotImplementedError("This finite state machine is "
                                      "not complete.")

        final_components = self.final_components()
        if len(final_components) != 1:
            raise NotImplementedError("asymptotic_moments is only "
                                      "implemented for finite state machines "
                                      "with one final component.")
        final_component = final_components[0]

        if not final_component.digraph().is_aperiodic():
            raise NotImplementedError("asymptotic_moments is only "
                                      "implemented for finite state machines "
                                      "whose unique final component is "
                                      "aperiodic.")

        def get_matrix(fsm, x, y):
            return fsm.adjacency_matrix(
                entry=lambda transition: x**sum(transition.word_in) *
                                         y**(sum(transition.word_out) +
                                         sum(transition.from_state.word_out)))

        K = len(self.input_alphabet)
        R = PolynomialRing(QQ, ("x", "y", "z"))
        (x, y, z) = R.gens()
        try:
            M = get_matrix(self, x, y)
        except TypeError:
            verbose("Non-integer output weights lead to "
                    "significant performance degradation.", level=0)
            # fall back to symbolic ring
            R = SR
            x = R.symbol()
            y = R.symbol()
            z = R.symbol()
            M = get_matrix(self, x, y)
            def substitute_one(g):
                return g.subs({x: 1, y: 1, z: 1})
        else:
            def substitute_one(g):
                # the result of the substitution shall live in QQ,
                # not in the polynomial ring R, so the method
                # subs does not achieve the result.
                # Therefore, we need this helper function.
                return g(1, 1, 1)

        f = (M.parent().identity_matrix() - z/K*M).det()
        f_x = substitute_one(derivative(f, x))
        f_y = substitute_one(derivative(f, y))
        f_z = substitute_one(derivative(f, z))
        f_xy = substitute_one(derivative(f, x, y))
        f_xz = substitute_one(derivative(f, x, z))
        f_yz = substitute_one(derivative(f, y, z))
        f_yy = substitute_one(derivative(f, y, y))
        f_zz = substitute_one(derivative(f, z, z))

        e_2 = f_y / f_z
        v_2 = (f_y**2 * (f_zz+f_z) + f_z**2 * (f_yy+f_y)
               - 2*f_y*f_z*f_yz) / f_z**3
        c = (f_x * f_y * (f_zz+f_z) + f_z**2 * f_xy - f_y*f_z*f_xz
             - f_x*f_z*f_yz) / f_z**3

        return {'expectation': e_2*variable + SR(1).Order(),
                'variance': v_2*variable + SR(1).Order(),
                'covariance': c*variable + SR(1).Order()}


#*****************************************************************************


def is_Automaton(FSM):
    """
    Tests whether or not ``FSM`` inherits from :class:`Automaton`.

    TESTS::

        sage: from sage.combinat.finite_state_machine import is_FiniteStateMachine, is_Automaton
        sage: is_Automaton(FiniteStateMachine())
        False
        sage: is_Automaton(Automaton())
        True
        sage: is_FiniteStateMachine(Automaton())
        True
    """
    return isinstance(FSM, Automaton)


class Automaton(FiniteStateMachine):
    """
    This creates an automaton, which is a finite state machine, whose
    transitions have input labels.

    An automaton has additional features like creating a deterministic
    and a minimized automaton.

    See class :class:`FiniteStateMachine` for more information.

    EXAMPLES:

    We can create an automaton recognizing even numbers (given in
    binary and read from left to right) in the following way::

        sage: A = Automaton([('P', 'Q', 0), ('P', 'P', 1),
        ....:                ('Q', 'P', 1), ('Q', 'Q', 0)],
        ....:               initial_states=['P'], final_states=['Q'])
        sage: A
        Automaton with 2 states
        sage: A([0])
        True
        sage: A([1, 1, 0])
        True
        sage: A([1, 0, 1])
        False

    Note that the full output of the commands can be obtained by
    calling :meth:`.process` and looks like this::

        sage: A.process([1, 0, 1])
        (False, 'P')

    TESTS::

        sage: Automaton()
        Automaton with 0 states
    """

    def _repr_(self):
        """
        Represents the finite state machine as "Automaton with n
        states" where n is the number of states.

        INPUT:

        Nothing.

        OUTPUT:

        A string.

        EXAMPLES::

            sage: Automaton()._repr_()
            'Automaton with 0 states'
        """
        return "Automaton with %s states" % len(self._states_)

    def _latex_transition_label_(self, transition, format_function=latex):
        r"""
        Returns the proper transition label.

        INPUT:

        - ``transition`` - a transition

        - ``format_function`` - a function formatting the labels

        OUTPUT:

        A string.

        EXAMPLES::

            sage: F = Automaton([('A', 'B', 1)])
            sage: print latex(F)  # indirect doctest
            \begin{tikzpicture}[auto, initial text=, >=latex]
            \node[state] (v0) at (3.000000, 0.000000) {$\text{\texttt{A}}$};
            \node[state] (v1) at (-3.000000, 0.000000) {$\text{\texttt{B}}$};
            \path[->] (v0) edge node[rotate=360.00, anchor=south] {$1$} (v1);
            \end{tikzpicture}

        TESTS::

            sage: F = Automaton([('A', 'B', 0, 1)])
            sage: t = F.transitions()[0]
            sage: F._latex_transition_label_(t)
            \left[0\right]
        """
        return format_function(transition.word_in)


    def intersection(self, other, only_accessible_components=True):
        """
        Returns a new automaton which accepts an input if it is
        accepted by both given automata.

        INPUT:

        - ``other`` -- an automaton

        - ``only_accessible_components`` -- If ``True`` (default), then
          the result is piped through ``accessible_components``. If no
          ``new_input_alphabet`` is given, it is determined by
          :meth:`.determine_alphabets`.

        OUTPUT:

        A new automaton which computes the intersection
        (see below) of the languages of ``self`` and ``other``.

        The set of states of the new automaton is the cartesian product of the
        set of states of both given automata. There is a transition `((A, B),
        (C, D), a)` in the new automaton if there are transitions `(A, C, a)`
        and `(B, D, a)` in the old automata.

        The methods :meth:`.intersection` and
        :meth:`.cartesian_product` are the same (for automata).

        EXAMPLES::

            sage: aut1 = Automaton([('1', '2', 1),
            ....:                   ('2', '2', 1),
            ....:                   ('2', '2', 0)],
            ....:                  initial_states=['1'],
            ....:                  final_states=['2'],
            ....:                  determine_alphabets=True)
            sage: aut2 = Automaton([('A', 'A', 1),
            ....:                   ('A', 'B', 0),
            ....:                   ('B', 'B', 0),
            ....:                   ('B', 'A', 1)],
            ....:                  initial_states=['A'],
            ....:                  final_states=['B'],
            ....:                  determine_alphabets=True)
            sage: res = aut1.intersection(aut2)
            sage: (aut1([1, 1]), aut2([1, 1]), res([1, 1]))
            (True, False, False)
            sage: (aut1([1, 0]), aut2([1, 0]), res([1, 0]))
            (True, True, True)
            sage: res.transitions()
            [Transition from ('1', 'A') to ('2', 'A'): 1|-,
             Transition from ('2', 'A') to ('2', 'B'): 0|-,
             Transition from ('2', 'A') to ('2', 'A'): 1|-,
             Transition from ('2', 'B') to ('2', 'B'): 0|-,
             Transition from ('2', 'B') to ('2', 'A'): 1|-]

        For automata with epsilon-transitions, intersection is not well
        defined. But for any finite state machine, epsilon-transitions can be
        removed by :meth:`.remove_epsilon_transitions`.

        ::

            sage: a1 = Automaton([(0, 0, 0),
            ....:                 (0, 1, None),
            ....:                 (1, 1, 1),
            ....:                 (1, 2, 1)],
            ....:                 initial_states=[0],
            ....:                 final_states=[1],
            ....:                 determine_alphabets=True)
            sage: a2 = Automaton([(0, 0, 0), (0, 1, 1), (1, 1, 1)],
            ....:                 initial_states=[0],
            ....:                 final_states=[1],
            ....:                 determine_alphabets=True)
            sage: a1.intersection(a2)
            Traceback (most recent call last):
            ...
            ValueError: An epsilon-transition (with empty input)
            was found.
            sage: a1.remove_epsilon_transitions()  # not tested (since not implemented yet)
            sage: a1.intersection(a2)  # not tested
        """
        if not is_Automaton(other):
            raise TypeError(
                 "Only an automaton can be intersected with an automaton.")

        def function(transition1, transition2):
            if not transition1.word_in or not transition2.word_in:
                raise ValueError(
                    "An epsilon-transition (with empty input) was found.")
            if transition1.word_in == transition2.word_in:
                return (transition1.word_in, None)
            else:
                raise LookupError

        return self.product_FiniteStateMachine(
            other,
            function,
            only_accessible_components=only_accessible_components)

    cartesian_product = intersection

    def determinisation(self):
        """
        Returns a deterministic automaton which accepts the same input
        words as the original one.

        INPUT:

        Nothing.

        OUTPUT:

        A new automaton, which is deterministic.

        The labels of the states of the new automaton are frozensets
        of states of ``self``. The color of a new state is the
        frozenset of colors of the constituent states of ``self``.
        Therefore, the colors of the constituent states have to be
        hashable.

        The input alphabet must be specified.

        EXAMPLES::

            sage: aut = Automaton([('A', 'A', 0), ('A', 'B', 1), ('B', 'B', 1)],
            ....:                 initial_states=['A'], final_states=['B'])
            sage: aut.determinisation().transitions()
            [Transition from frozenset(['A'])
                          to frozenset(['A']): 0|-,
             Transition from frozenset(['A'])
                          to frozenset(['B']): 1|-,
             Transition from frozenset(['B'])
                          to frozenset([]): 0|-,
             Transition from frozenset(['B'])
                          to frozenset(['B']): 1|-,
             Transition from frozenset([])
                          to frozenset([]): 0|-,
             Transition from frozenset([])
                          to frozenset([]): 1|-]

        ::

            sage: A = Automaton([('A', 'A', 1), ('A', 'A', 0), ('A', 'B', 1),
            ....:                ('B', 'C', 0), ('C', 'C', 1), ('C', 'C', 0)],
            ....:               initial_states=['A'], final_states=['C'])
            sage: A.determinisation().states()
            [frozenset(['A']), frozenset(['A', 'B']),
            frozenset(['A', 'C']), frozenset(['A', 'C', 'B'])]

        ::

            sage: A = Automaton([(0, 1, 1), (0, 2, [1, 1]), (0, 3, [1, 1, 1]),
            ....:                (1, 0, -1), (2, 0, -2), (3, 0, -3)],
            ....:               initial_states=[0], final_states=[0, 1, 2, 3])
            sage: B = A.determinisation().relabeled()
            sage: all(t.to_state.label() == 2 for t in
            ....:     B.state(2).transitions)
            True
            sage: B.state(2).is_final
            False
            sage: B.delete_state(2)  # this is a sink
            sage: sorted(B.transitions())
            [Transition from 0 to 1: 1|-,
             Transition from 1 to 0: -1|-,
             Transition from 1 to 3: 1|-,
             Transition from 3 to 0: -2|-,
             Transition from 3 to 4: 1|-,
             Transition from 4 to 0: -3|-]

        Note that colors of states have to be hashable::

            sage: A = Automaton([[0, 0, 0]], initial_states=[0])
            sage: A.state(0).color = []
            sage: A.determinisation()
            Traceback (most recent call last):
            ...
            TypeError: unhashable type: 'list'
            sage: A.state(0).color = ()
            sage: A.determinisation()
            Automaton with 1 states

        TESTS:

        This is from #15078, comment 13.

        ::

            sage: D = {'A': [('A', 'a'), ('B', 'a'), ('A', 'b')],
            ....:      'C': [], 'B': [('C', 'b')]}
            sage: auto = Automaton(D, initial_states=['A'], final_states=['C'])
            sage: auto.is_deterministic()
            False
            sage: auto.process(list('aaab'))
            (False, 'A')
            sage: auto.states()
            ['A', 'C', 'B']
            sage: Ddet = auto.determinisation()
            sage: Ddet
            Automaton with 3 states
            sage: Ddet.is_deterministic()
            True
            sage: sorted(Ddet.transitions())
            [Transition from frozenset(['A']) to frozenset(['A', 'B']): 'a'|-,
             Transition from frozenset(['A']) to frozenset(['A']): 'b'|-,
             Transition from frozenset(['A', 'B']) to frozenset(['A', 'B']): 'a'|-,
             Transition from frozenset(['A', 'B']) to frozenset(['A', 'C']): 'b'|-,
             Transition from frozenset(['A', 'C']) to frozenset(['A', 'B']): 'a'|-,
             Transition from frozenset(['A', 'C']) to frozenset(['A']): 'b'|-]
            sage: Ddet.initial_states()
            [frozenset(['A'])]
            sage: Ddet.final_states()
            [frozenset(['A', 'C'])]
        """
        if any(len(t.word_in) > 1 for t in self.iter_transitions()):
            return self.split_transitions().determinisation()

        epsilon_successors = {}
        direct_epsilon_successors = {}
        for state in self.iter_states():
            direct_epsilon_successors[state] = set(
                t.to_state
                for t in self.iter_transitions(state)
                if not t.word_in)
            epsilon_successors[state] = set([state])

        old_count_epsilon_successors = 0
        count_epsilon_successors = len(epsilon_successors)

        while old_count_epsilon_successors < count_epsilon_successors:
            old_count_epsilon_successors = count_epsilon_successors
            count_epsilon_successors = 0
            for state in self.iter_states():
                for direct_successor in direct_epsilon_successors[state]:
                    epsilon_successors[state] = epsilon_successors[state].union(epsilon_successors[direct_successor])
                count_epsilon_successors += len(epsilon_successors[state])

        def set_transition(states, letter):
            result = set()
            for state in states:
                for transition in self.iter_transitions(state):
                    if transition.word_in == [letter]:
                        result.add(transition.to_state)
            result = result.union(*(epsilon_successors[s] for s in result))
            return (frozenset(result), [])

        result = self.empty_copy()
        new_initial_states = [frozenset(self.iter_initial_states())]
        result.add_from_transition_function(set_transition,
                                            initial_states=new_initial_states)

        for state in result.iter_states():
            state.is_final = any(s.is_final for s in state.label())
            state.color = frozenset(s.color for s in state.label())

        return result


    def minimization(self, algorithm=None):
        """
        Returns the minimization of the input automaton as a new automaton.

        INPUT:

        - ``algorithm`` -- Either Moore's algorithm (by
          ``algorithm='Moore'`` or as default for deterministic
          automata) or Brzozowski's algorithm (when
          ``algorithm='Brzozowski'`` or when the automaton is not
          deterministic) is used.

        OUTPUT:

        A new automaton.

        The resulting automaton is deterministic and has a minimal
        number of states.

        EXAMPLES::

            sage: A = Automaton([('A', 'A', 1), ('A', 'A', 0), ('A', 'B', 1),
            ....:                ('B', 'C', 0), ('C', 'C', 1), ('C', 'C', 0)],
            ....:               initial_states=['A'], final_states=['C'])
            sage: B = A.minimization(algorithm='Brzozowski')
            sage: B.transitions(B.states()[1])
            [Transition from frozenset([frozenset(['A', 'C', 'B']),
            frozenset(['C', 'B']), frozenset(['A', 'C'])]) to
            frozenset([frozenset(['A', 'C', 'B']), frozenset(['C', 'B']),
            frozenset(['A', 'C']), frozenset(['C'])]): 0|-,
            Transition from frozenset([frozenset(['A', 'C', 'B']),
            frozenset(['C', 'B']), frozenset(['A', 'C'])]) to
            frozenset([frozenset(['A', 'C', 'B']), frozenset(['C', 'B']),
            frozenset(['A', 'C'])]): 1|-]
            sage: len(B.states())
            3
            sage: C = A.minimization(algorithm='Brzozowski')
            sage: C.transitions(C.states()[1])
            [Transition from frozenset([frozenset(['A', 'C', 'B']),
            frozenset(['C', 'B']), frozenset(['A', 'C'])]) to
            frozenset([frozenset(['A', 'C', 'B']), frozenset(['C', 'B']),
            frozenset(['A', 'C']), frozenset(['C'])]): 0|-,
            Transition from frozenset([frozenset(['A', 'C', 'B']),
            frozenset(['C', 'B']), frozenset(['A', 'C'])]) to
            frozenset([frozenset(['A', 'C', 'B']), frozenset(['C', 'B']),
            frozenset(['A', 'C'])]): 1|-]
            sage: len(C.states())
            3

        ::

            sage: aut = Automaton([('1', '2', 'a'), ('2', '3', 'b'),
            ....:                  ('3', '2', 'a'), ('2', '1', 'b'),
            ....:                  ('3', '4', 'a'), ('4', '3', 'b')],
            ....:                  initial_states=['1'], final_states=['1'])
            sage: min = aut.minimization(algorithm='Brzozowski')
            sage: [len(min.states()), len(aut.states())]
            [3, 4]
            sage: min = aut.minimization(algorithm='Moore')
            Traceback (most recent call last):
            ...
            NotImplementedError: Minimization via Moore's Algorithm is only
            implemented for deterministic finite state machines
        """
        deterministic = self.is_deterministic()

        if algorithm == "Moore" or (algorithm is None and deterministic):
            return self._minimization_Moore_()
        elif algorithm == "Brzozowski" or (algorithm is None and not deterministic):
            return self._minimization_Brzozowski_()
        else:
            raise NotImplementedError("Algorithm '%s' is not implemented. Choose 'Moore' or 'Brzozowski'" % algorithm)


    def _minimization_Brzozowski_(self):
        """
        Returns a minimized automaton by using Brzozowski's algorithm.

        See also :meth:`.minimization`.

        TESTS::

            sage: A = Automaton([('A', 'A', 1), ('A', 'A', 0), ('A', 'B', 1),
            ....:                ('B', 'C', 0), ('C', 'C', 1), ('C', 'C', 0)],
            ....:               initial_states=['A'], final_states=['C'])
            sage: B = A._minimization_Brzozowski_()
            sage: len(B.states())
            3
        """
        return self.transposition().determinisation().transposition().determinisation()


    def _minimization_Moore_(self):
        """
        Returns a minimized automaton by using Moore's algorithm.

        See also :meth:`.minimization`.

        TESTS::

            sage: aut = Automaton([('1', '2', 'a'), ('2', '3', 'b'),
            ....:                  ('3', '2', 'a'), ('2', '1', 'b'),
            ....:                  ('3', '4', 'a'), ('4', '3', 'b')],
            ....:                  initial_states=['1'], final_states=['1'])
            sage: min = aut._minimization_Moore_()
            Traceback (most recent call last):
            ...
            NotImplementedError: Minimization via Moore's Algorithm is only
            implemented for deterministic finite state machines
        """
        if self.is_deterministic():
            return self.quotient(self.equivalence_classes())
        else:
            raise NotImplementedError("Minimization via Moore's Algorithm is only " \
                                      "implemented for deterministic finite state machines")


    def process(self, *args, **kwargs):
        """
        .. WARNING::

            The default output of this method is scheduled to change.
            This docstring describes the new default behaviour, which can
            already be achieved by setting
            ``FSMOldProcessOutput`` to ``False``.

        Returns whether the automaton accepts the input and the state
        where the computation stops.

        INPUT:

        - ``input_tape`` -- The input tape can be a list with entries from
          the input alphabet.

        - ``initial_state`` -- (default: ``None``) The state in which
          to start. If this parameter is ``None`` and there is only
          one initial state in the machine, then this state is taken.

        - ``full_output`` -- (default: ``True``) If set, then the full
          output is given, otherwise only whether the sequence is accepted
          or not (the first entry below only).

        OUTPUT:

        The full output is a pair, where

        - the first entry is ``True`` if the input string is accepted and

        - the second gives the state reached after processing the
          input tape (This is a state with label ``None`` if the input
          could not be processed, i.e., when at one point no
          transition to go could be found.).

        Note that in the case the automaton is not
        deterministic, one possible path is gone. This means that in
        this case the output can be wrong. Use
        :meth:`.determinisation` to get a deterministic automaton
        machine and try again.

        By setting ``FSMOldProcessOutput`` to ``False``
        the new desired output is produced.

        EXAMPLES::

            sage: sage.combinat.finite_state_machine.FSMOldProcessOutput = False  # activate new output behavior
            sage: from sage.combinat.finite_state_machine import FSMState
            sage: NAF_ = FSMState('_', is_initial = True, is_final = True)
            sage: NAF1 = FSMState('1', is_final = True)
            sage: NAF = Automaton(
            ....:     {NAF_: [(NAF_, 0), (NAF1, 1)], NAF1: [(NAF_, 0)]})
            sage: [NAF.process(w) for w in [[0], [0, 1], [1, 1], [0, 1, 0, 1],
            ....:                           [0, 1, 1, 1, 0], [1, 0, 0, 1, 1]]]
            [(True, '_'), (True, '1'), (False, None),
             (True, '1'), (False, None), (False, None)]

        If we just want a condensed output, we use::

            sage: [NAF.process(w, full_output=False)
            ....:     for w in [[0], [0, 1], [1, 1], [0, 1, 0, 1],
            ....:               [0, 1, 1, 1, 0], [1, 0, 0, 1, 1]]]
            [True, True, False, True, False, False]

        It is equivalent to::

            sage: [NAF(w) for w in [[0], [0, 1], [1, 1], [0, 1, 0, 1],
            ....:                   [0, 1, 1, 1, 0], [1, 0, 0, 1, 1]]]
            [True, True, False, True, False, False]

        The following example illustrates the difference between
        non-existing paths and reaching a non-final state::

            sage: NAF.process([2])
            (False, None)
            sage: NAF.add_transition(('_', 's', 2))
            Transition from '_' to 's': 2|-
            sage: NAF.process([2])
            (False, 's')
        """
        if FSMOldProcessOutput:
            from sage.misc.superseded import deprecation
            deprecation(16132, "The output of Automaton.process "
                               "(and thus of Automaton.__call__) "
                               "will change. Please use the corresponding "
                               "functions from FiniteStateMachine "
                               "for the original output.")
            return super(Automaton, self).process(*args, **kwargs)

        if not kwargs.has_key('full_output'):
            kwargs['full_output'] = True

        it = self.iter_process(*args, **kwargs)
        for _ in it:
            pass

        # process output
        if kwargs['full_output']:
            return (it.accept_input, it.current_state)
        else:
            return it.accept_input


#*****************************************************************************


def is_Transducer(FSM):
    """
    Tests whether or not ``FSM`` inherits from :class:`Transducer`.

    TESTS::

        sage: from sage.combinat.finite_state_machine import is_FiniteStateMachine, is_Transducer
        sage: is_Transducer(FiniteStateMachine())
        False
        sage: is_Transducer(Transducer())
        True
        sage: is_FiniteStateMachine(Transducer())
        True
    """
    return isinstance(FSM, Transducer)


class Transducer(FiniteStateMachine):
    """
    This creates a transducer, which is a finite state machine, whose
    transitions have input and output labels.

    An transducer has additional features like creating a simplified
    transducer.

    See class :class:`FiniteStateMachine` for more information.

    EXAMPLES:

    We can create a transducer performing the addition of 1 (for
    numbers given in binary and read from right to left) in the
    following way::

        sage: T = Transducer([('C', 'C', 1, 0), ('C', 'N', 0, 1),
        ....:                 ('N', 'N', 0, 0), ('N', 'N', 1, 1)],
        ....:                initial_states=['C'], final_states=['N'])
        sage: T
        Transducer with 2 states
        sage: T([0])
        [1]
        sage: T([1,1,0])
        [0, 0, 1]
        sage: ZZ(T(15.digits(base=2)+[0]), base=2)
        16

    Note that we have padded the binary input sequence by a `0` so
    that the transducer can reach its final state.

    TESTS::

        sage: Transducer()
        Transducer with 0 states
    """

    def _repr_(self):
        """
        Represents the transducer as "Transducer with n states" where
        n is the number of states.

        INPUT:

        Nothing.

        OUTPUT:

        A string.

        EXAMPLES::

            sage: Transducer()._repr_()
            'Transducer with 0 states'
        """
        return "Transducer with %s states" % len(self._states_)

    def _latex_transition_label_(self, transition, format_function=latex):
        r"""
        Returns the proper transition label.

        INPUT:

        - ``transition`` - a transition

        - ``format_function`` - a function formatting the labels

        OUTPUT:

        A string.

        EXAMPLES::

            sage: F = Transducer([('A', 'B', 1, 2)])
            sage: print latex(F)  # indirect doctest
            \begin{tikzpicture}[auto, initial text=, >=latex]
            \node[state] (v0) at (3.000000, 0.000000) {$\text{\texttt{A}}$};
            \node[state] (v1) at (-3.000000, 0.000000) {$\text{\texttt{B}}$};
            \path[->] (v0) edge node[rotate=360.00, anchor=south] {$1\mid 2$} (v1);
            \end{tikzpicture}

        TESTS::

            sage: F = Transducer([('A', 'B', 0, 1)])
            sage: t = F.transitions()[0]
            sage: F._latex_transition_label_(t)
            \left[0\right] \mid \left[1\right]
        """
        return (format_function(transition.word_in) + "\\mid "
                + format_function(transition.word_out))


    def intersection(self, other, only_accessible_components=True):
        """
        Returns a new transducer which accepts an input if it is accepted by
        both given finite state machines producing the same output.

        INPUT:

        - ``other`` -- a transducer

        - ``only_accessible_components`` -- If ``True`` (default), then
          the result is piped through ``accessible_components``. If no
          ``new_input_alphabet`` is given, it is determined by
          :meth:`.determine_alphabets`.

        OUTPUT:

        A new transducer which computes the intersection
        (see below) of the languages of ``self`` and ``other``.

        The set of states of the transducer is the cartesian product of the
        set of states of both given transducer. There is a transition `((A,
        B), (C, D), a, b)` in the new transducer if there are
        transitions `(A, C, a, b)` and `(B, D, a, b)` in the old transducers.

        EXAMPLES::

            sage: transducer1 = Transducer([('1', '2', 1, 0),
            ....:                           ('2', '2', 1, 0),
            ....:                           ('2', '2', 0, 1)],
            ....:                          initial_states=['1'],
            ....:                          final_states=['2'],
            ....:                          determine_alphabets=True)
            sage: transducer2 = Transducer([('A', 'A', 1, 0),
            ....:                           ('A', 'B', 0, 0),
            ....:                           ('B', 'B', 0, 1),
            ....:                           ('B', 'A', 1, 1)],
            ....:                          initial_states=['A'],
            ....:                          final_states=['B'],
            ....:                          determine_alphabets=True)
            sage: res = transducer1.intersection(transducer2)
            sage: res.transitions()
            [Transition from ('1', 'A') to ('2', 'A'): 1|0,
             Transition from ('2', 'A') to ('2', 'A'): 1|0]

        In general, transducers are not closed under intersection. But
        for transducer which do not have epsilon-transitions, the
        intersection is well defined (cf. [BaWo2012]_). However, in
        the next example the intersection of the two transducers is
        not well defined. The intersection of the languages consists
        of `(a^n, b^n c^n)`. This set is not recognizable by a
        *finite* transducer.

        ::

            sage: t1 = Transducer([(0, 0, 'a', 'b'),
            ....:                  (0, 1, None, 'c'),
            ....:                  (1, 1, None, 'c')],
            ....:                 initial_states=[0],
            ....:                 final_states=[0, 1],
            ....:                 determine_alphabets=True)
            sage: t2 = Transducer([('A', 'A', None, 'b'),
            ....:                  ('A', 'B', 'a', 'c'),
            ....:                  ('B', 'B', 'a', 'c')],
            ....:                 initial_states=['A'],
            ....:                 final_states=['A', 'B'],
            ....:                 determine_alphabets=True)
            sage: t2.intersection(t1)
            Traceback (most recent call last):
            ...
            ValueError: An epsilon-transition (with empty input or output)
            was found.

        REFERENCES:

        .. [BaWo2012] Javier Baliosian and Dina Wonsever, *Finite State
           Transducers*, chapter in *Handbook of Finite State Based Models and
           Applications*, edited by Jiacun Wang, Chapman and Hall/CRC, 2012.
        """
        if not is_Transducer(other):
            raise TypeError(
                "Only a transducer can be intersected with a transducer.")

        def function(transition1, transition2):
            if not transition1.word_in or not transition2.word_in \
                    or not transition1.word_out or not transition2.word_out:
                raise ValueError("An epsilon-transition "
                                 "(with empty input or output) was found.")
            if transition1.word_in == transition2.word_in \
                    and transition1.word_out == transition2.word_out:
                return (transition1.word_in, transition1.word_out)
            else:
                raise LookupError

        return self.product_FiniteStateMachine(
            other,
            function,
            only_accessible_components=only_accessible_components)


    def cartesian_product(self, other, only_accessible_components=True):
        """
        .. WARNING::

            The default output of this method is scheduled to change.
            This docstring describes the new default behaviour, which can
            already be achieved by setting
            ``FSMOldCodeTransducerCartesianProduct`` to ``False``.

        Return a new transducer which can simultaneously process an input with
        ``self`` and ``other`` where the output labels are pairs of the
        original output labels.

        INPUT:

        - ``other`` - a finite state machine

        - ``only_accessible_components`` -- If ``True`` (default), then
          the result is piped through ``accessible_components``. If no
          ``new_input_alphabet`` is given, it is determined by
          :meth:`.determine_alphabets`.

        OUTPUT:

        A transducer which can simultaneously process an input with ``self``
        and ``other``.

        The set of states of the new transducer is the cartesian product of
        the set of states of ``self`` and ``other``.

        Let `(A, B, a, b)` be a transition of self and `(C, D, c, d)` be
        a transition of other. Then there is a transition `((A, C), (B,
        D), a, (b, d))` in the new transducer if `a = c`.

        EXAMPLES:

        Originally a different output was constructed by
        ``Transducer.cartesian_product``. This output is now produced by
        ``Transducer.intersection``.

        ::

            sage: transducer1 = Transducer([('A', 'A', 0, 0),
            ....:                           ('A', 'A', 1, 1)],
            ....:                          initial_states=['A'],
            ....:                          final_states=['A'],
            ....:                          determine_alphabets=True)
            sage: transducer2 = Transducer([(0, 1, 0, ['b', 'c']),
            ....:                           (0, 0, 1, 'b'),
            ....:                           (1, 1, 0, 'a')],
            ....:                          initial_states=[0],
            ....:                          final_states=[1],
            ....:                          determine_alphabets=True)
            sage: result = transducer1.cartesian_product(transducer2)
            doctest:1: DeprecationWarning: The output of
            Transducer.cartesian_product will change.
            Please use Transducer.intersection for the original output.
            See http://trac.sagemath.org/16061 for details.
            sage: result
            Transducer with 0 states

        By setting ``FSMOldCodeTransducerCartesianProduct`` to ``False``
        the new desired output is produced.

        ::

            sage: sage.combinat.finite_state_machine.FSMOldCodeTransducerCartesianProduct = False
            sage: result = transducer1.cartesian_product(transducer2)
            sage: result
            Transducer with 2 states
            sage: result.transitions()
            [Transition from ('A', 0) to ('A', 1): 0|(0, 'b'),(None, 'c'),
             Transition from ('A', 0) to ('A', 0): 1|(1, 'b'),
             Transition from ('A', 1) to ('A', 1): 0|(0, 'a')]
            sage: result([1, 0, 0])
            [(1, 'b'), (0, 'b'), (None, 'c'),  (0, 'a')]
            sage: (transducer1([1, 0, 0]), transducer2([1, 0, 0]))
            ([1, 0, 0], ['b', 'b', 'c', 'a'])

        Also final output words are correctly processed::

            sage: transducer1.state('A').final_word_out = 2
            sage: result = transducer1.cartesian_product(transducer2)
            sage: result.final_states()[0].final_word_out
            [(2, None)]

        The following transducer counts the number of 11 blocks minus
        the number of 10 blocks over the alphabet ``[0, 1]``.

        ::

            sage: count_11 = transducers.CountSubblockOccurrences(
            ....:     [1, 1],
            ....:     input_alphabet=[0, 1])
            sage: count_10 = transducers.CountSubblockOccurrences(
            ....:     [1, 0],
            ....:     input_alphabet=[0, 1])
            sage: count_11x10 = count_11.cartesian_product(count_10)
            sage: difference = transducers.sub([0, 1])(count_11x10)
            sage: T = difference.simplification().relabeled()
            sage: T.initial_states()
            [1]
            sage: sorted(T.transitions())
            [Transition from 0 to 1: 0|-1,
             Transition from 0 to 0: 1|1,
             Transition from 1 to 1: 0|0,
             Transition from 1 to 0: 1|0]
            sage: input =  [0, 1, 1,  0, 1,  0, 0, 0, 1, 1, 1,  0]
            sage: output = [0, 0, 1, -1, 0, -1, 0, 0, 0, 1, 1, -1]
            sage: T(input) == output
            True

        If ``other`` is an automaton, then :meth:`.cartesian_product` returns
        ``self`` where the input is restricted to the input accepted by
        ``other``.

        For example, if the transducer transforms the standard
        binary expansion into the non-adjacent form and the automaton
        recognizes the binary expansion without adjacent ones, then the
        cartesian product of these two is a transducer which does not change
        the input (except for changing ``a`` to ``(a, None)`` and ignoring a
        leading `0`).

        ::

            sage: NAF = Transducer([(0, 1, 0, None),
            ....:                   (0, 2, 1, None),
            ....:                   (1, 1, 0, 0),
            ....:                   (1, 2, 1, 0),
            ....:                   (2, 1, 0, 1),
            ....:                   (2, 3, 1, -1),
            ....:                   (3, 2, 0, 0),
            ....:                   (3, 3, 1, 0)],
            ....:                  initial_states=[0],
            ....:                  final_states=[1],
            ....:                  determine_alphabets=True)
            sage: aut11 = Automaton([(0, 0, 0), (0, 1, 1), (1, 0, 0)],
            ....:                   initial_states=[0],
            ....:                   final_states=[0, 1],
            ....:                   determine_alphabets=True)
            sage: res = NAF.cartesian_product(aut11)
            sage: res([1, 0, 0, 1, 0, 1, 0])
            [(1, None), (0, None), (0, None), (1, None), (0, None), (1, None)]

        This is obvious because if the standard binary expansion does not have
        adjacent ones, then it is the same as the non-adjacent form.

        Be aware that :meth:`.cartesian_product` is not commutative.

        ::

            sage: aut11.cartesian_product(NAF)
            Traceback (most recent call last):
            ...
            TypeError: Only an automaton can be intersected with an automaton.
        """
        if FSMOldCodeTransducerCartesianProduct:
            from sage.misc.superseded import deprecation
            deprecation(16061, "The output of Transducer.cartesian_product "
                               "will change. Please use "
                               "Transducer.intersection for the original "
                               "output.")
            return self.intersection(
                other,
                only_accessible_components=only_accessible_components)

        def function(transition1, transition2):
            if transition1.word_in == transition2.word_in:
                return (transition1.word_in,
                        list(itertools.izip_longest(transition1.word_out,
                                                    transition2.word_out)))
            else:
                raise LookupError

        def final_function(s1, s2):
            return list(itertools.izip_longest(s1.final_word_out,
                                               s2.final_word_out))

        return self.product_FiniteStateMachine(
            other,
            function,
            final_function=final_function,
            only_accessible_components=only_accessible_components)


    def simplification(self):
        """
        Returns a simplified transducer.

        INPUT:

        Nothing.

        OUTPUT:

        A new transducer.

        This function simplifies a transducer by Moore's algorithm,
        first moving common output labels of transitions leaving a
        state to output labels of transitions entering the state
        (cf. :meth:`.prepone_output`).

        The resulting transducer implements the same function as the
        original transducer.

        EXAMPLES::

            sage: fsm = Transducer([("A", "B", 0, 1), ("A", "B", 1, 0),
            ....:                           ("B", "C", 0, 0), ("B", "C", 1, 1),
            ....:                           ("C", "D", 0, 1), ("C", "D", 1, 0),
            ....:                           ("D", "A", 0, 0), ("D", "A", 1, 1)])
            sage: fsms = fsm.simplification()
            sage: fsms
            Transducer with 2 states
            sage: fsms.transitions()
            [Transition from ('A', 'C')
                          to ('B', 'D'): 0|1,
             Transition from ('A', 'C')
                          to ('B', 'D'): 1|0,
             Transition from ('B', 'D')
                          to ('A', 'C'): 0|0,
             Transition from ('B', 'D')
                          to ('A', 'C'): 1|1]
            sage: fsms.relabeled().transitions()
            [Transition from 0 to 1: 0|1,
             Transition from 0 to 1: 1|0,
             Transition from 1 to 0: 0|0,
             Transition from 1 to 0: 1|1]

        ::

            sage: fsm = Transducer([("A", "A", 0, 0),
            ....:                   ("A", "B", 1, 1),
            ....:                   ("A", "C", 1, -1),
            ....:                   ("B", "A", 2, 0),
            ....:                   ("C", "A", 2, 0)])
            sage: fsm_simplified = fsm.simplification()
            sage: fsm_simplified
            Transducer with 2 states
            sage: fsm_simplified.transitions()
            [Transition from ('A',) to ('A',): 0|0,
             Transition from ('A',) to ('B', 'C'): 1|1,0,
             Transition from ('A',) to ('B', 'C'): 1|-1,0,
             Transition from ('B', 'C') to ('A',): 2|-]

        ::

            sage: from sage.combinat.finite_state_machine import duplicate_transition_add_input
            sage: T = Transducer([('A', 'A', 1/2, 0),
            ....:                 ('A', 'B', 1/4, 1),
            ....:                 ('A', 'C', 1/4, 1),
            ....:                 ('B', 'A', 1, 0),
            ....:                 ('C', 'A', 1, 0)],
            ....:                initial_states=[0],
            ....:                final_states=['A', 'B', 'C'],
            ....:                on_duplicate_transition=duplicate_transition_add_input)
            sage: sorted(T.simplification().transitions())
            [Transition from ('A',) to ('A',): 1/2|0,
             Transition from ('A',) to ('B', 'C'): 1/2|1,
             Transition from ('B', 'C') to ('A',): 1|0]

        Illustrating the use of colors in order to avoid identification of states::

            sage: T = Transducer( [[0,0,0,0], [0,1,1,1],
            ....:                  [1,0,0,0], [1,1,1,1]],
            ....:                 initial_states=[0],
            ....:                 final_states=[0,1])
            sage: sorted(T.simplification().transitions())
            [Transition from (0, 1) to (0, 1): 0|0,
             Transition from (0, 1) to (0, 1): 1|1]
            sage: T.state(0).color = 0
            sage: T.state(0).color = 1
            sage: sorted(T.simplification().transitions())
            [Transition from (0,) to (0,): 0|0,
             Transition from (0,) to (1,): 1|1,
             Transition from (1,) to (0,): 0|0,
             Transition from (1,) to (1,): 1|1]
        """
        fsm = deepcopy(self)
        fsm.prepone_output()
        return fsm.quotient(fsm.equivalence_classes())


    def process(self, *args, **kwargs):
        """
        .. WARNING::

            The default output of this method is scheduled to change.
            This docstring describes the new default behaviour, which can
            already be achieved by setting
            ``FSMOldProcessOutput`` to ``False``.

        Returns whether the transducer accepts the input, the state
        where the computation stops and which output is generated.

        INPUT:

        - ``input_tape`` -- The input tape can be a list with entries from
          the input alphabet.

        - ``initial_state`` -- (default: ``None``) The state in which
          to start. If this parameter is ``None`` and there is only
          one initial state in the machine, then this state is taken.

        - ``full_output`` -- (default: ``True``) If set, then the full
          output is given, otherwise only the generated output (the
          third entry below only). If the input is not accepted, a
          ``ValueError`` is raised.

        OUTPUT:

        The full output is a triple, where

        - the first entry is ``True`` if the input string is accepted,

        - the second gives the reached state after processing the
          input tape (This is a state with label ``None`` if the input
          could not be processed, i.e., when at one point no
          transition to go could be found.), and

        - the third gives a list of the output labels used during
          processing.

        Note that in the case the transducer is not
        deterministic, one possible path is gone. This means, that in
        this case the output can be wrong.

        By setting ``FSMOldProcessOutput`` to ``False``
        the new desired output is produced.

        EXAMPLES::

            sage: sage.combinat.finite_state_machine.FSMOldProcessOutput = False  # activate new output behavior
            sage: from sage.combinat.finite_state_machine import FSMState
            sage: A = FSMState('A', is_initial = True, is_final = True)
            sage: binary_inverter = Transducer({A:[(A, 0, 1), (A, 1, 0)]})
            sage: binary_inverter.process([0, 1, 0, 0, 1, 1])
            (True, 'A', [1, 0, 1, 1, 0, 0])

        If we are only interested in the output, we can also use::

            sage: binary_inverter([0, 1, 0, 0, 1, 1])
            [1, 0, 1, 1, 0, 0]

        The following transducer transforms `0^n 1` to `1^n 2`::

            sage: T = Transducer([(0, 0, 0, 1), (0, 1, 1, 2)])
            sage: T.state(0).is_initial = True
            sage: T.state(1).is_final = True

        We can see the different possibilites of the output by::

            sage: [T.process(w) for w in [[1], [0, 1], [0, 0, 1], [0, 1, 1],
            ....:                         [0], [0, 0], [2, 0], [0, 1, 2]]]
            [(True, 1, [2]), (True, 1, [1, 2]),
             (True, 1, [1, 1, 2]), (False, None, None),
             (False, 0, [1]), (False, 0, [1, 1]),
             (False, None, None), (False, None, None)]

        If we just want a condensed output, we use::

            sage: [T.process(w, full_output=False)
            ....:      for w in [[1], [0, 1], [0, 0, 1]]]
            [[2], [1, 2], [1, 1, 2]]
            sage: T.process([0, 1, 2], full_output=False)
            Traceback (most recent call last):
            ...
            ValueError: Invalid input sequence.

        It is equivalent to::

            sage: [T(w) for w in [[1], [0, 1], [0, 0, 1]]]
            [[2], [1, 2], [1, 1, 2]]
            sage: T([0, 1, 2])
            Traceback (most recent call last):
            ...
            ValueError: Invalid input sequence.
        """
        if FSMOldProcessOutput:
            from sage.misc.superseded import deprecation
            deprecation(16132, "The output of Transducer.process "
                               "(and thus of Transducer.__call__) "
                               "will change. Please use the corresponding "
                               "functions from FiniteStateMachine "
                               "for the original output.")
            return super(Transducer, self).process(*args, **kwargs)

        if not kwargs.has_key('full_output'):
            kwargs['full_output'] = True

        it = self.iter_process(*args, **kwargs)
        for _ in it:
            pass

        # process output
        if kwargs['full_output']:
            if it.current_state.label() is None:
                return (it.accept_input, it.current_state, None)
            else:
                return (it.accept_input, it.current_state, it.output_tape)
        else:
            if not it.accept_input:
                raise ValueError("Invalid input sequence.")
            return it.output_tape



#*****************************************************************************


def is_FSMProcessIterator(PI):
    """
    Tests whether or not ``PI`` inherits from :class:`FSMProcessIterator`.

    TESTS::

        sage: from sage.combinat.finite_state_machine import is_FSMProcessIterator, FSMProcessIterator
        sage: is_FSMProcessIterator(FSMProcessIterator(FiniteStateMachine([[0, 0, 0, 0]], initial_states=[0])))
        True
    """
    return isinstance(PI, FSMProcessIterator)


class FSMProcessIterator(SageObject):
    """
    This class is for processing an input string on a finite state
    machine.

    An instance of this class is generated when
    :meth:`FiniteStateMachine.process` or
    :meth:`FiniteStateMachine.iter_process` of the finite state
    machine is invoked. It behaves like an iterator which, in each
    step, takes one letter of the input and runs (one step on) the
    finite state machine with this input. More precisely, in each
    step, the process iterator takes an outgoing transition of the
    current state, whose input label equals the input letter of the
    tape. The output label of the transition, if present, is written
    on the output tape.

    INPUT:

    - ``fsm`` -- The finite state machine on which the input should be
      processed.

    - ``input_tape`` -- The input tape. It can be anything that is
      iterable.

    - ``initial_state`` -- The initial state in which the machine
      starts. If this is ``None``, the unique inital state of the finite
      state machine is takes. If there are several, a ``ValueError`` is
      raised.

    The process (iteration) stops if there are no more input letters
    on the tape. In this case a StopIteration exception is thrown. As
    result the following attributes are available:

    - ``accept_input`` -- Is ``True`` if the reached state is a final state.

    - ``current_state`` -- The current/reached state in the process.

    - ``output_tape`` -- The written output.

    Current values of those attributes (except ``accept_input``) are
    (also) available during the iteration.

    OUTPUT:

    An iterator.

    EXAMPLES:

    The following transducer reads binary words and outputs a word,
    where blocks of ones are replaced by just a single one. Further
    only words that end with a zero are accepted.

    ::

        sage: T = Transducer({'A': [('A', 0, 0), ('B', 1, None)],
        ....:                 'B': [('B', 1, None), ('A', 0, [1, 0])]},
        ....:     initial_states=['A'], final_states=['A'])
        sage: input = [1, 1, 0, 0, 1, 0, 1, 1, 1, 0]
        sage: T.process(input)
        (True, 'A', [1, 0, 0, 1, 0, 1, 0])

    The function :meth:`FiniteStateMachine.process` created a new
    ``FSMProcessIterator``. We can do that manually, too, and get full
    access to the iteration process::

        sage: from sage.combinat.finite_state_machine import FSMProcessIterator
        sage: it = FSMProcessIterator(T, input_tape=input)
        sage: for _ in it:
        ....:     print (it.current_state, it.output_tape)
        ('B', [])
        ('B', [])
        ('A', [1, 0])
        ('A', [1, 0, 0])
        ('B', [1, 0, 0])
        ('A', [1, 0, 0, 1, 0])
        ('B', [1, 0, 0, 1, 0])
        ('B', [1, 0, 0, 1, 0])
        ('B', [1, 0, 0, 1, 0])
        ('A', [1, 0, 0, 1, 0, 1, 0])
        sage: it.accept_input
        True

    TESTS::

        sage: T = Transducer([[0, 0, 0, 0]])
        sage: T.process([])
        Traceback (most recent call last):
        ...
        ValueError: No state is initial.

    ::

        sage: T = Transducer([[0, 1, 0, 0]], initial_states=[0, 1])
        sage: T.process([])
        Traceback (most recent call last):
        ...
        ValueError: Several initial states.

    """
    def __init__(self, fsm, input_tape=None, initial_state=None, **kwargs):
        """
        See :class:`FSMProcessIterator` for more information.

        EXAMPLES::

            sage: from sage.combinat.finite_state_machine import FSMProcessIterator
            sage: inverter = Transducer({'A': [('A', 0, 1), ('A', 1, 0)]},
            ....:     initial_states=['A'], final_states=['A'])
            sage: it = FSMProcessIterator(inverter, input_tape=[0, 1])
            sage: for _ in it:
            ....:     pass
            sage: it.output_tape
            [1, 0]
        """
        self.fsm = fsm
        if initial_state is None:
            fsm_initial_states = self.fsm.initial_states()
            try:
                self.current_state = fsm_initial_states[0]
            except IndexError:
                raise ValueError("No state is initial.")
            if len(fsm_initial_states) > 1:
                raise ValueError("Several initial states.")
        else:
            self.current_state = initial_state
        self.output_tape = []
        if input_tape is None:
            self._input_tape_iter_ = iter([])
        else:
            if hasattr(input_tape, '__iter__'):
                self._input_tape_iter_ = iter(input_tape)
            else:
                raise ValueError("Given input tape is not iterable.")

    def __iter__(self):
        """
        Returns ``self``.

        TESTS::

            sage: from sage.combinat.finite_state_machine import FSMProcessIterator
            sage: inverter = Transducer({'A': [('A', 0, 1), ('A', 1, 0)]},
            ....:     initial_states=['A'], final_states=['A'])
            sage: it = FSMProcessIterator(inverter, input_tape=[0, 1])
            sage: id(it) == id(iter(it))
            True
        """
        return self

    def next(self):
        """
        Makes one step in processing the input tape.

        INPUT:

        Nothing.

        OUTPUT:

        It returns the taken transition. A ``StopIteration`` exception is
        thrown when there is nothing more to read.

        EXAMPLES::

            sage: from sage.combinat.finite_state_machine import FSMProcessIterator
            sage: inverter = Transducer({'A': [('A', 0, 1), ('A', 1, 0)]},
            ....:     initial_states=['A'], final_states=['A'])
            sage: it = FSMProcessIterator(inverter, input_tape=[0, 1])
            sage: it.next()
            Transition from 'A' to 'A': 0|1
            sage: it.next()
            Transition from 'A' to 'A': 1|0
            sage: it.next()
            Traceback (most recent call last):
            ...
            StopIteration

        TESTS::

            sage: Z = Transducer()
            sage: s = Z.add_state(0)
            sage: s.is_initial = True
            sage: s.is_final = True
            sage: s.final_word_out = [1, 2]
            sage: Z.process([])
            (True, 0, [1, 2])
        """
        if hasattr(self, 'accept_input'):
            raise StopIteration
        try:
            # process current state
            transition = None
            try:
                transition = self.current_state.hook(
                    self.current_state, self)
            except AttributeError:
                pass
            self.write_word(self.current_state.word_out)

            # get next
            if not isinstance(transition, FSMTransition):
                next_word = []
                found = False

                try:
                    while not found:
                        next_word.append(self.read_letter())
                        try:
                            transition = self.get_next_transition(
                                next_word)
                            found = True
                        except ValueError:
                            pass
                except StopIteration:
                    # this means input tape is finished
                    if len(next_word) > 0:
                        self.current_state = FSMState(None,
                                                      allow_label_None=True)
                    raise StopIteration

            # process transition
            try:
                transition.hook(transition, self)
            except AttributeError:
                pass
            self.write_word(transition.word_out)

            # go to next state
            self.current_state = transition.to_state

        except StopIteration:
            # this means, either input tape is finished or
            # someone has thrown StopIteration manually (in one
            # of the hooks)
            if self.current_state.label is None or not self.current_state.is_final:
                self.accept_input = False
            if not hasattr(self, 'accept_input'):
                self.accept_input = True
            if self.current_state.is_final:
                self.write_word(self.current_state.final_word_out)
            raise StopIteration

        # return
        return transition

    def read_letter(self):
        """
        Reads a letter from the input tape.

        INPUT:

        Nothing.

        OUTPUT:

        A letter.

        Exception ``StopIteration`` is thrown if tape has reached
        the end.

        EXAMPLES::

            sage: from sage.combinat.finite_state_machine import FSMProcessIterator
            sage: inverter = Transducer({'A': [('A', 0, 1), ('A', 1, 0)]},
            ....:     initial_states=['A'], final_states=['A'])
            sage: it = FSMProcessIterator(inverter, input_tape=[0, 1])
            sage: it.read_letter()
            0
        """
        return self._input_tape_iter_.next()

    def write_letter(self, letter):
        """
        Writes a letter on the output tape.

        INPUT:

        - ``letter`` -- the letter to be written.

        OUTPUT:

        Nothing.

        EXAMPLES::

            sage: from sage.combinat.finite_state_machine import FSMProcessIterator
            sage: inverter = Transducer({'A': [('A', 0, 1), ('A', 1, 0)]},
            ....:     initial_states=['A'], final_states=['A'])
            sage: it = FSMProcessIterator(inverter, input_tape=[0, 1])
            sage: it.write_letter(42)
            sage: it.output_tape
            [42]
        """
        self.output_tape.append(letter)

    def write_word(self, word):
        """
        Writes a word on the output tape.

        INPUT:

        - ``word`` -- the word to be written.

        OUTPUT:

        Nothing.

        EXAMPLES::

            sage: from sage.combinat.finite_state_machine import FSMProcessIterator
            sage: inverter = Transducer({'A': [('A', 0, 1), ('A', 1, 0)]},
            ....:     initial_states=['A'], final_states=['A'])
            sage: it = FSMProcessIterator(inverter, input_tape=[0, 1])
            sage: it.write_word([4, 2])
            sage: it.output_tape
            [4, 2]
        """
        for letter in word:
            self.write_letter(letter)

    def get_next_transition(self, word_in):
        """
        Returns the next transition according to ``word_in``. It is
        assumed that we are in state ``self.current_state``.

        INPUT:

        - ``word_in`` -- the input word.

        OUTPUT:

        The next transition according to ``word_in``. It is assumed
        that we are in state ``self.current_state``. If no transition
        matches, a ``ValueError`` is thrown.

        EXAMPLES::

            sage: from sage.combinat.finite_state_machine import FSMProcessIterator
            sage: inverter = Transducer({'A': [('A', 0, 1), ('A', 1, 0)]},
            ....:     initial_states=['A'], final_states=['A'])
            sage: it = FSMProcessIterator(inverter, input_tape=[0, 1])
            sage: it.get_next_transition([0])
            Transition from 'A' to 'A': 0|1
            sage: it.get_next_transition([2])
            Traceback (most recent call last):
            ...
            ValueError: No transition with input [2] found.
        """
        for transition in self.current_state.transitions:
            if transition.word_in == word_in:
                return transition
        raise ValueError("No transition with input %s found." % (word_in,))


#*****************************************************************************


def setup_latex_preamble():
    """
    This function adds the package ``tikz`` with support for automata
    to the preamble of Latex so that the finite state machines can be
    drawn nicely.

    INPUT:

    Nothing.

    OUTPUT:

    Nothing.

    In the Sage notebook, you probably want to use
    ``latex.mathjax_avoid_list('tikzpicture')`` such that
    :func:`~sage.misc.latex.view` actually shows the result.
    See the section on :ref:`finite_state_machine_LaTeX_output`
    in the introductory examples of this module.

    TESTS::

        sage: from sage.combinat.finite_state_machine import setup_latex_preamble
        sage: setup_latex_preamble()
        sage: latex.mathjax_avoid_list('tikzpicture')
    """
    latex.add_package_to_preamble_if_available('tikz')
    latex.add_to_preamble('\\usetikzlibrary{automata}')


#*****************************************************************************
