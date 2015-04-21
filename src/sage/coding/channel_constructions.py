r"""
Channels

Given an input space and an output space, a channel takes element from
the input space (the message) and transforms it into an element of the output space
(the transmitted message).

This file contains the following elements:

    - *AbstractChannel*, the abstract class for Channels
    - *ChannelStaticErrorRate*, which creates a specific number of errors in each
      transmitted message
    - *ChannelErrorErasure*, which creates a specific number of errors and a
      specific number of erasures in each transmitted message
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

from sage.rings.finite_rings.constructor import GF
from sage.misc.prandom import randint, random
from sage.modules.free_module_element import vector
from sage.misc.abstract_method import abstract_method
from sage.combinat.cartesian_product import CartesianProduct
from sage.modules.free_module import VectorSpace

def _random_error_position(n , number_errors):
    r"""
    Returns a list of exactly ``number_errors`` random numbers between 0 and ``n-1``
    This is a helper function, for internal use only.
    This function was taken from codinglib (https://bitbucket.org/jsrn/codinglib/)
    and was written by Johan Nielsen.

    INPUT:

    - ``number_errors`` -- the number of elements in the list

    - ``n`` -- upper bound for the elements of the list

    OUTPUT:

    - A list of integers

    EXAMPLES::

        sage: sage.coding.channel_constructions._random_error_position(6, 2) # random
        [1, 4]
    """
    error_position = []
    i = 0
    while i < n and number_errors > 0:
        if random() < number_errors/(n-i):
            error_position.append(i)
            number_errors -= 1
        i += 1
    return error_position

def _random_error_vector(n, F, error_positions):
    r"""
    Return a vector of length ``n`` over ``F`` filled with random non-zero coefficients
    at the positions given by ``error_positions``.
    This is a helper function, for internal use only.
    This function was taken from codinglib (https://bitbucket.org/jsrn/codinglib/)
    and was written by Johan Nielsen.

    INPUT:

    - ``n`` -- the length of the vector

    - ``F`` -- the field over which the vector is defined

    - ``error_positions`` -- the non-zero positions of the vector

    OUTPUT:

    - a vector of ``F``

    EXAMPLES::

        sage: sage.coding.channel_constructions._random_error_vector(5, GF(2), [1,3])
        (0, 1, 0, 1, 0)
    """
    vect = [F.zero()]*n
    for i in error_positions:
        while vect[i].is_zero():
            vect[i] = F.random_element()
    return vector(F, vect)

def _tuple_to_integer(value):
    r"""
    Returns an integer from ``value``. If ``value`` is a tuple, it will return a random
    integer between its bounds.
    This is a helper function, for internal use only.

    INPUT:

    - ``value`` -- an integer or a couple of integers

    OUTPUT:

    - an integer

    EXAMPLES::

        sage: sage.coding.channel_constructions._tuple_to_integer(4)
        4

        sage: sage.coding.channel_constructions._tuple_to_integer((1,5)) # random
        3
    """
    value = value if not hasattr(value, "__iter__") else randint(value[0], value[1])
    return value





class AbstractChannel(object):
    r"""
    Abstract top-class for Channel objects.

    All channel objects must heritate from this class. To implement a channel subclass, one should
    do the following:

    - heritate from this class,

    - call the super constructor,

    - override :meth:`transmit_unsafe`.

    While not being mandatory, it might be useful to reimplement representation methods (``__repr__`` and
    ``_latex_``), along with a comparison method (``__eq__``).

    This abstract class provides the following parameters:

        - ``input_space`` -- the space of the words to transmit

        - ``output_space`` -- the space of the transmitted words
    """

    def __init__(self, input_space, output_space):
        r"""
        Initializes parameters for a Channel object.

        This is a private method, which should be called by the constructor
        of every encoder, as it automatically initializes the mandatory
        parameters of a Channel object.

        INPUT:

        - ``input_space`` -- the space of the words to transmit

        - ``output_space`` -- the space of the transmitted words

        EXAMPLES:

        We first create a new Channel subclass::

            sage: class ChannelExample(sage.coding.channel_constructions.AbstractChannel):
            ....:   def __init__(self, input_space, output_space):
            ....:       super(ChannelExample, self).__init__(input_space, output_space)

        We now create a member of our newly made class::

            sage: input = VectorSpace(GF(7), 6)
            sage: output = VectorSpace(GF(7), 5)
            sage: Chan = ChannelExample(input, output)

        We can check its parameters::

            sage: Chan.input_space()
            Vector space of dimension 6 over Finite Field of size 7
            sage: Chan.output_space()
            Vector space of dimension 5 over Finite Field of size 7
        """
        self._input_space = input_space
        self._output_space = output_space

    def transmit(self, message):
        r"""
        Returns ``message``, modified accordingly with the algorithm of the channel it was
        transmitted through.
        Checks if ``message`` belongs to the input space, and returns an exception if not.

        INPUT:

        - ``message`` -- a vector

        OUTPUT:

        - a vector of the output space of ``self``

        EXAMPLES::

            sage: F = VectorSpace(GF(59), 6)
            sage: n_err = 2
            sage: Chan = channels.ChannelStaticErrorRate(F, n_err)
            sage: msg = F((4, 8, 15, 16, 23, 42))
            sage: Chan.transmit(msg) # random
            (4, 14, 15, 16, 17, 42)

        If we transmit a vector which is not in the input space of ``self``::

            sage: F = VectorSpace(GF(59), 6)
            sage: n_err = 2
            sage: Chan = channels.ChannelStaticErrorRate(F, n_err)
            sage: msg = (4, 8, 15, 16, 23, 42)
            sage: Chan.transmit(msg)
            Traceback (most recent call last):
            ...
            TypeError: Message must be an element of the input space for the given channel
        """
        if message in self.input_space():
            return self.transmit_unsafe(message)
        else :
            raise TypeError("Message must be an element of the input space for the given channel")

    def input_space(self):
        r"""
        Returns the input space of ``self``.

        EXAMPLES::

            sage: F = VectorSpace(GF(59), 6)
            sage: n_err = 2
            sage: Chan = channels.ChannelStaticErrorRate(F, n_err)
            sage: Chan.input_space()
            Vector space of dimension 6 over Finite Field of size 59

        """
        return self._input_space

    def output_space(self):
        r"""
        Returns the output space of ``self``.

        EXAMPLES::

            sage: F = VectorSpace(GF(59), 6)
            sage: n_err = 2
            sage: Chan = channels.ChannelStaticErrorRate(F, n_err)
            sage: Chan.output_space()
            Vector space of dimension 6 over Finite Field of size 59
        """
        return self._output_space
    
    @abstract_method
    def transmit_unsafe(self, message):
        r"""
        Returns ``message``, modified accordingly with the algorithm of the channel it was
        transmitted through.
        This method does not check if ``message`` belongs to the input space of``self``.

        This is an abstract method which should be reimplemented in all the subclasses of
        Channel.
        """
        raise NotImplementedError

