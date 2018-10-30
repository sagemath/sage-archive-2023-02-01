r"""
Groups of imaginary elements

.. NOTE::

    One main purpose of such groups is in an
    :doc:`asymptotic ring's <asymptotic_ring>`
    :doc:`growth group <growth_group>` when an element like `n^z`
    (for some constant `z`) is split into
    `n^{\Re z + I\Im z}`.
    (Note that the first summand in the exponent determines the growth,
    the second does not influence the growth.)

AUTHORS:

- Daniel Krenn (2018)

Classes and Methods
===================
"""
#*****************************************************************************
# Copyright (C) 2018 Daniel Krenn <dev@danielkrenn.at>
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 2 of the License, or
# (at your option) any later version.
#                  http://www.gnu.org/licenses/
#*****************************************************************************
from __future__ import absolute_import

