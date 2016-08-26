r"""
Verbosity System and Logging in SageMath

Howto: Logging
==============

Using Python's Logging Module
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Import it::

    sage: import logging
    sage: logging.basicConfig()  # only needed once

Setting the level::

    sage: logging.getLogger().setLevel(logging.INFO)

Log something::

    sage: logger = logging.getLogger(__name__)
    sage: logger.info('Hello. I am talking to you.')
    INFO:__main__:Hello. I am talking to you.

If we haven't set the logging level to ``logging.INFO``, then the previous
wouldn't have been shown.
::

    sage: logger.debug('Hello. I am really talking a lot.')

The latter is not shown as the current logging level is only
``logging.INFO`` and not ``logging.DEBUG``.

Reset the level::

    sage: logging.getLogger().setLevel(logging.WARNING)

Warnings are still shown at this default level (``logging.WARNING``)::

    sage: logger.warn('Hello. I am warning you.')
    WARNING:__main__:Hello. I am warning you.

And that's all.

There are a lot more features, see
:python:`Logging facility for Python<library/logging.html>`.


Using SageMath's Verbosity System
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Alternatively, this module provides
:func:`verbose`, :func:`set_verbose`, :func:`get_verbose` which can
be used as follows::

    sage: from sage.misc.verbose import verbose, set_verbose, get_verbose
    sage: set_verbose(1)
    sage: t = verbose("This is SageMath.", level=0)
    verbose 0 (<module>) This is SageMath.
    sage: t = verbose("This is SageMath.", level=1)
    verbose 1 (<module>) This is SageMath.
    sage: t = verbose("This is SageMath.", level=2)


Logging Levels of SageMath and Python
=====================================

.. csv-table::
    :class: contentstable
    :widths: 20, 20
    :delim: |

    SageMath | Python
    `-2` | ``logging.CRITICAL``
    `-1` | ``logging.ERROR``
    `0` | ``logging.WARNING``
    `1` | ``logging.INFO``
    `2` | ``logging.DEBUG``


Various
=======

AUTHORS:

- Daniel Krenn (2016)


Functions
=========
"""
#*****************************************************************************
#       Copyright (C) 2016 Daniel Krenn <dev@danielkrenn.at>
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#                  http://www.gnu.org/licenses/
#*****************************************************************************
