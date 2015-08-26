r"""
Parallelization utilities

This module defines

- the singleton class :class:`TensorParallelCompute` to gather the information
  relative to the parallelization of tensor algebra (basically the number of
  processes to be used)
- the global functions :func:`set_nproc` and :func:`get_nproc_tensor` to be used in
  a Sage session for managing the number of processes involved in the
  parallelization.

Some examples of parallelization are provided in the documentation
of :meth:`sage.tensor.modules.comp.Components.contract`.

AUTHORS:

- Marco Mancini, Eric Gourgoulhon, Michal Bejger (2015): initial version

"""

#******************************************************************************
#       Copyright (C) 2015 Marco Mancini <marco.mancini@obspm.fr>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#  as published by the Free Software Foundation; either version 2 of
#  the License, or (at your option) any later version.
#                  http://www.gnu.org/licenses/
#******************************************************************************

from sage.structure.sage_object import SageObject
from sage.misc.fast_methods import Singleton
from sage.parallel.ncpus import ncpus


class TensorParallelCompute(Singleton, SageObject):
    r"""
    Singleton class for managing the number of processes used in parallel
    computations involved in tensor algebra.

    EXAMPLES::

        sage: from sage.tensor.modules.parallel_utilities import TensorParallelCompute
        sage: TensorParallelCompute()
        Number of cpu used = 1
        sage: TensorParallelCompute().set(4)
        sage: TensorParallelCompute()
        Number of cpu used = 4
        sage: TensorParallelCompute().set(1)
        sage: TensorParallelCompute()
        Number of cpu used = 1

    """
    def __init__(self):
        r"""
        Only a single instance of this class is created (singleton model)

        TEST::

            sage: from sage.tensor.modules.parallel_utilities import TensorParallelCompute
            sage: TP = TensorParallelCompute()
            sage: TP
            Number of cpu used = 1

        Test of the singleton character::

            sage: TensorParallelCompute() is TP
            True

        The test suite is passed::

            sage: TestSuite(TP).run()

        """
        self._nproc = 1
        self._use_paral = False

    def _repr_(self):
        r"""
        String representation of the object.

        TEST::

            sage: from sage.tensor.modules.parallel_utilities import TensorParallelCompute
            sage: TensorParallelCompute()._repr_()
            'Number of cpu used = 1'

        """
        return "Number of cpu used = {}".format(self._nproc)

    def set(self, nproc=None):
        r"""
        Set the number of processes to be launched in parallel
        computations

        INPUT:

        - ``nproc`` -- (defaut: ``None``) number of processes; if ``None``, the
          number of processes will be set to the maximum of cores found on the
          computer.

        EXAMPLES:

        The default is a single processor (no parallelization)::

            sage: from sage.tensor.modules.parallel_utilities import TensorParallelCompute
            sage: TensorParallelCompute()
            Number of cpu used = 1

        Asking for parallelization on 4 cores::

            sage: TensorParallelCompute().set(4)
            sage: TensorParallelCompute()
            Number of cpu used = 4

        Using all the cores available on the computer::

            sage: TensorParallelCompute().set()
            sage: TensorParallelCompute()  # random
            Number of cpu used = 8

        Switching off the parallelization::

            sage: TensorParallelCompute().set(1)
            sage: TensorParallelCompute()
            Number of cpu used = 1

        """
        self._nproc = ncpus() if nproc is None else nproc
        self._use_paral = True if self._nproc!=1 else False


def set_nproc_tensor(nproc=None):
    r"""
    Set the number of processes for parallelizing computations
    relative to tensor algebra.

    INPUT:

    - ``nproc`` -- (defaut: ``None``) number of processes; if ``None``, the
      number of processes will be set to the maximum of cores found on the
      computer.

    EXAMPLES:

    The default is a single processor (no parallelization)::

        sage: get_nproc_tensor()
        1

    Asking for parallelization on 4 cores::

        sage: set_nproc_tensor(4)
        sage: get_nproc_tensor()
        4

    Using all the cores available on the computer::

        sage: set_nproc_tensor()
        sage: get_nproc_tensor()  # random
        8

    Switching off the parallelization::

        sage: set_nproc_tensor(1)
        sage: get_nproc_tensor()
        1

    See :meth:`sage.tensor.modules.comp.Components.contract` for a concrete
    example of use.

    """
    TensorParallelCompute().set(nproc)

def get_nproc_tensor():
    r"""
    Return the number of processes used in parallelized tensorial
    computations.

    EXAMPLES:

    The default is a single processor (no parallelization)::

        sage: get_nproc_tensor()
        1

    Asking for parallelization on 4 cores::

        sage: set_nproc_tensor(4)
        sage: get_nproc_tensor()
        4

    """
    return TensorParallelCompute()._nproc
