r"""
Parallelization utilities

This module defines

- the singleton class :class:`Parallelism` to gather the information
  relative to the parallelization  (basically the number of  processes to be used)
- the global functions :func:`set_nproc_tensor` and :func:`get_nproc_tensor` to be used in
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
from sage.rings.integer import Integer

class Parallelism(Singleton, SageObject):
    r"""
    Singleton class for managing the number of processes used in parallel
    computations involved in tensor algebra.

    EXAMPLES::

        sage: from sage.parallel.parallelism import Parallelism
        sage: Parallelism().get('tensor')
        1
        sage: Parallelism().set('tensor',4)
        sage: Parallelism().get('tensor')
        4
        sage: Parallelism().set('tensor',1)
        sage: Parallelism().get('tensor')
        1

    """
    def __init__(self):
        r"""
        Only a single instance of this class is created (singleton model)

        TEST::

            sage: from sage.parallel.parallelism import Parallelism
            sage: TP = Parallelism()
            sage: TP
            Number of cpu used = 1

        Test of the singleton character::

            sage: Parallelism() is TP
            True

        The test suite is passed::

            sage: TestSuite(TP).run()

        """
        self._default = None
        self._nproc = {'tensor' : 1}


    def _repr_(self):
        r"""
        String representation of the object.

        TEST::

            sage: from sage.parallel.parallelism import Parallelism
            sage: Parallelism()._repr_()
            'Number of cpu used = 1'

        """
        return str(self._nproc)

    def set(self, field, nproc=None):
        r"""
        Set the number of processes to be launched in parallel
        computations regarding some specific field in sage

        INPUT:
        -

        - ``field`` -- specify the part of sage for which set the
          number of processes for parallel computations. 

        - ``nproc`` -- (defaut: ``None``) number of processes; if ``None``, the
          number of processes will be set to the maximum of cores found on the
          computer.

        EXAMPLES:

        The default is a single processor (no parallelization)::

            sage: from sage.parallel.parallelism import Parallelism
            sage: Parallelism()
            Number of cpu used = 1

        Asking for parallelization on 4 cores in tensor algebra::

            sage: Parallelism().set('tensor',4)
            sage: Parallelism()
            Number of cpu used = 4

        Using all the cores available on the computer::

            sage: Parallelism().set('tensor')
            sage: Parallelism()  # random
            Number of cpu used = 8

        Switching off the parallelization::

            sage: Parallelism().set('tensor',1)
            sage: Parallelism()
            Number of cpu used = 1

        """
        if field not in self._nproc :
            raise KeyError(
                """entry for field {0} is not implemented in the Parallelism""".format(field)
            )
        if nproc is not None :
            if not isinstance(nproc,(int,Integer)):
                raise TypeError("nproc must be integer")

        if nproc is None:
            if self._default is None:
                self._nproc[field] = ncpus()  
            else: 
                self._nproc[field] = self._default
        else:
            self._nproc[field] = nproc

    def get(self,field):
        r"""
        Get the number of processes which wiil be used in parallel
        computations regarding some specific field in sage

        INPUT:
        -

        - ``field`` -- specify the part of sage for which get the
          number of processes for parallel computations. 

        EXAMPLES:

            The default is a single processor (no parallelization)::

            sage: from sage.parallel.parallelism import Parallelism
            sage: Parallelism().get('tensor')
            1

            Asking for parallelization on 4 cores::
            
            sage: Parallelism().set('tensor',4)
            sage: Parallelism().get('tensor')
            4

        """
        if field not in self._nproc :
            raise KeyError("""entry for field {} does not correspond
                              is not implemented in the Parallelism""".format(field))

        return self._nproc[field]


    def get_all(self):
        r"""
        Get the number of processes which wiil be used in parallel
        computations in any field in sage

        INPUT:
        -

        EXAMPLES:

            sage: from sage.parallel.parallelism import Parallelism
            sage: Parallelism().get_all()
            {
             'tensor': 1
            }

            Asking for parallelization on 4 cores::
            
            sage: Parallelism().set('tensor',4)
            sage: Parallelism().get_all()
            {
             'tensor': 4
            }

        """
        return self._nproc


    def set_default(self, nproc=None):
        r"""
        Set the default number of processes to be launched in parallel
        computations for all fields in sage

        INPUT:
        -

        - ``nproc`` -- (defaut: ``None``) default number of processes; if ``None``, the
          number of processes will be set to the maximum of cores found on the
          computer.

        EXAMPLES:

        The default is a single processor (no parallelization)::

            sage: from sage.parallel.parallelism import Parallelism
            sage: Parallelism()

        Asking for parallelization on 4 cores for all fields::

            sage: Parallelism().set_default(4)
            sage: Parallelism()
            Number of cpu used = 4

        Using all the cores available on the computer::

            sage: Parallelism().set_default('tensor')
            sage: Parallelism()  # random
            Number of cpu used = 8

        """
        if nproc is not None :
            if not isinstance(nproc,(int,Integer)):
                raise TypeError("nproc must be integer")

        self._default = nproc

    def get_default(self):
        r"""
        Get the default number of processes to be launched in parallel
        computations for all fields in sage

        INPUT:
        -

        EXAMPLES:

        The default is a single processor (no parallelization)::

            sage: from sage.parallel.parallelism import Parallelism
            sage: Parallelism()
            sage: Parallelism().get_default('tensor')
            1


        """
        return self._default

