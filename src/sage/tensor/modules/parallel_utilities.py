r"""
Parallelization utilities

This module defines a singleton class containing parallelization
informations when using tensor algebra.


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


class TensorParallelism(Singleton, SageObject):
    r"""
    Class governing Tensors Parallelism.

    Contains information about parallelism of Tensor algebra:
    
    * ``nproc`` -- (default: 1) number of processors to use in computations.   
    * ``use_paral`` -- (default: ``False``) flag, ``True`` if parallelism is used.
    
    """
    def __init__(self): 
        r"""
        Intializes the parallelism. Only an instance (Singleton) of this
        class is created.
        Default : 1 processor is used, so no parallelism.

        TEST ::

            sage: TP = TensorParallelism()
            sage: print TP
            Number of cpu used = 1

        Test of the singleton character::

            sage: TensorParallelism() is TP
            True
            
        """
        self.nproc = 1
        self.use_paral = False

    def _repr_(self):
        r"""
        String to print the number of cores.

        TEST ::

            sage: TensorParallelism()._repr_()
            'Number of cpu used = 1'

        """
        return "Number of cpu used = "+str(self.nproc)

    def set(self,nproc=None):
        r"""
        Changes the status of the parallelism in SageManifolds.
        If not argument is given, the number of processors will set
        to the maximum of cores found.

        EXAMPLE ::

            sage: print TensorParallelism()
            Number of cpu used = 1
            sage: TensorParallelism().set(4)
            sage: print TensorParallelism()
            Number of cpu used = 4

        """
        self.nproc = ncpus() if nproc is None else nproc
        self.use_paral = True if self.nproc!=1 else False

