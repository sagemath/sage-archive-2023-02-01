r"""
Parallelization utilities

This module defines a class containing parallelization informations

AUTHORS:

- Eric Gourgoulhon, Michal Bejger, Marco Mancini (2014-2015): initial version

"""

#******************************************************************************
#       Copyright (C) 2015 Eric Gourgoulhon <eric.gourgoulhon@obspm.fr>
#       Copyright (C) 2015 Michal Bejger <bejger@camk.edu.pl>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#  as published by the Free Software Foundation; either version 2 of
#  the License, or (at your option) any later version.
#                  http://www.gnu.org/licenses/
#******************************************************************************

from sage.parallel.ncpus import ncpus


class ManifoldParallelism():
    r"""
    Class governing SageManifolds Parallelism.

    Contains information about parallelism of SageManifolds:
    nproc = Number of processors (nproc) to use.
    use_paral = Flag if parallelism is used.
    Default is nproc=1 and use_paral=False
    
    """
    def __init__(self,nproc=1):
        r"""
        Intializes the parallelism. Normally only an instance of this
        class is created in a global (visible) variable.
        Default : 1 processor is used, so no parallelism.
        """
        self.nproc = nproc
        self.use_paral = True if self.nproc!=1 else False

    def __str__(self):
        r"""
        String to print the number of cores.
        """
        return "Number of cpu used = "+str(self.nproc)

    def set(self,nproc=None):
        r"""
        Changes the status of the parallelism in SageManifolds.
        If not argument is given, the number of processors will set
        to the maximum of cores found.
        """
        self.nproc = ncpus() if nproc is None else nproc
        self.use_paral = True if self.nproc!=1 else False



# initialisation of SageManifolds Parallelism
# probably to define in another file
manifoldPara = ManifoldParallelism(1)
