"""
The Sage ZMQ Kernel

Version of the IPython kernel when running Sage inside the IPython
notebook or remote IPython sessions.
"""

#*****************************************************************************
#       Copyright (C) 2015 Volker Braun <vbraun.name@gmail.com>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#  as published by the Free Software Foundation; either version 2 of
#  the License, or (at your option) any later version.
#                  http://www.gnu.org/licenses/
#*****************************************************************************


from IPython.kernel.zmq.ipkernel import Kernel
from IPython.kernel.zmq.zmqshell import ZMQInteractiveShell
from IPython.utils.traitlets import Type

from sage.repl.interpreter import SageNotebookInteractiveShell


class SageZMQInteractiveShell(SageNotebookInteractiveShell, ZMQInteractiveShell):
    pass


class SageKernel(Kernel):    
    shell_class = Type(SageZMQInteractiveShell)

