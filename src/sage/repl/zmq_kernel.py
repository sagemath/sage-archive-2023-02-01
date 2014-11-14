"""
The Sage ZMQ Kernel

Version of the IPython kernel when running Sage inside the IPython
notebook or remote IPython sessions.
"""

from IPython.kernel.zmq.ipkernel import Kernel
from IPython.kernel.zmq.zmqshell import ZMQInteractiveShell
from IPython.utils.traitlets import Type

from sage.repl.interpreter import SageNotebookInteractiveShell


class SageZMQInteractiveShell(SageNotebookInteractiveShell, ZMQInteractiveShell):
    pass


class SageKernel(Kernel):    
    shell_class = Type(SageZMQInteractiveShell)

