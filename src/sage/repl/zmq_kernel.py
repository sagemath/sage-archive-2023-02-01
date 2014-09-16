"""
The Sage ZMQ Kernel
"""

from IPython.kernel.zmq.ipkernel import Kernel
from IPython.kernel.zmq.zmqshell import ZMQInteractiveShell
from IPython.utils.traitlets import Type

from sage.repl.interpreter import SageInteractiveShell


class SageZMQInteractiveShell(SageInteractiveShell, ZMQInteractiveShell):
    pass


class SageKernel(Kernel):    
    shell_class = Type(SageZMQInteractiveShell)

