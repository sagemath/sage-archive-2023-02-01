#!/usr/bin/env python
############################################################################
#
#   DSAGE: Distributed SAGE
#
#       Copyright (C) 2006, 2007 Yi Qiang <yqiang@gmail.com>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#
#    This code is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
#    General Public License for more details.
#
#  The full text of the GPL is available at:
#
#                  http://www.gnu.org/licenses/
############################################################################


from sage.all import *

import datetime
import sys

from twisted.spread import pb
from zope.interface import implements
from twisted.cred import portal, credentials
from twisted.spread.interfaces import IJellyable
from twisted.spread.pb import IPerspective, AsReferenceable
from twisted.python import log
from twisted.internet import defer

from sage.dsage.interface.dsage_interface import DSage
from sage.dsage.dist_functions.all import *
from sage.dsage.database.job import Job
from sage.dsage.misc.hostinfo import HostInfo
from sage.dsage.server.hostinfo_tracker import hostinfo_list
from sage.dsage.errors.exceptions import BadJobError, BadTypeError

try:
    from IPython.Shell import MTInteractiveShell
    from IPython.ipmaker import make_IPython

    import threading


    def hijack_reactor():
        """Modifies Twisted's reactor with a dummy so user code does
        not block IPython.  This function returns the original
        'twisted.internet.reactor' that has been hijacked.

        NOTE: Make sure you call this *AFTER* you've installed
        the reactor of your choice.
        """
        from twisted import internet
        orig_reactor = internet.reactor

        class DummyReactor(object):
            def __init__(self):
                self.orig_reactor = internet.reactor
            def run(self):
                pass
            def __getattr__(self, name):
                return getattr(orig_reactor, name)
            def __setattr__(self, name, value):
                return setattr(orig_reactor, name, value)

        internet.reactor = DummyReactor()
        return orig_reactor


    class IPShellTwisted(threading.Thread):
        """
        Run a Twisted reactor while in an IPython session.

        Python commands can be passed to the thread where they will be
        executed.  This is implemented by periodically checking for
        passed code using a Twisted reactor callback.
        """

        TIMEOUT = 0.1 # Millisecond interval between reactor

        def __init__(self, argv=None, user_ns=None, debug=Integer(1),
                     shell_class=MTInteractiveShell):

            from twisted.internet import reactor
            self.reactor = hijack_reactor()

            mainquit = self.reactor.stop

            # Make sure IPython keeps going after reactor stop.
            def reactorstop():
                pass
            self.reactor.stop = reactorstop
            reactorrun_orig = self.reactor.run
            self.quitting = False
            def reactorrun():
                while True and not self.quitting:
                    reactorrun_orig()
            self.reactor.run = reactorrun

            self.IP = make_IPython(argv, user_ns=user_ns, debug=debug,
                                   shell_class=shell_class,
                                   on_kill=[mainquit])

            threading.Thread.__init__(self)

        def run(self):
            self.IP.mainloop()
            self.quitting = True
            self.IP.kill()

        def mainloop(self):
            self.reactor.callLater(self.TIMEOUT, self.on_timer)
            self.start()
            self.reactor.run()
            self.join()

        def on_timer(self):
            self.IP.runcode()
            self.reactor.callLater(self.TIMEOUT, self.on_timer)


    exists = True


except ImportError:
    exists = False


if __name__ == '__main__':
    # Sample usage.

    # Create the shell object. This steals twisted.internet.reactor
    # for its own purposes, to make sure you've already installed a
    # reactor of your choice.
    shell = IPShellTwisted(
        argv=[],
        user_ns=globals())

    # Run the mainloop.  This runs the actual reactor.run() method.
    # The twisted.internet.reactor object at this point is a dummy
    # object that passes through to the actual reactor, but prevents
    # run() from being called on it again.
    print 'Starting the Distributed SAGE console...'
    shell.mainloop()

    # You must exit IPython to terminate your program.
    print 'Goodbye!'


