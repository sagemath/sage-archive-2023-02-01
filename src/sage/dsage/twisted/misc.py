##############################################################################
#
#  DSAGE: Distributed SAGE
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
##############################################################################

import threading, sys
from twisted.internet import defer, reactor
from twisted.python.failure import Failure

# This code is from
# http://twistedmatrix.com/trac/ticket/1042
def blockingCallFromThread(func, *args, **kwargs):
    # print func
    # print args
    # print kwargs
    e = threading.Event()
    l = []
    def _got_result(result):
        print result
        l.append(result)
        e.set()
        return None
    def wrapped_func():
        d = defer.maybeDeferred(func, *args, **kwargs)
        d.addBoth(_got_result)
    reactor.callFromThread(wrapped_func)
    e.wait()
    result = l[0]
    if isinstance(result, Failure):
        # Whee!  Cross-thread exceptions!
        result.raiseException()
    else:
        return result