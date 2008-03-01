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

# This code is from http://twistedmatrix.com/trac/ticket/1042

def blockingCallFromThread(f, *a, **kw):
    """
    Run a function in the reactor from a thread, and wait for the result
    synchronously, i.e. until the callback chain returned by the function get a
    result.

    @param f: the callable to run in the reactor thread
    @type f: any callable.
    @param a: the arguments to pass to C{f}.
    @param kw: the keyword arguments to pass to C{f}.

    @return: the result of the callback chain.
    @raise: any error raised during the callback chain.
    """

    import Queue
    from twisted.internet import reactor
    from twisted.internet import defer
    from twisted.python import failure

    queue = Queue.Queue()
    def _callFromThread():
        result = defer.maybeDeferred(f, *a, **kw)
        result.addBoth(queue.put)
    reactor.callFromThread(_callFromThread)
    result = queue.get()
    if isinstance(result, failure.Failure):
        result.raiseException()
    return result