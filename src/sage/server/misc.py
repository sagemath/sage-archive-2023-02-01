"""
Miscellaneous Server Functions

.. WARNING::

    Parts of this file are duplicated in sagenb/misc/misc.py
"""

#*****************************************************************************
#       Copyright (C) 2006, 2007 William Stein <wstein@gmail.com>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#  as published by the Free Software Foundation; either version 2 of
#  the License, or (at your option) any later version.
#                  http://www.gnu.org/licenses/
#*****************************************************************************

import socket
from sage.misc.viewer import browser


def print_open_msg(address, port, secure=False, path=""):
    """
    Print a message on the screen suggesting that the user open their
    web browser to a certain URL.

    INPUT:

    - ``address`` -- a string; a computer address or name

    - ``port`` -- an int; a port number

    - ``secure`` -- a bool (default: False); whether to prefix the URL
      with 'http' or 'https'

    - ``path`` -- a string; the URL's path following the port.

    EXAMPLES::

        sage: from sage.server.misc import print_open_msg
        sage: print_open_msg('localhost', 8000, True)
        ****************************************************
        *                                                  *
        * Open your web browser to https://localhost:8000  *
        *                                                  *
        ****************************************************
        sage: print_open_msg('sagemath.org', 8000, False)
        ******************************************************
        *                                                    *
        * Open your web browser to http://sagemath.org:8000  *
        *                                                    *
        ******************************************************
        sage: print_open_msg('sagemath.org', 90, False)
        ****************************************************
        *                                                  *
        * Open your web browser to http://sagemath.org:90  *
        *                                                  *
        ****************************************************
        sage: print_open_msg('sagemath.org', 80, False)
        **************************************************
        *                                                *
        *  Open your web browser to http://sagemath.org  *
        *                                                *
        **************************************************
    """
    if port == 80:
        port = ''
    else:
        port = ':%s' % port
    s = "Open your web browser to http%s://%s%s%s" % ('s' if secure
                                                      else '', address,
                                                      port, path)
    t = len(s)
    if t % 2:
        t += 1
        s += ' '
    n = max(t + 4, 50)
    k = n - t - 1
    j = k/2
    print '*'*n
    print '*'+ ' '*(n - 2) + '*'
    print '*' + ' '*j + s + ' '*j + '*'
    print '*'+ ' '*(n - 2) + '*'
    print '*'*n


def find_next_available_port(host, start, max_tries=100, verbose=False):
    """
    Find the next available port on a given host, that is, a port for
    which a current connection attempt returns a 'Connection refused'
    error message.  If no port is found, raise a RuntimeError exception.

    INPUT:

    - ``host`` - address to check

    - ``start`` - an int; the starting port number for the scan

    - ``max_tries`` - an int (default: 100); how many ports to scan

    - ``verbose`` - a bool (default: True); whether to print information
      about the scan

    OUTPUT:

    - an int - the port number

    EXAMPLES::

        sage: from sage.server.misc import find_next_available_port
        sage: find_next_available_port('127.0.0.1', 9000, verbose=False)   # random output -- depends on network
        9002
    """
    from sage.misc.misc import alarm, cancel_alarm
    alarm_count = 0
    for port in range(start, start+max_tries+1):
        try:
            alarm(1)
            s = socket.socket(socket.AF_INET, socket.SOCK_STREAM)
            s.connect((host, port))
        except socket.error as msg:
            if msg[1] == 'Connection refused':
                if verbose:
                    print "Using port = %s" % port
                return port
        except KeyboardInterrupt:
            if verbose:
                print "alarm"
            alarm_count += 1
            if alarm_count >= 10:
                break
            pass
        finally:
            cancel_alarm()
        if verbose:
            print "Port %s is already in use." % port
            print "Trying next port..."
    raise RuntimeError("no available port.")


def open_page(address, port, secure, path=""):
    from os import system
    if secure:
        rsrc = 'https'
    else:
        rsrc = 'http'

    system('%s %s://%s:%s%s 1>&2 > /dev/null &' % (browser(),
                                                   rsrc, address,
                                                   port, path))
