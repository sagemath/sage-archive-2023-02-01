"""
Misc code useful for the notebook
"""

#############################################################################
#       Copyright (C) 2006, 2007 William Stein <wstein@gmail.com>
#  Distributed under the terms of the GNU General Public License (GPL)
#  The full text of the GPL is available at:
#                  http://www.gnu.org/licenses/
#############################################################################

import os

from   sage.misc.viewer     import browser

def print_open_msg(address, port, secure=False, path=""):
    s = "Open your web browser to http%s://%s:%s%s"%('s' if secure else '', address, port, path)
    t = len(s)
    if t%2:
        t += 1
        s += ' '
    n = max(t+4, 50)
    k = n - t  - 1
    j = k/2
    print '*'*n
    print '*'+ ' '*(n-2) + '*'
    print '*' + ' '*j + s + ' '*j + '*'
    print '*'+ ' '*(n-2) + '*'
    print '*'*n


import socket
def find_next_available_port(start, max_tries=100, verbose=True):
    s = socket.socket(socket.AF_INET, socket.SOCK_STREAM)
    for port in range(start, start+max_tries+1):
        try:
            s.connect(('', port))
        except socket.error, msg:
            return port
        else:
            if verbose:
                print "Port %s is already in use."%port
                print "Trying next port..."
            port += 1
    raise RuntimeError, "no available port."


def open_page(address, port, secure, path=""):
    if secure:
        rsrc = 'https'
    else:
        rsrc = 'http'

    os.system('%s %s://%s:%s%s 1>&2 > /dev/null &'%(browser(), rsrc, address, port, path))

