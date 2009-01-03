"""
Get resource usage of process

AUTHORS:
    -- William Stein (2006-03-04): initial versoin
"""

#############################################################################
#       Copyright (C) 2006 William Stein <wstein@gmail.com>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#
#                  http://www.gnu.org/licenses/
#############################################################################

import os

def top():
    """
    Return the top output line that contains this running SAGE
    process.
    """
    U = os.uname()[0].lower()
    pid = os.getpid()

    if U == 'linux':
        cmd = 'top -b -n 1 -p %s'%pid
    elif U == 'darwin':
        cmd = 'top -l 1 |grep "^ *%s "'%pid
    elif U == 'sunos':
        cmd = 'top -b -n 65635 |grep "^ *%s "'%pid
    else:
        raise NotImplementedError, "top not implemented on platform %s"%U

    r = os.popen(cmd).read()
    r = r.strip()
    i = r.rfind('\n')
    if i == -1:
        return r
    return r[i+1:]

def get_memory_usage(t=None):
    """
    Return memory usage.

    INPUT:
        t -- None or output of previous call; (only used on Linux)

    OUTPUT:
        * Linux -- Returns float number (in megabytes)
        * OS X -- returns string (VSIZE column of top)
        * other -- not implemented for any other operating systems
    """
    U = os.uname()[0].lower()
    if U == 'linux':
        if t is None:
            return linux_memory_usage()
        else:
            return linux_memory_usage() - t
    elif U == 'darwin':
        return top().split()[-1]
    elif U == 'sunos':
        return top().split()[-5]
    else:
        raise NotImplementedError, "memory usage not implemented on platform %s"%U



########################################################################
# The following is adapted from
#   http://aspn.activestate.com/ASPN/Cookbook/Python/Recipe/286222
# Python Cookbook, by Jean Brouwers
########################################################################

_proc_status = '/proc/%d/status' % os.getpid()

def VmB(VmKey):
    """
    Function used internally by this module.
    """
    global _proc_status, _scale
     # get pseudo file  /proc/<pid>/status
    try:
        t = open(_proc_status)
        v = t.read()
        t.close()
    except:
        return 0.0  # non-Linux?
     # get VmKey line e.g. 'VmRSS:  9999  kB\n ...'
    i = v.index(VmKey)
    v = v[i:].split(None, 3)  # whitespace
    if len(v) < 3:
        return 0.0  # invalid format?
    return float(v[1])/1024.0

def linux_memory_usage():
    """
    Return memory usage in megabytes.
    """
    return VmB('VmSize:')
