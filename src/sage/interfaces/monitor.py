###############################################################################
#   SAGE: System for Algebra and Geometry Experimentation
#       Copyright (C) 2005, 2006 William Stein <wstein@gmail.com>
#  Distributed under the terms of the GNU General Public License (GPL)
#  The full text of the GPL is available at:
#                  http://www.gnu.org/licenses/
###############################################################################

import os, time

PID = os.getpid()

import sage.misc.misc as misc
F = '%s/spawn'%misc.SAGE_TMP

def monitor(pid, interval, cmd):
    cmd = 'sage-monitor %s %s %s &'%(PID, pid, interval)
    os.system(cmd)

    ##################################################################
    #
    # NOTE: On Cygwin the os.system calls is instant and the os.spawnl
    # takes *SEVERAL SECONDS*.  so we do *not* use spawn -- also spawn
    # doesn't work right under linux for this app.
    # os.spawnl(os.P_NOWAIT, 'sage-monitor', PID, pid, interval)



