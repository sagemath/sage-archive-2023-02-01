###############################################################################
#   SAGE: System for Algebra and Geometry Experimentation
#       Copyright (C) 2005, 2006 William Stein <wstein@gmail.com>
#  Distributed under the terms of the GNU General Public License (GPL)
#  The full text of the GPL is available at:
#                  http://www.gnu.org/licenses/
###############################################################################


import os

PID = os.getpid()

def monitor(pid, interval):
    cmd = 'sage-monitor %s %s %s &'%(PID, pid, interval)
    os.system(cmd)

