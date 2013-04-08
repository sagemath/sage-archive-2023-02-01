"""nodoctest
"""

###############################################################################
#   Sage: System for Algebra and Geometry Experimentation
#       Copyright (C) 2007 William Stein <wstein@gmail.com>
#  Distributed under the terms of the GNU General Public License (GPL)
#  The full text of the GPL is available at:
#                  http://www.gnu.org/licenses/
###############################################################################
# Triva Note: For the name "sage-cleaner", think of the
# "The Cleaner" from Pulp Fiction:
#      http://www.frankjankowski.de/quiz/illus/keitel.jpg
###############################################################################
import os

import sage.misc.misc as misc

def cleaner(pid, cmd=''):
    if cmd != '':
        cmd = cmd.strip().split()[0]
    # This is safe, since only this process writes to this file.
    F = os.path.join(misc.SAGE_TMP, 'spawned_processes')
    if os.path.exists(F):
        o = open(F,'a')
    else:
        if not os.path.exists(str(misc.SAGE_TMP)):
            return
        o = open(F,'w')
    o.write('%s %s\n'%(pid, cmd))
    o.close()
    start_cleaner_if_not_running()

################

def start_cleaner_if_not_running():
    D = '%s/tmp_cleaner.pid'%misc.DOT_SAGE
    try:
        pid = int(open(D).read())
        os.kill(pid,0)
        return
    except (IOError, OSError, ValueError):
        os.system('sage-cleaner >/dev/null 2>/dev/null &')   # it has died



