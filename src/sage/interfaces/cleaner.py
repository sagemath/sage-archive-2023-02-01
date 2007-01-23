###############################################################################
#   SAGE: System for Algebra and Geometry Experimentation
#       Copyright (C) 2007 William Stein <wstein@gmail.com>
#  Distributed under the terms of the GNU General Public License (GPL)
#  The full text of the GPL is available at:
#                  http://www.gnu.org/licenses/
###############################################################################
# Triva Note: For the name "sage-cleaner", think of the
# "The Cleaner" from Pulp Fiction:
#      http://www.frankjankowski.de/quiz/illus/keitel.jpg
###############################################################################


import sage.misc.misc as misc
F = '%s/spawned_processes'%misc.SAGE_TMP

def cleaner(pid, cmd):
    cmd = cmd.strip().split()[0]
    # This is safe, since only this process writes to this file.
    o = open(F,'a')
    o.write('%s %s\n'%(pid, cmd))
    o.close()



