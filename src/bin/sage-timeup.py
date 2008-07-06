###########################################################
# Based on a script by Andrew Dalke:
#    http://projects.scipy.org/pipermail/numpy-discussion/2008-July/035415.html
#
# Included in Sage by Mike Hansen and William Stein
###########################################################


import time

seen = set()
import_order = []
elapsed_times = {}
level = 0
parent = None
children = {}

def new_import(name, globals={}, locals={}, fromlist=[]):
     global level, parent
     if name in seen:
         return old_import(name, globals, locals, fromlist)
     seen.add(name)
     import_order.append((name, level, parent))
     t1 = time.time()
     old_parent = parent
     parent = name
     level += 1
     module = old_import(name, globals, locals, fromlist)
     level -= 1
     parent = old_parent
     t2 = time.time()
     elapsed_times[name] = t2-t1
     return module

old_import = __builtins__.__import__

__builtins__.__import__ = new_import

from sage.all import *

parents = {}
for name, level, parent in import_order:
     parents[name] = parent

print "== Tree =="
for name, level,parent in import_order:
     try:
          print "%s%s: %.3f (%s)" % (" "*level, name, elapsed_times[name],
parent)
     except KeyError:
          pass

print "\n"
print "== Slowest (including children) =="
slowest = sorted((t, name) for (name, t) in elapsed_times.items())[-50:]
for elapsed_time, name in slowest[::-1]:
     print "%.3f %s (%s)" % (elapsed_time, name, parents[name])



