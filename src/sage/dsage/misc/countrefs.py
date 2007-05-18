#!/usr/bin/env python
# Copyright 2003-2005 Andrew Bennetts.

# Simple object counting for CPython.  Abuses the gc module to find most class
# and type objects, and then reports their reference counts.  Because instances
# hold a reference to their class/type, this is a good approximation of the
# number of instances, particularly as the number of instances gets large.
# Non-heaptype types probably aren't reported, but that's not usually
# significant.

# TODO:
#  * Make this a script wrapper, like profile.py, so that you can use:
#       countrefs.py --log=refs.log -t 5 -n 30 some-script.py
#    to run some-script.py, logging the top 30 reference counts to refs.log
#    every 5 seconds.  It should still be an importable module, of course.

import gc, sys, types
import threading, time

__all__ = ['logInThread', 'mostRefs', 'printCounts']

def mostRefs(n=30):
    d = {}
    for obj in gc.get_objects():
        if type(obj) is types.InstanceType:
            cls = obj.__class__
        else:
            cls = type(obj)
        d[cls] = d.get(cls, 0) + 1
    counts = [(x[1],x[0]) for x in d.items()]
    counts.sort()
    counts = counts[-n:]
    counts.reverse()
    return counts


def printCounts(counts, file=None):
    for c, obj in counts:
        if file is None:
            print c, obj
        else:
            file.write("%s %s\n" % (c, obj))


def logInThread(n=30):
    """Start a thread to log the top n ref counts every second."""
    reflog = file('/tmp/refs.log','w')
    t = threading.Thread(target=_logRefsEverySecond, args=(reflog, n))
    t.setDaemon(True) # allow process to exit without explicitly stopping thread
    t.start()
    return t


def _logRefsEverySecond(log, n):
    while True:
        printCounts(mostRefs(n=n), file=log)
        log.write('\n')
        log.flush()
        time.sleep(1)


if __name__ == '__main__':
    counts = mostRefs()
    printCounts(counts)

