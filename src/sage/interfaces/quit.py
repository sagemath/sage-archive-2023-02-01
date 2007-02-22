"""nodoctest"""

import os, time

expect_objects = []

def expect_quitall(verbose=False):
    for P in expect_objects:
        R = P()
        if not R is None:
            try:
                R.quit(verbose=verbose)
            except RuntimeError, msg:
                if verbose:
                    print msg
    kill_spawned_jobs()

def kill_spawned_jobs():
    file = '%s/tmp/%s/spawned_processes'%(os.environ['DOT_SAGE'], os.getpid())
    if not os.path.exists(file):
        return
    for L in open(file).readlines():
        i = L.find(' ')
        pid = L[:i].strip()
        cmd = L[i+1:].strip()
        j = 0
        while j < 3:
            j += 1
            if not is_running(pid):
                break
            try:
                os.killpg(int(pid), 9)
            except OSError, msg:
                pass
            else:
                j += 1
                if j > 5:
                    os.kill(int(pid), 9)
                    break

def is_running(pid):
    """
    Return True if and only if there is a process with id pid running.
    """
    try:
        os.kill(int(pid),0)
        return True
    except (OSError, ValueError):
        return False


