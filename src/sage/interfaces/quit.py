"""nodoctest"""

import os, time, pexpect

expect_objects = []

def expect_quitall(verbose=False):
    for P in expect_objects:
        R = P()
        if not R is None:
            try:
                R.quit(verbose=verbose)
                pass
            except RuntimeError:
                pass
    kill_spawned_jobs()

def kill_spawned_jobs():
    file = '%s/tmp/%s/spawned_processes'%(os.environ['DOT_SAGE'], os.getpid())
    if not os.path.exists(file):
        return
    for L in open(file).readlines():
        i = L.find(' ')
        pid = L[:i].strip()
        cmd = L[i+1:].strip()
        try:
            #s = "Killing spawned job %s"%pid
            #if len(cmd) > 0:
            #    s += ' (%s)'%cmd
            #print s
            os.killpg(int(pid), 9)
        except:
            pass

def is_running(pid):
    """
    Return True if and only if there is a process with id pid running.
    """
    try:
        os.kill(int(pid),0)
        return True
    except (OSError, ValueError):
        return False


