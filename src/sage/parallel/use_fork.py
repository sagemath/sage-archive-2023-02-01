class p_iter_fork:
    def __init__(self, ncpus, timeout=0):
        self.ncpus = ncpus
        self.timeout = timeout

    def __call__(self, f, inputs):
        """
        Parallel iterator using fork.

        INPUT:

            - `f` -- a Python function that need not be picklable or anything else!

            - ``inputs`` -- a list of pickleable pairs (args, kwds), where
                 args is a tuple and kwds is a dictionary.

        OUTPUT:

        EXAMPLES::

        """
        n = self.ncpus
        v = list(inputs)
        import os, sys, signal
        from sage.structure.sage_object import save, load
        from sage.misc.misc import tmp_dir, walltime
        dir = tmp_dir()
        timeout = self.timeout

        workers = {}
        try:
            while len(v) > 0 or len(workers) > 0:
                # Spawn up to n subprocesses
                while len(v)>0 and len(workers) < n:
                    pid = os.fork()
                    if pid:
                        # The parent master process
                        workers[pid] = [v[0], walltime(), '']
                        del v[0]
                    else:
                        # The subprocess
                        x = v[0]
                        try:
                            # Make it so all stdout is sent to a file so it can
                            # be displayed.
                            sys.stdout = open('%s/%s.out'%(dir, os.getpid()), 'w')

                            # Run some commands to tell sage that its pid has changed.
                            import sage.misc.misc
                            reload(sage.misc.misc)
                            from sage.interfaces.quit import expect_objects
                            for I in expect_objects:
                                I1 = I()
                                if I1:
                                    I1._expect = None  #invalidate this interface
                                    I1._Expect__initialized = False
                                    I1._session_number += 1

                            # Now evaluate the function
                            value = f(*x[0], **x[1])

                            # And save the result to disk.
                            save(value, '%s/%s.sobj'%(dir, os.getpid()))
                        except Exception, msg:
                            print msg
                        finally:
                            sys.stdout.flush()
                            os._exit(0)

                if len(workers) > 0:
                    # Now wait for one subprocess to finish and report the result.
                    # However, wait at most the time since the oldest process started
                    if timeout:
                        def mysig(a,b):
                            raise RuntimeError
                        oldest = min([X[1] for X in workers.values()])
                        signal.signal(signal.SIGALRM, mysig)
                        signal.alarm(int(walltime() - oldest)+1)
                    try:
                        pid = os.wait()[0]
                        signal.signal(signal.SIGALRM, signal.SIG_IGN)
                    except RuntimeError:
                        signal.signal(signal.SIGALRM, signal.SIG_IGN)
                        # Kill workers that are too old
                        for pid, X in workers.iteritems():
                            if walltime() - X[1] > timeout:
                                print "Killing subprocess %s with input %s which took too long"%(pid, X[0])
                                os.kill(pid,9)
                                X[-1] = ' (timed out)'
                    else:
                        # collect data from process that successfuly terminated
                        sobj = '%s/%s.sobj'%(dir, pid)
                        if not os.path.exists(sobj):
                            X = "NO DATA" + workers[pid][-1]  # the message field
                        else:
                            X = load(sobj)
                            os.unlink(sobj)
                        out = '%s/%s.out'%(dir, pid)
                        if not os.path.exists(out):
                            output = "NO OUTPUT"
                        else:
                            output = open(out).read()
                            os.unlink(out)

                        if output.strip():
                            print output,

                        yield (workers[pid][0], X)
                        del workers[pid]

        except Exception, msg:
            print msg

        finally:
            try:
                for X in os.listdir(dir):
                    os.unlink(dir + '/' + X)
                os.rmdir(dir)
            except OSError, msg:
                print msg

            # kill -9 everything in workers that is left.
            if len(workers) > 0:
                print "Killing any remaining workers..."
                for pid in workers.keys():
                    try:
                        os.kill(pid, 9)
                    except OSError, msg:
                        print msg
                os.wait()

