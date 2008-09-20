r"""nodoctest
WARNING:
The following examples will not work if you have not run \code{dsage.setup()}.
    sage: from sage.dsage.misc.misc import find_open_port
    sage: port = find_open_port().next()
    sage: dsage.server(blocking=False, port=port, verbose=False, ssl=False, log_level=3)
    Going into testing mode...
    sage: dsage.worker(blocking=False, port=port, verbose=False, ssl=False, log_level=3, poll=0.1, authenticate=False)
    sage: sleep(2.0)
    sage: d = dsage.connect(port=port, ssl=False)
    sage: sleep(2.0)
    sage: a = d('2 + 3')
    sage: a.wait(timeout=30)
    sage: a
    5

Kill the job

    sage: job_id = a.kill()
    sage: v = [d('%s^2'%i) for i in range(100,103)]

Set timeout to 30 seconds so it will not hang the doctests indefinitely.

    sage: _ = [x.wait(timeout=30) for x in v]    # long time
    sage: print v                                # long time
    [10000, 10201, 10404]

    sage: _ = [x.kill() for x in v]
    sage: a = 5
    sage: b = 5
    sage: j = d('a+b', user_vars={'a': a, 'b': b})
    sage: j.wait()
    sage: j
    10
    sage: t = DistributedFunctionTest(d, 5) # long time
    sage: t.wait(timeout=60) # long time
    sage: t.done # long time
    True
    sage: t.result # long time
    15

    The following code block makes sure that things exit cleanly

    sage: dsage.kill_all()
    sage: from twisted.internet import reactor
    sage: reactor.callFromThread(reactor.stop); sleep(1)
    [DSage] Closed connection to localhost
    sage: d._dsage_thread.join()
"""
