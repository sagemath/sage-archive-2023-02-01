"""
These test that DSage is *really* working for normal users locally
on their system:

WARNING: Currently these non-blocking startups leave processes
hanging around!
   sage: port = randint(8000, 9000)
   sage: dsage.server(blocking=False, port=port, verbose=False, ssl=False)
   sage: dsage.worker(blocking=False, port=port, verbose=False, ssl=False)
   sage: sleep(1.0)
   sage: d = DSage(port=port, ssl=False)
   sage: sleep(1.0)
   sage: a = d('2 + 3')
   sage: a.wait()
   sage: a
   5
   sage: v = [d('%s^2'%i) for i in range(100,103)]

Set timeout to 30 seconds so it will not hang the doctests indefinitely.

   sage: _ = [x.wait(timeout=30) for x in v]
   sage: print v
   [10000, 10201, 10404]
"""
