"""
These test that DSage is *really* working for normal users locally
on their system:

WARNING: Currently these non-blocking startups leave processes
hanging around!

   sage: dsage.server(blocking=False, verbose=False)
   sage: dsage.worker(blocking=False, verbose=False)
   sage: d = DSage()
   sage: sleep(0.5)
   sage: a = d('2 + 3')
   sage: a.wait()
   sage: a
   5
   sage: v = [d('%s^2'%i) for i in range(100,103)]

Set timeout to 10 seconds it it will not hang the doctests indefinitely.

   sage: _ = [x.wait(timeout=15) for x in v]
   sage: print v
   [10000, 10201, 10404]
"""
