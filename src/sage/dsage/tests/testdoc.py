"""
These test that DSage is *really* working for normal users locally
on their system:

WARNING: Currently these non-blocking startups leave processes
hanging around!

   sage: dsage.server(blocking=False)
   sage: dsage.worker(blocking=False)
   sage: d = DSage()
   sage: sleep(0.5)
   sage: a = d('2 + 3')
   sage: a.wait()
   sage: a
   5
   sage: v = [d('%s^2'%i) for i in range(100,103)]
   sage: _ = [x.wait() for x in v]
   sage: print v
   [10000, 10201, 10404]
"""
