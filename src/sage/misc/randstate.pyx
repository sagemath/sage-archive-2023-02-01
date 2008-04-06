r"""
Random Number States

AUTHORS:
    -- Carl Witty (2008-03): new file

This module manages all the available pseudo-random number generators
in Sage.  (For the rest of the documentation in this module, we will
drop the ``pseudo''.)

The goal is to allow algorithms using random numbers to be
reproducible from one run of \sage to the next, and (to the extent
possible) from one machine to the next (even across different
operating systems and architectures).

There are two parts to the API.  First we will describe the
command-line-oriented API, for setting random number generator seeds.
Then we will describe the library API, for people writing \sage
library code that uses random numbers.

We'll start with the simplest usage: setting fixed random number seeds
and showing that these lead to reproducible results.

sage: K.<x> = QQ[]
sage: G = PermutationGroup([[(1,2,3),(4,5)], [(1,2)]])
sage: rgp = Gp()
sage: def rtest():
...       current_randstate().set_seed_gp(rgp)
...       return (ZZ.random_element(1000), RR.random_element(),
...               K.random_element(), G.random_element(),
...               rgp.random(), ntl.ZZ_random(99999),
...               random())

The above test shows the results of five different random number
generators, in three different processes.  The random elements from
ZZ, RR, and K all derive from a single GMP-based random number
generator.  The random element from G comes from a GAP subprocess.
The random number from rgp is from a Pari/gp subprocess.  NTL's
ZZ_random uses a separate NTL random number generator in the main
\sage process.  And random() is from a Python \class{random.Random}
object.

Here we see that setting the random number seed really does make the
results of these random number generators reproducible.

sage: set_random_seed(0)
sage: rtest()
(303, -0.187141682737491, 1/2*x^2 - 1/95*x - 1/2, (1,2), 1698070399, 8045, 0.96619117347084138)
sage: set_random_seed(1)
sage: rtest()
(978, 0.184109262667515, -3*x^2 - 1/12, (2,3), 1046254370, 60359, 0.83350776541997362)
sage: set_random_seed(2)
sage: rtest()
(207, 0.505364206568040, 4*x^2 + 1/2, (1,2)(4,5), 2137873234, 27695, 0.19982565117278328)
sage: set_random_seed(0)
sage: rtest()
(303, -0.187141682737491, 1/2*x^2 - 1/95*x - 1/2, (1,2), 1698070399, 8045, 0.96619117347084138)
sage: set_random_seed(1)
sage: rtest()
(978, 0.184109262667515, -3*x^2 - 1/12, (2,3), 1046254370, 60359, 0.83350776541997362)
sage: set_random_seed(2)
sage: rtest()
(207, 0.505364206568040, 4*x^2 + 1/2, (1,2)(4,5), 2137873234, 27695, 0.19982565117278328)

Once we've set the random number seed, we can check what seed was used.
(This is not the current random number state; it does not change when
random numbers are generated.)

sage: set_random_seed(12345)
sage: initial_seed()
12345L
sage: rtest()
(720, 0.0216737401150802, x^2 - x, (), 1506569166, 14005, 0.92053315995181839)
sage: initial_seed()
12345L

If \code{set_random_seed()} is called with no argument, then a new
seed is automatically selected.  On operating systems that support it,
the new seed comes from \function{os.urandom}; this is intended to be
a truly random (not pseudo-random), cryptographically secure number.
(Whether it is actually cryptographically secure depends on operating
system details that are outside the control of \sage.)

If \function{os.urandom} is not supported, then the new seed comes
from the current time, which is definitely not cryptographically
secure.

sage: set_random_seed()
sage: r = rtest()
sage: r         # random
(909, -0.407373370020575, 6/7*x^2 + 1, (1,2,3)(4,5), 985329107, 21461, 0.30047071049504859)

After setting a new random number seed with \code{set_random_seed()},
we can use \function{initial_seed} to see what seed was automatically
selected, and call \function{set_random_seed} to restart the same
random number sequence.

sage: s = initial_seed()
sage: s         # random
336237747258024892084418842839280045662L
sage: set_random_seed(s)
sage: r2 = rtest()
sage: r == r2
True

Whenever \sage starts, \code{set_random_seed()} is called just before
command line interaction starts; so every \sage run starts with a
different random number seed.  This seed can be recovered with
\code{initial_seed()} (as long as the user has not set a different
seed with \function{set_random_seed}), so that the results of this run
can be reproduced in another run; or this automatically selected seed
can be overridden with, for instance, \code{set_random_seed(0)}.

We can demonstrate this startup behavior by running a new instance of
\sage as a subprocess.

sage: subsage = Sage()
sage: s = ZZ(subsage('initial_seed()'))
sage: r = ZZ(subsage('ZZ.random_element(2^200)'))
sage: s         # random
161165040149656168853863459174502758403
sage: r         # random
1273828861620427462924151488498075119241254209468761367941442
sage: set_random_seed(s)
sage: r == ZZ.random_element(2^200)
True

Note that wrappers of all the random number generation methods from
Python's \module{random} module are available at the \sage command
line, and these wrappers are properly affected by set_random_seed().

sage: set_random_seed(0)
sage: random(), getrandbits(20), uniform(5.0, 10.0), normalvariate(0, 1)
(0.111439293741037, 539332L, 8.26785106378383, 1.3893337539828183)
sage: set_random_seed(1)
sage: random(), getrandbits(20), uniform(5.0, 10.0), normalvariate(0, 1)
(0.82940228518742587, 624859L, 5.77894484361117, -0.42013668263087578)
sage: set_random_seed(0)
sage: random(), getrandbits(20), uniform(5.0, 10.0), normalvariate(0, 1)
(0.111439293741037, 539332L, 8.26785106378383, 1.3893337539828183)

That pretty much covers what you need to know for command-line use of
this module.  Now let's move to what authors of \sage library code
need to know about the module.

First, we'll cover doctesting.  Every docstring now has an implicit
\code{set_random_seed(0)} prepended.  Any uses of \code{# random} that
are based on random numbers under the control of this module should be
removed, and the reproducible answers inserted instead.

This practice has two potential drawbacks.  First, it increases the work
of maintaining doctests; for instance, in a long docstring that has
many doctests that depend on random numbers, a change near the beginning
(for instance, adding a new doctest) may invalidate all later doctests
in the docstring.  To reduce this downside, you may add calls to
\code{set_random_seed(0)} throughout the docstring (in the extreme case,
before every doctest).

Second, the \code{# random} in the doctest served as a signal to the
reader of the docstring that the result was unpredictable and that it
would not be surprising to get a different result when trying out the
examples in the doctest.  If a doctest specifically refers to
\code{ZZ.random_element()} (for instance), this is presumably enough
of a signal to render this function of \code{# random} unnecessary.
However, some doctests are not obviously (from the name) random, but
do depend on random numbers internally, such as the
\method{composition_series} method of a \class{PermutationGroup}.  In
these cases, the convention is to insert the following text at the
beginning of the EXAMPLES section.

\begin{verbatim}

        These computations use pseudo-random numbers, so we set the
        seed for reproducible testing.
            sage: set_random_seed(0)

\end{verbatim}

Note that this call to \code{set_random_seed(0)} is redundant, since
\code{set_random_seed(0)} is automatically inserted at the beginning
of every docstring; but it makes the example reproducible for somebody
who just types the lines from the doctest and doesn't know about the
automatic \code{set_random_seed(0)}.

Next, let's cover setting the random seed from library code.  The
first rule is that library code should never call
\function{set_random_seed}; this function is only for command-line
use.  Instead, if the library code wants to use a different random
seed, it should use \code{with seed(s):}; this will use the new seed
within the scope of the \code{with}, but will revert to the previous seed
once the \code{with} is completed.  (Or the library can use
\code{with seed():} to get a seed automatically selected using
\code{os.urandom()} or the current time, in the same way as described for
\code{set_random_seed()} above.)

Ideally, using \code{with seed(s):} should not affect the outer random
number sequence at all; we will call this property ``isolation.''  We
achieve isolation for most, but not all, of the random number generators
in \sage (we fail for generators, such as NTL, that do not provide an API
to retrieve the current random number state).

We'll demonstrate isolation.  First, we show the sequence of random numbers
that you get without intervening \code{with seed}.

sage: set_random_seed(0)
sage: r1 = rtest(); r1
(303, -0.187141682737491, 1/2*x^2 - 1/95*x - 1/2, (1,2), 1698070399, 8045, 0.96619117347084138)
sage: r2 = rtest(); r2
(105, -0.581229341007821, -x^2 - x - 6, (1,3,2)(4,5), 697338742, 1271, 0.001767155077382232)

We get slightly different results with an intervening \code{with seed}.

sage: set_random_seed(0)
sage: r1 == rtest()
True
sage: with seed(1): rtest()
(978, 0.184109262667515, -3*x^2 - 1/12, (2,3), 1046254370, 60359, 0.83350776541997362)
sage: r2m = rtest(); r2m
(105, -0.581229341007821, -x^2 - x - 6, (1,3,2)(4,5), 697338742, 19769, 0.001767155077382232)
sage: r2m == r2
False

We can see that \var{r2} and \var{r2m} are the same except for the
call to \function{ntl.ZZ_random}, which produces different results
with and without the \code{with seed}.

However, we do still get a partial form of isolation, even in this
case, as we see in this example:

sage: set_random_seed(0)
sage: r1 == rtest()
True
sage: with seed(1): (rtest(), rtest())
((978, 0.184109262667515, -3*x^2 - 1/12, (2,3), 1046254370, 60359, 0.83350776541997362), (138, -0.247578366457583, 2*x - 24, (), 1077419109, 10234, 0.0033332230808060803))
sage: r2m == rtest()
True

The NTL results after the \code{with seed} don't depend on how many
NTL random numbers were generated inside the \code{with seed}.

Unfortunately, Cython does not yet support the \code{with} statement.
Instead, you must use the equivalent \code{try}/\code{finally} statement:

sage: set_random_seed(0)
sage: r1 == rtest()
True
sage: ctx = seed(1)
sage: try:
...       ctx.__enter__()
...       rtest()
... finally:
...       ctx.__exit__(None, None, None)
<sage.misc.randstate.randstate object at 0x...>
(978, 0.184109262667515, -3*x^2 - 1/12, (2,3), 1046254370, 60359, 0.83350776541997362)
False
sage: r2m == rtest()
True

(In general, the above code is not exactly equivalent to the \code{with}
statement, because if an exception happens in the body, the real
\code{with} statement will pass the exception information as parameters to
the \method{__exit__} method.  However, our \method{__exit__} method
ignores the exception information anyway, so the above is equivalent in
our case.)

Now we come to the last part of the documentation: actually generating
random numbers in library code.  First, the easy case: if you generate
random numbers only by calling other \sage library code (such as
\method{random_element} methods on parents), you don't need to do
anything special; the other code presumably already interacts with
this module correctly.

Otherwise, it depends on what random number generator you want to use.

\begin{description}
\item[gmp_randstate_t] If you want to use some random number generator
that takes a \code{gmp_randstate_t} (like \code{mpz_urandomm} or
\code{mpfr_urandomb}), then use code like the following:
\begin{verbatim}
from sage.misc.randstate cimport randstate, current_randstate
...

    cdef randstate rstate = current_randstate()
\end{verbatim}
then a \code{gmp_randstate_t} is available as \code{rstate.gmp_state}.

Fetch the current \class{randstate} with \code{current_randstate()} in
every function that wants to use it; don't cache it globally or
in a class.  (Such caching would break \method{set_random_seed}).

\item[Python] If you want to use the random number generators from the
\module{random} module, you have two choices.  The slightly easier choice
is to import functions from \module{sage.misc.prandom}; for instance,
you can simply replace \code{from random import randrange} with
\code{from sage.misc.prandom import randrange}.  However, this is
slighly less efficient, because the wrappers in \module{sage.misc.prandom}
look up the current \class{randstate} on each call.  If you're generating
many random numbers in a row, it's faster to instead do
\begin{verbatim}
from sage.misc.randstate import current_randstate
...

    randrange = current_randstate().python_random().randrange
\end{verbatim}

Fetch the current \class{randstate} with \code{current_randstate()} in
every function that wants to use it; don't cache the
\class{randstate}, the \class{Random} object returned by
\method{python_random}, or the bound methods on that \class{Random}
object globally or in a class.  (Such caching would break
\method{set_random_seed}).

\item[GAP] If you are calling code in GAP that uses random numbers,
call \method{set_seed_gap} at the beginning of your function, like this:
\begin{verbatim}
from sage.misc.randstate import current_randstate
...

    current_randstate().set_seed_gap()
\end{verbatim}

Fetch the current \class{randstate} with \code{current_randstate()} in
every function that wants to use it; don't cache it globally or
in a class.  (Such caching would break \method{set_random_seed}).

\item[Pari] If you are calling code in the Pari library that uses
random numbers, call \method{set_seed_pari} at the beginning of your
function, like this:
\begin{verbatim}
from sage.misc.randstate import current_randstate
...

    current_randstate().set_seed_pari()
\end{verbatim}

Fetch the current \class{randstate} with \code{current_randstate()} in
every function that wants to use it; don't cache it globally or
in a class.  (Such caching would break \method{set_random_seed}).

\item[Pari/gp] If you are calling code in a Pari/gp subprocess that
uses random numbers, call \method{set_seed_gp} at the beginning of
your function, like this:
\begin{verbatim}
from sage.misc.randstate import current_randstate
...

    current_randstate().set_seed_gp()
\end{verbatim}
This will set the seed in the gp process in sage.interfaces.gp.gp.
If you have a different gp process, say in the variable \var{my_gp}, then
call \code{set_seed_gp(my_gp)} instead.

Fetch the current \class{randstate} with \code{current_randstate()} in
every function that wants to use it; don't cache it globally or
in a class.  (Such caching would break \method{set_random_seed}).

\item[libc] If you are writing code that calls the libc function
\code{random()}: don't!  The \code{random()} function does not give
reproducible results across different operating systems, so we can't
make portable doctests for the results.  Instead, do:
\begin{verbatim}
from sage.misc.randstate cimport random
\end{verbatim}
The \function{random} function in \module{sage.misc.randstate} gives a
31-bit random number (the same range as the \function{random} in libc,
on all systems I'm familiar with), but it uses the \code{gmp_randstate_t}
in the current \class{randstate}, so it is portable.

However, you may still need to set the libc random number state; for
instance, if you are wrapping a library that uses \code{random()}
internally and you don't want to change the library.  In that case,
call \method{set_seed_libc} at the beginning of your function, like this:
\begin{verbatim}
from sage.misc.randstate import current_randstate
...

    current_randstate().set_seed_libc()
\end{verbatim}

Fetch the current \class{randstate} with \code{current_randstate()} in
every function that wants to use it; don't cache it globally or
in a class.  (Such caching would break \method{set_random_seed}).

\end{description}
"""

cdef extern from "mpz_pylong.h":
    cdef int mpz_set_pylong(mpz_t dst, src) except -1

cdef extern from "stdlib.h":
    long RAND_MAX
    long c_libc_random "random"()
    void c_libc_srandom "srandom"(unsigned int seed)

import binascii
import os
import time
import weakref

use_urandom = False
# Check whether os.urandom() works.
try:
    os.urandom(1)
    use_urandom = True
except NotImplementedError:
    pass

# Holds the current randstate object.
cdef randstate _current_randstate

# For each kind of random number generator, keep track of which
# randstate object was the most recent one to seed it.  (Before
# the generator has been seeded, these will be None.)
cdef randstate _libc_seed_randstate
cdef randstate _ntl_seed_randstate
cdef randstate _gap_seed_randstate
cdef randstate _pari_seed_randstate
# For each gp subprocess that has been seeded, keep track of which
# randstate object was the most recent one to seed it.
_gp_seed_randstates = weakref.WeakKeyDictionary()

cpdef randstate current_randstate():
    r"""
    Return the current random number state.

    EXAMPLES:
        sage: current_randstate()
        <sage.misc.randstate.randstate object at 0x...>
        sage: current_randstate().python_random().random()
        0.111439293741037
    """
    return _current_randstate

# Keep track of the stack of randstates involved in "with seed(s):".
randstate_stack = []

# The following components of Sage use random numbers that I have not
# figured out how to seed: fpLLL, mwrank, mpfi
# (there are probably others; these are the ones I noticed while trying
# to remove "# random" from doctests)

cdef class randstate:
     r"""
     The randstate class.  This class keeps track of random number
     states and seeds.  Type \code{sage.misc.randstate?} for much more
     information on random numbers in \sage.
     """

     def __init__(self, seed=None):
         r"""
         Initialize a new \class{randstate} object with the given seed
         (which must be coercible to a Python long).

         If no seed is given, then a seed is automatically selected
         using \function{os.urandom} if it is available, or the current
         time otherwise.

         EXAMPLES:
             sage: from sage.misc.randstate import randstate
             sage: r = randstate(54321); r
             <sage.misc.randstate.randstate object at 0x...>
             sage: r.seed()
             54321L
             sage: r = randstate(); r
             <sage.misc.randstate.randstate object at 0x...>
             sage: r.seed()     # random
             305866218880103397618377824640007711767L

         Note that creating a \class{randstate} with a seed of 0
         is vastly faster than any other seed (over a thousand times
         faster in my test).

             sage: timeit('randstate(0)') # random
             625 loops, best of 3: 1.38 us per loop
             sage: timeit('randstate(1)') # random
             125 loops, best of 3: 3.59 ms per loop
         """
         gmp_randinit_default(self.gmp_state)

         cdef mpz_t mpz_seed

         if seed is None:
             if use_urandom:
                 seed = long(binascii.hexlify(os.urandom(16)), 16)
             else:
                 seed = long(time.time() * 256)
         else:
             seed = long(seed)

         # If seed==0, leave it at the default seed used by
         # gmp_randinit_default()
         if seed:
             mpz_init(mpz_seed)
             mpz_set_pylong(mpz_seed, seed)
             gmp_randseed(self.gmp_state, mpz_seed)
             mpz_clear(mpz_seed)

         self._seed = seed

     def seed(self):
         r"""
         Return the initial seed of a randstate object.  (This is not
         the current state; it does not change when you get random
         numbers.)

         EXAMPLES:
             sage: from sage.misc.randstate import randstate
             sage: r = randstate(314159)
             sage: r.seed()
             314159L
             sage: r.python_random().random()
             0.111439293741037
             sage: r.seed()
             314159L
         """
         return self._seed

     def python_random(self):
         r"""
         Return a \class{random.Random} object.  The first time it is
         called on a given \class{randstate}, a new \class{random.Random}
         is created (seeded from the \emph{current} \class{randstate});
         the same object is returned on subsequent calls.

         It is expected that \method{python_random} will only be
         called on the current \class{randstate}.

         EXAMPLES:
             sage: set_random_seed(5)
             sage: rnd = current_randstate().python_random()
             sage: rnd.random()
             0.013558022446944151
             sage: rnd.randrange(1000)
             544
         """
         if self._python_random is not None:
             return self._python_random

         import random
         from sage.rings.integer_ring import ZZ
         rand = random.Random()
         rand.seed(long(ZZ.random_element(long(1)<<128)))
         self._python_random = rand
         return rand

     cpdef ZZ_seed(self):
         r"""
         When called on the current randstate, returns a 128-bit Integer
         suitable for seeding another random number generator.

         EXAMPLES:
             sage: set_random_seed(1414)
             sage: current_randstate().ZZ_seed()
             48314508034782595865062786044921182484
         """
         from sage.rings.integer_ring import ZZ
         return ZZ.random_element(long(1)<<128)

     cpdef long_seed(self):
         r"""
         When called on the current randstate, returns a 128-bit
         Python long suitable for seeding another random number generator.

         EXAMPLES:
             sage: set_random_seed(1618)
             sage: current_randstate().long_seed()
             256056279774514099508607350947089272595L
         """
         from sage.rings.integer_ring import ZZ
         return long(ZZ.random_element(long(1)<<128))

     cpdef set_seed_libc(self, bint force):
         r"""
         Checks to see if \code{self} was the most recent randstate
         to seed the libc random number generator.  If not, seeds the
         libc random number generator.  (Do not use the libc random
         number generator if you have a choice; its randomness is poor,
         and the random number sequences it produces are not portable
         across operating systems.)

         If the argument \var{force} is \code{True}, seeds the generator
         unconditionally.

         EXAMPLES:
             sage: from sage.misc.randstate import _doctest_libc_random
             sage: set_random_seed(0xBAD)
             sage: current_randstate().set_seed_libc(False)
             sage: _doctest_libc_random()   # random
             1070075918
         """
         global _libc_seed_randstate
         if force or _libc_seed_randstate is not self:
             c_libc_srandom(gmp_urandomb_ui(self.gmp_state, sizeof(int)*8))
             _libc_seed_randstate = self

     cpdef set_seed_ntl(self, bint force):
         r"""
         Checks to see if \code{self} was the most recent randstate
         to seed the NTL random number generator.  If not, seeds
         the generator.  If the argument \var{force} is \code{True},
         seeds the generator unconditionally.

         EXAMPLES:
             sage: set_random_seed(2008)

         This call is actually redundant; \function{ntl.ZZ_random} will
         seed the generator itself.  However, we put the call in
         to make the coverage tester happy.
             sage: current_randstate().set_seed_ntl(False)
             sage: ntl.ZZ_random(10^40)
             9121810031060285085432621721962973171231
         """
         global _ntl_seed_randstate
         if force or _ntl_seed_randstate is not self:
             import sage.libs.ntl.ntl_ZZ as ntl_ZZ
             from sage.rings.integer_ring import ZZ
             ntl_ZZ.ntl_setSeed(ZZ.random_element(long(1)<<128))
             _ntl_seed_randstate = self

     def set_seed_gap(self):
         r"""
         Checks to see if \code{self} was the most recent randstate
         to seed the GAP random number generator.  If not, seeds
         the generator.

         EXAMPLES:
             sage: set_random_seed(99900000999)
             sage: current_randstate().set_seed_gap()
             sage: gap.Random(1, 10^50)
             24922966241004817547714978350347109473543873184564
         """
         global _gap_seed_randstate
         if _gap_seed_randstate is not self:
             from sage.interfaces.gap import gap

             if self._gap_saved_seed is not None:
                 seed = self._gap_saved_seed
             else:
                 import sage.rings.integer_ring as integer_ring
                 from sage.rings.integer_ring import ZZ
                 seed = ZZ.random_element(long(1)<<128)

             prev_seed = gap.Reset(gap.GlobalMersenneTwister, seed)

             if _gap_seed_randstate is not None:
                 _gap_seed_randstate._gap_saved_seed = prev_seed

             _gap_seed_randstate = self

     def set_seed_gp(self, gp=None):
         r"""
         Checks to see if \code{self} was the most recent randstate
         to seed the random number generator in the given instance
         of gp.  (If no instance is given, uses the one in
         sage.interfaces.gp.gp.)  If not, seeds the generator.

         EXAMPLES:
             sage: set_random_seed(987654321)
             sage: current_randstate().set_seed_gp()
             sage: gp.random()
             331225388
         """
         if gp is None:
             import sage.interfaces.gp
             gp = sage.interfaces.gp.gp

         cdef randstate prev

         try:
             prev = _gp_seed_randstates[gp]
         except KeyError:
             prev = None


         if prev is not self:
             if self._gp_saved_seeds is not None and self._gp_saved_seeds.has_key(gp):

                 seed = self._gp_saved_seeds[gp]
             else:
                 seed = self.c_random()

             prev_seed = gp.getrand()
             gp.setrand(seed)

             if prev is not None:
                 if prev._gp_saved_seeds == None:
                     prev._gp_saved_seeds = weakref.WeakKeyDictionary()
                 prev._gp_saved_seeds[gp] = prev_seed

             _gp_seed_randstates[gp] = self

     def set_seed_pari(self):
         r"""
         Checks to see if \code{self} was the most recent randstate to
         seed the Pari random number generator.  If not, seeds the
         generator.

         EXAMPLES:
             sage: set_random_seed(5551212)
             sage: current_randstate().set_seed_pari()
             sage: pari.pari_rand31()
             1530973455
         """
         global _pari_seed_randstate
         if _pari_seed_randstate is not self:
             from sage.libs.pari.gen import pari

             if self._pari_saved_seed is not None:
                 seed = self._pari_saved_seed
             else:
                 seed = self.c_random()

             prev_seed = pari.getrand()
             pari.setrand(seed)

             if _pari_seed_randstate is not None:
                 _pari_seed_randstate._pari_saved_seed = prev_seed

             _pari_seed_randstate = self

     cpdef int c_random(self):
         r"""
         Returns a 31-bit random number.  Intended for internal
         use only; instead of calling current_randstate().c_random(),
         it is equivalent (but probably faster) to call
         sage.misc.randstate.random().

         EXAMPLES:
             sage: set_random_seed(1207)
             sage: current_randstate().c_random()
             2008037228

         We verify the equivalence mentioned above.
             sage: from sage.misc.randstate import random
             sage: set_random_seed(1207)
             sage: random()
             2008037228
         """
         return gmp_urandomb_ui(self.gmp_state, 31)

     cpdef double c_rand_double(self):
         r"""
         Returns a random floating-point number between 0 and 1.

         EXAMPLES:
             sage: set_random_seed(2718281828)
             sage: current_randstate().c_rand_double()
             0.22437207488974298
         """
         cdef double a = gmp_urandomb_ui(self.gmp_state, 25) * (1.0 / 33554432.0) # divide by 2^25
         cdef double b = gmp_urandomb_ui(self.gmp_state, 28) * (1.0 / 9007199254740992.0) # divide by 2^53
         return a+b

     def __dealloc__(self):
         r"""
         Free up the memory from the gmp_randstate_t in a randstate.

         EXAMPLES:
             sage: from sage.misc.randstate import randstate
             sage: foo = randstate()
             sage: foo = None
         """
         gmp_randclear(self.gmp_state)

     def __enter__(self):
         r"""
         Use a randstate object as a with statement context manager;
         switches this randstate to be the current randstate,
         to be switched back on exit from the with statement.

         For this purpose, we usually use the seed alias for randstate.

         EXAMPLES:
             sage: from sage.misc.randstate import randstate
             sage: seed is randstate
             True
             sage: set_random_seed(-12345)
             sage: ZZ.random_element(10^30)
             197130468050826967386035500824
             sage: ZZ.random_element(10^30)
             601704412330400807050962541983
             sage: set_random_seed(-12345)
             sage: ZZ.random_element(10^30)
             197130468050826967386035500824
             sage: with seed(12345):
             ...       ZZ.random_element(10^30)
             197130468050826967386035500824
             sage: ZZ.random_element(10^30)
             601704412330400807050962541983
         """
         global _current_randstate
         randstate_stack.append(_current_randstate)
         _current_randstate = self
         return self

     def __exit__(self, ty, value, traceback):
         r"""
         Use a randstate object as a with statement context manager;
         restores the previous randstate as the current randstate.

         For this purpose, we usually use the seed alias for randstate.

         EXAMPLES:
             sage: from sage.misc.randstate import randstate
             sage: seed is randstate
             True
             sage: set_random_seed(-12345)
             sage: ZZ.random_element(10^30)
             197130468050826967386035500824
             sage: ZZ.random_element(10^30)
             601704412330400807050962541983
             sage: set_random_seed(-12345)
             sage: ZZ.random_element(10^30)
             197130468050826967386035500824
             sage: with seed(12345):
             ...       ZZ.random_element(10^30)
             197130468050826967386035500824
             sage: ZZ.random_element(10^30)
             601704412330400807050962541983
         """
         global _current_randstate
         _current_randstate = randstate_stack.pop()
         return False

cpdef set_random_seed(seed=None):
    r"""
    Set the current random number seed from the given \var{seed}
    (which must be coercible to a Python long).

    If no seed is given, then a seed is automatically selected
    using \function{os.urandom} if it is available, or the current
    time otherwise.

    Type \code{sage.misc.randstate?} for much more
    information on random numbers in \sage.

    This function is only intended for command line use.  Never call
    this from library code; instead, use "with seed(s):".

    Note that setting the random number seed to 0 is much faster than
    using any other number.

    EXAMPLES:
        sage: set_random_seed(5)
        sage: initial_seed()
        5L
    """
    global _current_randstate
    _current_randstate = randstate(seed)

set_random_seed()

# Create an alias for randstate to be used in context managers
seed = randstate

cpdef int random():
    r"""
    Returns a 31-bit random number.  Intended as a drop-in replacement for
    the libc \code{random()} function.

    EXAMPLES:
        sage: set_random_seed(31)
        sage: from sage.misc.randstate import random
        sage: random()
        32990711
    """
    return gmp_urandomb_ui(_current_randstate.gmp_state, 31)

def initial_seed():
    r"""
    Returns the initial seed used to create the current randstate.

    EXAMPLES:
        sage: set_random_seed(42)
        sage: initial_seed()
        42L

    If you set a random seed (by failing to specify the seed), this is how
    you retrieve the seed actually chosen by \sage.  This can also be
    used to retrieve the seed chosen for a new \sage run (if the user
    has not used set_random_seed()).

        sage: set_random_seed()
        sage: initial_seed()          # random
        121030915255244661507561642968348336774L
    """
    return _current_randstate._seed

def benchmark_libc():
    r"""
    This function was used to test whether moving from libc to GMP's
    Mersenne Twister for random numbers would be a significant slowdown.

    EXAMPLES:
        sage: from sage.misc.randstate import benchmark_libc, benchmark_mt
        sage: timeit('benchmark_libc()')  # random
        125 loops, best of 3: 1.95 ms per loop
        sage: timeit('benchmark_mt()')    # random
        125 loops, best of 3: 2.12 ms per loop
    """
    cdef int i
    cdef randstate rstate = _current_randstate
    for i from 0 <= i < 100000:
        c_libc_random()

def benchmark_mt():
    r"""
    This function was used to test whether moving from libc to GMP's
    Mersenne Twister for random numbers would be a significant slowdown.

    EXAMPLES:
        sage: from sage.misc.randstate import benchmark_libc, benchmark_mt
        sage: timeit('benchmark_libc()')  # random
        125 loops, best of 3: 1.95 ms per loop
        sage: timeit('benchmark_mt()')    # random
        125 loops, best of 3: 2.11 ms per loop
    """
    cdef int i
    cdef randstate rstate = _current_randstate
    for i from 0 <= i < 100000:
        gmp_urandomb_ui(rstate.gmp_state, 32)

cpdef int _doctest_libc_random():
    r"""
    Returns the result of \code{random()} from libc.

    Only for use in doctests; this should not actually be used in \sage,
    since the resulting random number stream is not portable across
    operating systems.

    EXAMPLES:
        sage: from sage.misc.randstate import _doctest_libc_random
        sage: _doctest_libc_random()     # random
        910236436
    """
    return c_libc_random()
