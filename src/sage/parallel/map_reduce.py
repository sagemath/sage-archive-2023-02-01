r"""
Parallel computations using RecursivelyEnumeratedSet and Map-Reduce

There exists an efficient way to distribute computations when you have a set
`S` of objects defined by :func:`RecursivelyEnumeratedSet` (see
:mod:`sage.sets.recursively_enumerated_set` for more details) over which you
would like to perform the following kind of operations :

* Compute the cardinality of a (very large) set defined recursively (through a
  call to
  :class:`RecursivelyEnumeratedSet of forest type<sage.combinat.backtrack.SearchForest>`)

* More generally, compute any kind of generating series over this set

* Test a conjecture : i.e. find an element of `S` satisfying a specific
  property; conversely, check that all of them do

* Count/list the elements of `S` having a specific property

* Apply any map/reduce kind of operation over the elements of `S`

AUTHORS :

- Florent Hivert -- code, documentation (2012-2016)

- Jean Baptiste Priez -- prototype, debugging help on MacOSX (2011-June, 2016)
- Nathann Cohen -- Some doc (2012)

Contents
--------

- :ref:`basic-usage`
- :ref:`advanced-use`
- :ref:`profiling`
- :ref:`logging`
- :ref:`protocol-description`
- :ref:`examples`

How is this different from usual MapReduce ?
--------------------------------------------

This implementation is specific to
:class:`RecursivelyEnumeratedSet of forest type<sage.combinat.backtrack.SearchForest>`,
and uses its properties to do its job. Not only mapping
and reducing is done on different processors but also **generating the elements
of** `S`.

.. _basic-usage:

How can I use all that stuff?
-----------------------------

First, you need the information necessary to describe a
:class:`RecursivelyEnumeratedSet of forest
type<sage.combinat.backtrack.SearchForest>` representing your set `S` (see
:mod:`sage.sets.recursively_enumerated_set`).  Then, you need to provide a Map
function as well as a Reduce function. Here are some examples :

* **Counting the number of elements**: In this situation, the map function
  can be set to ``lambda x : 1``, and the reduce function just adds the
  values together, i.e. ``lambda x,y : x+y``.

  Here's the Sage code for binary words of length `\leq 16` ::

      sage: seeds = [[]]
      sage: succ = lambda l: [l+[0], l+[1]] if len(l) <= 15 else []
      sage: S = RecursivelyEnumeratedSet(seeds, succ,
      ....:   structure='forest', enumeration='depth')
      sage: map_function = lambda x: 1
      sage: reduce_function = lambda x,y: x+y
      sage: reduce_init = 0
      sage: S.map_reduce(map_function, reduce_function, reduce_init)
      131071

  One can check that this is indeed the number of binary words of
  length `\leq 16` ::

      sage: factor(131071 + 1)
      2^17


  Note that the function mapped and reduced here are equivalent to the default
  values of the :meth:`sage.combinat.backtrack.SearchForest.map_reduce` method
  so that to compute the number of element you only need to call::

      sage: S.map_reduce()
      131071

  You don't need to use :func:`RecursivelyEnumeratedSet`, you can use directly
  :class:`RESetMapReduce`. This is needed if you want to have fine
  control over the parallel execution (see :ref:`advanced-use` below)::

      sage: from sage.parallel.map_reduce import RESetMapReduce
      sage: S = RESetMapReduce(
      ....:   roots = [[]],
      ....:   children = lambda l: [l+[0], l+[1]] if len(l) <= 15 else [],
      ....:   map_function = lambda x : 1,
      ....:   reduce_function = lambda x,y: x+y,
      ....:   reduce_init = 0 )
      sage: S.run()
      131071

* **Generating series**: In this situation, the map function associates a
  monomial to each element of `S`, while the Reduce function is still equal to
  ``lambda x,y : x+y``.

  Here's the Sage code for binary words of length `\leq 16` ::

      sage: S = RecursivelyEnumeratedSet(
      ....:   [[]], lambda l: [l+[0], l+[1]] if len(l) < 16 else [],
      ....:   structure='forest', enumeration='depth')
      sage: sp = S.map_reduce(
      ....:   map_function = lambda z: x**len(z),
      ....:   reduce_function = lambda x,y: x+y,
      ....:   reduce_init = 0 )
      sage: sp
      65536*x^16 + 32768*x^15 + 16384*x^14 + 8192*x^13 + 4096*x^12 + 2048*x^11 + 1024*x^10 + 512*x^9 + 256*x^8 + 128*x^7 + 64*x^6 + 32*x^5 + 16*x^4 + 8*x^3 + 4*x^2 + 2*x + 1

  This is of course `\sum_{i=0}^{i=16} (2x)^i`::

      sage: bool(sp == sum((2*x)^i for i in range(17)))
      True

  Here is another example where we count permutations of size `\leq 8` (here
  we use the default values)::

      sage: S = RecursivelyEnumeratedSet( [[]],
      ....:   lambda l: ([l[:i] + [len(l)] + l[i:] for i in range(len(l)+1)]
      ....:               if len(l) < 8 else []),
      ....:   structure='forest', enumeration='depth')
      sage: sp = S.map_reduce(lambda z: x**len(z)); sp
      40320*x^8 + 5040*x^7 + 720*x^6 + 120*x^5 + 24*x^4 + 6*x^3 + 2*x^2 + x + 1

  This is of course `\sum_{i=0}^{i=8} i! x^i`::

      sage: bool(sp == sum(factorial(i)*x^i for i in range(9)))
      True

* **Post Processing**: We now demonstrate the use of ``post_process``. We
  generate the permutation as previously, but we only perform the map/reduce
  computation on those of even ``len``. Of course we get the even part of the
  previous generating series::

      sage: S = RecursivelyEnumeratedSet( [[]],
      ....:   lambda l: ([l[:i] + [len(l)+1] + l[i:] for i in range(len(l)+1)]
      ....:               if len(l) < 8 else []),
      ....:   post_process = lambda l : l if len(l) % 2 == 0 else None,
      ....:   structure='forest', enumeration='depth')
      sage: sp = S.map_reduce(lambda z: x**len(z)); sp
      40320*x^8 + 720*x^6 + 24*x^4 + 2*x^2 + 1

  This is also useful for example to call a constructor on the generated
  elements::

      sage: S = RecursivelyEnumeratedSet( [[]],
      ....:   lambda l: ([l[:i] + [len(l)+1] + l[i:] for i in range(len(l)+1)]
      ....:               if len(l) < 5 else []),
      ....:   post_process = lambda l : Permutation(l) if len(l) == 5 else None,
      ....:   structure='forest', enumeration='depth')
      sage: sp = S.map_reduce(lambda z: x**(len(z.inversions()))); sp
      x^10 + 4*x^9 + 9*x^8 + 15*x^7 + 20*x^6 + 22*x^5 + 20*x^4 + 15*x^3 + 9*x^2 + 4*x + 1

  We get here a polynomial called the `x`-factorial of `5` that is
  `\prod_{i=1}^{i=5} \frac{1-x^i}{1-x}`::

      sage: (prod((1-x^i)/(1-x) for i in range(1,6))).simplify_rational()
      x^10 + 4*x^9 + 9*x^8 + 15*x^7 + 20*x^6 + 22*x^5 + 20*x^4 + 15*x^3 + 9*x^2 + 4*x + 1


* **Listing the objects**: One can also compute the list of objects in a
  :class:`RecursivelyEnumeratedSet of forest type<sage.combinat.backtrack.SearchForest>`
  using :class:`RESetMapReduce`. As an example, we compute the set of number
  beetween 1 and 63, generated by their binary expansion::

      sage: S = RecursivelyEnumeratedSet( [1],
      ....:    lambda l: [(l<<1)|0, (l<<1)|1] if l < 1<<5 else [],
      ....:    structure='forest', enumeration='depth')

  Here is the list computed without :class:`RESetMapReduce`::

      sage: serial = list(S)
      sage: serial
      [1, 2, 4, 8, 16, 32, 33, 17, 34, 35, 9, 18, 36, 37, 19, 38, 39, 5, 10, 20, 40, 41, 21, 42, 43, 11, 22, 44, 45, 23, 46, 47, 3, 6, 12, 24, 48, 49, 25, 50, 51, 13, 26, 52, 53, 27, 54, 55, 7, 14, 28, 56, 57, 29, 58, 59, 15, 30, 60, 61, 31, 62, 63]

  Here is how to perform the parallel computation. The order of the lists
  depends on the synchronisation of the various computation processes and
  therefore should be considered as random::

      sage: parall = S.map_reduce( lambda x: [x], lambda x,y: x+y, [] )
      sage: parall   # random
      [1, 3, 7, 15, 31, 63, 62, 30, 61, 60, 14, 29, 59, 58, 28, 57, 56, 6, 13, 27, 55, 54, 26, 53, 52, 12, 25, 51, 50, 24, 49, 48, 2, 5, 11, 23, 47, 46, 22, 45, 44, 10, 21, 43, 42, 20, 41, 40, 4, 9, 19, 39, 38, 18, 37, 36, 8, 17, 35, 34, 16, 33, 32]
      sage: sorted(serial) == sorted(parall)
      True


.. _advanced-use:

Advanced use
------------

Fine control of the execution of a map/reduce computations is obtained by
passing parameters to the :meth:`RESetMapReduce.run` method. One can use the
three following parameters:

- ``max_proc`` -- maximum number of process used.
  default: number of processor on the machine
- ``timeout`` -- a timeout on the computation (default: ``None``)
- ``reduce_locally`` -- whether the workers should reduce locally
  their work or sends results to the master as soon as possible.
  See :class:`RESetMapReduceWorker` for details.


Here is an example or how to deal with timeout::


    sage: from sage.parallel.map_reduce import RESetMPExample, AbortError
    sage: EX = RESetMPExample(maxl = 8)
    sage: try:
    ....:     res = EX.run(timeout=0.1)
    ....: except AbortError:
    ....:     print "Computation timeout"
    ....: else:
    ....:     print "Computation normally finished"
    ....:     res
    Computation timeout

The following should not timeout even on a very slow machine::

    sage: try:
    ....:     res = EX.run(timeout=60)
    ....: except AbortError:
    ....:     print "Computation Timeout"
    ....: else:
    ....:     print "Computation normally finished"
    ....:     res
    Computation normally finished
    40320*x^8 + 5040*x^7 + 720*x^6 + 120*x^5 + 24*x^4 + 6*x^3 + 2*x^2 + x + 1


As for ``reduce_locally``, one should not see any difference, except for speed
during normal usage. Most of the time the user should leave it to ``True``,
unless he sets up a mecanism to consume the partial results as soon as they
arrive. See :class:`RESetParallelIterator` and in particular the ``__iter__``
method for a example of consumer use.


.. _profiling:

Profiling
---------

It is possible the profile a map/reduce computation. First we create a
:class:`RESetMapReduce` object::

    sage: from sage.parallel.map_reduce import RESetMapReduce
    sage: S = RESetMapReduce(
    ....:     roots = [[]],
    ....:     children = lambda l: [l+[0], l+[1]] if len(l) <= 15 else [],
    ....:     map_function = lambda x : 1,
    ....:     reduce_function = lambda x,y: x+y,
    ....:     reduce_init = 0 )

The profiling is activated by the ``profile`` parameter. The value provided
should be a prefix (including a possible directory) for the profile dump::

    sage: prof = tmp_dir('RESetMR_profile')+'profcomp'
    sage: res = S.run(profile=prof) # random
    [RESetMapReduceWorker-1:58] (20:00:41.444) Profiling in /home/user/.sage/temp/mymachine.mysite/32414/RESetMR_profilewRCRAx/profcomp1 ...
    ...
    [RESetMapReduceWorker-1:57] (20:00:41.444) Profiling in /home/user/.sage/temp/mymachine.mysite/32414/RESetMR_profilewRCRAx/profcomp0 ...
    sage: res
    131071

In this example, the profile have been dumped in files such as
``profcomp0``. One can then load and print them as follows. See
:class:`profile.profile` for more details::

    sage: import cProfile, pstats
    sage: st = pstats.Stats(prof+'0')
    sage: st.strip_dirs().sort_stats('cumulative').print_stats() #random
    ...
       Ordered by: cumulative time

       ncalls  tottime  percall  cumtime  percall filename:lineno(function)
            1    0.023    0.023    0.432    0.432 map_reduce.py:1211(run_myself)
        11968    0.151    0.000    0.223    0.000 map_reduce.py:1292(walk_branch_locally)
    ...
    <pstats.Stats instance at 0x7fedea40c6c8>

.. SEEALSO::

    `The Python Profilers <https://docs.python.org/2/library/profile.html>`_
    for more detail on profiling in python.


.. _logging:

Logging
-------

The computation progress is logged through a :class:`logging.Logger` in
:obj:`sage.parallel.map_reduce.logger` together with :class:`logging.StreamHandler`
and a :class:`logging.Formatter`. They are currently configured to print
warning message on the console.

.. SEEALSO::

    `Logging facility for Python <https://docs.python.org/2/library/logging.html>`_
    for more detail on logging and log system configuration.

.. note::

    Calls to logger which involve printing the node are commented out in the
    code, because the printing (to a string) of the node can be very time
    consuming depending on the node and it happens before the decision whether
    the logger should record the string or drop it.


.. _protocol-description:

How does it work ?
------------------

The scheduling algorithm we use here is any adaptation of :wikipedia:`Work_stealing`:

    In a work stealing scheduler, each processor in a computer system has a
    queue of work items (computational tasks, threads) to perform. [...]. Each
    work items are initially put on the queue of the processor executing the
    work item. When a processor runs out of work, it looks at the queues of
    other processors and "steals" their work items. In effect, work stealing
    distributes the scheduling work over idle processors, and as long as all
    processors have work to do, no scheduling overhead occurs.

For communication we use Python's basic :mod:`multiprocessing` module. We
first describe the different actors and communications tools used by the
system. The work is done under the coordination of a **master** object (an
instance of :class:`RESetMapReduce`) by a bunch of **worker** objects
(instances of :class:`RESetMapReduceWorker`).

Each running map reduce instance work on a :class:`RecursivelyEnumeratedSet of
forest type<sage.combinat.backtrack.SearchForest>` called here `C` and is
coordinated by a :class:`RESetMapReduce` object called the **master**. The
master is in charge of lauching the work, gathering the results and cleaning
up at the end of the computation. It doesn't perform any computation
associated to the generation of the element `C` nor the computation of the
mapped function. It however occasionally perform a reduce, but most reducing
is by default done by the workers. Also thanks to the work-stealing algorithm,
the master is only involved in detecting the termination of the computation
but all the load balancing is done at the level of the worker.

Workers are instance of :class:`RESetMapReduceWorker`. They are responsible of
doing the actual computations: elements generation, mapping and reducing. They
are also responsible of the load balancing thanks to work-stealing.

Here is a description of the attribute of the **master** relevant to the
map-reduce protocol:

- ``master._results`` -- a :class:`~multiprocessing.queues.SimpleQueue` where
  the master gathers the results sent by the workers.
- ``master._active_tasks`` -- a :class:`~multiprocessing.Semaphore` recording
  the number of active task.  The work is done when it gets to 0.
- ``master._done`` -- a :class:`~multiprocessing.Lock` which ensures that
  shutdown is done only once.
- ``master._abort`` -- a :func:`~multiprocessing.Value` storing a shared
  :class:`ctypes.c_bool` which is ``True`` if the computation was aborted before
  all the worker runs out of work.
- ``master._workers`` -- a list of :class:`RESetMapReduceWorker` objects. Each worker is
  identified by its position in this list.

Each worker is a process (:class:`RESetMapReduceWorker` inherits from
:class:`~multiprocessing.Process`) which contains:

- ``worker._iproc`` -- the identifier of the worker that is its position in the
  master's list of workers
- ``worker._todo`` -- a :class:`collections.deque` storing of nodes of the
  worker. It is used as a stack by the worker. Thiefs steal from the bottom of
  this queue.
- ``worker._request`` -- a :class:`~multiprocessing.queues.SimpleQueue` storing
  steal request submitted to ``worker``.
- ``worker._read_task``, ``worker._write_task`` -- a
  :class:`~multiprocessing.queues.Pipe` used to transfert node during steal.
- ``worker._thief`` -- a :class:`~threading.Thread` which is in charge of stealing from
  ``worker._todo``.

Here is a schematic of the architecture:

.. _figure-map_reduce_arch:

.. figure:: ../../media/map_reduce_arch.png


How thefts are performed
------------------------

During normal time, that is when all worker are active) a worker ``W`` is
iterating though a loop inside
:meth:`RESetMapReduceWorker.walk_branch_locally`. Work nodes are taken from
and new nodes ``W._todo`` are appended to ``W._todo``. When a worker ``W`` is
running out of work, that is ``worker._todo`` is empty, then it tries to steal
some work (ie: a node) from another worker. This is performed in the
:meth:`RESetMapReduceWorker.steal` method.

From the point of view of ``W`` here is what happens:

- ``W`` signals to the master that it is idle :meth:`master._signal_task_done`;
- ``W`` chose a victim ``V`` at random;
- ``W`` sends a request to ``V`` : it puts its identifier into ``V._request``;
- ``W`` tries to read a node from ``W._read_task``. Then three things may happen:

  + a proper node is read. Then the theft was a success and ``W`` starts
    working locally on the received node.
  + ``None`` is received. This means that ``V`` was idle. Then ``W`` tries
    another victim.
  + ``AbortError`` is received. This means either that the computation was
    aborted or that it simply succeded and that no more work is required by
    ``W``. Therefore an ``AbortError`` exception is raised leading to ``W`` to
    shutdown.

We now describe the protocol on the victims side. Each worker process contains
a :class:`Thread` which we call ``T`` for thief which acts like some kinds of
Troyan horse during theft. It is normally blocked waiting for a steal request.

From the point of view of ``V`` and ``T``, here is what happens:

- during normal time ``T`` is blocked waiting on ``V._request``;
- upon steal request, ``T`` wakes up receiving the identification of ``W``;
- ``T`` signal to the master that a new task is starting by
  :meth:`master._signal_task_start`;
- Two things may happen depending if the queue ``V._todo`` is empty or not.
  Remark that due to the GIL, there is no parallel execution between the
  victim ``V`` and its thief tread ``T``.

  + If ``V._todo`` is empty, then ``None`` is answered on
    ``W._write_task``. The task is immediately signaled to end the the master
    through :meth:`master._signal_task_done`.
  + Otherwise, a node is removed from the bottom of ``V._todo``. The node is
    sent to ``W`` on ``W._write_task``. The task will be ended by ``W``, that
    is when finished working on the subtree rooted at the node, ``W`` will
    call :meth:`master._signal_task_done`.

The end of the computation
--------------------------

To detect when a computation is finished, we keep a synchronized integer which
count the number of active task. This is essentially a semaphore but semaphore
are broken on Darwin's OSes so we ship two implementations depending on the os
(see :class:`ActiveTaskCounter` and :class:`ActiveTaskCounterDarwin` and note
below).

When a worker finishes working on a task, it calls
:meth:`master._signal_task_done`. This decrease the task counter
``master._active_tasks``. When it reaches 0, it means that there are no more
nodes: the work is done. The worker executes :meth:`master._shutdown` which
sends ``AbortError`` on all :meth:`worker._request` and
:meth:`worker._write_task` Queues. Each worker or thief thread receiving such
a message raise the corresponding exception, stoping therefore its work. A
lock called ``master._done`` ensures that shutdown is only done once.

Finally, it is also possible to interrupt the computation before its ends
calling :meth:`master.abort()`. This is done by putting
``master._active_tasks`` to 0 and calling :meth:`master._shutdown`.

.. warning:: The MacOSX Semaphore bug

    Darwin's OSes do not correctly implement POSIX's semaphore semantic.
    Indeed, on this system, acquire may fail and return False not only because
    the semaphore is equal to zero but also **because someone else is trying to
    acquire** at the same time. This renders the usage of Semaphore impossible
    on MacOSX so that on this system we use a synchronized integer.


.. _examples:

Are there examples of classes ?
-------------------------------

Yes ! Here, there are:

- :class:`RESetMPExample` -- a simple basic example
- :class:`RESetParallelIterator` -- a more advanced example using non standard
  communication configuration.

Tests
-----

Generating series for sum of strictly decreassing list of integer smaller than
15::

    sage: y = var('y')
    sage: R = RESetMapReduce(
    ....:  roots = [([], 0, 0)] +[([i], i, i) for i in range(1,15)],
    ....:  children = lambda (list, sum, last):
    ....:      [(list + [i], sum + i, i) for i in range(1,last)],
    ....:  map_function = lambda (li, sum, unused): y**sum)
    sage: sg = R.run()
    sage: bool(sg == expand(prod((1+y^i) for i in range(1,15))))
    True


Classes and methods
-------------------
"""
from multiprocessing import Process, Value, Semaphore, Lock, cpu_count
from multiprocessing.queues import Pipe, SimpleQueue
from multiprocessing.sharedctypes import RawArray
from threading import Thread
from sage.sets.recursively_enumerated_set import RecursivelyEnumeratedSet # _generic
from sage.misc.lazy_attribute import lazy_attribute
import collections, copy, sys, random, ctypes


import logging
logger = logging.getLogger(__name__)
logger.__doc__ = """
A logger for :mod:`sage.parallel.map_reduce`

.. SEEALSO::

    `Logging facility for Python <https://docs.python.org/2/library/logging.html>`_
    for more detail on logging and log system configuration.
"""
logger.setLevel(logging.WARN)
#logger.setLevel(logging.INFO)
#logger.setLevel(logging.DEBUG)
ch = logging.StreamHandler()
ch.setLevel(logging.DEBUG)
formatter = logging.Formatter(
    '[%(processName)s-%(threadName)s] (%(asctime)s.%(msecs)03.f) %(message)s',
    datefmt='%H:%M:%S')
ch.setFormatter(formatter)
logger.addHandler(ch)



def proc_number(max_proc = None):
    r"""
    Computing the number of process used

    INPUT:

    - ``max_proc`` -- the maximum number of process used

    EXAMPLE::

        sage: from sage.parallel.map_reduce import proc_number
        sage: proc_number() # random
        8
        sage: proc_number(max_proc=1)
        1
        sage: proc_number(max_proc=2) in (1, 2)
        True
    """
    if max_proc is None:
        return max(cpu_count(), 1)
    else:
        return min(max_proc, max(cpu_count(), 1))


class AbortError(Exception):
    r"""
    Exception for aborting parallel computations

    This is used both as exception or as abort message

    TESTS::

        sage: from sage.parallel.map_reduce import AbortError
        sage: raise AbortError
        Traceback (most recent call last):
        ...
        AbortError
    """
    pass


class ActiveTaskCounterDarwin(object):
    r"""
    Handling the number of Active Tasks

    A class for handling the number of active task in distributed computation
    process. This is essentially a semaphore, but Darwin's OSes do not
    correctly implement POSIX's semaphore semantic. So we use a shared integer
    with a lock.
    """
    def __init__(self, task_number):
        r"""
        TESTS::

            sage: from sage.parallel.map_reduce import ActiveTaskCounterDarwin as ATC
            sage: t = ATC(4)
            sage: TestSuite(t).run(skip="_test_pickling", verbose=True)
        """
        self._active_tasks = Value(ctypes.c_int, task_number)
        self._lock = Lock()

    def __repr__(self):
        """
        TESTS::

            sage: from sage.parallel.map_reduce import ActiveTaskCounterDarwin as ATC
            sage: ATC(4)
            ActiveTaskCounter(value=4)
        """
        return "ActiveTaskCounter(value=%s)"%(self._active_tasks.value)

    def task_start(self):
        r"""
        Increment the task counter by one.

        OUTPUT:

        Calling :meth:`task_start` on a zero or negative counter returns 0,
        otherwise increment the counter and returns its value after the
        incrementation.

        EXAMPLES::

            sage: from sage.parallel.map_reduce import ActiveTaskCounterDarwin as ATC
            sage: c = ATC(4); c
            ActiveTaskCounter(value=4)
            sage: c.task_start()
            5
            sage: c
            ActiveTaskCounter(value=5)

        Calling :meth:`task_start` on a zero counter does nothing::

            sage: c = ATC(0)
            sage: c.task_start()
            0
            sage: c
            ActiveTaskCounter(value=0)
        """
        logger.debug("_signal_task_start called")
        with self._lock:
            # The following test is not necessary but is allows active thieves to
            # stop before receiving the poison pill.
            if self._active_tasks.value <= 0:
                return 0
            self._active_tasks.value += 1
            return self._active_tasks.value

    def task_done(self):
        r"""
        Decrement the task counter by one.

        OUTPUT:

        Calling :meth:`task_done` decrement the counter and returns its value
        after the decrementation.

        EXAMPLES::

            sage: from sage.parallel.map_reduce import ActiveTaskCounterDarwin as ATC
            sage: c = ATC(4); c
            ActiveTaskCounter(value=4)
            sage: c.task_done()
            3
            sage: c
            ActiveTaskCounter(value=3)

            sage: c = ATC(0)
            sage: c.task_done()
            -1
        """
        logger.debug("_signal_task_done called")
        with self._lock:
            self._active_tasks.value -= 1
            return self._active_tasks.value

    def abort(self):
        r"""
        Set the task counter to 0.

        EXAMPLES::

            sage: from sage.parallel.map_reduce import ActiveTaskCounterDarwin as ATC
            sage: c = ATC(4); c
            ActiveTaskCounter(value=4)
            sage: c.abort()
            sage: c
            ActiveTaskCounter(value=0)
        """
        with self._lock:
            self._active_tasks.value = 0


class ActiveTaskCounterPosix(object):
    r"""
    Handling the number of Active Tasks

    A class for handling the number of active task in distributed computation
    process. This is the standard implementation on POSIX compliant OSes. We
    essentially wrap a semaphore.

    .. note::

        A legitimate question is whether there is a need in keeping the two
        implementations. I ran the following experiment on my machine::

            S = RecursivelyEnumeratedSet( [[]],
                    lambda l: ([l[:i] + [len(l)] + l[i:] for i in range(len(l)+1)]
                      if len(l) < NNN else []),
                structure='forest', enumeration='depth')
            %time sp = S.map_reduce(lambda z: x**len(z)); sp

        For NNN = 10, averaging a dozen of runs, I got:

        - Posix complient implementation : 17.04 s
        - Darwin's implementation        : 18.26 s

        So there is a non negligible overhead. It will probably be worth if we
        tries to Cythonize the code. So I'm keeping both implementation.
    """
    def __init__(self, task_number):
        r"""
        TESTS::

            sage: from sage.parallel.map_reduce import ActiveTaskCounter as ATC
            sage: t = ATC(4)
            sage: TestSuite(t).run(skip="_test_pickling", verbose=True)
        """
        self._active_tasks = Semaphore(task_number)

    def __repr__(self):
        """
        TESTS::

            sage: from sage.parallel.map_reduce import ActiveTaskCounter as ATC
            sage: ATC(4)
            ActiveTaskCounter(value=4)
        """
        return "ActiveTaskCounter(value=%s)"%(self._active_tasks.get_value())

    def task_start(self):
        r"""
        Increment the task counter by one.

        OUTPUT:

        Calling :meth:`task_start` on a zero or negative counter returns 0,
        otherwise increment the counter and returns its value after the
        incrementation.

        EXAMPLES::

            sage: from sage.parallel.map_reduce import ActiveTaskCounter as ATC
            sage: c = ATC(4); c
            ActiveTaskCounter(value=4)
            sage: c.task_start()
            5
            sage: c
            ActiveTaskCounter(value=5)

        Calling :meth:`task_start` on a zero counter does nothing::

            sage: c = ATC(0)
            sage: c.task_start()
            0
            sage: c
            ActiveTaskCounter(value=0)
        """
        logger.debug("_signal_task_start called")
        # The following test is not necessary but is allows active thieves to
        # stop before receiving the poison pill.
        if self._active_tasks._semlock._is_zero():
            return 0
        self._active_tasks.release()
        return self._active_tasks.get_value()

    task_start.__doc__ = ActiveTaskCounterDarwin.task_start.__doc__

    def task_done(self):
        r"""
        Decrement the task counter by one.

        OUTPUT:

        Calling :meth:`task_done` decrement the counter and returns its value
        after the decrementation.

        EXAMPLES::

            sage: from sage.parallel.map_reduce import ActiveTaskCounter as ATC
            sage: c = ATC(4); c
            ActiveTaskCounter(value=4)
            sage: c.task_done()
            3
            sage: c
            ActiveTaskCounter(value=3)

            sage: c = ATC(0)
            sage: c.task_done()
            -1
        """
        logger.debug("_signal_task_done called")
        # We tests if the semaphore counting the number of active tasks is
        # becoming negative. This should not happen in normal
        # computations. However, in case of abort, we artificially put the
        # semaphore to 0 to stop the computation so it is needed.
        if not self._active_tasks.acquire(False):
            return -1
        return self._active_tasks.get_value()

    def abort(self):
        r"""
        Set the task counter to 0.

        EXAMPLES::

            sage: from sage.parallel.map_reduce import ActiveTaskCounter as ATC
            sage: c = ATC(4); c
            ActiveTaskCounter(value=4)
            sage: c.abort()
            sage: c
            ActiveTaskCounter(value=0)
        """
        while self._active_tasks.acquire(False):
            pass


ActiveTaskCounter = (ActiveTaskCounterDarwin if sys.platform == 'darwin'
                     else ActiveTaskCounterPosix)

# ActiveTaskCounter = ActiveTaskCounterDarwin # to debug DARWIN's implem



class RESetMapReduce(object):
    r"""
    Map-Reduce on recursively enumerated sets

    INPUT:

    Description of the set:

    - either ``forest=f`` -- where ``f`` is a
      :class:`RecursivelyEnumeratedSet of forest type<sage.combinat.backtrack.SearchForest>`

    - or a triple ``roots, children, post_process`` as follows

      - ``roots=r`` -- The root of the enumeration
      - ``children=c`` -- a function iterating through children node, given a parent nodes
      - ``post_process=p`` -- a post processing function

    The option ``post_process`` allows for customizing the nodes that
    are actually produced. Furthermore, if ``post_process(x)`` returns ``None``,
    then ``x`` won't be output at all.

    Decription of the map/reduce operation:

    - ``map_function=f`` -- (default to ``None``)
    - ``reduce_function=red`` -- (default to ``None``)
    - ``reduce_init=init`` -- (default to ``None``)

    .. seealso::

       :mod:`the Map/Reduce module <sage.parallel.map_reduce>` for
       details and examples.
    """
    def __init__(self, roots = None,
                 children = None,
                 post_process = None,
                 map_function = None,
                 reduce_function = None,
                 reduce_init = None,
                 forest = None):
        r"""
        TESTS::

            sage: from sage.parallel.map_reduce import RESetMapReduce
            sage: R = RESetMapReduce( [[]], lambda : [[]])
            sage: R
            <sage.parallel.map_reduce.RESetMapReduce object at 0x...>

        To silence the coverage checker::

            sage: TestSuite(R).run(skip=['_test_pickling'])
        """
        if forest is not None:
            if not all(x is None for x in (roots, children, post_process)):
                raise ValueError("forest arg is incompatible with roots, children and post_process")
            self._forest = forest
            self._roots = forest._roots
            self.children = forest.children
            if hasattr(forest, 'post_process'):
                self.post_process = forest.post_process
        else:
            if roots is not None: self._roots = roots
            if children is not None: self.children = children
            if post_process is not None: self.post_process = post_process
        if map_function is not None: self.map_function = map_function
        if reduce_function is not None: self.reduce_function = reduce_function
        if reduce_init is not None: self._reduce_init = reduce_init
        self._profile = None

    @lazy_attribute
    def _forest(self):
        r"""
        The forest underlying the map-reduce computation

        EXAMPLES::

            sage: from sage.parallel.map_reduce import RESetMPExample
            sage: EX = RESetMPExample()
            sage: f = EX._forest; f
            An enumerated set with a forest structure
            sage: f.an_element()
            []
        """
        return RecursivelyEnumeratedSet(
            self.roots(),
            self.children,
            post_process=self.post_process,
            structure='forest', enumeration='depth')


    def roots(self):
        r"""
        Return the roots of ``self``

        OUTPUT:

        an iterable of nodes

        .. note:: This should be overloaded in applications.

        EXAMPLES::

            sage: from sage.parallel.map_reduce import RESetMapReduce
            sage: S = RESetMapReduce(42)
            sage: S.roots()
            42
        """
        return self._roots

    def map_function(self, o):
        r"""
        Return the function mapped by ``self``

        INPUT:

        - ``o`` -- a node

        OUTPUT:

        By default ``1``.

        .. note:: This should be overloaded in applications.

        EXAMPLES::

            sage: from sage.parallel.map_reduce import RESetMapReduce
            sage: S = RESetMapReduce()
            sage: S.map_function(7)
            1
            sage: S = RESetMapReduce(map_function = lambda x: 3*x + 5)
            sage: S.map_function(7)
            26
         """
        return 1

    def reduce_function(self, a, b):
        r"""
        Return the reducer function for ``self``

        INPUT:

        - ``a``, ``b`` -- two value to be reduced

        OUTPUT:

        by default the sum of ``a`` and ``b``.

        .. note:: This should be overloaded in applications.

        EXAMPLES::

            sage: from sage.parallel.map_reduce import RESetMapReduce
            sage: S = RESetMapReduce()
            sage: S.reduce_function(4, 3)
            7
            sage: S = RESetMapReduce(reduce_function=lambda x,y: x*y)
            sage: S.reduce_function(4, 3)
            12
        """
        return a+b

    def post_process(self, a):
        r"""
        Return the post-processing function for ``self``

        INPUT: ``a`` -- a node

        By default, returns ``a`` itself

        .. note:: This should be overloaded in applications.

        EXAMPLES::

            sage: from sage.parallel.map_reduce import RESetMapReduce
            sage: S = RESetMapReduce()
            sage: S.post_process(4)
            4
            sage: S = RESetMapReduce(post_process=lambda x: x*x)
            sage: S.post_process(4)
            16
        """
        return a


    _reduce_init = 0

    def reduce_init(self):
        r"""
        Return the initial element for a reduction

        .. note:: This should be overloaded in applications.

        TESTS::

            sage: from sage.parallel.map_reduce import RESetMapReduce
            sage: S = RESetMapReduce()
            sage: S.reduce_init()
            0
            sage: S = RESetMapReduce(reduce_init = 2)
            sage: S.reduce_init()
            2
        """
        return copy.copy(self._reduce_init)


    def setup_workers(self, max_proc = None, reduce_locally=True):
        r"""
        Setup the communication channels

        INPUT:

        - ``mac_proc`` -- an integer: the maximum number of workers

        - ``reduce_locally`` -- whether the workers should reduce locally
          their work or sends results to the master as soon as possible.
          See :class:`RESetMapReduceWorker` for details.

        TESTS::

            sage: from sage.parallel.map_reduce import RESetMapReduce
            sage: S = RESetMapReduce()
            sage: S.setup_workers(2)
            sage: S._results
            <multiprocessing.queues.SimpleQueue object at 0x...>
            sage: len(S._workers)
            2
        """
        self._nprocess = proc_number(max_proc)
        self._results = SimpleQueue()
        self._active_tasks = ActiveTaskCounter(self._nprocess)
        self._done = Lock()
        self._abort = Value(ctypes.c_bool, False)
        sys.stdout.flush()
        sys.stderr.flush()
        self._workers = [RESetMapReduceWorker(self, i, reduce_locally)
                         for i in range(self._nprocess)]

    def start_workers(self):
        r"""
        Lauch the workers

        The worker should have been created using :meth:`setup_workers`.

        TESTS::

            sage: from sage.parallel.map_reduce import RESetMapReduce
            sage: S = RESetMapReduce(roots=[])
            sage: S.setup_workers(2)
            sage: S.start_workers()
            sage: all(w.is_alive() for w in S._workers)
            True

            sage: sleep(1)
            sage: all(not w.is_alive() for w in S._workers)
            True

        Cleanups::

            sage: S.finish()
        """
        if self._nprocess == 0:
            raise ValueError("No process connected")
        logger.info("Starting work with %s processes", self._nprocess)
        logger.debug("Distributing tasks")
        for i, task in enumerate(self.roots()):
            self._workers[i % len(self._workers)]._todo.append(task)
        logger.debug("Starting processes")
        sys.stdout.flush()
        sys.stderr.flush()
        for w in self._workers: w.start()

    def get_results(self):
        r"""
        Get the results from the queue

        OUTPUT:

        the reduction of the results of all the workers, that is the result of
        the map/reduce computation.

        EXAMPLES::

            sage: from sage.parallel.map_reduce import RESetMapReduce
            sage: S = RESetMapReduce()
            sage: S.setup_workers(2)
            sage: for v in [1, 2, None, 3, None]: S._results.put(v)
            sage: S.get_results()
            6

        Cleanups::

            sage: del S._results, S._active_tasks, S._done, S._workers
        """
        res = self.reduce_init()
        active_proc = self._nprocess
        while active_proc > 0:
            newres = self._results.get()
            if newres is not None:
                logger.debug("Got one result")
                res = self.reduce_function(res, newres)
            else:
                active_proc -= 1
        return res


    def finish(self):
        r"""
        Destroys the worker and all the communication objects.

        Also gathers the communication statistics before destroying the workers.

        TESTS::

            sage: from sage.parallel.map_reduce import RESetMPExample
            sage: S = RESetMPExample(maxl=5)
            sage: S.setup_workers(2) # indirect doctest
            sage: S._workers[0]._todo.append([])
            sage: for w in S._workers: w.start()
            sage: _ = S.get_results()
            sage: S._shutdown()
            sage: S.print_communication_statistics()
            Traceback (most recent call last):
            ...
            AttributeError: 'RESetMPExample' object has no attribute '_stats'

            sage: S.finish()

            sage: S.print_communication_statistics()
            #proc: ...
            ...

            sage: _ = S.run() # Cleanup

        .. seealso:: :meth:`print_communication_statistics`
        """
        self._abort = self._abort.value
        if not self._abort:
            logger.debug("Joining worker processes...")
            for worker in self._workers:
                logger.debug("Joining %s"%worker.name)
                worker.join()
            logger.debug("Joining done")
        else:
            logger.debug("Killing worker processes...")
            for worker in self._workers:
                logger.debug("Terminating %s"%worker.name)
                worker.terminate()
            logger.debug("Killing done")

        del self._results, self._active_tasks, self._done
        self._get_stats()
        del self._workers


    def abort(self):
        r"""
        Abort the current parallel computation

        EXAMPLES::

            sage: from sage.parallel.map_reduce import RESetParallelIterator
            sage: S = RESetParallelIterator( [[]],
            ....:   lambda l: [l+[0], l+[1]] if len(l) < 17 else [])
            sage: it = iter(S)
            sage: it.next()
            []
            sage: S.abort()
            sage: hasattr(S, 'work_queue')
            False

        Cleanups::

            sage: S.finish()
        """
        logger.info("Abort called")
        self._abort.value = True
        self._active_tasks.abort()
        self._shutdown()

    def _shutdown(self):
        r"""
        Called to shutdown the workers

        Sends a poison pill to all workers and their thief thread.

        EXAMPLES::

            sage: from sage.parallel.map_reduce import RESetParallelIterator
            sage: S = RESetParallelIterator( [[]],
            ....:   lambda l: [l+[0], l+[1]] if len(l) < 20 else [])
            sage: S.setup_workers(2)
            sage: for w in S._workers: w.start()
            sage: S._shutdown()

        Cleanups::

            sage: S.finish()
        """
        if self._done.acquire(False):
            logger.debug("***************** FINISHED ******************")
            logger.debug("Sending poison pills")
            for worker in self._workers:
                worker._request.put(AbortError)
            for worker in self._workers:
                worker._write_task.send(AbortError)

    def _signal_task_start(self):
        r"""
        Signal a starting task

        Used by the worker to signal that a new task is starting. As soon as
        there are no more active task, the work is done, in which case an
        :exc:`AbortError` is raised.

        EXAMPLES::

            sage: from sage.parallel.map_reduce import RESetParallelIterator
            sage: S = RESetParallelIterator( [[]],
            ....:   lambda l: [l+[0], l+[1]] if len(l) < 20 else [])
            sage: S.setup_workers(2)
            sage: S._active_tasks
            ActiveTaskCounter(value=2)

            sage: S._signal_task_start()
            sage: S._active_tasks
            ActiveTaskCounter(value=3)

        Signaling one time too many raise a ``AbortError``::

            sage: S._signal_task_done()
            sage: S._signal_task_done()
            sage: S._signal_task_done()
            Traceback (most recent call last):
            ...
            AbortError
        """
        if self._active_tasks.task_start() == 0:
            raise AbortError

    def _signal_task_done(self):
        r"""
        Signal a done task

        Used by the worker to signal that a task is done. As soon as
        there are no more active task, the work is done, in which case an
        :exc:`AbortError` is raised.

        EXAMPLES::

            sage: from sage.parallel.map_reduce import RESetParallelIterator
            sage: S = RESetParallelIterator( [[]],
            ....:   lambda l: [l+[0], l+[1]] if len(l) < 20 else [])
            sage: S.setup_workers(2)
            sage: S._active_tasks
            ActiveTaskCounter(value=2)

            sage: S._signal_task_done()
            sage: S._active_tasks
            ActiveTaskCounter(value=1)

            sage: S._signal_task_done()
            Traceback (most recent call last):
            ...
            AbortError

        Cleanups::

            sage: del S._results, S._active_tasks, S._done, S._workers
        """
        # We tests if the semaphore counting the number of active tasks is
        # becoming negative. This should not happen in normal
        # computations. However, in case of abort, we artificially put the
        # semaphore to 0 to stop the computation so that it is needed.
        if self._active_tasks.task_done() <= 0:
            logger.debug("raising AbortError")
            self._shutdown()
            raise AbortError

    def random_worker(self):
        r"""
        Returns a random workers

        OUTPUT:

        A worker for ``self`` chosed at random

        EXAMPLES::

            sage: from sage.parallel.map_reduce import RESetMPExample, RESetMapReduceWorker
            sage: from threading import Thread
            sage: EX = RESetMPExample(maxl=6)
            sage: EX.setup_workers(2)
            sage: EX.random_worker()
            <RESetMapReduceWorker(RESetMapReduceWorker-..., initial)>
            sage: EX.random_worker() in EX._workers
            True

        Cleanups::

            sage: del EX._results, EX._active_tasks, EX._done, EX._workers
        """
        victim = random.randint(0, len(self._workers)-1)
        return self._workers[victim]

    def run(self,
            max_proc = None,
            reduce_locally = True,
            timeout=None,
            profile=None):
        r"""
        Run the computations

        INPUT:

        - ``max_proc`` -- maximum number of process used.
          default: number of processor on the machine
        - ``reduce_locally`` -- See :class:`RESetMapReduceWorker` (default: ``True``)
        - ``timeout`` -- a timeout on the computation (default: ``None``)
        - ``profile`` -- directory/filename prefix for profiling, or ``None``
          for no profiling (default: ``None``)

        OUTPUT:

        the result of the map/reduce computation or an exception
        :exc:`AbortError` if the computation was interrupted or timeout.

        EXAMPLES::

            sage: from sage.parallel.map_reduce import RESetMPExample
            sage: EX = RESetMPExample(maxl = 8)
            sage: EX.run()
            40320*x^8 + 5040*x^7 + 720*x^6 + 120*x^5 + 24*x^4 + 6*x^3 + 2*x^2 + x + 1

        Here is an example or how to deal with timeout::

            sage: from sage.parallel.map_reduce import AbortError
            sage: try:
            ....:     res = EX.run(timeout=0.1)
            ....: except AbortError:
            ....:     print "Computation timeout"
            ....: else:
            ....:     print "Computation normally finished"
            ....:     res
            Computation timeout

        The following should not timeout even on a very slow machine::

            sage: from sage.parallel.map_reduce import AbortError
            sage: try:
            ....:     res = EX.run(timeout=60)
            ....: except AbortError:
            ....:     print "Computation Timeout"
            ....: else:
            ....:     print "Computation normally finished"
            ....:     res
            Computation normally finished
            40320*x^8 + 5040*x^7 + 720*x^6 + 120*x^5 + 24*x^4 + 6*x^3 + 2*x^2 + x + 1
        """
        self._profile=profile
        self.setup_workers(max_proc, reduce_locally)
        self.start_workers()
        if timeout is not None:
            from threading import Timer
            self._timer = Timer(timeout, self.abort)
            self._timer.start()
        self.result = self.get_results()
        self.finish()
        if timeout is not None:
            self._timer.cancel()
        logger.info("Returning")
        if self._abort:
            raise AbortError
        else:
            return self.result

    def _get_stats(self):
        r"""
        Gather the communication statistics and the end of a run

        EXAMPLES::

            sage: from sage.parallel.map_reduce import RESetMPExample
            sage: S = RESetMPExample(maxl=6)
            sage: S.run() # indirect doctest
            720*x^6 + 120*x^5 + 24*x^4 + 6*x^3 + 2*x^2 + x + 1
        """
        res = []
        for i in range(self._nprocess):
            res.append(tuple(self._workers[i]._stats))
        self._stats = res

    def print_communication_statistics(self, blocksize = 16):
        r"""
        Print the communication statistics in a nice way

        EXAMPLES::

            sage: from sage.parallel.map_reduce import RESetMPExample
            sage: S = RESetMPExample(maxl=6)
            sage: S.run()
            720*x^6 + 120*x^5 + 24*x^4 + 6*x^3 + 2*x^2 + x + 1

            sage: S.print_communication_statistics()    # random
            #proc:        0    1    2    3    4    5    6    7
            reqs sent:    5    2    3   11   21   19    1    0
            reqs rcvs:   10   10    9    5    1   11    9    2
            - thefs:      1    0    0    0    0    0    0    0
            + thefs:      0    0    1    0    0    0    0    0
        """
        res = [""] # classical trick to have a local variable shared with the
        # local function (see e.g:
        # http://stackoverflow.com/questions/2609518/python-nested-function-scopes).
        def pstat(name, start, end, ist):
            res[0] += "\n" + name
            res[0] += " ".join(
                "%4i"%(self._stats[i][ist]) for i in range(start, end))
        for start in range(0, self._nprocess, blocksize):
            end = min(start+blocksize, self._nprocess)
            res[0] = "#proc:     "+" ".join("%4i"%(i) for i in range(start, end))
            pstat("reqs sent: ", start, end, 0)
            pstat("reqs rcvs: ", start, end, 1)
            pstat("- thefs:   ", start, end, 2)
            pstat("+ thefs:   ", start, end, 3)
        print res[0]

    def run_serial(self):
        r"""
        Serial run of the computation (mostly for tests)

        EXAMPLES::

            sage: from sage.parallel.map_reduce import RESetMPExample
            sage: EX = RESetMPExample(maxl = 4)
            sage: EX.run_serial()
            24*x^4 + 6*x^3 + 2*x^2 + x + 1
        """
        import functools
        return functools.reduce(self.reduce_function,
                                (self.map_function(x) for x in self._forest),
                                self.reduce_init())


class RESetMapReduceWorker(Process):
    """
    Worker for generate-map-reduce

    This shouldn't be called directly, but instead created by
    :meth:`RESetMapReduce.setup_workers`.

    INPUT:

    - ``mapred`` -- the instance of :class:`RESetMapReduce` for which
      this process is working.

    - ``iproc`` -- the id of this worker.

    - ``reduce_locally`` -- when reducing the results. Three possible values
      are supported:

      * ``True`` -- means the reducing work is done all locally, the result is
        only sent back at the end of the work. This ensure the lowest level of
        communication.

      * ``False`` -- results are sent back after each finished branches, when
        the process is asking for more work.
    """
    def __init__(self, mapred, iproc, reduce_locally):
        r"""
        TESTS::

            sage: from sage.parallel.map_reduce import RESetMPExample, RESetMapReduceWorker
            sage: EX = RESetMPExample()
            sage: RESetMapReduceWorker(EX, 200, True)
            <RESetMapReduceWorker(RESetMapReduceWorker-..., initial)>
        """
        Process.__init__(self)
        self._iproc = iproc
        self._todo = collections.deque()
        self._request = SimpleQueue()  # Faster than Queue
        # currently this is not possible to have to simultaneous read or write
        # on the following Pipe. So there is no need to have a queue.
        self._read_task, self._write_task = Pipe(duplex=False)
        self._mapred = mapred
        self._stats  =  RawArray('i', 4)
        self._reduce_locally = reduce_locally

    def _thief(self):
        r"""
        The thief thread of a worker process

        EXAMPLES::

            sage: from sage.parallel.map_reduce import RESetMPExample, RESetMapReduceWorker
            sage: from threading import Thread
            sage: EX = RESetMPExample(maxl=6)
            sage: EX.setup_workers(2)

            sage: w0, w1 = EX._workers
            sage: w0._todo.append(42)
            sage: thief0 = Thread(target = w0._thief, name="Thief")
            sage: thief0.start()

            sage: w1.steal()
            42
            sage: w0._todo
            deque([])
        """
        logger.debug("Thief started")
        reqs = 0
        thefts = 0

        try:
            for ireq in iter(self._request.get, AbortError):
                reqs +=1
                target = self._mapred._workers[ireq]
                logger.debug("Got a Steal request from %s"%target.name)
                self._mapred._signal_task_start()
                try:
                    work = self._todo.popleft()
                except IndexError:
                    target._write_task.send(None)
                    logger.debug("Failed Steal %s"%target.name)
                    self._mapred._signal_task_done()
                else:
                    target._write_task.send(work)
                    logger.debug("Succesful Steal %s"%target.name)
                    thefts += 1
        except AbortError:
            logger.debug("Thief aborted")
        else:
            logger.debug("Thief received poison pill")
        if self._mapred._abort.value:  # Computation was aborted
            self._todo.clear()
        else: # Check that there is no remaining work
            assert len(self._todo) == 0, "Bad stop the result may be wrong"

        self._stats[1] = reqs
        self._stats[2] = thefts
        logger.debug("Thief Exiting")

    def steal(self):
        r"""
        Steal some node from another worker

        OUTPUT:

        a node stolen from another worker choosed at random

        EXAMPLES::

            sage: from sage.parallel.map_reduce import RESetMPExample, RESetMapReduceWorker
            sage: from threading import Thread
            sage: EX = RESetMPExample(maxl=6)
            sage: EX.setup_workers(2)

            sage: w0, w1 = EX._workers
            sage: w0._todo.append(42)
            sage: thief0 = Thread(target = w0._thief, name="Thief")
            sage: thief0.start()

            sage: w1.steal()
            42
        """
        self._mapred._signal_task_done()
        node = None
        while node is None:
            victim = self._mapred.random_worker()
            if victim is not self:
                logger.debug("Trying to steal from %s"%(victim.name))
                victim._request.put(self._iproc)
                self._stats[0] += 1
                logger.debug("waiting from steal answer from %s"%(victim.name))
                node = self._read_task.recv()
                # logger.debug("Request answer: %s"%(node,))
                if node is AbortError:
                    raise AbortError
        # logger.debug("Received a stolen node: %s"%(node,))
        self._stats[3] += 1
        return node

    def run(self):
        r"""
        The main function executed by the worker

        Calls :meth:`run_myself` after possibly setting up parallel profiling.

        EXAMPLES::

            sage: from sage.parallel.map_reduce import RESetMPExample, RESetMapReduceWorker
            sage: EX = RESetMPExample(maxl=6)
            sage: EX.setup_workers(1)

            sage: w = EX._workers[0]
            sage: w._todo.append(EX.roots()[0])

            sage: w.run()
            sage: sleep(1)
            sage: w._todo.append(None)

            sage: EX.get_results()
            720*x^6 + 120*x^5 + 24*x^4 + 6*x^3 + 2*x^2 + x + 1

        Cleanups::

            sage: del EX._results, EX._active_tasks, EX._done, EX._workers
        """
        profile = self._mapred._profile
        if profile is not None:
            from multiprocessing import current_process
            import cProfile
            PROFILER = cProfile.Profile()
            PROFILER.runcall(self.run_myself)

            output = profile + str(self._iproc)
            logger.warn("Profiling in %s ..."%output)
            PROFILER.dump_stats(output)
        else:
            self.run_myself()

    def run_myself(self):
        r"""
        The main function executed by the worker

        EXAMPLES::

            sage: from sage.parallel.map_reduce import RESetMPExample, RESetMapReduceWorker
            sage: EX = RESetMPExample(maxl=6)
            sage: EX.setup_workers(1)

            sage: w = EX._workers[0]
            sage: w._todo.append(EX.roots()[0])
            sage: w.run_myself()

            sage: sleep(1)
            sage: w._todo.append(None)

            sage: EX.get_results()
            720*x^6 + 120*x^5 + 24*x^4 + 6*x^3 + 2*x^2 + x + 1

        Cleanups::

            sage: del EX._results, EX._active_tasks, EX._done, EX._workers
        """
        logger.debug("Started")
        mapred = self._mapred
        reduce_init = mapred.reduce_init
        results = mapred._results

        self._stats[0] = 0
        self._stats[3] = 0
        logger.debug("Launching thief")
        self._thief = Thread(target = self._thief, name="Thief")
        self._thief.start()
        self._res = reduce_init()

        try:
            while True:
                try:
                    node = self._todo.pop()
                except IndexError:
                    node = self.steal()
                self.walk_branch_locally(node)
                if not self._reduce_locally:
                    self.send_partial_result()
        except AbortError:
            logger.debug("Worker Done !")
            results.put(self._res)
        results.put(None)
        self._thief.join()
        del self._request
        self._read_task.close()
        self._write_task.close()
        del self._read_task, self._write_task
        del self._mapred
        del self._stats
        logger.debug("Exiting")

    def send_partial_result(self):
        r"""
        Send results to the MapReduce process

        Send the result stored in ``self._res`` to the master an reinitialize it to
        ``master.reduce_init``.

        EXAMPLES::

            sage: from sage.parallel.map_reduce import RESetMPExample, RESetMapReduceWorker
            sage: EX = RESetMPExample(maxl=4)
            sage: EX.setup_workers(1)
            sage: w = EX._workers[0]
            sage: w._res = 4
            sage: w.send_partial_result()
            sage: w._res
            0
            sage: EX._results.get()
            4
        """
        self._mapred._results.put(self._res)
        self._res = self._mapred.reduce_init()

    def walk_branch_locally(self, node):
        r"""
        Work locally

        Performs the map/reduce computation on the subtrees rooted at ``node``.

        INPUT:

        - ``node`` -- the root of the subtree explored.

        OUTPUT:

        nothing, the result are stored in ``self._res``

        This is where the actual work is performed.

        EXAMPLES::

            sage: from sage.parallel.map_reduce import RESetMPExample, RESetMapReduceWorker
            sage: EX = RESetMPExample(maxl=4)
            sage: w = RESetMapReduceWorker(EX, 0, True)
            sage: def sync(): pass
            sage: w.synchronize = sync
            sage: w._res = 0

            sage: w.walk_branch_locally([])
            sage: w._res
            x^4 + x^3 + x^2 + x + 1

            sage: w.walk_branch_locally(w._todo.pop())
            sage: w._res
            2*x^4 + x^3 + x^2 + x + 1

            sage: while True: w.walk_branch_locally(w._todo.pop())
            Traceback (most recent call last):
            ...
            IndexError: pop from an empty deque
            sage: w._res
            24*x^4 + 6*x^3 + 2*x^2 + x + 1
        """
        mapred = self._mapred
        children = mapred.children
        post_process = mapred.post_process
        fun  = mapred.map_function
        reduc = mapred.reduce_function

        # logger.debug("Working on %s..."%(node,))
        while True:
            res = post_process(node)
            if res is not None:
                self._res = reduc(self._res, fun(res))
            newnodes = iter(children(node))
            try:
                node = newnodes.next()
            except StopIteration:
                return
            self._todo.extend(newnodes)


class RESetMPExample(RESetMapReduce):
    r"""
    An example of map reduce class

    INPUT:

    - ``maxl`` -- the maximum size of permutations generated (default to `9`).

    This compute the generating series of permutations counted by their size
    upto size ``maxl``.

    EXAMPLE::

        sage: from sage.parallel.map_reduce import RESetMPExample
        sage: EX = RESetMPExample()
        sage: EX.run()
        362880*x^9 + 40320*x^8 + 5040*x^7 + 720*x^6 + 120*x^5 + 24*x^4 + 6*x^3 + 2*x^2 + x + 1

    .. seealso:: This is an example of :class:`RESetMapReduce`

    """
    def __init__(self, maxl = 9):
        r"""
        TESTS::

            sage: from sage.parallel.map_reduce import RESetMPExample
            sage: RESetMPExample()
            <sage.parallel.map_reduce.RESetMPExample object at 0x...>
        """
        RESetMapReduce.__init__(self)
        from sage.calculus.var import var
        self.x = var('x')
        self.maxl = maxl

    def roots(self):
        r"""
        Return the empty permutation

        EXAMPLE::

            sage: from sage.parallel.map_reduce import RESetMPExample
            sage: RESetMPExample().roots()
            [[]]
        """
        return [[]]

    def children(self, l):
        r"""
        Return the chidren of the permutation `l`

        INPUT:

        - ``l`` -- a list containing a permutation

        OUTPUT:

        the lists of ``len(l)`` inserted at all possible positions into ``l``

        EXAMPLE::

            sage: from sage.parallel.map_reduce import RESetMPExample
            sage: RESetMPExample().children([1,0])
            [[2, 1, 0], [1, 2, 0], [1, 0, 2]]
        """
        return [ l[:i] + [len(l)] + l[i:]
                 for i in range(len(l)+1) ] if len(l) < self.maxl else []

    def map_function(self, l):
        r"""
        The monomial associated to the permutation `l`

        INPUT:

        - ``l`` -- a list containing a permutation

        OUTPUT:

        ``x^len(l)``.

        EXAMPLE::

            sage: from sage.parallel.map_reduce import RESetMPExample
            sage: RESetMPExample().map_function([1,0])
            x^2
        """
        return self.x**len(l)


class RESetParallelIterator(RESetMapReduce):
    r"""
    A parallel iterator for recursively enumerated sets

    This demonstrate how to use :class:`RESetMapReduce` to get an iterator on
    a recursively enumerated sets for which the computations are done in
    parallel.

    EXAMPLE::

        sage: from sage.parallel.map_reduce import RESetParallelIterator
        sage: S = RESetParallelIterator( [[]],
        ....:   lambda l: [l+[0], l+[1]] if len(l) < 15 else [])
        sage: sum(1 for _ in S)
        65535
    """
    def map_function(self, z):
        r"""
        Return a singleton tuple

        INPUT: ``z`` -- a node

        OUPUT: ``(z, )``

        EXAMPLE::

            sage: from sage.parallel.map_reduce import RESetParallelIterator
            sage: S = RESetParallelIterator( [[]],
            ....:   lambda l: [l+[0], l+[1]] if len(l) < 15 else [])
            sage: S.map_function([1, 0])
            ([1, 0],)
        """
        return (z,)

    reduce_init = tuple

    def __iter__(self):
        r"""
        EXAMPLE::

            sage: from sage.parallel.map_reduce import RESetParallelIterator
            sage: S = RESetParallelIterator( [[]],
            ....:   lambda l: [l+[0], l+[1]] if len(l) < 15 else [])
            sage: it = iter(S)
            sage: it.next() # random
            [1, 1, 0]
            sage: it.next() # random
            [1, 1, 0, 1]
            sage: sum(1 for _ in it)
            65533
        """
        self.setup_workers(reduce_locally=False)
        self.start_workers()
        active_proc = self._nprocess
        while True:
            newres = self._results.get()
            if newres is not None:
                logger.debug("Got some results")
                for r in newres:
                    yield r
            else:
                active_proc -= 1
                if active_proc == 0:
                    break
        self.finish()

