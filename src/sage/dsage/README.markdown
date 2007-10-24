META
====

Last Updated: 10/09/07
Contact Info: yqiang ATNOSPAM gmail.com


INTRODUCTION
============

Distributed SAGE is a framework that allows one to do distributed
computing from within SAGE. It includes a server, client and workers as
well as a set of classes that one can subclass from to write distributed
computation jobs.


REQUIREMENTS
============

* [**SAGE 2.6**](http://www.sagemath.org)
* [**OpenSSH**](http://www.openssh.org)


OVERVIEW
========

There are 3 distinct parts of Distributed SAGE:

-   Server
    The server is responsible for job distribution, submission and
    collection.

-   Client
    The client is responsible for submitting new jobs to the server
    and collecting the results.

-   Monitor
    The monitor controls workers which do the actual computation of
    jobs. Sometimes monitors are referred to just as 'workers'.


QUICK-START
===========

1.  Launch sage
2.  Run

    `sage: dsage.setup()`

    For a **really** quick start, just hit ENTER on all questions.
    This will create all the necessary supporting files to get **DSAGE**
    running. It will create the databases, set up a private/public key for
    authentication and create a SSL certificate for the server.
3.  Launch a server, monitor and get a connection to the server:

    `sage: D = dsage.start_all()`

    This will start 2 workers by default.  You can change it by passing in the
    `workers=N` argument where `N` is the number of workers you want.
4.  To do a computation, use D just like any other SAGE interface. For
    example:

    `sage: j = D('2+2')`
    `sage: j.wait()`
    `sage: j`
    `4`

    Explanation:
    Line 1 returns a JobWrapper object which is how one accesses the
    results of a computation.

    Line 2 `j.wait()` will block until the job is finished

    Line 3 calls `__repr__` on the job which in this case prints the result

    You can put **any valid** SAGE code in the quotes and it will be
    converted into a JobWrapper object and performed by a worker.


EXAMPLES
========

All examples below assume that you already have a connection to the DSAGE
server and have workers up and running. You can get to that stage by doing the
following from SAGE:

    sage: D = dsage.start_all()

Example 1:
----------

In this example we have written some code in a separate file, called
**dsage.sage** and would like to send that file as a job to a worker. The
contents of the file are simply:

    DSAGE_RESULT = factor(randint(10^20, 10^30))

Our session looks like this:

    sage: j = D.eval_file('dsage.sage', 'test_job')
    sage: j
    No output.
    sage: j.result
    2 * 5 * 7 * 23 * 743 * 877 * 1213 * 462270878024426767


Example 2:
----------
The above example is not very useful because rarely do you want to send just
**ONE** job to a worker. A more likely situation is if you want to factor a
list of 1000 integers over one hundred workers. We are going to use list
comprehension for the following example:


INCREASING YOUR COMPUTING POWER
===============================

At some point you will want to add more workers to your worker pool. To do
this is very easy, since authentication is done with public keys. You can find
your public/private key pair in:

    $HOME/.sage/dsage/dsage_key
    $HOME/.sage/dsage/dsage_key.pub

We have a DSAGE server *A* and a DSAGE worker *B* which are on separate
machines.

Copy your public/private key pair from *A* to *B* using the same
path. Now, in a session on *B* issue the following command to connect to *A*:

    sage: dsage.worker(server=**hostname of *A***, port=**port of server A**)

Alternatively, you can also connect as an **anonymous** worker by issuing the following command:

    sage: dsage.worker(server=**hostname of *A***, port=**port of server A**,
                       anonymous=True)

Please note that some jobs will not be processed by anonymous workers because
of the security setting of the job, which will be explained in the Jobs
sections.

FILES
=====

All **DSAGE** related files are in `$HOME/.sage/dsage`.

* `$HOME/.sage/dsage/db` contains the database
* `$HOME/.sage/dsage/tmp_worker_files` contains the jobs workers process
  (currently these are never cleaned up, but in the future they might be)


SETUP
=====

Setting up DSAGE is done through the configuration utility `dsage.setup()` and
the related utilities `dsage.setup_client()`, `dsage.setup_worker()` and
`dsage.setup_server`.


CONFIGURATION
=============

As of **SAGE 2.6**, **DSAGE** does not use configuration files anymore.
Everything can be configured on the fly by passing in parameters to the
various **DSAGE** methods. That means that if you want a permanent
configuration, you should write a **SAGE** script that hardcodes the
parameters you wish to use.
To find out the parameters that are available, you can simply look at the
method definitions of `dsage.server()`, `dsage.worker()`, and `DSage()`.


USAGE (NOT ACCURATE)
=====

1.  Client
2.  Monitor
3.  Server


ADVANCED USAGE
==============

More complicated usage of **DSAGE** will require the user to write some code.
We have tried to make this as simple as possible by providing most of the
functionality in a super class called DistributedFunction which is located in
`sage.dsage.dist_functions.DistributedFunction`.
There are several examples that illustrate how this is done.  You can find
these examples in

- `sage.dsage.dist_functions.dist_function`
- `sage.dsage.dist_functions.dist_povray`
- `sage.dsage.dist_functions.dist_factor`


BUGS
====

Send bug reports to yqiang atNOSPAM gmail.com

AUTHORS
=======

[Yi Qiang](http://www.yiqiang.org)
