META
----

Last Updated: 12/14/06
Contact Info: yqiang ATNOSPAM gmail.com


INTRODUCTION
------------

Distributed SAGE is a framework that allows one to do distributed
computing from within SAGE. It includes a server, client and workers as
well as a set of classes that one can subclass from to write distributed
computation jobs.


REQUIREMENTS
------------

* [**SAGE 2.6**](http://www.sagemath.org)
* [**OpenSSH**](http://www.openssh.org)


OVERVIEW
--------

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
-----------

1.  Launch sage

2.  Run 'dsage.setup()'. This is a one-stop-shop to get **DSAGE** running.
    It will create the databases and set up a private/public key for
    authentication

3.  Launch a server, monitor and get a connection to the server:

    `sage: D = dsage.start_all()`

    This will start 2 workers by default.

4.  To do a computation, use D just like any other SAGE interface. For
    example:

    `sage: j = D('2+2')`
    `sage: j.wait()`
    `sage: j`
    `4`

    Explanation:
    `D('2+2')` returns a JobWrapper object which is how one accesses the
    results of a computation.

    `j.wait()` will block until the job is finished

    j calls `__repr__` on the job which in this case prints the result

    You can put **any valid** SAGE code in the quotes and it will be
    converted into a JobWrapper object and performed by a worker.


FILES
-----

All **DSAGE** related files are in `$HOME/.sage/dsage`.

* `$HOME/.sage/dsage/db` contains the database
* `$HOME/.sage/dsage/tmp_worker_files` contains the jobs workers process
  (currently these are never cleaned up, but in the future they might be)


SETUP
-----

Setting up DSAGE is done through the configuration utility `dsage.setup()` and
the related utilities `dsage.setup_client()`, `dsage.setup_worker()` and
`dsage.setup_server`.


CONFIGURATION
-------------

As of **SAGE 2.6**, **DSAGE** does not use configuration files anymore.
Everything can be configured on the fly by passing in parameters to the
various **DSAGE** methods. That means that if you want a permanent
configuration, you should write a **SAGE** script that hardcodes the
parameters you wish to use.
To find out the parameters that are available, you can simply look at the
method definitions of `dsage.server()`, `dsage.worker()`, and `DSage()`.


USAGE (NOT ACCURATE)
-----

1.  Client
2.  Monitor
3.  Server


ADVANCED USAGE
--------------

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
----
Send bug reports to yqiang atNOSPAM gmail.com

AUTHORS
-------
[Yi Qiang](http://www.yiqiang.org)
