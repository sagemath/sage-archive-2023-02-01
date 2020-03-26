.. nodoctest

.. highlight:: shell-session

.. _chapter-portability_testing:

=============================
Testing on multiple platforms
=============================

Sage is intended to build and run on a variety of platforms,
including all major Linux distributions, as well as MacOS, and
Windows (with Cygwin).

There is considerable variation among these platforms.
To ensure that Sage continues to build correctly on users'
machines, it is crucial to test changes to Sage, in particular
when external packages are added or upgraded, on a wide
spectrum of platforms.

Sage patchbots
==============

The `Sage patchbots <https://wiki.sagemath.org/patchbot>`_ will automatically test your Trac ticket
by attempting an incremental build of Sage and running doctests.

Sage buildbots
==============

The `Sage Release buildbot <https://wiki.sagemath.org/buildbot>`_
builds entire tarballs (e.g., all the development releases) on a
variety of machines.

Developers' and users' tests on sage-release
============================================

Sage developers and users are encouraged to contribute to testing
releases that are announced on
`Sage Release <https://groups.google.com/forum/#!forum/sage-release>`_ on their machines
and to report test results (success and failures) by responding to the
announcements.

Testing Sage on a different platform using Docker
=================================================

`Docker <https://www.docker.com>`_ is a popular virtualization software,
running Linux operating system images in containers on a shared Linux
kernel.  Using Docker Desktop for Mac and Windows, the containers can
also be run on these platforms.

All major Linux distributions provide ready-to-use Docker images,
which are published via `Docker Hub <https://hub.docker.com>`_.  For example, to run
the current stable (LTS) version of Ubuntu interactively, you can use
the shell command::

  [mkoeppe@sage sage]$ docker run -it ubuntu:latest
  root@9f3398da43c2:/# 

Here ``ubuntu`` is referred to as the "image (name)" and ``latest`` as
the "tag".  Other releases of Ubuntu are available under different
tags, such as ``xenial`` or ``devel``.

The above command drops you in a root shell on the container::

  root@9f3398da43c2:/# uname -a
  Linux 9f3398da43c2 4.19.76-linuxkit #1 SMP Thu Oct 17 19:31:58 UTC 2019 x86_64 x86_64 x86_64 GNU/Linux
  root@9f3398da43c2:/# df -h
  Filesystem      Size  Used Avail Use% Mounted on
  overlay         181G  116G   56G  68% /
  tmpfs            64M     0   64M   0% /dev
  tmpfs           2.7G     0  2.7G   0% /sys/fs/cgroup
  shm              64M     0   64M   0% /dev/shm
  /dev/sda1       181G  116G   56G  68% /etc/hosts
  tmpfs           2.7G     0  2.7G   0% /proc/acpi
  tmpfs           2.7G     0  2.7G   0% /sys/firmware

Exiting the shell terminates the container::

  root@9f3398da43c2:/# ^D
  [mkoeppe@sage sage]$

Let us work with a distclean Sage source tree.  If you are using git,
a good way to get one (without losing a precious installation in
``SAGE_LOCAL``) is by creating a new worktree::

  [mkoeppe@sage sage] git worktree add worktree-ubuntu-latest
  [mkoeppe@sage sage] cd worktree-ubuntu-latest
  [mkoeppe@sage worktree-ubuntu-latest] ls
  COPYING.txt ... Makefile ... configure.ac ... src tox.ini

This is not bootstrapped (``configure`` is missing), so let's bootstrap it::

  [mkoeppe@sage worktree-ubuntu-latest] make configure
  ...

We can start a container again with same image, ``ubuntu:latest``, but this time let's
mount the current directory into it::

  [mkoeppe@sage worktree-ubuntu-latest]$ docker run -it --mount type=bind,source=$(pwd),target=/sage ubuntu:latest
  root@39d693b2a75d:/# mount | grep sage
  osxfs on /sage type fuse.osxfs (rw,nosuid,nodev,relatime,user_id=0,group_id=0,allow_other,max_read=1048576)
  root@39d693b2a75d:/# cd sage
  root@39d693b2a75d:/sage# ls
  COPYING.txt ... Makefile ... config configure configure.ac ... src tox.ini
  
Typical Docker images provide minimal installations of packages only::
  
  root@39d693b2a75d:/sage# command -v python
  root@39d693b2a75d:/sage# command -v gcc
  root@39d693b2a75d:/sage#

As you can see above, the image ``ubuntu:latest`` has neither a Python nor
a GCC installed, which are among the build prerequisites of Sage.  We
need to install them using the distribution's package manager first.
 
Sage facilitates testing various distributions on Docker as follows.

Discovering the system's package system
---------------------------------------

::

  root@39d693b2a75d:/sage# build/bin/sage-guess-package-system 
  debian

Let's install gcc, hoping that the Ubuntu package providing it is
simply named ``gcc``.  If we forgot what the package manager on
Debian-derived distributions is called, we can ask Sage for a
reminder::

  root@39d693b2a75d:/sage# build/bin/sage-print-system-package-command debian install gcc
  sudo apt-get install gcc

We are already root, so we can drop the ``sudo``, of course.
And we remember that we need to fetch the current package lists
from the server first::

  root@39d693b2a75d:/sage# apt-get update
  root@39d693b2a75d:/sage# apt-get install gcc

Using Sage's database of distribution prerequisites
---------------------------------------------------

The source code of the Sage distribution contains a database of
package names in various distributions' package managers.  For
example, the file ``build/pkgs/debian.txt`` contains the following

.. code-block:: yaml

  # This file, build/pkgs/debian.txt, contains names of Debian/Ubuntu packages
  # needed for installation of Sage from source.
  #
  # In addition, the files build/pkgs/SPKG/debian.txt contain the names
  # of packages that provide the equivalent of SPKG.
  #
  # Everything on a line after a # character is ignored.
  binutils
  make
  m4
  perl
  # python3-minimal is not enough on debian buster, ubuntu bionic - it does not have urllib
  python3    # system python for bootstrapping the build
  tar
  bc
  gcc
  # On debian buster, need C++ even to survive 'configure'. Otherwise:
  # checking how to run the C++ preprocessor... /lib/cpp
  # configure: error: in `/sage':
  # configure: error: C++ preprocessor "/lib/cpp" fails sanity check
  g++
  # Needed if we download some packages from a https upstream URL
  ca-certificates

From this information, we know that we can use the following command
on our container to install the necessary build prerequisites::

  root@39d693b2a75d:/sage# apt-get install binutils make m4 perl python3 tar bc gcc g++ ca-certificates
  Reading package lists... Done
  Building dependency tree       
  Reading state information... Done
  tar is already the newest version (1.29b-2ubuntu0.1).
  The following additional packages will be installed:
  ...
  Done.

Now we can start the build::

  root@39d693b2a75d:/sage# ./configure 
  checking for a BSD-compatible install... /usr/bin/install -c
  checking for root user... yes
  configure: error: You cannot build Sage as root, switch to an unprivileged user.  (If building in a container, use --enable-build-as-root.)

Let's just follow this helpful hint::

  root@39d693b2a75d:/sage# ./configure --enable-build-as-root
  checking for a BSD-compatible install... /usr/bin/install -c
  ...


Using Sage's database of equivalent distribution packages
---------------------------------------------------------

At the end of the ``./configure`` run, Sage issued a message like the
following::

  configure: notice: the following SPKGs did not find equivalent system packages: arb boost boost_cropped bzip2 ... yasm zeromq zlib
  checking for the package system in use... debian
  configure: hint: installing the following system packages is recommended and may avoid building some of the above SPKGs from source:
  configure:   $ sudo apt-get install libflint-arb-dev ... yasm libzmq3-dev libz-dev
  configure: After installation, re-run configure using:
  configure:   $ ./config.status --recheck && ./config.status

This information comes from Sage's database of equivalent distribution
packages.  For example::

  root@39d693b2a75d:/sage# ls build/pkgs/arb/distros/
  arch.txt	conda.txt	debian.txt	gentoo.txt
  root@39d693b2a75d:/sage# cat build/pkgs/arb/distros/debian.txt 
  libflint-arb-dev

Note that these package equivalencies are based on a current stable or
testing version of the distribution; the packages are not guaranteed
to exist in every release or derivative distribution.

The Sage distribution is intended to build correctly no matter what
superset of the set of packages forming the minimal build
prerequisites is installed on the system.  If it does not, this is a
bug of the Sage distribution and should be reported and fixed on a
ticket.  Crucial part of a bug report is the configuration of the
system, in particular a list of installed packages and their versions.

Let us install a subset of these packages::

  root@39d693b2a75d:/sage# apt-get install libbz2-dev bzip2 yasm libz-dev
  Reading package lists... Done
  ...
  Setting up zlib1g-dev:amd64 (1:1.2.11.dfsg-0ubuntu2) ...
  root@39d693b2a75d:/sage#

  
Committing a container to disk
------------------------------

After terminating the container, we can create a new image
corresponding to its current state::

  root@39d693b2a75d:/sage# ^D
  [mkoeppe@sage worktree-ubuntu-latest]$ docker ps -a | head -n3
  CONTAINER ID        IMAGE                   COMMAND                   CREATED             STATUS
  39d693b2a75d        ubuntu:latest           "/bin/bash"               8 minutes ago       Exited (0) 6 seconds ago
  9f3398da43c2        ubuntu:latest           "/bin/bash"               8 minutes ago       Exited (0) 8 minutes ago
  [mkoeppe@sage worktree-ubuntu-latest]$ docker commit 39d693b2a75d ubuntu-latest-minimal-17
  sha256:4151c5ca4476660f6181cdb13923da8fe44082222b984c377fb4fd6cc05415c1

Here, ``39d693b2a75d`` was the container id (which appeared in the
shell prompts and in the output of ``docker ps``), and
``ubuntu-latest-minimal-17`` is an arbitrary symbolic name for the new
image.  The output of the command is the id of the new image.  We can
use either the symbolic name or the id to refer to the new image.

We can run the image and get a new container with the same state as
the one that we terminated.  Again we want to mount our worktree into
it; otherwise, because we did not make a copy, the new container will
have no access to the worktree::

  [mkoeppe@sage worktree-ubuntu-latest]$ docker run -it \
    --mount type=bind,source=$(pwd),target=/sage ubuntu-latest-minimal-17
  root@73987568712c:/# cd sage
  root@73987568712c:/sage# command -v gcc
  /usr/bin/gcc
  root@73987568712c:/sage# command -v yasm
  /usr/bin/yasm
  root@73987568712c:/sage# ^D
  [mkoeppe@sage worktree-ubuntu-latest]$

The image ``ubuntu-latest-minimal-17`` can be run in as many
containers as we want and can also be shared with other users or
developers so that they can run it in a container on their machine.
(See the Docker documentation on how to do this.)

This facilitates collaboration on fixing portability bugs of the Sage
distribution.  After reproducing a portability bug on a container,
several developers can work on fixing the bug using containers running
on their respective machines.


Generating Dockerfiles
----------------------

Sage also provides a script for generating a ``Dockerfile``, which is
a recipe for automatically building a new image::

  [mkoeppe@sage sage]$ build/bin/write-dockerfile.sh debian "@(standard|optional)" > Dockerfile

(The second argument is a bash extglob pattern for the package type.)
  
.. this interface should be improved obviously. See #29146 - Refactor tox.ini and build/bin/write_dockerfile.sh

The ``Dockerfile`` instructs the command ``docker build`` to build a
new Docker image.  Let us take a quick look at the generated file;
this is slightly simplified::

  [mkoeppe@sage sage]$ cat Dockerfile 
  # Automatically generated by SAGE_ROOT/build/bin/write-dockerfile.sh
  # the :comments: separate the generated file into sections
  # to simplify writing scripts that customize this file
  ...

First, it instructs ``docker build`` to start from an existing base
image...::

  ...
  ARG BASE_IMAGE=ubuntu:latest
  FROM ${BASE_IMAGE}
  ...

Then, to install system packages...::

  ...
  RUN apt-get update && DEBIAN_FRONTEND=noninteractive apt-get install -qqq --no-install-recommends --yes binutils make m4 perl python3 ... yasm libzmq3-dev libz-dev && apt-get clean

Then, to bootstrap and configure...::

  RUN mkdir -p /sage
  WORKDIR /sage
  ADD Makefile VERSION.txt README.md bootstrap configure.ac sage ./
  ADD src/doc/bootstrap src/doc/bootstrap
  ADD m4 ./m4
  ADD build ./build
  ADD src/bin/sage-version.sh src/bin/sage-version.sh
  RUN ./bootstrap
  ADD src/bin src/bin
  ADD src/Makefile.in src/Makefile.in
  ARG EXTRA_CONFIGURE_ARGS=""
  RUN ./configure --enable-build-as-root --enable-option-checking ${EXTRA_CONFIGURE_ARGS} || (cat config.log; exit 1)

Finally, to build and test...::

  ARG NUMPROC=8
  ENV MAKE="make -j${NUMPROC}"
  ARG USE_MAKEFLAGS="-k"
  RUN make ${USE_MAKEFLAGS} base-toolchain
  ARG TARGETS_PRE="sagelib-build-deps"
  RUN make ${USE_MAKEFLAGS} ${TARGETS_PRE}
  ADD src src
  ARG TARGETS="build ptest"
  RUN make ${USE_MAKEFLAGS} ${TARGETS}

You can customize the image build process by passing build arguments to the
command ``docker build``.  For example::

  [mkoeppe@sage sage]$ docker build . -f Dockerfile \
    --build-arg BASE_IMAGE=ubuntu:latest \
    --build-arg NUMPROC=4 \
    --build-arg EXTRA_CONFIGURE_ARGS="--with-python=2"

These arguments (and their default values) are defined using ``ARG``
commands in the ``Dockerfile``.

The above command will build Sage from scratch and will therefore take
quite long.  Let us instead just do a partial build, consisting of one
small package, by setting the arguments ``TARGETS_PRE`` and
``TARGETS``.  We use a silent build (``make V=0``)::

  [mkoeppe@sage sage]$ docker build . -f Dockerfile \
    --build-arg TARGETS_PRE=ratpoints \
    --build-arg TARGETS=ratpoints \
    --build-arg USE_MAKEFLAGS="V=0"
  Sending build context to Docker daemon    285MB
  Step 1/28 : ARG BASE_IMAGE=ubuntu:latest
  ...
  Step 2/28 : FROM ${BASE_IMAGE}
   ---> 549b9b86cb8d
  ...
  Step 25/28 : RUN make SAGE_SPKG="sage-spkg -y -o" ${USE_MAKEFLAGS} ${TARGETS_PRE}
  ...
  make[1]: Entering directory '/sage/build/make'
  sage-logger -p 'sage-spkg -y -o  ratpoints-2.1.3.p5' '/sage/logs/pkgs/ratpoints-2.1.3.p5.log'
  [ratpoints-2.1.3.p5] installing. Log file: /sage/logs/pkgs/ratpoints-2.1.3.p5.log
    [ratpoints-2.1.3.p5] successfully installed.
  make[1]: Leaving directory '/sage/build/make'

  real	0m18.886s
  user	0m1.779s
  sys	0m0.314s
  Sage build/upgrade complete!
  ...
  ---> 2d06689d39fa
  Successfully built 2d06689d39fa

We can now start a container using the image id shown in the last step::

  [mkoeppe@sage sage]$ docker run -it 2d06689d39fa bash
  root@fab59e09a641:/sage# ls -l logs/pkgs/
  total 236
  -rw-r--r-- 1 root root 231169 Mar 26 22:07 config.log
  -rw-r--r-- 1 root root   6025 Mar 26 22:27 ratpoints-2.1.3.p5.log
  root@fab59e09a641:/sage# ls -l local/lib/*rat*
  -rw-r--r-- 1 root root 177256 Mar 26 22:27 local/lib/libratpoints.a
  
You can customize the image build process further by editing the
``Dockerfile``.  For example, by default, the generated ``Dockerfile``
configures, builds, and tests Sage.  By deleting or commenting out the
commands for the latter, you can adjust the Dockerfile to stop after
the ``configure`` phase, for example.

``Dockerfile`` is the default filename for Dockerfiles.  You can
change it to any other name, but it is recommended to use
``Dockerfile`` as a prefix, such as ``Dockerfile-debian-standard``.
It should be placed within the tree rooted at the current directory
(``.``); if you want to put it elsewhere, you need to learn about
details of "Docker build contexts".

Note that in contrast to the workflow described in the above sections,
the ``Dockerfile`` **copies** a snapshot of your Sage worktree into
the build container, using ``ADD`` commands, instead of mounting the
directory into it.  This copying is subject to the exclusions in the
``.gitignore`` file (via a symbolic link from ``.dockerignore``).
Therefore, only the sources are copied, but not your configuration
(such as the file ``config.status``), nor the ``$SAGE_LOCAL`` tree,
nor any other build artefacts.

Because of this, you can build a Docker image using the generated
``Dockerfile`` from your main Sage development tree.  It does not have
to be distclean to start, and the build will not write into it at all.
Hence, you can continue editing and compiling your Sage development
tree even while Docker builds are running.


Debugging a portability bug using Docker
----------------------------------------

Let us do another partial build.  We choose a package that we suspect
might not work on all platforms, ``surf``, which was marked as
"experimental" in 2017::

  [mkoeppe@sage sage]$ docker build . -f Dockerfile \
    --build-arg BASE_IMAGE=ubuntu:latest \
    --build-arg NUMPROC=4 \
    --build-arg TARGETS_PRE=surf \
    --build-arg TARGETS=surf
  Sending build context to Docker daemon    285MB
  Step 1/28 : ARG BASE_IMAGE=ubuntu:latest
  Step 2/28 : FROM ${BASE_IMAGE}
   ---> 549b9b86cb8d
  ...
  Step 24/28 : ARG TARGETS_PRE="sagelib-build-deps"
   ---> Running in 17d0ddb5ad7b
  Removing intermediate container 17d0ddb5ad7b
   ---> 7b51411520c3
  Step 25/28 : RUN make SAGE_SPKG="sage-spkg -y -o" ${USE_MAKEFLAGS} ${TARGETS_PRE}
   ---> Running in 61833bea6a6d
  make -j4 build/make/Makefile --stop
  ...
  [surf-1.0.6-gcc6] Attempting to download package surf-1.0.6-gcc6.tar.gz from mirrors
  ...
  [surf-1.0.6-gcc6] http://mirrors.mit.edu/sage/spkg/upstream/surf/surf-1.0.6-gcc6.tar.gz
  [surf-1.0.6-gcc6] [......................................................................]
  [surf-1.0.6-gcc6] surf-1.0.6-gcc6
  [surf-1.0.6-gcc6] ====================================================
  [surf-1.0.6-gcc6] Setting up build directory for surf-1.0.6-gcc6
  ...
  [surf-1.0.6-gcc6] /usr/bin/ld: cannot find -lfl
  [surf-1.0.6-gcc6] collect2: error: ld returned 1 exit status
  [surf-1.0.6-gcc6] Makefile:504: recipe for target 'surf' failed
  [surf-1.0.6-gcc6] make[3]: *** [surf] Error 1
  ...
  [surf-1.0.6-gcc6] ************************************************************************
  [surf-1.0.6-gcc6] Error installing package surf-1.0.6-gcc6
  [surf-1.0.6-gcc6] ************************************************************************
  ...
  Makefile:2088: recipe for target '/sage/local/var/lib/sage/installed/surf-1.0.6-gcc6' failed
  make[1]: *** [/sage/local/var/lib/sage/installed/surf-1.0.6-gcc6] Error 1
  make[1]: Target 'surf' not remade because of errors.
  make[1]: Leaving directory '/sage/build/make'
  ...
  Error building Sage.

  The following package(s) may have failed to build (not necessarily
  during this run of 'make surf'):

  * package:         surf-1.0.6-gcc6
    last build time: Mar 26 22:07
    log file:        /sage/logs/pkgs/surf-1.0.6-gcc6.log
    build directory: /sage/local/var/tmp/sage/build/surf-1.0.6-gcc6

  ...
  Makefile:31: recipe for target 'surf' failed
  make: *** [surf] Error 1
  The command '/bin/sh -c make SAGE_SPKG="sage-spkg -y -o" ${USE_MAKEFLAGS} ${TARGETS_PRE}' returned a non-zero code: 2

Note that no image id is shown at the end; the build failed, and no
image is created.  However, the container in which the last step of
the build was attempted exists::

  [mkoeppe@sage sage]$ docker ps -a |head -n3
  CONTAINER ID        IMAGE                      COMMAND                   CREATED             STATUS
  61833bea6a6d        7b51411520c3               "/bin/sh -c 'make SA…"    9 minutes ago       Exited (2) 1 minute ago
  73987568712c        ubuntu-latest-minimal-17   "/bin/bash"               24 hours ago        Exited (0) 23 hours ago

We can copy the build directory from the container for inspection::

  [mkoeppe@sage sage]$ docker cp 61833bea6a6d:/sage/local/var/tmp/sage/build ubuntu-build
  [mkoeppe@sage sage]$ ls ubuntu-build/surf*/src
  AUTHORS         TODO            curve           misc
  COPYING         acinclude.m4    debug           missing
  ChangeLog       aclocal.m4      dither          mkinstalldirs
  INSTALL         background.pic  docs            mt
  Makefile        config.guess    draw            src
  Makefile.am     config.log      drawfunc        surf.1
  Makefile.global config.status   examples        surf.xpm
  Makefile.in     config.sub      gtkgui          yaccsrc
  NEWS            configure       image-formats
  README          configure.in    install-sh

Alternatively, we can use ``docker commit`` as explained earlier to
create an image from the container::

  [mkoeppe@sage sage]$ docker commit 61833bea6a6d
  sha256:003fbd511016fe305bd8494bb1747f0fbf4cb2c788b4e755e9099d9f2014a60d
  [mkoeppe@sage sage]$ docker run -it 003fbd511 bash
  root@2d9ac65f4572:/sage# (cd /sage/local/var/tmp/sage/build/surf* && /sage/sage --buildsh)

  Starting subshell with Sage environment variables set.  Don't forget
  to exit when you are done.  Beware:
   * Do not do anything with other copies of Sage on your system.
   * Do not use this for installing Sage packages using "sage -i" or for
     running "make" at Sage's root directory.  These should be done
     outside the Sage shell.

  Bypassing shell configuration files...

  Note: SAGE_ROOT=/sage
  (sage-buildsh) root@2d9ac65f4572:surf-1.0.6-gcc6$ ls /usr/lib/libfl*
  /usr/lib/libflint-2.5.2.so  /usr/lib/libflint-2.5.2.so.13.5.2  /usr/lib/libflint.a  /usr/lib/libflint.so
  (sage-buildsh) root@2d9ac65f4572:surf-1.0.6-gcc6$ apt-get update && apt-get install apt-file    
  (sage-buildsh) root@2d9ac65f4572:surf-1.0.6-gcc6$ apt-file update
  (sage-buildsh) root@2d9ac65f4572:surf-1.0.6-gcc6$ apt-file search "/usr/lib/libfl.a"
  flex-old: /usr/lib/libfl.a
  freebsd-buildutils: /usr/lib/libfl.a
  (sage-buildsh) root@2d9ac65f4572:surf-1.0.6-gcc6$ apt-get install flex-old
  (sage-buildsh) root@2d9ac65f4572:surf-1.0.6-gcc6$ ./spkg-install 
  checking for a BSD-compatible install... /usr/bin/install -c
  checking whether build environment is sane... yes
  ...
    /usr/bin/install -c  surf /sage/local/bin/surf
   /usr/bin/install -c -m 644 ./surf.1 /sage/local/share/man/man1/surf.1
  make[3]: Leaving directory '/sage/local/var/tmp/sage/build/surf-1.0.6-gcc6/src'
  make[2]: Leaving directory '/sage/local/var/tmp/sage/build/surf-1.0.6-gcc6/src'
  make[1]: Leaving directory '/sage/local/var/tmp/sage/build/surf-1.0.6-gcc6/src'
  (sage-buildsh) root@2d9ac65f4572:surf-1.0.6-gcc6$ exit
  root@2d9ac65f4572:/sage# exit
  [mkoeppe@sage sage]$

A standard case of bitrot.
  

Symbolic names of system configurations
---------------------------------------

To facilitate communication about host system configurations and for
automating testing, Sage defines symbolic names composed of several
`Tox "factors" <https://tox.readthedocs.io/en/latest/config.html#complex-factor-conditions>`_ in the file ``$SAGE_ROOT/tox.ini``.

The **system factor** describes a base operating system
image.

- Examples are ``ubuntu-focal``, ``debian-buster``,
  ``archlinux-latest``, ``fedora-30``, ``slackware-14.2``,
  ``centos-7-i386``, and ``ubuntu-bionic-arm64``.

- See ``$SAGE_ROOT/tox.ini`` for a complete list, and to which images
  on Docker hub they correspond.

The **packages factor** describes a list of system packages to be
installed on the system before building Sage:

- ``minimal`` installs the system packages known to Sage to provide
  minimal prerequisites for bootstrapping and building the Sage
  distribution.
  
- ``standard`` additionally installs all known system packages that
  are equivalent to standard packages of the Sage distribution, for
  which the mechanism ``spkg-configure.m4`` is implemented.
  This corresponds to the type pattern ``@(standard)``.

- ``maximal`` does the same for all standard and optional packages.
  This corresponds to the type pattern ``@(standard|optional)``.

The factors are connected by a hyphen to name a system configuration,
such as ``debian-buster-standard`` and ``centos-7-i386-minimal``.


Automatic Docker-based build testing using tox
----------------------------------------------

``tox`` is Python package that is widely used for automating tests of
Python projects.

Install ``tox`` for use with your system Python, for example using::

  [mkoeppe@sage sage]$ pip install --user tox

A tox "environment" is identified by a name composed of several
factors.  In addition to the system factor and the packages factor
that we introduced above, there are the following factors:

The **technology** factor describes how the environment is run:

- ``docker`` builds a Docker image as described above

..
   - ``local`` runs testing on the host OS instead.  This is currently only
     implemented for ``local-homebrew-macos``; see below.

The **configuration** factor (allowed to be empty):

- ``python2`` adds the argument ``--with-python=2`` to the
  ``configure`` run.

The factors are connected by a hyphen to name a tox environment.  To
run an environment::

  [mkoeppe@sage sage]$ tox -e docker-slackware-14.2-minimal
  [mkoeppe@sage sage]$ tox -e docker-ubuntu-bionic-standard-python2
  
Extra arguments to ``docker build`` can be supplied through the
environment variable ``EXTRA_DOCKER_BUILD_ARGS``.  For example, for
a silent build (``make V=0``), use::
  
  [mkoeppe@sage sage]$ EXTRA_DOCKER_BUILD_ARGS="--build-arg USE_MAKEFLAGS=\"V=0\"" \
    tox -e docker-ubuntu-bionic-standard


Automatic parallel tox runs on GitHub Actions
---------------------------------------------

The Sage source tree includes a default configuration for GitHub
Actions that runs tox on a multitude of platforms on every pull
request to a repository for which GitHub Actions are enabled.

This is defined in the file ``$SAGE_ROOT/.github/workflows/tox.yml``.

An additional GitHub Actions workflow for testing on Cygwin, not based
on tox, is defined in the file
``$SAGE_ROOT/.github/workflows/ci-cygwin.yml``.

GitHub Actions runs these build jobs on 2-core machines with 7 GB of
RAM memory and 14 GB of SSD disk space, cf.
`here <https://help.github.com/en/actions/reference/virtual-environments-for-github-hosted-runners#supported-runners-and-hardware-resources>`_,
and has a time limit of 6h per job. This is just barely enough for a
typical ``minimal`` build followed by make ptest to succeed; and
plenty of time for a typical ``standard`` build to succeed.


..
   Workflow for responding to portability bug reports
   ==================================================

   1. Reproduce the platform and the bug.

   2. Fix the bug.

