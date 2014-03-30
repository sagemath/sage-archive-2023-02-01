.. _chapter-walkthrough:

========================
Sage Development Process
========================

This section is a concise overview of the Sage development process. In
it, we will see how to make changes to the Sage source code and
communicate these changes back to the Sage project. If you are a
beginner to Sage development, this introductory guide is here to help
you become familiar with the Sage development process.

Most of the commands in the following section will not work unless you
have an account on the bug tracker. If you want to contribute to Sage,
it is a good idea to get an account now (see
:ref:`section-trac-account`). Alternatively, see the
:ref:`section-walkthrough-readonly` section.


.. _section-walkthrough-setup:

Development Setup
=================

.. _section-walkthrough-setup-git:

Installing Git
--------------

First, open a shell and run::
    
    [user@localhost]$ git
    usage: git [--version] [--help] [-C <path>] [-c name=value]
               [--exec-path[=<path>]] [--html-path] [--man-path] [--info-path]
               [-p|--paginate|--no-pager] [--no-replace-objects] [--bare]
               [--git-dir=<path>] [--work-tree=<path>] [--namespace=<name>]
               <command> [<args>]
    
    The most commonly used git commands are:
       add        Add file contents to the index
       bisect     Find by binary search the change that introduced a bug
       branch     List, create, or delete branches
       checkout   Checkout a branch or paths to the working tree
       clone      Clone a repository into a new directory
       commit     Record changes to the repository
       diff       Show changes between commits, commit and working tree, etc
       fetch      Download objects and refs from another repository
       grep       Print lines matching a pattern
       init       Create an empty Git repository or reinitialize an existing one
       log        Show commit logs
       merge      Join two or more development histories together
       mv         Move or rename a file, a directory, or a symlink
       pull       Fetch from and integrate with another repository or a local branch
       push       Update remote refs along with associated objects
       rebase     Forward-port local commits to the updated upstream head
       reset      Reset current HEAD to the specified state
       rm         Remove files from the working tree and from the index
       show       Show various types of objects
       status     Show the working tree status
       tag        Create, list, delete or verify a tag object signed with GPG
    
    'git help -a' and 'git help -g' lists available subcommands and some
    concept guides. See 'git help <command>' or 'git help <concept>'
    to read about a specific subcommand or concept.

This is your basic entry point into git, and it explains how to get
help for the subcommands. If you get a "command not found" error, then
you don't have git installed. Now is the time to install it, see
:ref:`chapter-git-setup` for instructions.

While you are at it, tell git how you want to be known::

    [user@localhost]$ git config --global user.name "Your Name"
    [user@localhost]$ git config --global user.email you@yourdomain.example.com


.. _section-walkthrough-setup-git-trac:

Installing the Git-Trac Command
-------------------------------

Git is a separate project from trac, and the two do not know how to
talk to each other. To simplify the development, we have a special
``git trac`` subcommand of the git suite. Note that this really is
only to simplify interaction with our trac issue management, you can
perform every development task with just git and a web browser to
interact with trac. See :ref:`chapter-manual-git` instead if you
prefer to do everything by hand::

    [user@localhost]$ git clone https://github.com/sagemath/git-trac-command.git
    Cloning into 'git-trac-command'...
    [...]
    Checking connectivity... done.
    [user@localhost]$ source git-trac-command/enable.sh
    Prepending the git-trac command to your search PATH

This creates a directory ``git-trac-command``. Sourcing the
``enable.sh`` script in there is just a quick and dirty way to enable
it temporarily. You probably want a more permanent installation on
your system later, see `the README
<https://github.com/sagemath/git-trac-command>`_ for details.


.. _section-walkthrough-setup-source:

The Sage Source Code
--------------------

Obviously you need the Sage source code to develop, so we download it
from github (which is a public read-only mirror of our internal git
repository)::

    [user@localhost]$ git clone git://github.com/sagemath/sage.git
    Cloning into 'sage'...
    [...]
    Checking connectivity... done.
    
This creates a directory named ``sage`` containing the Sage
sources. Now go into that directory and tell ``git trac`` about your
trac account::

    $ git trac config --user USERNAME --pass PASSWORD
    Trac xmlrpc URL:
        http://trac.sagemath.org/xmlrpc (anonymous)
        http://trac.sagemath.org/login/xmlrpc (authenticated)
        realm sage.math.washington.edu
    Username: USERNAME
    Password: PASSWORD
    Retrieving SSH keys...
        1024 ab:1b:7c:c9:9b:48:fe:dd:59:56:1e:9d:a4:a6:51:9d  My SSH Key
    
where you have to replace USERNAME with your trac user name and
PASSWORD with your trac password, of course. The password is stored in
``.git/config``, so make sure that it is not readable by other users
on your system.

If there is no SSH key listed then then you haven't uploaded your SSH
public key to the trac server. You should do that now following the
instructions to :ref:`section-trac-ssh-key`. Whereas trac uses a
password, you must have a SSH key set up to upload anything to our git
repository.


.. _section-walkthrough-readonly:

Readonly Access
---------------

Note that the ``git trac config`` command will automatically add a
``trac`` remote git repository to your list of remotes if
necessary. Hence, following the above instructions you have two remote
repositories set up::

    [user@localhost]$ git remote -v
    origin      git://github.com/sagemath/sage.git (fetch)
    origin      git://github.com/sagemath/sage.git (push)
    trac        git@trac.sagemath.org:sage.git (fetch)
    trac        git@trac.sagemath.org:sage.git (push)

If you **do not have a trac account** you can setup ``trac`` as readonly::

    [user@localhost]$ git trac config --readonly
    [user@localhost]$ git remote -v
    origin      git://github.com/sagemath/sage.git (fetch)
    origin      git://github.com/sagemath/sage.git (push)
    trac        git://trac.sagemath.org/sage.git (fetch)
    trac        git://trac.sagemath.org/sage.git (push)

If you do not want to use the ``git trac`` subcommand then you can set
up the remote by hand as described in the section on
:ref:`section-git-trac`.
  

.. _section-walkthrough-add:

Contributing to the Sage Source Code
====================================

.. _section-walkthrough-add-create:

Create a Ticket
---------------

Suppose you have written an algorithm for calculating the last twin prime, and
want to add it to Sage. You would first open a ticket for that::

    [user@localhost]$ git trac create 'Last Twin Prime'
    Remote branch: u/user/last_twin_prime
    Newly-created ticket number: 12345
    Ticket URL: http://trac.sagemath.org/12345
    Local branch: t/12345/last_twin_prime

This will create a new trac ticket titled "Last Twin Prime" with a
*remote branch* ``u/user/last_twin_prime`` attached to it. The remote
branch name is automatically derived from the ticket title; If you
don't like this then you can use the ``-b`` switch to specify it
explicitly. See ``git trac create -h`` for details. This new branch is
automatically checked out for you with the *local branch* name
``t/12345/last_twin_prime``.

.. note::

    Only some trac fields are filled in automatically. See
    :ref:`section-trac-fields` for what trac fieds are available and
    how we use them.

Alternatively, you can use the `web interface to the Sage trac
development server <http://trac.sagemath.org>`_ to open a new ticket,
just log in and click on "Create Ticket". Or maybe somebody else
already opened a ticket. Then, to get a suitable local branch to make
your edits, you would just run::

    [user@localhost]$ git trac checkout 12345
    Loading ticket #12345...
    Checking out Trac #13744 remote branch u/user/last_twin_prime -> local branch t/12345/last_twin_prime...

The ``git trac checkout`` command downloads an existing branch (as
specified in the "Branch:" field on the trac ticket) or creates a new
one if there is none yet. Just like the create command, you can
specify the remote branch name explicitly using the ``-b`` switch if
you want.



.. _section-walkthrough-branch-names:

Note on Branch Names
--------------------

Trac tickets that are finished or in the process of being worked on
can have a git branch attached to them. This is the "Branch:" field in
the ticket description. The branch name is generally of the form
``u/user/description``, where ``user`` is the name of the user who
made the branch and ``description`` is some free-form short
description (and can include further slashes, but not whitespace). Our
git server implements the following access restrictions for **remote
branch names**:

* Only the developer with the ``user`` trac account can create
  branches starting with ``u/user/``.

* Everybody can write to branches named ``public/description``.

Depending on your style of collaboration, you can use one or the
other. The ``git trac`` subcommands defaults to the former.

As a convention, the ``git trac`` subcommand uses **local branch
names** of the form ``t/12345/description``, where the number is the
trac ticket number. The script uses this number to figure out the
ticket from the local branch name. You can rename the local branches
if you want, but if they don't contain the ticket number then you will
have to specify the ticket number manually when you are uploading your
changes.


.. _section-walkthrough-add-edit:

Editing the Source Code
-----------------------

A branch is a copy (except that it doesn't take up twice the space) of
the Sage source code where you can store your modifications to the
Sage source code and which you can upload to trac tickets. If you used
the ``git trac`` script to :ref:`section-walkthrough-add-create`, then
you have a local branch already. Otherwise, see the
:ref:`section-walkthrough-add-local` section. You can list all
branches using::

    [user@localhost]$ git branch
      master
    * t/12345/last_twin_prime

The star indicates the currently active branch. To switch between your
local branches, use ``git checkout``::

    [user@localhost]$ git checkout master
    Switched to branch 'master'
    Your branch is up-to-date with 'github/master'.
    [user@localhost]$ git branch
    * master
      t/12345/last_twin_prime
    [user@localhost]$ git checkout t/12345/last_twin_prime
    Switched to branch 't/12345/last_twin_prime'
    Your branch is up-to-date with 'trac/u/user/last_twin_prime'.

Note that, unless you explicitly upload ("push") a branch to remote
git repository, the local branch will only be on your computer and not
visible to anyone else.

At this point you can start editing the source code. The subsequent
chapters of this developer guide explain how your code should look
like to fit into Sage, and how we ensure high code quality
throughout. Whenever you have reached one of your goals, you should
make a *commit*. This takes a snapshot of the whole Sage source code
that you have been working on and records the changes into your local
branch::

    [user@localhost]$ git add src/sage/primes/last_pair.py
    [user@localhost]$ git status
    # On branch t/12345/last_twin_prime
    # Changes to be committed:
    #   (use "git reset HEAD <file>..." to unstage)
    #
    #   new file:   src/sage/primes/last_pair.py
    #
    [user@localhost]$ git commit -m 'found the last prime pair'

Note that you always have to explictly add changed files to the
staging area in order to be able to commit them. You can read more
about that in the :ref:`section-git-commit` section.

If you are working on a larger project, it can be useful to break up
your work into multiple commits: Each commit is saved, enabling you to
retrieve older versions of files from the repository. So, even if you
accidentally delete something, you can get it back later. Also, if you
find a mistake in one of your earlier commits, then you just correct
it in the Sage source code and then add another commit on top.


.. _section-walkthrough-add-push:

Uploading Changes to Trac
-------------------------

At some point, you may wish to share your changes with the rest of us:
maybe it is ready for review, or maybe you are collaborating with
someone and want to share your changes "up until now". This is simply
done by::

    [user@localhost]$ git trac push
    Pushing to Trac #12345...
    Guessed remote branch: u/user/last_twin_prime

    To git@trac.sagemath.org:sage.git
     * [new branch]      HEAD -> u/user/last_twin_prime

    Changing the trac "Branch:" field...

This uploads your changes to a remote branch on the `Sage git server
<http://git.sagemath.org/sage.git>`_. The ``git trac`` command uses
the following logic to find out the remote branch name:

* By default, the remote branch name will be whatever is already on
  the trac ticket.

* If there is no remote branch yet, the branch will be called
  ``u/user/description`` (``u/user/last_twin_prime`` in the example).
  
* You can use the ``--branch`` option to specify the remote branch
  name explicitly, but it needs to follow the naming convention from
  :ref:`section-walkthrough-branch-names` for you to have write
  permission.

It is common to go through some iterations of commits before you
upload, and you will probably also have pushed your changes a few
times before your changes are ready for review.

If you are happy with the changes you uploaded, you want somebody else
to review them, so they can be included into the next version of
Sage. If your ticket is ready for review, you should set it to
``needs_review`` on the trac server. Also, add yourself as an author
for that ticket by inserting the following as the first line::

    Authors: Your Real Name


.. _section-walkthrough-add-pull:

Downloading Changes from Trac
-----------------------------

If somebody else worked on a ticket, or if you just switched
computers, you'll want to get the latest version of the branch from a
ticket into your local branch. This is done with::

    [user@localhost]$ git trac pull

Technically, this does a "merge" (just like the standard ``git pull``)
command. See :ref:`section-git-merge` for more background information.


.. _section-walkthrough-add-local:

Starting Without a Ticket
-------------------------

You might not want to create a trac ticket for your changes. For
example, if you are only working on your own code or if you are making
experimental changes that you are likely to throw away if they do not
work out. In that case, you can also start a branch that only lives in
your local repository. To do this, you use checkout but specify a
branch name instead of the ticket number. For example, to create a new
branch ``my_branch``, you would run::

    [user@localhost]$ git branch my_branch master
    [user@localhost]$ git checkout my_branch

which is equivalent to the abbreviated version::

    [user@localhost]$ git checkout -b my_branch master

The newly created branch starts at ``master`` as specified, but you
can use any other starting point. 

You can upload your local branch later to an existing ticket. This
works exactly like in the case where you started with a ticket, except
that you have to specify the ticket number. That is::

    [user@localhost]$ git trac push TICKETNUM
    
where you have to replace ``TICKETNUM`` with the number of the trac
ticket.


.. _section-walkthrough-merge:

Merging
=======

As soon as you are working on a bigger project that spans multiple
tickets you will want to base your work on branches that have not been
merged into Sage yet. This is natural in collaborative development,
and in fact you are very much encouraged to split your work into
logically different parts. Ideally, each part that is useful on its
own and and can be reviewed independently should be a different
ticket, instead of a huge patch bomb.

For this purpose, you can incorporate branches from other tickets (or
just other local branches) into your current branch. This is called
merging, and all it does is include commits from other branches into
your current branch. In particular, this is done when a new Sage
release is made: the finished tickets are merged with the Sage master
and the result is the next Sage version. Git is smart enough to not
merge commits twice. In particular, it is possible to merge two
branches, one of which had already merged the other branch. The syntax
for merging is easy::

    [user@localhost]$ git merge other_branch

This creates a new "merge" commit, joining your current branch and
``other_branch``.

.. warning::

    You should avoid merging branches both ways. Once A merged B and B
    merged A, there is no way to distinguish commits that were
    originally made in A or B. Effectively, merging both ways combines
    the branches and makes individual review impossible.

    In practice, you should only merge when one of the following holds:

    * Either two tickets conflict, then you have to merge one into the
      other in order to resolve the merge conflict.

    * Or you definitely need a feature that has been developed as part
      of another branch.

A special case of merging is merging in the ``master`` branch. This
brings your local branch up to date with the newest Sage version. The
above warning against unnecessary merges still applies, though. Try to
do all of your development with the Sage version that you originally
started with. The only reason for merging in the master branch is if
you need a new feature or if your branch conflicts.


.. _section-walkthrough-review:

Reviewing
=========

This section gives an example how to review using the ``sage``
command.  For a detailed discussion of Sage's review process, see
:ref:`Reviewing Patches <section-review-patches>`. If you go to the
`web interface to the Sage trac development server
<http://trac.sagemath.org>`_ then you can click on the "Branch:" field
and see the code that is added by combining all commits of the
ticket. This is what needs to be reviewed.

The ``git trac`` command gives you two commands that might be handy
(replace ``12345`` with the actual ticket number) if you do not want
to use the web interface:

* ``git trac get 12345`` displays the trac ticket directly in your
  terminal.

* ``git trac review 12345`` downloads the branch from the ticket and
  shows you what is being added, analogous to clicking on the
  "Branch:" field.



.. _section-walkthrough-collaborate:

Collaboration
=============

It is very easy to collaborate by just going through the above steps any number of times::

    # Alice
    <EDIT EDIT>
    git add .
    git commit
    git trac push

    # Bob
    git trac pull
    <EDIT EDIT>
    git add .
    git commit 
    git trac push

    # Alice
    git trac pull
    <EDIT EDIT>
    git add .
    git commit 
    git trac push
    (etc)

The obvious problem is when you both work on the same ticket simultaneously::

    # Alice
    <EDIT EDIT>
    git add .
    git commit
    git trac push

    # Bob
    <EDIT EDIT>
    git add .
    git commit
    git trac push

Bob gets an error message since the remote changed outside of his
control. The impolite solution is to overwrite Alice's changes::
  
    # Bob should not do this
    git trac push --force

but this is probably not the correct solution. Instead, Bob should
download Alice's changes first::

    # Bob should do this instead
    git trac pull

This will try to merge the changes that Alice made into the ones that
Bob made. Then Bob should check whether all seems okay, and if so,
upload the changes::

    # Bob
    git trac push   # works now

It is possible that the changes cannot be automatically merged. In
that case, Bob will have to do some manual fixup after downloading and
before uploading::
  
    # Bob
    <EDIT EDIT FOR FIXUP>
    git add .
    git commit -m "Resolved the merge conflict"
    git trac push


