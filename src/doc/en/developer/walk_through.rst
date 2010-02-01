.. _chapter-walk-through:

=======================================
Walking Through the Development Process
=======================================

This section presents a general process for working on Sage. We will
walk through each step of the process, explaining specific steps and
commands to help you get started on developing Sage. Along the way,
other sources of information will be presented where you could find
more detailed information on some particular issue. If you are a
beginner to Sage development, this introductory guide is here to help
you become familiar with the Sage development process.


.. _section-modify-source:

Modifying Sage source code
--------------------------

The name ``SAGE_ROOT`` refers to the directory where Sage is
installed on your system, e.g. ``/sage/sage-4.1.2/``. Some people
install Sage under their home directory, in which case ``SAGE_ROOT``
would be something like ``/home/username/sage-4.3.1``. Software
packages that are shipped with Sage would then be located under
``/home/username/sage-4.3.1``. In other words, ``SAGE_ROOT`` is the
top-level directory of your Sage installation.

The majority of Sage code for development work is under
``SAGE_ROOT/devel``, where the symbolic link ``SAGE_ROOT/devel/sage``
points to the current branch that you are using. The default branch is
the branch named ``sage-main`` and is also referred to as the main
branch. By default, the symbolic link ``SAGE_ROOT/devel/sage`` points
to ``SAGE_ROOT/devel/sage-main``. Under each branch is a directory named
``sage``, which contains code for interfacing with third-party
packages and code that comprises the Sage library. You can think of
the directory ``SAGE_ROOT/devel/sage-main/sage`` (or
``SAGE_ROOT/devel/sage/sage``) as where the Sage library resides.

Suppose you have modified some part of the Sage library, e.g. add new
code or delete some code. How do you get Sage to know about your
changes? You need to rebuild the Sage library so that it is updated
with your changes. Navigate to ``SAGE_ROOT`` and rebuild the library
as follows::

    ./sage -b main

If ``SAGE_ROOT/devel/sage`` points to the main branch, you could use
the above command or more simply the following command

::

    ./sage -b

The switch ``-b`` is for (re)building the Sage library and the argument
``main`` is the branch you want to rebuild. In this case, you want to
rebuild the Sage library as contained in the main branch, i.e under
``SAGE_ROOT/devel/sage-main``. In fact, the command ``./sage -b``
rebuilds whatever branch that the symbolic link
``SAGE_ROOT/devel/sage`` points to. Any of the above two commands does
not rebuild everything in the Sage library, but only those files in the
library that your changes affect. During the rebuild process, the
affected source files are copied elsewhere and compiled, Cython files
get converted to C code and compiled, and so on. Let this process
manage itself and do not edit the copies of the source files at their
destination.

You may want to create a totally distinct installation of Sage
compiled from source, where you will work so as to not impact a
version of Sage that you actually use. This could be true even if you
make clones (or sandboxes) as described in the next section.


.. _section-create-sandbox:

Creating a sandbox
------------------

With ``SAGE_ROOT`` as the default directory, you can create a new
branch (or clone) named ``test`` with the command::

    ./sage -clone test

Note that a general way to create a clone is the syntax
``./sage -clone <branch-name>``, where you should replace
``<branch-name>`` with the name you want to name your new
clone. Creating a clone may take a couple of minutes to conclude. The
result is a whole new subdirectory ``SAGE_ROOT/devel/sage-test`` where
the name ``sage-`` is automatically prepended to ``test``, and the
symbolic link ``SAGE_ROOT/devel/sage`` will now point to this new
subdirectory. Running Sage, or rebuilding, will use this version of
the source.

Now you can safely experiment with Sage code as much as you would
like. A cloned version of the Sage library can be found under

::

    SAGE_ROOT/devel/sage-test/sage

The Sage library consists of Python files (with the file extension
``.py``) or Cython files (with extension ``.pyx``), organized into
subdirectories that mirror the module hierarchy reflected in the help
files, Python types, etc.  (With the ``sage`` symbolic link, you will
now see why sometimes you see directories that look like
```.../sage/sage/...``.)

Apart from the directory that holds the Sage library, another special
directory is

::

    SAGE_ROOT/devel/sage-test/doc/en/reference

You can think of that directory as containing "table of contents"
files for documentation, each such file having the extension ``.rst``
to indicate that they follow the ReST format for documenting source
code.  Most documentation comes from the docstrings in the source
files. When you build the HTML or PDF versions of the reference
manual, the generated output lands in subdirectories of
``SAGE_ROOT/devel/sage-test/doc/output``.

To make another clone of the main branch, you need to first switch
back to the untouched code in the main branch.  From ``SAGE_ROOT``,
type

::

    ./sage -b main

to switch back to the original version, i.e. the main branch.  Now
make a new clone as shown above.

If a clone, say ``test``, is a mess and has nothing of value to you,
switch back to the main branch, then delete the directory
``SAGE_ROOT/devel/sage-test`` and everything below it. To know which
branch (or clone) you are using, issue the command

::

    ./sage -branch

This will report the current branch that the symbolic link
``SAGE_ROOT/devel/sage`` points to. Such information is also reported
by the command ``./sage`` when the current branch is not the main
branch.

There are many more arguments you can pass to Sage. For a list of
basic arguments, execute

::

    ./sage -help

The command

::

    ./sage -advanced


will report a list of advanced arguments in addition to the list of
basic arguments as output by ``./sage -help``.


.. _section-review-patch-walkthrough:

Reviewing a patch
-----------------

See also the section :ref:`section-review-patches` for further
guidelines on reviewing patches. Before reviewing a patch, you can
choose to download a patch file from the Trac server (the little
download icon next to the file name is the easiest way).  Clicking on
the file name will show you a side-by-side comparison view that is
useful for previewing changes; red shading is deletions, green shading
is additions.

To apply a patch to the code in your sandbox (see
:ref:`section-create-sandbox` for information on creating a sandbox),
follow these steps:

#. Run Sage: from ``SAGE_ROOT``, type ``./sage``.
#. Apply the patch: at the Sage command line, type::

       hg_sage.apply("<full-path-and-filename.patch>")

#. Quit Sage: use the command ``exit``.
#. Rebuild Sage: use the command ``./sage -b`` to rebuild the affected
   files in the Sage library.

In step 2, you are using Sage's simplified interface to the
`Mercurial <http://mercurial.selenic.com>`_
revision control system.  This command will add the patch as a new
"changeset" and "commit" the changes.  At the Sage command line, you
can run ``hg_sage.log()`` to see before/after changes to the Sage
library. In step 4, you should only see a few files copied, modified,
etc.  Unaffected files should not be part of this step.  Look for
compilation errors in this output and modify your changes as
appropriate. Avoid producing patches that result in compilation
errors or errors in building the documentation. (You want a working
Sage installation, right?)

To actually test out a patch, do the following:

#. Experiment with the functionality proposed by the patch. Verify
   results are correct by hand computations, test bad input, outrageous
   situations, etc.
#. Run tests on the affected files. From ``SAGE_ROOT``, issue the
   command ``./sage -t devel/sage-test/path-to-directory-or-file`` to
   run doctests on the affected file(s). Failures should be reported
   on the ticket and are reason to move the ticket to "needs work".
#. If affected files pass tests, then run ``./sage -testall``. This
   will take a while to complete. No, it is not optional.  A reviewer
   or release manager could discover this step was skipped and request
   that you modify your patch to fix any resulting doctest failures.
#. Ensure that the documentation builds. From ``SAGE_ROOT``, run
   ``./sage -docbuild reference html``, which will build the HTML
   version of the documentation.  Check the "look" of affected files
   in the output directory for the documentation (see above).
#. Check for full doctest coverage. From ``SAGE_ROOT``, run
   ``./sage -coverage <file>``  which will provide a complete report.
   Less than 100% coverage is another reason to return a patch to
   "needs work" status.

For more information on doctesting the Sage library, see
:ref:`chapter-doctesting`.


Creating a change
-----------------

To make a change to Sage (fix a bug, add new functionality), proceed
as follows:

#. Make a fresh clone, as discussed in :ref:`section-create-sandbox`.
#. Apply any precursor patches not in your current version, as
   demonstrated in :ref:`section-review-patch-walkthrough`.
#. Edit source files (see :ref:`section-modify-source` for location),
   test building Sage, test functionality, and so on.
#. Once you have something you like, do everything suggested for
   reviewing a patch.  It is a waste of time for a reviewer to start
   on reviewing a patch and find that tests fail, documentation was
   not tested, etc.  It would save any reviewer a lot of time if your
   patches have been fully tested before you submit them for review.
   Everybody makes mistakes, everybody has bugs they did not
   anticipate, and everybody writes code that can be improved---that is
   why there are code reviews.  But do not cut corners.


Submitting a change
-------------------

Here is how to prepare a patch with your changes:

#. Register for a Trac account at the URL
   ``http://trac.sagemath.org/sage_trac/register``. If you have
   problems with the registration process, please refer to the page
   ``http://www.sagemath.org/contact.html`` for the relevant person to
   contact about your registration issues. Most people use some
   variant of their real name, especially if they already have a
   reputation within mathematics.  Edit the main Trac page where there
   is a list of developers and add yourself with a link to your web
   page. Make sure to sort your Trac username alphabetically.
#. If it does not already exist, make a Trac ticket for your changes.
   Provide a one-line summary and then a description of the problem.
   Include a link to a sage-devel discussion if appropriate.  Choose a
   component, if this is a defect or enhancement, set your real name in
   the author field.  It works well if you have your Trac settings such
   that you get an email every time the ticket changes.  Make a note of
   the ticket number.
#. Create a ``.hgrc``  Mercurial configuration file in your home
   directory.  Specify your name and email address here, so it will
   identify you as the author of a patch, in the form `` Bill Smith
   <bsmith@bigu.edu>``. Here is a template for your ``.hgrc`` file:

   ::

       [ui]
       username = Carl Friedrich Gauss <cfgauss@uni-goettingen.de>

       [extensions]
       # Enable the Mercurial queue extension.
       hgext.mq =

   The Mercurial project website ``http://mercurial.selenic.com``
   contains many tutorials on using Mercurial.
#. If necessary, first switch to the branch holding your changes. From
   the Sage command line interface, run ``hg_sage.status()``.  The
   output will be a list of modified files, preceded by a capital ``M``.
   Check that this is what you expect.  For explanation of other
   letters, see the Mercurial documentation on the ``hg status``
   command.
#. From the Sage command line, run ``hg_sage.diff()``. This will show
   you the changes you have made. A plus sign is new code being added,
   a minus sign is code being deleted.  This should look like the
   changes you have made.
#. Now run ``hg_sage.commit()`` from the Sage command line.  This will
   package your changes as a single Mercurial "changeset", allowing
   others (reviewers, release manager) to add your changes to their
   versions of Sage.  An editor window will pop up (set your favorite
   editor in the ``.hgrc`` file mentioned above) where you should
   enter a *one-line* message describing the patch. This message is
   known as the commit message for your patch.  You are encouraged to
   write commit messages of the form
   ``Trac XXXX: <description-goes-here>`` using the Trac ticket number
   and then have a concise description, e.g. "fix echelon form error"
   or "add echelon form over finite fields." Some people also write
   commit messages in the form ``#xxxx: <description-goes-here>``,
   which is also acceptable. A key information to provide in a commit
   message is the ticket number.
#. Run the command ``hg_sage.log()`` from the Sage command line.  The
   first entry should be your changeset.  Note the changeset number,
   which is probably 5 decimal digits.
#. Next, issue the command

   ::

       ``hg_sage.export(<changeset-number>, "/path-to-somewhere/trac_XXXX_short_descriptor.patch")``

   where ``short_descriptor`` is really short, like
   ``echelon_form_fix`` or at most ``finite_field_echelon_form``.
#. You can preview your patch using a "diff viewer". Some people use
   kompare on Linux, others use kdiff3.
#. Upload your patch to the Trac server.
#. Feel free to CC another developer (use their Trac username from the
   list on the main page) if you think they might be able to review your
   change.  If somebody else originated, or commented on the Trac ticket,
   they will be notified of your change if they have set Trac to email
   them of any changes.


Updating a change
-----------------

Your first patch would likely have a review that suggests
changes. Here is one way to update your patch.  (There is probably a
better way, but the following steps should be easy to follow.)

#. Make a new fresh clone.  Read :ref:`section-create-sandbox` to be
   sure you clone the right stuff (i.e. do not clone the branch you
   changed).  We will call this clone ``test2`` here.
#. Apply your patch, but not with ``hg_sage.apply()``.  You want to
   make the changes without doing a commit.  (There is a switch that
   will prevent a commit, but by doing this, you will see how to do
   this at the system level.)  First make
   ``SAGE_ROOT/devel/sage-test2/`` your working directory.  Then at
   the system command line, run::

       patch -p1 /path-to-somewhere/trac_XXXX_short_descriptor.patch

   which will be like you just edited the source files with all the
   changes from your original patch.  Now you can edit to reflect a
   reviewer's suggestions and prepare a new patch.
#. When you upload to Trac, you can replace the file with one of the
   same name. The comments will include an indication of when the
   upload happened, so nobody will be confused about when the
   replacement happened.


Being more efficient
--------------------

Making a new clone for every review and for each revision to a patch
seems rather inefficient.  If you agree, then learn about Mercurial
queues.  They allow you to have two stacks (queues) of
"works-in-progress".  One stack is applied to your Sage installation,
the other is off to the side, i.e. the unapplied stack.  You can
pop/push from one stack to the other.  You can work on a named set of
changes and create a patch file from those changes without making an
irreversible commit. So a reviewer's changes are trivial to make, and
you can easily remove them after a review, or leave them in place to
review something that depends on them.  Also, it is very easy to
update a patch of your own in response to a reviewer's comments
(without making a whole sequence of patches). The following URLs
contain introductory tutorials on using Mercurial queues:

* http://mercurial.selenic.com/wiki/MqExtension
* http://wiki.sagemath.org/MercurialQueues
* https://developer.mozilla.org/en/Mercurial_Queues

The online book
`Mercurial: The Definitive Guide <http://hgbook.red-bean.com>`_
by Bryan O'Sullivan contains numerous examples on using Mercurial. See
especially Chapters 12 and 13 for explanation on how to effectively
use Mercurial queue.
