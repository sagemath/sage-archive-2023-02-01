.. _chapter-trac:

=====================================================
The Sage Trac Server: Submitting Patches and Packages
=====================================================

What should you do with your Mercurial patches for Sage? You should
post them on the Sage trac server.

The Sage trac server, located at
http://trac.sagemath.org/sage_trac/, is where Sage bugs are listed
and patched, new code is posted and reviewed, and ideas for
extending and improving Sage are discussed. Thus if you find a bug
in Sage, or if you have new code to submit, or if you have
corrections for the documentation, you should post on the trac
server.

Items on the server are called "tickets", and anyone may browse the
tickets: just visit http://trac.sagemath.org/sage_trac/report. You
need to open an account, though, if you want to comment on a
ticket, submit a patch, or create a new ticket. See the
`trac server <http://trac.sagemath.org/sage_trac>`_
for more information about obtaining an account. This chapter contains
various guidelines on using the trac server.


Reporting bugs
==============

"The first step is admitting you have a problem."

If you think you have found a bug in Sage, you should first search
through the following Google groups for postings related to your
possible bug:

* ``sage-devel``: http://groups.google.com/group/sage-devel
* ``sage-support``: http://groups.google.com/group/sage-support

Maybe the problem you have encountered has already been discussed. You
should also search the trac server to see if anyone else has opened a
ticket about your bug.

If you do not find anything, and you are not sure that you have found
a bug, ask about it on ``sage-devel``. You might be asked to open a
new ticket on the trac server. As mentioned above, you need an account
to do this. To report a bug, login and click on the "New ticket"
button. Type a meaningful one-liner in the "Short summary" box, with
more information in the larger box below. You should include at least
one explicit, reproducible example illustrating your bug (and/or the
steps required to reproduce the buggy behavior). You should also
include the version of Sage (and any relevant packages) you are using,
and operating system information, being precise as possible (32-bit,
64-bit, ...).

Between the "Summary" and "Full description" boxes, there is a
place to choose the "Type" of the ticket: "Defect", "Enhancement",
or "Task". Use your best judgment here; a bug should probably be
reported as a "Defect".

Choose a priority for your bug, keeping in mind that the "blocker"
label should be used very sparingly. Also pick a component for your
bug; this is sometimes straightforward. If your bug deals with
Sage's calculus implementation, choose "calculus". If it is not
obvious, do your best. Choose a milestone; if you are not sure what
to choose, just choose the numbered version of Sage from the menu
("sage-4.3.3", for example). Type in some helpful keywords. In the
box labeled "Assign to", type "somebody" if you are not sure what
else to do.

Hit the "Preview" button to make sure everything looks okay, and
then hit "Submit ticket".

If you do not have an account on the trac system to report directly,
you are still encouraged to report any possible bug to the
``sage-devel`` mailing list at ``sage-devel@googlegroups.com``.
The list is moderated for new users and requires subscription.
In your bug report to ``sage-devel``, make sure to include the
following information:

- **operating system**: as precise as possible and architecture
  (32-bit, 64-bit, ...)

- affected version: the exact **version number** and the downloaded
  package (source, precompiled, virtual machine image, or an upgrade
  from a previous version (which one?))

- provide a **reproducible example** and/or define the steps to
  reproduce the erroneous behaviour.

Thank you in advance for reporting bugs to improve Sage in the future!


Guidelines for opening tickets
==============================

In addition to bug reports, you should also open a ticket if you
have some new code which extends Sage's capabilities. If you have a
feature request, start a discussion on ``sage-devel`` first,
and then if there seems to be general agreement that you have a
good idea, open a ticket describing the idea.

When you consider opening a new ticket, please bear the following
points in mind.

- Before opening a ticket, make sure that nobody else has opened a
  ticket about the same or closely related issue.

- It is much better to open several specific tickets than one that
  is very broad. Indeed, a single ticket which deals with lots of
  different issues can be quite problematic, and should be avoided.

- Be precise: If foo does not work on OS X but is fine on Linux,
  mention that in the title. Use the keyword option so that
  searches will pick up the issue.

- The problem described in the ticket must be solvable. For
  example, it would be silly to open a ticket whose purpose was
  "Make Sage the best mathematical software in the world". There is
  no metric to measure this properly and it is highly subjective.

- If appropriate, provide URLs to background information or email
  threads relevant to the problem you are reporting.


Patching bugs/working on tickets
================================

If you have code which fixes a bug or deals with some issue in
Sage, here is what to do. First, use Mercurial to create a patch
file. See :ref:`chapter-walk-through` for more information on using
Mercurial to produce/manage patches. If the issue has been reported as
a ticket on the trac server, attach your patch file to that ticket: go
to the ticket, click on the "Attach File" button, and follow the
directions. On the ticket page, you should add a comment explaining
your patch. Some relevant information include:

* The version of Sage you used to create the patch. If the patch is
  based on Sage x.y.z, ensure you include such information.

* If the ticket has more than one patch, explicitly specify which ones
  are to be used. Are all of the patches to be applied? Or only a
  subset of the patches on the ticket?

* If more than one patch is to be applied, state the order in which
  those patches are to be applied.

* Does the ticket depend on another ticket? Sometimes, a ticket
  requires that the patches on another ticket be applied first. Be
  sure to include such information if relevant.

If there is no trac ticket associated to this issue, create one (as
explained in the previous sections) describing the issue and your
solution, and attach your patch.

The following are some other relevant issues:

- Every bug fixed should result in a doctest.

- Cooperative debugging via IRC is faster by at least an order of
  magnitude. If you have not learned how to use IRC, please do so.
  If you have problems using IRC because of firewalls, but you do
  have an account on the machine ``sage.math``, you can use irssi via
  ssh there. If you have a flaky connection, you can use it together
  with the program screen.

- This is not an issue with defects, but there are many enhancements
  possible for Sage and too few developers to implement all the
  good ideas. The trac server is useful for keeping ideas
  in a central place because in the Google groups they tend to get
  lost once they drop off the first page.

- If you are a developer, be nice and try to solve a stale/old
  ticket every once in a while.

- Some people regularly do triage. Triage in this context means
  that we look at new bugs and classify them according to our
  perceived priority. It is very likely that different people will
  see priorities of bugs very differently from us, so please let
  us know if you see a problem with specific tickets.

- **Patches Preferred**: Patches are easier to review, edit and
  can be merged without affecting the history. So we greatly prefer
  patches over Mercurial bundles. If you do have a large number of
  patches, a bundle can still be better than patches. One
  alternative to bundles is to use Mercurial queues to flatten the
  history. That might or might not be desirable. See
  :ref:`chapter-walk-through` for further information on using
  Mercurial queues to produce/manage patches.


.. _section-review-patches:

Reviewing patches
=================

All code that goes into Sage is peer-reviewed, to ensure that the
conventions discussed in this manual are followed, to make sure that
there are sufficient examples and doctests in the documentation, and
to try to make sure that the code does, mathematically, what it is
supposed to.

If someone (other than you) has posted a patch for a ticket on the
trac server, you can review it. Look at the patch (by clicking on
the file name in the list of attachments) to see if it makes sense.
Download it (from the window displaying the patch, see the
"Download" option at the bottom of the page). Apply it (using
``hg_sage.patch('filename')``, for example) to your copy of
Sage, and build Sage with the new code by typing ``sage -b``.

Now ask yourself questions such as the following:

- Does the new source code make sense?

- When you run it in Sage, does it fix the problem reported on the
  ticket?

- Does it introduce any new problems?

- Is it documented sufficiently, including both explanation and
  doctests? This is **very** important: all code in Sage must have
  doctests, so even if the patch is for code which did not have a
  doctest before, the new version must include one. In particular,
  all new code must be **100% doctested**. Use the command
  ``sage -coverage <files>`` to see the coverage percentage of
  ``<files>``.

- In particular, is there a doctest illustrating that the bug has
  been fixed? If a function used to give the wrong answer and this
  patch fixes that, then if possible, it should include a doctest
  illustrating its new success.

- If the patch claims to speed up some computation, does the ticket
  contain code examples to illustrate the claim? The ticket should
  explain the speed efficiency before applying the patch. It should
  also explain the speed efficiency gained after applying the patch.
  In both the "before" and "after" explanation, there should be
  code samples to illustrate the claims. It is not sufficient to
  just mention that the patch results in a speed-up of up to x
  percent or y factor.

- Does the reference manual build without errors? You can test the
  reference manual using the command ``sage -docbuild reference html``
  to build the HTML version. The PDF version of the reference manual
  must also build without errors. Use the command
  ``sage -docbuild reference pdf`` to test it out. The latter command
  requires that you have LaTeX installed on your system.

- Do all doctests pass without errors? You can test the Sage
  library with ``make test`` or ``make ptest`` (edit the number
  of threads in ``$SAGE_ROOT/makefile`` before using ``ptest``). See
  :ref:`chapter-doctesting` for more information on doctesting the
  Sage library.

- Do the code and documentation follow conventions documented in the
  following sections?

  - :ref:`chapter-conventions`
  - :ref:`chapter-python`
  - :ref:`chapter-cython`

If the answers to these and other such reasonable questions are yes,
then you might want to give the patch a positive review. On the main
ticket page, write a comment in the box explaining your review. If you
think there are issues with the patch, explain them in the comment box
as well. Browse the tickets on the trac server to see how things are
done.


Closing tickets
===============

**Do not close tickets.** That is the job of the acting Sage release
manager. If you feel strongly that a ticket should be closed or
deleted, send an email to the current release manager explaining the
situation. You can also comment on the ticket, explaining why it
should be closed. A related issue is re-opening tickets. You should
refrain from re-opening a ticket that is already closed.


Reasons to invalidate tickets
=============================

**One Issue Per Ticket**: A ticket must cover only one issue
and should not be a laundry list of unrelated issues. If a ticket
covers more than one issue, we cannot close it and while some of
the patches have been applied to a given release, the ticket would
remain in limbo.

**No Patch Bombs**: Code that goes into Sage is peer-reviewed. If
you show up with an 80,000 lines of code bundle that completely
rips out a subsystem and replaces it with something else, you can
imagine that the review process will be a little tedious. These
huge patch bombs are problematic for several reasons and we prefer
small, gradual changes that are easy to review and apply. This is
not always possible (e.g. coercion rewrite), but it is still highly
recommended that you avoid this style of development unless there
is no way around it.

**Sage Specific**: Sage's philosophy is that we ship everything
(or close to it) in one source tarball to make debugging possible.
You can imagine the combinatorial explosion we would have to deal
with if you replaced only ten components of Sage with external
packages. Once you start replacing some of the more essential
components of Sage that are commonly packaged (e.g. Pari, GAP,
lisp, gmp), it is no longer a problem that belongs in our tracker.
If your distribution's Pari package is buggy for example, file a
bug report with them. We are usually willing and able to solve
the problem, but there are no guarantees that we will help you
out. Looking at the open number of tickets that are Sage specific,
you hopefully will understand why.

**No Support Discussions**: The trac installation is not meant to
be a system to track down problems when using Sage. Tickets should
be clearly a bug and not "I tried to do X and I couldn't get it to
work. How do I do this?" That is usually not a bug in Sage and it
is likely that ``sage-support`` can answer that question for you. If
it turns out that you did hit a bug, somebody will open a concise
and to-the-point ticket.

**No Closing Or Invalidating**: Unless you have admin powers in
trac (which includes all the people who have ever done releases of
Sage), do not close tickets unless you are explicitly told to do so.
If you think that a ticket is invalid or has been fixed, just comment
on it and the current release manager will take a look and close it if
appropriate.

**Solution Must Be Achievable**: Tickets must be achievable. Many
times, tickets that fall into this category usually ran afoul to
some of the other rules listed above. An example would be to
"Make Sage the best CAS in the world". There is no metric to
measure this properly and it is highly subjective.


Milestones vs. releases
=======================

Milestones are usually goals to be met while working toward a
release. In Sage's trac, we use milestones instead of releases, but
unless somebody volunteers to clean up all the old milestones, we
will stick with the current model. It does not make a whole lot of
difference if we use milestone instead of release.

Finely grained releases are good. Release early and often is the way
to go, especially as more and more patches are coming in.

It is a good idea to make a big release and schedule at least one
more bug fix release after that to sort out the inevitable
"doctest X is broken on distribution Y and compiler Z" problem.
Given the number of compilers and operating systems out there, one
has to be realistic to expect problems. A compile farm would
certainly help to catch issues early.


Assigning tickets
=================

- Each ticket must have a milestone assigned. If you are unsure,
  assign it to the current milestone.

- If a ticket has a patch or spkg that is ready to be reviewed,
  assign it against the current milestone.

- Defect vs. enhancement vs. task: this can be tricky, but a defect
  should be something that leads to an exception or a mathematically
  wrong result.

- If you are unsure to whom to assign the ticket, assign it to
  "somebody" or "tba", which stands for "to be assigned".

- Certain categories have default people who get assigned all
  issues. For example, Jane Smith might be the default person who gets
  assigned all tickets relating to calculus. This means that Jane
  looks after tickets in that category, but not necessarily the person
  who is to fix all open tickets relating to calculus.

- If you have been assigned a ticket, you should either accept it
  or assign it back to "somebody" or "tba". Many people do not accept
  pending tickets at the moment. You have accepted a ticket if your
  name has a star next to it.
