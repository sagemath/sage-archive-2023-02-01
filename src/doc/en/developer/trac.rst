.. _chapter-trac:

======================================================
The Sage Trac Server: Submitting Patches and Packages
======================================================

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
ticket, submit a patch, or create a new ticket. See the web page
http://wiki.sagemath.org/TracGuidelines for more information about
obtaining an account, as well as general guidelines for using the
trac server.

Reporting Bugs
==============

"The first step is admitting you have a problem."

If you think you've found a bug in Sage, you should first search
through the Google groups ``sage-devel`` and
``sage-support`` for postings related to your possible bug:
maybe it has already been discussed. You should also search the
trac server to see if anyone else has opened a ticket about your
bug.

If you don't find anything, and you're not sure that you've found a
bug, ask about it on ``sage-devel``.

If you don't find anything, and if you're positive you've found a
bug, open a new ticket on the trac server. As mentioned above, you
need an account to do this. To report a bug, log in and click on
the "New ticket" button. Type a meaningful one-liner in the
"Short summary" box, with more information in the larger box below.
You should include an explicit, reproducible example illustrating
your bug (and/or the steps required to reproduce the buggy
behavior). You should also include the version of Sage (and any
relevant packages) you are using, and operating system information,
being precise as possible (32-bit, 64-bit, ...).

Between the "Summary" and "Full description" boxes, there is a
place to choose the "Type" of the ticket: "Defect", "Enhancement",
or "Task". Use your best judgment here; a bug should probably be
reported as a "Defect".

Choose a priority for your bug, keeping in mind that the "blocker"
label should be used very sparingly. Also pick a component for your
bug; this is sometimes straightforward - if your bug deals with
Sage's calculus implementation, choose "calculus". If it is not
obvious, do your best. Choose a milestone; if you're not sure what
to choose, just choose the numbered version of sage from the menu
("sage-3.1.4", for example). Type in some helpful keywords. In the
box labeled "Assign to", type "somebody" if you're not sure what
else to do.

Hit the "Preview" button to make sure everything looks okay, and
then hit "Submit ticket".

Guidelines for Opening Tickets
==============================

In addition to bug reports, you should also open a ticket if you
have some new code which extends Sage's capabilities. If you have a
feature request, start a discussion on ``sage-devel`` first,
and then if there seems to be general agreement that you have a
good idea, open a ticket describing the idea.

Other comments:


-  Before opening a ticket, make sure that nobody else has opened a
   ticket about the same or closely related issue.

-  It is much better to open several specific tickets than one that
   is very broad. Indeed, a single ticket which deals with lots of
   different issues can be quite problematic, and should be avoided.

-  Be precise: If foo doesn't work on OS X but is fine on Linux,
   mention that in the title; also use the keyword option so that
   searches will pick up the issue.

-  The problem described in the ticket must be solvable. For
   example, it would be silly to open a ticket whose purpose was
   "Make Sage the best mathematical software in the world". There is
   no metric to measure this properly and it is highly subjective.


Patching Bugs/Working on Tickets
================================

If you have code which fixes a bug or deals with some issue in
Sage, here's what to do: first, use Mercurial to create a patch
file. If the issue has been reported as a ticket on the trac
server, attach your patch file to that ticket: go to the ticket,
click on the "Attach File" button, and follow the directions. On
the ticket page, you should add a comment explaining your patch,
and you should also change the summary for the ticket from
``description of bug`` to
``[with patch, needs review] description of bug``.

If there is no trac ticket associated to this issue, create one (as
explained in the previous sections) describing the issue and your
solution, attach your patch, and give it a summary of the form
``[with patch, needs review] description of bug goes here``.

Reviewing Patches
=================

All code that goes into Sage is peer-reviewed, to ensure that the
conventions discussed in this manual are followed, to make sure
that there are sufficient examples and doctests in the
documentation, and to try to make sure that the code does,
mathematically, what it is supposed to.

If someone (other than you) has posted a patch for a ticket on the
trac server, you can review it. Look at the patch (by clicking on
the file name in the list of attachments) to see if it makes sense.
Download it (from the window displaying the patch, see the
"Download" option at the bottom of the page). Apply it (using
``hg_sage.patch('filename')``, for example) to your copy of
Sage, and build Sage with the new code by typing
``sage -b``.

Now ask yourself questions like these:


-  Does the new source code make sense?

-  When you run it in Sage, does it fix the problem reported on the
   ticket?

-  Does it fail to introduce any new problems?

-  Is it documented sufficiently, including both explanations and
   doctests? (This is **very** important: all code in Sage must have
   doctests, so even if the patch is for code which didn't have a
   doctest before, the new version must include one.)

-  In particular, is there a doctest illustrating that the bug has
   been fixed? If a function used to give the wrong answer and this
   patch fixes that, then if possible, it should include a doctest
   illustrating its new success.


If the answers to these and other such reasonable questions are
yes, then you might want to give the patch a positive review. On
the main ticket page, write a comment in the box and change the
summary from ``[with patch, needs review] description of bug`` to
``[with patch, positive review] description of bug``. If you feel
there are issues with the patch, explain them in the comment box,
and change the summary to
{[with patch, negative review]
description of bug}, or
{[with patch, needs work]
description of bug}, or
{[with patch, positive review pending fixes]
description of bug},
or something similar. Browse the tickets on the trac server to see
how things are done.

By the way, if you review a patch which deals with the Sage
manuals, say, instead of the source code, then you need to use
``hg_doc.patch('filename')`` instead of
``hg_sage.patch('filename')`` to apply it, and you need to
follow the directions in :ref:`chapter-sage_manuals` to build the new
documentation.

Closing Tickets
===============

Don't close tickets. That is the job of the acting Sage release
manager. If you feel strongly that a ticket should be closed or
deleted, send email to the current release manager explaining the
situation.