"""
Attach a file to a running instance of Sage.
"""

class Attach:
    r"""
    Attach a file to a running instance of Sage.

    \note{\code{attach} is \strong{not} a function and is not part of the Python language.}

    \code{load} is exactly the same as attach, but doesn't
    automatically reload a file when it changes.

    You attach a file, e.g., \file{foo.sage} or \file{foo.py} or
    \file{foo.spyx}, to a running Sage session by typing

        sage: attach foo.sage   # or foo.py   or foo.spyx  (not tested)

    The contents of the file are then loaded, which means they are
    read into the running Sage session.  For example, if
    \file{foo.sage} contains \samp{x=5}, after attaching
    \file{foo.sage} the variable \var{x} will be set to 5.  Moreover,
    any time you change \file{foo.sage}, the attached file will be
    re-read automatically (with no intervention on your part).

    USAGE: \samp{attach file1 file2 ...}
        -- space-separated list of .py, .spyx, and .sage files.

    EFFECT:
           Each file is read in and added to an internal list of watched files.
           The meaning of reading a file in depends on the file type:
           \begin{itemize}
           \item[.py] -- read in with no preparsing (so, e.g., \code{2\^3} is 2 bit-xor 3),
           \item[.sage] -- preparsed then the result is read in
           \item[.spyx] -- \strong{not} preparsed.  Compiled to a module \module{m} then \code{from m import *}
                        is executed.
           \end{itemize}

    Type \code{attached_files()} for a list of all currently attached files.
    You can remove files from this list to stop them from being watched.

    NOTE: attach is exactly the same as load, except it keeps track of the
    loaded file and automatically reloads it when it changes.
    """
    def __repr__(self):
        return self.__doc__

    def __call__(self, *args, **kwds):
        raise RuntimeError, "Use 'attach filename' instead, where filename is a .py, .sage, or .spyx file."

attach = Attach()
