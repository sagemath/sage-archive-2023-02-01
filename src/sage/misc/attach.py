class Attach:
    r"""
    Attach a file to a running instance of SAGE.

    attach is *not* a function and is not part of the Python language.
    You attach a file, e.g., foo.sage or foo.py or foo.spyx, to a
    running SAGE session by typing

        sage.: attach foo.sage   # or foo.py   or foo.spyx

    The contents of the file are then loaded, which means they are
    read into the running SAGE session.  For example, if foo.sage
    contains 'x=5', after attaching foo.sage the variable x will be
    set to 5.  Moreover, any time you change foo.sage, the attached
    file will be re-read automatically (with no intervention on your
    part).

    INPUT: attach file1 file2 ...
        -- space-separated list of .py, .spyx, and .sage files.

    EFFECT:
        -- Each file is read in and added to an internal list of watched files.
           The meaning of reading a file in depends on the file type:
               .py   -- read in with no preparsing (so, e.g., 2^3 is 2 bit-or 3),
               .sage -- preparsed then the result is read in
               .spyx -- *not* preparsed.  Compiled to a module m then "from m import *"
                        is executed.

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
