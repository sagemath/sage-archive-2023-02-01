"""
Attach a file to a running instance of Sage.
"""

class Attach:
    """
    Attach a file to a running instance of Sage.
    """
    def __repr__(self):
        return self.__doc__

    def __call__(self, *args, **kwds):
        raise RuntimeError, "Use 'attach filename' instead, where filename is a .py, .sage, or .spyx file."

attach = Attach()
