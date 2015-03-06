"""
Old Displayhook

Just for compatibility with the notebook, you should not use this any
more. Look into ``sage/repl/display`` instead.
"""

from sage.repl.display.python_hook import DoctestDisplayHook
from sage.repl.display.formatter import SageNBTextFormatter


class DisplayHook(DoctestDisplayHook):

    def __init__(self):
        """
        Python constructor

        EXAMPLES::

            sage: from sage.repl.display.python_hook import DoctestDisplayHook
            sage: d = DoctestDisplayHook()
            sage: d(set([1, 2, 3]))       # Sage commandline output
            {1, 2, 3}
            sage: print(set([1, 2, 3]))   # Plain Python output
            set([1, 2, 3])
        """
        # do not call super, we set our own formatter
        self.formatter = SageNBTextFormatter()
