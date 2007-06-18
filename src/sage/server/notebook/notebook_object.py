#############################################################################
#       Copyright (C) 2007 William Stein <wstein@gmail.com>
#  Distributed under the terms of the GNU General Public License (GPL)
#  The full text of the GPL is available at:
#                  http://www.gnu.org/licenses/
#############################################################################


"""
Notebook control object

This is used for configuring and starting the SAGE notebook server.
"""

import run_notebook

class NotebookObject:
    """
    Start the notebook.

    Type notebook.notebook? for more help.
    """
    def __call__(self, *args, **kwds):
        return self.notebook(*args, **kwds)

    notebook = run_notebook.notebook_twisted
    setup    = run_notebook.notebook_setup

notebook = NotebookObject()
