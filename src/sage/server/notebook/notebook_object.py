import twist

class NotebookObject:
    """
    Start the notebook.

    Type notebook.notebook? for more help.
    """
    def __call__(self, *args, **kwds):
        return self.notebook(*args, **kwds)

    notebook = twist.notebook_twisted
    setup = twist.notebook_setup

notebook = NotebookObject()
