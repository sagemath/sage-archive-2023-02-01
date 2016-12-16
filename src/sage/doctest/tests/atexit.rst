Register an atexit function to remove the file given by the
``DOCTEST_DELETE_FILE`` environment variable::

    sage: import atexit
    sage: fn = os.environ['DOCTEST_DELETE_FILE']
    sage: atexit.register(os.unlink, fn)
    <built-in function unlink>
