def autogen_all(force=False):
    """
    Regenerate the automatically generated files of the Sage library.

    INPUT:

    - ``force`` -- whether we force rebuilding the files (default is ``False``)
    """
    import pari
    pari.rebuild(force=force)
