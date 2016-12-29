import os


def autogen_all():
    """
    Regenerate the automatically generated files of the Sage library.
    """
    from sage_setup.autogen import pari
    pari.rebuild()

    from sage_setup.autogen import interpreters
    interpreters.rebuild(os.path.join("sage", "ext", "interpreters"))

    # Return list of sub-packages that should be appended to the list of
    # packages built/installed by setup.py
    #
    # In the case of Pari it just adds new files to an existing package (rather
    # than autogenerating the entire sub-package) so it's omitted here.
    return ['sage.ext.interpreters']
