r"""
Some tools for developers
"""
#*****************************************************************************
#  Copyright (C) 2011 Nicolas M. Thiery <nthiery at users.sf.net>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#                  http://www.gnu.org/licenses/
#******************************************************************************

def runsnake(command):
    """
    Graphical profiling with ``runsnake``

    INPUT:

    - ``command`` -- the command to be run as a string.

    EXAMPLES::

        sage: runsnake("list(SymmetricGroup(3))")        # optional - requires runsnake

    ``command`` is first preparsed (see :func:`preparse`)::

        sage: runsnake('for x in range(1,4): print x^2') # optional - requires runsnake
        1
        4
        9

    :func:`runsnake` requires the program ``runsnake``. Due to non
    trivial dependencies (python-wxgtk, ...), installing it within the
    Sage distribution is unpractical. Hence, we recommend installing
    it with the system wide Python. On Ubuntu 10.10, this can be done
    with::

        > sudo apt-get install python-profiler python-wxgtk2.8 python-setuptools
        > sudo easy_install RunSnakeRun

    See the ``runsnake`` website for instructions for other platforms.

    :func:`runsnake` further assumes that the system wide Python is
    installed in ``/usr/bin/python``.

    .. seealso::

        - `The runsnake website <http://www.vrplumber.com/programming/runsnakerun/>`_
        - ``%prun``
        - :class:`Profiler`

    """
    import cProfile, os
    from sage.misc.misc import tmp_filename, get_main_globals
    from sage.misc.preparser import preparse
    tmpfile = tmp_filename()
    cProfile.runctx(preparse(command.lstrip().rstrip()), get_main_globals(), locals(), filename=tmpfile)
    os.system("/usr/bin/python -E `which runsnake` %s &"%tmpfile)

def import_statements(*objects, **options):
    """
    Display import statements for the given objects

    INPUT:

    - ``*objects`` -- a sequence of objects
    - ``lazy`` -- a boolean (default: True)

    EXAMPLES::

        sage: import_statements(WeylGroup, lazy_attribute)
        from sage.combinat.root_system.weyl_group import WeylGroup
        from sage.misc.lazy_attribute import lazy_attribute

    If ``lazy`` is True, then :func:`lazy_import` statements are
    displayed instead::

        sage: import_statements(WeylGroup, lazy_attribute, lazy=True)
        from sage.misc.lazy_import import lazy_import
        lazy_import('sage.combinat.root_system.weyl_group', 'WeylGroup')
        lazy_import('sage.misc.lazy_attribute', 'lazy_attribute')

    .. todo::

        This is not correct::

            sage: import_statements(ZZ)
            from sage.categories.euclidean_domains import EuclideanDomains.parent_class

        This should be::

            sage: import_statements(ZZ)      # todo: not implemented
            from sage.rings.integer_ring import ZZ
    """
    lazy = "lazy" in options and options["lazy"]
    if lazy:
        print "from sage.misc.lazy_import import lazy_import"
    for obj in objects:
        if hasattr(obj, '__name__'):
            name = obj.__name__
        else:
            name = obj
        if lazy:
            print "lazy_import('%s', '%s')"%(obj.__module__, name)
        else:
            print "from %s import %s"%(obj.__module__, name)

