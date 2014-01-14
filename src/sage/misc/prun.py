"""
Prun -- running notebook cell contents in a profiler

AUTHORS:

- Timo Kluck (2012-09-05): initial version
"""

#*****************************************************************************
#       Copyright (C) 2012 William Stein <wstein@gmail.com>
#                          Timo Kluck <tkluck@infty.nl>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#
#                  http://www.gnu.org/licenses/
#*****************************************************************************

from preparser import preparse_file
class Prun(object):
    r"""
    Run cell contents in a profiler

    When the first line in a notebook cell contains (only) ``%prun``, the
    cell contents is executed by a profiler, which means that once execution
    has finished, you will get a listing of all functions called and the time
    spent in each.

    You can influence the appearance of the listing by passing arguments to
    the directive. Use

        ``%prun(sort_args, print_args)``

    where ``sort_args`` is a list of arguments to be passed to ``pstats.Stats::sort_stats``
    and where ``print_args`` is a list of arguments to be passed to ``pstats.Stats::print_stats``.

    For example,

        ``%prun(['time'],[0.1])``

    will sort the functions by time, and will only print the first 10%. See the Python
    documentation for ``pstats.Stats`` for details.

    AUTHOR:

    - Timo Kluck (September 2012)
    """
    def eval(self, cmd, sage_globals, locals=None):
        r"""
        EXAMPLES::

            sage: from sage.misc.prun import prun
            sage: prun.eval("len([1,2,3])", globals(), locals())
            20 function calls in 0... seconds
            <BLANKLINE>
            Ordered by: standard name
            <BLANKLINE>
            ...
        """
        import cProfile
        cmd = preparse_file(cmd, sage_globals)
        return cProfile.runctx(cmd, sage_globals, locals)

    def __call__(self, sort_args, print_args):
        r"""
        EXAMPLES::

            sage: from sage.misc.prun import prun
            sage: prun(['time'],[0.1]).eval("len([1,2,3])", globals(), locals())
            20 function calls in 0... seconds
            <BLANKLINE>
            Ordered by: internal time
            <BLANKLINE>
            ...
        """
        class PrunWithOptions(object):
            def eval(self, cmd, sage_globals, locals=None):
                r"""
                EXAMPLES::

                    sage: from sage.misc.prun import prun
                    sage: prun(['time'],[0.1]).eval("len([1,2,3])", globals(), locals())
                    20 function calls in 0... seconds
                    <BLANKLINE>
                    Ordered by: internal time
                    <BLANKLINE>
                    ...
                """
                import cProfile, pstats
                cmd = preparse_file(cmd, sage_globals)
                stats =  pstats.Stats(cProfile.Profile().runctx(cmd, sage_globals, locals))
                stats.sort_stats(*sort_args).print_stats(*print_args)
        return PrunWithOptions()

prun = Prun()

