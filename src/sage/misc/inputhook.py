"""
The Sage Input Hook

This is a hook into the IPython input prompt and will be called
periodically (every 100ms) while Python is sitting idle. We use it to
reload attached files if they have changed.
"""

###########################################################################
#       Copyright (C) 2013 Volker Braun <vbraun.name@gmail.com>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#                  http://www.gnu.org/licenses/
###########################################################################


from sage.misc.attached_files import reload_attached_files_if_modified


def sage_inputhook():
    """
    The input hook.

    This function will be called every 100ms when IPython is idle at
    the command prompt.

    EXAMPLES::

        sage: from sage.misc.interpreter import get_test_shell
        sage: shell = get_test_shell()
        sage: tmp = tmp_filename(ext='.py')
        sage: f = open(tmp, 'w'); f.write('a = 2\n'); f.close()
        sage: shell.run_cell('%attach ' + tmp)
        sage: shell.run_cell('a')
        2
        sage: f = open(tmp, 'w'); f.write('a = 3\n'); f.close()

    Note that the doctests are never really at the command prompt, so
    we call the input hook manually::

        sage: shell.run_cell('from sage.misc.inputhook import sage_inputhook')
        sage: shell.run_cell('sage_inputhook()')
        ### reloading attached file tmp_....py modified at ... ###
        0

        sage: shell.run_cell('a')
        3
        sage: shell.run_cell('detach({0})'.format(repr(tmp)))
        sage: shell.run_cell('attached_files()')
        []
    """
    reload_attached_files_if_modified()
    return 0





