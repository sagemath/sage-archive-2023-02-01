import os
if 'SAGE_CLEAN' not in os.environ:
    import sage.misc.misc
    from sage.misc.interpreter import preparser, _ip
    preparser(True)

    import sage.all_cmdline
    sage.all_cmdline._init_cmdline(globals())

    _ip.ex('from sage.all import Integer, RealNumber')
    os.chdir(os.environ["CUR"])
    import sage.misc.interpreter

    from sage.misc.interpreter import attached_files

    branch = sage.misc.misc.branch_current_hg_notice(sage.misc.misc.branch_current_hg())
    if branch:
        print branch

    if not os.environ.has_key('SAGE_IMPORTALL') or os.environ['SAGE_IMPORTALL'] != "no":
        _ip.ex('from sage.all_cmdline import *')


    startup_file = os.environ.get('SAGE_STARTUP_FILE', '')
    if os.path.exists(startup_file):
        _ip.options.autoexec.append('load %s'%startup_file)

    from sage.misc.sage_timeit import sage_timeit
    _ip.expose_magic('timeit', lambda self, s: sage_timeit(s, _ip.user_ns))

    from sage.misc.preparser import preparse
    old_prun = _ip.IP.magic_prun
    _ip.expose_magic('prun', lambda self, s: old_prun(preparse(s)))
