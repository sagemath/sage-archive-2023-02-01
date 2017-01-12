def _get_pari_instance():
    # There are two constraints for the virtual stack size:
    # 1) on 32-bit systems, even virtual memory can be a scarce
    #    resource since it is limited by 4GB (of which the kernel
    #    needs a significant part)
    # 2) the system should actually be able to handle a stack size
    #    as large as the complete virtual stack.
    # As a simple heuristic, we set the virtual stack to 1/4 of the
    # virtual memory.
    from sage.misc.getusage import virtual_memory_limit

    sizemax = virtual_memory_limit() // 4

    from sage.env import CYGWIN_VERSION
    if CYGWIN_VERSION and CYGWIN_VERSION < (2, 5, 2):
        # Cygwin's mmap is broken for large NORESERVE mmaps (>~ 4GB) See
        # http://trac.sagemath.org/ticket/20463 So we set the max stack
        # size to a little below 4GB (putting it right on the margin proves
        # too fragile)
        #
        # The underlying issue is fixed in Cygwin v2.5.2
        sizemax = min(sizemax, 0xf0000000)

    from sage.libs.cypari2 import Pari
    P = Pari(1000000, sizemax)

    # pari_init_opts() overrides MPIR's memory allocation functions,
    # so we need to reset them.
    from sage.ext.memory import init_memory_functions
    init_memory_functions()

    return P

pari = _get_pari_instance()
