import sage.misc.misc
import sage.server.support

def cython(code,
          verbose=False, compile_message=False,
          make_c_file_nice=False, use_cache=False):
    """
    Given a block of Cython (Sage's variant of Pyrex) code (as a text
    string), this function compiles it using your C compiler, and
    includes it into the global scope of the module that called this
    function.

    WARNING: Only use this from Python code, not from extension code,
    since from extension code you would change the global scope (i.e.,
    of the Sage interpreter). And it would be stupid, since you're
    already writing Cython!

    Also, never use this in the standard Sage library.  Any code that
    uses this can only run on a system that has a C compiler
    installed, and we want to avoid making that assumption for casual
    Sage usage.  Also, any code that uses this in the library would
    greatly slow down startup time, since currently there is no
    caching.

    AUTHOR: William Stein, 2006-10-31

    TODO: Need to create a clever caching system so code only gets
    compiled once.
    """
    tmpfile = sage.misc.misc.tmp_filename() + ".spyx"
    open(tmpfile,'w').write(code)
    sage.server.support.cython_import_all(tmpfile, globals(),
                                         verbose=verbose, compile_message=compile_message,
                                         use_cache=use_cache,
                                         create_local_c_file=False)


