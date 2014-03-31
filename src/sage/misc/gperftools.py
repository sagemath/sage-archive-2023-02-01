"""
C Function Profiler Using Google Perftools

Note that the profiler samples 100x per second by default. In
particular, you cannot profile anything shorter than 10ms. You can
adjust the rate with the ``CPUPROFILE_FREQUENCY`` evironment variable
if you want to change it.

EXAMPLES::

    sage: from sage.misc.crun import Profiler
    sage: prof = Profiler()
    sage: prof.start()       # optional - gperftools
    sage: prof.stop()        # optional - gperftools
    PROFILE: interrupts/evictions/bytes = ...

REFERENCE:

Uses the `Google performance analysis tools
<https://code.google.com/p/gperftools>`_. Note that they are not
included in Sage, you have to install them yourself on your system.

AUTHORS:

- Volker Braun (2014-03-31): initial version
"""

#*****************************************************************************
#       Copyright (C) 2014 Volker Braun <vbraun.name@gmail.com>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#
#                  http://www.gnu.org/licenses/
#*****************************************************************************

from sage.structure.sage_object import SageObject


libc = None
libprofiler = None



class Profiler(SageObject):

    def __init__(self, filename=None):
        """
        Interface to the gperftools profiler
        
        INPUT:
    
        - ``filename`` -- string or ``None`` (default). The file name to log to. By default, a if 

        EXAMPLES::

            sage: from sage.misc.crun import Profiler
            sage: Profiler()
            Profiler logging to ...
        """        
        if filename is None:
            from sage.misc.temporary_file import tmp_filename
            self._filename = tmp_filename(ext='.perf')
        else:
            self._filename = filename

    def filename(self):
        """
        Return the file name

        OUTPUT:

        String.

        EXAMPLES::

            sage: from sage.misc.crun import Profiler
            sage: prof = Profiler()
            sage: prof.filename()
            '.../tmp_....perf'
        """
        return self._filename

    def _repr_(self):
        """
        Return string representation

        OUTPUT:

        String.

        EXAMPLES::

            sage: from sage.misc.crun import Profiler
            sage: Profiler()
            Profiler logging to .../tmp....perf
        """
        return 'Profiler logging to {0}'.format(self.filename())

    def _libc(self):
        """
        Return libc
        
        OUTPUT:

        A ctypes shared library handle.

        EXAMPLES::

            sage: from sage.misc.crun import Profiler
            sage: Profiler()._libc()
            <CDLL 'libc...', handle ... at ...>
        """
        global libc
        if libc is not None:
            return libc
        import ctypes, ctypes.util
        name = ctypes.util.find_library('c')
        if name:
            libc = ctypes.CDLL(name)
            return libc
        else:
            raise ImportError('failed to open libc')

    def _libprofiler(self):
        """
        Return libprofiler
        
        OUTPUT:

        A ctypes shared library handle.

        EXAMPLES::

            sage: from sage.misc.crun import Profiler
            sage: Profiler()._libprofiler()    # optional - gperftools
            <CDLL 'libprofiler...', handle ... at ...>
        """
        global libprofiler
        if libprofiler is not None:
            return libprofiler
        import ctypes, ctypes.util
        name = ctypes.util.find_library('profiler')
        if name:
            libprofiler = ctypes.CDLL(name)
            return libprofiler
        else:
            raise ImportError('failed to open libprofiler, make sure gperftools is installed')

    def start(self):
        """
        Start profiling
        
        EXAMPLES::

            sage: from sage.misc.crun import Profiler
            sage: prof = Profiler()
            sage: prof.start()    # optional - gperftools
            sage: # do something
            sage: prof.stop()     # optional - gperftools
            PROFILE: interrupts/evictions/bytes = ...
        """
        from signal import SIGPROF, SIG_DFL
        self._previous_sigprof_handler = self._libc().signal(SIGPROF, SIG_DFL)
        profiler = self._libprofiler()
        rc = profiler.ProfilerStart(self.filename())
        if rc < 0:
            raise ValueError('profiler failed to start')

    def stop(self):
        """
        Stop the CPU profiler
        
        EXAMPLES::

            sage: from sage.misc.crun import Profiler
            sage: prof = Profiler()
            sage: prof.start()    # optional - gperftools
            sage: # do something
            sage: prof.stop()     # optional - gperftools
            PROFILE: interrupts/evictions/bytes = ...
        """
        profiler = self._libprofiler()
        profiler.ProfilerStop()

    def _pprof(self):
        """
        Return the name of the ``pprof`` binary.

        OUTPUT:

        String.

        EXAMPLES::

            sage: from sage.misc.crun import Profiler
            sage: prof = Profiler()
            sage: prof._pprof()
            'pprof'
        """
        return 'pprof'
        
    def _executable(self):
        """
        Return the name of the Sage Python interpreter.

        OUTPUT:

        String.

        EXAMPLES::

            sage: from sage.misc.crun import Profiler
            sage: prof = Profiler()
            sage: prof._executable()
            '.../python'
        """
        import sys
        return sys.executable

    def _call_pprof(self, *args, **kwds):
        """
        Run the pprof binary

        INPUT:

        - ``args`` -- list of strings. The arguments to ``pprof``.

        - ``kwds`` -- keyword arguments passed to
          ``subprocess.check_call``.

        EXAMPLES::

            sage: from sage.misc.crun import Profiler
            sage: prof = Profiler()
            sage: prof._call_pprof('--help')   # optional - gperftools
            Usage:
            pprof [options] <program> <profiles>
            ...
        """
        from subprocess import check_call
        check_call([self._pprof()] + list(args), **kwds)

    def top(self, cumulative=False):
        """
        Print text report

        OUTPUT:

        Nothing. A textual report is printed to stdout.

        EXAMPLES::

            sage: from sage.misc.crun import Profiler
            sage: prof = Profiler()
            sage: prof.start()    # optional - gperftools
            sage: # do something
            sage: prof.stop()     # optional - gperftools
            PROFILE: interrupts/evictions/bytes = ...
            sage: prof.top()   # optional - gperftools
            Using local file ...
            Using local file ...
        """
        args = []
        if cumulative:
            args += ['--cum']
        args += ['--text', self._executable(), self.filename()]
        self._call_pprof(*args)
            
    def save(self, filename, cumulative=False, verbose=True):
        """
        Save report to disk.

        INPUT:

        - ``filename`` -- string. The filename to save at. Must end
          with one of ``.dot``, ``.ps``, ``.pdf``, ``.svg``, ``.gif``,
          or ``.txt`` to specify the output file format.
        
        - ``cumulative`` -- boolean (optional, default:
          ``False``). Whether to return cumulative timings.

        - ``verbose`` -- boolean (optional, default:
          ``True``). Whether to print informational messages.

        EXAMPLES::

            sage: from sage.misc.crun import Profiler
            sage: prof = Profiler()
            sage: prof.start()    # optional - gperftools
            sage: prof.stop()     # optional - gperftools
            PROFILE: interrupts/evictions/bytes = ...
            sage: f = tmp_filename(ext='.txt')
            sage: prof.save(f, verbose=False)       # optional - gperftools
        """
        args = []
        if cumulative:
            args += ['--cum']
        if filename.endswith('.dot'):
            args += ['--dot']
        elif filename.endswith('.ps'):
            args += ['--ps']
        elif filename.endswith('.pdf'):
            args += ['--pdf']
        elif filename.endswith('.svg'):
            args += ['--svg']
        elif filename.endswith('.txt'):
            args += ['--text']
        elif filename.endswith('.gif'):
            args += ['--gif']
        else:
            raise ValueError('unknown extension')
        args += [self._executable(), self.filename()]
        stderr = sys.stdout if verbose else False
        with open(filename, 'wb') as outfile:
            self._call_pprof(*args, stdout=outfile, stderr=stderr)


def crun(s, evaluator):
    """
    Profile single statement.

    - ``s`` -- string. Sage code to profile.

    - ``evaluator`` -- callable to evaluate.

    EXAMPLES::

        sage: import sage.misc.gperftools
        sage: ev = lambda ex:eval(ex, globals(), locals())
        sage: sage.misc.gperftools.crun('1+1', evaluator=ev)   # optional - gperftools
        PROFILE: interrupts/evictions/bytes = ...
        Using local file ...
        Using local file ...
    """
    prof = Profiler()
    from sage.misc.preparser import preparse
    py_s = preparse(s)
    prof.start()
    evaluator(py_s)
    prof.stop()
    prof.top()

