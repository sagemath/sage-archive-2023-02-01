"""
Fortran compiler
"""
import __builtin__
import os

from sage.misc.temporary_file import tmp_dir


class InlineFortran:
    def __init__(self,globals):
        self.globals=globals
        self.library_paths=[]
        self.libraries=[]
        self.verbose = False

    def __repr__(self):
        return "Interface to Fortran compiler"

    def __call__(self, *args, **kwds):
        return self.eval(*args, **kwds)

    def eval(self, x, globals=None, locals=None):
        """
        EXAMPLES::

            sage: from sage.misc.inline_fortran import InlineFortran, _example
            sage: fortran = InlineFortran(globals())
            sage: fortran(_example)
            sage: import numpy
            sage: n = numpy.array(range(10),dtype=float)
            sage: fib(n,int(10))
            sage: n
            array([  0.,   1.,   1.,   2.,   3.,   5.,   8.,  13.,  21.,  34.])

        TESTS::

            sage: os.chdir(SAGE_ROOT)
            sage: fortran.eval("SYNTAX ERROR !@#$")
            Traceback (most recent call last):
            ...
            RuntimeError: failed to compile Fortran code:...
            sage: os.getcwd() == SAGE_ROOT
            True
        """
        if len(x.splitlines()) == 1 and os.path.exists(x):
            filename = x
            x = open(x).read()
            if filename.lower().endswith('.f90'):
                x = '!f90\n' + x

        from numpy import f2py

        # Create everything in a temporary directory
        mytmpdir = tmp_dir()

        try:
            old_cwd = os.getcwd()
            os.chdir(mytmpdir)

            old_import_path = os.sys.path
            os.sys.path.append(mytmpdir)

            name = "fortran_module"  # Python module name
            # if the first line has !f90 as a comment, gfortran will
            # treat it as Fortran 90 code
            if x.startswith('!f90'):
                fortran_file = name + '.f90'
            else:
                fortran_file = name + '.f'

            s_lib_path = ""
            s_lib = ""
            for s in self.library_paths:
                s_lib_path = s_lib_path + "-L%s "

            for s in self.libraries:
                s_lib = s_lib + "-l%s "%s

            log = name + ".log"
            extra_args = '--quiet --f77exec=sage-inline-fortran --f90exec=sage-inline-fortran %s %s >"%s" 2>&1'%(
                s_lib_path, s_lib, log)

            f2py.compile(x, name, extra_args = extra_args, source_fn=fortran_file)
            log_string = open(log).read()

            # f2py.compile() doesn't raise any exception if it fails.
            # So we manually check whether the compiled file exists.
            # NOTE: the .so extension is used expect on Cygwin,
            # that is even on OS X where .dylib might be expected.
            soname = name
            uname = os.uname()[0].lower()
            if uname[:6] == "cygwin":
                soname += '.dll'
            else:
                soname += '.so'
            if not os.path.isfile(soname):
                raise RuntimeError("failed to compile Fortran code:\n" + log_string)

            if self.verbose:
                print(log_string)

            m = __builtin__.__import__(name)
        finally:
            os.sys.path = old_import_path
            os.chdir(old_cwd)
            try:
                import shutil
                shutil.rmtree(mytmpdir)
            except OSError:
                # This can fail for example over NFS
                pass

        for k, x in m.__dict__.iteritems():
            if k[0] != '_':
                self.globals[k] = x

    def add_library(self,s):
       self.libraries.append(s)

    def add_library_path(self,s):
       self.library_paths.append(s)


_example = """
C FILE: FIB1.F
      SUBROUTINE FIB(A,N)
C
C     CALCULATE FIRST N FIBONACCI NUMBERS
C
      INTEGER N
      REAL*8 A(N)
      DO I=1,N
         IF (I.EQ.1) THEN
            A(I) = 0.0D0
         ELSEIF (I.EQ.2) THEN
            A(I) = 1.0D0
         ELSE
            A(I) = A(I-1) + A(I-2)
         ENDIF
      ENDDO
      END
C END FILE FIB1.F
"""
