import __builtin__
import os

from sage.misc.misc import tmp_filename

count=0


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

    def eval(self,x,globals=None, locals=None):
        """
        EXAMPLES::

            sage: from sage.misc.inline_fortran import InlineFortran, _example
            sage: test_fortran = InlineFortran(globals())   # optional -- fortran
            sage: test_fortran(_example)                    # optional -- fortran
            sage: import numpy
            sage: n = numpy.array(range(10),dtype=float)
            sage: fib(n,int(10))                            # optional -- fortran
            sage: n                                         # optional -- fortran
            array([  0.,   1.,   1.,   2.,   3.,   5.,   8.,  13.,  21.,  34.])
        """
        if len(x.splitlines()) == 1 and os.path.exists(x):
            filename = x
            x = open(x).read()
            if filename.lower().endswith('.f90'):
                x = '!f90\n' + x
        global count
        # On linux g77_shared should be a script that runs sage_fortran -shared
        # On OS X it should be a script that runs gfortran -bundle -undefined dynamic_lookup
        path = os.environ['SAGE_LOCAL']+'/bin/sage-g77_shared'
        from numpy import f2py
        old_import_path=os.sys.path
        cwd=os.getcwd()
        os.sys.path.append(cwd)

        #name = tmp_dir() + '/fortran_module_%d'%count
        name = 'fortran_module_%d'%count
        if os.path.exists(name):
            os.unlink(name)
        s_lib_path=""
        s_lib=""
        for s in self.library_paths:
            s_lib_path=s_lib_path+"-L"+s+" "

        for s in self.libraries:
            s_lib=s_lib +"-l"+s + " "

        # if the first line has !f90 as a comment gfortran will treat it as
        # fortran 90 code
        if x.startswith('!f90'):
            fname = os.path.join(tmp_filename() +'.f90')
        else:
            fname = os.path.join(tmp_filename() +'.f')

        log = tmp_filename()
        extra_args = '--quiet --f77exec=%s --f90exec=%s %s %s  1>&2 >"%s"'%(
                    path, path, s_lib_path, s_lib, log)

        f2py.compile(x, name, extra_args = extra_args, source_fn=fname)

        log_string = open(log).read()

        os.unlink(log)
        os.unlink(fname)

        if self.verbose:
            print log_string

        count += 1
        try:
            m=__builtin__.__import__(name)
        except ImportError:
            if not self.verbose:
                print log_string
            return
        finally:
            os.sys.path=old_import_path
            os.unlink(name + '.so')

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
