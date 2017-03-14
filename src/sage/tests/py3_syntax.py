"""
Test that the Sage library uses Python 3 syntax

EXAMPLES::

    sage: from sage.tests.py3_syntax import Python3SyntaxTest
    sage: py3_syntax = Python3SyntaxTest('sage', 'sage_setup')
    sage: py3_syntax.run_tests('.py')   # long time
"""

from __future__ import print_function

import os
import itertools
import subprocess

from sage.env import SAGE_SRC


class SortedDirectoryWalkerABC(object):

    def __init__(self, *directories):
        """
        Walk directory tree in a reproducible manner

        INPUT:

        -- ``*directories`` - Roots of the directory trees to walk.

        EXAMPLES::

            sage: from sage.tests.py3_syntax import Python3SyntaxTest
            sage: Python3SyntaxTest('sage/tests')
            <sage.tests.py3_syntax.Python3SyntaxTest object at 0x...>
        """
        self._directories = tuple(
            os.path.join(SAGE_SRC, d) for d in directories
        )

    def __iter__(self):
        """
        Iterate over files

        OUTPUT:

        Iterator over triples consisting of path, filename, and file
        extension. The iteration order only depends on the filenames
        and not on filesystem layout.

        EXAMPLES::

            sage: from sage.tests.py3_syntax import Python3SyntaxTest
            sage: test = Python3SyntaxTest('sage/tests/french_book')
            sage: next(iter(test))
            ('src/sage/tests/french_book', 'README', '')
        """
        tree_walk = itertools.chain(*map(os.walk, self._directories))
        for path, _, files in tree_walk:
            path = os.path.relpath(path)
            files.sort()
            for filename in files:
                if filename.startswith('.'):
                    continue
                _, ext = os.path.splitext(filename)
                yield (path, filename, ext)

    def run_tests(self, *extensions):
        """
        Run tests

        The abstract :meth:`test` is called on each file with a
        matching extension.

        INPUT:

        -- ``*extensions`` - the file extensions (including the leading dot)
            to check.

        EXAMPLES::

            sage: from sage.tests.py3_syntax import SortedDirectoryWalkerABC
            sage: walker = SortedDirectoryWalkerABC('sage/tests/french_book')
            sage: walker.run_tests('.py')
            Traceback (most recent call last):
            ...
            NotImplementedError: to be implemented in derived class
        """
        buf = []
        for (path, filename, ext) in self:
            if ext in extensions:
                fqn = os.path.join(path, filename)
                buf.append(fqn)
                if len(buf) > 20:
                    self.test(*buf)
                    buf = []
        self.test(*buf)

    def test(self, *filenames):
        """
        Test the given filenames

        To be implemented in derived classes.

        INPUT:

        -- ``*filename`` -- string. The full qualified filenames to check.

        EXAMPLES::

            sage: from sage.tests.py3_syntax import SortedDirectoryWalkerABC
            sage: walker = SortedDirectoryWalkerABC()
            sage: from sage.env import SAGE_SRC
            sage: filename = os.path.join(SAGE_SRC, 'sage', 'tests', 'py3_syntax.py')
            sage: walker.test(filename)
            Traceback (most recent call last):
            ...
            NotImplementedError: to be implemented in derived class
        """
        raise NotImplementedError('to be implemented in derived class')


PYTHON3_ENV = dict(os.environ)
PYTHON3_ENV.pop('PYTHONPATH', None)


class Python3SyntaxTest(SortedDirectoryWalkerABC):

    def test(self, *filenames):
        r"""
        Test that the given filenames are valid Python 3 syntax

        INPUT:

        -- ``filename`` -- string. The full qualified filename to check.

        EXAMPLES::

            sage: import os, tempfile
            sage: src = tempfile.NamedTemporaryFile(suffix='.py', delete=False)
            sage: src.write('print "invalid print statement"')
            sage: src.close()
            sage: from sage.tests.py3_syntax import Python3SyntaxTest
            sage: py3_syntax = Python3SyntaxTest()
            sage: py3_syntax.test(src.name)
            Invalid Python 3 syntax found:
              File "...py", line 1
                print "invalid print statement"
                                              ^
            SyntaxError: Missing parentheses in call to 'print'
            sage: os.unlink(src.name)
        """
        cmd = [
            'python3',
            '-m', 'py_compile'
        ] + list(filenames)
        process = subprocess.Popen(
            cmd,
            env=PYTHON3_ENV,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE
        )
        stdout, stderr = process.communicate()
        if process.returncode == 0:
            return
        print('Invalid Python 3 syntax found:')
        if stdout:
            print(stdout)
        if stderr:
            print(stderr)
