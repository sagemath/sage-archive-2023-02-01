r"""
This file contains some tests that Sage command line options actually
do something.

We test the following command line options:

test.py
/path/to/test.py
test.sage
/path/to/test.sage
--advanced
--branch
-c
--cython
--ecl
--experimental
--gap
--gp
-h
--help
--hg
--ipython
--kash
--lisp
--maxima
--min
--mwrank
--optional
--preparse
--python
-q
--R
--root
--scons
--sh
--singular
--sqlite3
--standard
--startuptime
-t
-v
--zzfoobar (illegal option)


AUTHORS:

- Jeroen Demeyer (2010-11-20): initial version (#10300)

"""

from subprocess import *
import os, select

def test_executable(args, input="", timeout=50.0):
    """
    Run the program defined by ``args`` using the string ``input`` on
    the standard input.

    INPUT:

    - ``args`` -- a list of program arguments, the first being the
      executable.

    - ``input`` -- a string serving as standard input.  Usually, this
      should end with a newline.

    - ``timeout`` -- if the program produces no output for ``timeout``
      seconds, a RuntimeError is raised.

    OUTPUT: a tuple ``(out, err, ret)`` with the standard output,
    standard error and exitcode of the program run.


    EXAMPLES::

        sage: from sage.tests.cmdline import test_executable
        sage: (out, err, ret) = test_executable(["cat"], "Hello World!")
        sage: out
        'Hello World!'
        sage: err
        ''
        sage: ret
        0

    We test the timeout option::

        sage: (out, err, ret) = test_executable(["sleep", "1"], timeout=0.1)
        Traceback (most recent call last):
        ...
        RuntimeError: timeout in test_executable()

    TESTS:

    Run Sage itself with various options::

        sage: (out, err, ret) = test_executable(["sage"])
        sage: out.find(version()) >= 0
        True
        sage: err
        ''
        sage: ret
        0

        sage: (out, err, ret) = test_executable(["sage"], "3^33\n")
        sage: out.find(version()) >= 0
        True
        sage: out.find("5559060566555523") >= 0
        True
        sage: err
        ''
        sage: ret
        0

        sage: (out, err, ret) = test_executable(["sage", "-q"], "3^33\n")
        sage: out.find(version()) >= 0
        False
        sage: out.find("5559060566555523") >= 0
        True
        sage: err
        ''
        sage: ret
        0

        sage: (out, err, ret) = test_executable(["sage", "-c", "print 3^33"])
        sage: print out
        5559060566555523
        sage: err
        ''
        sage: ret
        0

        sage: (out, err, ret) = test_executable(["sage", "--min", "-c", "print 3^33"])
        sage: print out
        5559060566555523
        sage: err
        ''
        sage: ret
        0

        sage: (out, err, ret) = test_executable(["sage", "--startuptime"])
        sage: out.find("sage.all: ") >= 0
        True
        sage: err
        ''
        sage: ret
        0

    Test help::

        sage: (out, err, ret) = test_executable(["sage", "-h"])
        sage: out.find("Optional arguments:") >= 0
        True
        sage: err
        ''
        sage: ret
        0

        sage: (out, err, ret) = test_executable(["sage", "--help"])
        sage: out.find("Optional arguments:") >= 0
        True
        sage: err
        ''
        sage: ret
        0

        sage: (out, err, ret) = test_executable(["sage", "--advanced"])
        sage: out.find("run with no output prompts") >= 0
        True
        sage: err
        ''
        sage: ret
        0

    Basic information about the Sage installation::

        sage: (out, err, ret) = test_executable(["sage", "-v"])
        sage: out.find(version()) >= 0
        True
        sage: err
        ''
        sage: ret
        0

        sage: (out, err, ret) = test_executable(["sage", "--branch"])
        sage: len(out) >= 2   # at least one character + newline
        True
        sage: err
        ''
        sage: ret
        0

        sage: (out, err, ret) = test_executable(["sage", "--root"])
        sage: len(out) >= 2   # at least one character + newline
        True
        sage: err
        ''
        sage: ret
        0

    Test ``sage-run`` on a Python file, both with an absolute and with a relative path::

        sage: import tempfile
        sage: F = tempfile.NamedTemporaryFile(suffix=".py")
        sage: F.write("print 3^33\n")
        sage: F.flush()
        sage: (out, err, ret) = test_executable(["sage", F.name])
        sage: print out
        34
        sage: err
        ''
        sage: ret
        0
        sage: (dir,filename) = os.path.split(F.name)
        sage: os.chdir(dir)
        sage: (out, err, ret) = test_executable(["sage", filename])
        sage: print out
        34
        sage: err
        ''
        sage: ret
        0
        sage: del F   # Close and delete the file

    The same as above, but now with a ``.sage`` file.  This indirectly
    also tests the preparser::

        sage: import tempfile
        sage: F = tempfile.NamedTemporaryFile(suffix=".sage")
        sage: py_file = F.name[:-5] + ".py"  # Will be created by the preparser
        sage: F.write("k.<a> = GF(5^3); print a^124\n")
        sage: F.flush()
        sage: (out, err, ret) = test_executable(["sage", F.name])
        sage: os.unlink(py_file)
        sage: print out
        1
        sage: err
        ''
        sage: ret
        0
        sage: (dir,filename) = os.path.split(F.name)
        sage: os.chdir(dir)
        sage: (out, err, ret) = test_executable(["sage", filename])
        sage: os.unlink(py_file)
        sage: print out
        1
        sage: err
        ''
        sage: ret
        0
        sage: del F   # Close and delete the file

    Testing "sage --preparse FILE" and "sage -t FILE".  First create a file and preparse it::

        sage: import os
        sage: s = '\"\"\"\nThis is a test file.\n\"\"\"\ndef my_add(a,b):\n    \"\"\"\n    Add a to b.\n\n        EXAMPLES::\n\n            sage: my_add(2,2)\n            4\n        \"\"\"\n    return a+b\n'
        sage: script = os.path.join(SAGE_TMP, 'my_script.sage')
        sage: F = open(script, 'w')
        sage: F.write(s)
        sage: F.close()
        sage: os.chdir(SAGE_TMP)
        sage: (out, err, ret) = test_executable(["sage", "--preparse", script])
        sage: ret
        0
        sage: os.path.exists(os.path.join(SAGE_TMP, 'my_script.py'))
        True

    Now test my_script.sage and the preparsed version my_script.py::

        sage: (out, err, ret) = test_executable(["sage", "-t", script])
        sage: ret
        0
        sage: out.find("All tests passed!") >= 0
        True
        sage: (out, err, ret) = test_executable(["sage", "-t", os.path.join(SAGE_TMP, 'my_script.py')])
        sage: ret
        0
        sage: out.find("All tests passed!") >= 0
        True

    Now for a file which should fail tests.  Run this test with
    SAGE_TESTDIR equal to a temporary directory, because failed doctests
    leave files lying around in SAGE_TESTDIR::

        sage: s = s.replace('4', '5') # (2+2 != 5)
        sage: F = open(script, 'w')
        sage: F.write(s)
        sage: F.close()
        sage: OLD_TESTDIR = os.environ['SAGE_TESTDIR']
        sage: os.environ['SAGE_TESTDIR'] = SAGE_TMP
        sage: (out, err, ret) = test_executable(["sage", "-t", script])
        sage: ret
        128
        sage: out.find("1 items had failures:") >= 0
        True
        sage: os.environ['SAGE_TESTDIR'] = OLD_TESTDIR  # just in case

    Test external programs being called by Sage::

        sage: (out, err, ret) = test_executable(["sage", "--sh"], "echo Hello World\nexit 42\n")
        sage: out.find("Hello World\n") >= 0
        True
        sage: ret
        42

        sage: (out, err, ret) = test_executable(["sage", "--sh", "-c", "echo Hello World; exit 42"])
        sage: out.find("Hello World\n") >= 0
        True
        sage: ret
        42

    When testing IPython, use Sage's ipython directory, to avoid
    incompatibilities for config files among different versions of
    IPython. ::

        sage: os.environ['IPYTHONDIR'] = os.path.join(os.environ['DOT_SAGE'], 'ipython')
        sage: (out, err, ret) = test_executable(["sage", "--ipython"], "\n3**33\n")
        sage: out.find("5559060566555523") >= 0
        True
        sage: err
        ''
        sage: ret
        0

        sage: (out, err, ret) = test_executable(["sage", "--python"], "print 3^33\n")
        sage: out
        '34\n'
        sage: err
        ''
        sage: ret
        0

        sage: (out, err, ret) = test_executable(["sage", "--cython"])
        sage: print err
        Cython (http://cython.org) is a compiler for code written in the
        Cython language.  Cython is based on Pyrex by Greg Ewing.
        ...

        sage: (out, err, ret) = test_executable(["sage", "--ecl"], "(* 12345 54321)\n")
        sage: out.find("Embeddable Common-Lisp") >= 0
        True
        sage: out.find("670592745") >= 0
        True
        sage: err
        ''
        sage: ret
        0

        sage: (out, err, ret) = test_executable(["sage", "--lisp"], "(* 12345 54321)\n")
        sage: out.find("Embeddable Common-Lisp") >= 0
        True
        sage: out.find("670592745") >= 0
        True
        sage: err
        ''
        sage: ret
        0

        sage: (out, err, ret) = test_executable(["sage", "--gap", "-q"], "Size(SymmetricGroup(5));\n")
        sage: out
        '120\n'
        sage: err
        ''
        sage: ret
        0

        sage: (out, err, ret) = test_executable(["sage", "--kash", "-b", "3^33;\n"])  # optional - kash
        sage: out.find("5559060566555523") >= 0  # optional - kash
        True
        sage: err  # optional - kash
        ''
        sage: ret  # optional - kash
        0

        sage: (out, err, ret) = test_executable(["sage", "--mwrank", "-v0", "-q", "-o"], "0 0 1 -7 6 0 0 0 0 0\n")
        sage: out
        'Curve [0,0,1,-7,6] :\tRank = 3\n[[3],[[1,-1],[-2,3],[-7/4,25/8]]]\n\n\n'
        sage: err
        ''
        sage: ret
        0

        sage: (out, err, ret) = test_executable(["sage", "--singular"], "12345*54321;\n")
        sage: out.find("A Computer Algebra System for Polynomial Computations") >= 0
        True
        sage: out.find("670592745") >= 0
        True
        sage: err
        ''
        sage: ret
        0

    Check that ``sage-make_relative`` did its job.  We test it on the
    ``ipython`` script::

        sage: open(os.path.join(SAGE_ROOT, "local", "bin", "ipython")).readline()
        '#!/usr/bin/env python\n'

    Test GP using the ``-f`` option which prevents the reading of a ``.gprc``
    configuration file::

        sage: (out, err, ret) = test_executable(["sage", "--gp", "-f"], "3^33\nquit(42)\n")
        sage: out.find("PARI/GP") >= 0
        True
        sage: out.find("5559060566555523") >= 0
        True
        sage: err
        ''
        sage: ret
        42

    Some programs of which we check functionality using only ``--version``::

        sage: (out, err, ret) = test_executable(["sage", "--hg", "--version"])
        sage: out.find("Mercurial Distributed SCM") >= 0
        True
        sage: err
        ''
        sage: ret
        0

        sage: (out, err, ret) = test_executable(["sage", "--maxima", "--version"])
        sage: out.find("Maxima ") >= 0
        True
        sage: err
        ''
        sage: ret
        0

        sage: (out, err, ret) = test_executable(["sage", "--R", "--version"])
        sage: out.find("R version ") >= 0
        True
        sage: err
        ''
        sage: ret
        0

        sage: (out, err, ret) = test_executable(["sage", "--scons", "--version"])
        sage: out.find("SCons") >= 0
        True
        sage: err
        ''
        sage: ret
        0

        sage: (out, err, ret) = test_executable(["sage", "--sqlite3", "--version"])
        sage: out.startswith("3.")
        True
        sage: err
        ''
        sage: ret
        0

    Check some things requiring an internet connection::

        sage: (out, err, ret) = test_executable(["sage", "--standard"])  # optional - internet
        sage: out.find("atlas") >= 0  # optional - internet
        True
        sage: err  # optional - internet
        ''
        sage: ret  # optional - internet
        0

        sage: (out, err, ret) = test_executable(["sage", "--optional"])  # optional - internet
        sage: out.find("database_cremona_ellcurve") >= 0  # optional - internet
        True
        sage: err  # optional - internet
        ''
        sage: ret  # optional - internet
        0

        sage: (out, err, ret) = test_executable(["sage", "--experimental"])  # optional - internet
        sage: out.find("macaulay2") >= 0  # optional - internet
        True
        sage: err  # optional - internet
        ''
        sage: ret  # optional - internet
        0

    Check an illegal command line option.  This outputs an error to stdout,
    but we allow stderr in case this changes in the future::

        sage: (out, err, ret) = test_executable(["sage", "--zzfoobar"])
        sage: (out+err).find("unknown option: --zzfoobar") >= 0
        True
        sage: ret > 0
        True

    """
    p = Popen(args, stdin=PIPE, stdout=PIPE, stderr=PIPE)
    if input: p.stdin.write(input)
    p.stdin.close()
    fdout = p.stdout.fileno()
    fderr = p.stderr.fileno()
    out = ""; err = ""

    while True:
        # Try reading from fdout and fderr
        rfd = []
        if fdout: rfd.append(fdout)
        if fderr: rfd.append(fderr)
        if len(rfd) == 0: break
        rlist = select.select(rfd, [], [], timeout)[0]

        if len(rlist) == 0:
            # Timeout!
            p.terminate()
            raise RuntimeError("timeout in test_executable()")
        if fdout in rlist:
            s = os.read(fdout, 1024)
            if s == "": fdout = None   # EOF
            out += s
        if fderr in rlist:
            s = os.read(fderr, 1024)
            if s == "": fderr = None   # EOF
            err += s

    return (out, err, p.wait())
