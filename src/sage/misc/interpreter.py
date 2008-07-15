"""nodoctest
Preparses input from the interpreter

Modified input.

  -- All ^'s (not in strings) are replaced by **'s.

  -- If M is a variable and i an integer,
     then M.i is replaced by M.gen(i), so generators can be
     accessed as in MAGMA.

  -- quit alone on a line quits.

  -- load to load in scripts

  -- Most int literals n are replaced by ZZ(n) Thus 2/3 is a rational
     number.  If they are in []'s right after a valid identifier they
     aren't replaced.

  -- real literals get wrapped in "RR" (with a precision)

  -- the R.<x,y,z> = ... notation

TODO:
  I have no plans for anything further, except to improve the
  robustness of the above.  Preparsing may work incorrectly for
  multi-line input lines in some cases; this will be fixed.

All other input is processed normally.

It automatically converts *most* integer literals to SAGE Integer's and
decimal literals to SAGE Reals.  It does not convert indexes into
1-d arrays, since those have to be ints.

I also extended the load command so it *really* works exactly
like the SAGE interpreter, so e.g., ^ for exponentiation is
allowed.  Also files being loaded can themselves load other files.
Finally, I added an "attach" command, e.g.,
    attach 'file'
that works kind of like attach in MAGMA.  Whenever you enter a blank
line in the SAGE interpreter, *all* attached files that have changed
are automatically reloaded.  Moreover, the attached files work according
to the SAGE interpreter rules, i.e., ^ --> **, etc.

I also fixed it so ^ is not replaced by ** inside strings.

Finally, I added back the M.n notation for the n-th generator
of object M, again like in MAGMA.

EXAMPLE:
    sage: 2/3
          2/3
    sage: type(2/3)
          <type 'sage.rings.rational.Rational'>
    sage: a = 49928420832092
    sage: type(a)
          <type 'sage.rings.integer.Integer'>
    sage: a.factor()
          [(2, 2), (11, 1), (1134736837093, 1)]
    sage: v = [1,2,3]
    sage: type(v[0])
          <type 'sage.rings.integer.Integer'>

If we don't make potential list indices int's, then lots of stuff
breaks, or users have to type v[int(7)], which is insane.
A fix would be to only not make what's in the brackets an
Integer if what's before the bracket is a valid identifier,
so the w = [5] above would work right.

    sage: s = "x^3 + x + 1"
    sage: s
          'x^3 + x + 1'
    sage: pari(s)
          x^3 + x + 1
    sage: f = pari(s)
    sage: f^2
          x^6 + 2*x^4 + 2*x^3 + x^2 + 2*x + 1
    sage: V = VectorSp
    VectorSpace                      VectorSpace_subspace
    VectorSpace_ambient              VectorSpace_subspace_with_basis
    VectorSpace_generic
    sage: V = VectorSpace(Q,3)
    sage: V.0
          (1 0 0)
    sage: V.1
          (0 1 0)
    sage: s = "This. Is. It."
    sage: print s
    This. Is. It.
"""

#*****************************************************************************
#       Copyright (C) 2004 William Stein <wstein@gmail.com>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#
#                  http://www.gnu.org/licenses/
#*****************************************************************************
import IPython.ipapi
_ip = IPython.ipapi.get()

__author__ = 'William Stein <wstein@gmail.com> et al.'
__license__ = 'GPL'

import os
import log

import remote_file

from IPython.iplib import InteractiveShell

import preparser_ipython
from preparser import preparse_file

import cython

#import signal
#def sig(x,n):
#    raise KeyboardInterrupt
#def unsetsig():
#    signal.signal(signal.SIGINT, sig)

# IPython has a prefilter() function that analyzes each input line. We redefine
# it here to first pre-process certain forms of input

# The prototype of any alternate prefilter must be like this one (the name
# doesn't matter):
# - line is a string containing the user input line.
# - continuation is a parameter which tells us if we are processing a first line of
#   user input or the second or higher of a multi-line statement.



attached = { }

def attached_files():
    """
    Return a list of all files attached to the current session.
    """
    global attached
    X = attached.keys()
    X.sort()
    return X

def load_startup_file(file):
    if os.path.exists(file):
        X = do_prefilter_paste('load "%s"'%file,False)
        _ip.runlines(X)
    if os.path.exists('attach.sage'):
        X = do_prefilter_paste('attach "attach.sage"',False)
        _ip.runlines(X)


def do_prefilter_paste(line, continuation):
    """
    Alternate prefilter for input.

    INPUT:
        line -- a single line; must *not* have any newlines in it
        continuation -- whether the input line is really part
                     of the previous line, because of open parens or backslash.
    """
    if '\n' in line:
        raise RuntimeError, "bug in function that calls do_prefilter_paste -- there can be no newlines in the input"

    global attached

    # This is so it's OK to have lots of blank space at the
    # beginning of any non-continuation line.

    if continuation:
        # strip ...'s that appear in examples
        L = line.lstrip()
        if L[:3] == '...':
            line = L[3:]
    else:
        line = line.lstrip()

    line = line.rstrip()

    if not line.startswith('attach ') and not line.startswith('load ') and not line.startswith('%run '):
        for F in attached.keys():
            tm = attached[F]
            if os.path.exists(F) and os.path.getmtime(F) > tm:
                # Reload F.
                try:
                    if F.endswith('.py'):
                        _ip.runlines('%%run -i "%s"'%F)
                    elif F.endswith('.sage'):
                        _ip.runlines('%%run -i "%s"'%preparse_file_named(F))
                    elif F.endswith('.spyx') or F.endswith('.pyx'):
                        X = load_cython(F)
                        __IPYTHON__.push(X)
                    else:
                        line = 'load("%s")'%F
                    t = os.path.getmtime(F)
                    attached[F] = t
                except IOError:
                    del attached[F]


    # Get rid of leading sage: so that pasting of examples from the documentation
    # works.  This is like MAGMA's SetLinePrompt(false).
    for prompt in ['sage:', '>>>']:
        if not continuation:
            while True:
                strip = False
                if line[:3] == prompt:
                    line = line[3:].lstrip()
                    strip = True
                elif line[:5] == prompt:
                    line = line[5:].lstrip()
                    strip = True
                if not strip:
                    break
                else:
                    line = line.lstrip()

    # 'quit' alone on a line to quit.
    if line.lower() in ['quit', 'exit', 'quit;', 'exit;']:
        line = '%quit'

    #################################################################
    # An interactive load command, like iload in MAGMA.
    #################################################################
    if line[:6] == 'iload ':
        try:
            name = str(eval(line[6:]))
        except:
            name = str(line[6:].strip())
        try:
            F = open(name)
        except IOError:
            raise ImportError, 'Could not open file "%s"'%name

        print 'Interactively loading "%s"'%name
        n = len(__IPYTHON__.input_hist)
        for L in F.readlines():
            L = L.rstrip()
            Llstrip = L.lstrip()
            raw_input('sage: %s'%L.rstrip())
            __IPYTHON__.input_hist_raw.append(L)
            if Llstrip[:5] == 'load ' or Llstrip[:7] == 'attach ' \
                   or Llstrip[:6] == 'iload ':
                log.offset -= 1
                L = do_prefilter_paste(L, False)
                if len(L.strip()) > 0:
                    _ip.runlines(L)
                L = ''
            else:
                L = preparser_ipython.preparse_ipython(L, not continuation)
            __IPYTHON__.input_hist.append(L)
            __IPYTHON__.push(L)
        log.offset += 1
        return ''


    #################################################################
    # A "load" command, like \r file in PARI or load "file" in MAGMA
    #################################################################
    if line[:5] == 'load ':
        # The -i so the file is run with the same environment,
        # e.g., including the "from sage import *"
        try:
            name = str(eval(line[5:])).strip()
        except:
            name = str(line[5:].strip())
        if name.lower().startswith('http://'):
            name = remote_file.get_remote_file(name)
        if isinstance(name, str):
            if not os.path.exists(name):
                raise ImportError, "File '%s' not found (be sure to give .sage, .py, or .pyx extension)"%name
            elif name.endswith('.py'):
                try:
                    line = '%run -i "' + name + '"'
                except IOError, s:
                    print s
                    raise ImportError, "Error loading '%s'"%name
            elif name.endswith('.sage'):
                try:
                    line = '%run -i "' + preparse_file_named(name) + '"'
                except IOError, s:
                    print s
                    raise ImportError, "Error loading '%s'"%name
                    line = ""
            elif name.endswith('.spyx') or name.endswith('.pyx'):
                line = load_cython(name)
            else:
                line = 'load("%s")'%name
                #raise ImportError, "Loading of '%s' not implemented (load .py, .pyx, and .sage files)"%name
                #line = ''

    # This is an attach command like in MAGMA.  The file to attach is
    # any file that could be loaded.  At attach time it is loaded as
    # above.  It is put in a list that is a variable with scope this
    # interpreter.py file.  Each time prefilter_paste is called and
    # continuation is False, check all attached file names and if any
    # have changed, compile the file, then use %run to load them again
    # as above.  This is like the MAGMA attach, but in some ways
    # better.  It is very nice for interactive code development.

    if line.startswith('attach '):
        # The -i so the file is run with the same environment,
        # e.g., including the "from sage import *"
        try:
            name = str(eval(line[7:]))
        except:
            name = str(line[7:].strip())
        name = os.path.abspath(name)
        if not os.path.exists(name):
            raise ImportError, "File '%s' not found  (be sure to give .sage, .py, or .pyx extension)."%name
        elif name.endswith('.py'):
            try:
                line = '%run -i "' + name + '"'
                attached[name] = os.path.getmtime(name)
            except IOError, OSError:
                raise ImportError, "File '%s' not found."%name
        elif name.endswith('.sage'):
            try:
                line = '%run -i "' + preparse_file_named(name) + '"'
                attached[name] = os.path.getmtime(name)
            except IOError, OSError:
                raise ImportError, "File '%s' not found."%name
        elif name.endswith('.pyx') or name.endswith('.spyx'):
            try:
                line = load_cython(name)
                attached[name] = os.path.getmtime(name)
            except IOError, OSError:
                raise ImportError, "File '%s' not found."%name
                line = ''
        else:
            #line = 'load("%s")'%name
            raise ImportError, "Attaching of '%s' not implemented (load .py, .pyx, and .sage files)"%name

    if len(line) > 0:
        line = preparser_ipython.preparse_ipython(line, not continuation)
    return line

def load_cython(name):
    cur = os.path.abspath(os.curdir)
    try:
        mod, dir  = cython.cython(name, compile_message=True, use_cache=True)
    except (IOError, OSError, RuntimeError), msg:
        print "Error compiling cython file:\n%s"%msg
        return ''
    import sys
    sys.path.append(dir)
    return 'from %s import *'%mod

def preparse_file_named_to_stream(name, out):
    r"""
    Preparse file named \code{name} (presumably a .sage file), outputting to
    stream \code{out}.
    """
    name = os.path.abspath(name)
    dir, _ = os.path.split(name)
    cur = os.path.abspath(os.curdir)
    os.chdir(dir)
    contents = open(name).read()
    parsed = preparse_file(contents, attached, do_time=True)
    os.chdir(cur)
    out.write("# -*- encoding: utf-8 -*-\n")
    out.write('#'*70+'\n')
    out.write('# This file was *autogenerated* from the file %s.\n' % name)
    out.write('#'*70+'\n')
    out.write(parsed)

def preparse_file_named(name):
    r"""
    Preparse file named \code{name} (presumably a .sage file), outputting to a
    temporary file.  Returns name of temporary file.
    """
    import sage.misc.misc
    name = os.path.abspath(name)
    tmpfilename = os.path.abspath(sage.misc.misc.tmp_filename(name) + ".py")
    out = open(tmpfilename,'w')
    preparse_file_named_to_stream(name, out)
    out.close()
    return tmpfilename

def sage_prefilter(self, block, continuation):
    """
    SAGE's prefilter for input.  Given a string block (usually a
    line), return the preparsed version of it.

    INPUT:
        block -- string (usually a single line, but not always)
        continuation -- whether or not this line is a continuation.
    """
    try:
        block2 = ''
        first = True
        B = block.split('\n')
        for i in range(len(B)):
            L = B[i]
            M = do_prefilter_paste(L, continuation or (not first))
            first = False
            # The L[:len(L)-len(L.lstrip())]  business here preserves
            # the whitespace at the beginning of L.
            if block2 != '':
                block2 += '\n'
            lstrip = L.lstrip()
            if lstrip[:5] == 'sage:' or lstrip[:3] == '>>>' or i==0:
                block2 += M
            else:
                block2 += L[:len(L)-len(lstrip)] + M

    except None:

        print "WARNING: An error occured in the SAGE parser while"
        print "parsing the following block:"
        print block
        print "Please report this as a bug (include the output of typing '%hist')."
        block2 = block

    return InteractiveShell._prefilter(self, block2, continuation)


import sage.server.support
def embedded():
    return sage.server.support.EMBEDDED_MODE

ipython_prefilter = InteractiveShell.prefilter
do_preparse=True
def preparser(on=True):
    """
    Turn on or off the SAGE preparser.

    NOTE: This only works on the command line.  To turn off preparsing
    in the notebook, switch to python mode.

    INPUT:
        on -- bool (default: True) if True turn on preparsing; if False, turn it off.

    EXAMPLES:
        sage: 2/3
        2/3
        sage: preparser(False)
        sage: 2/3
        0
        sage: preparser(True)
        sage: 2^3
        8
    """
    global do_preparse
    if on:
        do_preparse = True
        InteractiveShell.prefilter = sage_prefilter
    else:
        do_preparse = False
        InteractiveShell.prefilter = ipython_prefilter


import sagedoc
import IPython.OInspect
IPython.OInspect.getdoc = sagedoc.my_getdoc
IPython.OInspect.getsource = sagedoc.my_getsource


import __builtin__
_prompt = 'sage'

def set_sage_prompt(s):
    global _prompt
    _prompt = str(s)

def sage_prompt():
    log.update()
    return '%s'%_prompt

__builtin__.sage_prompt = sage_prompt


