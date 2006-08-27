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
          <type 'rational.Rational'>
    sage: a = 49928420832092
    sage: type(a)
          <type 'integer.Integer'>
    sage: a.factor()
          [(2, 2), (11, 1), (1134736837093, 1)]
    sage: v = [1,2,3]
    sage: type(v[0])
          <type 'integer.Integer'>

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
#       Copyright (C) 2004 William Stein <wstein@ucsd.edu>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#
#                  http://www.gnu.org/licenses/
#*****************************************************************************

__author__ = 'William Stein <wstein@ucsd.edu>'
__version__= '1.3.4'
__license__ = 'GPL'
__date__   = 'July 2006'

import os
import log

from IPython.iplib import InteractiveShell

import preparser_ipython
from preparser import preparse_file

import pyrex

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
        ipmagic(X)
    if os.path.exists('attach.sage'):
        X = do_prefilter_paste('attach "attach.sage"',False)
        ipmagic(X)


def do_prefilter_paste(line, continuation):
    """
    Alternate prefilter for input.
    """
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

    if not line[:7] == 'attach ' and not line[:5] == 'load ':
        for F in attached.keys():
            tm = attached[F]
            if os.path.exists(F) and os.path.getmtime(F) > tm:
                # Reload F.
                try:
                    if F[-3:] == '.py':
                        ipmagic('run -i "%s"'%F)
                    elif F[-5:] == '.sage':
                        ipmagic('run -i "%s"'%process_file(F))
                    else:
                        X = load_pyrex(F)
                        __IPYTHON__.push(X)

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
            name = line[6:].strip()
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
                    ipmagic(L)
                L = ''
            else:
                L = preparser_ipython.preparse_ipython(L)
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
            name = eval(line[5:])
        except:
            name = line[5:].strip()
        if isinstance(name, str):
            if name[-3:] == '.py':
                try:
                    line = '%run -i "' + name + '"'
                except IOError, s:
                    print s
                    raise ImportError, "Error loading '%s'"%name
            elif name[-5:] == '.sage':
                try:
                    line = '%run -i "' + process_file(name) + '"'
                except IOError, s:
                    print s
                    raise ImportError, "Error loading '%s'"%name
                    line = ""
            elif name[-5:] == '.spyx':
                line = load_pyrex(name)
            else:
                try:
                    line = '%run -i "' + name + '"'
                except IOError, s:
                    print s
                    raise ImportError, "Maybe file '%s' does not exist"%name
                    line = ""


    elif line[:13] == 'save_session(':
        line = 'save_session(locals(), %s'%line[13:]

    elif line[:14] == 'save_session (':
        line = 'save_session (locals(), %s'%line[14:]

    elif line[:13] == 'load_session(':
        line = 'load_session(locals(), %s'%line[13:]

    elif line[:14] == 'load_session (':
        line = 'load_session (locals(), %s'%line[14:]

    elif ''.join(line.split()) == 'show_identifiers()':
        line = 'show_identifiers(locals())'

    # This is an attach command like in MAGMA.  The file to attach is
    # any file that could be loaded.  At attach time it is loaded as
    # above.  It is put in a list that is a variable with scope this
    # interpreter.py file.  Each time prefilter_paste is called and
    # continuation is False, check all attached file names and if any
    # have changed, compile the file, then use %run to load them again
    # as above.  This is like the MAGMA attach, but in some ways
    # better.  It is very nice for interactive code development.

    if line[:7] == 'attach ':
        # The -i so the file is run with the same environment,
        # e.g., including the "from sage import *"
        try:
            name = eval(line[7:])
        except:
            name = line[7:].strip()
        name = os.path.abspath(name)
        if not os.path.exists(name):
            raise ImportError, "File '%s' not found."%name
        elif name[-3:] == '.py':
            try:
                line = '%run -i "' + name + '"'
                attached[name] = os.path.getmtime(name)
            except IOError, OSError:
                raise ImportError, "File '%s' not found."%name
        elif name[-5:] == '.sage':
            try:
                line = '%run -i "' + process_file(name) + '"'
                attached[name] = os.path.getmtime(name)
            except IOError, OSError:
                raise ImportError, "File '%s' not found."%name
        elif name[-5:] == '.spyx':
            try:
                line = load_pyrex(name)
                attached[name] = os.path.getmtime(name)
            except IOError, OSError:
                raise ImportError, "File '%s' not found."%name
                line = ''
        else:
            raise ImportError, "load: file (=%s) must have extension .py, .sage or .spyx"%name
    if len(line) > 0:
        line = preparser_ipython.preparse_ipython(line)
    return line

def load_pyrex(name):
    cur = os.path.abspath(os.curdir)
    try:
        mod, dir  = pyrex.pyrex(name, compile_message=True)
    except (IOError, OSError, RuntimeError), msg:
        print "Error compiling pyrex file:\n%s"%msg
        return ''
    import sys
    sys.path.append(dir)
    return 'from %s import *'%mod

def process_file(name):
    cur = os.path.abspath(os.curdir)
    dir, name = os.path.split(os.path.abspath(name))
    name2 = "%s/%s.py"%(dir,name[:name.find('.')])
    os.chdir(dir)
    contents = open(name).read()
    parsed = preparse_file(contents, attached)
    os.chdir(cur)
    W = open(name2,'w')
    W.write('#'*70+'\n')
    W.write('# This file was *autogenerated* from the file %s.\n'%name)
    W.write('#'*70+'\n')
    W.write(parsed)
    W.close()
    return name2


def sage_prefilter(self, line, continuation):
    """
    Alternate prefilter for input.
    """
    try:
        line2 = do_prefilter_paste(line, continuation)

    except None:

        print "WARNING: An error occured in the SAGE parser while"
        print "parsing the following line:"
        print line
        print "Please report this as a bug (include the output of typing '%hist')."
        line2 = line

    from IPython.iplib import InteractiveShell
    return InteractiveShell._prefilter(self, line2, continuation)



# Rebind this to be the new IPython prefilter:
InteractiveShell.prefilter = sage_prefilter

# Clean up the namespace.
#del InteractiveShell, sage_prefilter

#######################################
## Add detex'ing of documentation
#######################################
import sagedoc
import inspect

def my_getdoc(obj):
    try:
        ds = obj._sage_doc_()
    except (AttributeError, TypeError):  # TypeError for interfaces
        try:
            ds = inspect.getdoc(obj)
        except:
            return None
    if ds is None:
        return None
    return sagedoc.format(ds)

import IPython.OInspect
IPython.OInspect.getdoc = my_getdoc



import log

import __builtin__
_prompt = 'sage'

def set_sage_prompt(s):
    global _prompt
    _prompt = str(s)

def sage_prompt():
    log.update()
    return '%s'%_prompt

__builtin__.sage_prompt = sage_prompt


