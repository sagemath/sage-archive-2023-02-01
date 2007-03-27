"""
SAGE pre-parser.

AUTHOR:
    -- William Stein (2006-02-19): fixed bug when loading .py files.
    -- William Stein (2006-03-09): * fixed crash in parsing exponentials
                                   * precision of real literals now determined
                                     by digits of input (like mathematica).
    -- Joe Wetherell (2006-04-14): * added MAGMA-style constructor preparsing.

EXAMPLES:

PREPARSE:
    sage: preparse('2/3')
    'Integer(2)/Integer(3)'
    sage: preparse('2.5')
    "RealNumber('2.5')"
    sage: preparse('2^3')
    'Integer(2)**Integer(3)'
    sage: preparse('a^b')            # exponent
    'a**b'
    sage: preparse('a**b')
    'a**b'
    sage: preparse('G.0')            # generator
    'G.gen(0)'
    sage: preparse('a = 939393R')    # raw
    'a = 939393'

In SAGE methods can also be called on integer and real literals (note
that in pure Python this would be a syntax error).
    sage: 16.sqrt()
    4
    sage: 87.factor()
    3 * 29
    sage: 15.10.sqrt()
    3.88587184554509
    sage: preparse('87.sqrt()')
    'Integer(87).sqrt()'
    sage: preparse('15.10.sqrt()')
    "RealNumber('15.10').sqrt()"

Note that calling methods on int literals in pure Python is a
syntax error, but SAGE allows this for SAGE integers and reals,
because users frequently request it.

    sage: eval('4.__add__(3)')
    Traceback (most recent call last):
    ...
    SyntaxError: invalid syntax


RAW LITERALS:

Raw literals are not preparsed, which can be useful from an efficiency
point of view.  Just like Python ints are denoted by an L, in SAGE raw
integer and floating literals are followed by an"r" (or "R") for raw,
meaning not preparsed.

We create a raw integer.
    sage: a = 393939r
    sage: a
    393939
    sage: type(a)
    <type 'int'>

We create a raw float:
    sage: z = 1.5949r
    sage: z
    1.5949
    sage: type(z)
    <type 'float'>

You can also use an upper case letter:
    sage: z = 3.1415R
    sage: z
    3.1415000000000002
    sage: type(z)
    <type 'float'>

This next example illustrates how raw literals can be very useful in
certain cases.  We make a list of even integers up to 10000.
    sage: v = [ 2*i for i in range(10000)]

This talkes a noticeable fraction of a second (e.g., 0.25 seconds).
After preparsing, what Python is really executing is the following:

    sage: preparse('v = [ 2*i for i in range(10000)]')
    'v = [ Integer(2)*i for i in range(Integer(10000))]'

If instead we use a raw 2 we get executation that is ``instant'' (0.00 seconds):
    sage: v = [ 2r * i for i in range(10000r)]

Behind the scenes what happens is the following:
    sage: preparse('v = [ 2r * i for i in range(10000r)]')
    'v = [ 2 * i for i in range(10000)]'

WARNING: The result of the above two expressions is different.  The
first one computes a list of SAGE integers, whereas the second creates
a list of Python integers.  Python integers are typically much more
efficient than SAGE integers when they are very small; large SAGE
integers are much more efficient than Python integers, since they are
implemented using the GMP C library.



"""



###########################################################################
#       Copyright (C) 2006 William Stein <wstein@gmail.com>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#                  http://www.gnu.org/licenses/
###########################################################################

import os

def isalphadigit_(s):
    return s.isalpha() or s.isdigit() or s=="_"


in_single_quote = False
in_double_quote = False
in_triple_quote = False
def in_quote():
    return in_single_quote or in_double_quote or in_triple_quote

def preparse(line, reset=True, do_time=False, ignore_prompts=False):
    global in_single_quote, in_double_quote, in_triple_quote
    line = line.rstrip()  # xreadlines leaves the '\n' at end of line
    L = line.lstrip()
    if len(L) > 0 and L[0] in ['#', '!']:
        return line

    if L.startswith('...'):
        i = line.find('...')
        return line[:i+3] + preparse(line[i+3:], reset=reset, do_time=do_time, ignore_prompts=ignore_prompts)

    # Wrap integers with ZZ() and reals with RR().
    def wrap_num(i, line, is_real, num_start):
        zz = line[num_start:i]
        if is_real or '.' in zz:
            if zz[-1] == '.' and i < len(line) and line[i].isalpha():
                # by popular demand -- this allows, e.g., 173.sqrt().
                if '.' in zz[:-1]:
                    O = "RealNumber('"; C="')."
                else:
                    O = "Integer("; C = ")."
                zz = zz[:-1]
            else:
                O = "RealNumber('"; C="')"
        else:
            O = "Integer("; C = ")"
        line = line[:num_start] + O + zz + C + line[i:]
        return line, len(O+C)

    i = 0
    num_start = -1
    in_number = False
    is_real = False
    if reset:
        in_single_quote = False
        in_double_quote = False
        in_triple_quote = False


    if ignore_prompts:
        # Get rid of leading sage: so that pasting of examples from
        # the documentation works.
        for prompt in ['sage:', '>>>']:
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

    while i < len(line):
        # Update quote parsing
        if line[i] == "'":
            if not in_quote():
                in_single_quote = True
                i += 1
                continue
            elif in_single_quote:
                in_single_quote = False
                i += 1
                continue
        elif line[i:i+3] == '"""':
            if not in_quote():
                in_triple_quote = True
                i += 3
                continue
            elif in_triple_quote:
                in_triple_quote = False
                i += 3
                continue
        elif line[i] == '"':
            if not in_quote():
                in_double_quote = True
                i += 1
                continue
            elif in_double_quote:
                in_double_quote = False
                i += 1
                continue


        # Decide if we should wrap a particular integer or real literal
        if in_number:
            if line[i] == ".":
                is_real = True
            elif not line[i].isdigit():
                # end of a number
                # Do we wrap?
                if in_quote():
                    # do not wrap
                    pass
                elif i < len(line) and line[i] in 'eE':
                    # Yes, in scientific notation, so will wrap later
                    is_real = True
                    i += 1
                    if i < len(line) and line[i] == '-':
                        i += 2
                    continue
                elif i < len(line) and line[i] in 'rR':
                    # Raw number so do not wrap; but have to get rid of the "r".
                    line = line[:i] + line[i+1:]
                else:
                    line, n = wrap_num(i, line, is_real, num_start)
                    i += n
                in_number = False
                is_real = False
                continue

        # Support for generator construction syntax:
        # "obj.<gen0,gen1,...,genN> = objConstructor(...)"
        # is converted into
        # "obj = objConstructor(..., names=("gen0", "gen1", ..., "genN")); \
        #  (gen0, gen1, ..., genN,) = obj.gens()"
        #
        # Also, obj.<gen0,gen1,...,genN> = R[...] is converted into
        # "obj = R['gen0,gen1,..., genN']; (gen0, gen1, ..., genN,) = obj.gens()"
        #
        # LIMITATIONS:
        #    - The entire constructor *must* be on one line.
        #
        # AUTHORS:
        #     -- 2006-04-14: Joe Wetherell (jlwether@alum.mit.edu)
        #     -- 2006-04-17: William Stein - improvements to allow multiple statements.
        #     -- 2006-05-01: William -- fix bug that Joe found
        #     -- 2006-10-31: William -- fix so obj doesn't have to be mutated
        elif (line[i:i+2] == ".<") and not in_quote():
            try:
                gen_end = line.index(">", i+2)
            except ValueError:
                # Syntax Error -- let Python notice and raise the error
                i += 2
                continue

            gen_begin = i
            while gen_begin > 0 and line[gen_begin-1] != ';':
                gen_begin -= 1

            # parse out the object name and the list of generator names
            gen_obj = line[gen_begin:i].strip()
            gen_list = [s.strip() for s in line[i+2:gen_end].split(',')]
            for g in gen_list:
                if not g.isalnum() or len(g) == 0 or not g[0].isalpha():
                    raise SyntaxError, "variable name (='%s') must be alpha-numeric and begin with a letter"%g

            # format names as a list of strings and a list of variables
            gen_names = tuple(gen_list)
            gen_vars  = ", ".join(gen_list)

            # find end of constructor:
            #    either end of line, next semicolon, or next #.
            line_after = line[gen_end:]
            c = line_after.find('#')
            if c==-1: c = len(line_after)
            s = line_after.find(';')
            if s==-1: s = len(line_after)
            c = min(c,s) + gen_end

            # Find where the parenthesis of the constructor ends
            if line[:c].rstrip()[-1] == ']':
                # brackets constructor
                c0 = line[:c].find(']')
                d0 = line[:c0].rfind('[')
                if c0 == -1:
                    raise SyntaxError, 'constructor must end with ) or ]'
                line_new = '%s"%s"%s; (%s,) = %s.gens()'%(
                    line[:i] + line[gen_end+1:d0+1], gen_vars,
                    line[c0:c], gen_vars, gen_obj)
            else:
                c0 = line[:c].rfind(')')
                # General constructor -- rewrite the input line as two commands
                # We have to determine whether or not to put a comma before
                # the list of names.  We do this only if there are already
                # arguments to the constructor.  Some constructors have no
                # arguments, e.g., "K.<a> = f.root_field(  )"
                c1 = line[:c0].rfind('(')
                if len(line[c1+1:c0].strip()) > 0:
                    sep = ','
                else:
                    sep = ''

                line_new = '%s%snames=%s); (%s,) = %s.gens()'%(
                    line[:i] + line[gen_end+1:c0], sep, gen_names,
                    gen_vars, gen_obj)

            line = line_new + line[c:]
            #i = len(line_new)
            i += 1

            continue

        # exponents can be either ^ or **
        elif line[i] == "^" and not in_quote():
            line = line[:i] + "**" + line[i+1:]
            i += 2
            continue

        elif line[i] == "." and i > 0 and i < len(line)-1 and not in_quote() and \
                 (isalphadigit_(line[i-1]) or line[i-1] == ")" or line[i-1] == ']') and line[i+1].isdigit():
            # Generators: replace all ".<number>" by ".gen(<number>)"
            # If . is preceeded by \, then replace "\." by ".".
            j = i+1
            while j < len(line) and line[j].isdigit():
                j += 1
            line = line[:i] + ".gen(" + line[i+1:j] + ")" + line[j:]
            i = j+4

        if     not in_number and \
               not in_quote()and \
               (line[i].isdigit() or \
                   (len(line)>i+1 and line[i] == '.' and line[i+1].isdigit())) and \
               (i == 0 or (i > 0 and not (isalphadigit_(line[i-1]) \
                                          or line[i-1] == ')'))):
            in_number = True
            num_start = i

        # Decide if we hit a comment, so we're done.
        if line[i] == '#' and not (in_single_quote or in_double_quote or in_triple_quote):
            i = len(line)
            break

        i += 1

    if in_number:
        line, _ = wrap_num(i, line, is_real, num_start)

    # Time command like in MAGMA: (commented out, since it's standard in IPython now)
    L = line.lstrip()

    if do_time:
        if L[:5] == "time ":
            # strip semicolon from end of line
             if line[-1:] == ";":
                 line = line[:-1]
             indent = ' '*(len(line) - len(L))
             line = indent + '__time__=misc.cputime(); __wall__=misc.walltime(); %s; print \
        "Time: CPU %%.2f s, Wall: %%.2f s"%%(misc.cputime(__time__), misc.walltime(__wall__))'%L[4:]

    return line



######################################################
## Apply the preparser to an entire file
######################################################

def preparse_file(contents, attached={}, magic=True,
                  do_time=False, ignore_prompts=False):
    if not isinstance(contents, str):
        raise TypeError, "contents must be a string"

    assert isinstance(contents, str)
    # We keep track of which files have been loaded so far
    # in order to avoid a recursive load that would result
    # in creating a massive file and crashing.
    loaded_files = []

    F = []
    A = contents.splitlines()
    i = 0
    while i < len(A):
        L = A[i].rstrip()
        if magic and L[:7] == "attach ":
            name = os.path.abspath(eval(L[7:]))
            try:
                if not attached.has_key(name):
                    t = os.path.getmtime(name)
                    attached[name] = t
            except IOError, OSError:
                pass
            L = 'load ' + L[7:]

        if magic and L[:5] == "load ":
            try:
                name_load = str(eval(L[5:]))
            except:
                name_load = L[5:].strip()
            if name_load in loaded_files:
                i += 1
                continue
            loaded_files.append(name_load)
            if name_load[-3:] == '.py':
                _ip.magic('run -i "%s"'%name_load)
                L = ''
            elif name_load[-5:] == '.sage':
                try:
                    G = open(name_load)
                except IOError:
                    print "File '%s' not found, so skipping load of %s"%(name_load, name_load)
                    i += 1
                    continue
                else:
                    A = A[:i] + G.readlines() + A[i+1:]
                    continue
            elif name_load[-5:] == '.spyx':
                import interpreter
                L = interpreter.load_sagex(name_load)
            else:
                print "Loading of '%s' not implemented (load .py, .spyx, and .sage files)"%name_load
                L = ''
                continue
        M = preparse(L, reset=(i==0), do_time=do_time, ignore_prompts=ignore_prompts)
        F.append(M)
        i += 1
    # end while

    return '\n'.join(F)

