"""
SAGE pre-parser.

AUTHOR:
    -- William Stein (2006-02-19): fixed bug when loading .py files.
    -- William Stein (2006-03-09): * fixed crash in parsing exponentials
                                   * precision of real literals now determined
                                     by digits of input (like mathematica).
"""
#EXAMPLES:
#These examples all illustrate input lines whose pre-parsing is subtle.
# (todo)



###########################################################################
#       Copyright (C) 2006 William Stein <wstein@ucsd.edu>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#                  http://www.gnu.org/licenses/
###########################################################################

import os

def isalphadigit_(s):
    return s.isalpha() or s.isdigit() or s=="_"

def last_bracket_is_after_identifier(line, i):
    """
    Return True if and only if the previous bracket before line[i]
    is after a valid identifier followed possibly by some space.
    """
    line = line[:i]
    j = line.rfind("[")
    if j == -1:
        return False
    line = line[:j].rstrip()
    if len(line) == 0:
        return False
    c = line[len(line)-1]
    return c == "]" or c == ")" or isalphadigit_(c)


in_single_quote = False
in_double_quote = False
in_triple_quote = False
def in_quote():
    return in_single_quote or in_double_quote or in_triple_quote

def preparse(line, reset=True):
    global in_single_quote, in_double_quote, in_triple_quote
    line = line.split("\n")[0]   # xreadlines leaves the '\n' at end of line
    L = line.lstrip()
    if len(L) > 0 and L[0] in ['#', '!', '%']:
        return line

    # Wrap integers with ZZ() and reals with RR().
    def wrap_num(i, line, is_real, num_start):
        if is_real:
            # make real field with slightly less precision than
            # number of digits of our number
            n = int(3.32192*(i-num_start-1))
            # we may want to change to something like this:
            #  add a few digits (how many)
            #n = int(3.32192*(i-num_start-1)) + 5
            O = "RealField(max(%s,RR.precision()))('"%n; C = "')"
        else:
            O = "ZZ("; C = ")"
        line = line[:num_start] + O + line[num_start:i]  \
               + C + line[i:]
        return line, len(O+C)

    i = 0
    num_start = -1
    bracket_depth = 0
    seen_comma = False
    in_number = False
    is_real = False
    if reset:
        in_single_quote = False
        in_double_quote = False
        in_triple_quote = False

    while i < len(line):
        if bracket_depth > 0 and line[i] == ",":
            seen_comma = True

        # Decide if we should wrap a particular integer or real literal
        if in_number:
            if line[i] == "." and not (i+1 < len(line) and line[i+1].isalpha()):
                is_real = True
            elif not line[i].isdigit():
                # end of a number
                # Do we wrap?
                if i < len(line) and line[i] in 'eE':
                    i += 1   # skip wrapping and parsing
                    if i < len(line) and line[i] == '-':
                        i += 2
                elif bracket_depth == 0 or (bracket_depth > 0 and \
                                          not last_bracket_is_after_identifier(line, i)):
                    line, n = wrap_num(i, line, is_real, num_start)
                    i += n
                in_number = False
                is_real = False
                continue

        # experimental support for generator construction
        # syntax:  "obj.<gen0,gen1,...,genN> = objConstructor(...)"
        # is converted into
        # "obj = objConstructor(...); \
        #  obj.assign_names(["gen0", "gen1", ..., "genN"]); \
        #  (gen0, gen1, ..., genN,) = obj.gens()"
        # LIMITATIONS:
        # - The entire constructor must be on one line.
        # - The line must contain no other statements.
        elif line[i] == "." and line[i+1] == "<" and not in_quote():
            try:
                gen_end = line.index(">", i+2)
            except ValueError:
                # Syntax Error -- let Python notice and raise the error
                i += 2
                continue
            # parse out the object name and the list of generator names
            gen_obj = line[:i].strip()
            gen_list = map(lambda s: s.strip(), line[i+2:gen_end].split(","))
            # format names as a list of strings and a list of variables
            gen_names = str(gen_list)
            gen_vars  = ", ".join(gen_list)
            # rewrite the input line as three commands
            line = "; ".join([line[:i] + line[gen_end+1:],
                              "%s.assign_names(%s)" % (gen_obj, gen_names),
                              "(%s,) = %s.gens()" % (gen_vars, gen_obj)])
            continue

        # exponents can be either ^ or **
        elif line[i] == "^" and not in_quote():
            line = line[:i] + "**" + line[i+1:]
            i += 2
            continue

        elif line[i] == "." and i > 0 and i < len(line)-1 and not in_quote() and \
                 (isalphadigit_(line[i-1]) or line[i-1] == ")" or line[i-1] == ']') and line[i+1].isdigit():
            # Generators like in MAGMA: replace all ".<number>" by ".gen(<number>)"
            # If . is preceeded by \, then replace "\." by ".".
            j = i+1
            while j < len(line) and line[j].isdigit():
                j += 1
            line = line[:i] + ".gen(" + line[i+1:j] + ")" + line[j:]
            i = j+4

        # Update bracket depth
        if line[i] == "[" and not in_quote():
            bracket_depth += 1

        elif line[i] == "]" and not in_quote():
            bracket_depth -= 1
            if bracket_depth == 0:
                seen_comma = False

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

        if line[i] == '"':
            if not in_quote():
                in_double_quote = True
                i += 1
                continue
            elif in_double_quote:
                in_double_quote = False
                i += 1
                continue

        if line[i:i+3] == '"""':
            if not in_quote():
                in_triple_quote = True
                i += 3
                continue
            elif in_triple_quote:
                in_triple_quote = False
                i += 3
                continue

        if not in_number and not in_quote() and line[i].isdigit() and \
               (i == 0 or (i > 0 and not isalphadigit_(line[i-1]))):
            in_number = True
            num_start = i
        i += 1

    if in_number:
        line, _ = wrap_num(i, line, is_real, num_start)

    # Time command like in MAGMA: (commented out, since it's standard in IPython now)
    L = line.lstrip()
    #if L[:5] == "time ":
    #    # strip semicolon from end of line
    #     if line[-1:] == ";":
    #         line = line[:-1]
    #     indent = ' '*(len(line) - len(L))
    #     line = indent + '__time__=misc.cputime(); __wall__=misc.walltime(); %s; print \
    #"Time: CPU %%.2f s, Wall: %%.2f s"%%(misc.cputime(__time__), misc.walltime(__wall__))'%L[4:]

    return line



######################################################
## Apply the preparser to an entire file
######################################################

def preparse_file(contents, attached={}):
    if not isinstance(contents, str):
        raise TypeError, "contents must be a string"

    assert isinstance(contents, str)
    # We keep track of which files have been loaded so far
    # in order to avoid a recursive load that would result
    # in creating a massive file and crashing.
    loaded_files = []

    F = []
    A = contents.split('\n')
    i = 0
    while i < len(A):
        L = A[i]
        if L[:7] == "attach ":
            name = eval(L[7:])
            try:
                if not attached.has_key(name):
                    t = os.path.getmtime(name)
                    attached[name] = t
            except IOError, OSError:
                pass
            L = 'load ' + L[7:]
        if L[:5] == "load ":
            try:
                name_load = eval(L[5:])
                if name_load in loaded_files:
                    i += 1
                    continue
                loaded_files.append(name_load)
            except:
                pass
            else:
                if isinstance(name_load, str):
                    if name_load[-3:] == '.py':
                        ipmagic('run -i "%s"'%name_load)
                        L = ''
                    elif name_load[-5:] == '.sage':
                        try:
                            G = open(name_load)
                        except IOError:
                            print "File '%s' not found, so skipping re-load"%name
                            i += 1
                            continue
                        else:
                            A = A[:i] + G.readlines() + A[i+1:]
                            continue
                    else:
                        import interpreter
                        L = interpreter.load_pyrex(name_load)
        F.append(preparse(L, reset=False))
        i += 1
    # end while

    return '\n'.join(F)

