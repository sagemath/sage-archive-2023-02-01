"""
SAGE pre-parser.

AUTHOR:
    -- William Stein (2006-02-19): fixed bug when loading .py files.
    -- William Stein (2006-03-09): * fixed crash in parsing exponentials
                                   * precision of real literals now determined
                                     by digits of input (like mathematica).
    -- Joe Wetherell (2006-04-14): * added MAGMA-style constructor preparsing.
    -- Bobby Moretti (2007-01-25): * added preliminary function assignment
                                     notation
    -- Robert Bradshaw (2007-09-19): * strip_string_literals, containing_block
                                       utility functions. Arrr!
                                     * Add [1,2,..,n] notation.
    -- Robert Bradshaw (2008-01-04): * Implicit multiplication

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
    sage: preparse('a b c in L')     # implicit multiplication
    'a*b*c in L'
    sage: preparse('2e3x + 3exp(y)')
    "RealNumber('2e3')*x + Integer(3)*exp(y)"

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

SYMBOLIC FUNCTIONAL NOTATION:
    sage: a=10; f(theta, beta) = theta + beta; b = x^2 + theta
    sage: f
    (theta, beta) |--> theta + beta
    sage: a
    10
    sage: b
    x^2 + theta
    sage: f(theta,theta)
    2*theta

    sage: a = 5; f(x,y) = x*y*sqrt(a)
    sage: f
    (x, y) |--> sqrt(5)*x*y

This involves an =-, but should still be turned into a symbolic expression:
    sage: preparse('a(x) =- 5')
    '_=var("x");a=symbolic_expression(- Integer(5)).function(x)'
    sage: f(x)=-x
    sage: f(10)
    -10

This involves -=, which should not be turned into a symbolic
expression (of course a(x) isn't an identifier, so this will never be
valid):
    sage: preparse('a(x) -= 5')
    'a(x) -= Integer(5)'

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
import os, re
import pdb

def isalphadigit_(s):
    return s.isalpha() or s.isdigit() or s=="_"

keywords = """
and       del       from      not       while
as        elif      global    or        with
assert    else      if        pass      yield
break     except    import    print
class     exec      in        raise
continue  finally   is        return
def       for       lambda    try
""".split()

in_single_quote = False
in_double_quote = False
in_triple_quote = False

def in_quote():
    return in_single_quote or in_double_quote or in_triple_quote


def strip_string_literals(code):
    r"""
    Returns a string with all literal quotes replaced with
    labels and a dict of labels for re-subsitution. This makes
    parsing much easier.

    EXAMPLES:
        sage: from sage.misc.preparser import strip_string_literals
        sage: s, literals = strip_string_literals(r'''['a', "b", 'c', "d\""]''')
        sage: s
        '[%(L1)s, %(L2)s, %(L3)s, %(L4)s]'
        sage: literals
        {'L4': '"d\\""', 'L2': '"b"', 'L3': "'c'", 'L1': "'a'"}
        sage: print s % literals
        ['a', "b", 'c', "d\""]
        sage: print strip_string_literals(r'-"\\\""-"\\"-')[0]
        -%(L1)s-%(L2)s-

    Triple-quotes are handled as well.
        sage: s, literals = strip_string_literals("[a, '''b''', c, '']")
        sage: s
        '[a, %(L1)s, c, %(L2)s]'
        sage: print s % literals
        [a, '''b''', c, '']

    Comments are subsituted too:
        sage: s, literals = strip_string_literals("code '#' # ccc 't'"); s
        'code %(L1)s #%(L2)s'
        sage: s % literals
        "code '#' # ccc 't'"
    """
    new_code = []
    literals = {}
    counter = 0
    start = q = 0
    in_quote = False
    raw = False
    while True:
        sig_q = code.find("'", q)
        dbl_q = code.find('"', q)
        hash_q = code.find('#', q)
        q = min(sig_q, dbl_q)
        if q == -1: q = max(sig_q, dbl_q)
        if not in_quote and hash_q != -1 and (q == -1 or hash_q < q):
            # it's a comment
            newline = code.find('\n', hash_q)
            if newline == -1: newline = len(code)
            counter += 1
            label = "L%s" % counter
            literals[label] = code[hash_q+1:newline]
            new_code.append(code[start:hash_q].replace('%','%%'))
            new_code.append("#%%(%s)s" % label)
            start = q = newline
        elif q == -1:
            new_code.append(code[start:].replace('%','%%'))
            break
        elif in_quote:
            if not raw and code[q-1] == '\\':
                k = 2
                while code[q-k] == '\\':
                    k += 1
                if k % 2 == 0:
                    q += 1
            if code[q:q+len(in_quote)] == in_quote:
                counter += 1
                label = "L%s" % counter
                literals[label] = code[start:q+len(in_quote)]
                new_code.append("%%(%s)s" % label)
                q += len(in_quote)
                start = q
                in_quote = False
            else:
                q += 1
        else:
            raw = q>0 and code[q-1] in ['r', 'R']
            if len(code) >= q+3 and (code[q+1] == code[q] == code[q+2]):
                in_quote = code[q]*3
            else:
                in_quote = code[q]
            new_code.append(code[start:q].replace('%', '%%'))
            start = q
            q += len(in_quote)

    return "".join(new_code), literals


def containing_block(code, ix, delimiters=['()','[]','{}'], require_delim=True):
    """
    Returns the smallest range (start,end) such that code[start,end]
    is delimited by balanced parentheses/brackets/braces.

    EXAMPLES:
        sage: from sage.misc.preparser import containing_block
        sage: s = "factor(next_prime(L[5]+1))"
        sage: s[22]
        '+'
        sage: start, end = containing_block(s, 22); print start, end
        17 25
        sage: s[start:end]
        '(L[5]+1)'
        sage: s[20]
        '5'
        sage: start, end = containing_block(s, 20); s[start:end]
        '[5]'
        sage: start, end = containing_block(s, 20, delimiters=['()']); s[start:end]
        '(L[5]+1)'
        sage: start, end = containing_block(s, 10); s[start:end]
        '(next_prime(L[5]+1))'
    """
    openings = "".join([d[0] for d in delimiters])
    closings = "".join([d[-1] for d in delimiters])
    levels = [0] * len(openings)
    p = 0
    start = ix
    while start >= 0:
        start -= 1
        if start == -1:
            if require_delim:
                raise SyntaxError, "Unbalanced or missing ()'s"
            else:
                break
        if code[start] in openings:
            p = openings.index(code[start])
            levels[p] -= 1
            if levels[p] == -1:
                break
        elif code[start] in closings:
            p = closings.index(code[start])
            levels[p] += 1
    if start == -1:
        return 0, len(code)
    end = ix
    level = 0
    while end < len(code):
        end += 1
        if end == len(code):
            raise SyntaxError, "Unbalanced or missing ()'s"
        if code[end] == openings[p]:
            level += 1
        elif code[end] == closings[p]:
            level -= 1
            if level == -1:
                break
    return start, end+1


def parse_ellipsis(code, preparse_step=True):
    """
    Preparse [0,2,..,n] notation.

    EXAMPLES:
        sage: from sage.misc.preparser import parse_ellipsis
        sage: parse_ellipsis("[1,2,..,n]")
        '(ellipsis_range(1,2,Ellipsis,n))'
        sage: parse_ellipsis("for i in (f(x) .. L[10]):")
        'for i in (ellipsis_iter(f(x) ,Ellipsis, L[10])):'
    """
    ix = code.find('..')
    while ix != -1:
        if ix == 0:
            raise SyntaxError, "Cannot start line with ellipsis."
        elif code[ix-1]=='.':
            # '...' be valid Python in index slices
            code = code[:ix-1] + "Ellipsis" + code[ix+2:]
        elif len(code) >= ix+3 and code[ix+2]=='.':
            # '...' be valid Python in index slices
            code = code[:ix] + "Ellipsis" + code[ix+3:]
        else:
            start_list, end_list = containing_block(code, ix, ['()','[]'])
            arguments = code[start_list+1:end_list-1].replace('...', ',Ellipsis,').replace('..', ',Ellipsis,')
            arguments = re.sub(r',\s*,', ',', arguments)
            if preparse_step:
                arguments = arguments.replace(';', ', step=')
            range_or_iter = 'range' if code[start_list]=='[' else 'iter'
            code = "%s(ellipsis_%s(%s))%s" %  (code[:start_list],
                                               range_or_iter,
                                               arguments,
                                               code[end_list:])
        ix = code.find('..')
    return code

def strip_prompts(line):
    r"""Get rid of leading sage: and >>> prompts so that pasting of examples from
    the documentation works.

    sage: from sage.misc.preparser import strip_prompts
    sage: strip_prompts("sage: 2 + 2")
    '2 + 2'
    sage: strip_prompts(">>>   3 + 2")
    '3 + 2'
    sage: strip_prompts("  2 + 4")
    '  2 + 4'
    """
    for prompt in ['sage:', '>>>']:
        if line.startswith(prompt):
            line = line[len(prompt):].lstrip()
            break
    return line

def parse_generators(line, start_index):
    r"""Parse R.<a, b> in line[start_index:].

    Returns (modified_line, continue_index); you should resume parsing
    modified_line[continue_index:].

    Support for generator construction syntax:
    "obj.<gen0,gen1,...,genN> = objConstructor(...)"
    is converted into
    "obj = objConstructor(..., names=("gen0", "gen1", ..., "genN")); \
     (gen0, gen1, ..., genN,) = obj.gens()"

    Also, obj.<gen0,gen1,...,genN> = R[interior] is converted into
    "obj = R[interior]; (gen0, gen1, ..., genN,) = obj.gens()"

    LIMITATIONS:
       - The entire constructor *must* be on one line.

    AUTHORS:
        -- 2006-04-14: Joe Wetherell (jlwether@alum.mit.edu)
        -- 2006-04-17: William Stein - improvements to allow multiple statements.
        -- 2006-05-01: William -- fix bug that Joe found
        -- 2006-10-31: William -- fix so obj doesn't have to be mutated

    TESTS:
        sage: from sage.misc.preparser import preparse

        Vanilla:

        sage: preparse("R.<x> = ZZ['x']")
        "R = ZZ['x']; (x,) = R._first_ngens(Integer(1))"
        sage: preparse("R.<x,y> = ZZ['x,y']")
        "R = ZZ['x,y']; (x, y,) = R._first_ngens(Integer(2))"

        No square brackets:

        sage: preparse("R.<x> = PolynomialRing(ZZ, 'x')")
        "R = PolynomialRing(ZZ, 'x',names=('x',)); (x,) = R._first_ngens(Integer(1))"
        sage: preparse("R.<x,y> = PolynomialRing(ZZ, 'x,y')")
        "R = PolynomialRing(ZZ, 'x,y',names=('x', 'y')); (x, y,) = R._first_ngens(Integer(2))"

        Names filled in:

        sage: preparse("R.<x> = ZZ[]")
        "R = ZZ['x']; (x,) = R._first_ngens(Integer(1))"
        sage: preparse("R.<x,y> = ZZ[]")
        "R = ZZ['x, y']; (x, y,) = R._first_ngens(Integer(2))"

        Names given not the same as generator names:

        sage: preparse("R.<x> = ZZ['y']")
        "R = ZZ['y']; (x,) = R._first_ngens(Integer(1))"
        sage: preparse("R.<x,y> = ZZ['u,v']")
        "R = ZZ['u,v']; (x, y,) = R._first_ngens(Integer(2))"

        Number fields:

        sage: preparse("K.<a> = QQ[2^(1/3)]")
        'K = QQ[Integer(2)**(Integer(1)/Integer(3))]; (a,) = K._first_ngens(Integer(1))'
        sage: preparse("K.<a, b> = QQ[2^(1/3), 2^(1/2)]")
        'K = QQ[Integer(2)**(Integer(1)/Integer(3)), Integer(2)**(Integer(1)/Integer(2))]; (a, b,) = K._first_ngens(Integer(2))'

        Just the .<> notation:

        sage: preparse("R.<x> = ZZx")
        'R = ZZx; (x,) = R._first_ngens(Integer(1))'
        sage: preparse("R.<x, y> = a+b")
        'R = a+b; (x, y,) = R._first_ngens(Integer(2))'

        Ensure we don't eat too much:

        sage: preparse("R.<x, y> = ZZ;2")
        'R = ZZ; (x, y,) = R._first_ngens(Integer(2));Integer(2)'
        sage: preparse("R.<x, y> = ZZ['x,y'];2")
        "R = ZZ['x,y']; (x, y,) = R._first_ngens(Integer(2));Integer(2)"
    """
    i = start_index
    if not line.startswith(".<", i):
        return (line, i)
    try:
        gen_end = line.index(">", i+2)
    except ValueError:
        # Syntax Error -- let Python notice and raise the error
        i += 2
        return (line, i)

    gen_begin = i
    while gen_begin > 0 and line[gen_begin-1] != ';':
        gen_begin -= 1

    # parse out the object name and the list of generator names
    gen_obj = line[gen_begin:i].strip()
    gen_list = [s.strip() for s in line[i+2:gen_end].split(',')]
    for g in gen_list:
        if (not g.isalnum() and not g.replace("_","").isalnum()) or len(g) == 0 or not g[0].isalpha():
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

    gens_assignment = '; (%s,) = %s._first_ngens(%s)' % (gen_vars, gen_obj, gen_vars.count(',')+1)
    ring_assignment = ""
    # Find where the parenthesis of the constructor ends
    if line[:c].rstrip()[-1] == ']':
        # brackets constructor
        c0 = line[:c].find(']')
        d0 = line[:c0].rfind('[')
        if c0 == -1:
            raise SyntaxError, 'constructor must end with ) or ]'

        in_square_brackets = line[d0+1:c0]
        if in_square_brackets.strip() == '':
            # as a convenience to the user, 'K.<a> = ZZ[]' -> 'K.<a> = ZZ["a"]'
            in_square_brackets = "'%s'" % gen_vars

        ring_assignment = '%s%s%s' % (line[:i] + line[gen_end+1:d0+1], in_square_brackets, line[c0:c])
    elif line[:c].rstrip()[-1] == ')':
        c0 = line[:c].rfind(')')
        # General constructor -- rewrite the input line as two commands
        # We have to determine whether or not to put a comma before
        # the list of names.  We do this only if there are already
        # arguments to the constructor.  Some constructors have no
        # arguments, e.g., "K.<a> = f.root_field(  )"
        c1 = line[:c0].rfind('(')
        in_parentheses = line[c1+1:c0].strip()
        if len(in_parentheses) > 0:
            sep = ','
        else:
            sep = ''
        ring_assignment = '%s%snames=%s)' % (line[:i] + line[gen_end+1:c0], sep, gen_names)
    else:
        ring_assignment = line[:i] + line[gen_end+1:c]

    line = ring_assignment + gens_assignment + line[c:]
    i += 1

    return (line, i)

eq_chars_pre = ["=", "!", ">", "<", "+", "-", "*", "/", "^"]

def preparse(line, reset=True, do_time=False, ignore_prompts=False):
    r"""
    sage: preparse("ZZ.<x> = ZZ['x']")
    "ZZ = ZZ['x']; (x,) = ZZ._first_ngens(Integer(1))"
    sage: preparse("ZZ.<x> = ZZ['y']")
    "ZZ = ZZ['y']; (x,) = ZZ._first_ngens(Integer(1))"
    sage: preparse("ZZ.<x,y> = ZZ[]")
    "ZZ = ZZ['x, y']; (x, y,) = ZZ._first_ngens(Integer(2))"
    sage: preparse("ZZ.<x,y> = ZZ['u,v']")
    "ZZ = ZZ['u,v']; (x, y,) = ZZ._first_ngens(Integer(2))"
    sage: preparse("ZZ.<x> = QQ[2^(1/3)]")
    'ZZ = QQ[Integer(2)**(Integer(1)/Integer(3))]; (x,) = ZZ._first_ngens(Integer(1))'
    sage: QQ[2^(1/3)]
    Number Field in a with defining polynomial x^3 - 2
    """
    try:
        # [1,2,..,n] notation
        L, literals = strip_string_literals(line)
        L = parse_ellipsis(L, preparse_step=False)
        line = L % literals
    except SyntaxError:
        pass

    line = implicit_mul(line)

    # find where the parens are for function assignment notation
    oparen_index = -1
    cparen_index = -1

    # for caclulus function notation:
    paren_level = 0
    first_eq_index = -1

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

    in_args = False

    if reset:
        in_single_quote = False
        in_double_quote = False
        in_triple_quote = False


    if ignore_prompts:
        # Get rid of leading sage: and >>> so that pasting of examples from
        # the documentation works.
        line = strip_prompts(line)

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

        elif line[i] == ";" and not in_quote():
            line = line[:i+1] + preparse(line[i+1:], reset, do_time, ignore_prompts)
            i = len(line)
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
        elif not in_quote() and (line[i:i+2] == ".<"):
            line, i = parse_generators(line, i)

        # Support for calculus-like function assignment, the line
        # "f(x,y,z) = sin(x^3 - 4*y) + y^x"
        # gets turnd into
        # '_=var("x,y,z");f=symbolic_expression(sin(x**Integer(3) - Integer(4)*y) + y**x).function(x,y,z)'
        # AUTHORS:
        #   - Bobby Moretti: initial version - 02/2007
        #   - William Stein: make variables become defined if they aren't already defined.

        elif (line[i] == "(") and not in_quote():
            paren_level += 1
            # we need to make sure that this is the first open paren we find
            if oparen_index == -1:
                oparen_index = i
            i += 1
            continue

        elif (line[i] == ")") and not in_quote():
            cparen_index = i
            paren_level -= 1
            i += 1
            continue

        elif (line[i] == "=") and paren_level == 0 and not in_quote():
            if first_eq_index == -1:
                first_eq_index = i
            eq = i

            if cparen_index == -1:
                i += 1
                continue

            # make sure the '=' sign is on its own, representing assignment
            if eq+1 < len(line) and (line[eq-1] in eq_chars_pre or line[eq+1] == '='):
                i += 1
                continue

            line_before = line[:oparen_index].strip()
            if line_before == "" or not line_before.isalnum():
                i += 1
                continue

            if eq != first_eq_index:
                i += 1
                continue

            vars_end = cparen_index
            vars_begin = oparen_index+1

            # figure out where the line ends
            line_after = line[vars_end+1:]
            try:
                a = line.index("#")
            except ValueError:
                a = len(line)

            try:
                b = line.index(";")
            except ValueError:
                b = len(line)

            a =  min(a,b)

            vars = line[vars_begin:vars_end].split(",")
            vars = [v.strip() for v in vars]
            b = []


            # construct the parsed line
            k = line[:vars_begin-1].rfind(';')
            if k == -1:
                k = len(line) - len(line.lstrip())
            b.append(line[:k])
            b.append('_=var("%s");'%(','.join(vars)))
            b.append(line[k:vars_begin-1])
            b.append('=')
            b.append('symbolic_expression(')
            b.append(line[eq+1:a].strip())
            b.append(').function(')
            b.append(','.join(vars))
            b.append(')')
            b.append(line[a:])

            line =  ''.join(b)

            # i should get set to the position of the first paren after the new
            # assignment operator
            n = len(line_before)
            i = line.find('=', n) + 2

            continue
        #
        ####### END CALCULUS ########

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
               not in_quote():

            if i < len(line)-1 and line[i] == '\\':
                j = i+1
                while j < len(line) and line[j].isspace():
                    j += 1

                while j < len(line) and not line[j] in '*/;:\\#\'"':
                    j += 1
                line = line[:i] + "._backslash_(" + line[i+1:j] + ')' + line[j:]

            elif (line[i].isdigit() or \
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
                L = interpreter.load_cython(name_load)
            else:
                #print "Loading of '%s' not implemented (load .py, .spyx, and .sage files)"%name_load
                L = 'load("%s")'%name_load

        M = preparse(L, reset=(i==0), do_time=do_time, ignore_prompts=ignore_prompts)
        F.append(M)
        i += 1
    # end while

    return '\n'.join(F)

def implicit_mul(code, level=5):
    """
    Insert explicit *'s for implicit multiplication.

    INPUT:
        code  -- the code with missing *'s
        level -- how agressive to be in placing *'s
                   0) Do nothing
                   1) numeric followed by alphanumeric
                   2) closing parentheses followed by alphanumeric
                   3) Spaces between alphanumeric
                  10) Adjacent parentheses (may mangle call statements)

    EXAMPLES:
        sage: from sage.misc.preparser import implicit_mul
        sage: implicit_mul('(2x^2-4x+3)a0')
        '(2*x^2-4*x+3)*a0'
        sage: implicit_mul('a b c in L')
        'a*b*c in L'
        sage: implicit_mul('1r + 1e3 + 5exp(2)')
        '1r + 1e3 + 5*exp(2)'
        sage: implicit_mul('f(a)(b)', level=10)
        'f(a)*(b)'
    """
    def re_no_keyword(pattern, code):
        for _ in range(2): # do it twice in because matches don't overlap
            for m in reversed(list(re.finditer(pattern, code))):
                left, right = m.groups()
                if left not in keywords and right not in keywords:
                    code = "%s%s*%s%s" % (code[:m.start()],
                                          left,
                                          right,
                                          code[m.end():])
        return code

    code, literals = strip_string_literals(code)
    if level >= 1:
        no_mul_token = " '''_no_mult_token_''' "
        code = re.sub(r'( *)time ', r'\1time %s' % no_mul_token, code)  # first word may be magic 'time'
        code = re.sub(r'\b(\d+(?:\.\d+)?(?:e\d+)?)([rR]\b)', r'\1%s\2' % no_mul_token, code)  # exclude such things as 10r
        code = re.sub(r'\b(\d+(?:\.\d+)?)e([-\d])', r'\1%se%s\2' % (no_mul_token, no_mul_token), code)  # exclude such things as 1e5
        code = re_no_keyword(r'\b(\d+(?:\.\d+)?) *([a-zA-Z_(]\w*)\b', code)
    if level >= 2:
        code = re.sub(r'(\%\(L\d+\))s', r'\1%ss%s' % (no_mul_token, no_mul_token), code) # literal strings
        code = re_no_keyword(r'(\)) *(\w+)', code)
    if level >= 3:
        code = re_no_keyword(r'(\w+) +(\w+)', code)
    if level >= 10:
        code = re.sub(r'\) *\(', ')*(', code)
    code = code.replace(no_mul_token, '')
    return code % literals
