r"""
Module that creates and modifies parse trees of well formed
boolean formulas.

AUTHORS:
    -- Chris Gorecki

EXAMPLES:
    sage: import sage.logic.logicparser as logicparser
    sage: s = 'a|b&c'
    sage: t = logicparser.parse(s)
    sage: t
    (['|', 'a', ['&', 'b', 'c']], ['a', 'b', 'c'])
"""

from types import *
import string

__symbols = '()&|~<->^'
__op_list = ['~', '&', '|', '^', '->', '<->']

def parse(s):
    r"""
    This function produces a parse tree from a boolean formula s.

    INPUT:
        s -- a string containing a boolean formula.

    OUTPUT:
        Returns the tuple (parse tree of s, variables in s).

    EXAMPLES:
        sage: import sage.logic.logicparser as logicparser
        sage: s = 'a|b&c'
        sage: t = logicparser.parse(s)
        sage: t
        (['|', 'a', ['&', 'b', 'c']], ['a', 'b', 'c'])
    """
    toks, vars_order = tokenize(s)
    tree = tree_parse(toks)
    # special case of tree == single variable
    if isinstance(tree, StringType):
        return ['&', tree, tree], vars_order
    return tree, vars_order

def polish_parse(s):
    r"""
    This function produces a full syntax parse tree from
    a boolean formlua s.

    INPUT:
        s -- a string containing a boolean expression

    OUTPUT:
        Returns a list representing the full syntax parse tree of s, with
        the symbols in s arranged from left to right as in polish notation

    EXAMPLES:
        sage: import sage.logic.logicparser as logicparser
        sage: s = 'a|~~b'
        sage: t = logicparser.polish_parse(s)
        sage: t
        ['|', 'a', ['~', ['~', 'b']]]
    """
    toks, vars_order = tokenize(s)
    tree = tree_parse(toks, polish = True)
    # special case where the formula s is a single variable
    if isinstance(tree, StringType):
        return vars_order
    return tree

def tokenize(s):
    r"""
    This function produces a string of tokens from s.

    INPUT:
        s -- a string containing a boolean formula.

    OUTPUT:
        Returns a tuple consisting of (tokens in s, variables in s).

    EXAMPLES:
        sage: import sage.logic.logicparser as logicparser
        sage: s = 'a|b&c'
        sage: t = logicparser.tokenize(s)
        sage: t
        (['(', 'a', '|', 'b', '&', 'c', ')'], ['a', 'b', 'c'])
    """
    i = 0
    toks = ['(']
    vars_order = []

    while i < len(s):
        tok = ""
        skip = valid = 1
        if s[i] in '()~&|^':
            tok = s[i]
        elif s[i:i + 2] == '->':
            tok = '->'
            skip = 2
        elif s[i:i + 3] == '<->':
            tok = '<->'
            skip = 3
        # check to see if '-', '<' or '>' are used incorretly
        elif s[i] in '<->':
            msg = "'%s' can only be used as part of the operators '<->' or '->'." % (s[i])
            raise SyntaxError, msg
        if len(tok) > 0:
            toks.append(tok)
            i += skip
            continue
        else:
            # token is a variable name
            if s[i] == ' ':
                i += 1
                continue

            while i < len(s) and s[i] not in __symbols and s[i] != ' ':
                tok += s[i]
                i += 1

            if len(tok) > 0:
                if tok[0] not in string.letters:
                    valid = 0
                for c in tok:
                    if c not in string.letters and c not in string.digits and c != '_':
                        valid = 0

            if valid == 1:
                toks.append(tok)
                if tok not in vars_order:
                    vars_order.append(tok)
            else:
                msg = 'invalid variable name ' + tok
                msg += ": identifiers must begin with a letter and contain only "
                msg += "alphanumerics and underscores"
                raise NameError, msg

    toks.append(')')
    return toks, vars_order

def tree_parse(toks, polish = False):
    r"""
    This function produces a parse tree from the tokens in toks.

    INPUT:
        toks -- a list of tokens.
        polish -- (default : False) A boolean value.
                  When true, causes the function to produce the full
                  syntax parse tree from toks

    OUTPUT:
        -- Returns a parse tree of the tokens toks.
        -- If polish = True is passed as an argument,
           returns the full syntax parse tree of toks


    EXAMPLES:
        sage: import sage.logic.logicparser as logicparser
        sage: t = ['(', 'a', '|', 'b', '&', 'c', ')']
        sage: logicparser.tree_parse(t)
        ['|', 'a', ['&', 'b', 'c']]

        note the difference when polish = True
        sage: t = ['(', 'a', '->', '~', '~', 'b', ')']
        sage: logicparser.tree_parse(t)
        ['->', 'a', 'b']
        sage: t = ['(', 'a', '->', '~', '~', 'b', ')']
        sage: logicparser.tree_parse(t, polish = True)
        ['->', 'a', ['~', ['~', 'b']]]
    """
    stack = []
    for tok in toks:
        stack.append(tok)
        if tok == ')':
            lrtoks = []
            while tok != '(':
                tok = stack.pop()
                lrtoks.insert(0, tok)
            branch = parse_ltor(lrtoks[1:-1], polish = polish)
            stack.append(branch)
    return stack[0]

def parse_ltor(toks, n = 0, polish = False):
    r"""
    This function produces a parse tree from the tokens in toks under
    the precondition that each token in toks is atomic.

    INPUT:
        toks -- a list of tokens.
        n -- an integer representing which order of operations
             are occurring.
        polish -- (default : False) A boolean value. If True,
                  double negations are not cancelled and negated statements
                  are turned into lists of length 2

    OUTPUT:
       -- Returns a parse tree of the tokens toks.
       -- If polish = True is passed as an argument, returns the full syntax parse
          tree of the tokens toks

    EXAMPLES:
        sage: import sage.logic.logicparser as logicparser
        sage: t = ['a', '|', 'b', '&', 'c']
        sage: logicparser.parse_ltor(t)
        ['|', 'a', ['&', 'b', 'c']]

        note the difference when polish = True is passed
        sage: t = ['a', '->', '~', '~', 'b']
        sage: logicparser.parse_ltor(t)
        ['->', 'a', 'b']
        sage: t = ['a', '->', '~', '~', 'b']
        sage: logicparser.parse_ltor(t, polish = True)
        ['->', 'a', ['~', ['~', 'b']]]

    """
    i = 0
    for tok in toks:
        if tok == __op_list[n]:
            if tok == '~':
                if not polish:
                    # cancel double negations
                    if toks[i] == '~' and toks[i + 1] == '~':
                        del toks[i]
                        del toks[i]
                        return parse_ltor(toks, n)
                    args = [toks[i], toks[i + 1], None]
                    toks[i] = args
                    del toks[i + 1]
                    return parse_ltor(toks, n)
                # This executes when creating the full syntax parse tree
                else:
                    j = i
                    while toks[j] == '~':
                        j += 1
                    while j > i:
                        args = [toks[j - 1], toks[j]]
                        toks[j - 1] = args
                        del toks[j]
                        j -= 1
                    return parse_ltor(toks, n = n, polish = polish)
            else:
                args = [toks[i - 1], toks[i], toks[i + 1]]
                toks[i - 1] = [args[1], args[0], args[2]]
                del toks[i]
                del toks[i]
                return parse_ltor(toks, n)
        i += 1
    if n + 1 < len(__op_list):
        return parse_ltor(toks, n + 1)
    if len(toks) > 1:
        raise SyntaxError
    return toks[0]

def apply_func(tree, func):
    r"""
    This function applies func to each node of tree.

    INPUT:
        tree -- a parse tree.
        func -- a function to be applied to tree.

    OUTPUT:
        Returns a parse tree after func has been applied
        to it.

    EXAMPLES:
        sage: import sage.logic.logicparser as logicparser
        sage: t = ['|', ['&', 'a', 'b'], ['&', 'a', 'c']]
        sage: f = lambda t: t
        sage: logicparser.apply_func(t, f)
        ['|', ['&', 'a', 'b'], ['&', 'a', 'c']]
    """
    if type(tree[1]) is ListType and type(tree[2]) is ListType:
        lval = apply_func(tree[1], func)
        rval = apply_func(tree[2], func)
    elif type(tree[1]) is ListType:
        lval = apply_func(tree[1], func)
        rval = tree[2]
    elif type(tree[2]) is ListType:
        lval = tree[1]
        rval = apply_func(tree[2], func)
    else:
        lval = tree[1]
        rval = tree[2]
    return func([tree[0], lval, rval])


