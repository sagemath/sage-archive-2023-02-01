r"""
Module that creates and modifies parse trees of well formed boolean formulas.

A parse tree of a boolean formula is a nested list, where each branch is either
a single variable, or a formula composed of either two variables and a binary
operator or one variable and a unary operator. The function parse produces
a parse tree that is simplified for the purposes of more efficient truth value
evaluation. The function polish_parse produces the full parse tree of a boolean
formula which is used in functions related to proof and inference.  That is,
parse is meant to be used with functions in the logic module that perform
semantic operations on a boolean formula, and polish_parse is to be used with
functions that perform syntactic operations on a boolean formula.

AUTHORS:

- Chris Gorecki (2007): initial version

- Paul Scurek (2013-08-01): added polish_parse, cleaned up python code,
  updated docstring formatting

- Paul Scurek (2013-08-06): added recover_formula, recover_formula_internal,
  prefix_to_infix, to_infix_internal

EXAMPLES:

Find the parse tree and variables of a string representation of a boolean formula::

    sage: import sage.logic.logicparser as logicparser
    sage: s = 'a|b&c'
    sage: t = logicparser.parse(s)
    sage: t
    (['|', 'a', ['&', 'b', 'c']], ['a', 'b', 'c'])

Find the full syntax parse tree of a string representation of a boolean formula::

    sage: import sage.logic.logicparser as logicparser
    sage: s = '(a&b)->~~c'
    sage: logicparser.polish_parse(s)
    ['->', ['&', 'a', 'b'], ['~', ['~', 'c']]]

Find the tokens and distinct variables of a boolean formula::

    sage: import sage.logic.logicparser as logicparser
    sage: s = '~(a|~b)<->(c->c)'
    sage: logicparser.tokenize(s)
    (['(', '~', '(', 'a', '|', '~', 'b', ')', '<->', '(', 'c', '->', 'c', ')', ')'], ['a', 'b', 'c'])

Find the parse tree of a boolean formula from a list of the formula's tokens::

    sage: import sage.logic.logicparser as logicparser
    sage: t = ['(', 'a', '->', '~', 'c', ')']
    sage: logicparser.tree_parse(t)
    ['->', 'a', ['~', 'c', None]]
    sage: r = ['(', '~', '~', 'a', '|', 'b', ')']
    sage: logicparser.tree_parse(r)
    ['|', 'a', 'b']

Find the full syntax parse tree of a boolean formula from a list of tokens::

    sage: import sage.logic.logicparser as logicparser
    sage: t = ['(', 'a', '->', '~', 'c', ')']
    sage: logicparser.tree_parse(t, polish = True)
    ['->', 'a', ['~', 'c']]
    sage: r = ['(', '~', '~', 'a', '|', 'b', ')']
    sage: logicparser.tree_parse(r, polish = True)
    ['|', ['~', ['~', 'a']], 'b']


"""
#*****************************************************************************
#       Copyright (C) 2007 Chris Gorecki <chris.k.gorecki@gmail.com>
#       Copyright (C) 2013 Paul Scurek <scurek86@gmail.com>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#  as published by the Free Software Foundation; either version 2 of
#  the License, or (at your option) any later version.
#                  http://www.gnu.org/licenses/
#*****************************************************************************

from types import *
import string

__symbols = '()&|~<->^'
__op_list = ['~', '&', '|', '^', '->', '<->']

def parse(s):
    r"""
    Return a parse tree from a boolean formula s.

    INPUT:

    - ``s`` -- a string containing a boolean formula.

    OUTPUT:

    A list containing the prase tree and a list containing the
    variables in a boolean formula in this order:

    1. the list containing the pase tree
    2. the list containing the variables

    EXAMPLES:

    This example illustrates how to produce the parse tree of a boolean formula s.

    ::

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
    Return the full syntax parse tree from a boolean formula s.

    INPUT:

    - ``s`` -- a string containing a boolean expression

    OUTPUT:

    The full syntax parse tree as a nested list

    EXAMPLES:

    This example illustrates how to find the full syntax parse tree of a boolean formula.

    ::

        sage: import sage.logic.logicparser as logicparser
        sage: s = 'a|~~b'
        sage: t = logicparser.polish_parse(s)
        sage: t
        ['|', 'a', ['~', ['~', 'b']]]

    AUTHORS:

    - Paul Scurek (2013-08-03)
    """
    toks, vars_order = tokenize(s)
    tree = tree_parse(toks, polish = True)
    # special case where the formula s is a single variable
    if isinstance(tree, StringType):
        return vars_order
    return tree

def recover_formula(prefix_tree):
    r"""
    Recover the formula from a parse tree in prefix form.

    INPUT:

    - ``prefix_tree`` -- a list. This is a full syntax parse
      tree in prefix form.

    OUTPUT:

    The formula as a string

    EXAMPLES:

    This example illustrates the recovery of a formula from a parse tree.

    ::

        sage: import sage.logic.propcalc as propcalc
        sage: import sage.logic.logicparser as logicparser
        sage: t = ['->', ['&', 'a', ['~', ['~', 'c']]], ['~', ['|', ['~', 'c'], 'd']]]
        sage: logicparser.recover_formula(t)
        '(a&~~c)->~(~c|d)'

    ::

        sage: f = propcalc.formula("a&(~~c|d)")
        sage: logicparser.recover_formula(f.full_tree())
        'a&(~~c|d)'

    ::

        sage: r = ['~', 'a']
        sage: logicparser.recover_formula(r)
        '~a'

    ::

        sage: s = ['d']
        sage: logicparser.recover_formula(s)
        'd'

    .. NOTES::

        The function :func:`polish_parse` may be passed as an argument,
        but :func:`tree_parse` may not unless the parameter ``polish``
        is set equal to True.

    AUTHORS:

    - Paul Scurek (2013-08-06)
    """
    formula = ''
    if not isinstance(prefix_tree, list):
        raise Exception, "recover_formula takes a parse tree as input"
    else:
        formula = apply_func(prefix_tree, recover_formula_internal)
        if prefix_tree[0] == '~' or len(prefix_tree) == 1:
            return formula
        else:
            return formula[1:-1]

def recover_formula_internal(prefix_tree):
    r"""
    Recover the formula from a parse tree in prefix form.

    INPUT:

    - ``prefix_tree`` -- a list. This is a simple tree
      with at most one operator in prefix form.

    OUTPUT:

    The formula as a string

    EXAMPLES:

    This example illustrates recovering the formula from a parse tree.

    ::

        sage: import sage.logic.logicparser as logicparser
        sage: import sage.logic.propcalc as propcalc
        sage: t = ['->', 'a', 'b']
        sage: logicparser.recover_formula_internal(t)
        '(a->b)'

    ::

        sage: r = ['~', 'c']
        sage: logicparser.recover_formula_internal(r)
        '~c'

    ::

        sage: s = ['d']
        sage: logicparser.recover_formula_internal(s)
        'd'

    We can pass :func:`recover_formula_internal` as an argument in :func:`apply_func`.

    ::

        sage: f = propcalc.formula("~(d|c)<->(a&~~~c)")
        sage: logicparser.apply_func(f.full_tree(), logicparser.recover_formula_internal)
        '(~(d|c)<->(a&~~~c))'

    .. NOTE::

        This function is for internal use by logicparser.py.  The function
        recovers the formula of a simple parse tree in prefix form.  A
        simple parse tree contains at most one operator.


        The function :func:`polish_parse` may be passed as an argument,
        but :func:`tree_parse` may not unless the parameter ``polish``
        is set equal to True.

    AUTHORS:

    - Paul Scurek (2013-08-06)
    """
    if len(prefix_tree) == 3:
        return '(' + prefix_tree[1] + prefix_tree[0] + prefix_tree[2] + ')'
    else:
        return ''.join(prefix_tree)

def prefix_to_infix(prefix_tree):
    r"""
    Convert a parse tree from prefix form to infix form.

    INPUT:

    - ``prefix_tree`` -- a list. This is a full syntax parse
      tree in prefix form.

    OUTPUT:

    A list containing the tree in infix form

    EXAMPLES:

    This example illustrates converting a prefix tree to an infix tree.

    ::

        sage: import sage.logic.logicparser as logicparser
        sage: import sage.logic.propcalc as propcalc
        sage: t = ['|', ['~', 'a'], ['&', 'b', 'c']]
        sage: logicparser.prefix_to_infix(t)
        [['~', 'a'], '|', ['b', '&', 'c']]

    ::

        sage: f = propcalc.formula("(a&~b)<->~~~(c|d)")
        sage: logicparser.prefix_to_infix(f.full_tree())
        [['a', '&', ['~', 'b']], '<->', ['~', ['~', ['~', ['c', '|', 'd']]]]]

    .. NOTES::

        The function :func:`polish_parse` may be passed as an argument,
        but :func:`tree_parse` may not unless the parameter ``polish``
        is set equal to True.

    AUTHORS:

    - Paul Scurek (2013-08-06)
    """
    if not isinstance(prefix_tree, list):
        raise Exception, "prefix_to_infix takes a parse tree as input"
    return apply_func(prefix_tree, to_infix_internal)

def to_infix_internal(prefix_tree):
    r"""
    Convert a simple parse tree from prefix form to infix form.

    INPUT:

    - ``prefix_tree`` -- a list. This is a simple parse tree
      in prefix form with at most one operator.

    OUTPUT:

    The tree in infix form as a list.

    EXAMPLES:

    This example illustrates converting a simple tree from prefix to infix form.

    ::

        sage: import sage.logic.logicparser as logicparser
        sage: import sage.logic.propcalc as propcalc
        sage: t = ['|', 'a', 'b']
        sage: logicparser.to_infix_internal(t)
        ['a', '|', 'b']

    We can pass :func:`to_infix_internal` as an argument in :func:`apply_func`.

    ::

        sage: f = propcalc.formula("(a&~b)<->~~~(c|d)")
        sage: logicparser.apply_func(f.full_tree(), logicparser.to_infix_internal)
        [['a', '&', ['~', 'b']], '<->', ['~', ['~', ['~', ['c', '|', 'd']]]]]

    .. NOTES::

        This function is for internal use by logicparser.py. It converts a simple
        parse tree from prefix form to infix form.  A simple parse tree contains
        at most one operator.

        The function :func:`polish_parse` may be passed as an argument,
        but :func:`tree_parse` may not unless the parameter ``polish``
        is set equal to True.

    AUTHORS:

    - Paul Scurek (2013-08-06)
    """
    if prefix_tree[0] != '~' and len(prefix_tree) == 3:
        return [prefix_tree[1], prefix_tree[0], prefix_tree[2]]
    else:
        return prefix_tree

def tokenize(s):
    r"""
    Return the tokens and the distinct variables appearing in a boolean formula s.

    INPUT:

    - ``s`` -- a string representation of a boolean formula

    OUTPUT:

    The tokens and variables as an ordered pair of lists in the following order:

    1. A list containing the tokens of s, in the order they appear in s
    2. A list containing the distinct variables in s, in the order they appearn in s

    EXAMPLES:

    This example illustrates how to tokenize a string representation of a boolean formula.

    ::

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
        # check to see if '-', '<' or '>' are used incorrectly
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
    Return a parse tree from the tokens in toks.

    INPUT:

    - ``toks`` -- a list of tokens from a boolean formula

    - ``polish`` -- (default: False) a boolean.  When true, tree_parse will return
      the full syntax parse tree.

    OUTPUT:

    A parse tree in the form of a nested list that depends on ``polish`` as follows:

    polish == False -- Return a simplified parse tree.

    polish == True -- Return the full syntax parse tree.

    EXAMPLES:

    This example illustrates the use of tree_parse when polish == False.

    ::

        sage: import sage.logic.logicparser as logicparser
        sage: t = ['(', 'a', '|', 'b', '&', 'c', ')']
        sage: logicparser.tree_parse(t)
        ['|', 'a', ['&', 'b', 'c']]

    We now demonstrate the use of tree_parse when polish == True.

    ::

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
    Return a parse tree from toks, where each token in toks is atomic.

    INPUT:

    - ``toks`` -- a list of tokens. Each token is atomic.

    - ``n`` -- (default: 0) an integer representing which order of
      operations are occurring

    - ``polish`` -- (default: False) a boolean.  When true, double negations
      are not cancelled and negated statements are turned into list of length two.

    OUTPUT:

    The parse tree as a nested list that depends on ``polish`` as follows:

    polish == False - Return a simplified parse tree.

    polish == True - Return the full syntax parse tree.

    EXAMPLES:

    This example illustrates the use of parse_ltor when polish == False.

    ::

        sage: import sage.logic.logicparser as logicparser
        sage: t = ['a', '|', 'b', '&', 'c']
        sage: logicparser.parse_ltor(t)
        ['|', 'a', ['&', 'b', 'c']]

    ::

        sage: import sage.logic.logicparser as logicparser
        sage: t = ['a', '->', '~', '~', 'b']
        sage: logicparser.parse_ltor(t)
        ['->', 'a', 'b']

    We now repeat the previous example, but with polish == True.

    ::

        sage: import sage.logic.logicparser as logicparser
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
    Apply func to each node of tree.  Return a new parse tree.

    INPUT:

    - ``tree`` -- a parse tree of a boolean formula

    - ``func`` -- a function to be applied to each node of tree.  This may
      be a function that comes from elsewhere in the logic module.

    OUTPUT:

    The new parse tree in the form of a nested list

    EXAMPLES:

    This example uses :func:`apply_func` where ``func`` switches two entries of tree.

    ::

        sage: import sage.logic.logicparser as logicparser
        sage: t = ['|', ['&', 'a', 'b'], ['&', 'a', 'c']]
        sage: f = lambda t: [t[0], t[2], t[1]]
        sage: logicparser.apply_func(t, f)
        ['|', ['&', 'c', 'a'], ['&', 'b', 'a']]
    """
    # used when full syntax parse tree is passed as argument
    if len(tree) == 1:
        return func(tree)
    # used when full syntax parse tree is passed as argument
    elif len(tree) == 2:
        rval = apply_func(tree[1], func)
        return func([tree[0], rval])
    elif isinstance(tree[1], list) and isinstance(tree[2], list):
        lval = apply_func(tree[1], func)
        rval = apply_func(tree[2], func)
    elif isinstance(tree[1], list):
        lval = apply_func(tree[1], func)
        rval = tree[2]
    elif isinstance(tree[2], list):
        lval = tree[1]
        rval = apply_func(tree[2], func)
    else:
        lval = tree[1]
        rval = tree[2]
    return func([tree[0], lval, rval])


