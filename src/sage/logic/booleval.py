r"""
Module associated with logicparser to do evaluations of boolean formulas.

AUTHORS:
    -- Chris Gorecki

EXAMPLES:
    sage: import sage.logic.booleval as booleval
    sage: t = ['|', ['&', 'a', 'b'], ['&', 'a', 'c']]
    sage: d = {'a' : True, 'b' : False, 'c' : True}
    sage: booleval.eval_formula(t, d)
    True
    sage: d['a'] = False
    sage: booleval.eval_formula(t, d)
    False
"""

from types import *
import logicparser

__vars = {}

def eval_formula(tree, vars):
    r"""
    Evaluates the tree using the boolean values contained in dictionary
    and returns a single boolean value.

    INPUT:
        tree -- a list of three elements corresponding to a branch of a
                parse tree.
        vars -- a dictionary containing variable keys and boolean values.

    OUTPUT:
        Returns the boolean evaluation of a boolean formula.

    EXAMPLES:
        sage: import sage.logic.booleval as booleval
        sage: t = ['|', ['&', 'a', 'b'], ['&', 'a', 'c']]
        sage: d = {'a' : True, 'b' : False, 'c' : True}
        sage: booleval.eval_formula(t, d)
        True
        sage: d['a'] = False
        sage: booleval.eval_formula(t, d)
        False
    """
    global __vars
    __vars = vars
    b = logicparser.apply_func(tree, eval_f)
    return b

def eval_f(tree):
    r"""
    This function can be applied to a parse tree to evaluate it.

    INPUT:
        tree -- a list of three elements corresponding to a branch of a
                parse tree.

    OUTPUT:
         Returns a boolean evaluation of the tree.

    EXAMPLES:
         sage: import sage.logic.booleval as booleval
         sage: booleval.eval_f(['&', True, False])
         False
         sage: booleval.eval_f(['^', True, True])
         False
         sage: booleval.eval_f(['|', False, True])
         True
    """
    return eval_op(tree[0], tree[1], tree[2])

def eval_op(op, lv, rv):
    r"""
    This function evaluates lv and rv according to the operator op.

    INPUT:
        op -- a string/char representing a boolean operator.
        lv -- a boolean or variable.
        rv -- a boolean or variable.

    OUTPUT:
        Returns the evaluation of lv op rv.

    EXAMPLES:
        sage: import sage.logic.booleval as booleval
        sage: booleval.eval_op('&', True, False)
        False
        sage: booleval.eval_op('^', True, True)
        False
        sage: booleval.eval_op('|', False, True)
        True
    """
    lval = rval = None
    if(lv == False):
        lval = False
    elif(lv == True):
        lval = True
    elif(lv != None):
        lval = __vars[lv]

    if(rv == False):
        rval = False
    elif(rv == True):
        rval = True
    elif(rv != None):
        rval = __vars[rv]

    if(op == '~'):
        return not lval
    elif(op == '&'):
        return lval and rval
    elif(op == '|'):
        return lval or rval
    elif(op == '^'):
        return lval ^ rval
    elif(op == '->'):
        return (not lval) or rval
    elif(op == '<->'):
        return (not lval or rval) and (not rval or lval)
    else:  #one variable
        return __vars[op]