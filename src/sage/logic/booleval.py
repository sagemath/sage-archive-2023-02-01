r"""
Evaluation of Boolean Formulas

AUTHORS:

- Chris Gorecki (2006): initial version

- Paul Scurek (2013-08-05): updated docstring formatting

EXAMPLES:

We can assign values to the variables and evaluate a formula::

    sage: import sage.logic.booleval as booleval
    sage: t = ['|', ['&', 'a', 'b'], ['&', 'a', 'c']]
    sage: d = {'a' : True, 'b' : False, 'c' : True}
    sage: booleval.eval_formula(t, d)
    True

We can change our assignment of values by modifying the dictionary::

    sage: d['a'] = False
    sage: booleval.eval_formula(t, d)
    False
"""
#*****************************************************************************
#       Copyright (C) 2006 Chris Gorecki <chris.k.gorecki@gmail.com>
#       Copyright (C) 2013 Paul Scurek <scurek86@gmail.com>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#  as published by the Free Software Foundation; either version 2 of
#  the License, or (at your option) any later version.
#                  http://www.gnu.org/licenses/
#*****************************************************************************

import logicparser

# dictionary containing variable keys and boolean values
__vars = {}


def eval_formula(tree, vdict):
    r"""
    Evaluate the tree and return a boolean value.

    INPUT:

    - ``tree`` -- a list of three elements corresponding to a branch of a
      parse tree

    - ``vdict`` -- a dictionary containing variable keys and boolean values

    OUTPUT:

    The result of the evaluation as a boolean value.

    EXAMPLES:

    This example illustrates evaluating a boolean formula::

        sage: import sage.logic.booleval as booleval
        sage: t = ['|', ['&', 'a', 'b'], ['&', 'a', 'c']]
        sage: d = {'a' : True, 'b' : False, 'c' : True}
        sage: booleval.eval_formula(t, d)
        True

    ::

        sage: d['a'] = False
        sage: booleval.eval_formula(t, d)
        False
    """
    global __vars
    __vars = vdict
    b = logicparser.apply_func(tree, eval_f)
    return b

def eval_f(tree):
    r"""
    Evaluate the tree.

    INPUT:

    - ``tree`` -- a list of three elements corresponding to a branch of a
      parse tree

    OUTPUT:

    The result of the evaluation as a boolean value.

    EXAMPLES:

    This example illustrates how to evaluate a parse tree::

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
    Evaluate ``lv`` and ``rv`` according to the operator ``op``.

    INPUT:

    - ``op`` -- a string or character representing a boolean operator

    - ``lv`` -- a boolean or variable

    - ``rv`` -- a boolean or variable

    OUTPUT:

    The evaluation of ``lv op rv`` as a boolean value.

    EXAMPLES:

    We can evaluate an operator given the values on either side::

        sage: import sage.logic.booleval as booleval
        sage: booleval.eval_op('&', True, False)
        False

        sage: booleval.eval_op('^', True, True)
        False

        sage: booleval.eval_op('|', False, True)
        True
    """
    lval = rval = None
    if lv == False:
        lval = False
    elif lv == True:
        lval = True
    elif lv is not None:
        lval = __vars[lv]

    if rv == False:
        rval = False
    elif rv == True:
        rval = True
    elif rv is not None:
        rval = __vars[rv]

    if op == '~':
        return not lval
    elif op == '&':
        return lval and rval
    elif op == '|':
        return lval or rval
    elif op == '^':
        return lval ^ rval
    elif op == '->':
        return (not lval) or rval
    elif op == '<->':
        return (not lval or rval) and (not rval or lval)
    else:  # one variable
        return __vars[op]
