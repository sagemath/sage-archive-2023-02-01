r"""
    Module associated with booleval, logictable, and logicparser, to do
    boolean evaluation of boolean formulas.
    Formulas consist of the operators &, |, ~, ^, ->, <->, corresponding
    to and, or, not, xor, if...then, if and only if.  Operators can
    be applied to variables that consist of a leading letter and trailing
    underscores and alphanumerics.  Parentheses may be used to
    explicitly show order of operation.

    AUTHORS:
        -- Chris Gorecki

    EXAMPLES:
        sage: import sage.logic.propcalc as propcalc
        sage: f = propcalc.formula("a&((b|c)^a->c)<->b")
        sage: g = propcalc.formula("boolean<->algebra")
        sage: (f&~g).ifthen(f)
        ((a&((b|c)^a->c)<->b)&(~(boolean<->algebra)))->(a&((b|c)^a->c)<->b)

     EXAMPLES:
        sage: import sage.logic.propcalc as propcalc
        sage: f = propcalc.formula("a&((b|c)^a->c)<->b")
        sage: g = propcalc.formula("boolean<->algebra")
        sage: (f&~g).ifthen(f)
        ((a&((b|c)^a->c)<->b)&(~(boolean<->algebra)))->(a&((b|c)^a->c)<->b)

    We can create a truth table from a formula.
        sage: f.truthtable()
        a      b      c      value
        False  False  False  True
        False  False  True   True
        False  True   False  False
        False  True   True   False
        True   False  False  True
        True   False  True   False
        True   True   False  True
        True   True   True   True
        sage: f.truthtable(end=3)
        a      b      c      value
        False  False  False  True
        False  False  True   True
        False  True   False  False
        sage: f.truthtable(start=4)
        a      b      c      value
        True   False  False  True
        True   False  True   False
        True   True   False  True
        True   True   True   True
        sage: propcalc.formula("a").truthtable()
        a      value
        False  False
        True   True

    Now we can evaluate the formula for a given set of inputs.
        sage: f.evaluate({'a':True, 'b':False, 'c':True})
        False
        sage: f.evaluate({'a':False, 'b':False, 'c':True})
        True

    And we can convert a boolean formula to conjunctive normal form.
        sage: f.convert_cnf_table()
        sage: f
        (a|~b|c)&(a|~b|~c)&(~a|b|~c)
        sage: f.convert_cnf_recur()
        sage: f
        (a|~b|c)&(a|~b|~c)&(~a|b|~c)

    We can also simplify an expression.
        sage: f = propcalc.formula("a&((b|c)^a->c)<->b")
        sage: f.truthtable()
        a      b      c      value
        False  False  False  True
        False  False  True   True
        False  True   False  False
        False  True   True   False
        True   False  False  True
        True   False  True   False
        True   True   False  True
        True   True   True   True
        sage: f.simplify()
        (~a&~b)|(a&~b&~c)|(a&b)
        sage: f.truthtable()
        a      b      c      value
        False  False  False  True
        False  False  True   True
        False  True   False  False
        False  True   True   False
        True   False  False  True
        True   False  True   False
        True   True   False  True
        True   True   True   True

    Or determine if an epression is satisfiable, a contradiction, or a tautology.
        sage: f = propcalc.formula("a|b")
        sage: f.is_satisfiable()
        True
        sage: f = f & ~f
        sage: f.is_satisfiable()
        False
        sage: f.is_contradiction()
        True
        sage: f = f | ~f
        sage: f.is_tautology()
        True

    The equlaity operator compares semantic equivlance.
        sage: f = propcalc.formula("(a|b)&c")
        sage: g = propcalc.formula("c&(b|a)")
        sage: f == g
        True
        sage: g = propcalc.formula("a|b&c")
        sage: f == g
        False

    It is an error to create a formula with bad syntax.
        sage: propcalc.formula("")
        Traceback (most recent call last):
        ...
        SyntaxError: malformed statement
        sage: propcalc.formula("a&b~(c|(d)")
        Traceback (most recent call last):
        ...
        SyntaxError: malformed statement
        sage: propcalc.formula("a&&b")
        Traceback (most recent call last):
        ...
        SyntaxError: malformed statement
        sage: propcalc.formula("a&b a")
        Traceback (most recent call last):
        ...
        SyntaxError: malformed statement

        It is also an error to not abide by the naming conventions.
        sage: propcalc.formula("~a&9b")
        Traceback (most recent call last):
        ...
        NameError: invalid variable name 9b: identifiers must begin with a letter and contain only alphanumerics and underscores

"""
#*******************************************************************************************
#copyright (C) 2006 William Stein <wstein@gmail.com>
#copyright (C) 2006 Chris Gorecki <chris.k.gorecki@gmail.com>
#Distributed under the terms of the GNU General Public License (GPL)
#http://www.gnu.org/licenses/
#*******************************************************************************************
import booleval
import logictable
import logicparser
import boolopt
from types import *

latex_operators = [('&', '\\wedge '),
                   ('|', '\\vee '),
                   ('~', '\\neg '),
                   ('^', '\\oplus '),
                   ('<->', '\\leftrightarrow '),
                   ('->', '\\rightarrow ')]

class BooleanFormula:
    __expression = ""
    __tree = []
    __vars_order = []

    def __init__(self, exp, tree, vo):
        r"""
        This function initializes the data fields and is called when a
        new statement is created.

        INPUT:
            self -- the calling object.
            exp -- a string containing the logic expression to be manipulated.
            tree -- a list containing the parse tree of the expression.
            vo -- a list of the variables in the expression in order,
                  with each variable occurring only once.

        OUTPUT:
            Effectively returns an instance of this class.

        EXAMPLES:
        This example illustrates the creation of a statement.
            sage: import sage.logic.propcalc as propcalc
            sage: s = propcalc.formula("a&b|~(c|a)")
            sage: s
            a&b|~(c|a)
        """
        self.__expression = ""
        #remove white space from expression
        for c in exp:
            if(c != ' '):
                self.__expression += c
        self.__tree = tree
        self.__vars_order = vo

    def __repr__(self):
        r"""
        Returns a string representation of this statement.

        INPUT:
            self -- the calling object.

        OUTPUT:
            Returns the string representation of this statement.

        EXAMPLES:
            sage: import sage.logic.propcalc as propcalc
            sage: propcalc.formula("man->monkey&human")
            man->monkey&human
        """
        return self.__expression

    def _latex_(self):
        r"""
        Returns a LaTeX  representation of this statement.

        INPUT:
            self -- the calling object.

        OUTPUT:
            Returns the latex representation of this statement.

        EXAMPLES:
            sage: import sage.logic.propcalc as propcalc
            sage: s = propcalc.formula("man->monkey&human")
            sage: latex(s)
            man\rightarrow monkey\wedge human

            sage: f = propcalc.formula("a & ((~b | c) ^ a -> c) <-> ~b")
            sage: latex(f)
            a\wedge ((\neg b\vee c)\oplus a\rightarrow c)\leftrightarrow \neg b
        """
        latex_expression = self.__expression
        for old, new in latex_operators:
            latex_expression = latex_expression.replace(old, new)
        return latex_expression

    def tree(self):
        return self.__tree

    def __or__(self, other):
        r"""
        Overloads the | operator to 'or' two statements together.

        INPUT:
            self -- the right hand side statement.
            other -- the left hand side statement.

        OUTPUT:
            Returns a new statement that is the first statement logically
                or'ed together.
        EXAMPLES:
            sage: import sage.logic.propcalc as propcalc
            sage: s = propcalc.formula("a&b")
            sage: f = propcalc.formula("c^d")
            sage: s | f
            (a&b)|(c^d)
        """
        return self.add_statement(other, '|')

    def __and__(self, other):
        r"""
        Overloads the & operator to 'and' two statements together.

        INPUT:
            self -- the right hand side statement.
            other -- the left hand side statement.

        OUTPUT:
            Returns a new statement that is the first statement logically
            and'ed together.

        EXAMPLES:
            sage: import sage.logic.propcalc as propcalc
            sage: s = propcalc.formula("a&b")
            sage: f = propcalc.formula("c^d")
            sage: s & f
            (a&b)&(c^d)
        """
        return self.add_statement(other, '&')

    def __xor__(self, other):
        r"""
        Overloads the ^ operator to xor two statements together.

        INPUT:
            self -- the right hand side statement.
            other -- the left hand side statement.

        OUTPUT:
            Returns a new statement that is the first statement logically
            xor'ed together.

        EXAMPLES:
            sage: import sage.logic.propcalc as propcalc
            sage: s = propcalc.formula("a&b")
            sage: f = propcalc.formula("c^d")
            sage: s ^ f
            (a&b)^(c^d)
        """
        return self.add_statement(other, '^')

    def __pow__(self, other):
        r"""
        Overloads the ^ operator to xor two statements together.

        INPUT:
            self -- the right hand side statement.
            other -- the left hand side statement.

        OUTPUT:
            Returns a new statement that is the first statement logically
            xor'ed together.

        EXAMPLES:
            sage: import sage.logic.propcalc as propcalc
            sage: s = propcalc.formula("a&b")
            sage: f = propcalc.formula("c^d")
            sage: s ^ f
            (a&b)^(c^d)
        """
        return self.add_statement(other, '^')

    def __invert__(self):
        r"""
        Overloads the ~ operator to not a statement.

        INPUT:
            self -- the right hand side statement.
            other -- the left hand side statement.

        OUTPUT:
            Returns a new statement that is the first statement logically
            not'ed.

        EXAMPLES:
            sage: import sage.logic.propcalc as propcalc
            sage: s = propcalc.formula("a&b")
            sage: ~s
            ~(a&b)
        """
        exp = '~(' + self.__expression + ')'
        parse_tree, vars_order = logicparser.parse(exp)
        return BooleanFormula(exp, parse_tree, vars_order)

    def ifthen(self, other):
        r"""
        Returns two statements attached by the -> operator.

        INPUT:
            self -- the right hand side statement.
            other -- the left hand side statement.

        OUTPUT:
            Returns a new statement that is the first statement logically
            ifthen'ed together.

       EXAMPLES:
           sage: import sage.logic.propcalc as propcalc
           sage: s = propcalc.formula("a&b")
           sage: f = propcalc.formula("c^d")
           sage: s.ifthen(f)
           (a&b)->(c^d)
        """
        return self.add_statement(other, '->')

    def iff(self, other):
        r"""
        Returns two statements atatched by the <-> operator.

        INPUT:
            self -- the right hand side statement.
            other -- the left hand side statement.

        OUTPUT:
            Returns a new statement that is the first statement logically
            ifandonlyif'ed together.

        EXAMPLES:
            sage: import sage.logic.propcalc as propcalc
            sage: s = propcalc.formula("a&b")
            sage: f = propcalc.formula("c^d")
            sage: s.iff(f)
            (a&b)<->(c^d)
        """
        return self.add_statement(other, '<->')

    def __eq__(self, other):
        r"""
        Overloads the == operator to determine if two expressions
        are logically equivalent.

        INPUT:
            self -- the right hand side statement.
            other -- the left hand side statement.

        OUTPUT:
            Returns true if the left hand side is equivlant to the
            right hand side, and false otherwise.

        EXAMPLES:
            sage: import sage.logic.propcalc as propcalc
            sage: f = propcalc.formula("(a|b)&c")
            sage: g = propcalc.formula("c&(b|a)")
            sage: f == g
            True
            sage: g = propcalc.formula("a|b&c")
            sage: f == g
            False
        """
        return self.equivlant(other)

    def truthtable(self, start = 0, end = -1):
        r"""
        This function returns a truthtable object corresponding to the given
        statement.  Each row of the table corresponds to a binary number, with
        each variable associated to a column of the number, and taking on
        a true value if that column has a value of 1.  Please see the
        logictable module for details.  The function returns a table that
        start inclusive and end exclusive so truthtable(0, 2) will include
        row 0, but not row 2.

        INPUT:
            self -- the calling object.
            start -- an integer representing the row of the truth
                     table from which to start intilized to 0, which
                     is the first row when all the variables are
                     false.
            end -- an integer representing the last row of the truthtable
                   to be created.  It is initialized to the last row of the
                   full table.

        OUTPUT:
            Returns the truthtable (a 2-D array with the creating statement
            tacked on the front) corresponding to the statement.

        EXAMPLES:
        This example illustrates the creation of a statement.
            sage: import sage.logic.propcalc as propcalc
            sage: s = propcalc.formula("a&b|~(c|a)")
            sage: s.truthtable()
            a      b      c      value
            False  False  False  True
            False  False  True   False
            False  True   False  True
            False  True   True   False
            True   False  False  False
            True   False  True   False
            True   True   False  True
            True   True   True   True

        We can now create truthtable of rows 1 to 4, inclusive
            sage: s.truthtable(1, 5)
            a      b      c      value
            False  False  True   False
            False  True   False  True
            False  True   True   False
            True   False  False  False

        There should be no errors.

        NOTES:
            When sent with no start or end paramaters, this is an
            exponential time function requiring O(2**n) time, where
            n is the number of variables in the expression.
        """
        max = 2 ** len(self.__vars_order)
        if(end < 0):
            end = max
        if(end > max):
            end = max
        if(start < 0):
            start = 0
        if(start > max):
            start = max
        keys, table = [], []
        vars = {}
        for var in self.__vars_order:
            vars[var] = False
            keys.insert(0, var)
        keys = list(keys)
        for i in range(start, end):
            j = 0
            row = []
            for key in keys:
                bit = self.get_bit(i, j)
                vars[key] = bit
                j += 1
                row.insert(0, bit)
            row.append(booleval.eval_formula(self.__tree, vars))
            table.append(row)
        keys.reverse()
        table = logictable.Truthtable(table, keys)
        return table

    def evaluate(self, var_values):
        r"""
        Evaluates a formula for the given input values.

    	INPUT:
    	    self -- the calling object.
    	    var_values -- a dictionary containing pairs of
                              variables and their boolean values.
                              All variable must be present.

    	OUTPUT:
    	    Rerturn the evaluation of the formula with the given
            inputs, either True or False.

        EXAMPLES:
            sage: import sage.logic.propcalc as propcalc
            sage: f = propcalc.formula("a&b|c")
            sage: f.evaluate({'a':False, 'b':False, 'c':True})
            True
            sage: f.evaluate({'a':True, 'b':False, 'c':False})
            False
        """
        return booleval.eval_formula(self.__tree, var_values)

    def is_satisfiable(self):
        r"""
        Is_satisfiable determines if there is some assingment of
        variables for which the formula will be true.

        INPUT:
            self -- the calling object.

        OUTPUT:
            True if the formula can be satisfied, False otherwise.

        EXAMPLES:
            sage: import sage.logic.propcalc as propcalc
            sage: f = propcalc.formula("a|b")
            sage: f.is_satisfiable()
            True
            sage: g = f & (~f)
            sage: g.is_satisfiable()
            False
        """
        table = self.truthtable().get_table_list()
        for row in table[1:]:
            if row[-1] == True:
	   	return True
        return False

    def is_tautology(self):
        r"""
        Is_tautology determines if the formula is always
        true.

        INPUT:
            self -- the calling object.

        OUTPUT:
            True if the formula is a tautology, False otherwise.

        EXAMPLES:
            sage: import sage.logic.propcalc as propcalc
            sage: f = propcalc.formula("a|~a")
            sage: f.is_tautology()
            True
            sage: f = propcalc.formula("a&~a")
            sage: f.is_tautology()
            False
            sage: f = propcalc.formula("a&b")
            sage: f.is_tautology()
            False

        """
        return not (~self).is_satisfiable()

    def is_contradiction(self):
        r"""
        Is_contradiction determines if the formula is always
        false.

    	INPUT:
     	    self -- the calling object.

        OUTPUT:
            True if the formula is a contradiction, False otherwise.

    	EXAMPLES:
    	    sage: import sage.logic.propcalc as propcalc
            sage: f = propcalc.formula("a&~a")
            sage: f.is_contradiction()
            True
            sage: f = propcalc.formula("a|~a")
            sage: f.is_contradiction()
            False
            sage: f = propcalc.formula("a|b")
            sage: f.is_contradiction()
            False
        """
        return not self.is_satisfiable()

    def equivlant(self, other):
        r"""
        This function determines if two formulas are semantically
        equivlant.

        INPUT:
            self -- the calling object.
            other -- a boolformula instance.

        OUTPUT:
            True if the two fromulas are logically equivlant, False
            otherwise.

        EXAMPLES:
    	    sage: import sage.logic.propcalc as propcalc
            sage: f = propcalc.formula("(a|b)&c")
            sage: g = propcalc.formula("c&(a|b)")
            sage: f.equivlant(g)
            True
            sage: g = propcalc.formula("a|b&c")
            sage: f.equivlant(g)
            False
        """
        return self.iff(other).is_tautology()

    def convert_cnf(self):
        r"""
        This function calls convert_cnf_table().
        Please refer to that documentation.
        """
        self.convert_cnf_table()

    def convert_cnf_table(self):
        r"""
        This function converts an instance of boolformula to conjunctive normal form.
        It does this by examining the truthtable of the formula, and thus takes O(2^n)
        time, where n is the number of variables.

        INPUT:
             self -- the calling object.

        OUTPUT:
            An instance of boolformula with an identical truth table that is in
            conjunctive normal form.

        EXAMPLES:
            sage: import sage.logic.propcalc as propcalc
            sage: s = propcalc.formula("a^b<->c")
            sage: s.convert_cnf()
    	    sage: s
            (a|b|~c)&(a|~b|c)&(~a|b|c)&(~a|~b|~c)

        NOTES:
            This method creates the cnf parse tree by examining the logic
            table of the formula.  Creating the table requires O(2^n) time
            where n is the number of variables in the formula.
        """
        str = ''
        t = self.truthtable()
        table = t.get_table_list()
        vars = table[0]
        for row in table[1:]:
            if(row[-1] == False):
                str += '('
                for i in range(len(row) - 1):
                    if(row[i] == True):
                        str += '~'
                    str += vars[i]
                    str += '|'
                str = str[:-1] + ')&'
        self.__expression = str[:-1]
        #in case of tautology
        if(len(self.__expression) == 0):
            self.__expression = '(' + self.__vars_order[0] + '|~' + self.__vars_order[0] + ')'
        self.__tree, self.__vars_order = logicparser.parse(self.__expression)

    def convert_cnf_recur(self):
        r"""
        This function converts an instance of boolformula to conjunctive normal form.
        It does this by applying a set of rules that are gaurenteed to convert the
        formula.  Worst case the converted expression has an O(2^n) increase in
        size (and time as well), but if the formula is already in CNF (or close to)
        it is only O(n).

        INPUT:
            self -- the calling object.

        OUTPUT:
            An instance of boolformula with an identical truth table that is in
            conjunctive normal form.

        EXAMPLES:
            sage: import sage.logic.propcalc as propcalc
            sage: s = propcalc.formula("a^b<->c")
            sage: s.convert_cnf_recur()
            sage: s
            (~a|a|c)&(~b|a|c)&(~a|b|c)&(~b|b|c)&(~c|a|b)&(~c|~a|~b)

        NOTES:
            This function can require an exponential blow up in space from the
            original expression.  This in turn can require large amounts of time.
            Unless a formula is already in (or close to) being in cnf convert_cnf()
            is typically preferred, but results can vary.
        """
        self.__tree = logicparser.apply_func(self.__tree, self.reduce_op)
        self.__tree = logicparser.apply_func(self.__tree, self.dist_not)
        self.__tree = logicparser.apply_func(self.__tree, self.dist_ors)
        self.convert_expression()

    def satformat(self):
        r"""
        This function returns the satformat representation of a boolean formula.
        See www.cs.ubc.ca/~hoos/SATLIB/Benchmarks/SAT/satformat.ps for a
        description of satformat.

        INPUT:
             self -- the calling object.

        OUTPUT:
            A string representing the satformat represetation of this object.

        EXAMPLES:
            sage: import sage.logic.propcalc as propcalc
            sage: f = propcalc.formula("a&((b|c)^a->c)<->b")
            sage: f.convert_cnf()
            sage: f
            (a|~b|c)&(a|~b|~c)&(~a|b|~c)
            sage: f.satformat()
            'p cnf 3 0\n1 -2 3 0 1 -2 -3 \n0 -1 2 -3'

        NOTES:
            If the instance of boolean formula has not been converted to
            CNF form by a call to convert_cnf() or convert_cnf_recur()
            satformat() will call convert_cnf().  Please see the notes for
            convert_cnf() and convert_cnf_recur() for performance issues.
        """
        self.convert_cnf_table()
        s = ''
        vars_num = {}
        i = 0
        clauses = 0
        for e in self.__vars_order:
            vars_num[e] = str(i + 1)
            i += 1
        i = 0
        w = 1
        while i < len(self.__expression):
            c = self.__expression[i]
            if(c == ')'):
                clauses += 1
            if(c in '()|'):
                i += 1
                continue
            if(c == '~'):
                s += '-'
            elif(c == '&'):
                s += '0 '
            else:
                varname = ''
                while(i < self.__expression[i] not in '|) '):
                    varname += self.__expression[i]
                    i += 1
                s += vars_num[varname] + ' '
            if(len(s) >= (w * 15) and s[-1] != '-'):
                s += '\n'
                w += 1
            i += 1
        s = 'p cnf ' + str(len(self.__vars_order)) + ' ' + str(clauses) + '\n' + s
        return s[:-1]

    def simplify(self):
        r"""
        This function uses the propcalc package to simplify an expression to
        its minimal form.

        INPUT:
             self -- the calling object.

        OUTPUT:
            A simplified expression.

        EXAMPLES:
            sage: import sage.logic.propcalc as propcalc
            sage: f = propcalc.formula("a&((b|c)^a->c)<->b")
            sage: f.truthtable()
            a      b      c      value
            False  False  False  True
            False  False  True   True
            False  True   False  False
            False  True   True   False
            True   False  False  True
            True   False  True   False
            True   True   False  True
            True   True   True   True
            sage: f.simplify()
            (~a&~b)|(a&~b&~c)|(a&b)
            sage: f.truthtable()
            a      b      c      value
            False  False  False  True
            False  False  True   True
            False  True   False  False
            False  True   True   False
            True   False  False  True
            True   False  True   False
            True   True   False  True
            True   True   True   True

        NOTES:
            If the instance of boolean formula has not been converted to
            cnf form by a call to convert_cnf() or convert_cnf_recur()
            satformat() will call convert_cnf().  Please see the notes for
            convert_cnf() and convert_cnf_recur() for performance issues.
        """
        exp = ''
        self.__tree = logicparser.apply_func(self.__tree, self.reduce_op)
        plf = logicparser.apply_func(self.__tree, self.convert_opt)
        wff = boolopt.PLFtoWFF()(plf) # convert to positive-normal form
        wtd = boolopt.WFFtoDNF()
        dnf = wtd(wff)
        dnf = wtd.clean(dnf)
        if(dnf == [] or dnf == [[]]):
            exp = self.__vars_order[0] + '&~' + self.__vars_order[0] + ' '
        opt = boolopt.optimize(dnf)
        if(exp == '' and (opt == [] or opt == [[]])):
            exp = self.__vars_order[0] + '|~' + self.__vars_order[0] + ' '
        if(exp == ''):
            for con in opt:
                s = '('
                for prop in con:
                    if(prop[0] == 'notprop'):
                       s += '~'
                    s += prop[1] + '&'
                exp += s[:-1] + ')|'
        self.__expression = exp[:-1]
        self.__tree, self.__vars_order = logicparser.parse(self.__expression)
        return BooleanFormula(self.__expression, self.__tree, self.__vars_order)

    def convert_opt(self, tree):
        r"""
        This function can be applied to a parse tree to convert to the tuple
        form used by bool opt.  The expression must only contain '&', '|', and
        '~' operators.

        INPUT:
            self -- the calling object.
            tree -- a list of three elements corrospsponding to a branch of a
                    parse tree.
        OUTPUT:
            A tree branch that does not contain ^, ->, or <-> operators.

        EXAMPLES:
            sage: import sage.logic.propcalc as propcalc, sage.logic.logicparser as logicparser
            sage: s = propcalc.formula("a&(b|~c)")
            sage: tree = ['&', 'a', ['|', 'b', ['~', 'c', None]]]
            sage: logicparser.apply_func(tree, s.convert_opt)
            ('and', ('prop', 'a'), ('or', ('prop', 'b'), ('not', ('prop', 'c'))))
        """
        if(type(tree[1]) is not TupleType and tree[1] != None):
            lval = ('prop', tree[1])
        else:
            lval = tree[1]
        if(type(tree[2]) is not TupleType and tree[2] != None):
            rval = ('prop', tree[2])
        else:
            rval = tree[2]
        if(tree[0] == '~'):
            return ('not', lval)
        if(tree[0] == '&'):
            op = 'and'
        if(tree[0] == '|'):
            op = 'or'
        return (op, lval, rval)

    def add_statement(self, other, op):
        r"""
        This function takes two statements and combines them
        together with the given operator.

        INPUT:
            self -- the left hand side statement object.
            other -- the right hand side statement object.
            op -- the character of the operation to be performed.

        OUTPUT:
            Returns a new statement that is the first statement attached to
            the second statement by the operator op.

        EXAMPLES:
            sage: import sage.logic.propcalc as propcalc
            sage: s = propcalc.formula("a&b")
            sage: f = propcalc.formula("c^d")
            sage: s.add_statement(f, '|')
            (a&b)|(c^d)
            sage: s.add_statement(f, '->')
            (a&b)->(c^d)
        """
        exp = '(' + self.__expression + ')' + op + '(' + other.__expression + ')'
        parse_tree, vars_order = logicparser.parse(exp)
        return BooleanFormula(exp, parse_tree, vars_order)

    def get_bit(self, x, c):
        r"""
        This function returns bit c of the number x.  The 0 bit is the
        low order bit.  Errors should be handled gracefully by a return
        of false, and negative numbers x always return false while a
        negative c will index from the high order bit.

        INPUT:
            self -- the calling object.
            x -- an integer, the number from which to take the bit.
            c -- an integer, the bit number to be taken, where 0 is
                 the low order bit.

        OUTPUT:
            returns True if bit c of number x is 1, False otherwise.

        EXAMPLES:
            sage: import sage.logic.propcalc as propcalc
            sage: s = propcalc.formula("a&b")
            sage: s.get_bit(2, 1)
            True
            sage: s.get_bit(8, 0)
            False

        It is not an error to have a bit out of range.
            sage: s.get_bit(64, 7)
            False

        Nor is it an error to use a negative number.
            sage: s.get_bit(-1, 3)
            False
            sage: s.get_bit(64, -1)
            True
            sage: s.get_bit(64, -2)
            False
        """
        bits = []
        while(x > 0):
            if(x % 2 == 0):
                b = False
            else:
                b = True
            x = int(x / 2)
            bits.append(b)
        if(c > len(bits) - 1):
            return False
        else:
            return bits[c]

    def reduce_op(self, tree):
        r"""
        This function can be applied to a parse tree to convert if-and-only-if,
        if-then, and xor operations to operations only involving and/or operations.

        INPUT:
            self -- the calling object.
            tree -- a list of three elements corresponding to a branch of a
                    parse tree.
        OUTPUT:
            A tree branch that does not contain ^, ->, or <-> operators.

        EXAMPLES:
            sage: import sage.logic.propcalc as propcalc, sage.logic.logicparser as logicparser
            sage: s = propcalc.formula("a->b^c")
            sage: tree = ['->', 'a', ['^', 'b', 'c']]
            sage: logicparser.apply_func(tree, s.reduce_op)
            ['|', ['~', 'a', None], ['&', ['|', 'b', 'c'], ['~', ['&', 'b', 'c'], None]]]
        """
        if(tree[0] == '<->'):
            #parse tree for (~tree[1]|tree[2])&(~tree[2]|tree[1])
            new_tree = ['&', ['|', ['~', tree[1], None], tree[2]], \
                       ['|', ['~', tree[2], None], tree[1]]]
        elif(tree[0] == '^'):
            #parse tree for (tree[1]|tree[2])&~(tree[1]&tree[2])
            new_tree = ['&', ['|', tree[1], tree[2]], \
                       ['~', ['&', tree[1], tree[2]], None]]
        elif(tree[0] == '->'):
            #parse tree for ~tree[1]|tree[2]
            new_tree = ['|', ['~', tree[1], None], tree[2]]
        else:
            new_tree = tree
        return new_tree

    def dist_not(self, tree):
        r"""
        This function can be applied to a parse tree to distribute not operators
        over other operators.

        INPUT:
            self -- the calling object.
            tree -- a list of three elements corrospsponding to a branch of a
                    parse tree.
        OUTPUT:
            A tree branch that does not contain un-distributed nots.

        EXAMPLES:
            sage: import sage.logic.propcalc as propcalc, sage.logic.logicparser as logicparser
            sage: s = propcalc.formula("~(a&b)")
            sage: tree = ['~', ['&', 'a', 'b'], None]
            sage: logicparser.apply_func(tree, s.dist_not) #long time
            ['|', ['~', 'a', None], ['~', 'b', None]]
        """
        if(tree[0] == '~' and type(tree[1]) is ListType):
            op = tree[1][0]
            if(op != '~'):
                if(op == '&'):
                    op = '|'
                else:
                    op = '&'
                new_tree = [op, ['~', tree[1][1], None], ['~', tree[1][2], None]]
                return logicparser.apply_func(new_tree, self.dist_not)
            else:
                #cancel double negative
                return tree[1][1]
        else:
            return tree

    def dist_ors(self, tree):
        r"""
        This function can be applied to a parse tree to distribute or over and.

        INPUT:
            self -- the calling object.
            tree -- a list of three elements corrospsponding to a branch of a
                    parse tree.
        OUTPUT:
            A tree branch that does not contain un-distributed ors.

        EXAMPLES:
            sage: import sage.logic.propcalc as propcalc, sage.logic.logicparser as logicparser
            sage: s = propcalc.formula("(a&b)|(a&c)")
            sage: tree = ['|', ['&', 'a', 'b'], ['&', 'a', 'c']]
            sage: logicparser.apply_func(tree, s.dist_ors) #long time
            ['&', ['&', ['|', 'a', 'a'], ['|', 'b', 'a']], ['&', ['|', 'a', 'c'], ['|', 'b', 'c']]]

        """
        if(tree[0] == '|' and type(tree[2]) is ListType and tree[2][0] == '&'):
            new_tree = ['&', ['|', tree[1], tree[2][1]], ['|', tree[1], tree[2][2]]]
            return logicparser.apply_func(new_tree, self.dist_ors)
        if(tree[0] == '|' and type(tree[1]) is ListType and tree[1][0] == '&'):
           new_tree = ['&', ['|', tree[1][1], tree[2]], ['|', tree[1][2], tree[2]]]
           return logicparser.apply_func(new_tree, self.dist_ors)
        return tree

    def to_infix(self, tree):
        r"""
        This function can be applied to a parse tree to convert it from prefix to
        infix.

        INPUT:
            self -- the calling object.
            tree -- a list of three elements corrospsponding to a branch of a
                    parse tree.

        OUTPUT:
            A tree branch in infix form.

        EXAMPLES:
            sage: import sage.logic.propcalc as propcalc, sage.logic.logicparser as logicparser
            sage: s = propcalc.formula("(a&b)|(a&c)")
            sage: tree = ['|', ['&', 'a', 'b'], ['&', 'a', 'c']]
            sage: logicparser.apply_func(tree, s.to_infix)
            [['a', '&', 'b'], '|', ['a', '&', 'c']]
        """
        if(tree[0] != '~'):
            return [tree[1], tree[0], tree[2]]
        return tree

    def convert_expression(self):
        r"""
        This function converts the string expression associated with an instance
        of boolformula to match with its tree reprsentation after being converted
        to conjunctive normal form.

        INPUT:
            self -- the calling object.

        OUTPUT:
            None.

        EXAMPLES:
            sage: import sage.logic.propcalc as propcalc
            sage: s = propcalc.formula("a^b<->c")
            sage: s.convert_cnf_recur(); s  #long time
            (~a|a|c)&(~b|a|c)&(~a|b|c)&(~b|b|c)&(~c|a|b)&(~c|~a|~b)

        """
        ttree = self.__tree[:]
        ttree = logicparser.apply_func(ttree, self.to_infix)
        self.__expression = ''
        str_tree = str(ttree)
        open_flag = False
        put_flag = False
        i = 0
        for c in str_tree:
            if(i < len(str_tree) - 1):
                op = self.get_next_op(str_tree[i:])
            if(op == '|' and not open_flag):
                self.__expression += '('
                open_flag = True
            if(i < len(str_tree) - 2 and str_tree[i + 1] == '&' and open_flag):
                open_flag = False
                self.__expression += ')'
            if(str_tree[i:i + 4] == 'None'):
                i += 4
            if(i < len(str_tree) and str_tree[i] not in ' \',[]'):
                self.__expression += str_tree[i]
            i += 1
        if(open_flag == True):
            self.__expression += ')'

    def get_next_op(self, str):
        r"""
        This function returns the next operator in a string.

        INPUT:
            self -- the calling object.
            tree -- a string containing a logical expression.

        OUTPUT:
            The next operator in the string.

        EXAMPLES:
            sage: import sage.logic.propcalc as propcalc
            sage: s = propcalc.formula("f&p")
            sage: s.get_next_op("abra|cadabra")
            '|'
        """
        i = 0
        while(i < len(str) - 1 and str[i] != '&' and str[i] != '|'):
            i += 1
        return str[i]
