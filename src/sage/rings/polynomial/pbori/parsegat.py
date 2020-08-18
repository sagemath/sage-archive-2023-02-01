from __future__ import print_function
#import pathadjuster

if __name__ == '__main__':
    from sys import path as search_path
    from os import path as file_path
    search_path.append(file_path.join(file_path.dirname(__file__), '..'))


def _exists():
    """PolyBoRi convention: checking optional components for prerequisites here

    sage: _exists()
    True
    """
    try:
        import pyparsing
    except ImportError:
        return False

    return True

#from .gbrefs import declare_ring
from . import *
from .ll import ll_encode
from optparse import OptionParser
from .statistics import used_vars_set
from itertools import chain
try:
    from cStringIO import StringIO
except ImportError:
    from io import StringIO
from sys import stderr
parser = OptionParser()

gat_max = 20000
next_max = 2000
output_max = 2000
inter_max = 20000
input_max = 2000
state_max = 2000

# declare_ring([Block("x",gat_max),Block("xnext",next_max,reverse =
# True),Block("xoutput",output_max,reverse =
# True),Block("xinter",inter_max,reverse =
# True),Block("xinput",input_max,reverse =
# True),Block("xstate",state_max,reverse = True)],globals()) input order is
# antitopological we reverse, when mapping to final variables
parser.add_option("--forward-iteration",
                  action="store_true", dest="forward", default=True,
                  help="switch between forward and backward iteration")

parser.add_option("--backward-iteration",
                  action="store_false", dest="forward", default=True,
                  help="switch between forward and backward iteration")

parser.add_option("--include-outputs",
                  action="store_true", dest="include_outputs", default=True,
                  help="include circuit outputs")

parser.add_option("--initialize-state", action="store", dest="initialize",
    default="noinit", type="choice", choices=["noinit", "random", "zero"],
    help="initialize state variables/only useful with forward iteration")


def parse_identifier(str, log, tokens):
    return int(tokens[0])


def throw_away(str, log, tokens):
    return "None"

assigned = set()
from random import Random


def add_negated(str, log, tokens):
    p = tokens[1]
    l = p.lead()
    if not l in assigned and r.randint(0, 200) == 0:
        assigned.add(l)
        print("NEG", p + 1)
        return p + 1
    else:
        return "NONE"


class FamiliarityException(Exception):
    """docstring for FamiliarityException"""

    def __init__(self):
        super(FamiliarityException, self).__init__()


def fix_symbol_name(str, log, tokens):
    return (tokens[0].replace("-", "_").replace("]", "").replace("[", "").
        replace("/", "_"))


class DeterminingEquation(object):
    """docstring for DeterminingEquation"""

    def __init__(self, variable, mapped_to):
        super(DeterminingEquation, self).__init__()
        self.variable = variable
        self.mapped_to = mapped_to

    def _get_equation(self):
        return self.variable + self.mapped_to
    equation = property(_get_equation)


# be careful: for next state/output we directly generate mapped_to + value
# instead of introducing a variable and mapping it later
class VariableManager(object):
    """docstring for VariableManager"""

    def __init__(self, ring, prefix="", initialize="noinit", **kwd):
        super(VariableManager, self).__init__()
        self.ring = ring
        self.output = []
        self.input = []
        self.state = []
        self.next = []
        self.deletion_candidates = []
        self.intermediate = []

        self.next_state_equations = []
        self.intermediate_equations = []
        self.__map = dict()
        self.prefix = prefix
        self.initialize = initialize
        self.ands = dict()
        self.xors = dict()
        self.tails = []
        for (k, v) in kwd.items():
            setattr(self, k, v)

    def gauss(self, determining_equations):
        l = []
        res = []
        output_set = set(self.output)
        for p in determining_equations:
            eq = p.equation
            if not p.variable in output_set and eq.deg() == 1:
                assert eq.lead() == p.variable
                l.append(eq)
                self.deletion_candidates.append(p.variable)
            else:
                res.append(p)
        l = gauss_on_polys(l)
        for p in l:
            res.append(DeterminingEquation(p.lead(), p + p.lead()))
        return res

    def ideals(self):
        def sort_key(p):
            return p.navigation().value()

        def my_sort(l):
            return sorted(l, key=sort_key)
        #self.intermediate_equations = self.gauss(self.intermediate_equations)
        tail_variables = set(used_vars_set([e.mapped_to for e in chain(self.
            next_state_equations, self.intermediate_equations)]).variables())
        to_delete = []
        for v in self.deletion_candidates:
            if not v in tail_variables:
                to_delete.append(v)
        to_delete = set(to_delete)
        self.intermediate_equations = [e for e in self.intermediate_equations
            if not e.variable in to_delete]

        # we don't delete the variable itself, as we will confuse
        # gluemultipliers,that looks on the lengths of self.intermediate...
        def equations(determining_equations):
            return [e.equation for e in determining_equations]
        ideal_state = []
        ideal_next_state = equations(self.next_state_equations)
        ideal_intermediate = equations(self.intermediate_equations)

        if self.initialize != "noinit":
            if self.initialize == "zero":
                initializer = zero_fun
            else:
                initializer = random_bit
            for v in variables:
                if str(v)[:1] == "s":
                    ideal_state.append(v + initializer())

        return [
            my_sort(self.apply_map(i)) for i in
            (ideal_state, ideal_intermediate, ideal_next_state)]

    def set_variable_name(self, i, name):
        _set_variable_name(self.ring, i, self.prefix + name)

    def parse_output_action(self, str, log, tokens):
        p = Polynomial(tokens[1])
        #self.output.append(p)
        #mapped_to = xoutput(len(self.output)-1)
        mapped_to = self.xoutput(len(self.output))
        self.output.append(mapped_to)
        self.set_variable_name(mapped_to.index(), "out_" + tokens[2])
        self.tails.append(p)
        if self.include_outputs:
            self.intermediate_equations.append(DeterminingEquation(mapped_to,
                p))
            return mapped_to + p
        else:
            return "NONE"

    def parse_buffer_action(self, str, log, tokens):
        return "NONE"

    def parse_next_state_action(self, str, log, tokens):
        p = Polynomial(tokens[1])
        #self.next.append(p)
        #mapped_to = xnext(len(self.next)-1)
        mapped_to = self.xnext(len(self.next))
        self.next.append(mapped_to)
        self.set_variable_name(mapped_to.index(), "nstate_" + tokens[2])
        self.next_state_equations.append(DeterminingEquation(mapped_to, p))
        self.tails.append(p)
        return mapped_to + p

    def parse_state_action(self, str, log, tokens):
        p = self.x(tokens[0])
        #self.state.append(p)
        #mapped_to = xstate(len(self.state)-1)
        mapped_to = self.xstate(len(self.state))
        self.state.append(mapped_to)
        self.set_variable_name(mapped_to.index(), "state_" + tokens[2])
        self.map(p, mapped_to)
        return "NONE"

    def parse_input_action(self, str, log, tokens):
        p = self.x(tokens[0])
        #self.input.append(p)
        #mapped_to = xinput(len(self.input)-1)
        mapped_to = self.xinput(len(self.input))
        self.input.append(mapped_to)
        self.set_variable_name(mapped_to.index(), "in_" + tokens[2])
        self.map(p, mapped_to)
        return "NONE"

    def evaluates_to_variables(self, signal):
        ands = self.ands[signal]
        return list(
            used_vars_set(ands).variables())

    def eval_and_gat(self, lit, inputs=None):
        if inputs is None:
            inputs = self.ands[lit.lex_lead()]
        assert len(inputs) > 0
        res = inputs[0]
        for i in inputs[1:]:
            res = res * i
        if lit.has_constant_part():
            res = res + 1
        return res

    def simplified_product(self, out, op1, op2):

        if op1.deg() == 1 and op2.deg() == 1:
            op1v = op1.lex_lead()
            op2v = op2.lex_lead()

        if op1v not in self.xors \
            and op2v not in self.xors \
            and op1v in self.ands \
            and op2v in self.ands:
            mapped_to = [p for p in chain(self.ands[op1v], self.ands[op2v])]
            for p in mapped_to:
                if p.deg() != 1:
                    return op1 * op2
            mapped_to = set([p.lex_lead() for p in mapped_to])
            #print >>stderr, op1, op2, mapped_to
            if len(mapped_to) <= 2:
                #xor pattern
                self.deletion_candidates.extend([op1v, op2v])
                assert op1.vars_as_monomial().deg() == 1
                assert op2.vars_as_monomial().deg() == 1
                op1 = self.eval_and_gat(op1)
                op2 = self.eval_and_gat(op2)
                outer_res = self.eval_and_gat(out, [op1, op2])
                self.xors[out] = outer_res
            else:
                try:
                    for (op1, op1v, op2, op2v) in \
                        [(op1, op1v, op2, op2v), (op2, op2v, op1, op1v)]:
                        (op1e, op2e) = [
                            self.evaluates_to_variables(left_parent)
                            for left_parent in (op1v, op2v)]
                        for (i, s) in enumerate(op2e):
                            if s in self.xors or s not in self.ands:
                                continue

                            if self.evaluates_to_variables(s) == op1e:
                                raise FamiliarityException()
                except FamiliarityException:
                    op1 = self.eval_and_gat(op1)
                    op2_inputs = list(self.ands[op2v])
                    for (i, s2) in enumerate(op2_inputs):
                        if s2.lex_lead() == s:
                            op2_inputs[i] = self.eval_and_gat(s2)
                    op2 = self.eval_and_gat(op2, inputs=op2_inputs)
                    self.deletion_candidates.extend((s, op1v, op2v))
                else:
                    # pattern
                    # 34 AND ~16 ~15
                    # 35 AND 14 34
                    # 36 AND 16 15
                    # 37 AND ~14 36
                    # 38 AND ~37 ~35

                    try:
                        (op1a, op2a) = [self.ands.get(v)
                            for v in (op1v, op2v)]
                        assert len(op1a) == 2
                        assert len(op2a) == 2
                        print("+", file=stderr)
                        for op1a in [[op1a[0], op1a[1]], op1a]:
                            for op2a in [[op2a[0], op2a[1]], op2a]:
                                print(op1a[0], op2a[0], file=stderr)
                                if op1a[0] == op2a[0] + 1:
                                    print("-", file=stderr)
                                    if op1a[1] in self.ands and \
                                        op2a[1] in self.ands and \
                                        op1a[1] not in self.xors \
                                        and op2a[1] not in self.xors:
                                        if set(self.ands[op1a[1]]) == set([1
                                            + s for s in self.ands[op2a[1]]]):
                                            raise FamiliarityException
                    except FamiliarityException:
                        self.deletion_candidates.extend([op1v, op2v, op2a[1].
                            lex_lead(), op1a[1].lex_lead()])
                        for group in (op2a, op1a):
                            group[1] = self.eval_and_gat(group[1])
                        (op1, op2) = [
                            self.eval_and_gat(op, inputs) for (op, inputs)
                            in [(op1, op1a), (op2, op2a)]]

        return op1 * op2

    def parse_and_action(self, string, log, tokens):

        inter = self.x(tokens[0])
        op1 = tokens[2]
        op2 = tokens[3]
        #self.intermediate.append(inter)
        #mapped_to = xinter(len(self.intermediate)-1)
        mapped_to = self.xinter(len(self.intermediate))
        self.ands[inter] = (op1, op2)
        self.intermediate.append(mapped_to)
        self.set_variable_name(mapped_to.index(), "y" + str(tokens[0]))
        self.map(inter, mapped_to)
        tail = self.simplified_product(inter, op1, op2)
        self.tails.append(tail)
        eq = inter + tail
        self.intermediate_equations.append(DeterminingEquation(inter, tail))
        return eq

    def map(self, from_, to):
        self.__map[from_] = to

    def apply_map(self, eq):
        encoded = ll_encode((k + v for (k, v) in self.__map.items()))
        return [ll_red_nf_noredsb(p, encoded) for p in eq]

    def parse_sign(self, str, log, tokens):
        if list(tokens) != [0]:
            return [1]
        else:
            return tokens

    def parse_idenfifier_ref(self, str, log, tokens):
        if tokens[1] in (0, 1):
            #from sys import stderr
            #stderr.write("TOKENS:"+repr(tokens))
            #stderr.flush()
            return Polynomial(tokens[0] + tokens[1])
        return tokens[0] + self.x(tokens[1])


def generate_gat_bnf(manager):

    from pyparsing import (Literal, CaselessLiteral, Word, Combine, Group,
        Optional, ZeroOrMore, Forward, nums, alphas, Or, restOfLine,
        OneOrMore, restOfLine, alphanums)

    identifier = Word(nums).setParseAction(parse_identifier)
    symbolic_name = (Word(alphanums + "_-[]/", alphanums + "_-[]/")). \
        setParseAction(fix_symbol_name)
    meta = ZeroOrMore(
    Or(
    [
    CaselessLiteral("WRITER"),
    CaselessLiteral("DESIGN"),
    CaselessLiteral("ORIGIN")]) + restOfLine)
    end = CaselessLiteral("END")
    and_ = CaselessLiteral("AND")
    input_ = CaselessLiteral("INPUT")
    state = CaselessLiteral("STATE")
    nstate = CaselessLiteral("NSTATE")
    output = CaselessLiteral("OUTPUT")
    buffer_ = CaselessLiteral("BUFFER")
    identifier_ref = \
        (Optional(Literal("~"), default=0).setParseAction(
        manager.parse_sign) + identifier).setParseAction(
            manager.parse_idenfifier_ref)

    input_line = (identifier + input_ + symbolic_name).setParseAction(manager.
        parse_input_action)
    and_line = (identifier + and_ + identifier_ref + identifier_ref). \
        setParseAction(manager.parse_and_action)
    buffer_line = (buffer_ + identifier_ref + symbolic_name).setParseAction(
        manager.parse_buffer_action)
    output_line = (output + identifier_ref + symbolic_name).setParseAction(
        manager.parse_output_action)
    state_line = (identifier + state + symbolic_name).setParseAction(manager. \
        parse_state_action)
    nstate_line = (nstate + identifier_ref + symbolic_name).setParseAction(
        manager.parse_next_state_action)
    assignment = Or([output_line, and_line, input_line, state_line,
        nstate_line, buffer_line])

    gat = meta + OneOrMore(assignment) + end
    return gat
generator = Random(123)


def parse(f, manager):
    f = open(f)
    content = f.read()
    bnf = generate_gat_bnf(manager)
    parsed = bnf.parseString(content)
    f.close()


def zero_fun():
    return 0


def random_bit():
    return generator.randint(0, 1)


def format_grouped(l, group_size=10, indent=0):
    s = StringIO()
    s.write("[")
    last = len(l) - 1
    for (i, e) in enumerate(l):
        if i % group_size == 0:
            s.write("\n" + indent * " ")
        s.write(e)
        if i != last:
            s.write(", ")
    s.write("]")
    return s.getvalue()


def generate_three_ideal_output(ideal_state, ideal_intermediate,
    ideal_next_state, variables):
    print("declare_ring(" + format_grouped([repr(str(v)) for v in variables],
        indent=4) + ")")
    print("block_start_hints=" + repr(block_starts))
    print("ideal_intermediate=[")
    print(",\n".join((str(p) for p in ideal_intermediate)))
    print("]")
    print("ideal_state=[")
    print(",\n".join((str(p) for p in ideal_state)))
    print("]")
    print("ideal_next_state=[")
    print(",\n".join((str(p) for p in ideal_next_state)))
    print("]")
    print("ideal = ideal_state+ideal_next_state+ideal_intermediate")

if __name__ == '__main__':
    (options, args) = parser.parse_args()
    kwd_args = dict()
    ring = declare_ring([
        "t",
        Block("x", gat_max, reverse=True),
        Block("xnext", next_max, reverse=True),
        Block("xoutput", output_max, reverse=True),
        Block("xinter", inter_max, reverse=True),
        Block("xinput", input_max, reverse=True),
        Block("xstate", state_max, reverse=True)], kwd_args)

    manager = VariableManager(ring=ring, include_outputs=options.
        include_outputs, initialize=options.initialize, **kwd_args)

    from sys import argv
    f = args[0]

    parse(f, manager)
    (ideal_state, ideal_intermediate, ideal_next_state) = manager.ideals()
    ideal = ideal_state + ideal_intermediate + ideal_next_state

    variables = []
    used_vars = set(used_vars_set(ideal).variables())
    if not options.forward:
        variables = list(chain(reversed(manager.next), reversed(manager.output
            ), reversed(manager.intermediate), reversed(manager.input),
            reversed(manager.state)))
    else:
        variables = list(
            chain(reversed(manager.output),
            reversed(manager.intermediate),
            reversed(manager.input),
            reversed(manager.state),
            reversed(manager.next)))
    variables = ["t"] + variables
    beginnings = [str(v)[:1] for v in variables]
    block_starts = []
    last = beginnings[0]

    for (i, s) in enumerate(beginnings):
        if s != last:
            block_starts.append(i)
        last = s

    generate_three_ideal_output(ideal_state, ideal_intermediate,
        ideal_next_state, variables)
