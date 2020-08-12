from . import BooleSet, interpolate_smallest_lex


class PartialFunction(object):
    """docstring for PartialFunction"""

    def __init__(self, zeros, ones):
        super(PartialFunction, self).__init__()
        self.zeros = zeros.set()
        self.ones = ones.set()

    def interpolate_smallest_lex(self):
        return interpolate_smallest_lex(self.zeros, self.ones)

    def __str__(self):
        return "PartialFunction(zeros=" + str(self.zeros) + ", ones=" + str(
            self.ones) + ")"

    def definedOn(self):
        return self.zeros.union(self.ones)

    def __add__(self, other):
        domain = self.definedOn().intersect(other.definedOn())
        zeros = self.zeros.intersect(other.zeros).union(self.ones.intersect(
            other.ones))
        ones = self.zeros.intersect(other.ones).union(self.ones.intersect(
            other.zeros))
        assert zeros.diff(domain).empty()
        assert ones.diff(domain).empty()
        return PartialFunction(zeros, ones)

    def __repr__(self):
        return str(self)

    def __mul__(self, other):
        zeros = self.zeros.union(other.zeros)
        ones = self.ones.intersect(other.ones)
        return PartialFunction(zeros, ones)

    def __or__(self, other):
        zeros = self.zeros.intersect(other.zeros)
        ones = self.ones.union(other.ones)
        return PartialFunction(zeros, ones)

    def __xor__(self, other):
        return self + other

    def __and__(self, other):
        return self * other
