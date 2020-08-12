# -*- python -*-
# encoding: utf-8
"""
cluster.py

Created by Michael Brickenstein on 2011-08-05.
Copyright 2011 The PolyBoRi Team. See LICENSE file.
"""

import sys
import os
from .statistics import used_vars, used_vars_set
from .PyPolyBoRi import Variable


def main():
    pass


class ClusterAlgorithmFailed(Exception):
    pass


class ClusterAlgorithm(object):
    def __init__(self, ideal, determination_modifier=1):
        if len(ideal) == 0:
            raise ValueError('ideal generators list should be non empty')

        self.ideal = ideal

        self.determination_modifier = determination_modifier

        self.used_variables_ideal = used_vars_set(ideal)

        self.number_of_used_variables_in_ideal = len(self.used_variables_ideal)
        self.used_variables_of_polynomial = dict(
            [(p, set(p.vars_as_monomial().variables())) for p in ideal])
        self.variables_introduction_mapping = dict()

        self.cluster = set()
        self.used_variables_cluster = set()
        self.build_variables_usage()
        self.initialize_variables_introduction_mapping()

    def build_variables_usage(self):
        self.variables_usage = dict()
        for (p, variables) in self.used_variables_of_polynomial.items():
            for v in variables:
                self.variables_usage.setdefault(v, []).append(p)

    def initialize_variables_introduction_mapping(self):
        self._build_variables_introduction_mapping()

    def _build_variables_introduction_mapping(self):
        def var_set_to_tuple(s):
            return tuple(sorted(s, key=Variable.index))
        self.variables_introduction_mapping.clear()
        for (p, var_set) in self.used_variables_of_polynomial.items():
            if p in self.cluster:
                continue
            as_tuple = var_set_to_tuple(var_set.difference(self.
                used_variables_cluster))
            self.variables_introduction_mapping.setdefault(
                as_tuple, []).append(p)

    def adjust_variables_introduction_mapping(self, introduced_variables):
        self._build_variables_introduction_mapping()

    def determined_enough(self):
        return (self.number_of_used_variables_in_cluster + self.
            determination_modifier <= len(self.cluster))

    def find_cluster(self):
        p = self.initial_choice()
        self.cluster = set()
        self.add_polynomial_to_cluster(p)
        while not self.determined_enough():
            self.increase_cluster()
        return list(self.cluster)

    def add_polynomial_to_cluster(self, p):
        self.cluster.add(p)
        self.used_variables_cluster = set(used_vars_set(self.cluster).
            variables())
        self.number_of_used_variables_in_cluster = len(self.
            used_variables_cluster)
        self.adjust_variables_introduction_mapping(self.
            used_variables_of_polynomial[p])

    def initial_choice(self):
        def max_key(entry):
            (entry_variable, entry_polynomials) = entry
            return len(entry_polynomials)
        (variable, polynomials) = max(self.variables_usage.items(),
            key=max_key)

        def min_key(p):
            return len(self.used_variables_of_polynomial[p])
        return min(polynomials, key=min_key)

    def increase_cluster(self):
        introduced_variables_possibilities = (list(self.
            variables_introduction_mapping.keys()))
        introduced_variables = min(introduced_variables_possibilities, key=len)
        polynomials = self.variables_introduction_mapping[introduced_variables]
        assert len(polynomials) > 0
        for p in polynomials:
            self.add_polynomial_to_cluster(p)
        if len(self.cluster) == len(self.ideal):
            raise ClusterAlgorithmFailed
        self.adjust_variables_introduction_mapping(introduced_variables)


def find_cluster(ideal):
    algorithm = ClusterAlgorithm(ideal)
    return algorithm.find_cluster()


if __name__ == '__main__':
    main()
