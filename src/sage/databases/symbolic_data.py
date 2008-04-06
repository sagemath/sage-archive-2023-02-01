"""nodoctest
Very thin wrapper for the SymbolicData benchmark and test ideals as
published on http://www.symbolicdata.org . From that website:

    For different purposes algorithms and implementations are tested on
    certified and reliable data. The development of tools and data for
    such tests is usually 'orthogonal' to the main implementational
    efforts, it requires different skills and technologies and is not
    loved by programmers. On the other hand, in many cases tools and data
    could easily be reused - with slight modifications - across similar
    projects. The SymbolicData Project is set out to coordinate such
    efforts within the Computer Algebra Community.  Commonly collected
    certified and reliable data can also be used to compare otherwise
    incomparable approaches, algorithms, and implementations. Benchmark
    suites and Challenges for symbolic computations are not as well
    established as in other areas of computer science. This is probably
    due to the fact that there are not yet well agreed aims of such a
    benchmarking. Nevertheless various (often high quality) special
    benchmarks are scattered through the literature.  During the last
    years efforts toward collection of test data for symbolic
    computations were intensified. They focused mainly on the creation of
    general benchmarks for different areas of symbolic computation and
    the collection of such activities on different Web site.  For further
    qualification of these efforts it would be of great benefit to create
    a commonly available digital archive of these special benchmark data
    scattered through the literature. This would provide the community
    with an electronic repository of certified data that could be
    addressed and extended during further development.

EXAMPLE:
    sage: sd = SymbolicData(); sd
    SymbolicData with 372 ideals

    sage: sd.ZeroDim__example_1
    Ideal (-10 + x2^2 + x1^2, -16 + 2*x2^2 + x1*x2 + x1^2) of Polynomial Ring in x1, x2 over Rational Field

    sage: sd.Katsura_3
    Ideal (-1 + 2*u3 + 2*u2 + 2*u1 + u0, -1*u2 + 2*u1*u3 + u1^2 + 2*u0*u2, 2*u2*u3 - u1 + 2*u1*u2 + 2*u0*u1, 2*u3^2 + 2*u2^2 + 2*u1^2 - u0 + u0^2) of Polynomial Ring in u0, u1, u2, u3 over Rational Field


    sage: sd.get_ideal('Katsura_3',GF(127),'degrevlex')
    Ideal (126 + 2*u3 + 2*u2 + 2*u1 + u0, 126*u2 + 2*u1*u3 + u1^2 + 2*u0*u2, 2*u2*u3 + 126*u1 + 2*u1*u2 + 2*u0*u1, 2*u3^2 + 2*u2^2 + 2*u1^2 + 126*u0 + u0^2) of Polynomial Ring in u0, u1, u2, u3 over Finite Field of size 127

AUTHOR: Martin Albrecht <malb@informatik.uni-bremen.de>

"""
import os
from xml.dom.minidom import parse
from sage.rings.rational_field import QQ
from sage.rings.ideal import Ideal
from sage.rings.polynomial.polynomial_ring_constructor import PolynomialRing


class SymbolicData:
    """
    Access to the benchmark and test ideals of the SymbolicData
    suite. This class needs the optional database-symbolicdata package
    to be installed.
    """
    def __init__(self):

        path=os.environ["SAGE_ROOT"]+"/data/symbolic_data"
        self.__intpath = path + "/Data/XMLResources/INTPS/"
        self.__genpath = path + "/Data/XMLResources/GenPS/"


    def get_ideal(self, name, base_ring=QQ, term_order="degrevlex"):
        """
        Returns the ideal given by 'name' over the base ring given by
        'base_ring' in a polynomial ring with the term order given by
        'term_order'.

        INPUT:
            name -- name as in the SymbolicData package
            base_ring -- base ring for the polynomial ring (default: QQ)
            term_order -- term order for the polynomial ring (default: degrevlex)

        OUTPUT:
            ideal as given by name in MPolynomialRing(base_ring,vars,term_order)

        """

        def _getTextFromNode(node):
            t = ""
            for n in node.childNodes:
                if n.nodeType == n.TEXT_NODE:
                    t += str(n.nodeValue)
                else:
                    raise NotTextNodeError
            return t

        def _dom2ideal(node):
            """
            """
            l = []

            if str(node.nodeName) in ['vars','poly']:
                l.append(_getTextFromNode(node))

            for c in node.childNodes:
                l += _dom2ideal(c)

            return l

        name = name.replace('__','.')

        try:
            name = self.__intpath + name + ".xml"
            open(name)
        except IOError:
            try:
                name = self.__genpath + name + ".xml"
                open(name)
            except IOError:
                raise AttributeError, "Ideal not found on disk"


        dom = parse(name)
        res = _dom2ideal(dom)
        vars,polys = res[0].replace("_",""),[p.replace("_","") for p in res[1:]]

        return Ideal(PolynomialRing(base_ring, len(vars.split(",")), vars), polys)


    def __repr__(self):
        try:
            l = len(self.trait_names())
        except AttributeError:
            l = 0
        return "SymbolicData with %d ideals"%l

    def __getattr__(self, name):
        return self.get_ideal(name)

    def trait_names(self):
        """
        """
        if hasattr(self,"__ideals"): return self.__ideals
        try:
            __ideals = [s.replace('.xml','') for s in  os.listdir(self.__intpath)]
            __ideals += [s.replace('.xml','') for s in  os.listdir(self.__genpath)]
            self.__ideals = [s.replace('.', '__') for s in __ideals]
            return self.__ideals
        except OSError:
            raise AttributeError, "Could not find symbolic data, you should perhaps install the optional package"
