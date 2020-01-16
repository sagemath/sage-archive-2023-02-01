# -*- coding: utf-8 -*-
r"""
FindStat - the Combinatorial Statistic Finder.

The FindStat database can be found at::

    sage: findstat()
    The Combinatorial Statistic Finder (https://www.findstat.org/)

Fix the following three notions:

- A *combinatorial collection* is a set `S` with interesting combinatorial properties,
- a *combinatorial map* is a combinatorially interesting map `f: S \to S'` between combinatorial collections, and
- a *combinatorial statistic* is a combinatorially interesting map `s: S \to \ZZ`.

You can use the sage interface to FindStat to:

- identify a combinatorial statistic from the values on a few small objects,
- obtain more terms, formulae, references, etc. for a given statistic,
- edit statistics and maps and submit new statistics.

AUTHORS:

- Martin Rubey (2015): initial version.
- Martin Rubey (2020): rewrite, adapt to new FindStat API

A guided tour
-------------

Retrieving information
^^^^^^^^^^^^^^^^^^^^^^

The most straightforward application of the FindStat interface is to
gather information about a combinatorial statistic.  To do this, we
supply :class:`findstat<FindStat>` with a list of ``(object, value)``
pairs.  For example::

    sage: PM = PerfectMatchings
    sage: r = findstat([(m, m.number_of_nestings()) for n in range(6) for m in PM(2*n)]); r    # optional -- internet
    0: St000042oMp00116 (quality [100, 100])
    1: St000041 (quality [20, 100])
    ...

The result of this query is a list (presented as a
:class:`sage.databases.oeis.FancyTuple`) of matches.  Each match
consists of a :class:`FindStatCompoundStatistic` `s: S \to \ZZ` and
an indication of the quality of the match.

The precise meaning of the result is as follows:

    The composition `f_n \circ ... \circ f_2 \circ f_1` applied to
    the objects sent to FindStat agrees with all ``(object, value)``
    pairs of `s` in the database.

    Suppose that the quality of the match is `(q_a, q_d)`.  Then
    `q_a` is the percentage of ``(object, value)`` pairs that are in
    the database among those which were sent to FindStat, and `q_d`
    is the percentage of ``(object, value)`` pairs with distinct
    values in the database among those which were sent to FindStat.

Put differently, if ``quality`` is not too small it is likely that
the statistic sent to FindStat equals `s \circ f_n \circ ... \circ
f_2 \circ f_1`.  If `q_a` is large, but `q_b` is small, then there
were many matches, but while the sought for statistic attains many
distinct values, the match found by FindStat covers only ``(object,
value)`` pairs for few values.

In the case at hand, for the match ``St000041``, the list of maps is
empty.  We can retrieve the description of the statistic from the
database as follows::

    sage: print(r[1].statistic().description())                                 # optional -- internet
    The number of nestings of a perfect matching.
    <BLANKLINE>
    <BLANKLINE>
    This is the number of pairs of edges $((a,b), (c,d))$ such that $a\le c\le d\le b$. i.e., the edge $(c,d)$ is nested inside $(a,b)$...

We can check the references::

    sage: r[1].statistic().references()                                         # optional -- internet
    0: [1]  de Médicis, A., Viennot, X. G., Moments des $q$-polynômes de Laguerre et la bijection de Foata-Zeilberger [[MathSciNet:1288802]]
    1: [2]  Simion, R., Stanton, D., Octabasic Laguerre polynomials and permutation statistics [[MathSciNet:1418763]]...

If you prefer, you can look at this information also in your browser::

    sage: r[1].statistic().browse()                                             # optional -- webbrowser

Another interesting possibility is to look for equidistributed
statistics.  Instead of submitting a list of ``(object, value)``
pairs, we pass a list of pairs ``(objects, values)``::

    sage: r = findstat([(PM(2*n), [m.number_of_nestings() for m in PM(2*n)]) for n in range(5)], depth=0); r # optional -- internet
    0: St000041 (quality [99, 100])
    1: St000042 (quality [99, 100])
    ...

This results tells us that the database contains another entry that
is equidistributed with the number of nestings on perfect matchings
of size at most `10`, namely the number of crossings.  Note that
there is a limit on the number of elements FindStat accepts for a
query, which is currently `1000`.  Queries with more than `1000`
elements are truncated.

Let us now look at a slightly more complicated example, where the
submitted statistic is the composition of a sequence of combinatorial
maps and a statistic known to FindStat.  We use the occasion to
advertise yet another way to pass values to FindStat::

    sage: r = findstat(Permutations, lambda pi: pi.saliances()[0], depth=2)     # optional -- internet

Let us pick one particular result::

    sage: s = next(s for s in r if s.statistic().id() == 51); s                 # optional -- internet
    St000051oMp00061oMp00069 (quality [87, 86])

To obtain the value of the statistic sent to FindStat on a given
object, apply the maps in the list in the given order to this object,
and evaluate the statistic on the result.  For example, let us check
that the result given by FindStat agrees with our statistic on the
following permutation::

    sage: pi = Permutation([3,1,4,5,2]); pi.saliances()[0]
    3

We first have to find out, what the maps and the statistic actually do::

    sage: print(s.statistic().description())                                    # optional -- internet
    The size of the left subtree of a binary tree.

    sage: print(s.statistic().sage_code())                                      # optional -- internet
    def statistic(T):
        return T[0].node_number()

    sage: print("\n\n".join(m.sage_code() for m in s.maps()))                   # optional -- internet, random
    def complement(sigma):
        return sigma.complement()
    <BLANKLINE>
    def increasing_tree_shape(sigma):
        return sigma.increasing_tree_shape()

So, the following should coincide with what we sent FindStat::

    sage: pi.complement().increasing_tree_shape()[0].node_number()
    3

Editing and submitting statistics
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Of course, often a statistic will not be in the database::

    sage: s = findstat([(d, randint(1,1000)) for d in DyckWords(4)]); s         # optional -- internet
    0: St000000: a new statistic on Cc0005: Dyck paths

In this case, and if the statistic might be "interesting", please
consider submitting it to the database using
:meth:`FindStatStatistic.submit`.

Also, you may notice omissions, typos or even mistakes in the
description, the code and the references.  In this case, simply
replace the value by using :meth:`FindStatStatistic.set_description`,
:meth:`FindStatStatistic.set_code` or
:meth:`FindStatStatistic.set_references`, and then
:meth:`FindStatStatistic.submit` your changes for review by the
FindStat team.

Classes and methods
-------------------

"""
#*****************************************************************************
#       Copyright (C) 2015 Martin Rubey <martin.rubey@tuwien.ac.at>,
#
#  Distributed under the terms of the GNU General Public License (GPL)
#  as published by the Free Software Foundation; either version 2 of
#  the License, or (at your option) any later version.
#                  http://www.gnu.org/licenses/
#*****************************************************************************
from __future__ import print_function
from six.moves import range
from six import iteritems, add_metaclass, string_types

from sage.misc.lazy_list import lazy_list
from sage.misc.inherit_comparison import InheritComparisonClasscallMetaclass
from sage.structure.element import Element
from sage.structure.parent import Parent
from sage.structure.unique_representation import UniqueRepresentation

from sage.categories.sets_cat import Sets
from sage.structure.sage_object import SageObject
from sage.structure.richcmp import richcmp

from sage.misc.misc import verbose
from sage.rings.integer import Integer
from sage.databases.oeis import FancyTuple

from ast import literal_eval
from collections import OrderedDict
from copy import deepcopy
import re
import webbrowser
import tempfile
import inspect
import json
import cgi
import requests

# import compatible with py2 and py3
from six.moves.urllib.parse import urlencode
from six.moves.urllib.request import Request, urlopen
from six.moves.urllib.error import HTTPError

# Combinatorial collections
from sage.combinat.alternating_sign_matrix import AlternatingSignMatrix, AlternatingSignMatrices
from sage.combinat.binary_tree import BinaryTree, BinaryTrees
from sage.combinat.core import Core, Cores
from sage.combinat.dyck_word import DyckWord, DyckWords
from sage.combinat.root_system.cartan_type import CartanType_abstract, CartanType
from sage.combinat.gelfand_tsetlin_patterns import GelfandTsetlinPattern, GelfandTsetlinPatterns
from sage.graphs.graph import Graph
from sage.combinat.composition import Composition, Compositions
from sage.combinat.partition import Partition, Partitions
from sage.combinat.ordered_tree import OrderedTree, OrderedTrees
from sage.combinat.parking_functions import ParkingFunction, ParkingFunction_class, ParkingFunctions
from sage.combinat.perfect_matching import PerfectMatching, PerfectMatchings
from sage.combinat.permutation import Permutation, Permutations
from sage.combinat.posets.posets import Poset, FinitePoset
from sage.combinat.posets.poset_examples import Posets
from sage.combinat.tableau import SemistandardTableau, SemistandardTableaux, StandardTableau, StandardTableaux
from sage.combinat.set_partition import SetPartition, SetPartitions
from sage.combinat.skew_partition import SkewPartition, SkewPartitions
from sage.graphs.graph_generators import graphs
from sage.combinat.words.word import Word
from sage.combinat.words.words import Words
from sage.combinat.words.abstract_word import Word_class

######################################################################
# the FindStat URLs
FINDSTAT_URL                = 'https://www.findstat.org/'
FINDSTAT_API                = FINDSTAT_URL + "api/"
FINDSTAT_API_COLLECTIONS    = FINDSTAT_API + 'CollectionsDatabase/'
FINDSTAT_API_STATISTICS     = FINDSTAT_API + 'StatisticsDatabase/'
FINDSTAT_API_MAPS           = FINDSTAT_API + 'MapsDatabase/'

FINDSTAT_URL_LOGIN          = FINDSTAT_URL + "?action=login"
FINDSTAT_URL_COLLECTIONS    = FINDSTAT_URL + 'CollectionsDatabase/'
FINDSTAT_URL_STATISTICS     = FINDSTAT_URL + 'StatisticsDatabase/'
FINDSTAT_URL_EDIT_STATISTIC = FINDSTAT_URL + 'EditStatistic/'
FINDSTAT_URL_NEW_STATISTIC  = FINDSTAT_URL + 'NewStatistic/'
FINDSTAT_URL_MAPS           = FINDSTAT_URL + 'MapsDatabase/'
FINDSTAT_URL_EDIT_MAP       = FINDSTAT_URL + 'EditMap/'
FINDSTAT_URL_NEW_MAP        = FINDSTAT_URL + 'NewMap/'


######################################################################
# the number of values FindStat allows to search for at most
FINDSTAT_MAX_VALUES = 1000
# the number of values FindStat needs at least to search for
FINDSTAT_MIN_VALUES = 3
# the number of maps that FindStat should compose at most to find a match
FINDSTAT_MAX_DEPTH = 4
# the maximal number of maps that FindStat composes by default to find a match
FINDSTAT_DEFAULT_DEPTH = 2
# the number of values FindStat allows to submit at most
FINDSTAT_MAX_SUBMISSION_VALUES = 1200

# separates name from description
FINDSTAT_SEPARATOR_NAME = "\n"
# separates references
FINDSTAT_SEPARATOR_REFERENCES = "\n"
# separates values
FINDSTAT_VALUE_SEPARATOR = ";"
FINDSTAT_MAP_SEPARATOR = "o"
# regexp to recognize a statistic identifier
FINDSTAT_STATISTIC_REGEXP = '^St[0-9]{6}$'
FINDSTAT_MAP_REGEXP = '^Mp[0-9]{5}$'
FINDSTAT_COLLECTION_REGEXP = '^Cc[0-9]{4}$'
FINDSTAT_STATISTIC_PADDED_IDENTIFIER  = "St%06d"
FINDSTAT_MAP_PADDED_IDENTIFIER  = "Mp%05d"
FINDSTAT_COLLECTION_PADDED_IDENTIFIER  = "Cc%04d"

######################################################################

# the format string for using POST
# WARNING: we use cgi.escape to avoid injection problems, thus we expect double quotes as field delimiters.
FINDSTAT_POST_HEADER = """
<script src="https://www.google.com/jsapi"></script>
<script>
    google.load("jquery", "1.3.2");
</script>

<script>
    $(document).ready(function() {$("#form").submit(); });
</script>
"""
FINDSTAT_NEWSTATISTIC_FORM_HEADER = '<form id="form" name="NewStatistic" action="%s" enctype="multipart/form-data" method="post" />'
FINDSTAT_NEWMAP_FORM_HEADER = '<form id="form" name="NewMap" action="%s" enctype="multipart/form-data" method="post" />'
FINDSTAT_FORM_FORMAT = '<input type="hidden" name="%s" value="%s" />'
FINDSTAT_FORM_FOOTER = '</form>'

######################################################################
class FindStat(UniqueRepresentation, SageObject):
    r"""
    The Combinatorial Statistic Finder.

    :class:`FindStat` is a class representing results of queries to
    the FindStat database.  This class is also the entry point to
    edit statistics and new submissions.  Use the shorthand
    :class:`findstat<FindStat>` to call it.

    """
    def __init__(self):
        r"""
        Initialize the database.

        TESTS::

            sage: findstat()
            The Combinatorial Statistic Finder (https://www.findstat.org/)
        """
        # user credentials if provided
        self._user_name  = ""
        self._user_email = ""

    def _repr_(self):
        r"""
        Return the representation of ``self``.

        EXAMPLES::

            sage: findstat()
            The Combinatorial Statistic Finder (https://www.findstat.org/)
        """
        return "The Combinatorial Statistic Finder (%s)" % FINDSTAT_URL

    def browse(self):
        r"""
        Open the FindStat web page in a browser.

        EXAMPLES::

            sage: findstat().browse()                                           # optional -- webbrowser
        """
        webbrowser.open(FINDSTAT_URL)

    def set_user(self, name=None, email=None):
        r"""
        Set the user for this session.

        INPUT:

        - ``name`` -- the name of the user.

        - ``email`` -- an email address of the user.

        This information is used when submitting a statistic with
        :meth:`FindStatStatistic.submit`.

        EXAMPLES::

            sage: findstat().set_user(name="Anonymous", email="invalid@org")

        .. NOTE::

            It is usually more convenient to login into the FindStat
            web page using the :meth:`login` method.
        """
        if not isinstance(name, string_types):
            raise ValueError("The given name is not a string.")
        if not isinstance(email, string_types):
            raise ValueError("The given email address is not a string.")
        self._user_name  = name
        self._user_email = email

    def login(self):
        r"""
        Open the FindStat login page in a browser.

        EXAMPLES::

            sage: findstat().login()                                            # optional -- webbrowser
        """
        webbrowser.open(FINDSTAT_URL_LOGIN)

def _data_from_iterable(iterable, mapping=False, domain=None,
                        codomain=None, check=True):
    """
    Return a list of pairs of lists of the same size, domain, and if
    applicable, codomain.

    INPUT:

    - iterable, a pair of lists of the same size, or an iterable of
      pairs, such that every pair consists of two iterables of the
      same size, or a single element and a single value.  Every
      object must be a :class:`SageObject`.

    - mapping -- (default: ``False``), ``False``, if the codomain is
      ``Integer`` and ``True`` if it is a FindStat collection

    - domain -- (optional), the domain, if ``None`` it is guessed
      from the iterable

    - codomain -- (optional), the codomain, if ``None`` it is guessed
      from the iterable

    """
    if isinstance(iterable, dict):
        iterator = iter(iterable.items())
    else:
        iterator = iter(iterable)
    query0 = next(iterator)
    try:
        query1 = next(iterator)
    except StopIteration:
        # query0 == (elts, vals), which we interpret as [(elts, vals)]
        pre_data = [(query0[0], query0[1])]
        # next(iterator) will raise StopIteration
    try:
        query2 = next(iterator)
        pre_data = [query0, query1, query2]
    except StopIteration:
        # (query0, query1) == (elts, vals) or
        # (query0, query1) == [(elts, vals), (elts, vals)]
        # in the former case, elts has at least 3 elements
        # in both cases, query0 and query1 are not objects
        # if query0 is a long list, we have to raise an error anyway
        elts, vals = list(query0), list(query1)
        if len(elts) == 2:
            if len(vals) != 2:
                raise ValueError("cannot interpret the given argument as a FindStat query")
            pre_data = [elts, vals]
        else:
            if len(elts) != len(vals):
                raise ValueError("FindStat expects the same number of objects (got %s) as values (got %s)" % (len(elts), len(vals)))
            pre_data = [(elts, vals)]

    # pre_data is a list of all elements of the iterator accessed so
    # far, for each of its elements and also the remainder ot the
    # iterator, each element is either a pair ``(object, value)`` or
    # a pair ``(objects, values)``
    elts, vals = pre_data[0]
    if domain is None:
        domain = FindStatCollection(elts)
    if mapping and codomain is None:
        codomain = FindStatCollection(vals)

    all_elements = set()
    def sanitize_pair(elts, vals):
        if domain.is_element(elts):
            elts = [elts]
            if mapping:
                vals = [vals]
            else:
                vals = [Integer(vals)]
        else:
            elts = list(elts)
            if check:
                bad = [elt for elt in elts if not domain.is_element(elt)]
                assert not bad, "%s are not elements of %s" % (bad, domain)
            if mapping:
                vals = list(vals)
            else:
                vals = list(map(Integer, vals))
        if len(elts) != len(vals):
            raise ValueError("FindStat expects the same number of objects as values in each pair")
        if check and mapping:
            bad = [elt for elt in vals if not codomain.is_element(elt)]
            assert not bad, "%s are not elements of %s" % (bad, codomain)
        for elt in elts:
            if elt in all_elements:
                raise ValueError("FindStat expects that every object occurs at most once: %s" % elt)
            all_elements.add(elt)

        return elts, vals

    lazy_data = lazy_list((sanitize_pair(elts, vals)
                           for elts, vals in iterator),
                          initial_values=[sanitize_pair(elts, vals)
                                          for elts, vals in pre_data])
    if mapping:
        return lazy_data, domain, codomain
    return lazy_data, domain


def _data_from_function(function, domain):
    return lazy_list(([elt], [value])
                     for elt, value in domain.first_terms(function))


def _data_from_data(data, max_values):
    """Return the first few pairs (of lists of the same size) with a
    total of at most ``max_values`` objects in the range of the
    collection.

    INPUT:

    - ``data``, an iterable over pairs of lists of the same size

    - ``max_values``, the maximal number of objects (and values) to
      return

    We assume that the number of elements in each pair weakly
    increases, to decide when to stop.
    """
    query = []
    total = min(max_values, FINDSTAT_MAX_VALUES)
    iterator = iter(data)
    while total > 0:
        try:
            elts, vals = next(iterator)
        except StopIteration:
            break
        if total >= len(elts):
            query.append((elts, vals))
            total -= len(elts)
        else:
            break # assuming that the next pair is even larger

    return query


def _distribution_from_data(data, domain, max_values):
    """Return the first few pairs (of lists of the same size) with a
    total of at most ``max_values`` objects in the range of the
    collection, combined by level.

    INPUT:

    - ``data``, an iterable over pairs of lists of the same size

    - ``domain``, a :class:`FindStatCollection`

    - ``max_values``, the maximal number of objects (and values) to
      return

    TESTS::

        sage: from sage.databases.findstat import _distribution_from_data, FindStatCollection
        sage: cc = FindStatCollection(1)                                        # optional -- internet
        sage: n = 3; l = lambda i: [pi for pi in Permutations(n) if pi(1) == i] # optional -- internet
        sage: data = [([pi for pi in l(i)], [pi(1) for pi in l(i)]) for i in range(1,n+1)] # optional -- internet
        sage: _distribution_from_data(data, cc, 10)                             # optional -- internet
        [([[1, 2, 3], [1, 3, 2], [2, 1, 3], [2, 3, 1], [3, 1, 2], [3, 2, 1]],
          [1, 1, 2, 2, 3, 3])]
    """
    lvl_dict = {} # lvl: elts, vals
    total = min(max_values, FINDSTAT_MAX_VALUES)
    iterator = iter(data)
    levels_with_sizes = domain.levels_with_sizes()
    while total > 0:
        try:
            elts, vals = next(iterator)
        except StopIteration:
            break
        if total < len(elts):
            break
        else:
            lvl = domain.element_level(elts[0])
            if levels_with_sizes[lvl] > total:
                # we assume that from now on levels become even larger
                break
            if not all(domain.element_level(elt) == lvl for elt in elts[1:]):
                raise ValueError("cannot combine %s into a distribution" % elts)
            lvl_elts, lvl_vals = lvl_dict.get(lvl, [[], []])
            lvl_dict[lvl] = (lvl_elts + elts, lvl_vals + vals)
            if levels_with_sizes[lvl] == len(lvl_dict[lvl][0]):
                total -= levels_with_sizes[lvl]

    return [(elts, vals) for lvl, (elts, vals) in lvl_dict.items()
            if levels_with_sizes[lvl] == len(elts)]

def _generating_functions_from_dict(gfs, style):
    if style == "dictionary":
        return gfs
    elif style == "list":
        return { key : [ gfs[key][deg] if deg in gfs[key] else 0
                         for deg in range(min(gfs[key]),max(gfs[key])+1) ]
                 for key in gfs }
    elif style == "polynomial":
        from sage.rings.polynomial.polynomial_ring_constructor import PolynomialRing
        from sage.rings.integer_ring import ZZ
        P = PolynomialRing(ZZ,"q")
        q = P.gen()
        return { level : sum( coefficient * q**exponent
                              for exponent,coefficient in iteritems(gen_dict) )
                 for level, gen_dict in iteritems(gfs)}
    else:
        raise ValueError("The argument 'style' (='%s') must be 'dictionary', 'polynomial', or 'list'." % style)

def findstat(query=None, values=None, distribution=None, domain=None,
             depth=FINDSTAT_DEFAULT_DEPTH, max_values=FINDSTAT_MAX_VALUES):
    r"""
    Return matching statistics.

    INPUT:

    One of the following:

    - an integer or a string representing a valid FindStat identifier
      (e.g. 45 or 'St000045').  The keyword arguments ``depth`` and
      ``max_values`` are ignored, ``values`` and ``distribution``
      must be ``None``.

    - a list of pairs of the form ``(object, value)``, or a
      dictionary from sage objects to integer values.  The keyword
      arguments ``depth`` and ``max_values`` are passed to the
      finder, ``values`` and ``distribution`` must be ``None``.

    - a list of pairs of the form (list of objects, list of values),
      or a single pair of the form (list of objects, list of values).
      In each pair there should be as many objects as values.  The
      keyword arguments ``depth`` and ``max_values`` are passed to
      the finder.

    - a collection and a list of pairs of the form (string, value),
      or a dictionary from strings to integer values.  The keyword
      arguments ``depth`` and ``max_values`` are passed to the
      finder.  This should only be used if the collection is not yet
      supported.

    - a collection and a callable.  The callable is used to generate
      ``max_values`` ``(object, value)`` pairs.  The number of terms
      generated may also be controlled by passing an iterable
      collection, such as ``Permutations(3)``.  The keyword arguments
      ``depth`` and ``max_values`` are passed to the finder.

    OUTPUT:

    An instance of a :class:`FindStatStatistic`, represented by

    - the FindStat identifier together with its name, or

    - a list of triples, each consisting of

        - the statistic

        - a list of strings naming certain maps

        - a number which says how many of the values submitted agree
          with the values in the database, when applying the maps in
          the given order to the object and then computing the
          statistic on the result.

    EXAMPLES:

    A particular statistic can be retrieved by its St-identifier or
    number::

        sage: findstat('St000041')                                              # optional -- internet
        St000041: The number of nestings of a perfect matching.

        sage: findstat(51)                                                      # optional -- internet
        St000051: The size of the left subtree of a binary tree.

        sage: findstat('St000042oMp00116')                                      # optional -- internet
        St000042oMp00116

    The database can be searched by providing a list of pairs::

        sage: l = [pi for n in range(7) for pi in Permutations(n)]
        sage: q = findstat([(pi, pi.length()) for pi in l], depth=0); q         # optional -- internet, random
        0: St000018 (quality [100, 100])
        1: St000000: a new statistic on Cc0001: Permutations

    or a dictionary::

        sage: p = findstat({pi: pi.length() for pi in l}, depth=0); p           # optional -- internet, random
        0: St000018 (quality [100, 100])
        1: St000000: a new statistic on Cc0001: Permutations

    Note however, that the results of these two queries need not
    compare equal, because we compare queries by the data
    sent, and the ordering of the data might be different.

    Another possibility is to send a collection and a function.  In
    this case, the function is applied to the first few objects of
    the collection::

        sage: findstat("Permutations", lambda pi: pi.length(), depth=0)         # optional -- internet
        0: St000018 (quality [87, 100])
        1: St000000: a new statistic on Cc0001: Permutations

    To search for a distribution, send a list of lists, or a single pair::

        sage: S = PerfectMatchings(10); findstat((S, [pi.number_of_nestings() for pi in S]), depth=0) # optional -- internet
        0: St000042 (quality [100, 100])
        1: St000041 (quality [9, 100])
        2: St000000: a new statistic on a subset of Cc0012: Perfect matchings

    Alternatively, specify the ``distribution`` parameter::

        sage: findstat(12, distribution=lambda pi: pi.number_of_nestings(), depth=0) # optional -- internet
        0: St000041 (quality [100, 100])
        1: St000042 (quality [100, 100])
        2: St000000: a new statistic on Cc0012: Perfect matchings

    Note that there is a limit, ``FINDSTAT_MAX_VALUES``, on the number
    of elements that may be submitted to FindStat, which is currently
    1000.  Therefore, the interface tries to truncate queries
    appropriately, but this may be impossible, especially with
    distribution searches::

        sage: S = Permutations(7); S.cardinality()                              # optional -- internet
        5040
        sage: findstat((S, [1 for a in S]))                                     # optional -- internet
        Traceback (most recent call last):
        ...
        ValueError: E016: You passed too few elements (0 < 3) to FindStat!

    Finally, we can also retrieve all statistics with a given domain::

        sage: findstat("Cc0024")                                                # optional -- internet
        Set of combinatorial statistics with domain Cc0024: Binary words in FindStat

        sage: findstat(domain="Cores")                                          # optional -- internet
        Set of combinatorial statistics with domain Cc0013: Cores in FindStat

    TESTS::

        sage: findstat("Permutations", lambda x: 1, depth="x")
        Traceback (most recent call last):
        ...
        ValueError: the depth of a FindStat query must be a non-negative integer less than or equal to 4

        sage: findstat("Permutations", lambda x: 1, depth=100)
        Traceback (most recent call last):
        ...
        ValueError: the depth of a FindStat query must be a non-negative integer less than or equal to 4

        sage: S = Permutation
        sage: findstat([(S([1,2]), 1), ([S([1,3,2]), S([1,2])], [2,3])])        # optional -- internet
        Traceback (most recent call last):
        ...
        ValueError: FindStat expects that every object occurs at most once: [1, 2]

    Check that values which can be converted to integers are supported::

        sage: findstat([(la, la[0]/1) for la in Partitions(10)], depth=0)       # optional -- internet
        0: St000147 (quality [100, 100])
        1: St000000: a new statistic on Cc0002: Integer partitions

    """
    try:
        depth = int(depth)
        assert 0 <= depth <= FINDSTAT_MAX_DEPTH
    except (ValueError, AssertionError):
        raise ValueError("the depth of a FindStat query must be a non-negative integer less than or equal to %i" % FINDSTAT_MAX_DEPTH)

    try:
        max_values = int(max_values)
        assert 0 <= max_values <= FINDSTAT_MAX_VALUES
    except (ValueError, AssertionError):
        raise ValueError("the maximal number of values for a FindStat query must be a non-negative integer less than or equal to %i" % FINDSTAT_MAX_VALUES)

    check_collection = True
    def get_data(raw, domain=None):
        all_data, domain = _data_from_iterable(raw, domain=domain,
                                               mapping=False,
                                               check=check_collection)
        data = _data_from_data(all_data, max_values)
        return all_data, data, domain, None

    def get_values(raw, domain=None):
        if callable(raw):
            all_data = _data_from_function(raw, domain)
            function = raw
        else:
            all_data, domain = _data_from_iterable(raw, domain=domain,
                                                   mapping=False,
                                                   check=check_collection)
            function = None
        data = _data_from_data(all_data, max_values)
        return all_data, data, domain, function

    def get_distribution(raw, domain=None):
        if callable(raw):
            all_data = _data_from_function(raw, domain)
            function = raw
        else:
            all_data, domain = _data_from_iterable(raw, domain=domain,
                                                   mapping=False,
                                                   check=check_collection)
            function = None
        data = _distribution_from_data(all_data, domain, max_values)
        return all_data, data, domain, function

    ######################################################################
    if query is None and values is None and distribution is None and domain is None:
        return FindStat()

    if values is not None and distribution is not None:
        raise ValueError("not both of `values` and `distribution` may be given for a FindStat query")

    if values is None and distribution is None:
        if query is None:
            return FindStatStatistics(domain=domain)

        if domain is None:
            if isinstance(query, Integer):
                return FindStatStatistic(query)

            if isinstance(query, string_types):
                if re.match(FINDSTAT_COLLECTION_REGEXP, query):
                    return FindStatStatistics(domain=query)

                return FindStatStatistic(query)

        values, query = query, None

    if query is not None:
        if domain is None:
            domain, query = query, None
        else:
            raise ValueError("the given arguments cannot be used for a FindStat search")

    if domain is not None:
        domain = FindStatCollection(domain)

    if values is not None:
        if isinstance(values, (string_types, Integer)):
            if domain is not None:
                raise ValueError("the domain must not be provided if a statistic identifier is given")
            return FindStatStatisticQuery(values_of=values, depth=depth)

        all_data, data, domain, function = get_values(values, domain)
        return FindStatStatisticQuery(data=data, domain=domain, depth=depth,
                                      all_data=all_data, function=function)

    if distribution is not None:
        if isinstance(distribution, (string_types, Integer)):
            if domain is not None:
                raise ValueError("the domain must not be provided if a statistic identifier is given")
            return FindStatStatisticQuery(distribution_of=distribution, depth=depth)

        all_data, data, domain, function = get_distribution(distribution, domain)
        return FindStatStatisticQuery(data=data, domain=domain, depth=depth,
                                      all_data=all_data, function=function)

    raise ValueError("the given arguments cannot be used for a FindStat search")

######################################################################

def findmap(*args, **kwargs):
    r"""
    Return matching maps.

    INPUT:

    One of the following:

    - an integer or a string representing a valid FindStat map identifier
      (e.g. 45 or 'Mp00045').

    - a list of pairs of the form ``(object, value)``, or a
      dictionary from sage objects to integer values.  The keyword
      arguments ``depth`` and ``max_values`` are passed to the
      finder.

    - a list of pairs of the form ``(list of objects, list of
      values)``, or a single pair of the form ``(list of objects,
      list of values)``.  In each pair there should be as many
      objects as values.  The keyword arguments ``depth`` and
      ``max_values`` are passed to the finder.

    - a collection and a list of pairs of the form (string, value),
      or a dictionary from strings to integer values.  The keyword
      arguments ``depth`` and ``max_values`` are passed to the
      finder.  This should only be used if the collection is not yet
      supported.

    - a collection and a callable.  The callable is used to generate
      ``max_values`` ``(object, value)`` pairs.  The number of terms
      generated may also be controlled by passing an iterable
      collection, such as ``Permutations(3)``.  The keyword arguments
      ``depth`` and ``max_values`` are passed to the finder.

    OUTPUT:

    An instance of a :class:`FindStatMap`, :class:`FindStatMapQuery`
    or :class:`FindStatMaps`.

    EXAMPLES:

    A particular map can be retrieved by its Mp-identifier or
    number::

        sage: findmap('Mp00001')                                                # optional -- internet
        Mp00001: to semistandard tableau via monotone triangles

        sage: findmap(1)                                                        # optional -- internet
        Mp00001: to semistandard tableau via monotone triangles

        sage: findmap("Mp00087oMp00058")                                        # optional -- internet
        Mp00087oMp00058

    The database can be searched by providing a list of pairs::

        sage: l = [pi for n in range(7) for pi in Permutations(n)]
        sage: q = findmap([(pi, pi.inverse().complement()) for pi in l], depth=2); q      # optional -- internet, random
        0: Mp00066oMp00064 (quality [100])

    or a dictionary::

        sage: p = findmap({pi: pi.inverse().complement() for pi in l}, depth=2); p        # optional -- internet, random
        0: Mp00066oMp00064 (quality [100])

    Note however, that the results of these two queries need not
    compare equal, because we compare queries by the data
    sent, and the ordering of the data might be different.

    Another possibility is to send a collection and a function.  In
    this case, the function is applied to the first few objects of
    the collection::

        sage: findmap("Permutations", lambda pi: pi.inverse().complement(), depth=2)      # optional -- internet
        0: Mp00066oMp00064 (quality [100])

    In rare cases, it may not be possible to guess the codomain of a
    map, in which case it can be provided as second argument or
    keyword argument::

        sage: findmap("Dyck paths", "Perfect matchings", lambda D: [(a+1, b) for a,b in D.tunnels()]) # optional -- internet
        0: Mp00146 (quality [100])

        sage: findmap("Dyck paths", "Set partitions", lambda D: [(a+1, b) for a,b in D.tunnels()]) # optional -- internet
        0: Mp00092oMp00146 (quality [50])

    Finally, we can also retrieve all maps with a given domain or codomain::

        sage: findmap("Cc0024")                                                 # optional -- internet
        Set of combinatorial maps with domain Cc0024: Binary words used by FindStat

        sage: findmap(codomain="Cores")                                         # optional -- internet
        Set of combinatorial maps with codomain Cc0013: Cores used by FindStat

    """
    if len(args) > 3:
        raise TypeError("findmap takes at most 3 positional arguments (%s given)" % len(args))

    bad_args = set(kwargs).difference(["values", "distribution",
                                       "domain", "codomain",
                                       "depth", "max_values"])
    if bad_args:
        raise TypeError("findmap got unexpected keyword arguments '%s'" % bad_args)

    max_values = kwargs.get("max_values", FINDSTAT_MAX_VALUES)
    depth = kwargs.get("depth", FINDSTAT_DEFAULT_DEPTH)
    values = kwargs.get("values", None)
    distribution = kwargs.get("distribution", None)
    domain = kwargs.get("domain", None)
    codomain = kwargs.get("codomain", None)

    try:
        depth = int(depth)
        assert 0 <= depth <= FINDSTAT_MAX_DEPTH
    except (ValueError, AssertionError):
        raise ValueError("The depth of a FindStat query must be a non-negative integer less than or equal to %i." % FINDSTAT_MAX_DEPTH)

    try:
        max_values = int(max_values)
        assert 0 <= max_values <= FINDSTAT_MAX_VALUES
    except (ValueError, AssertionError):
        raise ValueError("The maximal number of values for a FindStat query must be a non-negative integer less than or equal to %i." % FINDSTAT_MAX_VALUES)

    check_collection = True
    def get_data(raw, domain=None, codomain=None):
        all_data, domain, codomain = _data_from_iterable(raw, domain=domain,
                                                         codomain=codomain,
                                                         mapping=True,
                                                         check=check_collection)
        data = _data_from_data(all_data, max_values)
        return all_data, data, domain, codomain, None

    def get_values(raw, domain=None, codomain=None):
        if callable(raw):
            all_data = _data_from_function(raw, domain)
            if codomain is None:
                codomain = FindStatCollection(all_data[0][1][0])
            function = raw
        else:
            all_data, domain, codomain = _data_from_iterable(raw, domain=domain,
                                                             codomain=codomain,
                                                             mapping=True,
                                                             check=check_collection)
            function = None
        data = _data_from_data(all_data, max_values)
        return all_data, data, domain, codomain, function

    def get_distribution(raw, domain=None, codomain=None):
        if callable(raw):
            all_data = _data_from_function(raw, domain)
            function = raw
        else:
            all_data, domain, codomain = _data_from_iterable(raw, domain=domain,
                                                             codomain=codomain,
                                                             mapping=True,
                                                             check=check_collection)
            function = None
        data = _distribution_from_data(all_data, domain, max_values)
        return all_data, data, domain, codomain, function

    def is_collection(arg):
        try:
            FindStatCollection(arg)
            return True
        except ValueError:
            return False

    def check_domain(arg, domain):
        if domain is not None:
            raise TypeError("the domain was specified twice, as positional argument (%s) and as keyword domain=%s" % (arg, domain))
        return arg

    def check_codomain(arg, codomain):
        if codomain is not None:
            raise TypeError("the codomain was specified twice, as positional argument (%s) and as keyword codomain=%s" % (arg, codomain))
        return arg

    def check_values(arg, values):
        if values is not None:
            raise TypeError("values were specified twice, as positional argument (%s) and as keyword values=%s" % (arg, values))
        return arg

    ######################################################################
    if values is not None and distribution is not None:
        raise ValueError("Not both of `values` and `distribution` may be given for a FindStat query.")

    if len(args) == 1:
        if (values is None and distribution is None
            and domain is None and codomain is None
            and (isinstance(args[0], Integer)
                 or (isinstance(args[0], string_types)
                     and not is_collection(args[0])))):
            return FindStatMap(args[0])

        elif (isinstance(args[0], string_types) and
              is_collection(args[0])):
            domain = check_domain(args[0], domain)

        else:
            values = check_values(args[0], values)

    elif len(args) == 2:
        domain = check_domain(args[0], domain)
        if isinstance(args[1], (Integer, string_types)):
            codomain = check_codomain(args[1], codomain)
        else:
            values = check_values(args[1], values)

    elif len(args) == 3:
        domain = check_domain(args[0], domain)
        codomain = check_codomain(args[1], codomain)
        values = check_values(args[2], values)

    if domain is not None:
        domain = FindStatCollection(domain)
    if codomain is not None:
        codomain = FindStatCollection(codomain)

    if (values is None and distribution is None
        and (domain is not None or codomain is not None)):
        return FindStatMaps(domain=domain, codomain=codomain)

    if values is not None:
        if isinstance(values, (string_types, Integer)):
            if domain is not None or codomain is not None:
                raise ValueError("domain and codomain must not be provided if a map identifier is given")
            return FindStatMapQuery(values_of=values, depth=depth)

        all_data, data, domain, codomain, function = get_values(values, domain, codomain)
        return FindStatMapQuery(data=data, domain=domain, codomain=codomain, depth=depth,
                                all_data=all_data, function=function)

    if distribution is not None:
        if isinstance(distribution, (string_types, Integer)):
            if domain is not None or codomain is not None:
                raise ValueError("domain and codomain must not be provided if a map identifier is given")
            return FindStatMapQuery(distribution_of=values, depth=depth)

        all_data, data, domain, function = get_distribution(distribution, domain)
        return FindStatMapQuery(data=data, domain=domain, codomain=codomain, depth=depth,
                                all_data=all_data, function=function)

    raise ValueError("The given arguments cannot be used for a FindStat search.")

######################################################################

class FindStatFunction(SageObject):
    """
    A class containing the common methods of :class:`FindStatMap` and
    :class:`FindStatStatistic`.
    """
    def __call__(self, elt):
        if self._function is False:
            raise ValueError("execution of code provided by FindStat is not enabled for %s" % self)
        if self._function is True:
            if not self.sage_code():
                raise ValueError("there is no code available for %s" % self)
            from sage.repl.preparse import preparse
            try:
                l = {}
                code = "from sage.all import *\n" + preparse(self.sage_code())
                exec(code, l)
            except SyntaxError:
                raise ValueError("could not execute code for %s" % self)
            if isinstance(self, FindStatStatistic):
                self._function = eval("statistic", l)
            elif isinstance(self, FindStatMap):
                self._function = eval(self.code_name(), l)
            else:
                raise ValueError("cannot execute code for %s" % self)

        return self._function(elt)

    def __repr__(self):
        r"""
        Return the representation of the FindStat statistic or map.

        OUTPUT:

        A string, the identifier and the name of the statistic.  If
        the statistic was modified (see :meth:`modified`) this is
        also indicated.

        EXAMPLES::

            sage: findstat(1)                                                   # optional -- internet
            St000001: The number of reduced words for a permutation.

            sage: findstat(914)                                                 # optional -- internet
            St000914: The sum of the values of the Möbius function of a poset.
        """
        if self._modified:
            return "%s(modified): %s" % (self.id_str(), self.name())
        else:
            return "%s: %s" % (self.id_str(), self.name())

    def is_new(self):
        """
        Return ``True`` if and only if the statistic or map is new.
        """
        return self.id() == 0

    def _data(self):
        # initializes self._modified_data on first call
        if not self.is_new():
            self._fetch_data()
        # some of the data are lists, so we need to deepcopy
        return deepcopy(self._modified_data)

    def id(self):
        r"""
        Return the FindStat identifier of the statistic or map.

        OUTPUT:

        The FindStat identifier of the statistic or map, as an integer.

        EXAMPLES::

            sage: findstat(1).id()                                              # optional -- internet
            1
        """
        return int(self._id[2:])

    def id_str(self):
        r"""
        Return the FindStat identifier of the statistic or map.

        OUTPUT:

        The FindStat identifier of the statistic or map, as a string.

        EXAMPLES::

            sage: findstat(1).id_str()                                          # optional -- internet
            'St000001'
        """
        return self._id

    def description(self):
        r"""
        Return the description of the statistic or map.

        OUTPUT:

        A string.  For statistics, the first line is used as name.

        EXAMPLES::

            sage: print(findstat(1).description())                              # optional -- internet
            The number of reduced words for a permutation.
            <BLANKLINE>
            This is...
        """
        return self._data()["Description"]

    def set_description(self, value):
        r"""
        Set the description of the statistic or map.

        INPUT:

        - a string -- for statistics, this is the name of the
          statistic followed by its description on a separate line.

        This information is used when submitting the statistic or map with
        :meth:`submit`.

        EXAMPLES::

            sage: q = findstat([(d, randint(1,1000)) for d in DyckWords(4)]); q[-1]  # optional -- internet
            St000000: a new statistic on Cc0005: Dyck paths
            sage: q[-1].set_description("Random values on Dyck paths.\nNot for submission.")   # optional -- internet
            sage: q[-1]                                                              # optional -- internet
            St000000(modified): a new statistic on Cc0005: Dyck paths
            sage: print(q[-1].description())                                         # optional -- internet
            Random values on Dyck paths.
            Not for submission.

        """
        if value != self.description():
            self._modified = True
            self._modified_data["Description"] = value

    def name(self):
        r"""
        Return the name of the statistic or map.

        OUTPUT:

        A string.  For statistics, this is just the first line of the
        description.

        EXAMPLES::

            sage: findstat(1).name()                                            # optional -- internet
            'The number of reduced words for a permutation.'

        """
        import sys
        if sys.version_info[0] < 3:
            return self._data()["Name"].encode("utf-8")
        return self._data()["Name"]

    def references(self):
        r"""
        Return the references associated with the statistic or map.

        OUTPUT:

        An instance of :class:`sage.databases.oeis.FancyTuple`, each
        item corresponds to a reference.

        EXAMPLES::

            sage: findstat(1).references()                                      # optional -- internet
            0: [1]  Edelman, P., Greene, C., Balanced tableaux [[MathSciNet:0871081]]
            1: [2]  Number of simple allowable sequences on 1..n containing the permutation 12...n. [[OEIS:A005118]]
            2: [3]  Total number of reduced decompositions for all permutations in S_n. [[OEIS:A246865]]
        """
        result = []
        refs = self.references_raw()
        if refs:
            refs = refs.split(FINDSTAT_SEPARATOR_REFERENCES)
        else:
            return FancyTuple([])
        bibs = self._data()["Bibliography"]
        for ref in refs:
            parts = ref.partition("[[")
            parts = parts[:-1] + parts[2].partition("]]")
            comment = parts[0]
            link = parts[2]
            if link == "":
                result.append(ref)
            else:
                try:
                    bibitem = bibs[link]
                except KeyError:
                    # this means that the link is unhandled
                    result.append(ref)
                else:
                    author_title = ", ".join(e for e in [bibitem["Author"], bibitem["Title"]]
                                             if e)
                    result.append(comment + author_title + " " + "".join(parts[1:]))

        import sys
        if sys.version_info[0] < 3:
            return FancyTuple([ref.encode("utf-8") for ref in result])

        return FancyTuple([ref for ref in result])

    def references_raw(self):
        r"""
        Return the unrendered references associated with the statistic or map.
        """
        return self._data()["References"]

    def set_references_raw(self, value):
        r"""
        Set the references associated with the statistic or map.

        INPUT:

        - a string -- the individual references should be separated
          by FINDSTAT_SEPARATOR_REFERENCES, which is "\\r\\n".

        OUTPUT:

        This information is used when submitting the statistic with
        :meth:`submit`.

        EXAMPLES::

            sage: q = findstat([(d, randint(1, 1000)) for d in DyckWords(4)]); q[-1]      # optional -- internet
            St000000: a new statistic on Cc0005: Dyck paths
            sage: q[-1].set_references_raw("[1] The wonders of random Dyck paths, Anonymous Coward, [[arXiv:1102.4226]].\r\n[2] [[oeis:A000001]]")  # optional -- internet
            sage: q[-1].references()                                                      # optional -- internet
            0: [1] The wonders of random Dyck paths, Anonymous Coward, [[arXiv:1102.4226]].
            1: [2] [[oeis:A000001]]

        """
        if value != self.references_raw():
            self._modified = True
            self._modified_data["References"] = value

    def sage_code(self):
        r"""
        Return the sage code associated with the statistic or map.

        OUTPUT:

        An empty string or a string of the form::

            def statistic(x):
                ...

        or::

            def mapping(x):
                ...

        EXAMPLES::

            sage: print(findstat(1).code())                                     # optional -- internet
            def statistic(x):
                return sum(1 for _ in x.reduced_words_iterator())

        """
        return self._data()["SageCode"]

    def set_sage_code(self, value):
        r"""
        Set the code associated with the statistic or map.

        INPUT:

        - a string -- contributors are encouraged to submit sage code
          in the form::

            def statistic(x):
                ...

        or::

            def mapping(x):
                ...

        This information is used when submitting the statistic with
        :meth:`submit`.

        EXAMPLES::

            sage: q = findstat([(d, randint(1,1000)) for d in DyckWords(4)]); q[-1]       # optional -- internet
            St000000: a new statistic on Cc0005: Dyck paths
            sage: q[-1].set_sage_code("def statistic(x):\r\n    return randint(1,1000)")  # optional -- internet
            sage: print(q[-1].sage_code())                                                # optional -- internet
            def statistic(x):
                return randint(1,1000)
        """
        if value != self.sage_code():
            self._modified = True
            self._modified_data["SageCode"] = value

@add_metaclass(InheritComparisonClasscallMetaclass)
class FindStatStatistic(Element, FindStatFunction):
    r"""
    A FindStat statistic.

    :class:`FindStatStatistic` is a class representing a
    combinatorial statistic available in the FindStat database.

    This class provides methods to inspect and update various
    properties of these statistics.

    INPUT:

    - a string or an integer representing its FindStat id.

    EXAMPLES::

        sage: from sage.databases.findstat import FindStatStatistic
        sage: FindStatStatistic(1)                                              # optional -- internet
        St000001: The number of reduced words for a permutation.

    .. SEEALSO::

        :class:`FindStatStatistics`

    """
    @staticmethod
    def __classcall_private__(cls, entry):
        """
        Retrieve a statistic from the database.

        TESTS::

            sage: from sage.databases.findstat import FindStatStatistic
            sage: FindStatStatistic("abcdefgh")                                 # optional -- internet
            Traceback (most recent call last):
            ...
            ValueError: the value 'abcdefgh' is not a valid FindStat statistic identifier
        """
        return FindStatStatistics()(entry)

    def __init__(self, parent, id, new_data=None):
        self._function = False # determines whether FindStat code may be executed
        self._modified_data = None
        self._modified_first_terms = None
        self._modified_first_terms_raw = None
        self._modified = False # set in every method modifying the data
        self._id = id
        if self.is_new():
            self._initialize_data(new_data)
            self._initialize_first_terms(new_data)
            self._function = new_data.function()
        Element.__init__(self, parent)

    def __reduce__(self):
        """Return a function and its arguments needed to create this
        statistic.

        TESTS::

            sage: from sage.databases.findstat import FindStatStatistic
            sage: c = FindStatStatistic(62)                                     # optional -- internet
            sage: loads(dumps(c)) == c                                          # optional -- internet
            True

        """
        return FindStatStatistic, (self.id(),)

    def _fetch_data(self):
        if self._modified_data is not None:
            return
        fields = "Bibliography,Code,Description,Domain,Name,References,SageCode"
        fields_Bibliography = "Author,Title"
        url = (FINDSTAT_API_STATISTICS + self._id
               + "?fields=" + fields
               + "&fields[Bibliography]=" + fields_Bibliography)
        included = json.load(urlopen(url))["included"]
        # slightly simplify the representation
        data = {key: val for key, val in included["Statistics"][self._id].items()}
        # we replace the list of identifiers in Bibliography with the dictionary
        data["Bibliography"] = included["References"]
        self._modified_data = data

    def _initialize_data(self, data):
        self._modified_data = {"Bibliography": {},
                               "Code": data.code(),
                               "Description" : "",
                               "Domain": data.domain(),
                               "Name": "a new statistic on %s" % data.domain(),
                               "References": "",
                               "SageCode": ""}

    # WARNING: before modifying self._modified_first_terms_raw this
    # method or self._fetch_first_terms, or initialize_first_terms
    # has to be called; it is cached because it needs to be executed
    # only once
    def _fetch_first_terms_raw(self):
        r"""
        Initialize the first terms of the statistic, as (string, value)
        pairs.

        This fetches the data from FindStat, and sets
        ``self._modified_first_terms_raw``.
        """
        if self._modified_first_terms_raw is not None:
            return
        fields = "Values"
        url = FINDSTAT_API_STATISTICS + self._id + "?fields=" + fields
        values = json.load(urlopen(url))["included"]["Statistics"][self._id]["Values"]
        values = [tuple(pair) for pair in values]
        self._modified_first_terms_raw = values

    # WARNING: before modifying self._modified_first_terms this
    # method or initialize_first_terms has to be called; it is cached because it needs to be
    # executed only once
    def _fetch_first_terms(self):
        r"""
        Initialize the first terms of the statistic, as ``(object,
        value)`` pairs.

        This sets ``self._modified_first_terms``.

        """
        if self._modified_first_terms is not None:
            return
        from_str = self.collection().from_string()
        self._fetch_first_terms_raw()
        values = [(from_str(obj), Integer(val))
                  for obj, val in self._modified_first_terms_raw]
        self._modified_first_terms = values

    def _initialize_first_terms(self, data):
        r"""
        Initialize the first terms of the statistic, as ``(object,
        value)`` pairs.

        This sets ``self._modified_first_terms``.

        """
        to_str = self.collection().to_string()
        self._modified_first_terms = [(objs[0], vals[0])
                                      for objs, vals in data.data() if len(vals) == 1]
        self._modified_first_terms_raw = [(to_str(obj), Integer(val))
                                          for obj, val in self._modified_first_terms]

    def first_terms(self):
        r"""
        Return the first terms of the statistic.

        OUTPUT:

        A dictionary from sage objects representing an element of the
        appropriate collection to integers.  If the statistic is in
        the FindStat database, the dictionary contains exactly the
        pairs in the database.

        EXAMPLES::

            sage: findstat(1).first_terms()[Permutation([1,4,3,2])]             # optional -- internet
            2
        """
        # initialize self._modified_first_terms and
        # self._modified_first_terms_raw on first call
        if not self.is_new():
            self._fetch_first_terms()
        # a shallow copy suffices - tuples are immutable
        return dict(self._modified_first_terms)

    def _first_terms_raw(self):
        # initialize self._modified_first_terms_raw on first call
        if not self.is_new():
            self._fetch_first_terms_raw()
        # a shallow copy suffices - tuples are immutable
        return self._modified_first_terms_raw[:]

    def first_terms_str(self):
        r"""
        Return the first terms of the statistic in the format needed
        for a FindStat query.

        OUTPUT:

        A string, where each line is of the form ``object => value``,
        where ``object`` is the string representation of an element
        of the appropriate collection as used by FindStat and value
        is an integer.

        EXAMPLES::

            sage: findstat(1).first_terms_str()[:9]                             # optional -- internet
            u'[] => 1\r\n'

        """
        return "\r\n".join(key + " => " + str(val)
                           for (key, val) in self._first_terms_raw())

    def set_first_terms(self, values):
        r"""
        Update the first terms of the statistic.

        INPUT:

        - a list of pairs of the form ``(object, value)`` where
          ``object`` is a sage object representing an element of the
          appropriate collection and ``value`` is an integer.

        This information is used when submitting the statistic with
        :meth:`submit`.

        """
        to_str = self.collection().to_string()
        new = [tuple(to_str(obj), value) for obj, value in values]
        if sorted(new) != sorted(self.first_terms_str()):
            # initialize self._modified_first_term, because
            # self.first_terms_str() initializes only
            # self._modified_first_term_raw
            if not self.is_new():
                self._fetch_first_terms()
            self._modified = True
            self._modified_first_terms_raw = new
            self._modified_first_terms = values

    def collection(self):
        r"""
        Return the FindStat collection of the statistic.

        OUTPUT:

        The FindStat collection of the statistic as an instance of
        :class:`FindStatCollection`.

        EXAMPLES::

            sage: findstat(1).collection()                                      # optional -- internet
            Cc0001: Permutations
        """
        return FindStatCollection(self._data()["Domain"])

    domain = collection

    def code(self):
        r"""
        Return the code associated with the statistic or map.

        OUTPUT:

        A string.  Contributors are encouraged to submit sage code in the form::

            def statistic(x):
                ...

        but the string may also contain code for other computer
        algebra systems.

        EXAMPLES::

            sage: print(findstat(1).code())                                     # optional -- internet
            def statistic(x):
                return sum(1 for _ in x.reduced_words_iterator())

            sage: print(findstat(118).code())                                   # optional -- internet, random
            (* in Mathematica *)
            tree = {{{{}, {}}, {{}, {}}}, {{{}, {}}, {{}, {}}}};
            Count[tree, {{___}, {{___}, {{___}, {___}}}}, {0, Infinity}]
        """
        return self._data()["Code"]

    def set_code(self, value):
        r"""
        Set the code associated with the statistic.

        INPUT:

        - a string -- contributors are encouraged to submit sage code
          in the form::

            def statistic(x):
                ...

        This information is used when submitting the statistic with
        :meth:`submit`.

        EXAMPLES::

            sage: q = findstat([(d, randint(1,1000)) for d in DyckWords(4)]); q[-1] # optional -- internet
            St000000: a new statistic on Cc0005: Dyck paths
            sage: q[-1].set_code("def statistic(x):\r\n    return randint(1,1000)") # optional -- internet
            sage: print(q[-1].code())                                               # optional -- internet
            def statistic(x):
                return randint(1,1000)
        """
        if value != self.code():
            self._modified = True
            self._modified_data["Code"] = value

    def _generating_functions_dict(self):
        r""" Return the generating functions of ``self`` in a dictionary,
        computed from ``self.first_terms``.

        TESTS:

            sage: q = findstat((DyckWords(4), range(14)))                       # optional -- internet, indirect doctest
            sage: q.generating_functions()                                      # optional -- internet, indirect doctest
            {4: q^13 + q^12 + q^11 + q^10 + q^9 + q^8 + q^7 + q^6 + q^5 + q^4 + q^3 + q^2 + q + 1}

            sage: C = AlternatingSignMatrices(4)
            sage: q = findstat((C, range(C.cardinality())))                     # optional -- internet, indirect doctest
            sage: q.generating_functions()                                      # optional -- internet, indirect doctest
            {4: q^41 + q^40 + q^39 + q^38 + q^37 + q^36 + q^35 + q^34 + q^33 + q^32 + q^31 + q^30 + q^29 + q^28 + q^27 + q^26 + q^25 + q^24 + q^23 + q^22 + q^21 + q^20 + q^19 + q^18 + q^17 + q^16 + q^15 + q^14 + q^13 + q^12 + q^11 + q^10 + q^9 + q^8 + q^7 + q^6 + q^5 + q^4 + q^3 + q^2 + q + 1}
        """
        gfs = {}
        lvls = {}
        domain = self.collection()
        levels_with_sizes = domain.levels_with_sizes()
        for elt, val in self.first_terms().items():
            lvl = domain.element_level(elt)
            if lvl not in levels_with_sizes:
                continue
            if lvl not in gfs:
                gfs[lvl] = {}
            gfs[lvl][val] = gfs[lvl].get(val, 0) + 1
            lvls[lvl] = lvls.get(lvl, 0) + 1

        for lvl, size in lvls.items():
            if size < levels_with_sizes[lvl]:
                del gfs[lvl]
        return gfs

    def generating_functions(self, style="polynomial"):
        r"""
        Return the generating functions of ``self`` in a dictionary.

        The keys of this dictionary are the levels for which the
        generating function of ``self`` can be computed from the data
        of this statistic, and each value represents a generating
        function for one level, as a polynomial, as a dictionary, or as
        a list of coefficients.

        INPUT:

        - a string -- (default:"polynomial") can be
          "polynomial", "dictionary", or "list".

        OUTPUT:

        - if ``style`` is ``"polynomial"``, the generating function is
          returned as a polynomial.

        - if ``style`` is ``"dictionary"``, the generating function is
          returned as a dictionary representing the monomials of the
          generating function.

        - if ``style`` is ``"list"``, the generating function is
          returned as a list of coefficients of the generating function.

        EXAMPLES::

            sage: st = findstat(18)                                             # optional -- internet

            sage: st.generating_functions()                                     # optional -- internet
            {1: 1,
             2: q + 1,
             3: q^3 + 2*q^2 + 2*q + 1,
             4: q^6 + 3*q^5 + 5*q^4 + 6*q^3 + 5*q^2 + 3*q + 1,
             5: q^10 + 4*q^9 + 9*q^8 + 15*q^7 + 20*q^6 + 22*q^5 + 20*q^4 + 15*q^3 + 9*q^2 + 4*q + 1,
             6: q^15 + 5*q^14 + 14*q^13 + 29*q^12 + 49*q^11 + 71*q^10 + 90*q^9 + 101*q^8 + 101*q^7 + 90*q^6 + 71*q^5 + 49*q^4 + 29*q^3 + 14*q^2 + 5*q + 1}

            sage: st.generating_functions(style="dictionary")                   # optional -- internet
            {1: {0: 1},
             2: {0: 1, 1: 1},
             3: {0: 1, 1: 2, 2: 2, 3: 1},
             4: {0: 1, 1: 3, 2: 5, 3: 6, 4: 5, 5: 3, 6: 1},
             5: {0: 1, 1: 4, 2: 9, 3: 15, 4: 20, 5: 22, 6: 20, 7: 15, 8: 9, 9: 4, 10: 1},
             6: {0: 1,
              1: 5,
              2: 14,
              3: 29,
              4: 49,
              5: 71,
              6: 90,
              7: 101,
              8: 101,
              9: 90,
              10: 71,
              11: 49,
              12: 29,
              13: 14,
              14: 5,
              15: 1}}

            sage: st.generating_functions(style="list")                         # optional -- internet
            {1: [1],
             2: [1, 1],
             3: [1, 2, 2, 1],
             4: [1, 3, 5, 6, 5, 3, 1],
             5: [1, 4, 9, 15, 20, 22, 20, 15, 9, 4, 1],
             6: [1, 5, 14, 29, 49, 71, 90, 101, 101, 90, 71, 49, 29, 14, 5, 1]}
        """
        return _generating_functions_from_dict(self._generating_functions_dict(),
                                               style=style)

    def oeis_search(self, search_size=32, verbose=True):
        r"""
        Search the OEIS for the generating function of the statistic.

        INPUT:

        - ``search_size`` (default:32) the number of integers in the
          sequence. If too big, the OEIS result is corrupted.

        - ``verbose`` (default:True) if true, some information about
          the search are printed.

        OUTPUT:

        - a tuple of OEIS sequences, see
          :meth:`sage.databases.oeis.OEIS.find_by_description` for more
          information.

        EXAMPLES::

            sage: st = findstat(18)                                             # optional -- internet

            sage: st.oeis_search()                                              # optional -- internet
            Searching the OEIS for "1  1,1  1,2,2,1  1,3,5,6,5,3,1  1,4,9,15,20,22,20,15,9,4,1  ..."
            0: A008302: Triangle of Mahonian numbers T(n,k)...

            sage: st.oeis_search(search_size=13)                                # optional -- internet
            Searching the OEIS for "1  1,1  1,2,2,1  ..."
            0: A008302: Triangle of Mahonian numbers T(n,k)...
        """
        from sage.databases.oeis import oeis

        gen_funcs = self.generating_functions(style="list")

        OEIS_string = ""
        keys = sorted(gen_funcs.keys())
        counter = 0
        for key in keys:
            gen_func = gen_funcs[key]
            while gen_func[0] == 0:
                gen_func.pop(0)
            # we strip the result according to the search size. -- stumpc5, 2015-09-27
            gen_func = gen_func[:search_size]
            counter += len(gen_func)
            if search_size > 0:
                search_size -= len(gen_func)
            OEIS_func_string     = ",".join( str(coefficient) for coefficient in gen_func )
            OEIS_string         += OEIS_func_string + "  "
        OEIS_string = OEIS_string.strip()
        if counter >= 4:
            if verbose:
                print('Searching the OEIS for "%s"' % OEIS_string)
            return oeis( OEIS_string )
        else:
            if verbose:
                print("Too little information to search the OEIS for this statistic (only %s values given)." % counter)
            return

    ######################################################################
    # browse current statistic
    ######################################################################

    def browse(self):
        r"""
        Open the FindStat web page of the statistic in a browser.

        EXAMPLES::

            sage: findstat(41).browse()                                         # optional -- webbrowser
        """
        if self.is_new():
            self.submit()
        else:
            webbrowser.open(FINDSTAT_URL_STATISTICS + self.id_str())

    def submit(self):
        r"""
        Open the FindStat web page for editing the statistic in a browser.
        """
        args = dict()
        args["OriginalStatistic"] = self.id_str()
        args["Domain"]            = str(self.collection().id_str())
        args["Values"]            = self.first_terms_str()
        args["Description"]       = self.description()
        args["References"]        = self.references_raw()
        args["Code"]              = self.code()
        args["SageCode"]          = self.sage_code()
        args["CurrentAuthor"]     = FindStat()._user_name
        args["CurrentEmail"]      = FindStat()._user_email

        f = tempfile.NamedTemporaryFile(mode='w', suffix='.html', delete=False)
        verbose("Created temporary file %s" % f.name, caller_name='FindStat')
        f.write(FINDSTAT_POST_HEADER)
        if self.id() == 0:
            f.write(FINDSTAT_NEWSTATISTIC_FORM_HEADER % FINDSTAT_URL_NEW_STATISTIC)
        else:
            f.write(FINDSTAT_NEWSTATISTIC_FORM_HEADER %(FINDSTAT_URL_EDIT_STATISTIC+self.id_str()))
        for key, value in iteritems(args):
            if value:
                verbose("writing argument %s" % key, caller_name='FindStat')
                value_encoded = cgi.escape(str(value), quote=True)
                verbose("%s" % value_encoded, caller_name='FindStat')
                f.write((FINDSTAT_FORM_FORMAT %(key, value_encoded)))
            else:
                verbose("skipping argument %s because it is empty" % key, caller_name='FindStat')
        f.write(FINDSTAT_FORM_FOOTER)
        f.close()
        verbose("Opening file with webbrowser", caller_name='FindStat')
        webbrowser.open(f.name)

    # editing and submitting is really the same thing
    edit = submit

_all_statistics = {}

class FindStatStatistics(UniqueRepresentation, Parent):
    r"""
    The class of FindStat statistics.

    The elements of this class are combinatorial statistics currently
    in FindStat.

    EXAMPLES:

    We can print a nice list of the first few statistics currently in
    FindStat, sorted by domain::

        sage: from sage.databases.findstat import FindStatCollections, FindStatStatistic, FindStatStatistics
        sage: for cc, _ in zip(FindStatCollections(), range(3)):                # optional -- internet, random
        ....:     print(cc.name(style="plural"))
        ....:     for st, _ in zip(FindStatStatistics(cc), range(3)):
        ....:         print("    " + st.name())
        Semistandard tableaux
            The cocharge of a semistandard tableau.
            The charge of a semistandard tableau.
            The sum of the entries of a semistandard tableau.
        Gelfand-Tsetlin patterns
            The number of circled entries.
            The number of boxed entries.
            The number of special entries.
        Cores
            The length of a core.
            The size of a core.
            The number of strong covers of a core.

    """
    def __init__(self, domain=None):
        """
        TESTS::

            sage: from sage.databases.findstat import FindStatStatistics
            sage: M = FindStatStatistics()                                      # optional -- internet
            sage: TestSuite(M).run()                                            # optional -- internet
        """
        if domain is None:
            self._domain = None
        else:
            self._domain = FindStatCollection(domain)
        self._identifiers = None
        Parent.__init__(self, category=Sets())

    def _element_constructor_(self, id):
        """Initialize a FindStat statistic.

        INPUT:

        - ``id`` -- a string containing the FindStat identifier of
          the statistic, or the corresponding integer

        EXAMPLES::

            sage: from sage.databases.findstat import FindStatStatistic
            sage: FindStatStatistic(1)                                          # optional -- internet
            St000001: The number of reduced words for a permutation.

        """
        if isinstance(id, self.Element):
            return id
        if isinstance(id, FindStatStatisticQuery):
            if id.submittable():
                return self.element_class(self, FINDSTAT_STATISTIC_PADDED_IDENTIFIER % 0, id)
            raise ValueError("the given query cannot be submitted")
        if isinstance(id, (int, Integer)):
            id = FINDSTAT_STATISTIC_PADDED_IDENTIFIER % id
        if not isinstance(id, string_types):
            raise TypeError("the value '%s' is not a valid FindStat statistic identifier, nor a FindStat statistic query" % id)
        if FINDSTAT_MAP_SEPARATOR in id:
            return FindStatCompoundStatistic(id)
        if not re.match(FINDSTAT_STATISTIC_REGEXP, id) or id == FINDSTAT_STATISTIC_PADDED_IDENTIFIER % 0:
            raise ValueError("the value '%s' is not a valid FindStat statistic identifier" % id)
        if id not in _all_statistics or _all_statistics[id] is None:
            _all_statistics[id] = self.element_class(self, id)

        return _all_statistics[id]

    def _repr_(self):
        """
        Return the representation of the set of FindStat maps.

        EXAMPLES::

            sage: from sage.databases.findstat import FindStatStatistics
            sage: FindStatStatistics()                                          # optional -- internet
            Set of combinatorial statistics in FindStat

            sage: FindStatStatistics(12)                                        # optional -- internet
            Set of combinatorial statistics with domain Cc0012: Perfect matchings in FindStat
        """
        if self._domain is None:
            return "Set of combinatorial statistics in FindStat"
        return "Set of combinatorial statistics with domain %s in FindStat" % self._domain

    def __iter__(self):
        """
        Return an iterator over all FindStat statistics.

        EXAMPLES::

            sage: from sage.databases.findstat import FindStatStatistics
            sage: next(iter(FindStatStatistics("Permutations")))                # optional -- internet
            St000001: The number of reduced words for a permutation.

        """
        if self._identifiers is None:
            if self._domain is None:
                url = FINDSTAT_API_STATISTICS
            else:
                url = FINDSTAT_API_STATISTICS + "?Domain=%s" % self._domain.id_str()

            self._identifiers = json.load(urlopen(url))["data"]

        for st in self._identifiers:
            yield FindStatStatistic(st)

    Element = FindStatStatistic

class FindStatStatisticQuery(SageObject):
    """
    A class representing a query for FindStat (compound) statistics.
    """
    @staticmethod
    def _data_to_str(data, domain):
        to_str = domain.to_string()
        return "\n".join("\n".join(to_str(element) for element in elements)
                         + "\n====> "
                         + FINDSTAT_VALUE_SEPARATOR.join(str(value)
                                                         for value in values)
                         for elements, values in data)

    def __init__(self, data=None, values_of=None, distribution_of=None,
                 domain=None, all_data=None, function=None,
                 depth=FINDSTAT_DEFAULT_DEPTH,
                 debug=False):
        """Initialize a query for FindStat (compound) statistics.

        INPUT::

        - ``data`` -- (optional), a list of pairs ``(objects,
          values)``, where ``objects`` and ``values`` are all lists
          of the same length, the former are elements in the FindStat
          collection, the latter are integers

        - ``all_data`` -- (optional), a callable or a list in the
          same format as ``data``, which agrees with ``data``, and
          may be used for submission

        - ``values_of`` -- (optional), anything accepted by
          :class:`FindStatCompoundStatistic`

        - ``distribution_of`` -- (optional), anything accepted by
          :class:`FindStatCompoundStatistic`

        - ``domain`` -- (optional), anything accepted by
          :class:`FindStatCollection`

        - ``depth`` -- (optional), the number of maps to apply before
          applying the statistic


        Only one of ``data``, ``values_of`` and ``distribution_of``
        may be provided.  The parameter ``domain`` must be provided
        if and only if ``data`` is provided, or ``values_of`` or
        ``distribution_of`` are given as a function.

        The parameter ``all_data`` is only allowed, if ``data`` is
        provided.

        EXAMPLES::

            sage: from sage.databases.findstat import FindStatStatisticQuery
            sage: data = [[[pi],[pi[0]]] for n in range(1,7) for pi in Permutations(n)] # optional -- internet
            sage: FindStatStatisticQuery(domain=1, data=data, depth=1)          # optional -- internet
            0: St000054 (quality [100, 100])
            1: St000025oMp00127 (quality [100, 100])
            2: St000051oMp00072 with offset -1 (quality [100, 100])
            ...

        """
        self._data = data
        self._all_data = all_data
        self._function = function
        self._values_of = None
        self._distribution_of = None
        self._depth = depth

        if data is not None:
            assert all(param is None for param in [distribution_of, values_of])

            self._domain = FindStatCollection(domain)
            query = {"Domain": self._domain.id_str(),
                     "Data": self._data_to_str(self._data, self._domain)}

        elif distribution_of is not None:
            assert all(param is None for param in [data, all_data, values_of])

            self._distribution_of = FindStatCompoundStatistic(distribution_of)
            query = {"DistributionOf": distribution_of}

        elif values_of is not None:
            assert all(param is None for param in [data, all_data, distribution_of])

            self._values_of = FindStatCompoundStatistic(values_of)
            query = {"ValuesOf": values_of}

        else:
            raise ValueError("incompatible set of parameters: data: %s, distribution_of: %s, values_of: %s" % ((data, distribution_of, values_of)))

        if depth is not None:
            query["Depth"] = depth

        query["fields"] = "MatchingStatistic,Offset,Quality"
        if debug:
            print(query)
        response = requests.post(FINDSTAT_API_STATISTICS, data=query).json()

        if debug:
            print(response)
        if "data" not in response:
            raise ValueError(response["error"])
        result = []
        for match in response["data"]:
            entry = response["included"]["MatchingStatistics"][match]
            result.append(FindStatMatchingStatistic(entry["MatchingStatistic"],
                                                    entry["Offset"],
                                                    entry["Quality"]))
        if self.submittable():
            result.append(FindStatStatistic(self))
        self._result = FancyTuple(result)

    def submittable(self):
        return self._all_data is not None

    def domain(self):
        return self._domain

    def function(self):
        return self._function

    def code(self):
        if self._function is None:
            return ""
        try:
            code = inspect.getsource(self._function)
        except IOError:
            verbose("inspect.getsource could not get code from function provided", caller_name='FindStat')
            code = ""
        return code

    def data(self):
        return self._all_data

    def generating_functions(self, style="polynomial"):
        return _generating_functions_from_dict(self._generating_functions_dict(),
                                               style=style)

    def _generating_functions_dict(self, max_values=FINDSTAT_MAX_VALUES):
        gfs = {} # lvl: vals
        total = min(max_values, FINDSTAT_MAX_VALUES)
        iterator = iter(self.data())
        domain = self.domain()
        levels_with_sizes = domain.levels_with_sizes()
        while total > 0:
            try:
                elts, vals = next(iterator)
            except StopIteration:
                break
            if total < len(elts):
                break
            else:
                lvl = domain.element_level(elts[0])
                if lvl not in levels_with_sizes:
                    continue
                if levels_with_sizes[lvl] > total:
                    # we assume that from now on levels become even larger
                    break
                if not all(domain.element_level(elt) == lvl for elt in elts[1:]):
                    raise ValueError("cannot combine %s into a distribution" % elts)
                gfs[lvl] = gfs.get(lvl, []) + vals
                if levels_with_sizes[lvl] == len(gfs[lvl]):
                    total -= levels_with_sizes[lvl]

        return {lvl: {val: vals.count(val) for val in set(vals)}
                for lvl, vals in gfs.items()
                if levels_with_sizes[lvl] == len(vals)}

    def _repr_(self):
        return repr(self._result)

    def __getitem__(self, idx):
        return self._result[idx]

class FindStatCompoundStatistic(Element):
    def __init__(self, id_str):
        """
        TESTS::

            sage: findstat("St000001oMp00127")                                  # optional -- internet
            Traceback (most recent call last):
            ...
            ValueError: the statistic St000001: The number of reduced words for a permutation. cannot be composed with the map Mp00127

        """
        if isinstance(id_str, Integer):
            id_str = FINDSTAT_STATISTIC_PADDED_IDENTIFIER % id_str
        self._id_str = id_str
        composition = id_str.partition("o")
        self._statistic = FindStatStatistic(composition[0])
        self._maps = FindStatCompoundMap(composition[2])
        if self._maps.codomain() is not None and self._maps.codomain() != self._statistic.domain():
            raise ValueError("the statistic %s cannot be composed with the map %s" % (self._statistic, self._maps))
        Element.__init__(self, FindStatStatistics())

    def domain(self):
        return self._maps.domain()

    def __call__(self, elt):
        return self.statistic()(self.maps()(elt))

    def id_str(self):
        return self._id_str

    def _repr_(self):
        return self.id_str()

    def statistic(self):
        return self._statistic

    def maps(self):
        return self._maps

    def browse(self):
        r"""
        Open the FindStat web page of the compound statistic in a browser.

        EXAMPLES::

            sage: FindStatCompoundStatistic("Mp00087oMp00058").browse()         # optional -- webbrowser
        """
        webbrowser.open(FINDSTAT_URL_STATISTICS + self.id_str())


class FindStatMatchingStatistic(FindStatCompoundStatistic):
    def __init__(self, matching_statistic, offset, quality):
        self._quality = quality
        self._offset = offset
        super(FindStatMatchingStatistic, self).__init__(matching_statistic)

    def _repr_(self):
        if self._offset:
            return "%s with offset %s (quality %s)" % (self.id_str(), self._offset, self._quality)
        return "%s (quality %s)" % (self.id_str(), self._quality)

    def offset(self):
        return self._offset

    def quality(self):
        return self._quality[:]


class FindStatMapQuery(SageObject):
    @staticmethod
    def _data_to_str(data, domain, codomain):
        to_str_dom = domain.to_string()
        to_str_codom = codomain.to_string()
        return "\n".join("\n".join(to_str_dom(element) for element in elements)
                         + "\n====> "
                         + FINDSTAT_VALUE_SEPARATOR.join(to_str_codom(value)
                                                         for value in values)
                         for elements, values in data)

    def __init__(self, data=None, values_of=None, distribution_of=None,
                 domain=None, codomain=None, all_data=None, function=None,
                 depth=FINDSTAT_DEFAULT_DEPTH,
                 debug=False):
        """
        EXAMPLES::

            sage: from sage.databases.findstat import FindStatMapQuery
            sage: FindStatMapQuery(domain=1, codomain=1, data=[[[pi],[pi.reverse().inverse()]] for pi in Permutations(4)]) # optional -- internet
            0: Mp00066oMp00064 (quality [100])
        """
        self._data = data
        self._all_data = all_data
        self._function = function
        self._values_of = None
        self._distribution_of = None
        self._depth = depth

        if data is not None:
            assert all(param is None for param in [distribution_of, values_of])

            self._domain = FindStatCollection(domain)
            self._codomain = FindStatCollection(codomain)
            query = {"Domain": self._domain.id_str(),
                     "Codomain": self._codomain.id_str(),
                     "Data": self._data_to_str(self._data, self._domain, self._codomain)}

        elif distribution_of is not None:
            assert all(param is None for param in [data, all_data, values_of])

            self._distribution_of = FindStatCompoundMap(distribution_of)
            query = {"DistributionOf": distribution_of}

        elif values_of is not None:
            assert all(param is None for param in [data, all_data, distribution_of])

            self._values_of = FindStatCompoundMap(values_of)
            query = {"ValuesOf": values_of}

        else:
            raise ValueError("incompatible set of parameters: data: %s, distribution_of: %s, values_of: %s" % ((data, distribution_of, values_of)))

        if depth is not None:
            query["Depth"] = depth

        query["fields"] = "MatchingMap,Quality"
        if debug:
            print(query)
        response = requests.post(FINDSTAT_API_MAPS, data=query).json()

        if debug:
            print(response)
        result = []
        if "data" not in response:
            raise ValueError(response["error"])
        for match in response["data"]:
            entry = response["included"]["MatchingMaps"][match]
            result.append(FindStatMatchingMap(entry["MatchingMap"],
                                              entry["Quality"]))
        self._result = FancyTuple(result)

    def domain(self):
        return self._domain

    def codomain(self):
        return self._codomain

    def function(self):
        return self._function

    def _repr_(self):
        return repr(self._result)

    def __getitem__(self, idx):
        return self._result[idx]

class FindStatCompoundMap(Element):
    def __init__(self, id_str):
        """
        TESTS::

            sage: findmap("Mp00146oMp00127").domain()                           # optional -- internet
            Cc0001: Permutations

            sage: findmap("Mp00146oMp00127").codomain()                         # optional -- internet
            Cc0012: Perfect matchings

            sage: findmap("Mp00127oMp00146")                                    # optional -- internet
            Traceback (most recent call last):
            ...
            ValueError: the sequence of maps [Mp00146: to tunnel matching, Mp00127: left-to-right-maxima to Dyck path] cannot be composed
        """

        if isinstance(id_str, Integer):
            id_str = FINDSTAT_MAP_PADDED_IDENTIFIER % id_str
        if id_str == "":
            self._id_str = "id"
            self._maps = []
        else:
            self._id_str = id_str
            self._maps = [FindStatMap(m) for m in id_str.split("o")][::-1]
            if not all(self._maps[i].codomain() == self._maps[i+1].domain()
                       for i in range(len(self._maps)-1)):
                raise ValueError("the sequence of maps %s cannot be composed" % self._maps)
        Element.__init__(self, FindStatMaps())

    def __getitem__(self, i):
        return self._maps[i]

    def domain(self):
        if self._maps:
            return self._maps[0].domain()

    def codomain(self):
        if self._maps:
            return self._maps[-1].codomain()

    def __call__(self, elt):
        for m in self.maps():
            elt = m(elt)
        return elt

    def id_str(self):
        return self._id_str

    def _repr_(self):
        return self.id_str()

    def maps(self):
        return self._maps

    def browse(self):
        r"""
        Open the FindStat web page of the compound map in a browser.

        EXAMPLES::

            sage: findmap(41).browse()                                          # optional -- webbrowser
        """
        webbrowser.open(FINDSTAT_URL_MAPS + self.id_str())


class FindStatMatchingMap(FindStatCompoundMap):
    def __init__(self, matching_map, quality):
        self._quality = quality
        super(FindStatMatchingMap, self).__init__(matching_map)

    def _repr_(self):
        return "%s (quality %s)" % (self.id_str(), self._quality)

    def quality(self):
        return self._quality[:]

@add_metaclass(InheritComparisonClasscallMetaclass)
class FindStatMap(Element, FindStatFunction):
    r"""
    A FindStat map.

    :class:`FindStatMap` is a class representing a combinatorial
    map available in the FindStat database.

    The result of a :class:`findstat<FindStat>` query contains a
    (possibly empty) list of such maps.  This class provides methods
    to inspect various properties of these maps, in particular
    :meth:`code`.

    INPUT:

    - a string or an integer representing the FindStat identifier of
      the map.

    EXAMPLES::

        sage: from sage.databases.findstat import FindStatMap
        sage: FindStatMap(71)                                                   # optional -- internet
        Mp00071: descent composition

    .. SEEALSO::

        :class:`FindStatMaps`

    """
    @staticmethod
    def __classcall_private__(cls, entry):
        """
        Retrieve a map from the database.

        TESTS::

            sage: from sage.databases.findstat import FindStatMap
            sage: FindStatMap("abcdefgh")                                       # optional -- internet
            Traceback (most recent call last):
            ...
            ValueError: the value 'abcdefgh' is not a valid FindStat map identifier
        """
        return FindStatMaps()(entry)

    def __init__(self, parent, id):
        """
        Initialize the map.

        This should only be called in
        :meth:`FindStatMaps()._element_constructor_` via
        :meth:`FindStatMaps().element_class`.

        INPUT:

        - ``parent`` -- :class:`FindStatMaps`.

        - ``entry`` -- a dictionary containing the properties of the
          map, such as its name, code, and so on.

        TESTS::

            sage: from sage.databases.findstat import FindStatMap
            sage: FindStatMap(62).parent()                                      # optional -- internet
            Set of combinatorial maps used by FindStat

        """
        self._function = False # determines whether FindStat code may be executed
        self._modified_data = None
        self._modified = False
        self._id = id
        Element.__init__(self, parent)

    def __reduce__(self):
        """Return a function and its arguments needed to create this
        map.

        TESTS::

            sage: from sage.databases.findstat import FindStatMap
            sage: c = FindStatMap(62)                                           # optional -- internet
            sage: loads(dumps(c)) == c                                          # optional -- internet
            True

        """
        return (FindStatMap, (self.id(),))

    def _richcmp_(self, other, op):
        return richcmp(self.id(), other.id(), op)

    # WARNING: before modifying self._modified_data this method has
    # to be called, it is cached because it needs to be executed only
    # once
    def _fetch_data(self):
        if self._modified_data is not None:
            return
        fields = "Bibliography,Codomain,Description,Domain,Name,Properties,References,SageCode,SageName"
        fields_Bibliography = "Author,Title"
        url = (FINDSTAT_API_MAPS + self._id
               + "?fields=" + fields
               + "&fields[Bibliography]=" + fields_Bibliography)
        included = json.load(urlopen(url))["included"]
        # slightly simplify the representation
        data = {key: val for key, val in included["Maps"][self._id].items()}
        # we replace the list of identifiers in Bibliography with the dictionary
        data["Bibliography"] = included["References"]
        self._modified_data = data

    def browse(self):
        r"""
        Open the FindStat web page of the map in a browser.

        EXAMPLES::

            sage: from sage.databases.findstat import FindStatMap
            sage: FindStatMap(62).browse()                                      # optional -- webbrowser
        """
        if self.is_new():
            self.submit()
        else:
            webbrowser.open(FINDSTAT_URL_MAPS + self.id_str())

    def submit(self):
        r"""
        Open the FindStat web page for editing the map in a browser.
        """
        args = dict()
        args["OriginalMap"]        = self.id_str()
        args["Domain"]             = str(self.domain().id_str())
        args["Codomain"]           = str(self.codomain().id_str())
        args["Name"]               = self.name()
        args["Description"]        = self.description()
        args["References"]         = self.references_raw()
        args["Properties"]         = self.properties_raw()
        args["SageCode"]           = self.sage_code()
        args["SageName"]           = self.code_name()
        args["CurrentAuthor"]      = findstat._user_name
        args["CurrentEmail"]       = findstat._user_email

        f = tempfile.NamedTemporaryFile(suffix='.html', delete=False)
        verbose("Created temporary file %s" % f.name, caller_name='FindStat')
        f.write(FINDSTAT_POST_HEADER)
        if self.id() == 0:
            f.write(FINDSTAT_NEWMAP_FORM_HEADER % FINDSTAT_URL_NEW_MAP)
        else:
            f.write(FINDSTAT_NEWMAP_FORM_HEADER % (FINDSTAT_URL_EDIT_MAP+self.id_str()))
        for key, value in iteritems(args):
            if value:
                verbose("writing argument %s" % key, caller_name='FindStat')
                value_encoded = cgi.escape(str(value), quote=True)
                verbose("%s" % value_encoded, caller_name='FindStat')
                f.write((FINDSTAT_FORM_FORMAT %(key, value_encoded)))
            else:
                verbose("skipping argument %s because it is empty" % key, caller_name='FindStat')
        f.write(FINDSTAT_FORM_FOOTER)
        f.close()
        verbose("Opening file with webbrowser", caller_name='FindStat')
        webbrowser.open(f.name)

    # editing and submitting is really the same thing
    edit = submit

    def __hash__(self):
        """
        Returns a hash value for the map.

        EXAMPLES::

            sage: from sage.databases.findstat import FindStatMaps
            sage: sorted(list(FindStatMaps(codomain=17)))                       # optional -- internet
            [Mp00003: rotate counterclockwise,
             Mp00004: rotate clockwise,
             Mp00005: transpose,
             Mp00006: gyration,
             Mp00035: to alternating sign matrix,
             Mp00063: to alternating sign matrix,
             Mp00137: to symmetric ASM]
        """
        return self.id()

    def domain(self):
        r"""
        Return the FindStat collection which is the domain of the map.

        OUTPUT:

        The domain of the map as a :class:`FindStatCollection`.

        EXAMPLES::

            sage: from sage.databases.findstat import FindStatMap               # optional -- internet
            sage: FindStatMap(71).domain()                                      # optional -- internet
            Cc0001: Permutations
        """
        return FindStatCollection(self._data()["Domain"])

    def codomain(self):
        r"""
        Return the FindStat collection which is the codomain of the map.

        OUTPUT:

        The codomain of the map as a :class:`FindStatCollection`.

        EXAMPLES::

            sage: from sage.databases.findstat import FindStatMap               # optional -- internet
            sage: FindStatMap(71).codomain()                                    # optional -- internet
            Cc0006: Integer compositions
        """
        return FindStatCollection(self._data()["Codomain"])

    def properties_raw(self):
        r"""
        Return the properties of the map.

        OUTPUT:

        The properties as a string.

        EXAMPLES::

            sage: from sage.databases.findstat import FindStatMap               # optional -- internet
            sage: FindStatMap(71).properties_raw()                              # optional -- internet
            u'surjective, graded'
        """
        return self._data()["Properties"]

    def set_properties_raw(self, value):
        r"""
        Set the properties of the map.

        EXAMPLES::

            sage: from sage.databases.findstat import FindStatMap               # optional -- internet
            sage: FindStatMap(71).set_properties_raw(u'surjective')             # optional -- internet
            sage: FindStatMap(71).properties_raw()                              # optional -- internet
            u'surjective'
            sage: FindStatMap(71)                                               # optional -- internet
            Mp00071(modified): descent composition
            sage: from sage.databases.findstat import _all_maps                 # optional -- internet
            sage: del _all_maps['Mp00071']                                      # optional -- internet
            sage: FindStatMap(71)                                               # optional -- internet
            Mp00071: descent composition
        """
        if value != self.properties_raw():
            self._modified = True
            self._modified_data["Properties"] = value

    def code_name(self):
        r"""
        Return the name of the function defined by :meth:`code`.

        OUTPUT:

        A string.

        EXAMPLES::

            sage: from sage.databases.findstat import FindStatMap               # optional -- internet
            sage: print(FindStatMap(71).code_name())                            # optional -- internet
            descents_composition
        """
        return self._data()["SageName"]

    def set_code_name(self, value):
        r"""
        Set the sage code name of the map.

        INPUT:

        - a string -- for statistics, this is the name of the
          statistic followed by its description on a separate line.

        This information is used when submitting the statistic or map with
        :meth:`submit`.

        """
        if value != self.code_name():
            self._modified = True
            self._modified_data["SageName"] = value

_all_maps = {}

class FindStatMaps(UniqueRepresentation, Parent):
    r"""
    The class of FindStat maps.

    The elements of this class are combinatorial maps currently in
    FindStat.

    EXAMPLES:

    We can print a sample map, for each domain and codomain::

        sage: from sage.databases.findstat import FindStatCollections, FindStatMap, FindStatMaps
        sage: ccs = sorted(FindStatCollections())[:3]                           # optional -- internet
        sage: for cc_dom in ccs:                                                # optional -- internet, random
        ....:     for cc_codom in ccs:
        ....:         print(cc_dom.name(style="plural") + " -> " + cc_codom.name(style="plural"))
        ....:         try:
        ....:             print("    " + next(iter(FindStatMaps(cc_dom, cc_codom))).name())
        ....:         except StopIteration:
        ....:             pass
        Permutations -> Permutations
            Lehmer-code to major-code bijection
        Permutations -> Integer partitions
            Robinson-Schensted tableau shape
        Permutations -> Dyck paths
            left-to-right-maxima to Dyck path
        Integer partitions -> Permutations
        Integer partitions -> Integer partitions
            conjugate
        Integer partitions -> Dyck paths
            to Dyck path
        Dyck paths -> Permutations
            to non-crossing permutation
        Dyck paths -> Integer partitions
            to partition
        Dyck paths -> Dyck paths
            reverse

    """
    def __init__(self, domain=None, codomain=None):
        """
        TESTS::

            sage: from sage.databases.findstat import FindStatMaps
            sage: M = FindStatMaps()                                            # optional -- internet
            sage: TestSuite(M).run()                                            # optional -- internet
        """
        if domain is None:
            self._domain = None
        else:
            self._domain = FindStatCollection(domain)
        if codomain is None:
            self._codomain = None
        else:
            self._codomain = FindStatCollection(codomain)
        self._identifiers = None
        Parent.__init__(self, category=Sets())

    def _element_constructor_(self, id):
        """Initialize a FindStat map.

        INPUT:

        - ``id`` -- a string containing the FindStat identifier of
          the map, or an integer giving its id

        EXAMPLES::

            sage: from sage.databases.findstat import FindStatMap
            sage: FindStatMap(71)                                               # optional -- internet
            Mp00071: descent composition

        """
        if isinstance(id, self.Element):
            return id
        if isinstance(id, (int, Integer)):
            id = FINDSTAT_MAP_PADDED_IDENTIFIER % id
        if not isinstance(id, string_types):
            raise TypeError("the value %s is neither an integer nor a string" % id)
        if FINDSTAT_MAP_SEPARATOR in id:
            return FindStatCompoundMap(id)
        if not re.match(FINDSTAT_MAP_REGEXP, id) or id == FINDSTAT_MAP_PADDED_IDENTIFIER % 0:
            raise ValueError("the value '%s' is not a valid FindStat map identifier" % id)
        if id not in _all_maps or _all_maps[id] is None:
            _all_maps[id] = self.element_class(self, id)

        return _all_maps[id]

    def _repr_(self):
        """
        Return the representation of the set of FindStat maps.

        EXAMPLES::

            sage: from sage.databases.findstat import FindStatMaps
            sage: FindStatMaps()                                                # optional -- internet
            Set of combinatorial maps used by FindStat
        """
        if self._domain is None:
            text = []
        else:
            text = ["domain %s" % self._domain]
        if self._codomain is not None:
            text.append("codomain %s" % self._codomain)

        if text:
            return "Set of combinatorial maps with " + " and ".join(text) + " used by FindStat"
        return "Set of combinatorial maps used by FindStat"

    def __iter__(self):
        """
        Return an iterator over all FindStat maps.

        EXAMPLES::

            sage: from sage.databases.findstat import FindStatMaps
            sage: next(iter(FindStatMaps(domain=5)))                            # optional -- internet
            Mp00023: to non-crossing permutation

        """
        if self._identifiers is None:
            query = []
            if self._domain is not None:
                query.append("Domain=%s" % self._domain.id_str())
            if self._codomain is not None:
                query.append("Codomain=%s" % self._codomain.id_str())
            if query:
                url = FINDSTAT_API_MAPS + "?" + "&".join(query)
            else:
                url = FINDSTAT_API_MAPS
            self._identifiers = json.load(urlopen(url))["data"]

        for mp in self._identifiers:
            yield FindStatMap(mp)

    Element = FindStatMap

# helper for generation of CartanTypes
def _finite_irreducible_cartan_types_by_rank(n):
    """
    Return the Cartan types of rank n.

    INPUT:

    - n -- an integer.

    OUTPUT:

    The list of Cartan types of rank n.

    TESTS::

        sage: from sage.databases.findstat import _finite_irreducible_cartan_types_by_rank
        sage: _finite_irreducible_cartan_types_by_rank(2)
        [['A', 2], ['B', 2], ['G', 2]]
    """
    cartan_types      = [ CartanType(['A',n]) ]
    if n >= 2:
        cartan_types += [ CartanType(['B',n]) ]
    if n >= 3:
        cartan_types += [ CartanType(['C',n]) ]
    if n >= 4:
        cartan_types += [ CartanType(['D',n]) ]
    if n in [6,7,8]:
        cartan_types += [ CartanType(['E',n]) ]
    if n == 4:
        cartan_types += [ CartanType(['F',n]) ]
    if n == 2:
        cartan_types += [ CartanType(['G',n]) ]
    return cartan_types


@add_metaclass(InheritComparisonClasscallMetaclass)
class FindStatCollection(Element):
    r"""
    A FindStat collection.

    :class:`FindStatCollection` is a class representing a
    combinatorial collection available in the FindStat database.

    Its main use is to allow easy specification of the combinatorial
    collection when using :class:`findstat<FindStat>`.  It also
    provides methods to quickly access its FindStat web page
    (:meth:`browse`), check whether a particular element is actually
    in the range considered by FindStat (:meth:`in_range`), etc.

    INPUT:

    One of the following:

    - a string eg. 'Dyck paths' or 'DyckPaths', case-insensitive, or

    - an integer designating the FindStat id of the collection, or

    - a sage object belonging to a collection, or

    - an iterable producing a sage object belonging to a collection.

    EXAMPLES::

        sage: from sage.databases.findstat import FindStatCollection
        sage: FindStatCollection("Dyck paths")                                  # optional -- internet
        Cc0005: Dyck paths

        sage: FindStatCollection(5)                                             # optional -- internet
        Cc0005: Dyck paths

        sage: FindStatCollection(DyckWord([1,0,1,0]))                           # optional -- internet
        Cc0005: Dyck paths

        sage: FindStatCollection(DyckWords(2))                                  # optional -- internet
        a subset of Cc0005: Dyck paths

        sage: FindStatCollection(DyckWords)                                     # optional -- internet
        Cc0005: Dyck paths

    .. SEEALSO::

        :class:`FindStatCollections`
    """
    @staticmethod
    def __classcall_private__(cls, entry):
        """
        Retrieve a collection from the database.

        TESTS::

            sage: from sage.databases.findstat import FindStatCollection
            sage: FindStatCollection(0)                                         # optional -- internet
            Traceback (most recent call last):
            ...
            ValueError: Could not find FindStat collection for 0.
        """
        return FindStatCollections()(entry)

    def __init__(self, parent, id, data, sageconstructor_overridden):
        """Initialize the collection.

        This should only be called in
        :meth:`FindStatCollections()._element_constructor_` via
        :meth:`FindStatCollections().element_class`.

        INPUT:

        - ``parent`` -- :class:`FindStatCollections`.

        - ``id`` -- the (padded) FindStat identifier of the collection.

        - ``data`` -- a dictionary containing the properties of the
          collection, such as its name, the corresponding class in
          sage, and so on.

        - ``sageconstructor_overridden`` -- either ``None`` or an
          iterable which yields a subset of the elements of the
          collection.

        TESTS::

            sage: from sage.databases.findstat import FindStatCollection
            sage: FindStatCollection(5).parent()                                # optional -- internet
            Set of combinatorial collections used by FindStat

        """
        self._id = id
        self._data = data
        self._sageconstructor_overridden = sageconstructor_overridden

        Element.__init__(self, parent)

    def __reduce__(self):
        """Return a function and its arguments needed to create this
        collection.

        TESTS::

            sage: from sage.databases.findstat import FindStatCollection
            sage: c = FindStatCollection("Permutations")                        # optional -- internet
            sage: loads(dumps(c)) == c                                          # optional -- internet
            True

        """
        return FindStatCollection, (self.id(),)

    def _richcmp_(self, other, op):
        """
        TESTS::

            sage: from sage.databases.findstat import FindStatCollection, FindStatCollections
            sage: FindStatCollection("Permutations") == FindStatCollection("Permutations")          # optional -- internet
            True

            sage: FindStatCollection("Permutations") == FindStatCollection("Integer Partitions")    # optional -- internet
            False

            sage: FindStatCollection("Permutations") != FindStatCollection("Permutations")          # optional -- internet
            False

            sage: FindStatCollection("Permutations") != FindStatCollection("Integer Partitions")    # optional -- internet
            True

            sage: FindStatCollection("Permutations") == 1                                           # optional -- internet
            False

            sage: FindStatCollection("Permutations") != 1                                           # optional -- internet
            True

            sage: sorted(c for c in FindStatCollections())[0]                                       # optional -- internet
            Cc0001: Permutations

        It is not clear what we want here, but equality may be better::

            sage: FindStatCollection(graphs(3)) == FindStatCollection("Graphs")                     # optional -- internet
            True
        """
        return richcmp(self.id(), other.id(), op)

    def __hash__(self):
        """
        Returns a hash value for the collection.

        EXAMPLES::

            sage: from sage.databases.findstat import FindStatCollections
            sage: set(c for c in FindStatCollections())                         # optional -- internet
            {Cc0001: Permutations,
             Cc0002: Integer partitions,
             ...

        """
        return self.id()

    def is_supported(self):
        """
        Check whether the collection is fully supported by the interface.

        EXAMPLES::

            sage: from sage.databases.findstat import FindStatCollection
            sage: FindStatCollection(1).is_supported()                          # optional -- internet
            True

            sage: FindStatCollection(24).is_supported()                         # optional -- internet, random
            False

        """
        return "Code" in self._data

    def elements_on_level(self, level):
        """
        Return an iterable over the elements on the given level.
        """
        return self._data["Code"].elements_on_level(level)

    def element_level(self, element):
        """
        Return the level of an element.
        """
        return self._data["Code"].element_level(element)

    def is_element(self, element):
        """
        Return whether the element belongs to this collection.
        """
        return self._data["Code"].is_element(element)

    def levels_with_sizes(self):
        """
        Return a dictionary from levels to level sizes.
        """
        return self._data["LevelsWithSizes"]

    def in_range(self, element):
        r"""
        Check whether an element of the collection is in FindStat's precomputed range.

        INPUT:

        - ``element`` -- a sage object that belongs to the collection.

        OUTPUT:

        ``True``, if ``element`` is used by the FindStat search
        engine, and ``False`` if it is ignored.

        EXAMPLES::

            sage: from sage.databases.findstat import FindStatCollection
            sage: c = FindStatCollection("GelfandTsetlinPatterns")              # optional -- internet
            sage: c.in_range(GelfandTsetlinPattern([[2, 1], [1]]))              # optional -- internet
            True
            sage: c.in_range(GelfandTsetlinPattern([[3, 1], [1]]))              # optional -- internet
            True
            sage: c.in_range(GelfandTsetlinPattern([[7, 1], [1]]))              # optional -- internet, random
            False


        TESTS::

            sage: from sage.databases.findstat import FindStatCollections
            sage: l = FindStatCollections()                                     # optional -- internet
            sage: long = [9, 12, 14, 20]
            sage: for c in sorted(l):                                           # optional -- internet, random
            ....:     if c.id() not in long and c.is_supported():
            ....:         f = c.first_terms(lambda x: 1)
            ....:         print("{} {} {}".format(c, len(list(f)), all(c.in_range(e) for e, _ in f)))
            ....:
            Cc0001: Permutations 5913 True
            Cc0002: Integer partitions 1211 True
            Cc0005: Dyck paths 2055 True
            Cc0006: Integer compositions 1023 True
            Cc0007: Standard tableaux 1115 True
            Cc0010: Binary trees 2055 True
            Cc0013: Cores 100 True
            Cc0017: Alternating sign matrices 7917 True
            Cc0018: Gelfand-Tsetlin patterns 1409 True
            Cc0019: Semistandard tableaux 2374 True
            Cc0021: Ordered trees 2056 True
            Cc0022: Finite Cartan types 31 True
            Cc0023: Parking functions 18248 True
            Cc0024: Binary words 1022 True
            Cc0028: Skew partitions 1250 True

        """
        return self._data["Code"].element_level(element) in self._data["LevelsWithSizes"]

    def first_terms(self, function, level=None):
        r"""
        Compute the first few terms of the given statistic or map,
        restricted to the given level.

        INPUT:

        - ``function`` -- a callable

        - ``level`` -- (optional), the level to restrict to

        OUTPUT:

        A lazy list of pairs of the form ``(object, value)``.

        EXAMPLES::

            sage: from sage.databases.findstat import FindStatCollection
            sage: c = FindStatCollection("GelfandTsetlinPatterns")              # optional -- internet
            sage: c.first_terms(lambda x: 1)[:10].list()                        # optional -- internet, random
            [([[1, 0], [0]], 1),
             ([[1, 0], [1]], 1),
             ([[2, 0], [0]], 1),
             ([[2, 0], [1]], 1),
             ([[2, 0], [2]], 1),
             ([[1, 1], [1]], 1),
             ([[1, 0, 0], [0, 0], [0]], 1),
             ([[1, 0, 0], [1, 0], [0]], 1),
             ([[1, 0, 0], [1, 0], [1]], 1),
             ([[3, 0], [0]], 1)]

        """
        if self._sageconstructor_overridden is None:
            if level is None:
                g = (x
                     for level in self._data["LevelsWithSizes"]
                     for x in self._data["Code"].elements_on_level(level))
            else:
                g = (x for x in self._data["Code"].elements_on_level(level))
        else:
            if level is None:
                g = self._sageconstructor_overridden
            else:
                g = (x for x in self._sageconstructor_overridden
                     if self.element_level(x) == level)

        return lazy_list(((x, function(x)) for x in g))

    def id(self):
        r"""
        Return the FindStat identifier of the collection.

        OUTPUT:

        The FindStat identifier of the collection as an integer.

        EXAMPLES::

            sage: from sage.databases.findstat import FindStatCollection
            sage: c = FindStatCollection("GelfandTsetlinPatterns")              # optional -- internet
            sage: c.id()                                                        # optional -- internet
            18
        """
        return int(self._id[2:])

    def id_str(self):
        r"""
        Return the FindStat identifier of the collection.

        OUTPUT:

        The FindStat identifier of the collection as a string.

        EXAMPLES::

            sage: from sage.databases.findstat import FindStatCollection
            sage: c = FindStatCollection("GelfandTsetlinPatterns")              # optional -- internet
            sage: c.id_str()                                                    # optional -- internet
            u'Cc0018'
        """
        return self._id

    def browse(self):
        r"""
        Open the FindStat web page of the collection in a browser.

        EXAMPLES::

            sage: from sage.databases.findstat import FindStatCollection
            sage: FindStatCollection("Permutations").browse()                   # optional -- webbrowser
        """
        webbrowser.open(FINDSTAT_URL_COLLECTIONS + self._id)

    def to_string(self):
        r"""
        Return a function that returns a FindStat representation given an
        object.

        OUTPUT:

        The function that produces the string representation as
        needed by the FindStat search webpage.

        EXAMPLES::

            sage: from sage.databases.findstat import FindStatCollection
            sage: p = Poset((range(3), [[0, 1], [1, 2]]))                       # optional -- internet
            sage: c = FindStatCollection("Posets")                              # optional -- internet
            sage: c.to_string()(p)                                              # optional -- internet
            '([(0, 1), (1, 2)], 3)'

        """
        return self._data["Code"].element_to_string

    def from_string(self):
        r"""
        Return a function that returns the object given a FindStat
        representation.

        OUTPUT:

        The function that produces the sage object given its FindStat
        representation as a string.

        EXAMPLES::

            sage: from sage.databases.findstat import FindStatCollection
            sage: c = FindStatCollection("Posets")                              # optional -- internet
            sage: p = c.from_string()('([(0, 2), (2, 1)], 3)')                  # optional -- internet
            sage: p.cover_relations()                                           # optional -- internet
            [[0, 2], [2, 1]]

            sage: c = FindStatCollection("Binary Words")                        # optional -- internet
            sage: w = c.from_string()('010101')                                 # optional -- internet
            sage: w in c._data["Code"].elements_on_level(6)                     # optional -- internet
            True
        """
        return self._data["Code"].string_to_element

    def _repr_(self):
        r"""
        Return the representation of the FindStat collection.

        OUTPUT:

        The representation, including the identifier and the name.

        EXAMPLES::

            sage: from sage.databases.findstat import FindStatCollection
            sage: FindStatCollection("Binary trees")                            # optional -- internet
            Cc0010: Binary trees
        """
        if self._sageconstructor_overridden:
            return "a subset of %s: %s" % (self.id_str(), self._data["NamePlural"])
        return "%s: %s" % (self.id_str(), self._data["NamePlural"])

    def name(self, style="singular"):
        r"""
        Return the name of the FindStat collection.

        INPUT:

        - a string -- (default:"singular") can be
          "singular", or "plural".

        OUTPUT:

        The name of the FindStat collection, in singular or in plural.

        EXAMPLES::

            sage: from sage.databases.findstat import FindStatCollection
            sage: FindStatCollection("Binary trees").name()                     # optional -- internet
            u'Binary tree'

            sage: FindStatCollection("Binary trees").name(style="plural")       # optional -- internet
            u'Binary trees'
        """
        if style == "singular":
            return self._data["Name"]
        elif style == "plural":
            return self._data["NamePlural"]
        else:
            raise ValueError("Argument 'style' (=%s) must be 'singular' or 'plural'."%style)

from collections import namedtuple
SupportedFindStatCollection = namedtuple("SupportedFindStatCollection",
                                         ["string_to_element",
                                          "element_to_string",
                                          "elements_on_level", # return all elements on given level
                                          "element_level",     # return level of a given element
                                          "is_element"]) # return whether element is member of this collection (and, ideally, of no other collection)

SupportedFindStatCollections = {
    "Permutations":
    SupportedFindStatCollection(lambda x: Permutation(literal_eval(x)),
                                str,
                                Permutations,
                                lambda x: x.size(),
                                lambda x: isinstance(x, Permutation)),
    "BinaryWords":
    SupportedFindStatCollection(lambda x: Word((int(e) for e in str(x)), alphabet=[0,1]),
                                str,
                                lambda x: Words([0,1], length=x),
                                lambda x: x.length(),
                                lambda x: isinstance(x, Word_class)),

    "AlternatingSignMatrices":
    SupportedFindStatCollection(lambda x: AlternatingSignMatrix(literal_eval(x)),
                                lambda x: str(list(map(list, x.to_matrix().rows()))),
                                AlternatingSignMatrices,
                                lambda x: x.to_matrix().nrows(),
                                lambda x: isinstance(x, AlternatingSignMatrix)),
    "BinaryTrees":
    SupportedFindStatCollection(lambda x: BinaryTree(str(x)),
                                str,
                                BinaryTrees,
                                lambda x: x.node_number(),
                                lambda x: isinstance(x, BinaryTree)),
    "Cores":
    SupportedFindStatCollection(lambda x: (lambda pi, k: Core(pi, k))(*literal_eval(x)),
                                lambda X: "( " + X._repr_() + ", " + str(X.k()) + " )",
                                lambda x: Cores(x[1], x[0]),
                                lambda x: (x.length(), x.k()),
                                lambda x: isinstance(x, Core)),
    "DyckPaths":
    SupportedFindStatCollection(lambda x: DyckWord(literal_eval(x)),
                                lambda x: str(list(DyckWord(x))),
                                DyckWords,
                                lambda x: x.semilength(),
                                lambda x: isinstance(x, DyckWord)),
    "FiniteCartanTypes":
    SupportedFindStatCollection(lambda x: CartanType(*literal_eval(str(x))),
                                str,
                                _finite_irreducible_cartan_types_by_rank,
                                lambda x: x.rank(),
                                lambda x: isinstance(x, CartanType_abstract)),
    "GelfandTsetlinPatterns":
    SupportedFindStatCollection(lambda x: GelfandTsetlinPattern(literal_eval(x)),
                                str,
                                lambda x: (P
                                           for la in Partitions(x[1], max_length=x[0])
                                           for P in GelfandTsetlinPatterns(top_row=la + [0]*(x[0]-len(la)))),
                                lambda x: (len(x[0]), sum(x[0])),
                                lambda x: (x == GelfandTsetlinPatterns
                                           or isinstance(x, GelfandTsetlinPattern))),
    "Graphs":
    SupportedFindStatCollection(lambda x: (lambda E, V: Graph([list(range(V)),
                                                               lambda i,j: (i,j) in E or (j,i) in E],
                                                              immutable=True))(*literal_eval(x)),
                                lambda X: str((sorted(X.edges(False)), X.num_verts())),
                                graphs,
                                lambda x: x.num_verts(),
                                lambda x: isinstance(x, Graph)),
    "IntegerCompositions":
    SupportedFindStatCollection(lambda x: Composition(literal_eval(x)),
                                str,
                                Compositions,
                                lambda x: x.size(),
                                lambda x: isinstance(x, Composition)),
    "IntegerPartitions":
    SupportedFindStatCollection(lambda x: Partition(literal_eval(x)),
                                str,
                                Partitions,
                                lambda x: x.size(),
                                lambda x: isinstance(x, Partition)),
    "OrderedTrees":
    SupportedFindStatCollection(lambda x: OrderedTree(literal_eval(x)),
                                str,
                                OrderedTrees,
                                lambda x: x.node_number(),
                                lambda x: isinstance(x, OrderedTree)),
    "ParkingFunctions":
    SupportedFindStatCollection(lambda x: ParkingFunction(literal_eval(x)),
                                str,
                                ParkingFunctions,
                                lambda x: len(x),
                                lambda x: isinstance(x, ParkingFunction_class)),
    "PerfectMatchings":
    SupportedFindStatCollection(lambda x: PerfectMatching(literal_eval(x)),
                                str,
                                PerfectMatchings,
                                lambda x: x.size(),
                                lambda x: isinstance(x, PerfectMatching)),
    "Posets":
    SupportedFindStatCollection(lambda x: (lambda R, E: Poset((list(range(E)), R)))(*literal_eval(x)),
                                lambda X: str((sorted(X._hasse_diagram.cover_relations()),
                                               len(X._hasse_diagram.vertices()))),
                                Posets,
                                lambda x: x.cardinality(),
                                lambda x: isinstance(x, FinitePoset)),
    "SemistandardTableaux":
    SupportedFindStatCollection(lambda x: SemistandardTableau(literal_eval(x)),
                                str,
                                lambda x: (T for T in SemistandardTableaux(size=x[0], max_entry=x[1])
                                           if max(T.entries()) == x[1]),
                                lambda x: (x.size(), max(x.entries())),
                                lambda x: isinstance(x, SemistandardTableau)),
    "SetPartitions":
    SupportedFindStatCollection(lambda x: SetPartition(literal_eval(x.replace('{','[').replace('}',']'))),
                                str,
                                SetPartitions,
                                lambda x: x.size(),
                                lambda x: isinstance(x, SetPartition)),
    "StandardTableaux":
    SupportedFindStatCollection(lambda x: StandardTableau(literal_eval(x)),
                                str,
                                StandardTableaux,
                                lambda x: x.size(),
                                lambda x: isinstance(x, StandardTableau)),
    "SkewPartitions":
    SupportedFindStatCollection(lambda x: SkewPartition(literal_eval(x)),
                                str,
                                SkewPartitions,
                                lambda x: x.size(),
                                lambda x: isinstance(x, SkewPartition))}

class FindStatCollections(UniqueRepresentation, Parent):
    r"""
    The class of FindStat collections.

    The elements of this class are combinatorial collections in
    FindStat as of January 2020.  If a new collection was added to the
    web service since then, the dictionary ``SupportedFindStatCollections``
    in this module has to be updated accordingly.

    EXAMPLES::

        sage: from sage.databases.findstat import FindStatCollections
        sage: sorted(c for c in FindStatCollections())                          # optional -- internet
        [Cc0001: Permutations,
         Cc0002: Integer partitions,
         Cc0005: Dyck paths,
         Cc0006: Integer compositions,
         Cc0007: Standard tableaux,
         Cc0009: Set partitions,
         Cc0010: Binary trees,
         Cc0012: Perfect matchings,
         Cc0013: Cores,
         Cc0014: Posets,
         Cc0017: Alternating sign matrices,
         Cc0018: Gelfand-Tsetlin patterns,
         Cc0019: Semistandard tableaux,
         Cc0020: Graphs,
         Cc0021: Ordered trees,
         Cc0022: Finite Cartan types,
         Cc0023: Parking functions,
         Cc0024: Binary words,
         Cc0025: Plane partitions,
         Cc0026: Decorated permutations,
         Cc0027: Signed permutations,
         Cc0028: Skew partitions]
    """
    def _raise_unsupported_error(*args):
        """A placeholder function for unsupported collections that raises an
        error.

        sage: from sage.databases.findstat import FindStatCollection
        sage: FindStatCollection(24).first_terms(lambda x: 1);                  # optional -- internet, random, indirect doctest
        Traceback (most recent call last):
        ...
        NotImplementedError: This FindStatCollection is not yet supported.
        """
        raise NotImplementedError("This FindStatCollection is not yet supported.")

    # The following is used for unknown collections, with the
    # intention to make as much as possible working.
    #
    # In particular, we keep the elements provided by findstat.org as
    # strings.  Therefore, to_string and from_string should do
    # nothing.  Furthermore, in_range always returns True to make
    # submission of sequences possible.
    _findstat_unsupported_collection_default = [None, None, None,
                                                _raise_unsupported_error,
                                                _raise_unsupported_error,
                                                None,
                                                _raise_unsupported_error,
                                                lambda x, l: True,
                                                lambda x: x, lambda x: x]

    def __init__(self):
        """
        Fetch the collections from FindStat.

        TESTS::

            sage: from sage.databases.findstat import FindStatCollections
            sage: C = FindStatCollections()                                     # optional -- internet
            sage: TestSuite(C).run()                                            # optional -- internet
        """
        fields = "LevelsWithSizes,Name,NamePlural,NameWiki"
        url = FINDSTAT_API_COLLECTIONS + "?fields=" + fields
        response = urlopen(url)
        d = json.load(response, object_pairs_hook=OrderedDict)
        self._findstat_collections = d["included"]["Collections"]
        for id, data in self._findstat_collections.items():
            data["LevelsWithSizes"] = OrderedDict((literal_eval(level), size)
                                                  for level, size in data["LevelsWithSizes"].items())
            if data["NameWiki"] in SupportedFindStatCollections:
                data["Code"] = SupportedFindStatCollections[data["NameWiki"]]
            else:
                print("There is a new collection available at `%s`: %s." % (findstat(), data["NamePlural"]))
                print("To use it with this interface, it has to be added to the dictionary SupportedFindStatCollections in src/sage/databases/findstat.py of the SageMath distribution.  Please open a ticket on trac!")
#                print("Very likely, the following code would work:")
#                fields = "SageCodeElementToString,SageCodeElementsOnLevel,SageCodeStringToElement"
#                url = FINDSTAT_API_COLLECTIONS + id + "?fields=" + fields
#                print(json.load(urlopen(url))["included"]["Collections"][id])

        Parent.__init__(self, category=Sets())

    def _element_constructor_(self, entry):
        """Initialize a FindStat collection.

        INPUT:

        see :class:`FindStatCollection`.

        TESTS:

        Create an object and find its collection::

            sage: from sage.databases.findstat import FindStatCollection, FindStatCollections
            sage: sorted([FindStatCollection(c.first_terms(lambda x: 0)[0][0]) for c in FindStatCollections() if c.is_supported()])   # optional -- internet, random
            [Cc0001: Permutations,
             Cc0002: Integer partitions,
             Cc0005: Dyck paths,
             Cc0006: Integer compositions,
             Cc0007: Standard tableaux,
             Cc0009: Set partitions,
             Cc0010: Binary trees,
             Cc0012: Perfect matchings,
             Cc0013: Cores,
             Cc0014: Posets,
             Cc0017: Alternating sign matrices,
             Cc0018: Gelfand-Tsetlin patterns,
             Cc0019: Semistandard tableaux,
             Cc0020: Graphs,
             Cc0021: Ordered trees,
             Cc0022: Finite Cartan types,
             Cc0023: Parking functions,
             Cc0024: Binary words]

            sage: FindStatCollection(Permutation([1,2,3]))                      # optional -- internet
            Cc0001: Permutations

            sage: FindStatCollection(Permutations(3))                           # optional -- internet
            a subset of Cc0001: Permutations

            sage: FindStatCollection(PerfectMatching([[1,2]]))                  # optional -- internet
            Cc0012: Perfect matchings

            sage: FindStatCollection(PerfectMatchings(4))                       # optional -- internet
            a subset of Cc0012: Perfect matchings

            sage: FindStatCollection(SetPartition([[1,2]]))                     # optional -- internet
            Cc0009: Set partitions

            sage: FindStatCollection(CartanType("A3"))                          # optional -- internet
            Cc0022: Finite Cartan types

            sage: cc = FindStatCollection(graphs(3)); cc                        # optional -- internet
            a subset of Cc0020: Graphs
            sage: cc.first_terms(lambda x: x.edges(False)).list()               # optional -- internet
            [(Graph on 3 vertices, []),
             (Graph on 3 vertices, [(0, 2)]),
             (Graph on 3 vertices, [(0, 2), (1, 2)]),
             (Graph on 3 vertices, [(0, 1), (0, 2), (1, 2)])]

            sage: len(cc.first_terms(lambda x: x.edges(False)).list())          # optional -- internet
            4

        """
        if isinstance(entry, self.Element):
            return entry

        if isinstance(entry, string_types):
            # find by name in self._findstat_collections (ignoring case and spaces)
            def normalize(e):
                return "".join(e.split()).upper()

            for id, data in self._findstat_collections.items():
                if normalize(entry) in (normalize(id),
                                        normalize(data["NameWiki"]),
                                        normalize(data["NamePlural"]),
                                        normalize(data["Name"])):
                    return self.element_class(self, id, data, None)

        elif isinstance(entry, (int, Integer)):
            # find by id in _findstat_collections
            for id, data in self._findstat_collections.items():
                if entry == int(id[2:]):
                    return self.element_class(self, id, data, None)

        else:
            # find collection given an object or a constructor

            # it would be good to first check whether entry is an
            # element rather than a type - unfortunately, we cannot
            # test with isinstance(_, SageObject), since this is True
            # for CartanType.

            # first check whether the class fits:
            for id, data in self._findstat_collections.items():
                if ("Code" in data
                    and (data["Code"].is_element(entry)
                         # elements_on_level is rarely equal to entry
                         # (it may be a function), but it is
                         # convenient for some types
                         or data["Code"].elements_on_level == entry)):
                    return self.element_class(self, id, data, None)

            # check whether entry is iterable (it's not a string!)
            # producing a subset of objects

            # WARNING: SkewPartition is iterable and produces two
            # elements of Partition

            # WARNING: we have to remember all elements of the
            # generator, because it is used again, for example in
            # self.first_terms
            try:
                entries = lazy_list(entry)
                obj = entries[0]
            except (TypeError, IndexError):
                pass
            else:
                for id, data in self._findstat_collections.items():
                    if "Code" in data and data["Code"].is_element(obj):
                        return self.element_class(self, id, data, entries)

        raise ValueError("Could not find FindStat collection for %s." %str(entry))

    def _repr_(self):
        """
        Return the representation of the set of FindStat collections.

        EXAMPLES::

            sage: from sage.databases.findstat import FindStatCollections
            sage: FindStatCollections()                                         # optional -- internet
            Set of combinatorial collections used by FindStat
        """
        return "Set of combinatorial collections used by FindStat"

    def __iter__(self):
        """
        Return an iterator over all FindStat collections.

        EXAMPLES::

            sage: from sage.databases.findstat import FindStatCollections
            sage: sorted(FindStatCollections())[0]                              # optional -- internet
            Cc0001: Permutations
        """
        for c in self._findstat_collections:
            yield FindStatCollection(c)

    Element = FindStatCollection
