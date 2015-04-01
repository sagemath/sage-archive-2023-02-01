r"""
FindStat - the Combinatorial Statistic Finder.

You can use the sage interface to FindStat to:

- identify a combinatorial statistic from the values on a few small objects,
- obtain more terms, formulae, references, etc. for a given statistic,
- edit statistics and submit new statistics.

To access the database, use :class:`findstat<FindStat>`::

    sage: findstat
    The Combinatorial Statistic Finder (http://www.findstat.org/)

A guided tour
-------------

Introduction
^^^^^^^^^^^^

Let us fix three notions:

- a combinatorial collection is a set `S` with interesting combinatorial properties,
- a combinatorial map is a combinatorially interesting map `f: S \to S'`,
- a combinatorial statistic is a combinatorially interesting map `s: S \to \ZZ`.

Retrieving information
^^^^^^^^^^^^^^^^^^^^^^

The most straightforward application of the FindStat interface is to
gather information about a combinatorial statistic.  To do this, we
supply :class:`findstat<FindStat>` with a list of `(object, value)`
pairs.  For example::

    sage: PM8 = PerfectMatchings(8)
    sage: r = findstat([(m, m.number_of_nestings()) for m in PM8]); r           # optional -- internet,random
    0: (St000041: The number of nestings of a perfect matching. , [], 105)
    ...

The result of this query is a list (presented as a
:class:`sage.databases.oeis.FancyTuple`) of triples.  The first
element of each triple is a :class:`FindStatStatistic` `s: S \to
\ZZ`, the second element a list of :class:`FindStatMap`'s `f: S \to
S'`, and the third element is an integer::

    sage: (s, list_f, quality) = r[0]                                           # optional -- internet

In the case at hand, the list of maps is empty and the integer
`quality` equals the number of `(object, value)` pairs passed to
FindStat.  This means, that the statistic `s` matches the data
perfectly.  We can now retrieve the description from the database::

    sage: print s.description()                                                 # optional -- internet,random
    The number of nestings of a perfect matching.
    <BLANKLINE>
    <BLANKLINE>
    This is the number of pairs of edges $((a,b), (c,d))$ such that $a\le c\le d\le b$. i.e., the edge $(c,d)$ is nested inside $(a,b)$.

and check the references::

    sage: s.references()                                                        # optional -- internet,random
    0: [MV] combinatorics of orthogonal polynomials (A. de Medicis et X.Viennot, Moments des q-polynomes de Laguerre et la bijection de Foata-Zeilberger, Adv. Appl. Math., 15 (1994), 262-304)
    1: [SS] R. Simion and D. Stanton. Octabasic Laguerre polynomials and permutation statistics. J. Comput. Appl. Math., ...

If you prefer, you can look at this information also in your browser::

    sage: findstat(41).browse()                                                 # optional -- webbrowser

Another interesting possibility is to look for equidistributed
statistics.  Instead of submitting a list of pairs, we pass a pair of
lists::

    sage: r = findstat((PM8, [m.number_of_nestings() for m in PM8])); r         # optional -- internet,random
    0: (St000041: The number of nestings of a perfect matching. , [], 105)
    1: (St000042: The number of crossings of a perfect matching. , [], 105)
    ...


Let us now look at a slightly more complicated example, where the
submitted statistic is the composition of a sequence of combinatorial
maps and a statistic known to FindStat.  We use the occasion to
advertise yet another way to pass values to FindStat::

    sage: r = findstat(lambda pi: pi.saliances()[0], Permutations(4)); r        # optional -- internet,random
    0: (St000051: The size of the left subtree. , [Mp00069: complement, Mp00061: to increasing tree], 24)
    ...
    sage: (s, list_f, quality) = r[0]                                           # optional -- internet

To obtain the value of the statistic sent to FindStat on a given
object, apply the maps in the list in the given order to this object,
and evaluate the statistic on the result.  For example, let us check
that the result given by FindStat agrees with our statistic on the
following permutation::

    sage: pi = Permutation([3,1,4,5,2]); pi.saliances()[0]
    3

We first have to find out, what the maps and the statistic actually do::

    sage: print s.description()                                                 # optional -- internet,random
    The size of the left subtree.

    sage: print s.code()                                                        # optional -- internet,random
    def statistic(T):
        return T[0].node_number()

    sage: print list_f[0].code() + "\r\n" + list_f[1].code()                    # optional -- internet,random
    def complement(elt):
        n = len(elt)
        return elt.__class__(elt.parent(), map(lambda x: n - x + 1, elt) )
    <BLANKLINE>
    def increasing_tree_shape(elt, compare=min):
        return elt.increasing_tree(compare).shape()

So, the following should coincide with what we sent FindStat::

    sage: pi.complement().increasing_tree_shape()[0].node_number()
    3

Editing and submitting statistics
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Of course, often a statistic will not be in the database::

    sage: findstat([(d, randint(1,1000)) for d in DyckWords(4)])                # optional -- internet
    a new statistic on Cc0005: Dyck paths

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

AUTHORS:

- Martin Rubey (2015): initial version.

Classes and methods
-------------------
"""
from sage.misc.cachefunc import cached_method
from sage.misc.misc import verbose
from sage.structure.sage_object import SageObject
from sage.rings.integer import Integer
from sage.databases.oeis import FancyTuple

from string import join
from ast import literal_eval
from collections import OrderedDict
from urllib import urlencode
from urllib2 import Request, urlopen, HTTPError
import re
import webbrowser
import tempfile
import time
import inspect
import json
import cgi


# Combinatoral collections
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
from sage.combinat.posets.poset_examples import posets
from sage.combinat.tableau import SemistandardTableau, SemistandardTableaux, StandardTableau, StandardTableaux
from sage.combinat.set_partition import SetPartition, SetPartitions
from sage.graphs.graph_generators import graphs

######################################################################
# the FindStat URLs
FINDSTAT_URL                        = 'http://www.findstat.org/'
FINDSTAT_URL_RESULT                 = FINDSTAT_URL + "StatisticFinder/Result/"
FINDSTAT_URL_LOGIN                  = FINDSTAT_URL + "StatisticFinder?action=login"
FINDSTAT_URL_NEW                    = FINDSTAT_URL + 'StatisticsDatabase/NewStatistic/'
FINDSTAT_URL_EDIT                   = FINDSTAT_URL + 'StatisticsDatabase/EditStatistic/'
FINDSTAT_URL_BROWSE                 = FINDSTAT_URL + 'StatisticsDatabase/'

FINDSTAT_URL_DOWNLOADS              = 'http://downloads.findstat.org/'
FINDSTAT_URL_DOWNLOADS_STATISTICS   = FINDSTAT_URL_DOWNLOADS + "statistics/%s.json"
FINDSTAT_URL_DOWNLOADS_COLLECTIONS  = FINDSTAT_URL_DOWNLOADS + "collections.json"
FINDSTAT_URL_DOWNLOADS_MAPS         = FINDSTAT_URL_DOWNLOADS + "maps.json"

######################################################################
# the number of values FindStat allows to search for at most
FINDSTAT_MAX_VALUES = 200
# the number of maps that FindStat should compose at most to find a match
FINDSTAT_MAX_DEPTH = 5
# the number of values FindStat allows to submit at most
FINDSTAT_MAX_SUBMISSION_VALUES = 1200

# the fields of the FindStat database we expect
FINDSTAT_STATISTIC_IDENTIFIER      = 'StatisticIdentifier'
FINDSTAT_STATISTIC_COLLECTION      = 'StatisticCollection'
FINDSTAT_STATISTIC_DATA            = 'StatisticData'
FINDSTAT_STATISTIC_DESCRIPTION     = 'StatisticDescription'
FINDSTAT_STATISTIC_REFERENCES      = 'StatisticReferences'
FINDSTAT_STATISTIC_CODE            = 'StatisticCode'
FINDSTAT_STATISTIC_ORIGINAL_AUTHOR = 'StatisticOriginalAuthor' # unused, designates a dictionary with Name, Email, Time
FINDSTAT_STATISTIC_UPDATE_AUTHOR   = 'StatisticUpdateAuthor'   # unused, designates a dictionary with Name, Email, Time

FINDSTAT_POST_AUTHOR               = 'StatisticAuthor' # designates the name of the author
FINDSTAT_POST_EMAIL                = 'StatisticEmail'
FINDSTAT_POST_SAGE_CELL            = 'SageCellField'   # currently only used as post key
FINDSTAT_POST_EDIT                 = 'EDIT'            # only used as post key

FINDSTAT_COLLECTION_IDENTIFIER                = 'CollectionIdentifier'
FINDSTAT_COLLECTION_NAME                      = 'CollectionName'
FINDSTAT_COLLECTION_NAME_PLURAL               = 'CollectionNamePlural'
FINDSTAT_COLLECTION_NAME_WIKI                 = 'CollectionNameWiki'
FINDSTAT_COLLECTION_PARENT_LEVELS_PRECOMPUTED = 'CollectionLevelsPrecomputed'

FINDSTAT_MAP_IDENTIFIER  = 'MapIdentifier' # should be identical to FINDSTAT_MAP_IDENTIFIER
FINDSTAT_MAP_NAME        = 'MapName'
FINDSTAT_MAP_DESCRIPTION = 'MapDescription'
FINDSTAT_MAP_DOMAIN      = 'MapDomain'
FINDSTAT_MAP_CODOMAIN    = 'MapCodomain'
FINDSTAT_MAP_CODE        = 'MapCode'
FINDSTAT_MAP_CODE_NAME   = 'MapSageName'

FINDSTAT_QUERY_MATCHES       = 'QueryMatches'
FINDSTAT_QUERY_MATCHING_DATA = 'QueryMatchingData'
FINDSTAT_QUERY_MAPS          = 'QueryMaps'

# the entries of this list are required as post arguments for submitting or editing a statistic
FINDSTAT_EDIT_FIELDS = set([FINDSTAT_STATISTIC_IDENTIFIER,
                            FINDSTAT_STATISTIC_COLLECTION,
                            FINDSTAT_STATISTIC_DATA,
                            FINDSTAT_STATISTIC_DESCRIPTION,
                            FINDSTAT_STATISTIC_REFERENCES,
                            FINDSTAT_STATISTIC_CODE,
                            FINDSTAT_POST_AUTHOR,
                            FINDSTAT_POST_EMAIL,
                            FINDSTAT_POST_SAGE_CELL,
                            FINDSTAT_POST_EDIT])

# separates name from description
FINDSTAT_SEPARATOR_NAME = "\r\n"
# separates references
FINDSTAT_SEPARATOR_REFERENCES = "\r\n"

######################################################################

# the format string for using POST
# WARNING: we use cgi.escape to avoid injection problems, thus we expect double quotes as field delimiters.
FINDSTAT_POST_HEADER = """
<script src="http://www.google.com/jsapi"></script>
<script>
    google.load("jquery", "1.3.2");
</script>

<script>
    $(document).ready(function() {$("#form").submit(); });
</script>
"""
FINDSTAT_NEWSTATISTIC_FORM_HEADER = '<form id="form" name="NewStatistic" action="%s" enctype="multipart/form-data" method="post" />'
FINDSTAT_NEWSTATISTIC_FORM_FORMAT = '<input type="hidden" name="%s" value="%s" />'
FINDSTAT_NEWSTATISTIC_FORM_FOOTER = '</form>'

######################################################################
class FindStat():
    r"""
    The Combinatorial Statistic Finder.

    :class:`FindStat` is a class representing results of queries to
    the FindStat database.  This class is also the entry point to
    edit statistics and new submissions.  Use the shorthand
    :class:`findstat<FindStat>` to call it.

    INPUT:

    - an integer or a string representing a valid FindStat identifier
      (e.g. 45 or 'St000045').  The optional argument ``collection``
      should be ``None``, the optional arguments ``depth`` and
      ``max_values`` are ignored.

    - a list of pairs of the form (object, value), or a dictionary
      from sage objects to integer values.  The optional argument
      ``collection`` should be ``None``, the optional arguments
      ``depth`` and ``max_values`` are passed to the finder.

    - a list of pairs of the form (list of objects, list of values),
      or a single pair of the form (list of objects, list of values).
      In each pair there should be as many objects as values.  The
      optional argument ``collection`` should be ``None``, the
      optional arguments ``depth`` and ``max_values`` are passed to
      the finder.

    - a callable and a collection.  The callable is used to generate
      ``max_values`` (object, value) pairs.  The number of terms
      generated may also be controlled by passing an iterable
      collection, such as Permutations(3).  The optional arguments
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

        sage: findstat('St000045')                                              # optional -- internet,random
        St000045: The number of linear extensions of the tree.

        sage: findstat(3)                                                       # optional -- internet,random
        St000003: The number of [[/StandardTableaux|standard Young tableaux]] of the partition.

    The database can be searched by providing a list of pairs::

        sage: q = findstat([(pi, pi.length()) for pi in Permutations(4)]); q    # optional -- internet,random
        0: (St000018: The [[/Permutations/Inversions|number of inversions]] of a permutation., [], 24)
        1: (St000004: The [[/Permutations/Descents-Major|major index]] of a permutation., [Mp00062: inversion-number to major-index bijection], 24)
        ...

    or a dictionary::

        sage: p = findstat({pi: pi.length() for pi in Permutations(4)}); p      # optional -- internet,random
        0: (St000018: The [[/Permutations/Inversions|number of inversions]] of a permutation., [], 24)
        1: (St000004: The [[/Permutations/Descents-Major|major index]] of a permutation., [Mp00062: inversion-number to major-index bijection], 24)
        ...

    Note however, that the results of these two queries are not
    necessarily the same, because we compare queries by the data
    sent, and the ordering of the data might be different::

        sage: p == q                                                            # optional -- internet
        False

    Another possibility is to send a function and a collection.  In
    this case, the function is applied to the first few objects of
    the collection::

        sage: findstat(lambda pi: pi.length(), "Permutations")                  # optional -- internet,random
        0: (St000018: The [[/Permutations/Inversions|number of inversions]] of a permutation., [], 200)
        ...

    To search for a distribution, send a list of lists, or a single pair::

        sage: S4 = Permutations(4); findstat((S4, [pi.length() for pi in S4]))  # optional -- internet,random
        0: (St000004: The [[/Permutations/Descents-Major|major index]] of a permutation., [], 24)
        1: (St000018: The [[/Permutations/Inversions|number of inversions]] of a permutation., [], 24)
        ...
    """
    def __init__(self):
        r"""
        Initialize a cache from integers to
        :class:`FindStatStatistic` instances that avoids retrieving
        the same statistic over and over again.
        """

        self._statistic_cache = dict()

        # user credentials if provided
        self._user_name  = ""
        self._user_email = ""

    def __call__(self, query, collection=None, depth=2, max_values=FINDSTAT_MAX_VALUES):
        r"""
        Return an instance of a :class:`FindStatStatistic`.

        This should be the only way to access
        :class:`FindStatStatistic`.  We do the preprocessing of the
        data here, and call the appropriate method of
        :class:`FindStatStatistic` to launch the query.
        """
        try:
            depth = int(depth)
            assert 0 <= depth <= FINDSTAT_MAX_DEPTH
        except:
            raise ValueError("The depth must be a non-negative integer less than or equal to %i." %FINDSTAT_MAX_DEPTH)

        if collection is None:
            if isinstance(query, str):
                if re.match('^St[0-9]{6}$', query):
                    return self._statistic(Integer(query[2:].lstrip("0")))
                else:
                    raise ValueError("The value %s is not a valid statistic identifier." %query)

            elif isinstance(query, (int, Integer)):
                return self._statistic(query)

            elif isinstance(query, dict):
                # we expect a dictionary from objects to integers
                data = [([key], [value]) for (key, value) in query.iteritems()]
                collection = FindStatCollection(data[0][0][0])
                return FindStatStatistic(id=0, data=data,
                                         first_terms=query,
                                         collection=collection,
                                         depth=depth)._find_by_values(max_values=max_values)

            elif isinstance(query, (list, tuple)):
                # either a pair (list of objects, list of integers)
                # or a list of such or (object, integer) pairs

                # values must always be lists because otherwise we
                # get a trailing comma when printing
                if (len(query) == 2 and
                    isinstance(query[1], (list, tuple)) and
                    len(query[1]) != 0 and
                    isinstance(query[1][0], (int, Integer))):
                    # just a pair
                    if len(query[0]) != len(query[1]):
                           raise ValueError, "FindStat expects the same number of objects as values!"

                    data = [(query[0], list(query[1]))]
                    collection = FindStatCollection(data[0][0][0])
                    return FindStatStatistic(id=0, data=data,
                                             collection=collection,
                                             depth=depth)._find_by_values(max_values=max_values)
                else:
                    is_statistic = True
                    data = []
                    for (key, value) in query:
                        if isinstance(value, (list, tuple)):
                            data += [(key, list(value))]
                            if len(key) != len(value):
                                raise ValueError, "FindStat expects the same number of objects as values!"
                            is_statistic = False
                        else:
                            data += [([key], [value])]
                    collection = FindStatCollection(data[0][0][0])
                    if is_statistic:
                        return FindStatStatistic(id=0, data=data,
                                                 collection=collection,
                                                 first_terms=query,
                                                 depth=depth)._find_by_values(max_values=max_values)
                    else:
                        return FindStatStatistic(id=0, data=data,
                                                 collection=collection,
                                                 depth=depth)._find_by_values(max_values=max_values)

            else:
                raise ValueError("The given query, %s, cannot be used for a FindStat search." %query)

        else:
            if callable(query):
                if not isinstance(collection, FindStatCollection):
                    collection = FindStatCollection(collection)
                first_terms = collection.first_terms(query, max_values=max_values)
                data = [([key], [value]) for (key, value) in first_terms]
                try:
                    code = inspect.getsource(query)
                except IOError:
                    _ = verbose("inspect.getsource could not get code from function provided", caller_name='FindStat')
                    code = ""
                return FindStatStatistic(id=0, first_terms=first_terms,
                                         data=data, function=query, code=code,
                                         collection=collection,
                                         depth=depth)._find_by_values(max_values=max_values)
            else:
                raise ValueError("The given query, %s, cannot be used for a FindStat search." %query)

    def __repr__(self):
        r"""
        Return the representation of ``self``.

        EXAMPLES::

            sage: findstat
            The Combinatorial Statistic Finder (http://www.findstat.org/)
        """
        return "The Combinatorial Statistic Finder (%s)" % FINDSTAT_URL

    def browse(self):
        r"""
        Open the FindStat web page in a browser.

        EXAMPLES::

            sage: findstat.browse()                                             # optional -- webbrowser
        """
        webbrowser.open(FINDSTAT_URL)

    def set_user(self, name=None, email=None):
        r"""
        Set the user for this session.

        INPUT:

        - ``name`` -- the name of the user

        - ``email`` -- an email address of the user

        .. NOTE::

            It is usually more convenient to login into the FindStat
            web page using the :meth:`login` method.
        """
        if not isinstance(name, str):
            raise ValueError("The given name is not a string")
        if not isinstance(email, str):
            raise ValueError("The given email address is not a string")
        self._user_name  = name
        self._user_email = email

    def login(self):
        r"""
        Open the FindStat login page in a browser.
        """
        webbrowser.open(FINDSTAT_URL_LOGIN)

    ######################################################################

    def _statistic(self, id):
        r"""
        INPUT:

        - ``id``, an integer designating the FindStat id of a statistic

        OUTPUT:

        - A :class:`FindStatStatistic` instance.

        .. TODO::

            this method caches the statistics.  It may make sense to
            provide a method that clears the cache, or reloads a
            single statistic.
        """
        if id > 0:
            if id not in self._statistic_cache.keys():
                self._statistic_cache[id] = FindStatStatistic(id)._find_by_id()
            return self._statistic_cache[id]
        else:
            raise ValueError("The statistic identifier must be at least 1.")

######################################################################

class FindStatStatistic(SageObject):
    r"""
    The class of FindStat statistics.

    Do not instantiate this class directly.  Instead, use
    :class:`findstat<FindStat>`.
    """
    def __init__(self, id, first_terms=None, data=None, function=None, code="", collection=None, depth=None):
        self._depth = depth
        self._query = None
        self._modified = False

        self._id = id
        self._result = None

        self._first_terms = first_terms
        self._data = data
        self._function = function
        self._code = code
        self._collection = collection

        self._description = ""
        self._references = ""

    def __repr__(self):
        if self._query == "ID":
            if self._modified:
                return "%s(modified): %s" % (self.id_str(), self.name())
            else:
                return "%s: %s" % (self.id_str(), self.name())

        elif self._query == "data":
            if len(self._result) == 0:
                return "a new statistic on " + self._collection.__repr__()
            else:
                return self._result.__repr__()

        else:
            raise ValueError("self._query should be either 'ID' or 'data', but is %s" %self._query)

    def __eq__(self, other):
        """
        .. TODO::

            this is *very* rudimentary
        """
        if self._query == "ID" and other._query == "ID":
            if self._modified or other._modified:
                return False
            else:
                return self._id == other._id
        elif self._query == "data" and other._query == "data":
            if self._modified or other._modified:
                return False
            else:
                return self._data == other._data
        else:
            return False

    def __ne__(self, other):
        return not self.__eq__(other)

    ######################################################################
    # query = "ID"
    ######################################################################

    def _find_by_id(self):
        r"""
        Expects that `_id` is a valid identifier.

        OUTPUT:

        - self

        TESTS::

            sage: findstat(999999)
            Traceback (most recent call last):
            ...
            ValueError: St999999 is not a FindStat statistic identifier.
        """
        self._query = "ID"

        # get the database entry from FindStat
        url = FINDSTAT_URL_DOWNLOADS_STATISTICS %self.id_str()
        _ = verbose("Fetching URL %s ..." %url, caller_name='FindStat')
        try:
            self._raw = json.load(urlopen(url), object_pairs_hook=OrderedDict)
        except HTTPError, error:
            if error.code == 404:
                raise ValueError("%s is not a FindStat statistic identifier." %self.id_str())
            else:
                raise

        self._description = self._raw[FINDSTAT_STATISTIC_DESCRIPTION].encode("utf-8")
        self._references = self._raw[FINDSTAT_STATISTIC_REFERENCES].encode("utf-8")
        self._collection = FindStatCollection(self._raw[FINDSTAT_STATISTIC_COLLECTION])
        self._code = self._raw[FINDSTAT_STATISTIC_CODE]

        from_str = self._collection.from_string()
        # we want to keep FindStat's ordering here!
        if from_str is None:
            self._first_terms = self._raw[FINDSTAT_STATISTIC_DATA]
        else:
            self._first_terms = [(from_str(obj), Integer(val)) for (obj, val) in self._raw[FINDSTAT_STATISTIC_DATA].iteritems()]
        return self

    ######################################################################
    # query = "data"
    ######################################################################

    def _find_by_values(self, max_values=FINDSTAT_MAX_VALUES):
        r"""
        Expects that data is a list of pairs of the form (list of
        objects, list of values), each containing as many values as
        objects, and that collection is appropriately set.

        OUTPUT:

        - self

        TESTS::

            sage: findstat(lambda x: 1, "Permutations", depth=100)
            Traceback (most recent call last):
            ...
            ValueError: The depth must be a non-negative integer less than or equal to 5.
        """
        self._query = "data"

        # FindStat allows to search for at most FINDSTAT_MAX_VALUES
        # values.  For the user's convenience, from data, we take the
        # first min(max_values, FINDSTAT_MAX_VALUES) such that all
        # elements are in the precomputed range
        data = []
        total = min(max_values, FINDSTAT_MAX_VALUES)
        for (elements, values) in self._data:
            if total >= len(elements):
                if all(self._collection.in_range(e) for e in elements):
                    data += [(elements, values)]
                    total -= len(elements)

        # this might go wrong:
        try:
            assert data != []
        except:
            raise ValueError("after discarding elements not in the range, and keeping less than %s values, nothing remained to send to FindStat." %FINDSTAT_MAX_VALUES)

        url = FINDSTAT_URL_RESULT + self._collection._url_name + "/"

        to_str = self._collection.to_string()
        stat = [(map(to_str, keys), str(values)[1:-1]) for (keys, values) in data]

        stat_str = join([join(keys, "\n") + "\n====> " + values for (keys, values) in stat], "\n")
        _ = verbose("Sending data to FindStat %s" %stat_str, caller_name='FindStat')

        values = urlencode({"freedata": stat_str, "depth": str(self._depth), "caller": "Sage"})
        _ = verbose("Fetching URL %s with encoded data %s" %(url, values), caller_name='FindStat')

        request = Request(url, data=values)
        _ = verbose("Requesting %s" %request, caller_name='FindStat')

        response = urlopen(request)
        _ = verbose("Response was %s" %response.info(), caller_name='FindStat')

        try:
            result = json.load(response)
            self._result = FancyTuple((findstat(match[FINDSTAT_STATISTIC_IDENTIFIER]),
                                       [FindStatMap(mp[FINDSTAT_MAP_IDENTIFIER]) for mp in match[FINDSTAT_QUERY_MAPS]],
                                       len(match[FINDSTAT_QUERY_MATCHING_DATA]))
                                      for match in result[FINDSTAT_QUERY_MATCHES])
            return self
        except:
            raise IOError, "FindStat did not answer with a json response."

    ######################################################################

    def __getitem__(self, key):
        if self._query == "ID":
            raise TypeError("Use 'first_terms' to access the values of the statistic.")

        elif self._query == "data":
            return self._result[key]

        else:
            raise ValueError("self._query should be either 'ID' or 'data', but is %s" %self._query)

    def id(self):
        r"""
        Return the FindStat identifier of the statistic.

        OUTPUT:

        The FindStat identifier of the statistic (or 0), as an integer.

        EXAMPLES::

            sage: findstat(1).id()                                              # optional -- internet
            1
        """
        return self._id

    def id_str(self):
        r"""
        Return the FindStat identifier of the statistic.

        OUTPUT:

        The FindStat identifier of the statistic (or 'St000000'), as a string.

        EXAMPLES::

            sage: findstat(1).id_str()                                          # optional -- internet
            'St000001'
        """
        id = str(self._id)
        return 'St000000'[:-len(id)] + id

    def data(self):
        r"""
        Return the data used for querying the FindStat database.

        OUTPUT:

        The data provided by the user to query the FindStat database.
        When the database was searched using an identifier, ``data``
        is ``None``.

        EXAMPLES::

            sage: S4 = Permutations(4); findstat((S4, [pi.length() for pi in S4])).data()           # optional -- internet
            [(Standard permutations of 4,
              [0, 1, 1, 2, 2, 3, 1, 2, 2, 3, 3, 4, 2, 3, 3, 4, 4, 5, 3, 4, 4, 5, 5, 6])]
        """
        return self._data

    def modified(self):
        r"""
        Return whether the statistic was modified.

        OUTPUT:

        True, if the statistic was modified.
        """
        return self._modified

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
        return self._collection

    def function(self):
        r"""
        Return the function used to compute the values of the statistic.

        OUTPUT:

        The function used to compute the values of the statistic, or
        ``None``.

        EXAMPLES::

            sage: findstat(lambda pi: pi.length(), "Permutations").function()   # optional -- internet
            ...
            <function <lambda> at ...>
        """
        return self._function

    def first_terms(self):
        r"""
        Return the first terms of the statistic.

        OUTPUT:

        A list of pairs of the form ``(object, value)`` where object
        is an object is a sage object representing an element of the
        appropriate collection and value is an integer.  The list
        contains exactly the pairs in the database.

        EXAMPLES::

            sage: findstat(1).first_terms()                                     # optional -- internet,random
            [([1], 1),
             ([1, 2], 1),
             ([2, 1], 1),
             ([1, 2, 3], 1),
             ([1, 3, 2], 1),
             ([2, 1, 3], 1),
             ...
        """
        return self._first_terms

    def first_terms_str(self):
        r"""
        Return the first terms of the statistic in the format needed
        for a FindStat query.

        OUTPUT:

        A string, where each line is of the form `object => value`,
        where `object` is the string representation of an element of
        the appropriate collection as used by FindStat and value is
        an integer.

        EXAMPLES::

            sage: findstat(1).first_terms_str()[:10]                            # optional -- internet,random
            '[1] => 1\r\n'
        """
        if self._first_terms != None:
            to_str = self._collection.to_string()
            return join([to_str(key) + " => " + str(val)
                         for (key, val) in self._first_terms], "\r\n")
        else:
            return ""

    def description(self):
        r"""
        Return the description of the statistic.

        OUTPUT:

        A string, whose first line is used as the name of the
        statistic.

        EXAMPLES::

            sage: print findstat(1).description()                               # optional -- internet,random
            The number of ways to write a permutation as a minimal length product of simple transpositions.
            <BLANKLINE>
            That is, the number of reduced words for the permutation.  E.g., there are two reduced words for $[3,2,1] = (1,2)(2,3)(1,2) = (2,3)(1,2)(2,3)$.
        """
        return self._description

    def set_description(self, value):
        r"""
        Set the description of the statistic.

        INPUT:

        A string, whose first line is used as the name of the
        statistic.
        """
        if value != self._description:
            self._modified = True
            self._description = value

    def name(self):
        r"""
        Return the name of the statistic.

        OUTPUT:

        A string, which is just the first line of the description of
        the statistic.

        EXAMPLES::

            sage: findstat(1).name()                                            # optional -- internet,random
            u'The number of ways to write a permutation as a minimal length product of simple transpositions.'
        """
        return self._description.partition(FINDSTAT_SEPARATOR_NAME)[0]

    def references(self):
        r"""
        Return the references associated with the statistic.

        OUTPUT:

        An instance of :class:`sage.databases.oeis.FancyTuple`, each
        item corresponds to a reference.

        .. TODO::

            Since the references in the database are sometimes not
            formatted properly, this method is unreliable.  The
            string representation can be obtained via
            :attr:`_references`.

        EXAMPLES::

            sage: findstat(1).references()                                      # optional -- internet,random
            0: P. Edelman and C. Greene, Balanced tableaux, Adv. in Math., 63 (1987), pp. 42-99.
            1: [[OEIS:A005118]]
            2: [[oeis:A246865]]
        """
        l = [ref.strip() for ref in self._references.split(FINDSTAT_SEPARATOR_REFERENCES)]
        return FancyTuple([ref for ref in l if ref != ""])

    def set_references(self, value):
        r"""
        Set the references associated with the statistic.

        INPUT:

        A string.  The individual references should be separated by
        FINDSTAT_SEPARATOR_REFERENCES, which is "\\r\\n".
        """
        if value != self._references:
            self._modified = True
            self._references = value

    def code(self):
        r"""
        Return the code associated with the statistic.

        OUTPUT:

        A string.  Contributors are encouraged to submit sage code in the form::

            def statistic(x):
                ...

        but the string may also contain code for other computer
        algebra systems.

        EXAMPLES::

            sage: print findstat(1).code()                                      # optional -- internet,random
            def statistic(x):
                return len(x.reduced_words())

            sage: print findstat(118).code()                                    # optional -- internet,random
            (* in Mathematica *)
            tree = {{{{}, {}}, {{}, {}}}, {{{}, {}}, {{}, {}}}};
            Count[tree, {{___}, {{___}, {{___}, {___}}}}, {0, Infinity}]
        """
        return self._code

    def set_code(self, value):
        r"""
        Set the code associated with the statistic.

        INPUT:

        A string.  Contributors are encouraged to submit sage code in
        the form::

            def statistic(x):
                ...
        """
        if value != self._code:
            self._modified = True
            self._code = value

    ######################################################################
    # browse current statistic
    ######################################################################

    def browse(self):
        r"""
        Open the FindStat web page of the statistic in a browser.

        EXAMPLES::

            sage: findstat(45).browse()                                         # optional -- webbrowser
        """
        if self._query == "ID":
            webbrowser.open(FINDSTAT_URL_BROWSE + self.id_str())
        else:
            raise NotImplementedError("Would be nice to show the result of the query in the webbrowser.")

    ######################################################################
    # submit current (possibly incompletely defined) statistic
    ######################################################################

    def submit(self, max_values=FINDSTAT_MAX_SUBMISSION_VALUES):
        r"""
        Open the FindStat web page for editing the statistic ``self`` in a browser.

        .. TODO::

            decide whether we want to somehow take into account when
            there is a statistic that matches the data with depth 0.
        """
        # if the statistic is given as a function, and we have a new
        # statistic then update first_terms

        # it is not clear whether we want to do this also for old statistics.
        if self.function() and self.id() == 0:
            self._first_terms = self.collection().first_terms(self.function(),
                                                              max_values=max_values)

        args = dict()
        args[FINDSTAT_STATISTIC_IDENTIFIER]  = self._id
        args[FINDSTAT_STATISTIC_COLLECTION]  = str(self.collection().id())
        args[FINDSTAT_STATISTIC_DATA]        = self.first_terms_str()
        args[FINDSTAT_STATISTIC_DESCRIPTION] = self._description
        args[FINDSTAT_STATISTIC_REFERENCES]  = self._references
        args[FINDSTAT_STATISTIC_CODE]        = self.code()
        args[FINDSTAT_POST_SAGE_CELL]        = ""
        args[FINDSTAT_POST_EDIT]             = ""
        args[FINDSTAT_POST_AUTHOR]           = findstat._user_name
        args[FINDSTAT_POST_EMAIL]            = findstat._user_email

        assert set(args.keys()) == FINDSTAT_EDIT_FIELDS, "It appears that the list of required post variables for editing a statistic has changed.  Please update FindStatStatistic.submit()."

        # write the file
        f = tempfile.NamedTemporaryFile(delete=False)
        _ = verbose("Created temporary file %s" %f.name, caller_name='FindStat')
        f.write(FINDSTAT_POST_HEADER)
        if self.id() == 0:
            f.write(FINDSTAT_NEWSTATISTIC_FORM_HEADER %FINDSTAT_URL_NEW)
        else:
            f.write(FINDSTAT_NEWSTATISTIC_FORM_HEADER %(FINDSTAT_URL_EDIT+self.id_str()))
        for key, value in args.iteritems():
            _ = verbose("writing argument %s" %key, caller_name='FindStat')
            value_encoded = cgi.escape(str(value), quote=True)
            _ = verbose("%s" %value_encoded, caller_name='FindStat')
            f.write((FINDSTAT_NEWSTATISTIC_FORM_FORMAT %(key, value_encoded)))
        f.write(FINDSTAT_NEWSTATISTIC_FORM_FOOTER)
        f.close()
        _ = verbose("Opening file with webbrowser", caller_name='FindStat')
        _ = webbrowser.open(f.name)

        _ = verbose("Waiting a little before deleting the temporary file", caller_name='FindStat')
        time.sleep(1)

        f.unlink(f.name)

    # editing and submitting is really the same thing
    edit = submit

class FindStatCollection(SageObject):
    r"""
    The class of FindStat collections.
    """
    # helper for generation of CartanTypes
    def _finite_irreducible_cartan_types_by_rank(n):
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

    # we set up a dict of FindStat collections containing, with key
    # being the FINDSTAT_COLLECTION_IDENTIFIER to tuples, containing in this order:

    # * the FindStat name                                        (FINDSTAT_COLLECTION_NAME)
    # * the FindStat name plural                                 (FINDSTAT_COLLECTION_NAME_PLURAL)
    # * url's as needed by FindStat                              (FINDSTAT_COLLECTION_NAME_WIKI)
    # * sage element constructor
    # * sage constructor                                         (would be parent_initializer)
    # * list of arguments for constructor                        (FINDSTAT_COLLECTION_PARENT_LEVELS_PRECOMPUTED)
    # * a method to get the size of the sage object
    # * the (FindStat) string representations of the sage object (would be element_repr)
    # * sage constructors of the FindStat string representation  (would be element_constructor)

    # several fields are initialised with 'None', they are updated upon the first call to this class
    _findstat_collections = {
        17: [None, None, None, AlternatingSignMatrix, AlternatingSignMatrices, None,
             lambda x: x.to_matrix().nrows(),
             lambda x: str(map(list, list(x._matrix))),
             lambda x: AlternatingSignMatrix(literal_eval(x))],
        10: [None, None, None, BinaryTree,            BinaryTrees,             None,
             lambda x: x.node_number(),
             str,
             lambda x: BinaryTree(str(x))],
        13: [None, None, None, Core,                  lambda x: Cores(x[1], x[0]),
             None,
             lambda x: (x.length(), x.k()),
             lambda X: "( "+X._repr_()+", "+str(X.k())+" )",
             lambda x: (lambda pi, k: Core(pi, k))(*literal_eval(x))],
        5:  [None, None, None, DyckWord,              DyckWords,               None,
             lambda x: x.length()/2,
             lambda x: str(list(DyckWord(x))),
             lambda x: DyckWord(literal_eval(x))],
        22: [None, None, None, CartanType_abstract,   _finite_irreducible_cartan_types_by_rank,
             None,
             lambda x: x.rank(),
             str,
             lambda x: CartanType(*literal_eval(str(x)))],
        18: [None, None, None, GelfandTsetlinPattern, lambda x: GelfandTsetlinPatterns(*x),
             None,
             lambda x: (len(x), max(max(row) for row in x)),
             str,
             lambda x: GelfandTsetlinPattern(literal_eval(x))],
        20: [None, None, None, Graph,                 graphs,
             None,
             lambda x: x.num_verts(),
             lambda X: str((sorted(X.canonical_label().edges(False)), X.num_verts())),
             lambda x: (lambda E, V: Graph([range(V), lambda i,j: (i,j) in E or (j,i) in E], immutable=True))(*literal_eval(x))],
        6:  [None, None, None, Composition,           Compositions,            None,
             lambda x: x.size(),
             str,
             lambda x: Composition(literal_eval(x))],
        2:  [None, None, None, Partition,             Partitions,              None,
             lambda x: x.size(),
             str,
             lambda x: Partition(literal_eval(x))],
        21: [None, None, None, OrderedTree,           OrderedTrees,            None,
             lambda x: x.node_number(),
             str,
             lambda x: OrderedTree(literal_eval(x))],
        23: [None, None, None, ParkingFunction_class, ParkingFunctions,        None,
             len,
             str,
             lambda x: ParkingFunction(literal_eval(x))],
        12: [None, None, None, PerfectMatching,       PerfectMatchings,        None,
             lambda x: x.size(),
             str,
             lambda x: PerfectMatching(literal_eval(x))],
        1:  [None, None, None, Permutation,           Permutations,            None,
             lambda x: x.size(),
             str,
             lambda x: Permutation(literal_eval(x))],
        14: [None, None, None, FinitePoset,           posets,                  None,
             lambda x: x.cardinality(),
             lambda X: str((sorted(X._hasse_diagram.canonical_label().cover_relations()), len(X._hasse_diagram.vertices()))),
             lambda x: (lambda R, E: Poset((range(E), R)))(*literal_eval(x))],
        19: [None, None, None, SemistandardTableau,   lambda x: SemistandardTableaux(x, max_entry=4),
             None,
             lambda x: x.size(),
             str,
             lambda x: SemistandardTableau(literal_eval(x))],
        9:  [None, None, None, SetPartition,          SetPartitions,           None,
             lambda x: x.size(),
             str,
             lambda x: SetPartition(literal_eval(x.replace('{','[').replace('}',']')))],
        7:  [None, None, None, StandardTableau,       StandardTableaux,        None,
             lambda x: x.size(),
             str,
             lambda x: StandardTableau(literal_eval(x))]}

    r"""
    Objects are normalized using the method :meth:`to_string`.  This
    method should apply to objects produced by :meth:`first_terms` as
    well as to objects produced by :meth:`from_string`.

    EXAMPLES::

        sage: from sage.databases.findstat import FindStatCollection
        sage: FindStatCollection("Permutations")                                # optional -- internet
        Cc0001: Permutations

    TESTS::

        # when one of these tests does not pass, there is probably a new collection to be added.
        sage: cdata = FindStatCollection._findstat_collections.values()         # optional -- internet
        sage: cl = [FindStatCollection(x[0]) for x in cdata]
        # create an object and find its collection
        sage: [FindStatCollection(c.first_terms(lambda x: 0, max_values=1)[0][0]) for c in cl]
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
         Cc0023: Parking functions]
    """


    def __init__(self, entry):
        r"""
        Initialize a FindStat collection.

        INPUT:

        - a string, eg. 'Dyck paths' or "DyckPaths", case-insensitve
        - an integer designating the FindStat id of the collection
        - a sage object belonging to a collection
        - an iterable producing a sage object belonging to a collection

        EXAMPLES::

            sage: from sage.databases.findstat import FindStatCollection
            sage: FindStatCollection("Dyck paths")                              # optional -- internet
            Cc0005: Dyck paths

            sage: FindStatCollection(5)                                         # optional -- internet
            Cc0005: Dyck paths

            sage: FindStatCollection(DyckWord([1,0,1,0]))                       # optional -- internet
            Cc0005: Dyck paths

            sage: FindStatCollection(DyckWords(2))                              # optional -- internet
            Cc0005: Dyck paths
        """
        if not FindStatCollection._findstat_collections.values()[0][0]:
            for j in json.load(urlopen(FINDSTAT_URL_DOWNLOADS_COLLECTIONS)):
                c = FindStatCollection._findstat_collections[j[FINDSTAT_COLLECTION_IDENTIFIER]]
                c[0] = j[FINDSTAT_COLLECTION_NAME]
                c[1] = j[FINDSTAT_COLLECTION_NAME_PLURAL]
                c[2] = j[FINDSTAT_COLLECTION_NAME_WIKI]
                c[5] = literal_eval(j[FINDSTAT_COLLECTION_PARENT_LEVELS_PRECOMPUTED])

        self._sageconstructor_overridden = None
        bad = True
        def initialize_with(id, c):
            self._id = id
            (self._name, self._name_plural, self._url_name,
             self._sageclass, self._sageconstructor, self._range,
             self._to_size, self._to_str, self._from_str) = c

        if isinstance(entry, (str, unicode)):
            # find by name in _findstat_collections
            for (id, c) in FindStatCollection._findstat_collections.iteritems():
                if entry.upper() in (c[0].upper(), c[1].upper(), c[2].upper()):
                    initialize_with(id, c)
                    bad = False
                    break

        elif isinstance(entry, (int, Integer)):
            # find by id in _findstat_collections
            for (id, c) in FindStatCollection._findstat_collections.iteritems():
                if entry == id:
                    initialize_with(id, c)
                    bad = False
                    break
        else:
            # find collection given an object or a constructor

            # unfortunately, we cannot test with
            # isinstance(_, SageObject), since this is True for
            # CartanType.

            # TODO: entry == c[4] will work rarely because c[4] might be a function!
            # also, the error handling is only necessary because of this...
            for (id, c) in FindStatCollection._findstat_collections.iteritems():
                try:
                    if isinstance(entry, c[3]) or entry == c[4]:
                        initialize_with(id, c)
                        bad = False
                        break
                except:
                    # examples are
                    # graphs:
                    # TypeError: cannot compare graph to non-graph (<class 'sage.combinat.permutation.Permutations'>)
                    # perfect matchings:
                    # TypeError: descriptor 'parent' of 'sage.structure.sage_object.SageObject' object needs an argument
                    pass

            if bad:
                # check whether entry is iterable (it's not a string!)
                try:
                    obj = iter(entry).next()
                    self._sageconstructor_overridden = entry
                    for (id, c) in FindStatCollection._findstat_collections.iteritems():
                        if isinstance(obj, c[3]):
                            initialize_with(id, c)
                            bad = False
                            break

                except TypeError:
                    pass
        if bad:
            raise ValueError, "Could not find FindStat collection for " + str(entry)

    def __eq__(self, other):
        """
        TESTS::

            sage: from sage.databases.findstat import FindStatCollection
            sage: FindStatCollection("Permutations") == FindStatCollection("Permutations")          # optional -- internet
            True

            sage: FindStatCollection("Permutations") == FindStatCollection("Integer Partitions")    # optional -- internet
            False
        """
        return self.id() == other.id()

    def __ne__(self, other):
        """
        TESTS::

            sage: from sage.databases.findstat import FindStatCollection
            sage: FindStatCollection("Permutations") != FindStatCollection("Permutations")          # optional -- internet
            False

            sage: FindStatCollection("Permutations") != FindStatCollection("Integer Partitions")    # optional -- internet
            True
        """
        return not self.__eq__(other)

    def in_range(self, element):
        r"""
        Check whether an element of the collection is in FindStat's precomputed range.

        INPUT:

        - ``element`` -- a sage object that belongs to the collection

        OUTPUT:

        - ``True``, if ``element`` is used by the FindStat search
          engine, and ``False`` if it is ignored.

        EXAMPLES::

            sage: from sage.databases.findstat import FindStatCollection
            sage: c = FindStatCollection("GelfandTsetlinPatterns")              # optional -- internet
            sage: c.in_range(GelfandTsetlinPattern([[3, 1], [1]]))              # optional -- internet
            True
            sage: c.in_range(GelfandTsetlinPattern([[4, 1], [1]]))              # optional -- internet,random
            False
        """
        n = self.to_size()(element)
        return (n in self._range and element in self._sageconstructor(n))

    def first_terms(self, statistic, max_values=FINDSTAT_MAX_SUBMISSION_VALUES):
        r"""
        Compute the first few terms of the given statistic.

        INPUT:

        - ``statistic`` -- a callable.

        - ``max_values`` -- the number of terms to compute at most.

        OUTPUT:

        - a list of pairs of the form (object, value).

        EXAMPLES::

            sage: from sage.databases.findstat import FindStatCollection
            sage: c = FindStatCollection("GelfandTsetlinPatterns")              # optional -- internet
            sage: c.first_terms(lambda x: 1, max_values=10)                     # optional -- internet,random
            [([[0]], 1),
             ([[1]], 1),
             ([[2]], 1),
             ([[3]], 1),
             ([[0, 0], [0]], 1),
             ([[1, 0], [0]], 1),
             ([[1, 0], [1]], 1),
             ([[1, 1], [1]], 1),
             ([[2, 0], [0]], 1),
             ([[2, 0], [1]], 1)]
        """
        if self._sageconstructor_overridden is None:
            g = (x for n in self._range for x in self._sageconstructor(n))
        else:
            g = self._sageconstructor_overridden

        return [(x, statistic(x)) for (x,_) in zip(g, xrange(max_values))]

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
        return self._id

    def id_str(self):
        r"""
        Return the FindStat identifier of the collection.

        OUTPUT:

        The FindStat identifier of the collection as a string.

        EXAMPLES::

            sage: from sage.databases.findstat import FindStatCollection
            sage: c = FindStatCollection("GelfandTsetlinPatterns")              # optional -- internet
            sage: c.id_str()                                                    # optional -- internet
            'Cc0018'
        """
        id = str(self.id())
        return 'Cc0000'[:-len(id)] + id

    def browse(self):
        r"""
        Open the FindStat web page of the collection in a browser.

        EXAMPLES::

            sage: from sage.databases.findstat import FindStatCollection
            sage: FindStatCollection("Permutations").browse()                   # optional -- webbrowser
        """
        webbrowser.open(FINDSTAT_URL + self._url_name)

    def to_size(self):
        r"""
        Return a function that returns the FindStat size of an object.

        OUTPUT:

        The function that produces the size as needed by the
        constructor of the collection.

        EXAMPLES::

            sage: from sage.databases.findstat import FindStatCollection
            sage: c = FindStatCollection("GelfandTsetlinPatterns")              # optional -- internet
            sage: c.to_size()(GelfandTsetlinPattern([[4, 1], [1]]))             # optional -- internet
            (2, 4)
        """
        return self._to_size


    def to_string(self):
        r"""
        Return a function that returns the FindStat normal
        representation given an object.

        OUTPUT:

        The function that produces the string representation as
        needed by the FindStat search webpage.

        EXAMPLES::

            sage: from sage.databases.findstat import FindStatCollection
            sage: p = Poset((range(3), [[0, 1], [1, 2]]))                       # optional -- internet
            sage: c = FindStatCollection("Posets")                              # optional -- internet
            sage: c.to_string()(p)                                              # optional -- internet
            '([(0, 2), (2, 1)], 3)'
        """
        return self._to_str

    def from_string(self):
        r"""
        Return a function that returns the object given the FindStat
        normal representation.

        OUTPUT:

        The function that produces the sage object given its FindStat
        normal representation as a string.

        EXAMPLES::

            sage: from sage.databases.findstat import FindStatCollection
            sage: c = FindStatCollection("Posets")                              # optional -- internet
            sage: p = c.from_string()('([(0, 2), (2, 1)], 3)')                  # optional -- internet
            sage: p.cover_relations()                                           # optional -- internet
            [[0, 2], [2, 1]]
        """
        return self._from_str

    def __repr__(self):
        r"""
        Return the representation of the FindStat collection.

        OUTPUT:

        The representation, including the identifier and the name.

        EXAMPLES::

            sage: from sage.databases.findstat import FindStatCollection
            sage: FindStatCollection("Binary trees")                            # optional -- internet
            Cc0010: Binary trees
        """
        return "%s: %s" %(self.id_str(), self._name_plural)

    def name(self):
        r"""
        Return the name of the FindStat collection.

        OUTPUT:

        The name of the FindStat collection, in singular.

        EXAMPLES::

            sage: from sage.databases.findstat import FindStatCollection
            sage: FindStatCollection("Binary trees").name()                     # optional -- internet
            u'Binary tree'
        """
        return self._name

class FindStatMap(SageObject):
    r"""
    The class of FindStat maps.

    INPUT:

    - ``entry`` -- a string giving the FindStat name of the map, or
      an integer giving its id.
    """
    _findstat_maps = []

    def __init__(self, entry):
        if not self._findstat_maps:
            self._findstat_maps += json.load(urlopen(FINDSTAT_URL_DOWNLOADS_MAPS))

        bad = True
        if isinstance(entry, (str, unicode)):
            # find by name in _findstat_maps
            for c in FindStatMap._findstat_maps:
                if entry.upper() == c[FINDSTAT_MAP_NAME].upper():
                    self._map = c
                    bad = False
                    break

        elif isinstance(entry, (int, Integer)):
            # find by id in _findstat_maps
            for c in FindStatMap._findstat_maps:
                if entry == c[FINDSTAT_MAP_IDENTIFIER]:
                    self._map = c
                    bad = False
                    break

        if bad:
            raise ValueError, "Could not find FindStat map for " + str(entry)

    def id(self):
        r"""
        Return the FindStat identifier of the map.

        OUTPUT:

        The FindStat identifier of the map as an integer.

        EXAMPLES::

            sage: m = findstat(lambda pi: pi.length(), "Permutations")[1][1][0] # optional -- internet
            sage: m.id()                                                        # optional -- internet
            62
        """
        return self._map[FINDSTAT_MAP_IDENTIFIER]

    def id_str(self):
        r"""
        Return the FindStat identifier of the map.

        OUTPUT:

        The FindStat identifier of the map as a string.

        EXAMPLES::

            sage: m = findstat(lambda pi: pi.length(), "Permutations")[1][1][0] # optional -- internet
            sage: m.id_str()                                                    # optional -- internet
            'Mp00062'
        """

        id = str(self.id())
        return 'Mp00000'[:-len(id)] + id

    def __repr__(self):
        r"""
        Return the representation of the FindStat map.

        EXAMPLES::

            sage: from sage.databases.findstat import FindStatMap
            sage: FindStatMap(71)                                               # optional -- internet
            Mp00071: descent composition
        """
        return "%s: %s" %(self.id_str(), self._map[FINDSTAT_MAP_NAME])

    def __eq__(self, other):
        """
        TESTS::

            sage: from sage.databases.findstat import FindStatMap
            sage: FindStatMap(71) == FindStatMap(71)                            # optional -- internet
            True

            sage: FindStatMap(62) == FindStatMap(71)                            # optional -- internet
            False
        """
        return self.id() == other.id()

    def __ne__(self, other):
        """
        TESTS::

            sage: from sage.databases.findstat import FindStatMap
            sage: FindStatMap(71) != FindStatMap(71)                            # optional -- internet
            False

            sage: FindStatMap(62) != FindStatMap(71)                            # optional -- internet
            True
        """
        return not self.__eq__(other)


    def name(self):
        r"""
        Return the FindStat name of the map.

        OUTPUT:

        The name of the map as a string, as used by FindStat.

        EXAMPLES::

            sage: m = findstat(lambda pi: pi.length(), "Permutations")[1][1][0] # optional -- internet
            sage: m.name()                                                      # optional -- internet
            u'inversion-number to major-index bijection'
        """
        return self._map[FINDSTAT_MAP_NAME]

    def description(self):
        r"""
        Return the FindStat description of the map.

        OUTPUT:

        The description as a string.

        EXAMPLES::

            sage: m = findstat(lambda pi: pi.length(), "Permutations")[1][1][0] # optional -- internet
            sage: print m.description()                                         # optional -- internet,random
            Let $\sigma \in \mathcal{S}_n$ be a permutation.
            <BLANKLINE>
            Maps $\sigma$ to the permutation $\tau$ such that the major code of $\tau$ is given by the Lehmer code of $\sigma$.
            <BLANKLINE>
            In particular, the number of inversions of $\sigma$ equals the major index of $\tau$.
            <BLANKLINE>
            EXAMPLES:
            <BLANKLINE>
            $[3,4,1,2] \mapsto [3,1,4,2]$
        """
        return self._map[FINDSTAT_MAP_DESCRIPTION]

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
        return FindStatCollection(self._map[FINDSTAT_MAP_DOMAIN])

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
        return FindStatCollection(self._map[FINDSTAT_MAP_CODOMAIN])

    def code(self):
        r"""
        Return the code associated with the map.

        OUTPUT:

        A string.

        EXAMPLES::

            sage: from sage.databases.findstat import FindStatMap               # optional -- internet
            sage: print FindStatMap(71).code()                                  # optional -- internet
            def descents_composition(elt):
                if len(elt) == 0:
                    return Composition([])
                d = [-1] + elt.descents() + [len(elt)-1]
                return Composition([ d[i+1]-d[i] for i in range(len(d)-1)])
        """
        return self._map[FINDSTAT_MAP_CODE]

    def code_name(self):
        r"""
        Return the name of the function defined by :meth:`code`.

        OUTPUT:

        A string.

        EXAMPLES::

            sage: from sage.databases.findstat import FindStatMap               # optional -- internet
            sage: print FindStatMap(71).code_name()                             # optional -- internet
            descents_composition
        """
        return self._map[FINDSTAT_MAP_CODE_NAME]

findstat = FindStat()
