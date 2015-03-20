r"""
This is the Sage interface for the FindStat project at www.findstat.org.

TODO:

- in FindStat, it would be nice if the post keys would coincide with
  the json keys.  Below is a long list of constants, a good scheme
  should be worked out.  Possibly it is best to use the database
  field names themselves everywhere - also in the json templates.

To run all doctests, run

sage -t --optional=sage,internet path/to/interface.sage
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

FINDSTAT_URL = 'http://www.findstat.org'
FINDSTAT_URL_DOWNLOADS = 'http://downloads.findstat.org'

FINDSTAT_URL_RESULT    = FINDSTAT_URL + "/StatisticFinder/Result/"
FINDSTAT_URL_LOGIN     = FINDSTAT_URL + "/StatisticFinder?action=login"
FINDSTAT_URL_NEW       = FINDSTAT_URL + '/StatisticsDatabase/NewStatistic/'
FINDSTAT_URL_EDIT      = FINDSTAT_URL + '/StatisticsDatabase/EditStatistic/'
FINDSTAT_URL_BROWSE    = FINDSTAT_URL + '/StatisticsDatabase/'

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
# we might want to use these, too: MapSageName, MapCode

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
FINDSTAT_SEPARATOR_REFERENCES = "\r\n\r\n"

######################################################################

# the format string for using POST
FINDSTAT_POST_HEADER = """
<script src='http://www.google.com/jsapi'></script>
<script>
    google.load('jquery', '1.3.2');
</script>

<script>
    $(document).ready(function() {$("#form").submit(); });
</script>
"""

FINDSTAT_NEWSTATISTIC_FORM_HEADER = "<form id='form' name='NewStatistic' action='%s' enctype='multipart/form-data' method='post' />"
FINDSTAT_NEWSTATISTIC_FORM_FORMAT = "<input type='hidden' name='%s' value='%s' />"
FINDSTAT_NEWSTATISTIC_FORM_FOOTER = "</form>"

######################################################################
class FindStat():
    r"""
    The Combinatorial Statistic Finder.

    ``FindStat`` is a class representing results of queries, edits
    and new submissions to the FindStat database.

    In principle there are two different results of a query:

    - a FindStat ID yields a single statistic

    - statistic data yield a list of triples (statistic, maps, quality)

    We want a single entry point for both and we want to be able to
    prepare a new statistic for submission by using the second type
    of query.

    """
    def __init__(self):
        # a cache from integers to ``FindStatStatistic`` instances
        # that avoids retrieving the same statistic over and over
        # again
        self._statistic_cache = dict()

        # user credentials if provided
        self._user_name  = ""
        self._user_email = ""

    def __call__(self, query, collection=None, depth=2, max_values=FINDSTAT_MAX_VALUES):
        r"""
        Returns an instance of a ``FindStatStatistic``.  This should
        be the only way to access ``FindStatStatistic``.  We do the
        preprocessing of the data here, and call the appropriate
        method of ``FindStatStatistic`` to launch the query.

        INPUT:

        - an integer or a string representing a valid FindStat ID
        (e.g. 45 or 'St000045').  ``collection`` should be ``None``,
        ``depth`` and ``max_values`` are ignored.

        - a list of pairs of the form (object, value), or a
        dictionary from sage objects to integer values.
        ``collection`` should be ``None``, ``depth`` and
        ``max_values`` are passed to the finder.

        - a list of pairs of the form (list of objects, list of
        values), or a single pair of the form (list of objects, list
        of values).  In each pair there should be as many objects as
        values.  ``collection`` should be ``None``, ``depth`` and
        ``max_values`` are passed to the finder.

        - a callable and a collection.  The callable is used to
        generate ``max_values`` (object, value) pairs.  The number of
        terms generated may also be controlled by passing an iterable
        collection, such as Permutations(3).  ``depth`` and
        ``max_values`` are passed to the finder.

        OUTPUT:

        - the FindStat ID together with its name

        - a list of triples, each consisting of

          - the statistic

          - a list of strings naming certain maps

          - a number which says how many of the values submitted agree
            with the values in the database, when applying the maps in
            the given order to the object and then computing the
            statistic on the result.

        EXAMPLES::

        A particular statistic can be called by its St-number or number::

            sage: findstat('St000045')                   # optional -- internet
            St000045: The number of linear extensions of the tree.

            sage: findstat(3)                            # optional -- internet
            St000003: The number of [[/StandardTableaux|standard Young tableaux]] of the partition.

        The database can be searched by providing a dictionary::

            sage: stat = {pi: pi.length() for pi in Permutations(4)}
            sage: findstat(stat)                         # optional -- internet
            0: (St000018: The [[/Permutations/Inversions|number of inversions]] of a permutation., [], 24)
            1: (St000004: The [[/Permutations/Descents-Major|major index]] of a permutation., [Mp00062: inversion-number to major-index bijection], 24)
            2: (St000067: The inversion number of the alternating sign matrix., [Mp00063: to alternating sign matrix], 24)
            3: (St000008: The major index of the composition., [Mp00062: inversion-number to major-index bijection, Mp00071: descent composition], 24)
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
                                         collection = collection,
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
            The Combinatorial Statistics Finder (http://www.findstat.org)

        """
        return "The Combinatorial Statistics Finder (%s)" % FINDSTAT_URL

    def browse(self):
        r"""
        Open the FindStat web page in a browser.
        """
        webbrowser.open(FINDSTAT_URL)

    def set_user(self, name=None, email=None):
        r"""
        Sets the user for this session.
        """
        if not isinstance(name, str):
            raise ValueError("The given name is not a string")
        if not isinstance(email, str):
            raise ValueError("The given email address is not a string")
        self._user_name  = name
        self._user_email = email

    def login(self):
        r"""
        Open the login page in a browser.
        """
        webbrowser.open(FINDSTAT_URL_LOGIN)

    ######################################################################

    def _statistic(self, id):
        r"""

        INPUT:
        - id, an integer designating the FindStat id of a statistic

        OUTUT:
        - a FindStatStatistic instance

        TODO:
        - this method caches the statistics.  It may make sense to
        provide a method that clears the cache, or reloads a single
        statistic
        """
        if id > 0:
            if id not in self._statistic_cache.keys():
                self._statistic_cache[id] = FindStatStatistic(id)._find_by_id()
            return self._statistic_cache[id]
        else:
            raise ValueError("The statistic identifier must be at least 1.")

######################################################################

# we use read-only properties, so we need a new-style class here.
class FindStatStatistic(SageObject):
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
                return "%s(modified): %s" % (self.id_str, self.name)
            else:
                return "%s: %s" % (self.id_str, self.name)

        elif self._query == "data":
            if len(self._result) == 0:
                return "a new statistic on " + self.collection.__repr__()
            else:
                return self._result.__repr__()

        else:
            raise ValueError("self._query should be either 'ID' or 'data', but is %s" %self._query)

    def __eq__(self, other):
        """
        TODO: this is *very* rudimentary
        """
        if self._query == "ID" and other._query == "ID":
            if self._modified or other._modified:
                return False
            else:
                return self.id == other.id
        elif self._query == "data" and other._query == "data":
            if self._modified or other._modified:
                return False
            else:
                return self.data == other.data
        else:
            return False

    ######################################################################
    # query = "ID"
    ######################################################################

    def _find_by_id(self):
        r"""
        expects that _id is a valid identifier

        OUTPUT: self

        TESTS::

        sage: findstat(999999)
        Traceback (most recent call last):
        ...
        ValueError: St999999 is not a FindStat statistic identifier.
        """
        self._query = "ID"

        # get the database entry from FindStat
        filename = "/statistics/" + self.id_str + ".json"
        url = FINDSTAT_URL_DOWNLOADS + filename
        _ = verbose("Fetching URL %s ..." %url, caller_name='FindStat')
        try:
            self._raw = json.load(urlopen(url), object_pairs_hook=OrderedDict)
        except HTTPError, error:
            if error.code == 404:
                raise ValueError("%s is not a FindStat statistic identifier." %self.id_str)
            else:
                raise

        self._description = self._raw[FINDSTAT_STATISTIC_DESCRIPTION]
        self._references = FancyTuple(self._raw[FINDSTAT_STATISTIC_REFERENCES].split(FINDSTAT_SEPARATOR_REFERENCES))
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

        expects that data is a list of pairs of the form (list of
        objects, list of values), each containing as many values as
        objects, and that collection is appropriately set

        OUTPUT: self

        TESTS::

        sage: findstat(lambda x: 1, "Permutations", depth=100)
        Traceback (most recent call last):
        ...
        ValueError: The depth must be a non-negative integer less than or equal to 5.
        """
        self._query = "data"
        # FindStat allows to search for at most FINDSTAT_MAX_VALUES values.
        data = self._data[:min(max_values, FINDSTAT_MAX_VALUES)]

        # for distribution searches, we cannot simply restrict to the
        # first values, so we raise an error
        try:
            assert sum(len(objects) for objects, values in data) <= FINDSTAT_MAX_VALUES
        except:
            raise NotImplementedError("It is unclear how to restrict distribution searches to at most %s values" %FINDSTAT_MAX_VALUES)

        url = FINDSTAT_URL_RESULT + self._collection.url_name() + "/"

        to_str = self._collection.to_string()
        stat = [(map(to_str, keys), str(values)[1:-1]) for (keys, values) in data]

        stat_str = join([join(keys, "\n") + "\n====> " + values for (keys, values) in stat], "\n")
        _ = verbose("Sending data to Findstat %s" %stat_str, caller_name='FindStat')

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

    ######################################################################
    # properties
    ######################################################################

    # read-only:
    @property
    def id(self):
        return self._id

    @property
    def id_str(self):
        id = str(self.id)
        return 'St000000'[:-len(id)] + id

    @property
    def data(self):
        return self._data

    @property
    def modified(self):
        return self._modified

    @property
    def collection(self):
        return self._collection

    @property
    def function(self):
        return self._function

    @property
    def first_terms(self):
        return self._first_terms

    @property
    def first_terms_str(self):
        if self._first_terms != None:
            to_str = self._collection.to_string()
            return join([to_str(key) + " => " + str(val)
                         for (key, val) in self._first_terms], "\r\n")
        else:
            return ""

    # writable:
    @property
    def description(self):
        return self._description

    @description.setter
    def description(self, value):
        if value != self._description:
            self._modified = True
            self._description = value

    @property
    def name(self):
        return self.description.partition(FINDSTAT_SEPARATOR_NAME)[0]

    @name.setter
    def name(self, value):
        raise ValueError("The name is the first line of the description.  Please modify the description instead.")

    @property
    def references(self):
        return self._references

    @references.setter
    def references(self, value):
        if value != self._references:
            self._modified = True
            self._references = value

    @property
    def code(self):
        return self._code

    @code.setter
    def code(self, value):
        if value != self._code:
            self._modified = True
            self._code = value

    ######################################################################
    # browse current statistic
    ######################################################################

    def browse(self):
        r"""
        Open the FindStat web page associated to the statistic ``self`` in a browser.

        EXAMPLES::

            sage: findstat(45).browse()                            # optional -- internet, webbrowser
        """
        if self._query == "ID":
            webbrowser.open(FINDSTAT_URL_BROWSE + self.id_str)
        else:
            raise NotImplementedError("Would be nice to show the result of the query in the webbrowser.")

    ######################################################################
    # submit current (possibly incompletely defined) statistic
    ######################################################################

    def submit(self, max_values=FINDSTAT_MAX_SUBMISSION_VALUES):
        r"""
        Open the FindStat web page for edits in a browser.

        TODO:

        - decide whether we want to somehow take into account when
          there is a statistic that matches the data with depth 0

        """
        # if the statistic is given as a function, and we have a new
        # statistic then update first_terms

        # it is not clear whether we want to do this also for old statistics.
        if self.function and self.id == 0:
            self._first_terms = self.collection.first_terms(self.function,
                                                            max_values=max_values)

        args = dict()
        args[FINDSTAT_STATISTIC_IDENTIFIER]  = self.id
        args[FINDSTAT_STATISTIC_COLLECTION]  = str(self.collection.id())
        args[FINDSTAT_STATISTIC_DATA]        = self.first_terms_str
        args[FINDSTAT_STATISTIC_DESCRIPTION] = self.description
        args[FINDSTAT_STATISTIC_REFERENCES]  = join(self.references, FINDSTAT_SEPARATOR_REFERENCES)
        args[FINDSTAT_STATISTIC_CODE]        = self.code
        args[FINDSTAT_POST_SAGE_CELL]        = ""
        args[FINDSTAT_POST_EDIT]             = ""
        args[FINDSTAT_POST_AUTHOR]           = findstat._user_name
        args[FINDSTAT_POST_EMAIL]            = findstat._user_email

        assert set(args.keys()) == FINDSTAT_EDIT_FIELDS, "It appears that the list of required post variables for editing a statistic has changed.  Please update FindStatStatistic.submit()."

        _ = verbose("Submitting arguments %s" %args, caller_name='FindStat')

        # write the file
        f = tempfile.NamedTemporaryFile(delete=False)
        _ = verbose("Created temporary file %s" %f.name, caller_name='FindStat')
        f.write(FINDSTAT_POST_HEADER)
        if self.id == 0:
            f.write(FINDSTAT_NEWSTATISTIC_FORM_HEADER %FINDSTAT_URL_NEW)
        else:
            f.write(FINDSTAT_NEWSTATISTIC_FORM_HEADER %(FINDSTAT_URL_EDIT+self.id_str))
        for key, value in args.iteritems():
            f.write(FINDSTAT_NEWSTATISTIC_FORM_FORMAT %(key, value))
        f.write(FINDSTAT_NEWSTATISTIC_FORM_FOOTER)
        f.close()
        _ = verbose("Opening file with webbrowser", caller_name='FindStat')
        _ = webbrowser.open(f.name)

        _ = verbose("Waiting a little before deleting the temporary file", caller_name='FindStat')
        time.sleep(1)

        f.unlink(f.name)

    # editing and submitting is really the same thing
    edit = submit

class FindStatCollection():
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
        13: [None, None, None, Core,                  lambda x: Cores(*x),     None,
             lambda x: x.k(),
             lambda X: "( "+X._repr_()+", "+str(X.k())+" )",
             lambda x: (lambda pi, k: Core(pi, k))(*literal_eval(x))],
        5:  [None, None, None, DyckWord,              DyckWords,               None,
             lambda x: x.length(),
             lambda x: str(list(DyckWord(x))),
             lambda x: DyckWord(literal_eval(x))],
        22: [None, None, None, CartanType_abstract,   _finite_irreducible_cartan_types_by_rank,
             None,
             lambda x: x.rank(),
             str,
             lambda x: CartanType(*literal_eval(str(x)))],
        18: [None, None, None, GelfandTsetlinPattern, lambda x: GelfandTsetlinPatterns(*x),
             None,
             len,
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

    def __init__(self, entry):
        r"""
        Initializes a FindStat collection from

        * a string, eg. 'Dyck paths' or "DyckPaths"
        * a sage object

        TODO:

        * a FindStat collection should also work.

        TESTS::

            sage: FindStatCollection('Dyck paths')                 # optional -- internet
            Cc0005: Dyck paths

        """
        if not FindStatCollection._findstat_collections.values()[0][0]:
            for j in json.load(urlopen(FINDSTAT_URL_DOWNLOADS + "/collections.json")):
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
            for (id, c) in FindStatCollection._findstat_collections.iteritems():
                if isinstance(entry, c[3]) or entry == c[4]:
                    initialize_with(id, c)
                    bad = False
                    break

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
        return self.id == other.id

    def first_terms(self, statistic, max_values=FINDSTAT_MAX_SUBMISSION_VALUES):
        r"""

        given a statistic, compute the first few terms and return them as a list of pairs
        """
        if self._sageconstructor_overridden is None:
            g = (x for n in self._range for x in self._sageconstructor(n))
        else:
            g = self._sageconstructor_overridden

        return [(x, statistic(x))
                for (x,_) in zip(g, xrange(min(max_values, FINDSTAT_MAX_SUBMISSION_VALUES)))]

    def id(self):
        return self._id

    def id_str(self):
        id = str(self.id())
        return 'Cc0000'[:-len(id)] + id

    def url_name(self):
        return self._url_name

    def to_size(self):
        r"""

        returns a function that returns the FindStat size of a sage object
        """
        return self._to_size


    def to_string(self):
        r"""

        returns a function that converts a sage object into a FindStat string
        """
        return self._to_str

    def from_string(self):
        r"""
        Returns the sage constructor corresponding to the combinatorial collection, or `None`.

        OUTPUT:

        - string or `None`.

        """
        return self._from_str

    def __repr__(self):
        r"""
        Prints the name of the collection

        OUTPUT:

        - string.

        EXAMPLES::

            sage: FindStatCollection("Binary trees")                      # optional -- internet
            Cc0010: Binary trees
        """
        return "%s: %s" %(self.id_str(), self._name_plural)

    def name(self):
        r"""
        returns the name of the collection

        OUTPUT:

        - string.

        EXAMPLES::

            sage: FindStatCollection("Binary trees")                      # optional -- internet
            Cc0010: Binary trees
        """
        return self._name

class FindStatMap():
    r"""
    The class of FindStat maps.

    """
    _findstat_maps = []

    def __init__(self, entry):
        r"""
        Initializes a FindStat map
        """
        if not self._findstat_maps:
            self._findstat_maps += json.load(urlopen(FINDSTAT_URL_DOWNLOADS + "/maps.json"))

        bad = True
        if isinstance(entry, (str, unicode)):
            # find by name in _findstat_collections
            for c in FindStatMap._findstat_maps:
                if entry.upper() == c[FINDSTAT_MAP_NAME].upper():
                    self._map = c
                    bad = False
                    break

        elif isinstance(entry, (int, Integer)):
            # find by id in _findstat_collections
            for c in FindStatMap._findstat_maps:
                if entry == c[FINDSTAT_MAP_IDENTIFIER]:
                    self._map = c
                    bad = False
                    break

        if bad:
            raise ValueError, "Could not find FindStat map for " + str(entry)

    def id(self):
        return self._map[FINDSTAT_MAP_IDENTIFIER]

    def id_str(self):
        id = str(self.id())
        return 'Mp00000'[:-len(id)] + id


    def __repr__(self):
        r"""
        Prints the name of the map
        """
        return "%s: %s" %(self.id_str(), self._map[FINDSTAT_MAP_NAME])

    def name(self):
        return self._map[FINDSTAT_MAP_NAME]

    def description(self):
        return self._map[FINDSTAT_MAP_DESCRIPTION]

    def domain(self):
        return FindStatCollection(self._map[FINDSTAT_MAP_DOMAIN])

    def codomain(self):
        return FindStatCollection(self._map[FINDSTAT_MAP_CODOMAIN])


findstat = FindStat()
