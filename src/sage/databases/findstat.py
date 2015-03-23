r"""
FindStat - the Combinatorial Statistic Finder.

You can use the sage interface to FindStat to:

    - identify a combinatorial statistic from the values on a few small objects.
    - obtain more terms, formulae, references, etc. for a given statistic.
    - edit statistics and submit new statistics

To access the database, use :class:`findstat<FindStat>`::

    sage: findstat
    The Combinatorial Statistic Finder (http://www.findstat.org/)

AUTHORS:

- Martin Rubey (2015): initial version.


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

FINDSTAT_URL = 'http://www.findstat.org/'
FINDSTAT_URL_DOWNLOADS = 'http://downloads.findstat.org'

FINDSTAT_URL_RESULT    = FINDSTAT_URL + "StatisticFinder/Result/"
FINDSTAT_URL_LOGIN     = FINDSTAT_URL + "StatisticFinder?action=login"
FINDSTAT_URL_NEW       = FINDSTAT_URL + 'StatisticsDatabase/NewStatistic/'
FINDSTAT_URL_EDIT      = FINDSTAT_URL + 'StatisticsDatabase/EditStatistic/'
FINDSTAT_URL_BROWSE    = FINDSTAT_URL + 'StatisticsDatabase/'

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

    :class:`FindStat` is a class representing results of queries to
    the FindStat database.  This class is also the entry point to
    edit statistics and new submissions.  Use the shorthand
    :class:`findstat<FindStat>` to call it.

    INPUT:

    - an integer or a string representing a valid FindStat ID
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

    an instance of a :class:`FindStatStatistic`, represented by

    - the FindStat ID together with its name, or

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

        sage: findstat('St000045')                   # optional -- internet
        St000045: The number of linear extensions of the tree.

        sage: findstat(3)                            # optional -- internet
        St000003: The number of [[/StandardTableaux|standard Young tableaux]] of the partition.

    The database can be searched by providing a list of pairs::

        sage: q = findstat([(pi, pi.length()) for pi in Permutations(4)]); q # optional -- internet
        0: (St000018: The [[/Permutations/Inversions|number of inversions]] of a permutation., [], 24)
        1: (St000004: The [[/Permutations/Descents-Major|major index]] of a permutation., [Mp00062: inversion-number to major-index bijection], 24)
        ...

    or a dictionary::

        sage: p = findstat({pi: pi.length() for pi in Permutations(4)}); p # optional -- internet
        0: (St000018: The [[/Permutations/Inversions|number of inversions]] of a permutation., [], 24)
        1: (St000004: The [[/Permutations/Descents-Major|major index]] of a permutation., [Mp00062: inversion-number to major-index bijection], 24)
        ...

    Note however, that the results of these two queries are not
    necessarily the same, because we compare queries by the data
    sent, and the ordering of the data might be different::

        sage: p == q                                           # optional -- internet
        False

    Another possibility is to send a function and a collection.  In
    this case, the function is applied to the first few objects of
    the collection::

        sage: findstat(lambda pi: pi.length(), "Permutations") # optional -- internet
        0: (St000018: The [[/Permutations/Inversions|number of inversions]] of a permutation., [], 200)
        ...

    To search for a distribution, send a list of lists, or a single pair::

        sage: findstat((Permutations(4), [pi.length() for pi in Permutations(4)])) # optional -- internet
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
            The Combinatorial Statistic Finder (http://www.findstat.org/)

        """
        return "The Combinatorial Statistic Finder (%s)" % FINDSTAT_URL

    def browse(self):
        r"""
        Open the FindStat web page in a browser.

        EXAMPLES::

            sage: findstat.browse()                              # optional -- internet, webbrowser
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

        - a :class:`FindStatStatistic` instance.

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

# we use read-only properties, so we need a new-style class here.
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
        .. TODO::

            this is *very* rudimentary
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
        Open the FindStat web page of the statistic in a browser.

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
        Open the FindStat web page for editing the statistic ``self`` in a browser.

        .. TODO::

            decide whether we want to somehow take into account when
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

        sage: FindStatCollection("Permutations")                         # optional -- internet
        Cc0001: Permutations

    TESTS::

        # when one of these tests does not pass, there is probably a new collection to be added.
        sage: cdata = FindStatCollection._findstat_collections.values()  # optional -- internet
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

        .. TODO::

            a FindStat collection should also work.

        EXAMPLES::

            sage: from sage.databases.findstat import FindStatCollection
            sage: FindStatCollection("Dyck paths")                 # optional -- internet
            Cc0005: Dyck paths

            sage: FindStatCollection(5)                            # optional -- internet
            Cc0005: Dyck paths

            sage: FindStatCollection(DyckWord([1,0,1,0]))          # optional -- internet
            Cc0005: Dyck paths

            sage: FindStatCollection(DyckWords(2))                 # optional -- internet
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
            sage: FindStatCollection("Permutations") == FindStatCollection("Permutations")
            True

            sage: FindStatCollection("Permutations") != FindStatCollection("Integer Partitions")
            True
        """
        return self.id() == other.id()

    def in_range(self, element):
        r"""
        Check whether an element of the collection is in FindStat's precomputed range.

        INPUT:

        - ``element`` -- a sage object that belongs to the collection

        OUTPUT:

        - ``True``, if ``element`` is used by the FindStat search
          engine, and ``False`` if it is ignored.

        """
        n = self.to_size()(element)
        return element in self._sageconstructor(n)

    def first_terms(self, statistic, max_values=FINDSTAT_MAX_SUBMISSION_VALUES):
        r"""
        Compute the first few terms of the given statistic.

        INPUT:

        - ``statistic`` -- a callable.

        - ``max_values`` -- the number of terms to compute at most.

        OUTPUT:

        - a list of pairs of the form (object, value).

        """
        if self._sageconstructor_overridden is None:
            g = (x for n in self._range for x in self._sageconstructor(n))
        else:
            g = self._sageconstructor_overridden

        return [(x, statistic(x)) for (x,_) in zip(g, xrange(max_values))]

    def id(self):
        r"""
        Return the FindStat identifier of the collection as an integer.
        """
        return self._id

    def id_str(self):
        r"""
        Return the FindStat identifier of the collection as a string.
        """
        id = str(self.id())
        return 'Cc0000'[:-len(id)] + id

    def browse(self):
        r"""
        Open the FindStat web page of the collection in a browser.

        EXAMPLES::

            sage: from sage.databases.findstat import FindStatCollection
            sage: FindStatCollection("Permutations").browse()            # optional -- internet, webbrowser
        """
        webbrowser.open(FINDSTAT_URL + self._url_name)

    def to_size(self):
        r"""
        Return a function that returns the FindStat size of an object.

        OUTPUT:

        The function that produces the size as needed by the
        constructor of the collection.

        TESTS::

            sage: from sage.databases.findstat import FindStatCollection
            sage: c = FindStatCollection("GelfandTsetlinPatterns")                 # optional -- internet
            sage: c.to_size()(GelfandTsetlinPattern([[4, 1], [1]]))                # optional -- internet
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

        TESTS::

            sage: from sage.databases.findstat import FindStatCollection
            sage: p = Poset((range(3), [[0, 1], [1, 2]]))
            sage: c = FindStatCollection("Posets")                                 # optional -- internet
            sage: c.to_string()(p)                                                 # optional -- internet
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

        TESTS::

            sage: from sage.databases.findstat import FindStatCollection
            sage: c = FindStatCollection("Posets")                                 # optional -- internet
            sage: p = c.from_string()('([(0, 2), (2, 1)], 3)')                     # optional -- internet
            sage: p.cover_relations()                                              # optional -- internet
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
            sage: FindStatCollection("Binary trees")                      # optional -- internet
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
            sage: FindStatCollection("Binary trees").name()               # optional -- internet
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
            self._findstat_maps += json.load(urlopen(FINDSTAT_URL_DOWNLOADS + "/maps.json"))

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

        The identifier as an integer.
        """
        return self._map[FINDSTAT_MAP_IDENTIFIER]

    def id_str(self):
        r"""
        Return the FindStat identifier of the map.

        OUTPUT:

        The identifier as a string.
        """

        id = str(self.id())
        return 'Mp00000'[:-len(id)] + id

    def __repr__(self):
        r"""
        Return the representation of the FindStat map.
        """
        return "%s: %s" %(self.id_str(), self._map[FINDSTAT_MAP_NAME])

    def name(self):
        r"""
        Return the FindStat name of the map.

        OUTPUT:

        The name as a string.
        """

        return self._map[FINDSTAT_MAP_NAME]

    def description(self):
        r"""
        Return the FindStat description of the map.

        OUTPUT:

        The description as a string.
        """
        return self._map[FINDSTAT_MAP_DESCRIPTION]

    def domain(self):
        r"""
        Return the FindStat collection which is the domain of the map.

        OUTPUT:

        The domain of the map as a :class:`FindStatCollection`.

        """
        return FindStatCollection(self._map[FINDSTAT_MAP_DOMAIN])

    def codomain(self):
        r"""
        Return the FindStat collection which is the codomain of the map.

        OUTPUT:

        The codomain of the map as a :class:`FindStatCollection`.

        """
        return FindStatCollection(self._map[FINDSTAT_MAP_CODOMAIN])


findstat = FindStat()
